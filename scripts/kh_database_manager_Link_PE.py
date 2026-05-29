import sqlite3
import pandas as pd
import numpy as np
import time
import json
import ast
import argparse
from spectral_kh import KhovanovEngine

# --- CONFIGURATION ---
DB_PATH = "databases/link_research.db"
CSV_PATH = "data/linkinfo.csv"
DENSE_THRESHOLD = 4000


def init_db(conn):
    """Initialize the link database schema if it does not already exist."""
    cursor = conn.cursor()

    cursor.execute(
        '''
        CREATE TABLE IF NOT EXISTS links (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            pd_code TEXT,
            khovanov_polynomial TEXT
        )
        '''
    )

    cursor.execute(
        '''
        CREATE TABLE IF NOT EXISTS bidegrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            link_id INTEGER,
            h INTEGER,
            q INTEGER,
            dimension INTEGER,
            betti INTEGER,
            smallest_nonzero REAL,
            log_torsion REAL,
            is_truncated BOOLEAN,
            compute_time REAL,
            FOREIGN KEY (link_id) REFERENCES links (id),
            UNIQUE(link_id, h, q)
        )
        '''
    )

    cursor.execute(
        '''
        CREATE INDEX IF NOT EXISTS idx_bidegrees_link_id
        ON bidegrees(link_id)
        '''
    )

    conn.commit()


def parse_pd_code(raw_str):
    """Clean and parse PD codes from the CSV."""
    try:
        if pd.isna(raw_str):
            return None
        clean = str(raw_str).strip().replace(';', ',').replace('{', '[').replace('}', ']')
        clean = clean.replace(',]', ']').replace(', ]', ']')
        return ast.literal_eval(clean)
    except Exception:
        return None


def format_kh_poly(betti_dict):
    """Convert {(h, q): betti} into a polynomial string."""
    terms = []
    for (h, q), betti in sorted(betti_dict.items()):
        if betti == 1:
            terms.append(f"q^{q}t^{h}")
        else:
            terms.append(f"{betti}q^{q}t^{h}")
    return " + ".join(terms) if terms else "0"


def process_batch_parallel(args):
    conn = sqlite3.connect(DB_PATH, timeout=120)
    init_db(conn)

    try:
        df = pd.read_csv(CSV_PATH)
    except Exception as e:
        print(f"CRITICAL: Failed to read {CSV_PATH}. Error: {e}")
        conn.close()
        return

    # Partition rows across array tasks when launched in parallel.
    df = df[df.index % args.total_bins == args.array_id]

    for _, row in df.iterrows():
        link_name = str(row.iloc[0]).strip()
        pd_raw = row.iloc[1]

        pd_list = parse_pd_code(pd_raw)
        if not pd_list:
            continue

        cursor = conn.cursor()
        link_id = None

        try:
            # Ensure the master record exists.
            cursor.execute(
                "INSERT OR IGNORE INTO links (name, pd_code) VALUES (?, ?)",
                (link_name, json.dumps(pd_list)),
            )
            conn.commit()

            result = cursor.execute(
                "SELECT id, khovanov_polynomial FROM links WHERE name=?",
                (link_name,),
            ).fetchone()
            if result is None:
                raise RuntimeError(f"Could not retrieve database row for {link_name}")

            link_id, kh_poly = result

            # True resume logic: skip only if the final polynomial already exists.
            if kh_poly is not None:
                print(f"SKIP: {link_name} already completed.")
                continue

            # If a previous run left partial data, clear it and recompute cleanly.
            existing_rows = cursor.execute(
                "SELECT COUNT(*) FROM bidegrees WHERE link_id=?",
                (link_id,),
            ).fetchone()[0]
            if existing_rows > 0:
                print(f"RESET: clearing partial data for {link_name}")
                cursor.execute("DELETE FROM bidegrees WHERE link_id=?", (link_id,))
                conn.commit()

            print(f">>> NODE {args.array_id} | STARTING: {link_name}")

            engine = KhovanovEngine(pd_list)
            active_degrees = engine.get_active_bidegrees()
            betti_accumulator = {}

            for h, q in active_degrees:
                start_time = time.time()
                dim = len(engine.chain_groups.get((h, q), []))
                if dim == 0:
                    continue

                if dim < DENSE_THRESHOLD:
                    spec = engine.get_spectrum(h, q)
                    betti = int(np.sum(spec < 1e-9))
                    non_zero_eigs = spec[spec > 1e-9]
                    gap = float(non_zero_eigs[0]) if len(non_zero_eigs) > 0 else None
                    log_torsion = (
                        float(np.sum(np.log(non_zero_eigs))) if len(non_zero_eigs) > 0 else 0.0
                    )
                    is_truncated = False
                else:
                    betti, gaps = engine.get_low_spectrum_adaptive(h, q, tol=1e-9)
                    gap = float(gaps[0]) if len(gaps) > 0 else None
                    log_torsion = None
                    is_truncated = True

                if betti > 0:
                    betti_accumulator[(h, q)] = betti

                cursor.execute(
                    '''
                    INSERT OR REPLACE INTO bidegrees
                    (link_id, h, q, dimension, betti, smallest_nonzero, log_torsion, is_truncated, compute_time)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''',
                    (
                        link_id,
                        h,
                        q,
                        dim,
                        betti,
                        gap,
                        log_torsion,
                        is_truncated,
                        time.time() - start_time,
                    ),
                )
                conn.commit()

            poly_str = format_kh_poly(betti_accumulator)
            cursor.execute(
                "UPDATE links SET khovanov_polynomial = ? WHERE id = ?",
                (poly_str, link_id),
            )
            conn.commit()
            print(f"SUCCESS: {link_name}")

        except Exception as e:
            print(f"FAILED {link_name}: {e}")
            if link_id is not None:
                cursor.execute("DELETE FROM bidegrees WHERE link_id=?", (link_id,))
                cursor.execute(
                    "UPDATE links SET khovanov_polynomial = NULL WHERE id = ?",
                    (link_id,),
                )
                conn.commit()

    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--total-bins", type=int, default=1)
    args = parser.parse_args()
    process_batch_parallel(args)
