import sqlite3
import pandas as pd
import numpy as np
import time
import ast
import argparse
import re
from spectral_kh import KhovanovEngine

# =========================================================
# CONFIGURATION
# =========================================================
DB_PATH = "databases/twisted_unknot_research.db"
CSV_PATH = "data/TU_knots.csv"
DENSE_THRESHOLD = 4000
MAX_CROSSINGS = 13

# =========================================================
# Utilities
# =========================================================
def parse_pd_code(raw_str):
    """Clean and parse PD codes from CSV."""
    try:
        if pd.isna(raw_str):
            return None
        clean = str(raw_str).strip()
        clean = clean.replace(';', ',')
        clean = clean.replace('{', '[').replace('}', ']')
        clean = clean.replace(',]', ']').replace(', ]', ']')
        return ast.literal_eval(clean)
    except Exception:
        return None

def format_kh_poly(betti_dict):
    """Convert {(h,q): betti} dictionary into a polynomial string."""
    terms = []
    for (h, q), betti in sorted(betti_dict.items()):
        if betti == 1:
            terms.append(f"q^{q}t^{h}")
        else:
            terms.append(f"{betti}q^{q}t^{h}")
    return " + ".join(terms) if terms else "0"

def get_crossings_from_name(name):
    """
    Extract crossing number from names like TU_1, TU_2, ..., TU_13.
    """
    m = re.search(r'(\d+)$', str(name).strip())
    return int(m.group(1)) if m else None

def init_database():
    """Create database schema if it does not already exist."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS twisted_unknots (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE NOT NULL,
            crossings INTEGER,
            pd_code TEXT,
            khovanov_polynomial TEXT
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS bidegrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            twisted_unknot_id INTEGER NOT NULL,
            h INTEGER NOT NULL,
            q INTEGER NOT NULL,
            dimension INTEGER,
            betti INTEGER,
            smallest_nonzero REAL,
            log_torsion REAL,
            is_truncated BOOLEAN,
            compute_time REAL,
            FOREIGN KEY (twisted_unknot_id) REFERENCES twisted_unknots(id)
        )
    """)

    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_twisted_unknots_name
        ON twisted_unknots(name)
    """)

    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_bidegrees_twisted_unknot_id
        ON bidegrees(twisted_unknot_id)
    """)

    conn.commit()
    conn.close()

def load_twisted_unknots_from_csv():
    """
    Expect CSV columns:
      Name
      PD Notation
    """
    df = pd.read_csv(CSV_PATH)

    required_cols = {"Name", "PD Notation"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"{CSV_PATH} must contain columns {required_cols}. "
            f"Found columns: {list(df.columns)}"
        )

    df = df.copy()
    df["crossings"] = df["Name"].apply(get_crossings_from_name)
    df["pd_parsed"] = df["PD Notation"].apply(parse_pd_code)

    records = []
    for _, row in df.iterrows():
        name = str(row["Name"]).strip()
        crossings = row["crossings"]
        pd_code = row["pd_parsed"]

        if pd_code is None:
            print(f"Skipping {name}: could not parse PD notation")
            continue

        if crossings is None:
            print(f"Skipping {name}: could not determine crossing number from name")
            continue

        if crossings > MAX_CROSSINGS:
            continue

        records.append({
            "name": name,
            "crossings": int(crossings),
            "pd": pd_code
        })

    return records

# =========================================================
# Main processing logic
# =========================================================
def process_batch_parallel(args):
    conn = sqlite3.connect(DB_PATH, timeout=120)
    cur = conn.cursor()

    records = load_twisted_unknots_from_csv()

    # Partitioning logic for array jobs
    records = [rec for i, rec in enumerate(records) if i % args.total_bins == args.array_id]

    for rec in records:
        name = rec["name"]
        crossings = rec["crossings"]
        pd_list = rec["pd"]

        if not pd_list:
            continue

        # Ensure record exists
        cur.execute("""
            INSERT OR IGNORE INTO twisted_unknots (name, crossings, pd_code)
            VALUES (?, ?, ?)
        """, (name, crossings, str(pd_list)))
        conn.commit()

        twisted_id = cur.execute(
            "SELECT id FROM twisted_unknots WHERE name=?",
            (name,)
        ).fetchone()[0]

        # Resumption logic: skip if already computed
        existing_count = cur.execute(
            "SELECT COUNT(*) FROM bidegrees WHERE twisted_unknot_id=?",
            (twisted_id,)
        ).fetchone()[0]

        if existing_count > 0:
            continue

        print(f">>> NODE {args.array_id} | STARTING: {name}")

        try:
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
                    log_torsion = float(np.sum(np.log(non_zero_eigs))) if len(non_zero_eigs) > 0 else 0.0
                    is_truncated = False
                else:
                    betti, gaps = engine.get_low_spectrum_adaptive(h, q, tol=1e-9)
                    gap = float(gaps[0]) if len(gaps) > 0 else None
                    log_torsion = None
                    is_truncated = True

                if betti > 0:
                    betti_accumulator[(h, q)] = betti

                cur.execute("""
                    INSERT INTO bidegrees
                    (twisted_unknot_id, h, q, dimension, betti, smallest_nonzero, log_torsion, is_truncated, compute_time)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    twisted_id, h, q, dim, betti, gap, log_torsion,
                    is_truncated, time.time() - start_time
                ))
                conn.commit()

            poly_str = format_kh_poly(betti_accumulator)
            cur.execute("""
                UPDATE twisted_unknots
                SET khovanov_polynomial = ?
                WHERE id = ?
            """, (poly_str, twisted_id))
            conn.commit()

            print(f"SUCCESS: {name}")

        except Exception as e:
            print(f"FAILED {name}: {e}")
            try:
                cur.execute("DELETE FROM bidegrees WHERE twisted_unknot_id=?", (twisted_id,))
                conn.commit()
            except Exception as cleanup_error:
                print(f"CLEANUP FAILED for {name}: {cleanup_error}")

    conn.close()

# =========================================================
# Entry point
# =========================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--total-bins", type=int, default=1)
    args = parser.parse_args()

    init_database()
    process_batch_parallel(args)
