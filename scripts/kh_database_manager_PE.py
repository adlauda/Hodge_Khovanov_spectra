import sqlite3
import pandas as pd
import numpy as np
import time
import ast
import argparse
import os
from spectral_kh import KhovanovEngine

# --- CONFIGURATION ---
DB_PATH = "databases/knot_research.db"
CSV_PATH = "data/knotinfo.csv"

def init_db(conn):
    """Initialize the knot database with high-concurrency settings."""
    cursor = conn.cursor()
    # Use standard DELETE mode because NFS does not support WAL shared memory
    cursor.execute("PRAGMA journal_mode=DELETE;")
    cursor.execute("PRAGMA synchronous=NORMAL;")
    cursor.execute("PRAGMA busy_timeout=30000;") # Wait 30s if locked
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS knots (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            pd_code TEXT,
            khovanov_polynomial TEXT,
            turaev_genus INTEGER
        )
    ''')
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS bidegrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            knot_id INTEGER,
            h INTEGER,
            q INTEGER,
            dimension INTEGER,
            betti INTEGER,
            smallest_nonzero REAL,
            compute_time REAL,
            FOREIGN KEY (knot_id) REFERENCES knots (id),
            UNIQUE(knot_id, h, q)
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_bidegrees_knot_id ON bidegrees(knot_id)')
    conn.commit()

def parse_pd_code(raw_str):
    try:
        if pd.isna(raw_str): return None
        clean = str(raw_str).strip().replace(';', ',').replace('{', '[').replace('}', ']')
        clean = clean.replace(',]', ']').replace(', ]', ']')
        return ast.literal_eval(clean)
    except Exception: return None

def format_kh_poly(betti_dict):
    terms = []
    for (h, q), betti in sorted(betti_dict.items()):
        term = f"q^{q}t^{h}"
        if betti > 1: term = f"{betti}{term}"
        terms.append(term)
    return " + ".join(terms) if terms else "0"

def process_batch_parallel(args):
    # Use a longer timeout for the connection
    conn = sqlite3.connect(DB_PATH, timeout=60)
    init_db(conn)
    cursor = conn.cursor()
    
    try:
        df = pd.read_csv(CSV_PATH)
        
        # --- NEW FAST-FORWARD LOGIC ---
        # 1. Do a single bulk query to get all finished knots
        finished_df = pd.read_sql_query("SELECT name FROM knots WHERE khovanov_polynomial IS NOT NULL", conn)
        finished_set = set(finished_df['name'].tolist())
        
        # 2. Instantly drop all finished knots from the Pandas dataframe
        df = df[~df['Name'].isin(finished_set)].reset_index(drop=True)
        print(f"Node {args.array_id}: Fast-forwarded past {len(finished_set)} knots. {len(df)} left to distribute.")
        # ------------------------------
        
    except Exception as e:
        print(f"CSV/DB Read Error: {e}")
        return

    # Assign the remaining uncomputed knots to this specific worker node
    my_batch = df[df.index % args.total_bins == args.array_id]

    for _, row in my_batch.iterrows():
        knot_name = row['Name']
        pd_code = parse_pd_code(row['PD Notation'])
        if pd_code is None: continue

        # --- ROBUST DATABASE INITIALIZATION ---
        max_retries = 50
        already_done = None
        knot_id = None
        
        for attempt in range(max_retries):
            try:
                # 1. Initialize the knot entry
                cursor.execute("INSERT OR IGNORE INTO knots (name, pd_code) VALUES (?, ?)", 
                               (knot_name, str(pd_code)))
                conn.commit()
                
                knot_id = cursor.execute("SELECT id FROM knots WHERE name=?", (knot_name,)).fetchone()[0]
                
                # 2. Check if already finished
                already_done = cursor.execute("SELECT khovanov_polynomial FROM knots WHERE id=?", (knot_id,)).fetchone()[0]
                break # Success! Break out of the retry loop
                
            except sqlite3.OperationalError as e:
                if "locked" in str(e).lower() and attempt < max_retries - 1:
                    print(f"Node {args.array_id} hit NFS lock. Cycling connection (Attempt {attempt+1}/{max_retries})...")
                    time.sleep(10 + (os.getpid() % 15)) # Sleep 10 to 24 seconds
                    
                    # Force drop the ghost lock and reconnect
                    conn.close()
                    conn = sqlite3.connect(DB_PATH, timeout=60)
                    cursor = conn.cursor()
                    cursor.execute("PRAGMA journal_mode=DELETE;")
                    cursor.execute("PRAGMA busy_timeout=30000;")
                    continue
                else:
                    raise e # Let it crash if it fails 50 times

        if already_done:
            continue
            
        # --- DECOUPLED MATH ENGINE ---
        try:
            engine = KhovanovEngine(pd_code)
            t_genus = int(engine.get_turaev_genus())
            betti_accumulator = {}
            bidegree_results = [] # Local list to hold results
            
            # 3. Calculate all bidegrees (100% Math, NO Database Locks here!)
            for (h, q), dim in engine.chain_groups.items():
                start_time = time.time()
                
                betti, gaps = engine.get_low_spectrum_adaptive(h, q, tol=1e-7)
                gap = float(gaps[0]) if len(gaps) > 0 else None
                
                if betti > 0:
                    betti_accumulator[(h, q)] = betti
                
                # Append to our local list instead of writing to the DB
                bidegree_results.append((knot_id, h, q, len(dim), betti, gap, time.time() - start_time))
            
            # 4. LIGHTNING FAST BULK COMMIT
            poly_str = format_kh_poly(betti_accumulator)
            
            # Wrap the final write in our robust retry logic just in case
            for attempt in range(max_retries):
                try:
                    # Write all bidegrees at once in milliseconds
                    cursor.executemany('''
                        INSERT OR REPLACE INTO bidegrees 
                        (knot_id, h, q, dimension, betti, smallest_nonzero, compute_time)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', bidegree_results)
                    
                    # Old line:
                    # cursor.execute("UPDATE knots SET khovanov_polynomial = ? WHERE id = ?", (poly_str, knot_id))
                    
                    # New line:
                    cursor.execute("UPDATE knots SET khovanov_polynomial = ?, turaev_genus = ? WHERE id = ?", 
                                   (poly_str, t_genus, knot_id))
                    
                    conn.commit()
                    print(f"SUCCESS: {knot_name}")
                    break # Break out of retry loop on success
                    
                except sqlite3.OperationalError as e:
                    if "locked" in str(e).lower() and attempt < max_retries - 1:
                        time.sleep(2 + (os.getpid() % 5))
                        # Quick cycle for the bulk write
                        conn.close()
                        conn = sqlite3.connect(DB_PATH, timeout=60)
                        cursor = conn.cursor()
                        cursor.execute("PRAGMA journal_mode=DELETE;")
                        cursor.execute("PRAGMA busy_timeout=30000;")
                        continue
                    else:
                        raise e
                        
        except Exception as e:
            print(f"FAILED {knot_name}: {e}")
            try:
                cursor.execute("DELETE FROM bidegrees WHERE knot_id=?", (knot_id,))
                conn.commit()
            except: 
                pass
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--total-bins", type=int, default=1)
    args = parser.parse_args()
    process_batch_parallel(args)
