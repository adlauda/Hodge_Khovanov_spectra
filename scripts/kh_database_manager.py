"""
KHOVANOV HOMOLOGY DATABASE MANAGER
----------------------------------
This script manages the systematic calculation and storage of Khovanov homology 
invariants for a large-scale knot dataset. It interfaces with the KhovanovEngine
to compute Betti numbers, spectral gaps, and Reidemeister torsion.

The manager is designed for high-performance computing (HPC) environments,
featuring automatic resumption of interrupted jobs and an adaptive solver 
selection mechanism based on bidegree dimensionality.

Primary Metrics:
- Khovanov Polynomial: Human-readable representation of homology.
- Spectral Gap: The smallest non-zero eigenvalue (torsion candidate).
- Log-Torsion: The sum of logarithms of non-zero eigenvalues (volume invariant).
"""

import sqlite3
import pandas as pd
import numpy as np
import time
import json
import ast
import re
from spectral_kh import KhovanovEngine

# --- CONFIGURATION & THRESHOLDS ---
DB_PATH = "databases/knot_research.db"
CSV_PATH = "data/knotinfo.csv"

# Dimensions exceeding this limit trigger the Sparse/Iterative solver.
# This prevents out-of-memory (OOM) errors on large 12-crossing bidegrees.
DENSE_THRESHOLD = 4000 

def init_db():
    """
    Initializes the SQLite relational database.
    Schema includes a master table for Knot metadata and a detailed table 
    for individual (h, q) bidegree results.
    """
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # KNOTS TABLE: Anchors the knot diagram and its global invariants
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS knots (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            pd_code TEXT,
            khovanov_polynomial TEXT
        )
    ''')
    
    # BIDEGREES TABLE: Stores high-resolution spectral data for each grading
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS bidegrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            knot_id INTEGER,
            h INTEGER,
            q INTEGER,
            dimension INTEGER,
            betti INTEGER,
            smallest_nonzero REAL,
            log_torsion REAL,      
            is_truncated BOOLEAN,
            compute_time REAL,
            FOREIGN KEY (knot_id) REFERENCES knots (id)
        )
    ''')
    conn.commit()
    return conn

def format_kh_poly(betti_dict):
    """
    Converts a mapping of {(h, q): betti} into a standard Khovanov Polynomial string.
    Example: { (0,1): 1, (2,5): 1 } -> "q^1 + t^2q^5"
    """
    terms = []
    # Sort by h and then q for standard mathematical ordering
    for (h, q) in sorted(betti_dict.keys()):
        coeff = betti_dict[(h, q)]
        c_str = f"{coeff}" if coeff > 1 else ""
        t_str = f"t^{h}" if h != 0 else ""
        q_str = f"q^{q}" if q != 0 else ""
        
        term = f"{c_str}{t_str}{q_str}"
        if not term: term = f"{coeff}" # Case for constant term
        terms.append(term)
    
    return " + ".join(terms) if terms else "0"

def parse_pd_code(pd_string):
    """
    Sanitizes and converts Mathematica-style Planar Diagram strings into Python lists.
    Handles semicolon separators and varying bracket spacing.
    """
    if not isinstance(pd_string, str): return None
    # Normalize separators
    cleaned = pd_string.replace(';', ',').replace(' ', '')
    cleaned = re.sub(r'[^0-9,\[\]]', '', cleaned)
    try:
        return ast.literal_eval(cleaned)
    except Exception:
        return None

def process_batch():
    """
    Main processing loop. Iterates through the CSV dataset, performs
    homology calculations, and commits results to the database.
    """
    conn = init_db()
    
    try:
        df = pd.read_csv(CSV_PATH, sep=',')
    except Exception as e:
        print(f"CRITICAL: Failed to read {CSV_PATH}. Error: {e}")
        return

    for _, row in df.iterrows():
        # Positional indexing ensures compatibility with various CSV headers
        knot_name = str(row.iloc[0]).strip()
        pd_raw = row.iloc[1]
        
        if pd.isna(pd_raw): continue
        pd_list = parse_pd_code(str(pd_raw))
        if not pd_list: continue

        # --- STEP 1: RESUME LOGIC ---
        cursor = conn.cursor()
        cursor.execute("INSERT OR IGNORE INTO knots (name, pd_code) VALUES (?, ?)", 
                       (knot_name, json.dumps(pd_list)))
        conn.commit()
        
        knot_id = cursor.execute("SELECT id FROM knots WHERE name=?", (knot_name,)).fetchone()[0]
        
        # If bidegree entries exist for this knot, we assume it's finished or in progress
        if cursor.execute("SELECT count(*) FROM bidegrees WHERE knot_id=?", (knot_id,)).fetchone()[0] > 0:
            print(f"SKIP: {knot_name} is already present in database.")
            continue

        print(f"\n>>> PROCESSING: {knot_name}")
        
        try:
            engine = KhovanovEngine(pd_list)
            active_degrees = engine.get_active_bidegrees()
        except Exception as e:
            print(f"ERROR: Engine initialization failed for {knot_name}: {e}")
            continue

        betti_accumulator = {}

        # --- STEP 2: SPECTRAL ANALYSIS ---
        for h, q in active_degrees:
            start_time = time.time()
            dim = len(engine.chain_groups.get((h, q), []))
            
            if dim == 0:
                continue

            # Switch logic between exact full-spectrum and sparse-approximation
            if dim < DENSE_THRESHOLD:
                # FULL SPECTRUM: Allows for exact torsion (log_torsion) calculation
                spec = engine.get_spectrum(h, q)
                betti = int(np.sum(spec < 1e-9))
                non_zero_eigs = spec[spec > 1e-9]
                
                gap = float(non_zero_eigs[0]) if len(non_zero_eigs) > 0 else None
                log_torsion = float(np.sum(np.log(non_zero_eigs))) if len(non_zero_eigs) > 0 else 0.0
                is_truncated = False
            else:
                # SPARSE SPECTRUM: Efficiently find Betti numbers for large matrices
                betti, gaps = engine.get_low_spectrum_adaptive(h, q, tol=1e-9)
                gap = float(gaps[0]) if len(gaps) > 0 else None
                log_torsion = None # Not calculable without full spectrum
                is_truncated = True
            
            # Save Betti for final polynomial construction
            if betti > 0:
                betti_accumulator[(h, q)] = betti
            
            elapsed = time.time() - start_time
            
            # Log results for this specific bidegree
            cursor.execute('''
                INSERT INTO bidegrees 
                (knot_id, h, q, dimension, betti, smallest_nonzero, log_torsion, is_truncated, compute_time)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (knot_id, h, q, dim, betti, gap, log_torsion, is_truncated, elapsed))
        
        # --- STEP 3: FINALIZATION ---
        poly_str = format_kh_poly(betti_accumulator)
        cursor.execute("UPDATE knots SET khovanov_polynomial = ? WHERE id = ?", (poly_str, knot_id))
        conn.commit()
        
        print(f"SUCCESS: {knot_name} | Poly: {poly_str}")

if __name__ == "__main__":
    process_batch()