"""Exact rational Betti number of single bidegrees, by integer rank (FLINT).

For integer differentials with ``d_out @ d_in == 0`` (checked exactly),

    beta_Q = dim - rank_Q(d_out) - rank_Q(d_in).

``flint.fmpz_mat.rank`` computes the rank over the rationals exactly, so the
result is a proof for the block, not a modular bound.  Intended for blocks of
dimension up to a few thousand, such as ``13a_4874`` at ``(10, 31)``, where
two modular computations had only established ``beta_Q <= 1``.

    python checks/database/exact_block_betti.py --db databases/knot_research.db \
        --block 13a_4874 10 31 --out checks/database/exact_blocks.json
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import sqlite3
import time

import flint
import numpy as np

from spectral_kh.khovanov_core import KhovanovEngine


def exact_rank(matrix):
    dense = matrix.toarray()
    if dense.size == 0 or not dense.any():
        return 0
    return flint.fmpz_mat(dense.shape[0], dense.shape[1],
                          [int(x) for x in dense.ravel()]).rank()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--block", nargs=3, action="append", required=True,
                    metavar=("NAME", "H", "Q"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--pd-csv", default="data/knotinfo.csv",
                    help="input PD table, used when the database stores no PD code")
    args = ap.parse_args()

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    tables = [r[0] for r in con.execute("select name from sqlite_master where type='table'")]
    table = next(t for t in ("knots", "links", "twisted_unknots") if t in tables)
    results, engines = [], {}
    for name, h, q in args.block:
        h, q = int(h), int(q)
        start = time.perf_counter()
        if name not in engines:
            (pd,) = con.execute(f"select pd_code from {table} where name = ?", (name,)).fetchone()
            if pd is None:
                # Many 13-crossing rows carry no PD code; use the input table.
                with open(args.pd_csv, newline="") as f:
                    pd = next(r["PD Notation"] for r in csv.DictReader(f)
                              if r["Name"] == name).replace(";", ",")
            engines[name] = KhovanovEngine(ast.literal_eval(pd))
        engine = engines[name]
        d_out = engine.get_differential_matrix(h, q).astype(np.int64)
        d_in = engine.get_differential_matrix(h - 1, q).astype(np.int64)
        product = d_out @ d_in
        d2_zero = bool(product.nnz == 0 or not product.data.any())
        dim = len(engine.chain_groups.get((h, q), []))
        r_out, r_in = exact_rank(d_out), exact_rank(d_in)
        results.append({"name": name, "h": h, "q": q, "dimension": dim,
                        "d2_exactly_zero": d2_zero,
                        "rank_out_over_Q": r_out, "rank_in_over_Q": r_in,
                        "betti_over_Q": dim - r_out - r_in if d2_zero else None,
                        "seconds": round(time.perf_counter() - start, 2)})
        print(results[-1], flush=True)
    with open(args.out, "w") as f:
        json.dump({"method": "flint.fmpz_mat.rank (exact over Q)",
                   "flint_version": flint.__version__, "blocks": results}, f, indent=1)


if __name__ == "__main__":
    main()
