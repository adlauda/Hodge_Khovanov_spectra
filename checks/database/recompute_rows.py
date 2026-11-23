"""Recompute stored bidegrees with exact ranks and rank-indexed spectral gaps.

Why this differs from the historical drivers
--------------------------------------------
The historical drivers classified an eigenvalue of the Hodge Laplacian as zero
when it fell below a floating threshold, and the earliest driver asked a
Lanczos solver for the 50 smallest eigenvalues.  Lanczos does not reliably
return every copy of a repeated eigenvalue, so the kernel could be undercounted.

Here nothing is classified by a threshold.

* Betti numbers come from ranks of the integer differentials,
  ``beta = dim - rank(d_out) - rank(d_in)``, after an exact integer check that
  ``d_out @ d_in == 0``.  Ranks are computed modulo the prime 65521 by a blocked elimination
  whose products are exact in float64 (``spectral_kh.exact_rank``); a
  self-test compares it with FLINT.
  A modular rank cannot exceed the rational rank, so the value reported is a
  rigorous *upper* bound on the rational Betti number; it is compared
  elsewhere with KnotInfo's independently computed integral homology.
* The positive spectrum of the Laplacian in degree h is the union of the
  positive spectra of ``d_out^T d_out`` and ``d_in d_in^T`` (Hodge
  decomposition), and each of those equals the squared nonzero singular values
  of one differential.  For a differential of exact rank r, the smallest
  positive eigenvalue of its k-by-k Gram matrix is therefore the eigenvalue of
  index ``k - r`` in ascending order.  The index is fixed by the exact rank.
  The largest eigenvalue below that index is recorded as ``zero_leak``; it
  measures the floating-point error and should be many orders below the gap.
* As a cross-check the same gap is read off the full Laplacian, again at the
  index given by the exact Betti number (dense for moderate dimensions,
  shift-invert Lanczos asking for ``beta + 3`` eigenvalues otherwise).

The database is opened read-only; results go to a JSON-lines file.
"""

from __future__ import annotations

import argparse
import ast
import glob
import json
import os
import sqlite3
import struct
import sys
import time

import numpy as np
import scipy.linalg
from scipy import sparse
from scipy.sparse.linalg import eigsh

from spectral_kh.khovanov_core import KhovanovEngine
# The rank and Gram-matrix routines were written in this script, which is the
# form that produced the recorded results; they now live in the library
# unchanged so the database drivers use the same code.
from spectral_kh.exact_rank import gram_floor, rank_mod_p, self_test


def decode_betti(value):
    if value is None:
        return None
    if isinstance(value, bytes):
        return struct.unpack("<q", value)[0]
    return int(value)


def entity_table(con):
    tables = [r[0] for r in con.execute("select name from sqlite_master where type='table'")]
    for table, fk in (("knots", "knot_id"), ("links", "link_id"),
                      ("twisted_unknots", "twisted_unknot_id")):
        if table in tables:
            return table, fk
    raise SystemExit("no recognised entity table")


def laplacian_gap(engine, h, q, betti, dense_max):
    """Gap of the full Laplacian at the index given by the exact Betti number."""
    lap = engine.compute_laplacian(h, q)
    dim = lap.shape[0]
    if betti >= dim:
        return None, "no positive spectrum"
    if dim <= dense_max or betti + 3 >= dim - 1:
        lo = max(betti - 1, 0)
        w = scipy.linalg.eigvalsh(lap.toarray().astype(float), subset_by_index=[lo, betti])
        return float(w[-1]), "dense"
    w = np.sort(eigsh(lap.astype(float), k=betti + 3, sigma=-1e-5,
                      return_eigenvectors=False))
    return float(w[betti]), "shift-invert"


def large_gap(engine, h, q, betti, seeds=(1, 2)):
    """Gap for blocks too large for dense work, with the Betti number given.

    Matrix-free Lanczos for the smallest eigenvalues (no factorisation), run
    from two independent random start vectors and asking for ``betti + 4``
    eigenvalues.  A converged Lanczos run cannot skip the smallest eigenvalue
    of a symmetric matrix unless the start vector is orthogonal to its
    eigenspace, so the smallest positive value is reliable even though
    multiplicities are not; that is why the Betti number is supplied rather
    than counted.  The caller compares the two runs, and ``zero_leak`` reports
    the largest returned eigenvalue at an index below ``betti``.
    """
    lap = engine.compute_laplacian(h, q).astype(float)
    dim = lap.shape[0]
    if betti >= dim:
        return None, None, None
    runs = []
    for seed in seeds:
        v0 = np.random.default_rng(seed).standard_normal(dim)
        k = min(betti + 4, dim - 1)
        runs.append(np.sort(eigsh(lap, k=k, which="SA", v0=v0, tol=1e-11,
                                  ncv=min(dim, max(4 * k, 80)), maxiter=200000,
                                  return_eigenvectors=False)))
    gaps = [float(w[betti]) for w in runs]
    leak = max(float(abs(w[betti - 1])) if betti else 0.0 for w in runs)
    return gaps[0], gaps[1], leak


def label_shift(engine, stored_dims):
    """Offset (dh, dq) taking the database's bidegree labels to the engine's.

    Some historical tables were written with a different grading normalisation
    from the current engine.  The offset is accepted only if it carries the
    whole stored table of chain dimensions onto the engine's table exactly.
    """
    current = {k: len(v) for k, v in engine.chain_groups.items() if len(v)}
    stored = {k: d for k, d in stored_dims.items() if d}
    if not stored or not current:
        raise ValueError("empty dimension table")
    dh = min(h for h, _ in current) - min(h for h, _ in stored)
    dq = min(q for _, q in current) - min(q for _, q in stored)
    if {(h + dh, q + dq): d for (h, q), d in stored.items()} != current:
        raise ValueError("no constant label offset matches the stored dimensions")
    return dh, dq


def process(name, pd_code, rows, dense_max, rank_max, stored_dims=None):
    t0 = time.perf_counter()
    engine = KhovanovEngine(pd_code)
    build = time.perf_counter() - t0
    dh, dq = label_shift(engine, stored_dims) if stored_dims is not None else (0, 0)
    mats, ranks, floors = {}, {}, {}

    def mat(h, q):
        if (h, q) not in mats:
            mats[(h, q)] = engine.get_differential_matrix(h, q).astype(np.int64)
        return mats[(h, q)]

    def rank(h, q):
        if (h, q) not in ranks:
            ranks[(h, q)] = rank_mod_p(mat(h, q))
        return ranks[(h, q)]

    def floor(h, q):
        if (h, q) not in floors:
            floors[(h, q)] = gram_floor(mat(h, q), rank(h, q))
        return floors[(h, q)]

    out = []
    for row in rows:
        t1 = time.perf_counter()
        h, q = row["h"] + dh, row["q"] + dq  # engine labels
        dim = len(engine.chain_groups.get((h, q), []))
        d_out, d_in = mat(h, q), mat(h - 1, q)
        if max(d_out.shape + d_in.shape) > rank_max:
            # Too large for dense exact ranks: take the stored Betti number
            # (to be justified separately) and read the gap at that index.
            product = d_out @ d_in
            betti = row["betti"]
            gap, cross, leak = large_gap(engine, h, q, betti)
            out.append({
                "name": name, "h": row["h"], "q": row["q"], "row_id": row["id"],
                "engine_h": h, "engine_q": q,
                "dimension": dim, "stored_dimension": row["dimension"],
                "stored_betti": betti, "stored_gap": row["gap"],
                "stored_is_truncated": row["is_truncated"],
                "d2_exactly_zero": bool(product.nnz == 0 or not product.data.any()),
                "rank_out": None, "rank_in": None, "betti_mod_p": None,
                "gap": gap, "zero_leak": leak,
                "gap_from_out": None, "gap_from_in": None,
                "laplacian_gap": cross, "laplacian_method": "matrix-free Lanczos, second start vector",
                "seconds": round(time.perf_counter() - t1, 3),
            })
            mats.clear()
            continue
        product = d_out @ d_in
        d2_zero = bool(product.nnz == 0 or not product.data.any())
        r_out, r_in = rank(h, q), rank(h - 1, q)
        betti = dim - r_out - r_in
        f_out, f_in = floor(h, q), floor(h - 1, q)
        candidates = [f for f in (f_out, f_in) if f is not None]
        gap = min(f["floor"] for f in candidates) if candidates else None
        leak = max(f["zero_leak"] for f in candidates) if candidates else None
        cross, method = laplacian_gap(engine, h, q, betti, dense_max)
        out.append({
            "name": name, "h": row["h"], "q": row["q"], "row_id": row["id"],
            "engine_h": h, "engine_q": q,
            "dimension": dim, "stored_dimension": row["dimension"],
            "stored_betti": row["betti"], "stored_gap": row["gap"],
            "stored_is_truncated": row["is_truncated"],
            "d2_exactly_zero": d2_zero,
            "rank_out": r_out, "rank_in": r_in, "betti_mod_p": betti,
            "gap": gap, "zero_leak": leak,
            "gap_from_out": f_out["floor"] if f_out else None,
            "gap_from_in": f_in["floor"] if f_in else None,
            "laplacian_gap": cross, "laplacian_method": method,
            "seconds": round(time.perf_counter() - t1, 3),
        })
    return out, build


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db")
    ap.add_argument("--out")
    ap.add_argument("--select", choices=["truncated", "all"], default="truncated")
    ap.add_argument("--names", nargs="*", help="restrict to these diagram names")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--dense-max", type=int, default=3000,
                    help="largest Laplacian cross-checked with a dense solver")
    ap.add_argument("--rank-max", type=int, default=25000,
                    help="largest differential handled with dense exact ranks")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--auto-shift", action="store_true",
                    help="match the database's bidegree labels to the engine's by "
                         "the unique constant offset that reproduces all stored dimensions")
    ap.add_argument("--reverse", action="store_true",
                    help="work through the shard from the end of the name list")
    ap.add_argument("--done-glob",
                    help="other workers' output files; diagrams recorded there are skipped")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    table, fk = entity_table(con)
    where = "where b.is_truncated = 1" if args.select == "truncated" else ""
    work = {}
    for name, pd, rid, h, q, dim, betti, gap, trunc in con.execute(
        f"select e.name, e.pd_code, b.id, b.h, b.q, b.dimension, b.betti, "
        f"b.smallest_nonzero, b.is_truncated from bidegrees b join {table} e "
        f"on e.id = b.{fk} {where} order by e.name, b.h, b.q"
    ):
        if args.names and name not in args.names:
            continue
        work.setdefault(name, (pd, []))[1].append(
            {"id": rid, "h": h, "q": q, "dimension": dim, "betti": decode_betti(betti),
             "gap": gap, "is_truncated": trunc})
    names = sorted(work)[args.shard::args.shards]
    if args.limit:
        names = names[:args.limit]

    if args.reverse:
        names = names[::-1]

    def recorded():
        """Names already written by this worker or by any worker in --done-glob."""
        paths = set(glob.glob(args.done_glob)) if args.done_glob else set()
        paths.add(args.out)
        found = set()
        for path in paths:
            if os.path.exists(path):
                with open(path) as f:
                    for line in f:
                        try:
                            found.add(json.loads(line)["name"])
                        except (ValueError, KeyError):
                            pass  # a line another worker is still writing
        return found

    print(f"{len(names)} diagrams in shard, {len(recorded() & set(names))} already done",
          flush=True)

    for name in names:
        # Re-read before every diagram so workers sharing a list do not repeat
        # each other's finished work.
        if name in recorded():
            continue
        pd_code, rows = work[name]
        start = time.perf_counter()
        stored_dims = None
        if args.auto_shift:
            stored_dims = {(h, q): d for h, q, d in con.execute(
                f"select b.h, b.q, b.dimension from bidegrees b join {table} e "
                f"on e.id = b.{fk} where e.name = ?", (name,))}
        results, build = process(name, ast.literal_eval(pd_code), rows, args.dense_max,
                                 args.rank_max, stored_dims)
        # One diagram is written in a single append so a resumed run never
        # sees a partially recorded diagram as complete.
        with open(args.out, "a") as f:
            f.write("".join(json.dumps(r) + "\n" for r in results))
        print(f"{name}: {len(rows)} rows, build {build:.1f}s, "
              f"total {time.perf_counter() - start:.1f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
