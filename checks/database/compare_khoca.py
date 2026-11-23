"""Compare stored Betti numbers with an independent Khovanov homology program.

For every diagram in a knot database, the rational Khovanov homology of the
stored planar diagram is computed with khoca (L. Lewark and A. Lobb), which
uses Bar-Natan's tangle algorithm and shares no code with `spectral_kh`.  The
result is compared with

  * every stored bidegree row whose Betti number is not null, and
  * the stored Khovanov polynomial, for diagrams that have one.

khoca writes the homological degree with the opposite sign, so its term
t^a q^b is compared with the stored bidegree (h, q) = (-a, b).  The same
convention is applied to every diagram; no per-diagram mirror is allowed.

This checks Betti numbers only.  It says nothing about the spectral gaps.

Needs khoca (and its dependency py).  The database is opened read-only.

    python checks/database/compare_khoca.py --db databases/knot_research_14.db \
        --out checks/database/khoca_knot_research_14.json --workers 8
"""

from __future__ import annotations

import argparse
import ast
import collections
import hashlib
import json
import multiprocessing
import re
import sqlite3
import struct

TERM = re.compile(r"(\d+)?t\^(-?\d+)q\^(-?\d+)")
STORED_TERM = re.compile(r"(\d+)?q\^(-?\d+)t\^(-?\d+)")
_calculator = None


def decode_betti(value):
    if value is None:
        return None
    if isinstance(value, bytes):
        return struct.unpack("<q", value)[0]
    return int(value)


def khoca_betti(pd):
    """Rational Betti numbers {(h, q): rank} in the database's convention."""
    global _calculator
    if _calculator is None:
        import khoca
        _calculator = khoca.InteractiveCalculator(1)  # 1 = rational coefficients
    _, messages = _calculator(pd, print_messages=True)
    lines = messages[messages.index("Unreduced Homology:") + 1:]
    polynomial = next(line for line in lines if "t^" in line)
    betti = collections.Counter()
    for term in polynomial.split(" + "):
        m = TERM.fullmatch(term.strip())
        if not m:
            raise ValueError(f"cannot read khoca term {term!r}")
        betti[-int(m.group(2)), int(m.group(3))] += int(m.group(1) or 1)
    return dict(betti)


def stored_polynomial(text):
    betti = collections.Counter()
    for term in text.split(" + "):
        m = STORED_TERM.fullmatch(term.strip())
        if not m:
            raise ValueError(f"cannot read stored term {term!r}")
        betti[int(m.group(3)), int(m.group(2))] += int(m.group(1) or 1)
    return dict(betti)


def check(job):
    name, pd_text, polynomial, rows = job
    record = {"name": name, "complete": polynomial is not None,
              "rows_compared": 0, "row_mismatches": [], "null_betti_rows": 0}
    try:
        reference = khoca_betti(ast.literal_eval(pd_text))
    except Exception as exc:  # recorded, never silently skipped
        record["error"] = repr(exc)
        return record
    record["total_rank"] = sum(reference.values())
    seen = set()
    for h, q, betti in rows:
        seen.add((h, q))
        if betti is None:
            record["null_betti_rows"] += 1
            continue
        record["rows_compared"] += 1
        if betti != reference.get((h, q), 0):
            record["row_mismatches"].append(
                {"h": h, "q": q, "stored": betti, "khoca": reference.get((h, q), 0)})
    missing = sorted(k for k in reference if k not in seen)
    record["homology_without_row"] = [{"h": h, "q": q, "khoca": reference[h, q]}
                                      for h, q in missing]
    if polynomial is not None:
        try:
            record["polynomial_agrees"] = stored_polynomial(polynomial) == reference
        except ValueError as exc:
            record["polynomial_agrees"] = False
            record["polynomial_error"] = repr(exc)
    return record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--limit", type=int, help="check only the first N diagrams")
    args = ap.parse_args()

    with open(args.db, "rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    rows = collections.defaultdict(list)
    for knot_id, h, q, betti in con.execute("select knot_id, h, q, betti from bidegrees"):
        rows[knot_id].append((h, q, decode_betti(betti)))
    jobs = [(name, pd, poly, rows[knot_id]) for knot_id, name, pd, poly in con.execute(
        "select id, name, pd_code, khovanov_polynomial from knots order by name")]
    if args.limit:
        jobs = jobs[:args.limit]

    with multiprocessing.Pool(args.workers) as pool:
        records = pool.map(check, jobs, chunksize=8)

    def bad(r):
        return ("error" in r or r["row_mismatches"]
                or (r["complete"] and (r["homology_without_row"] or not r["polynomial_agrees"])))

    complete = [r for r in records if r["complete"]]
    partial = [r for r in records if not r["complete"]]
    import khoca
    from importlib.metadata import version
    summary = {
        "database": args.db,
        "sha256": digest,
        "reference": f"khoca {version('khoca')}, rational coefficients",
        "convention": "khoca term t^a q^b compared with stored bidegree (h, q) = (-a, b)",
        "diagrams": len(records),
        "complete_diagrams": len(complete),
        "complete_diagrams_agreeing": sum(not bad(r) for r in complete),
        "partial_diagrams": len(partial),
        "partial_diagrams_with_all_stored_rows_agreeing": sum(not bad(r) for r in partial),
        "rows_compared": sum(r["rows_compared"] for r in records),
        "rows_disagreeing": sum(len(r["row_mismatches"]) for r in records),
        "rows_with_null_betti": sum(r["null_betti_rows"] for r in records),
        "errors": sum("error" in r for r in records),
        "disagreeing_diagrams": [r for r in records if bad(r)],
        "partial_diagram_records": partial,
    }
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=1)
        f.write("\n")
    for key, value in summary.items():
        if not isinstance(value, list):
            print(f"{key}: {value}")
    print("disagreeing diagrams:", [r["name"] for r in summary["disagreeing_diagrams"]][:30])


if __name__ == "__main__":
    main()
