"""Graded Euler characteristic (Jones) consistency check for a result database.

For each diagram and quantum degree q, homology and chains have the same
Euler characteristic:

    sum_h (-1)^h betti(h, q)  ==  sum_h (-1)^h dim C^{h, q}.

The right side is combinatorial (it is the coefficient of the unnormalised
Jones polynomial), so this test does not depend on any eigensolver.  It is a
necessary condition only: errors in adjacent homological degrees can cancel.

Standard library only.  The database is opened read-only.
"""

from __future__ import annotations

import argparse
import collections
import json
import sqlite3
import struct


def decode_betti(value):
    if value is None:
        return None
    if isinstance(value, bytes):
        return struct.unpack("<q", value)[0]
    return int(value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out")
    args = ap.parse_args()
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    tables = [r[0] for r in con.execute("select name from sqlite_master where type='table'")]
    table, fk = next((t, f) for t, f in (("knots", "knot_id"), ("links", "link_id"),
                                         ("twisted_unknots", "twisted_unknot_id"))
                     if t in tables)
    names = dict(con.execute(f"select id, name from {table}"))
    chi_h = collections.defaultdict(int)
    chi_c = collections.defaultdict(int)
    null_betti = collections.Counter()
    rows = 0
    for eid, h, q, dim, betti in con.execute(
            f"select {fk}, h, q, dimension, betti from bidegrees"):
        rows += 1
        betti = decode_betti(betti)
        if betti is None:
            null_betti[eid] += 1
            continue
        sign = -1 if h % 2 else 1
        chi_h[(eid, q)] += sign * betti
        chi_c[(eid, q)] += sign * dim
    bad = collections.defaultdict(list)
    for key in chi_c:
        if chi_h[key] != chi_c[key]:
            bad[names[key[0]]].append({"q": key[1], "homology": chi_h[key],
                                       "chains": chi_c[key]})
    report = {
        "database": args.db, "table": table, "diagrams": len(names), "rows": rows,
        "diagram_q_pairs": len(chi_c),
        "diagrams_with_null_betti": len(null_betti),
        "diagrams_failing": len(bad),
        "diagram_q_pairs_failing": sum(len(v) for v in bad.values()),
        "failures": dict(sorted(bad.items())),
    }
    if args.out:
        with open(args.out, "w") as f:
            json.dump(report, f, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k != "failures"}))
    for name in list(report["failures"])[:8]:
        print(" ", name, report["failures"][name][:3])


if __name__ == "__main__":
    main()
