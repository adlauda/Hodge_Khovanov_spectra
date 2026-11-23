"""Compare stored link Betti numbers with LinkInfo's Khovanov polynomials.

LinkInfo lists one rational Khovanov Poincare polynomial per orientation
(``L6a4{0}``, ``L6a4{1}``, ...), while the link database stores one diagram per
unoriented name.  A stored table is accepted when it equals the LinkInfo table
of some orientation, or its mirror ``(h, q) -> (-h, -q)``.

Standard library only.  The database is opened read-only.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import sqlite3
import struct

FACTOR = re.compile(r"^([qt])(?:\^\(?(-?\d+)\)?)?$")


def decode_betti(value):
    if value is None:
        return None
    if isinstance(value, bytes):
        return struct.unpack("<q", value)[0]
    return int(value)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_side(text):
    """Parse a product such as ``2*q^4*t^2`` into (coefficient, q-exp, t-exp)."""
    text = text.strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    coeff, exps = 1, {"q": 0, "t": 0}
    for factor in text.split("*"):
        factor = factor.strip()
        if re.fullmatch(r"\d+", factor):
            coeff *= int(factor)
            continue
        m = FACTOR.match(factor)
        if not m:
            raise ValueError(f"cannot parse factor {factor!r}")
        exps[m.group(1)] += int(m.group(2) or 1)
    return coeff, exps["q"], exps["t"]


def parse_polynomial(text):
    table = collections.Counter()
    for term in text.split(" + "):
        num, _, den = term.partition("/")
        c, q, t = parse_side(num)
        if den:
            dc, dq, dt = parse_side(den)
            if dc != 1:
                raise ValueError(f"unexpected denominator in {term!r}")
            q, t = q - dq, t - dt
        table[(t, q)] += c
    return dict(table)


def mirror(table):
    return {(-h, -q): b for (h, q), b in table.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--linkinfo", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    csv.field_size_limit(10**9)
    ref = collections.defaultdict(dict)
    with open(args.linkinfo, newline="") as f:
        first = f.readline()
        delim = "|" if first.count("|") > first.count(",") else ","
        f.seek(0)
        for row in csv.DictReader(f, delimiter=delim):
            name, poly = row["name"], row.get("khovanov_polynomial", "").strip()
            if "{" not in name or not poly:
                continue
            base, orient = name.split("{")
            ref[base]["{" + orient] = parse_polynomial(poly)

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    names = dict(con.execute("select id, name from links"))
    stored = collections.defaultdict(dict)
    truncated = collections.defaultdict(set)
    for lid, h, q, betti, trunc in con.execute(
            "select link_id, h, q, betti, is_truncated from bidegrees order by id"):
        stored[lid][(h, q)] = decode_betti(betti)
        if trunc == 1:
            truncated[lid].add((h, q))

    summary = collections.Counter()
    mismatches, matched_as = [], {}
    for lid, name in sorted(names.items()):
        table = stored.get(lid, {})
        nonzero = {k: b for k, b in table.items() if b}
        if name not in ref:
            summary["links_without_reference"] += 1
            continue
        best = None
        for orient, poly in sorted(ref[name].items()):
            for tag, cand in (("direct", poly), ("mirror", mirror(poly))):
                keys = set(cand) | set(table)
                bad = [k for k in keys if cand.get(k, 0) != table.get(k, 0)]
                if best is None or len(bad) < len(best[0]):
                    best = (bad, orient, tag, cand)
        bad, orient, tag, cand = best
        matched_as[name] = f"{orient} {tag}"
        summary["rows_compared"] += len(table)
        summary["truncated_rows_compared"] += len(truncated[lid])
        if bad:
            summary["links_disagreeing"] += 1
            summary["rows_disagreeing"] += len(bad)
            summary["truncated_rows_disagreeing"] += sum(k in truncated[lid] for k in bad)
            mismatches.append({
                "name": name, "closest": matched_as[name],
                "bidegrees": [{"h": k[0], "q": k[1], "stored": table.get(k),
                               "linkinfo": cand.get(k, 0),
                               "truncated": k in truncated[lid]} for k in sorted(bad)]})
        else:
            summary["links_agreeing"] += 1

    report = {"database": args.db, "database_sha256": sha256(args.db),
              "linkinfo_sha256": sha256(args.linkinfo),
              "links_in_database": len(names),
              "summary": dict(sorted(summary.items())),
              "matched_orientation": matched_as, "mismatches": mismatches}
    with open(args.out, "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("matched_orientation", "mismatches")}, indent=1))
    for m in mismatches[:12]:
        print(m["name"], m["closest"], len(m["bidegrees"]), m["bidegrees"][:3])


if __name__ == "__main__":
    main()
