"""Compare every stored Betti number with KnotInfo's integral Khovanov homology.

KnotInfo's ``khovanov_unreduced_integral_vector`` lists entries
``[torsion, multiplicity, h, q]``; ``torsion == 0`` is the free part, so the
rational Betti number in bidegree (h, q) is the sum of those multiplicities.
That table was produced by programs unrelated to this repository's numerical
engine, so agreement is an independent check of the stored ``betti`` column.

Standard library only.  The database is opened read-only.

    python checks/database/compare_knotinfo.py \
        --db databases/knot_research.db \
        --knotinfo /path/to/knotinfo_data_complete.csv \
        --cleanup Paper/Current/reproducibility/DATABASE-CLEANUP.json \
        --out checks/database/knotinfo_comparison.json
"""

from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import sqlite3
import struct


def decode_betti(value):
    """Historical rows store Betti numbers as 8-byte little-endian blobs."""
    if value is None:
        return None
    if isinstance(value, bytes):
        if len(value) != 8:
            raise ValueError(f"unexpected Betti blob length {len(value)}")
        return struct.unpack("<q", value)[0]
    return int(value)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_knotinfo(path):
    csv.field_size_limit(10**9)
    with open(path, newline="") as f:
        first = f.readline()
        delim = "|" if first.count("|") > first.count(",") else ","
        f.seek(0)
        ref = {}
        for row in csv.DictReader(f, delimiter=delim):
            vec = row.get("khovanov_unreduced_integral_vector", "").strip()
            if not vec.startswith("["):
                continue
            betti = collections.Counter()
            for torsion, mult, h, q in ast.literal_eval(vec):
                if torsion == 0:
                    betti[(h, q)] += mult
            ref[row["name"]] = dict(betti)
    return ref


def mirror(table):
    return {(-h, -q): b for (h, q), b in table.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--knotinfo", required=True)
    ap.add_argument("--cleanup")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ref = load_knotinfo(args.knotinfo)
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    names = dict(con.execute("select id, name from knots"))
    stored = collections.defaultdict(dict)
    truncated = collections.defaultdict(set)
    duplicates = 0
    for kid, h, q, betti, trunc in con.execute(
        "select knot_id, h, q, betti, is_truncated from bidegrees order by id"
    ):
        if (h, q) in stored[kid]:
            duplicates += 1
        stored[kid][(h, q)] = decode_betti(betti)
        if trunc == 1:
            truncated[kid].add((h, q))

    summary = collections.Counter()
    mismatches = []
    orientation = {}
    for kid, name in sorted(names.items()):
        table = stored.get(kid, {})
        if name not in ref:
            summary["knots_without_reference"] += 1
            continue
        nonzero = {k: b for k, b in table.items() if b}
        direct, mirrored = ref[name], mirror(ref[name])
        # Choose the chirality that agrees in more bidegrees; an amphichiral
        # table makes the two identical.
        score = lambda r: sum(1 for k in set(r) | set(nonzero) if r.get(k, 0) == nonzero.get(k, 0))
        use, tag = (direct, "direct") if score(direct) >= score(mirrored) else (mirrored, "mirror")
        orientation[name] = tag
        summary[f"orientation_{tag}"] += 1
        bad = []
        for key in sorted(set(use) | set(table)):
            want = use.get(key, 0)
            if key not in table:
                bad.append({"h": key[0], "q": key[1], "stored": None, "knotinfo": want,
                            "truncated": False})
                summary["reference_bidegrees_missing_from_db"] += 1
                continue
            summary["rows_compared"] += 1
            if key in truncated[kid]:
                summary["truncated_rows_compared"] += 1
            if table[key] != want:
                is_trunc = key in truncated[kid]
                bad.append({"h": key[0], "q": key[1], "stored": table[key],
                            "knotinfo": want, "truncated": is_trunc})
                summary["rows_disagreeing"] += 1
                summary["truncated_rows_disagreeing" if is_trunc
                        else "untruncated_rows_disagreeing"] += 1
        if bad:
            summary["knots_disagreeing"] += 1
            mismatches.append({"name": name, "orientation": tag, "bidegrees": bad})
        else:
            summary["knots_agreeing"] += 1

    report = {
        "database": args.db,
        "database_sha256": sha256(args.db),
        "knotinfo_sha256": sha256(args.knotinfo),
        "knots_in_database": len(names),
        "duplicate_bidegree_rows_seen": duplicates,
        "summary": dict(sorted(summary.items())),
        "mismatches": mismatches,
    }

    if args.cleanup:
        with open(args.cleanup) as f:
            cleanup = json.load(f)
        verdict = collections.Counter()
        detail = []
        for group in cleanup["conflicting_betti_groups"]:
            name, key = group["name"], (group["h"], group["q"])
            table = ref[name] if orientation.get(name) != "mirror" else mirror(ref[name])
            want = table.get(key, 0)
            kept = [r for r in group["rows"] if r["id"] == group["retained_row_id"]][0]
            dropped = [r for r in group["rows"] if r["id"] != group["retained_row_id"]]
            kept_ok = kept["betti"] == want
            dropped_ok = any(r["betti"] == want for r in dropped)
            verdict["retained_matches_knotinfo" if kept_ok else "retained_differs"] += 1
            verdict["removed_matches_knotinfo" if dropped_ok else "removed_differs"] += 1
            detail.append({"name": name, "h": key[0], "q": key[1], "knotinfo": want,
                           "retained": kept["betti"],
                           "removed": [r["betti"] for r in dropped],
                           "removed_is_truncated": [r["is_truncated"] for r in dropped]})
        report["cleanup_conflicts"] = {"groups": len(detail), "verdict": dict(verdict),
                                       "detail": detail}

    with open(args.out, "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("mismatches", "cleanup_conflicts")}, indent=1))
    if args.cleanup:
        print("cleanup conflicts:", report["cleanup_conflicts"]["verdict"])
    for m in mismatches[:15]:
        print(m["name"], m["orientation"], m["bidegrees"][:4])


if __name__ == "__main__":
    main()
