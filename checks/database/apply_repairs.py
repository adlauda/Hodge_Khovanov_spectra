"""Build a repaired copy of a result database from verified recomputations.

Non-destructive: the source database is opened read-only and copied; the script
refuses to overwrite an existing output or report.

What changes in the copy
------------------------
1. Every Betti number stored as an 8-byte little-endian blob is rewritten as
   the same value with integer storage.
2. For each recomputed row that passes the acceptance tests below, ``betti``
   and ``smallest_nonzero`` are replaced by the verified values and
   ``is_truncated`` is cleared.  A row that fails any test is left untouched
   and listed in the report.
3. The stored Khovanov polynomial of every diagram whose Betti numbers changed
   is rebuilt in that diagram's existing string style.
4. With ``--pd-csv``, diagrams that store no PD code receive the code from
   the input table, after checking that every PD code already stored equals
   the table's entry for that diagram.
5. Tables ``row_repairs``, ``polynomial_repairs`` and ``pd_code_fills`` record
   every old and new value.

Acceptance tests for a recomputed row
-------------------------------------
* the recomputed chain dimension equals the stored dimension;
* ``d_out @ d_in`` is exactly zero in integer arithmetic;
* the Betti number equals the reference value (KnotInfo's integral Khovanov
  homology for knots; for links and twisted unknots the stored value, which
  the comparison scripts check against LinkInfo and against the unknot), and
  when exact modular ranks were computed they give the same number;
* the gap read from the differentials' Gram matrices and the gap read from the
  full Laplacian agree to a relative 1e-6;
* the largest eigenvalue classified as zero is below 1e-6 of the gap.
"""

from __future__ import annotations

import argparse
import ast
import collections
import csv
import glob
import hashlib
import json
import os
import re
import shutil
import sqlite3
import struct
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_knotinfo import load_knotinfo  # noqa: E402

TERM = re.compile(r"^(\d*)((?:[qt]\^-?\d+)*)$")
FACTOR = re.compile(r"([qt])\^(-?\d+)")
ENTITIES = (("knots", "knot_id"), ("links", "link_id"),
            ("twisted_unknots", "twisted_unknot_id"))


def decode_betti(value):
    if isinstance(value, bytes):
        if len(value) != 8:
            raise ValueError("unexpected Betti blob length")
        return struct.unpack("<q", value)[0]
    return value


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_polynomial(text):
    table = {}
    if text in (None, "", "0"):
        return table
    for term in text.split(" + "):
        m = TERM.match(term)
        if not m:
            raise ValueError(f"cannot parse term {term!r}")
        exps = {"q": 0, "t": 0}
        for var, exp in FACTOR.findall(m.group(2)):
            exps[var] += int(exp)
        key = (exps["t"], exps["q"])
        if key in table:
            raise ValueError(f"repeated bidegree in {text!r}")
        table[key] = int(m.group(1) or 1)
    return table


def format_qt(table):
    """Style of the parallel drivers: ``2q^-5t^-2``."""
    terms = []
    for (h, q), b in sorted(table.items()):
        terms.append(f"{b if b > 1 else ''}q^{q}t^{h}")
    return " + ".join(terms) if terms else "0"


def format_tq(table):
    """Style of the original driver: ``2t^-8q^-21``, omitting ``t^0``."""
    terms = []
    for (h, q), b in sorted(table.items()):
        terms.append(f"{b if b > 1 else ''}{f't^{h}' if h else ''}q^{q}")
    return " + ".join(terms)


def style_of(text):
    table = parse_polynomial(text)
    for formatter in (format_qt, format_tq):
        if formatter(table) == text:
            return formatter
    return None


def close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--recomputed", nargs="+", required=True,
                    help="JSON-lines files (globs allowed) from recompute_rows.py")
    ap.add_argument("--knotinfo", help="KnotInfo CSV; required for the knot database")
    ap.add_argument("--pd-csv", help="input table with columns Name and 'PD Notation'")
    args = ap.parse_args()

    for path in (args.output, args.report):
        if os.path.exists(path):
            raise SystemExit(f"refusing to overwrite {path}")
    source_digest = sha256(args.source)

    records = {}
    for pattern in args.recomputed:
        for path in sorted(glob.glob(pattern)):
            with open(path) as f:
                for line in f:
                    if line.strip():
                        r = json.loads(line)
                        records[r["row_id"]] = r

    shutil.copyfile(args.source, args.output)
    con = sqlite3.connect(args.output)
    con.execute("pragma foreign_keys = on")
    tables = [r[0] for r in con.execute("select name from sqlite_master where type='table'")]
    table, fk = next((t, f) for t, f in ENTITIES if t in tables)
    reference = None
    if table == "knots":
        if not args.knotinfo:
            raise SystemExit("--knotinfo is required for the knot database")
        reference = load_knotinfo(args.knotinfo)

    before = {rid: (eid, h, q, dim, betti, gap, trunc) for rid, eid, h, q, dim, betti, gap, trunc
              in con.execute(f"select id, {fk}, h, q, dimension, betti, smallest_nonzero, "
                             f"is_truncated from bidegrees")}
    names = dict(con.execute(f"select id, name from {table}"))
    flagged = {rid for rid, row in before.items() if row[6] == 1}

    con.execute("begin immediate")
    con.execute("""create table row_repairs (
        row_id integer primary key, old_betti integer, new_betti integer,
        old_smallest_nonzero real, new_smallest_nonzero real,
        old_is_truncated integer, betti_source text, gap_method text,
        zero_leak real, laplacian_cross_check real)""")
    con.execute("create table polynomial_repairs (name text primary key, "
                "old_polynomial text, new_polynomial text)")

    blobs = [(decode_betti(row[4]), rid) for rid, row in before.items()
             if isinstance(row[4], bytes)]
    con.executemany("update bidegrees set betti = ? where id = ?", blobs)

    rejected, accepted, touched = [], 0, set()
    counts = collections.Counter()
    for rid in sorted(flagged):
        eid, h, q, dim, betti, gap, trunc = before[rid]
        betti = decode_betti(betti)
        name = names[eid]
        r = records.get(rid)
        why = None
        if r is None:
            why = "not recomputed"
        elif (r["name"], r["h"], r["q"]) != (name, h, q):
            why = "record does not match row"
        elif r["dimension"] != dim:
            why = "dimension differs"
        elif not r["d2_exactly_zero"]:
            why = "d^2 != 0"
        if why is None:
            want = reference[name].get((h, q), 0) if reference is not None else betti
            if r["betti_mod_p"] is not None and r["betti_mod_p"] != want:
                why = f"modular Betti {r['betti_mod_p']} != reference {want}"
            elif reference is not None and r["betti_mod_p"] is None:
                why = "no exact rank"
        if why is None:
            new_gap, cross, leak = r["gap"], r["laplacian_gap"], r["zero_leak"]
            if new_gap is None:
                if cross is not None or want != dim:
                    why = "no gap but positive spectrum expected"
            elif cross is None or not close(new_gap, cross, 1e-6):
                why = "gap cross-check disagrees"
            elif leak is None or leak > 1e-6 * new_gap:
                why = "zero eigenvalues not separated from the gap"
        if why is not None:
            rejected.append({"row_id": rid, "name": name, "h": h, "q": q, "reason": why})
            continue
        source = ("reference table and exact modular ranks" if r["betti_mod_p"] is not None
                  else "reference table (block too large for dense ranks)")
        method = ("Gram matrices of the differentials, index from exact rank"
                  if r["rank_out"] is not None
                  else "Laplacian, stored Betti index, two start vectors: " + r["laplacian_method"])
        con.execute("update bidegrees set betti = ?, smallest_nonzero = ?, is_truncated = 0 "
                    "where id = ?", (want, new_gap, rid))
        con.execute("insert into row_repairs values (?,?,?,?,?,?,?,?,?,?)",
                    (rid, betti, want, gap, new_gap, trunc, source, method, leak, cross))
        accepted += 1
        counts["betti_changed"] += betti != want
        if betti != want:
            touched.add(eid)
        if (gap is None) != (new_gap is None):
            counts["gap_null_status_changed"] += 1
        elif gap is not None:
            counts["gap_changed_beyond_1e-6_relative"] += not close(gap, new_gap, 1e-6)
            counts["gap_changed_beyond_1e-9_relative"] += not close(gap, new_gap, 1e-9)

    # Rebuild polynomials of diagrams whose Betti numbers changed, and check
    # every other stored polynomial against the stored Betti numbers.
    after = collections.defaultdict(dict)
    for eid, h, q, betti in con.execute(f"select {fk}, h, q, betti from bidegrees"):
        if betti:
            after[eid][(h, q)] = betti
    inconsistent = []
    for eid, name, poly in con.execute(
            f"select id, name, khovanov_polynomial from {table}").fetchall():
        if poly is None:
            continue
        if eid in touched:
            formatter = style_of(poly)
            if formatter is None:
                raise SystemExit(f"unrecognised polynomial style for {name}")
            new = formatter(after[eid])
            con.execute(f"update {table} set khovanov_polynomial = ? where id = ?", (new, eid))
            con.execute("insert into polynomial_repairs values (?,?,?)", (name, poly, new))
        elif parse_polynomial(poly) != after[eid]:
            inconsistent.append(name)
    pd_filled = pd_mismatch = 0
    if args.pd_csv:
        with open(args.pd_csv, newline="") as f:
            table_pd = {r["Name"]: ast.literal_eval(r["PD Notation"].replace(";", ","))
                        for r in csv.DictReader(f)}
        con.execute("create table pd_code_fills (name text primary key, pd_code text)")
        current = con.execute(f"select id, name, pd_code from {table}").fetchall()
        for eid, name, pd in current:
            if pd is not None and (name not in table_pd or
                                   [list(x) for x in ast.literal_eval(pd)] != table_pd[name]):
                pd_mismatch += 1
        if pd_mismatch:
            raise SystemExit(f"{pd_mismatch} stored PD codes differ from {args.pd_csv}")
        for eid, name, pd in current:
            if pd is None and name in table_pd:
                text = str(table_pd[name])
                con.execute(f"update {table} set pd_code = ? where id = ?", (text, eid))
                con.execute("insert into pd_code_fills values (?, ?)", (name, text))
                pd_filled += 1
    con.commit()

    # Post-conditions.
    changed_other = 0
    repaired_ids = {r[0] for r in con.execute("select row_id from row_repairs")}
    for rid, eid, h, q, dim, betti, gap, trunc in con.execute(
            f"select id, {fk}, h, q, dimension, betti, smallest_nonzero, is_truncated "
            f"from bidegrees"):
        old = before[rid]
        same = (eid, h, q, dim, betti, gap, trunc) == old[:4] + (decode_betti(old[4]),) + old[5:]
        if not same and rid not in repaired_ids:
            changed_other += 1
    storage = dict(con.execute("select typeof(betti), count(*) from bidegrees group by 1"))
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": args.source, "source_sha256": source_digest,
        "entity_table": table,
        "bidegree_rows": len(before),
        "blob_betti_values_rewritten_as_integers": len(blobs),
        "betti_storage_types_after": storage,
        "rows_flagged_truncated_in_source": len(flagged),
        "rows_repaired": accepted,
        "rows_left_untouched": rejected,
        "rows_still_flagged_truncated": con.execute(
            "select count(*) from bidegrees where is_truncated = 1").fetchone()[0],
        "repaired_rows_with_betti_changed": counts["betti_changed"],
        "repaired_rows_with_gap_changed_beyond_1e-9_relative":
            counts["gap_changed_beyond_1e-9_relative"],
        "repaired_rows_with_gap_changed_beyond_1e-6_relative":
            counts["gap_changed_beyond_1e-6_relative"],
        "repaired_rows_with_gap_null_status_changed": counts["gap_null_status_changed"],
        "polynomials_rebuilt": len(touched),
        "pd_codes_filled_from_input_table": pd_filled,
        "diagrams_without_pd_code_after": con.execute(
            f"select count(*) from {table} where pd_code is null").fetchone()[0],
        "polynomials_inconsistent_with_stored_betti": inconsistent,
        "rows_changed_outside_repair_log": changed_other,
        "integrity_check": con.execute("pragma integrity_check").fetchone()[0],
        "foreign_key_violations": len(con.execute("pragma foreign_key_check").fetchall()),
    }
    con.execute("vacuum")
    con.close()
    report["output"] = args.output
    report["output_sha256"] = sha256(args.output)
    with open(args.report, "w") as f:
        json.dump(report, f, indent=1)
    shown = dict(report)
    shown["rows_left_untouched"] = len(rejected)
    shown["polynomials_inconsistent_with_stored_betti"] = len(inconsistent)
    print(json.dumps(shown, indent=1))
    for r in rejected[:10]:
        print("  left untouched:", r)


if __name__ == "__main__":
    main()
