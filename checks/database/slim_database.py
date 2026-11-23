"""Write a smaller copy of a result database without selected bidegree columns.

The published databases omit ``log_torsion`` and ``compute_time``, which no
figure, table or check in this repository reads.  Every other column, table
and index is kept.

Non-destructive: the source is opened read-only and copied; the script refuses
to overwrite an existing output or report.  After dropping the columns it
checks that every remaining column of every table is identical in source and
output, and records both checksums.

    python checks/database/slim_database.py --source full.db \
        --output slim.db --report slim.SLIM.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
from datetime import datetime, timezone

DEFAULT_DROP = ("log_torsion", "compute_time")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def columns(con, schema, table):
    return [r[1] for r in con.execute(f"pragma {schema}.table_info({table})")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--drop", nargs="+", default=list(DEFAULT_DROP),
                    help="columns of the bidegrees table to omit")
    args = ap.parse_args()
    for path in (args.output, args.report):
        if os.path.exists(path):
            raise SystemExit(f"refusing to overwrite {path}")

    source_digest = sha256(args.source)
    shutil.copyfile(args.source, args.output)
    con = sqlite3.connect(f"file:{os.path.abspath(args.output)}", uri=True, isolation_level=None)
    dropped = {}
    for column in args.drop:
        if column in columns(con, "main", "bidegrees"):
            dropped[column] = con.execute(
                f"select count(*) from bidegrees where {column} is not null").fetchone()[0]
            con.execute(f"alter table bidegrees drop column {column}")
    con.execute("vacuum")

    # Every remaining column of every table must be unchanged.
    con.execute("attach ? as src", (f"file:{os.path.abspath(args.source)}?mode=ro",))
    tables = [r[0] for r in con.execute(
        "select name from main.sqlite_master where type='table' and name not like 'sqlite_%'")]
    source_tables = [r[0] for r in con.execute(
        "select name from src.sqlite_master where type='table' and name not like 'sqlite_%'")]
    if sorted(tables) != sorted(source_tables):
        raise SystemExit("table list changed")
    rows = {}
    for table in tables:
        kept = columns(con, "main", table)
        expected = [c for c in columns(con, "src", table)
                    if not (table == "bidegrees" and c in dropped)]
        if kept != expected:
            raise SystemExit(f"unexpected columns in {table}")
        cols = ", ".join(f'"{c}"' for c in kept)
        for a, b in (("main", "src"), ("src", "main")):
            extra = con.execute(f"select count(*) from (select {cols} from {a}.{table} "
                                f"except select {cols} from {b}.{table})").fetchone()[0]
            if extra:
                raise SystemExit(f"{table}: {extra} rows of {a} not in {b}")
        n_out = con.execute(f"select count(*) from main.{table}").fetchone()[0]
        if n_out != con.execute(f"select count(*) from src.{table}").fetchone()[0]:
            raise SystemExit(f"{table}: row count changed")
        rows[table] = n_out
    indexes = [r[0] for r in con.execute(
        "select name from main.sqlite_master where type='index' order by name")]
    if indexes != [r[0] for r in con.execute(
            "select name from src.sqlite_master where type='index' order by name")]:
        raise SystemExit("index list changed")
    integrity = con.execute("pragma main.integrity_check").fetchone()[0]
    con.execute("detach src")
    con.close()
    if integrity != "ok":
        raise SystemExit(f"integrity check: {integrity}")

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": source_digest,
        "source_bytes": os.path.getsize(args.source),
        "columns_dropped_from_bidegrees": list(dropped),
        "non_null_values_dropped": dropped,
        "rows_per_table": rows,
        "indexes_kept": indexes,
        "remaining_columns_identical_to_source": True,
        "integrity_check": integrity,
        "output_sha256": sha256(args.output),
        "output_bytes": os.path.getsize(args.output),
    }
    with open(args.report, "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
