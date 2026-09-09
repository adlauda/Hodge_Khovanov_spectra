"""Merge complete independent missing-knot result shards into the primary database."""

import argparse
import sqlite3
import time
from pathlib import Path


def init_merge_table(connection):
    connection.execute("PRAGMA journal_mode=DELETE")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA busy_timeout=300000")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS merged_knot_results (
            knot TEXT PRIMARY KEY,
            source_db TEXT NOT NULL,
            merged_at REAL NOT NULL
        )
        """
    )
    connection.commit()


def complete_knots(shard):
    return shard.execute(
        """
        SELECT name, pd_code, khovanov_polynomial, turaev_genus
        FROM knot_results
        WHERE complete=1 AND khovanov_polynomial IS NOT NULL
        ORDER BY name
        """
    ).fetchall()


def merge_knot(primary, shard, source_path, result, replace):
    name, pd_code, polynomial, turaev_genus = result
    merged = primary.execute(
        "SELECT source_db FROM merged_knot_results WHERE knot=?", (name,)
    ).fetchone()
    if merged is not None and not replace:
        print(f"SKIP {name}: already merged from {merged[0]}", flush=True)
        return False

    knot_row = primary.execute("SELECT id FROM knots WHERE name=?", (name,)).fetchone()
    if knot_row is None:
        raise KeyError(f"Primary database has no knot row for {name}")
    knot_id = knot_row[0]
    rows = shard.execute(
        """
        SELECT h, q, dimension, betti, smallest_nonzero, compute_time
        FROM bidegree_results
        WHERE knot=?
        ORDER BY h, q
        """,
        (name,),
    ).fetchall()
    if not rows or any(betti is None for _, _, _, betti, _, _ in rows):
        raise ValueError(f"Shard {source_path} marks {name} complete with incomplete bidegrees")

    primary.execute("BEGIN IMMEDIATE")
    try:
        primary.execute(
            """
            UPDATE knots
            SET pd_code=COALESCE(pd_code, ?),
                turaev_genus=COALESCE(turaev_genus, ?),
                khovanov_polynomial=?
            WHERE id=?
            """,
            (pd_code, turaev_genus, polynomial, knot_id),
        )
        # Legacy databases have no UNIQUE(knot_id, h, q) constraint.
        # Update then insert under BEGIN IMMEDIATE to preserve existing IDs
        # and fields not supplied by the shard, without racing other writers.
        primary.executemany(
            """
            UPDATE bidegrees
            SET dimension=?, betti=?, smallest_nonzero=?, compute_time=?
            WHERE knot_id=? AND h=? AND q=?
            """,
            [(dim, betti, gap, seconds, knot_id, h, q)
             for h, q, dim, betti, gap, seconds in rows],
        )
        primary.executemany(
            """
            INSERT INTO bidegrees
                (knot_id, h, q, dimension, betti, smallest_nonzero, compute_time)
            SELECT ?, ?, ?, ?, ?, ?, ?
            WHERE NOT EXISTS (
                SELECT 1 FROM bidegrees WHERE knot_id=? AND h=? AND q=?
            )
            """,
            [(knot_id, *row, knot_id, row[0], row[1]) for row in rows],
        )
        primary.execute(
            """
            INSERT INTO merged_knot_results (knot, source_db, merged_at)
            VALUES (?, ?, ?)
            ON CONFLICT(knot) DO UPDATE SET
                source_db=excluded.source_db,
                merged_at=excluded.merged_at
            """,
            (name, str(source_path), time.time()),
        )
        primary.commit()
    except Exception:
        primary.rollback()
        raise
    print(f"MERGED {name}: {len(rows)} bidegrees", flush=True)
    return True


def shard_paths(values):
    paths = []
    for value in values:
        path = Path(value)
        paths.extend(sorted(path.glob("*.sqlite")) if path.is_dir() else [path])
    return paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-db", default="databases/knot_research.db")
    parser.add_argument("--results", nargs="+", required=True)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    primary = sqlite3.connect(args.primary_db, timeout=300)
    try:
        init_merge_table(primary)
        for path in shard_paths(args.results):
            if not path.exists():
                print(f"SKIP {path}: no result shard", flush=True)
                continue
            shard = sqlite3.connect(path)
            try:
                for result in complete_knots(shard):
                    merge_knot(primary, shard, path.resolve(), result, args.replace)
            finally:
                shard.close()
    finally:
        primary.close()


if __name__ == "__main__":
    main()
