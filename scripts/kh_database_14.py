import sqlite3
import pandas as pd
import time
import ast
import argparse
import os
import functools
import traceback

print = functools.partial(print, flush=True)

from spectral_kh import KhovanovEngine

DB_PATH = "databases/knot_research_14.db"
CSV_PATH = "data/knots14_pd.csv"

# A knot whose last_heartbeat is older than this many seconds is considered
# abandoned (the owner node likely died) and can be rescued by any other node.
STALE_SECONDS = 600  # 10 minutes


def init_db(conn):
    cursor = conn.cursor()
    cursor.execute("PRAGMA journal_mode=DELETE;")
    cursor.execute("PRAGMA synchronous=NORMAL;")
    cursor.execute("PRAGMA busy_timeout=30000;")

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS knots (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            pd_code TEXT,
            khovanov_polynomial TEXT,
            last_heartbeat REAL,
            claimed_by INTEGER
        )
    """)

    # Backwards-compatible migration for existing databases.
    for col, typ in [("last_heartbeat", "REAL"), ("claimed_by", "INTEGER")]:
        try:
            cursor.execute(f"ALTER TABLE knots ADD COLUMN {col} {typ}")
        except sqlite3.OperationalError:
            pass  # Column already exists

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS bidegrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            knot_id INTEGER,
            h INTEGER,
            q INTEGER,
            dimension INTEGER,
            betti INTEGER,
            smallest_nonzero REAL,
            compute_time REAL,
            FOREIGN KEY (knot_id) REFERENCES knots (id),
            UNIQUE(knot_id, h, q)
        )
    """)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_bidegrees_knot_id ON bidegrees(knot_id)")
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_knots_active ON knots(khovanov_polynomial, last_heartbeat)"
    )
    conn.commit()


def open_conn():
    conn = sqlite3.connect(DB_PATH, timeout=60)
    init_db(conn)
    return conn


def parse_pd_code(raw_str):
    try:
        if pd.isna(raw_str):
            return None
        clean = str(raw_str).strip().replace(";", ",").replace("{", "[").replace("}", "]")
        clean = clean.replace(",]", "]").replace(", ]", "]")
        return ast.literal_eval(clean)
    except Exception:
        return None


def format_kh_poly(betti_dict):
    terms = []
    for (h, q), betti in sorted(betti_dict.items()):
        term = f"q^{q}t^{h}"
        if betti > 1:
            term = f"{betti}{term}"
        terms.append(term)
    return " + ".join(terms) if terms else "0"


def reconnect(conn):
    try:
        conn.close()
    except Exception:
        pass
    return open_conn()


def update_heartbeat(conn, knot_id, array_id):
    """Stamp a fresh heartbeat on a knot we are actively computing."""
    cursor = conn.cursor()
    for attempt in range(5):
        try:
            cursor.execute(
                "UPDATE knots SET last_heartbeat = ?, claimed_by = ? WHERE id = ?",
                (time.time(), array_id, knot_id),
            )
            conn.commit()
            return
        except sqlite3.OperationalError as e:
            if "locked" in str(e).lower() and attempt < 4:
                time.sleep(1)
            else:
                return


def try_claim_partial(conn, array_id):
    """Find and claim an abandoned partial knot from anywhere in the DB."""
    cursor = conn.cursor()
    now = time.time()
    stale_before = now - STALE_SECONDS

    try:
        row = cursor.execute(
            """
            SELECT id, name, pd_code
            FROM knots
            WHERE khovanov_polynomial IS NULL
              AND (last_heartbeat IS NULL OR last_heartbeat < ?)
            ORDER BY (last_heartbeat IS NULL) DESC, last_heartbeat ASC
            LIMIT 1
            """,
            (stale_before,),
        ).fetchone()
    except sqlite3.OperationalError:
        return None

    if row is None:
        return None

    knot_id, name, pd_code_str = row

    try:
        cursor.execute(
            "UPDATE knots SET last_heartbeat = ?, claimed_by = ? WHERE id = ?",
            (now, array_id, knot_id),
        )
        conn.commit()
    except sqlite3.OperationalError:
        return None

    pd_code = parse_pd_code(pd_code_str)

    if pd_code is None:
        return None

    return (name, pd_code)


def process_one_knot(conn, knot_name, pd_code, array_id=0, heartbeat_every=5):
    cursor = conn.cursor()
    max_retries = 20
    knot_id = None
    already_done = None

    for attempt in range(max_retries):
        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO knots (name, pd_code, last_heartbeat, claimed_by)
                VALUES (?, ?, ?, ?)
                """,
                (knot_name, str(pd_code), time.time(), array_id),
            )
            conn.commit()

            row = cursor.execute(
                "SELECT id, khovanov_polynomial FROM knots WHERE name=?",
                (knot_name,),
            ).fetchone()

            if row is None:
                return conn, "failed", "could not fetch knot row"

            knot_id, already_done = row
            break

        except sqlite3.OperationalError as e:
            if "locked" in str(e).lower() and attempt < max_retries - 1:
                sleep_time = 2 + (os.getpid() % 5)
                print(f"[Node {array_id}] DB lock on init for {knot_name}; sleeping {sleep_time}s")
                time.sleep(sleep_time)
                conn = reconnect(conn)
                cursor = conn.cursor()
            else:
                return conn, "failed", f"db init error: {e}"

    if already_done is not None:
        return conn, "skip_done", "already completed"

    update_heartbeat(conn, knot_id, array_id)

    try:
        # CHECKPOINT RESUME
        existing_rows = cursor.execute(
            "SELECT h, q, betti FROM bidegrees WHERE knot_id = ?",
            (knot_id,),
        ).fetchall()

        existing_bidegrees = {(h, q): betti for h, q, betti in existing_rows}
        betti_accumulator = {k: v for k, v in existing_bidegrees.items() if v > 0}

        engine = KhovanovEngine(pd_code)
        total_bidegrees = len(engine.chain_groups)

        if len(existing_bidegrees) > 0:
            print(
                f"[Node {array_id}] {knot_name}: RESUMING! "
                f"{len(existing_bidegrees)}/{total_bidegrees} bidegrees already banked."
            )
        else:
            print(f"[Node {array_id}] {knot_name}: {total_bidegrees} bidegrees to compute")

        for j, ((h, q), dim) in enumerate(engine.chain_groups.items(), start=1):

            if (h, q) in existing_bidegrees:
                if j % heartbeat_every == 0 or j == total_bidegrees:
                    print(
                        f"[Node {array_id}] {knot_name}: skipped {j}/{total_bidegrees} "
                        f"(Already banked)"
                    )
                continue

            t0 = time.time()
            betti, gaps = engine.get_low_spectrum_adaptive(h, q, tol=1e-7)
            gap = float(gaps[0]) if len(gaps) > 0 else None
            elapsed = time.time() - t0

            if betti > 0:
                betti_accumulator[(h, q)] = betti

            for attempt in range(max_retries):
                try:
                    cursor.execute(
                        """
                        INSERT OR REPLACE INTO bidegrees
                        (knot_id, h, q, dimension, betti, smallest_nonzero, compute_time)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (knot_id, h, q, len(dim), betti, gap, elapsed),
                    )
                    conn.commit()
                    break
                except sqlite3.OperationalError as e:
                    if "locked" in str(e).lower() and attempt < max_retries - 1:
                        sleep_time = 2 + (os.getpid() % 5)
                        time.sleep(sleep_time)
                        conn = reconnect(conn)
                        cursor = conn.cursor()
                    else:
                        raise e

            update_heartbeat(conn, knot_id, array_id)

            if j % heartbeat_every == 0 or j == total_bidegrees:
                print(
                    f"[Node {array_id}] {knot_name}: completed {j}/{total_bidegrees} bidegrees"
                )

        poly_str = format_kh_poly(betti_accumulator)

        for attempt in range(max_retries):
            try:
                cursor.execute(
                    """
                    UPDATE knots
                    SET khovanov_polynomial = ?, claimed_by = NULL
                    WHERE id = ?
                    """,
                    (poly_str, knot_id),
                )
                conn.commit()
                return conn, "success", f"completed all {total_bidegrees} bidegrees"
            except sqlite3.OperationalError as e:
                if "locked" in str(e).lower() and attempt < max_retries - 1:
                    sleep_time = 2 + (os.getpid() % 5)
                    time.sleep(sleep_time)
                    conn = reconnect(conn)
                    cursor = conn.cursor()
                else:
                    raise e

    except Exception as e:
        return conn, "failed", f"compute error: {e}\n{traceback.format_exc()}"


def main(args):
    conn = open_conn()

    df = pd.read_csv(CSV_PATH)
    df.columns = [c.strip() for c in df.columns]

    if "name" not in df.columns or "pd_code" not in df.columns:
        raise ValueError(f"Expected columns ['name', 'pd_code'], got {list(df.columns)}")

    total_csv = len(df)

    finished_df = pd.read_sql_query(
        "SELECT name FROM knots WHERE khovanov_polynomial IS NOT NULL",
        conn
    )
    finished_set = set(finished_df["name"].tolist())

    df = df[~df["name"].isin(finished_set)].reset_index(drop=True)
    unfinished = len(df)

    if unfinished == 0:
        print("No remaining knots to process.")
        conn.close()
        return

    sample_n = min(args.sample_size, unfinished)
    df = df.sample(n=sample_n, random_state=args.seed).reset_index(drop=True)

    my_df = df[df.index % args.total_bins == args.array_id].reset_index(drop=True)
    
    # Sort in-bin partials to the top of this node's queue.
    partial_df = pd.read_sql_query(
        "SELECT name FROM knots WHERE khovanov_polynomial IS NULL", conn
    )
    partial_set = set(partial_df["name"].tolist())
    my_df['is_partial'] = my_df['name'].isin(partial_set)
    my_df = my_df.sort_values(by='is_partial', ascending=False).reset_index(drop=True)
    my_df = my_df.drop(columns=['is_partial'])

    print("=" * 70)
    print(f"Node {args.array_id} summary")
    print(f"Database           : {DB_PATH}")
    print(f"CSV                : {CSV_PATH}")
    print(f"Total in CSV       : {total_csv}")
    print(f"Already done       : {len(finished_set)}")
    print(f"Currently unfinished: {unfinished}")
    print(f"Requested sample   : {args.sample_size}")
    print(f"Actual sample      : {sample_n}")
    print(f"Total bins         : {args.total_bins}")
    print(f"This array id      : {args.array_id}")
    print(f"Assigned to node   : {len(my_df)}")
    print(f"Seed               : {args.seed}")
    print(f"Heartbeat every    : {args.heartbeat_every} bidegrees")
    print(f"Stale threshold    : {STALE_SECONDS}s (rescue eligibility)")
    print("=" * 70)

    if len(my_df) == 0:
        print(f"Node {args.array_id}: nothing assigned.")
        conn.close()
        return

    success = 0
    failed = 0
    skipped = 0
    rescued = 0
    t_start = time.time()

    assigned_queue = list(my_df.itertuples(index=False))
    total_assigned = len(assigned_queue)
    assigned_done = 0

    while True:
        rescue = try_claim_partial(conn, args.array_id)

        if rescue is not None:
            knot_name, pd_code = rescue
            print(f"[Node {args.array_id}] RESCUING abandoned partial: {knot_name}")
            t0 = time.time()
            conn, status, msg = process_one_knot(
                conn, knot_name, pd_code, args.array_id, args.heartbeat_every
            )
            elapsed = time.time() - t0

            if status == "success":
                success += 1
                rescued += 1
                print(f"[Node {args.array_id}] RESCUE SUCCESS {knot_name} ({elapsed:.1f}s) | {msg}")
            elif status == "skip_done":
                skipped += 1
                print(f"[Node {args.array_id}] RESCUE SKIP {knot_name} ({elapsed:.1f}s) | {msg}")
            else:
                failed += 1
                print(f"[Node {args.array_id}] RESCUE FAILED {knot_name} ({elapsed:.1f}s) | {msg}")

            continue

        if not assigned_queue:
            break

        row = assigned_queue.pop(0)
        knot_name = row.name
        pd_code = parse_pd_code(row.pd_code)
        assigned_done += 1

        print(f"[Node {args.array_id}] {assigned_done}/{total_assigned} starting {knot_name}")
        t0 = time.time()

        conn, status, msg = process_one_knot(
            conn, knot_name, pd_code, args.array_id, args.heartbeat_every
        )

        elapsed = time.time() - t0

        if status == "success":
            success += 1
            print(f"[Node {args.array_id}] SUCCESS {knot_name} ({elapsed:.1f}s) | {msg}")
        elif status == "skip_done":
            skipped += 1
            print(f"[Node {args.array_id}] SKIP {knot_name} ({elapsed:.1f}s) | {msg}")
        else:
            failed += 1
            print(f"[Node {args.array_id}] FAILED {knot_name} ({elapsed:.1f}s) | {msg}")

        if assigned_done % args.report_every == 0 or assigned_done == total_assigned:
            total_elapsed = time.time() - t_start
            done_count = success + failed + skipped
            rate = done_count / total_elapsed if total_elapsed > 0 else 0.0
            print(
                f"[Node {args.array_id}] Progress {assigned_done}/{total_assigned} assigned | "
                f"success={success}, failed={failed}, skipped={skipped}, rescued={rescued} | "
                f"elapsed={total_elapsed/60:.1f} min | rate={rate:.3f} knots/s"
            )

    total_elapsed = time.time() - t_start
    print("=" * 70)
    print(f"Node {args.array_id} finished")
    print(f"Assigned   : {total_assigned}")
    print(f"Rescued    : {rescued}")
    print(f"Success    : {success}")
    print(f"Failed     : {failed}")
    print(f"Skipped    : {skipped}")
    print(f"Wall time  : {total_elapsed/60:.1f} min")
    print("=" * 70)

    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--total-bins", type=int, default=1)
    parser.add_argument("--sample-size", type=int, default=100)
    parser.add_argument("--seed", type=int, default=14)
    parser.add_argument("--report-every", type=int, default=5)
    parser.add_argument("--heartbeat-every", type=int, default=5)
    args = parser.parse_args()

    if args.array_id < 0 or args.array_id >= args.total_bins:
        raise ValueError("--array-id must satisfy 0 <= array-id < total-bins")

    main(args)