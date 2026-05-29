#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

try:
    from spherogram import Link
except ImportError:
    print("Error: spherogram is not installed.", file=sys.stderr)
    sys.exit(1)

INPUT_FILES = ["data/crossings/14/14a-hyp.csv", "data/crossings/14/14n-hyp.csv"]
OUTPUT_FILE = "data/crossings/14/knotinfo_14_crossing_pd.csv"


def pd_to_string(pd_code):
    return "PD[" + ",".join(
        "X[" + ",".join(str(x) for x in crossing) + "]"
        for crossing in pd_code
    ) + "]"


def convert_dt_to_pd(dt_code):
    # Older spherogram versions often use the constructor directly.
    L = Link(dt_code)
    return pd_to_string(L.PD_code())


def process_file(input_file, writer):
    count = 0
    failures = 0

    with open(input_file, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            name = row["name"].strip()
            knot_sig = row["knot_sig"].strip()
            dt_code = row["dt_code"].strip()

            try:
                pd_code = convert_dt_to_pd(dt_code)
                writer.writerow({
                    "name": name,
                    "knot_sig": knot_sig,
                    "dt_code": dt_code,
                    "pd_code": pd_code,
                })
                count += 1
            except Exception as e:
                failures += 1
                print(f"Failed on {name} ({dt_code}): {e}", file=sys.stderr)

    return count, failures


def main():
    for fn in INPUT_FILES:
        if not Path(fn).exists():
            print(f"Missing input file: {fn}", file=sys.stderr)
            sys.exit(1)

    total_count = 0
    total_failures = 0

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["name", "knot_sig", "dt_code", "pd_code"]
        )
        writer.writeheader()

        for input_file in INPUT_FILES:
            count, failures = process_file(input_file, writer)
            total_count += count
            total_failures += failures
            print(f"{input_file}: wrote {count} rows, {failures} failures")

    print(f"Done. Wrote {total_count} rows to {OUTPUT_FILE}")
    print(f"Failures: {total_failures}")


if __name__ == "__main__":
    main()