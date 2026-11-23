#!/usr/bin/env python3
"""Check data/knots14_pd.csv against the Dowker-Thistlethwaite tables.

For every name in data/crossings/14/14a-hyp.csv and 14n-hyp.csv, the planar
diagram stored under that name in data/knots14_pd.csv is read off as a
Dowker-Thistlethwaite code from each of its 28 base points, in both
directions, and the table's code must be among the results up to an overall
change of sign.  A change of sign is a mirror image, which the tables do not
distinguish.

The test is exact and combinatorial.  Standard library only.

Usage, from the repository root:

    python scripts/crossings/14/check_14_pd_table.py
"""

import ast
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PD_TABLE = ROOT / "data/knots14_pd.csv"
DT_TABLES = [ROOT / "data/crossings/14/14a-hyp.csv",
             ROOT / "data/crossings/14/14n-hyp.csv"]


def dt_from_letters(code):
    """'bhN...' -> (4, 16, -28, ...): a=2, b=4, ..., capitals negative."""
    return tuple(2 * (ord(c.lower()) - 96) * (1 if c.islower() else -1) for c in code)


def passings(pd):
    """Pair the two passings of each crossing.

    A crossing [a, b, c, d] has incoming under-strand a and the other strands
    counterclockwise.  Strand k is followed by strand k + 1 (cyclically), so
    the over-strand enters along b if d follows b, and along d otherwise.
    Passing k is the arrival of strand k at a crossing.  Returns
    partner[k] and over[k] (True if passing k is an overpass), 0-based.
    """
    m = 2 * len(pd)
    partner = [None] * m
    over = [None] * m
    for a, b, c, d in pd:
        if c % m != (a + 1) % m:
            raise ValueError(f"under-strand {a} -> {c} is not consecutive")
        if d % m == (b + 1) % m:
            o = b
        elif b % m == (d + 1) % m:
            o = d
        else:
            raise ValueError(f"over-strand {b}, {d} is not consecutive")
        u, o = (a - 1) % m, (o - 1) % m
        if partner[u] is not None or partner[o] is not None:
            raise ValueError("a strand arrives at two crossings")
        partner[u], partner[o] = o, u
        over[u], over[o] = False, True
    if None in partner:
        raise ValueError("a strand arrives at no crossing")
    return partner, over


def dt_codes(pd):
    """All Dowker-Thistlethwaite codes of the diagram, one sign convention."""
    partner, over = passings(pd)
    m = len(partner)
    codes = set()
    for start in range(m):
        for step in (1, -1):
            # new label (1-based) of old passing p, and the inverse
            label = {(start + step * t) % m: t + 1 for t in range(m)}
            old = {v: k for k, v in label.items()}
            code = []
            for odd in range(1, m, 2):
                p = old[odd]
                even = label[partner[p]]
                if even % 2:
                    raise ValueError("odd passing paired with odd passing")
                code.append(-even if over[partner[p]] else even)
            codes.add(tuple(code))
    return codes


def main():
    with open(PD_TABLE, newline="") as f:
        pd_rows = list(csv.DictReader(f))
    pds = {r["name"]: r["pd_code"] for r in pd_rows}
    if len(pds) != len(pd_rows):
        sys.exit("duplicate names in the planar-diagram table")

    checked = 0
    bad = []
    names = set()
    for table in DT_TABLES:
        with open(table, newline="") as f:
            for row in csv.DictReader(f):
                name = row["name"].strip()
                names.add(name)
                want = dt_from_letters(row["dt_code"].strip())
                try:
                    codes = dt_codes(ast.literal_eval(pds[name]))
                except (KeyError, ValueError) as exc:
                    bad.append((name, repr(exc)))
                    continue
                if want not in codes and tuple(-x for x in want) not in codes:
                    bad.append((name, "code not found"))
                checked += 1

    extra = sorted(set(pds) - names)
    print(f"diagrams in {PD_TABLE.name}: {len(pds)}")
    print(f"names in the Dowker-Thistlethwaite tables: {len(names)}")
    print(f"diagrams compared: {checked}")
    print(f"diagrams with no table entry: {len(extra)}")
    print(f"mismatches: {len(bad)}")
    for item in bad[:20]:
        print("  ", *item)
    sys.exit(1 if bad or extra or len(names) != len(pds) else 0)


if __name__ == "__main__":
    main()
