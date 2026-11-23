# Hodge Khovanov Spectra

Manuscript, code, and data for
"Quantum computing and Khovanov homology."

The [September 2026 manuscript snapshot](Paper/Current/main.pdf) is archived
in `Paper/Current/`. The numerical data and verification programs are listed below.

## Contents

- `Paper/Current/` - September 2026 manuscript snapshot and supporting files.
- [`checks/`](checks/) - exact checks of the alternating and split-union formulas and two spectral gaps for the knot diagram `9_46`.
- [`checks/database/`](checks/database/) - November 2026 verification of the stored Betti numbers and gaps against
  KnotInfo, LinkInfo and exact ranks, with the programs, batch files and results.
- [`Paper/Current/reproducibility/`](Paper/Current/reproducibility/) - plotting script, numerical tables,
  database cleanup script and recovery record, and computational diagnostics.
- `src/spectral_kh/` - Khovanov complex and Hodge Laplacian implementation, with exact-rank Betti numbers in `exact_rank.py`.
- `tests/` - regression tests for the Betti numbers and gaps of small knots.
- `scripts/` - database generation scripts for knots, links, 14-crossing knots,
  and twisted unknots.
- `visualization/` - earlier plotting scripts, retained for reference.
- `data/` - input planar-diagram tables.
- `databases/` - SQLite result databases used by the plotting script. The knot, link and twisted-unknot databases
  were verified and repaired in November 2026; each has a `.REPAIR.json` record beside it.
- `outputs/plots/` - earlier figure files, retained for reference.

## Current numerical figures

The four numerical figures are in `Paper/Current/figures/`.

- `twisted_unknot_bidegree_gap.png`
- `observed_knots_gap_distribution.png`
- `observed_links_gap_distribution.png`
- `alternating_distinguished_bidegrees.png`

The paper folder also includes the other three diagram files needed to compile
the manuscript.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

The database-generation scripts import the installed `spectral_kh` package
directly. No repository-specific environment setup script is required.

## Generate the current numerical figures

```bash
python Paper/Current/reproducibility/plot_submission_numerics.py
```

This requires Python 3.11 or later and uses the included databases without
rerunning the eigensolvers. It writes the four figures and their summary
tables, checking the input database hashes recorded with the paper.

## Run the exact checks

Run these commands from the repository root with Python 3.10 or later and
assertions enabled. They use only the Python standard library.

```bash
python checks/check_alternating.py
python checks/check_split.py
python checks/check_gap_9_46.py
```

The regression tests for the engine's Betti numbers and gaps run with

```bash
python tests/test_exact_rank.py
```

The first script compares the spanning-tree and orientation formulas with a
direct rational Khovanov cube computation on eleven connected alternating
diagrams. The second checks the split-union formulas on three examples.
Their expected outputs are [output.json](checks/output.json) and
[split_output.json](checks/split_output.json). The third script uses exact
integer arithmetic to bound the first positive eigenvalue in bidegrees
`(-1,-1)` and `(0,-1)` of the knot diagram `9_46` between 0.0289984 and
0.0289985.

## Data Notes

The SQLite databases are included so readers can inspect the numerical data
without rerunning the full computation. The Python scripts in `scripts/` are the
database-generation entry points retained for reference and reruns.

`databases/knot_research.db` contains 12,965 knots and 1,401,192 bidegree rows,
with one row per knot and bidegree. It was made in two steps. A cleanup removed
duplicate rows, preserving retained row IDs and adding a unique index; its
selection rule, checksums, and removed-row recovery record are in
`Paper/Current/reproducibility/DATABASE-CLEANUP.json`. A repair in November
2026 then corrected 1,774 Betti numbers on 810 twelve-crossing knots that an
early eigensolver had undercounted, rebuilt those knots' Khovanov polynomials,
and stored every Betti number as an integer. Every Betti number in the file
now agrees with KnotInfo. No spectral gap changed. `link_research.db` and
`twisted_unknot_research.db` were checked the same way; their values were
already correct and only flags and storage types changed. Each database lists
its changes in the tables `row_repairs` and `polynomial_repairs`, and the full
account is in
[`NUMERICAL-REVIEW.md`](Paper/Current/reproducibility/NUMERICAL-REVIEW.md).

The files before the repair are at commit
`9a06f279dc323ec97d996e21b98ee915020f32a6`.

To keep the files small, the `log_torsion` and `compute_time` columns have
been removed from every database in `databases/`. Bidegree labels, chain
dimensions, Betti numbers, spectral gaps and all other tables are unchanged,
and each file has a `.SLIM.json` record with its checksums. The files with
those columns are at commit `5d61309a76522d5e8c59197de1d6e72fe861979c`.

The original database remains available at commit
`ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea`. To repeat the cleanup, provide that
original file to `Paper/Current/reproducibility/clean_knot_database.py` with
`--source`, `--output`, and `--report`, using new output paths.

`databases/knot_research_14_20261007.db` is a later snapshot of the
14-crossing computation, taken on 2026-10-07: 2,155 knots, of which 2,153 are
complete, and 277,468 bidegree rows, with one row per knot and bidegree and
Betti numbers stored as integers. `14ah_00002` and `14ah_00010` are still
partial. The paper's figures and tables do not use this file; they use
`databases/knot_research_14.db`, whose values are unchanged.
