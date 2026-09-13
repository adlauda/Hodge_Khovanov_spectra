# Hodge Khovanov Spectra

Manuscript, code, and data for
"Quantum computing and Khovanov homology."

The current paper is [Paper/Current/main.pdf](Paper/Current/main.pdf).
Edit [Paper/Current/main.tex](Paper/Current/main.tex), which contains the full
manuscript in one file. The same folder contains the bibliography, figures,
and instructions for compiling the paper and reproducing its numerical plots.

## Contents

- `Paper/Current/` - current manuscript, PDF, and supporting files.
- `Paper/Current/reproducibility/` - current plotting script, numerical tables,
  database cleanup script and recovery record, and computational diagnostics.
- `src/spectral_kh/` - Khovanov complex and Hodge Laplacian implementation.
- `scripts/` - database generation scripts for knots, links, 14-crossing knots,
  and twisted unknots.
- `visualization/` - earlier plotting scripts, retained for reference.
- `data/` - input planar-diagram tables.
- `databases/` - SQLite result databases used by the plotting scripts.
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

## Data Notes

The SQLite databases are included so readers can inspect the numerical data
without rerunning the full computation. The Python scripts in `scripts/` are the
database-generation entry points retained for reference and reruns.

`databases/knot_research.db` contains 12,965 knots and 1,401,192 bidegree rows,
with one row per knot and bidegree. The cleanup preserves retained row IDs
and adds a unique index to prevent duplicate bidegrees. The selection rule,
checksums, and removed-row recovery record are in
`Paper/Current/reproducibility/DATABASE-CLEANUP.json`.

The original database remains available at commit
`ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea`. To repeat the cleanup, provide that
original file to `Paper/Current/reproducibility/clean_knot_database.py` with
`--source`, `--output`, and `--report`, using new output paths.
