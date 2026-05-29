# Hodge Khovanov Spectra

Code and data for reproducing the numerical spectral-gap figures in
"Quantum computing and Khovanov homology."

This repository is intentionally narrow: it contains only the public artifacts
needed for the paper figures listed below. Exploratory studies, local CARC job
scripts, diagnostic notebooks, manuscript PDFs, and unrelated derived outputs
are not included.

## Contents

- `src/spectral_kh/` - Khovanov complex and Hodge Laplacian implementation.
- `scripts/` - database generation scripts for knots, links, 14-crossing knots,
  and twisted unknots.
- `visualization/` - plotting scripts for the retained paper figures.
- `data/` - input planar-diagram tables.
- `databases/` - SQLite result databases used by the plotting scripts.
- `outputs/plots/` - retained figure files.

## Retained Figures

| Figure file | Generation script | Source database |
| --- | --- | --- |
| `outputs/plots/twisted_unknot_gap.png` | `visualization/plot_twisted_unknot_gap.py` | `databases/twisted_unknot_research.db` |
| `outputs/plots/unified_plot_knots.png` | `visualization/plot_unified_knots.py` | `databases/knot_research.db`, `databases/knot_research_14.db` |
| `outputs/plots/unified_plot_links.png` | `visualization/plot_unified_links.py` | `databases/link_research.db` |

The paper also uses `10_DoS_Khovanov.png`; the script and source data for that
figure are not included in this repository.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

The database-generation scripts import the installed `spectral_kh` package
directly. No repository-specific environment setup script is required.

## Generate Retained Figures

```bash
python visualization/plot_twisted_unknot_gap.py
python visualization/plot_unified_knots.py
python visualization/plot_unified_links.py
```

These generate:

- `outputs/plots/twisted_unknot_gap.png`
- `outputs/plots/unified_plot_knots.png`
- `outputs/plots/unified_plot_links.png`

## Data Notes

The SQLite databases are included so readers can inspect the numerical data
without rerunning the full computation. The Python scripts in `scripts/` are the
database-generation entry points retained for reference and reruns.

`databases/knot_research.db` is about 89 MB, which is below GitHub's hard file
limit but above its recommended 50 MB threshold.
