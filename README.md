# Hodge Khovanov Spectra

Code and data for reproducing the numerical spectral-gap figures in
"Quantum computing and Khovanov homology."

## Contents

- `src/spectral_kh/` - Khovanov complex and Hodge Laplacian implementation.
- `scripts/` - database generation scripts for knots, links, 14-crossing knots,
  and twisted unknots.
- `visualization/` - plotting scripts for the retained paper figures.
- `data/` - input planar-diagram tables.
- `databases/` - SQLite result databases used by the plotting scripts.
- `outputs/plots/` - retained figure files.

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

The paper also uses `10_DoS_Khovanov.png`; the script and source data for that
figure are not included in this repository.
