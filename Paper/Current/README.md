# Manuscript

"Quantum computing and Khovanov homology", September 2026 version.

- [`main.pdf`](main.pdf) is the compiled manuscript.
- [`main.tex`](main.tex) is the complete source: preamble, sections, diagram
  macros and appendix.

## Compiling

Keep `figures/`, the two `.bib` files, `quantumarticle.cls` and
`schmidhuber.sty` beside `main.tex`. From this folder run

```
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

On Windows, `build.ps1` runs the same steps in a `build/` subfolder, stops on
undefined references or overfull boxes, and copies the result to `main.pdf`.

The databases are not needed to compile the paper.

## Numerical figures and data

The four numerical figures in `figures/` are produced from the databases in
`../../databases/` by

```
python reproducibility/plot_submission_numerics.py
```

[`reproducibility/NUMERICAL-REVIEW.md`](reproducibility/NUMERICAL-REVIEW.md)
describes where the data come from, how they were verified, and their limits.

The manuscript cites the repository at commit `ccfce9a3eaa3`. The knot
database has since had duplicate rows removed and 1,774 Betti numbers
corrected. Neither step changed a spectral gap, and the figures are
pixel-identical before and after the correction.
