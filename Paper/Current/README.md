# Current manuscript

Edit [main.tex](main.tex). It contains the entire paper, including the preamble,
authors, sections, diagram macros, and appendix.

Read [main.pdf](main.pdf). Run `build.ps1` in PowerShell after editing to rebuild
this PDF and its bibliography. In a LaTeX editor, select `main.tex` and use
pdfLaTeX + BibTeX. The `build` folder holds generated compilation files.

On another computer, let Dropbox finish downloading the whole `Current` folder
before opening `main.tex`. In the editor, compile with pdfLaTeX, run BibTeX,
then run pdfLaTeX twice. The equivalent terminal commands, run from this folder,
are `pdflatex main.tex`, `bibtex main`, and `pdflatex main.tex` twice.

The 12 September edit addresses the probability-command error reported by
the synced TeX Live 2024 log. If the error persists after syncing, retain the
new log so it can be checked against the updated source.

Keep `figures`, the two `.bib` files, `quantumarticle.cls`, and
`schmidhuber.sty` alongside the manuscript. The numerical summaries and scripts
are in `reproducibility`. In the GitHub checkout, the cleaned knot database is
at `../../databases/knot_research.db`; in the Dropbox working folder it is at
`reproducibility/data/knot_research.db`. Its cleanup record is
`reproducibility/DATABASE-CLEANUP.json`. Run
`python reproducibility/plot_submission_numerics.py` to regenerate the four
numerical figures and their summary tables. The script supports both layouts.
The databases are not needed to compile the paper.

This is the single-file R5 manuscript, promoted to the stable current location
on 11 September 2026. Its text was unchanged during the folder reorganization.
Sascha's opening and the mathematical corrections are retained.

Earlier drafts and exports remain in the Dropbox `Paper/Backup` folder.
GitHub carries this current version; future edits belong here.
