# Numerical review and reproducibility

Updated 2026-09-12. The numerical inputs start from repository snapshot
`ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea`. All four archived database files
match the Git blob hashes at that commit. The figures and summaries use the
cleaned knot database described below, with the other three databases unchanged.

## Local database cleanup

The raw knot database has 12,965 distinct knots. Its repeated bidegree rows
affect 552 knots, not additional knot types. There are 47,634 repeated
`(knot_id,h,q)` groups and 52,368 extra rows. All repeated chain dimensions
agree, and positive gaps differ by at most `7.925841538636291e-9`.

The Betti numbers disagree in 110 repeated bidegrees across 50 knots.
For example, `12a_3` at `(-1,-1)` has stored values 4 and 11. For every
bidegree of every knot, the highest-ID (last-inserted) row agrees with the
stored completed Khovanov polynomial, including all zero coefficients.
The working copy therefore retains that row. This reconciles the stored
records with their completion polynomial; it does not independently establish
which homology computation is mathematically correct. The 110 conflicts
remain explicitly recorded for possible independent recomputation.

The cleaned `knot_research.db` has 1,401,192 bidegree rows and no duplicate groups.
Every retained row, knot metadata row, and merge-log entry is unchanged.
A unique index on `(knot_id,h,q)` prevents new duplicate rows; it replaces
the redundant single-column `knot_id` index. The original published database
remains available at the pinned commit and in the Dropbox backup snapshots.
`DATABASE-CLEANUP.json` records the source and output checksums, every removed
row ID and its retained replacement, and the conflicting Betti records.
The cleanup passes SQLite integrity and foreign-key checks. It changes
per-diagram positive minima by at most `1.4337531162311734e-10`, with no change
to the rounded gap values quoted in the manuscript. Plots and summaries have
been regenerated from the working copy.

The cleaned copy preserves the existing Betti storage types (1,183,182 blobs
and 218,010 integers). Decode the eight-byte little-endian signed-integer blobs
before numerical Betti queries; the plotting script does not use this column.

## Code checks prompted by the manuscript comments

- Grading shifts are present in `src/spectral_kh/khovanov_core.py`, lines
  196-201. They change the bidegree labels, not the boundary matrices or their
  positive eigenvalues. `compute_laplacian`, lines 504-522, returns the sum of
  the upward and downward Laplacians without operator-norm rescaling.
- The base adaptive solver defaults to `tol=1e-9` (line 561). The original
  knot manager, link manager, and twisted-unknot manager use `1e-9`.
  `scripts/kh_database_manager_PE.py:142` and `scripts/kh_database_14.py:239`
  explicitly pass `tol=1e-7`. These are zero-classification thresholds, distinct
  from the sparse eigensolver's convergence tolerance. The published database
  does not record the exact driver and threshold for every historical run.
- `scripts/kh_database_14.py` randomly samples unfinished entries of
  `data/knots14_pd.csv`, using `DataFrame.sample` with a specified seed. The
  script performs no hyperbolicity test. The published table has 46,969
  diagrams (19,536 `14ah` and 27,433 `14nh` identifiers). Every one of the
  2,154 database PD codes matches its entry in that table. The input's Git
  blob hash is `a464594bd7d793e34e62a5b8135249b3bce9371b`.
  `scripts/crossings/14/make_14_pd_csv.py` identifies `14a-hyp.csv` and
  `14n-hyp.csv` as the upstream DT-code tables. Random selection from this
  supplied population is distinct from sampling uniformly from all knot types.
  The published snapshot still contains 2,133 completed and 21 partial records;
  there is no published fourteen-crossing completion update. The completed
  subset need not retain the inclusion probabilities of the initial draws.

All code references above are at the pinned public commit, available at
[Hodge Khovanov Spectra](https://github.com/adlauda/Hodge_Khovanov_spectra/tree/ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea).

## What is supplied

- `../main.tex`: the single active manuscript, including the numerical section.
- `clean_knot_database.py`: checked, non-destructive creation of the local
  working database; refuses to overwrite an existing output or cleanup record.
- Repository `databases/knot_research.db` and `DATABASE-CLEANUP.json`: cleaned
  working data and recovery record. The Dropbox layout instead keeps the
  working database at `data/knot_research.db` relative to this folder.
- `plot_submission_numerics.py`: deterministic plotting and summary script.
- `NUMERICS-MANIFEST.json`: pinned database sources, checksums, coverage,
  duplicate counts, storage issues, and figure checksums.
- `*_diagram_summary.csv`: one row per registered diagram, including completion
  state, counts of stored blocks, and smallest reported positive value.
- `knots_gap_summary.csv`, `links_gap_summary.csv`, and
  `alternating_gap_summary.csv`: exact statistics used in the figures.
- `frontier_diagnostics.json`: conservatively relabeled derivative of the
  earlier recomputation record. Its source SHA-256 is
  `73736d545eb63e15be35f9ed82cd246040260990f88570b19155f4c6be901747`.
- Repository `databases/`: the cleaned knot database and three unchanged
  input databases. The original files and KnotInfo diagram table are available
  at commit `ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea` and in the Dropbox
  `Paper/Backup/2026-09-08-R4/reproducibility/snapshots/` folder.
- `COMPLETION-VERIFICATION.json`: read-only comparison with the preceding
  snapshot, including preservation of all existing bidegree rows.

Commit `ccfce9a` contains the preceding database update and
`scripts/merge_missing_knot_results.py`. The cleanup documented here is a
subsequent derivative; `DATABASE-CLEANUP.json` identifies both database hashes.

Run `python reproducibility/plot_submission_numerics.py` from the submission
folder, or invoke that script by its path from any directory. The script uses
the supplied cleaned knot database and verifies its cleanup checksum. To
recreate the cleaned database in an empty output directory, run
`clean_knot_database.py` with the original database from the pinned commit;
`--source`, `--output`, and `--report` select explicit paths. Never overwrite
the archived source. Requires Python
3.11+, NumPy, pandas, and Matplotlib; revision 4 uses Python 3.12, NumPy
2.3.5, pandas 2.3.3, Matplotlib 3.10.6. The script verifies the working database
hashes before reading them with SQLite's read-only mode, and checks the
recorded cleanup source hash against the pinned original hash. Figure generation
only summarizes existing data. The CSV diagram table is used for the separate
orientation check below, not for plotting.

## Exact figures to use

| File under `../figures/` | Meaning |
|---|---|
| `observed_knots_gap_distribution.png` | One stored positive minimum per completed knot; medians, central 50%/80%, observed minimum; omitted-record counts stated. |
| `observed_links_gap_distribution.png` | Same summaries, each link weighted equally. |
| `twisted_unknot_bidegree_gap.png` | Exact gap in the specified unshifted unreduced bidegree compared with observed full-diagram minima. |
| `alternating_distinguished_bidegrees.png` | Completed alternating records compared with exact distinguished-bidegree formulas; absent diagrams marked as formulas. |

All four PNGs were visually inspected. The old unified plots should not be used
in the submission: their raster titles and legends still say “Global Worst-Case,”
and their means pool repeated bidegree rows. The new figures fit no asymptotic
curves and make no confidence-interval claim for their quantile bands.

## Published-snapshot coverage and historical data issues

The following counts describe the unchanged published snapshot, before the
local cleanup described above. Direct read-only SQLite queries found:

- Knots through 13 crossings: all 12,965 registered records are marked complete
  and have stored bidegrees. The 37 additions comprise 3 knots at 10 crossings,
  8 at 11, 9 at 12, and 17 at 13. The database has 1,453,560 bidegree rows,
  exactly 3,689 more than before. Every preceding bidegree row is unchanged,
  as is every previously completed knot row. There are 37 merge-log entries,
  no missing polynomial values, and `PRAGMA integrity_check` returns `ok`.
  These checks do not independently certify the spectra or Betti numbers.
- Fourteen-crossing subset: 2,154 registered; 2,133 marked complete; 21 partial.
  Plotting uses the completed records only. Completion is a stored-polynomial
  flag, not an independently validated enumeration of all expected blocks.
- Links: all 1,424 registered records are marked complete.
- Twisted unknots: all 13 records, n=1,...,13, are marked complete.
- Knot database: 47,634 duplicated `(knot_id,h,q)` groups and 52,368 extra rows.
  Duplicate positive gaps differ by at most `7.925841538636291e-9`; duplicated
  groups have identical dimensions. There are no duplicate bidegree groups in
  the other three databases. Per-diagram minima avoid giving repeated rows
  additional statistical weight.
- `betti` is stored as a blob in 1,199,591 knot rows, all 275,896 fourteen-crossing
  rows, 88 link rows, and 141 twisted-unknot rows. The plots do not read `betti`.
  SQL numerical filters or sums on that column are unsafe without type decoding.
- Null gaps occur in all datasets and cannot in general be equated with solver
  failure; some blocks legitimately have no positive spectrum. The raw database
  lacks complete per-row failure and solver records.

The overall fourteen-crossing minimum remains `0.0006119790255380953` when
partial records are removed. The alternating minimum changes: including partial
records gives `0.08259883661804181`, whereas the 886 completed alternating
records have minimum `0.08654902144568334`, attained by `14ah_00033`.
Adding `T(2,13)=13a_4878` lowers the 13-crossing alternating minimum from
`0.06797347509123662` at `13a_4874` to `0.05811636514789631`, agreeing
numerically with its distinguished-block formula. The overall nonalternating
frontier values in the paper remain unchanged.

The user reports that source-shard matching, repeated merges, and replacement
merges were tested successfully in an in-memory copy before the production
merge. Those tests were not rerun here. The published merge script updates
existing rows and inserts missing rows within `BEGIN IMMEDIATE`, preserving
row IDs and fields absent from the shards. This does not remove historical
duplicate groups, whose count is unchanged. Optional fourteen-crossing
completion, a whole-database Jones consistency check, and historical Betti-blob
repair remain undone. The new rows store integer Betti values, bringing that
count to 253,969; the 1,199,591 historical blob values are unchanged.

## What the earlier “certification” actually proves

For integer differentials A=d_out, B=d_in with AB=0, one has
`beta_Q = dim - rank_Q(A) - rank_Q(B)`. Ranks can only drop modulo a prime.
Thus a modular rank sum equal to the chain dimension proves beta_Q=0. All ten
displayed nonalternating frontier blocks have this pattern, and one prime is
already enough for that conclusion. The earlier script did not explicitly
record a separate exact check of AB=0, so its reliance on the engine's chain
construction should be acknowledged.

Positive modular nullity is different: the two-prime nullity 1 reported for
`13a_4874` establishes only beta_Q <= 1. Two primes agreeing does not prove
beta_Q=1. A rational kernel vector or another exact lower bound would be needed.

A small residual bounds the distance of an approximate eigenvalue to *some*
true eigenvalue when the residual is evaluated exactly. It does not show the
eigenvalue is the smallest positive one: for `L=diag(epsilon,1)`, `v=(0,1)` and
`lambda=1` have zero residual even though lambda_min=epsilon. Floating residuals
also lack rigorous outward rounding. Independent eigensolvers agreeing is good
numerical evidence, but the earlier JSON field `certified_interval` is not a
validated interval for the gap. The derivative supplied here removes that field
and labels the residuals as diagnostics. No certified gap or phase-gap promise
is inferred from them.

Also fix the phrase “both primes exceeding 2^61”: the first prime is `2^61-1`.

## Scope of the exact companion formulas

[Grlj and Lauda, Spectral geometry of Khovanov Laplacians](https://arxiv.org/html/2608.24298v1)
proves bidegree statements in Theorems 4.26, 4.29, and 4.32. Corollary 4.28,
Remark 4.30, and the conclusion explicitly leave the whole-complex minimum
question open. The formulas must therefore be labeled as distinguished-bidegree
gaps, even when observed minima agree.

The original CSV at the pinned commit has SHA-256
`f5d0a015177423c67e5fa2c428e2b6524ad04a4f9302c2c3f8ad2737d221340b`
(2,713,960 bytes). A lightweight Spherogram orientation check on that table
found all positive crossings for the torus representatives `3_1`, `5_1`, `7_1`,
`9_1`, `11a_367`, and `13a_4878`, consistent with the positive-braid formula.
The familiar even twist representatives `4_1`, `6_1`, `8_1`, and `10_1` have
`n-2` positive and 2 negative crossings, consistent with the companion's stated
shift. This sign check is not a proof that every tied diagram is identical to
the companion's standard representative. No family bidegree is assigned by
name to a fourteen-crossing census diagram. In the replacement section the
exact formulas explicitly refer to the companion's oriented representatives;
the empirical statistics minimize over the database's stored bidegrees.

## Earlier support scripts need repair before reuse

The following observations concern earlier support scripts reviewed on
2026-09-07, not the newly published merge script. The completed merge does not
establish that these separate historical concerns have been repaired:

- `certify_frontier.py` describes two-prime agreement as exact characteristic-zero
  rank for arbitrary positive nullity, and calls residual intervals certified.
  Its resume code expects `results.json` to be a list, but the supplied file is
  a dictionary with `blocks` and `cross_check_14nh_03671` keys, so resume fails.
- `complete_missing_knots.py` skips previously inserted rows even if they record
  a failed solve; on resumption its accumulated polynomial includes only newly
  processed rows and can overwrite the full polynomial incorrectly. LOBPCG
  fallback does not inspect residual/convergence evidence before accepting a
  result. Do not launch this script as a correctness-preserving repair workflow
  until these issues are fixed.
- The CARC note's proposed Euler-characteristic check is useful but insufficient
  to verify every Betti number: adjacent-degree errors can cancel. The public
  engine's `compute_jones_polynomial` itself uses numerical ranks, so it is not
  an independent exact Jones calculation.

The original engine, databases, and plotting scripts remain available at the
[pinned commit](https://github.com/adlauda/Hodge_Khovanov_spectra/tree/ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea).
This reproducibility folder contains the plotting script, summaries, cleanup
record, and diagnostics used by the current manuscript.
The archived data do not substantiate allocated CARC core-hours; those require
job-accounting records, so the replacement section omits those numbers.
