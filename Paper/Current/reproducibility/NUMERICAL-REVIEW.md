# Numerical review and reproducibility

Updated 2026-11-23. The numerical inputs start from repository snapshot
`ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea`. All four archived database files
match the Git blob hashes at that commit. The figures and summaries use the
repaired knot, link and twisted-unknot databases now in `databases/`,
described in the next section, and the fourteen-crossing database. On
2026-11-23 the `log_torsion` and `compute_time` columns were removed from all
five published files; see "Column removal" below. No other value changed.

The sections after the next one were written on 2026-09-12 and describe the
state at that date. Statements the November work supersedes are marked
"Resolved 2026-11-22" in place.

## Verification and repair of 2026-11-22

Programs, batch files and result files are in `checks/database/`, indexed in
its `README.md`. The repaired files replace `knot_research.db`,
`link_research.db` and `twisted_unknot_research.db` in `databases/`. The files
they were derived from are those of commit
`9a06f279dc323ec97d996e21b98ee915020f32a6`; "before the repair" below
refers to that commit.

### What was wrong

The earliest database driver computed the low spectrum of bidegrees of
dimension 4,000 and above with `eigsh(k=50, which='SM')` and counted the
returned eigenvalues below a threshold. Lanczos iteration does not reliably
return every copy of a repeated eigenvalue, so the kernel of the Laplacian
could be undercounted. Rows written by that path carry `is_truncated=1`.

- The cleaned knot database had 7,439 such rows on 1,389 knots (55 rows at
  `11a`, 1 at `11n`, 4,738 at `12a`, 2,645 at `12n`). In 1,774 of them, on 810
  knots, the stored Betti number was too small, by between 1 and 13. The stored
  Khovanov polynomials of those 810 knots were built from the same rows and
  were wrong accordingly.
- The spectral gaps in those rows were correct. No recomputed gap differs from
  the stored one by more than a relative `1.6e-12`.
- The link database had 88 such rows and the twisted-unknot database 141. All
  of their Betti numbers and gaps were correct.
- The 110 conflicting duplicate records were instances of the same fault: in
  each, the discarded row was a truncated row with the smaller Betti number.

### Independent references

- **KnotInfo.** Every stored knot Betti number was compared with the free part
  of KnotInfo's integral Khovanov homology, which comes from programs unrelated
  to this repository. Before the repair 1,399,418 of 1,401,192 rows agreed;
  the 1,774 that did not were all truncated rows. In all 110 conflicts the
  retained value equals KnotInfo's and the removed value does not, so the
  cleanup rule chose correctly every time.
- **LinkInfo.** All 1,424 stored link tables equal LinkInfo's Khovanov
  polynomial for some orientation, up to mirror.
- **Unknot.** Each of the 13 twisted-unknot diagrams has total Betti number 2
  in one homological degree. The database labels its bidegrees with a constant
  offset from the current engine's normalisation, (5, 15) for `TU_10` and
  `TU_11` and (6, 18) for `TU_12` and `TU_13`; in the engine's labels the
  homology sits at h = 0, q = ±1. The stored labels were kept.
- **Euler characteristic.** For each diagram and quantum degree the
  alternating sum of Betti numbers must equal that of the chain dimensions.
  This failed for 804 knots before the repair and for none after. It is a
  necessary condition only; adjacent errors can cancel, as they did for 6 of
  the 810 knots.

### Recomputation

Every truncated row was recomputed by `checks/database/recompute_rows.py`
without classifying any eigenvalue by a threshold.

- `d_out d_in = 0` was checked in integer arithmetic for every block.
- The Betti number is `dim - rank(d_out) - rank(d_in)` with ranks modulo 65521.
  A modular rank cannot exceed the rational rank, so this is an upper bound on
  the rational Betti number; it equalled the KnotInfo or LinkInfo value in
  every row. The rank routine agrees with FLINT on a self-test.
- The gap is the smallest positive eigenvalue of the two differentials' Gram
  matrices, read at the index fixed by the exact rank. It was recomputed from
  the full Laplacian by shift-invert Lanczos at the index fixed by the Betti
  number.

| Database | Rows | Largest dimension | Exact ranks | Betti changed | Gap vs stored | Two gap methods | Zero leak / gap |
|---|---|---|---|---|---|---|---|
| knots | 7,439 | 18,592 | 7,439 | 1,774 | `1.6e-12` | `4.4e-9` | `1.1e-11` |
| links | 88 | 7,000 | 88 | 0 | `2.7e-14` | `6.0e-12` | `7.6e-14` |
| twisted unknots | 141 | 180,180 | 28 | 0 | `1.6e-14` | `2.9e-14` | `7.4e-14` |

The last three columns are the largest relative values over the rows. In 113
twisted-unknot rows a differential has more than 12,000 rows or columns and
dense exact ranks were not attempted. Those rows use the stored Betti number,
which the unknot argument above fixes, and two Lanczos runs on the Laplacian
from independent start vectors.

A first pass over the twisted unknots was invalid: it looked up blocks by the
database's labels in the engine's normalisation and so examined the wrong
blocks. The chain-dimension test rejected all 141 rows. The pass was repeated
with the labels matched by the unique constant offset that reproduces every
stored chain dimension, and the invalid output was discarded.

The exact rational Betti numbers of nine small frontier blocks were also
computed with FLINT integer ranks (`exact_blocks_knots.json`,
`exact_blocks_links.json`). In particular `13a_4874` at (10, 31) has Betti
number exactly 1.

### Repaired databases

`checks/database/apply_repairs.py` wrote the three repaired files. A recomputed row was accepted only if its chain
dimension matched, `d^2 = 0` held, its Betti number equalled the reference,
the two gap computations agreed to `1e-6` and the zero leak was below `1e-6`
of the gap. Every row passed.

| File | SHA-256 | Rows repaired | Betti changed | Polynomials rebuilt | Blobs made integers |
|---|---|---|---|---|---|
| `knot_research.db` | `37c52286…7a431` | 7,439 | 1,774 | 810 | 1,183,182 |
| `link_research.db` | `acf404fd…46f83` | 88 | 0 | 0 | 88 |
| `twisted_unknot_research.db` | `82d58cc6…c3c40` | 141 | 0 | 0 | 141 |

- No row carries `is_truncated=1`, and every `betti` value has integer storage.
  The flag now means that a row's Betti number and gap have not been verified;
  `row_repairs.old_is_truncated` keeps the old flag.
- The 7,930 knots that stored no PD code now have the code from
  `data/knotinfo.csv`; the 5,035 codes already stored equal that table's.
- Tables `row_repairs`, `polynomial_repairs` and `pd_code_fills` in each file
  list every old and new value, and the `.REPAIR.json` file beside it gives
  the full checksums and post-conditions. No row changed outside that log.
- After the repair all 12,965 knots agree with KnotInfo in all 1,401,192 rows,
  all 1,424 links agree with LinkInfo, and no diagram fails the Euler check
  (`checks/database/after_repair/`, run on the files as they now stand).
- Rebuilding the three files from `checks/database/recompute/` reproduces the
  recorded checksums, and removing the two columns from the rebuilt files
  reproduces the checksums of the published ones.

### Column removal of 2026-11-23

To reduce the size of the published files,
`checks/database/slim_database.py` removed two columns of the `bidegrees`
table from every database in `databases/`: `compute_time`, and `log_torsion`
where present (the fourteen-crossing databases never had it). No figure,
table or check in this repository reads either column. The script verified
that every remaining column of every table, and the list of indexes, is
identical before and after. Each file has a `.SLIM.json` record beside it
with both checksums.

| File | Before | After | SHA-256 after |
|---|---|---|---|
| `knot_research.db` | 84.2 MB | 72.1 MB | `3cfd6618…` |
| `knot_research_14.db` | 21.5 MB | 18.0 MB | `191fc2be…` |
| `knot_research_14_20261007.db` | 21.7 MB | 16.3 MB | `01b7d6ef…` |
| `link_research.db` | 8.1 MB | 6.5 MB | `c8255715…` |
| `twisted_unknot_research.db` | 0.10 MB | 0.09 MB | `5536e339…` |

- The files with both columns are at commit
  `5d61309a76522d5e8c59197de1d6e72fe861979c`.
- The checksums in the "Repaired databases" table above are those of the
  files before the columns were removed; each `.SLIM.json` record names that
  checksum as its source. For `knot_research_14.db` the source is the
  snapshot checksum `2c3e9d55…` used by the manuscript.
- The September cleanup kept, for each duplicated bidegree, the newest row.
  For 35,959 bidegrees on 443 knots that row had no `log_torsion` while a
  removed row did. Those values were therefore already absent from the
  cleaned knot database; they remain in the file at commit `ccfce9a3eaa3`.
- The KnotInfo, LinkInfo and Euler checks in `checks/database/after_repair/`
  were rerun on the files after the column removal, with the same results.
- The database drivers in `scripts/` still create and write both columns.
  `scripts/merge_missing_knot_results.py` writes `compute_time` and needs a
  database that has it.

### Effect on the manuscript's numbers

None. The figures regenerated from the repaired databases are pixel-identical
to the previous ones, and the gap summary tables agree to a relative `1.3e-12`.
The plots do not read the `betti` column.

### Engine and driver changes

- `src/spectral_kh/exact_rank.py` holds the modular rank, the rank-indexed gap
  and the self-test.
- `KhovanovEngine.get_low_spectrum_adaptive`, which every database driver
  calls, now takes its Betti number from exact ranks and reads the gap at that
  index. Its return format is unchanged. A bidegree whose differentials exceed
  25,000 rows or columns falls back to the previous threshold method with a
  `RuntimeWarning`.
- `KhovanovEngine.compute_khovanov_polynomial` called two functions that do
  not exist in the module and could not run; it now works.
- `tests/test_exact_rank.py` checks the full Betti tables and minimum gaps of
  `3_1`, `4_1`, `8_19` and `9_42`, and optionally a 4,000-dimensional block of
  `12n_803` that the truncated solver recorded as 5 instead of 7.

### Still open

- The fourteen-crossing databases were not recomputed. KnotInfo has no
  Khovanov data at 14 crossings, so there is no independent reference. All
  complete records pass the Euler check; their Betti values are still stored
  as blobs. The 2026-10-07 snapshot has 2 partial records, `14ah_00002` and
  `14ah_00010`.
- Rows the drivers wrote from a dense full spectrum (dimension below 4,000)
  were not recomputed. Their Betti numbers are covered by the KnotInfo and
  LinkInfo comparisons; their gaps rest on the dense eigensolver.
- The gaps remain floating-point values without validated error bounds.
- The support scripts criticised under "Earlier support scripts need repair
  before reuse" are not in this repository and were not changed.
- The manuscript was not edited. It cites the snapshot at commit
  `ccfce9a3eaa3` with a cleaned copy, and every gap it reports is unchanged.
  The phrase "both primes exceeding 2^61" noted below no longer occurs in
  `main.tex`.

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
*Resolved 2026-11-22: the retained value equals KnotInfo's in all 110 cases.*

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
*Resolved 2026-11-22: the repaired databases store every Betti value as an
integer.*

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
- `DATABASE-CLEANUP.json`: recovery record of the cleaned knot database, from
  which the repaired `databases/knot_research.db` was derived. The Dropbox layout instead keeps the
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
- Repository `databases/`: the repaired knot, link and twisted-unknot
  databases with their `.REPAIR.json` records, and the unchanged
  fourteen-crossing databases. The cleaned knot database and unrepaired link
  and twisted-unknot databases are at commit
  `9a06f279dc323ec97d996e21b98ee915020f32a6`. The original files and KnotInfo diagram table are available
  at commit `ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea` and in the Dropbox
  `Paper/Backup/2026-09-08-R4/reproducibility/snapshots/` folder.
- `COMPLETION-VERIFICATION.json`: read-only comparison with the preceding
  snapshot, including preservation of all existing bidegree rows.

Commit `ccfce9a` contains the preceding database update and
`scripts/merge_missing_knot_results.py`. The cleanup documented here is a
subsequent derivative; `DATABASE-CLEANUP.json` identifies both database hashes.

Run `python reproducibility/plot_submission_numerics.py` from the submission
folder, or invoke that script by its path from any directory. The script uses
the repaired databases, verifying that each derives from the expected file
(for knots, the cleaned database with the recorded cleanup checksum) and that
its repair record reports no incomplete row. To
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
repair remain undone. *Resolved 2026-11-22 except for the fourteen-crossing
database: see "Verification and repair" above.* The new rows store integer Betti values, bringing that
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
*Resolved 2026-11-22: exact integer ranks give beta_Q=1 for this block.*

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
