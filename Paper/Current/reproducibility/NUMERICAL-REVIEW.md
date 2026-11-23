# Numerical data: provenance, verification and limits

Last updated 2026-11-23.

This document describes the numerical data behind the manuscript: where each
database comes from, what has been checked and how, what was corrected, and
what the data do not establish. Commands are run from the repository root
unless stated otherwise.

## Summary

- The manuscript's figures and tables use one number per diagram, the smallest
  stored positive Laplacian eigenvalue. They are regenerated from the databases
  in `databases/` by `plot_submission_numerics.py`, which checks every input
  file against recorded checksums.
- Every Betti number in the knot (through 13 crossings), link and
  twisted-unknot databases agrees with an independent reference: KnotInfo,
  LinkInfo, or the homology of the unknot.
- In November 2026, 1,774 Betti numbers on 810 twelve-crossing knots were
  found to be too small and were corrected. No spectral gap changed, and the
  figures regenerated afterwards are pixel-identical to the earlier ones.
- The spectral gaps are double-precision values without validated error
  bounds. The fourteen-crossing data have no independent reference.

## Provenance of the databases

The manuscript cites the repository at commit
`ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea`. The files now in `databases/`
derive from the files at that commit in at most three recorded steps. Each
step has a record with the checksum of its input and output, and the plotting
script follows the chain.

| File | At `ccfce9a` | After cleanup | After repair | Published |
|---|---|---|---|---|
| `knot_research.db` | `918832fc…` | `f06d7982…` | `37c52286…` | `b8362ff5…` |
| `link_research.db` | `6c0f6ee5…` | | `acf404fd…` | `c8255715…` |
| `twisted_unknot_research.db` | `72656193…` | | `82d58cc6…` | `5536e339…` |
| `knot_research_14.db` | `2c3e9d55…` | | | `191fc2be…` |

The entries are the first eight digits of SHA-256 checksums; the records hold
them in full.

| Step | Date | Program | Record |
|---|---|---|---|
| Cleanup: remove duplicate rows | 2026-09 | `clean_knot_database.py` | `DATABASE-CLEANUP.json` |
| Repair: correct Betti numbers, clear flags | 2026-11-22 | `checks/database/apply_repairs.py` | `databases/*.REPAIR.json` |
| Column removal: reduce file size | 2026-11-23 | `checks/database/slim_database.py` | `databases/*.SLIM.json` |

Intermediate files are in the Git history: after cleanup at commit
`9a06f279dc323ec97d996e21b98ee915020f32a6`, and after repair at commit
`5d61309a76522d5e8c59197de1d6e72fe861979c`.

`databases/knot_research_14_20261007.db` is a later snapshot of the
fourteen-crossing computation. The manuscript does not use it.

## Coverage

| Database | Diagrams | Complete | Bidegree rows |
|---|---|---|---|
| Knots through 13 crossings | 12,965 | 12,965 | 1,401,192 |
| Fourteen-crossing knots | 2,154 | 2,133 | 275,896 |
| Links | 1,424 | 1,424 | 110,646 |
| Twisted unknots, n = 1, …, 13 | 13 | 13 | 663 |

- A diagram is "complete" when its Khovanov polynomial is stored. The figures
  use complete diagrams only.
- The fourteen-crossing knots were drawn by `scripts/kh_database_14.py` from
  `data/knots14_pd.csv`, a table of 46,969 planar diagrams (19,536 `14ah` and
  27,433 `14nh`), using `DataFrame.sample` with a fixed seed over the entries
  not yet finished. Every PD code in the database equals its entry in that
  table. This is random selection from the supplied table, which is not the
  same as uniform sampling over knot types, and the completed subset need not
  keep the inclusion probabilities of the draws. The script performs no
  hyperbolicity test.
- A null gap does not in general indicate a solver failure: some bidegrees
  have no positive spectrum. The databases do not record which driver and
  threshold produced each row.
- With the 21 partial fourteen-crossing records removed, the overall
  fourteen-crossing minimum is `0.0006119790255380953`. The minimum over the
  886 complete alternating records is `0.08654902144568334`, at `14ah_00033`.
- The 13-crossing alternating minimum is `0.05811636514789631`, at
  `T(2,13) = 13a_4878`, agreeing numerically with its distinguished-bidegree
  formula.

## Cleanup of duplicate rows (September 2026)

The knot database at `ccfce9a` has 1,453,560 bidegree rows for 12,965 knots.
It contains 47,634 repeated `(knot_id, h, q)` groups, amounting to 52,368
extra rows on 552 knots, left by overlapping runs of different drivers. The
other databases have no repeated groups.

- All repeated rows agree on the chain dimension, and their positive gaps
  differ by at most `7.925841538636291e-9`.
- The Betti numbers disagree in 110 repeated bidegrees on 50 knots. For
  example, `12a_3` at (-1, -1) has stored values 4 and 11.
- For every bidegree, the row with the highest ID agrees with the knot's
  stored Khovanov polynomial. The cleanup keeps that row. The November
  comparison with KnotInfo confirmed this choice in all 110 conflicts.
- The cleaned file has 1,401,192 rows. Every kept row is unchanged.
  Per-diagram minima change by at most `1.4337531162311734e-10`, which does
  not alter any gap quoted in the manuscript.
- `DATABASE-CLEANUP.json` lists every removed row ID with the row kept in its
  place, and the conflicting Betti records.

To repeat the cleanup, extract the file at `ccfce9a` and run
`clean_knot_database.py` with `--source`, `--output` and `--report`. The
script refuses to overwrite existing files.

## Verification and repair (November 2026)

Programs, batch files and result files are in `checks/database/`, indexed in
its `README.md`. "Before the repair" below refers to the files at commit
`9a06f279dc323ec97d996e21b98ee915020f32a6`.

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
  the 1,774 that did not were all truncated rows.
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

The last three columns are the largest relative values over the rows. "Zero
leak" is the largest computed eigenvalue at an index the exact rank says is
zero. In 113 twisted-unknot rows a differential has more than 12,000 rows or
columns and dense exact ranks were not attempted. Those rows use the stored
Betti number, which the unknot argument above fixes, and two Lanczos runs on
the Laplacian from independent start vectors.

A first pass over the twisted unknots was invalid: it looked up blocks by the
database's labels in the engine's normalisation and so examined the wrong
blocks. The chain-dimension test rejected all 141 rows. The pass was repeated
with the labels matched by the unique constant offset that reproduces every
stored chain dimension, and the invalid output was discarded.

### Repaired databases

A recomputed row was accepted only if its chain dimension matched, `d^2 = 0`
held, its Betti number equalled the reference, the two gap computations
agreed to `1e-6` and the zero leak was below `1e-6` of the gap. Every row
passed.

| File | Rows repaired | Betti changed | Polynomials rebuilt | Blobs made integers |
|---|---|---|---|---|
| `knot_research.db` | 7,439 | 1,774 | 810 | 1,183,182 |
| `link_research.db` | 88 | 0 | 0 | 88 |
| `twisted_unknot_research.db` | 141 | 0 | 0 | 141 |

- No row carries `is_truncated=1`. The flag now means that a row's Betti
  number and gap have not been verified; `row_repairs.old_is_truncated` keeps
  the old flag.
- Every `betti` value has integer storage. Before the repair most were stored
  as eight-byte little-endian blobs, on which SQL comparisons and sums do not
  work.
- The 7,930 knots that stored no PD code now have the code from
  `data/knotinfo.csv`; the 5,035 codes already stored equal that table's.
- Tables `row_repairs`, `polynomial_repairs` and `pd_code_fills` in each file
  list every old and new value. No row changed outside that log.
- After the repair all 12,965 knots agree with KnotInfo in all 1,401,192 rows,
  all 1,424 links agree with LinkInfo, and no diagram fails the Euler check
  (`checks/database/after_repair/`, run on the published files).
- Rebuilding the three files from `checks/database/recompute/` reproduces the
  recorded checksums, and removing the columns and index from the rebuilt
  files reproduces the checksums of the published ones.

### Engine changes

- `src/spectral_kh/exact_rank.py` holds the modular rank, the rank-indexed gap
  and the self-test.
- `KhovanovEngine.get_low_spectrum_adaptive`, which every database driver
  calls, now takes its Betti number from exact ranks and reads the gap at that
  index. Its return format is unchanged. A bidegree whose differentials exceed
  25,000 rows or columns falls back to the previous threshold method with a
  `RuntimeWarning`.
- `KhovanovEngine.compute_khovanov_polynomial` called two functions that did
  not exist in the module; it now works.
- `tests/test_exact_rank.py` checks the full Betti tables and minimum gaps of
  `3_1`, `4_1`, `8_19` and `9_42`, and optionally a 4,000-dimensional block of
  `12n_803` that the truncated solver recorded as 5 instead of 7.

## Column and index removal (2026-11-23)

To reduce the size of the published files,
`checks/database/slim_database.py` removed two columns of the `bidegrees`
table from every database in `databases/`: `compute_time`, and `log_torsion`
where present (the fourteen-crossing databases never had it). It also removed
the index `unique_knot_bidegree` from `knot_research.db`. No figure, table or
check in this repository reads either column, and the index holds no data of
its own. The script verified that every remaining column of every table is
identical before and after, and that no other index changed.

| File | Before | After |
|---|---|---|
| `knot_research.db` | 84.2 MB | 51.5 MB |
| `knot_research_14.db` | 21.5 MB | 18.0 MB |
| `knot_research_14_20261007.db` | 21.7 MB | 16.3 MB |
| `link_research.db` | 8.1 MB | 6.5 MB |
| `twisted_unknot_research.db` | 0.10 MB | 0.09 MB |

- The files with both columns and the index are at commit
  `5d61309a76522d5e8c59197de1d6e72fe861979c`.
- `knot_research.db` still has exactly one row per knot and bidegree, but
  nothing now prevents a duplicate from being inserted, and lookups by knot
  scan the table. Before adding rows or querying knot by knot, rebuild the
  index with
  `CREATE UNIQUE INDEX unique_knot_bidegree ON bidegrees(knot_id, h, q);`.
- The drivers in `scripts/` create and write both columns, so they need a
  new database or one from the commit above, not the published files.
- The September cleanup kept the newest row of each duplicated bidegree. For
  35,959 bidegrees on 443 knots that row had no `log_torsion` while a removed
  row did. Those values were therefore absent from the cleaned knot database
  already; they are in the file at `ccfce9a`.

## What the checks establish, and what they do not

For integer differentials `A = d_out` and `B = d_in` with `AB = 0`, the
rational Betti number is `dim - rank_Q(A) - rank_Q(B)`.

- **Modular ranks give an upper bound.** A rank can only drop modulo a prime,
  so `dim - rank_p(A) - rank_p(B)` is at least the rational Betti number. A
  modular value of zero therefore proves the rational Betti number is zero,
  and one prime suffices. A positive modular value proves only an upper
  bound; agreement between two primes does not make it exact.
- **Exact integer ranks give the Betti number.** For nine small blocks at the
  low end of the gap distribution, ranks were computed over the rationals with
  FLINT, with `AB = 0` checked exactly (`checks/database/exact_blocks_*.json`).
  Seven have Betti number 0. `9_46` at (0, -1) and `13a_4874` at (10, 31) have
  Betti number exactly 1; for the latter, modular ranks alone had shown only
  that it is at most 1.
- **The stored Betti numbers** of knots and links are confirmed by agreement
  between the modular upper bound, or the dense spectrum, and the independent
  KnotInfo and LinkInfo tables.
- **The gaps are not certified.** A small residual bounds the distance from
  an approximate eigenvalue to some true eigenvalue. It does not show that the
  eigenvalue is the smallest positive one: for `L = diag(ε, 1)`, the pair
  `v = (0, 1)`, `λ = 1` has zero residual although the smallest eigenvalue is
  `ε`. Floating-point residuals also lack outward rounding. Agreement between
  independent solvers is good numerical evidence and is all that is claimed.
  `frontier_diagnostics.json` records residuals for eleven blocks as diagnostics,
  not as enclosures.
- **One exception:** `checks/check_gap_9_46.py` bounds the first positive
  eigenvalue of `9_46` in bidegrees (-1, -1) and (0, -1) between 0.0289984 and
  0.0289985 in exact integer arithmetic.
- **Rows below dimension 4,000** were written from a dense full spectrum and
  were not recomputed. Their Betti numbers are covered by the KnotInfo and
  LinkInfo comparisons; their gaps rest on the dense eigensolver.
- **The fourteen-crossing databases** were not recomputed. KnotInfo has no
  Khovanov data at 14 crossings, so there is no independent reference. All
  complete records pass the Euler check, and only the partial records fail it:
  21 in `knot_research_14.db`, and 2 (`14ah_00002`, `14ah_00010`) in the
  2026-10-07 snapshot. Their Betti values are stored as blobs in
  `knot_research_14.db` and as integers in the snapshot.
- **The engine's `compute_jones_polynomial`** uses numerical ranks and is not
  an independent exact Jones calculation.

## Notes on the engine

Line numbers refer to `src/spectral_kh/khovanov_core.py` at commit `ccfce9a`.

- The grading shifts (lines 196-201) change the bidegree labels, not the
  boundary matrices or their positive eigenvalues.
- `compute_laplacian` (lines 504-522) returns the sum of the upward and
  downward Laplacians without rescaling by an operator norm.
- At that commit the adaptive solver classified eigenvalues below a tolerance
  as zero, with default `1e-9`; `scripts/kh_database_manager_PE.py` and
  `scripts/kh_database_14.py` passed `1e-7`. This is a classification
  threshold, distinct from the eigensolver's convergence tolerance. Since the
  November change the tolerance is used only in the large-block fallback.

## Scope of the exact companion formulas

[Grlj and Lauda, Spectral geometry of Khovanov Laplacians](https://arxiv.org/html/2608.24298v1)
proves bidegree statements in Theorems 4.26, 4.29 and 4.32. Corollary 4.28,
Remark 4.30 and the conclusion leave open whether those bidegrees attain the
minimum over the whole complex. The formulas are therefore labelled as
distinguished-bidegree gaps in the figures, even where the observed minima
agree with them.

`data/knotinfo.csv` has SHA-256
`f5d0a015177423c67e5fa2c428e2b6524ad04a4f9302c2c3f8ad2737d221340b`. An
orientation check of that table with Spherogram found all crossings positive
for the torus representatives `3_1`, `5_1`, `7_1`, `9_1`, `11a_367` and
`13a_4878`, consistent with the positive-braid formula. The even twist
representatives `4_1`, `6_1`, `8_1` and `10_1` have `n - 2` positive and 2
negative crossings, consistent with the companion's stated shift. This sign
check does not prove that every such diagram is identical to the companion's
standard representative. No family bidegree is assigned by name to a
fourteen-crossing diagram. The exact formulas refer to the companion's
oriented representatives; the empirical statistics minimise over the stored
bidegrees.

## Figures

| File under `../figures/` | Content |
|---|---|
| `observed_knots_gap_distribution.png` | One stored positive minimum per complete knot: medians, central 50% and 80% bands, observed minimum, with counts of omitted records. |
| `observed_links_gap_distribution.png` | The same summaries, each link weighted equally. |
| `twisted_unknot_bidegree_gap.png` | Exact gap in the specified unshifted unreduced bidegree compared with observed full-diagram minima. |
| `alternating_distinguished_bidegrees.png` | Complete alternating records compared with exact distinguished-bidegree formulas; diagrams absent from the data are marked as formulas. |

The figures fit no asymptotic curves and make no confidence-interval claim
for their quantile bands. The earlier figures in `outputs/plots/` are kept for
reference only: their titles say "Global Worst-Case", and their means pool
the repeated bidegree rows that the cleanup later removed.

## Files in this folder

| File | Content |
|---|---|
| `plot_submission_numerics.py` | Deterministic plotting and summary script. Opens the databases read-only after verifying their checksum chain. |
| `NUMERICS-MANIFEST.json` | Database sources, checksums, coverage and figure checksums, written by the plotting script. |
| `*_diagram_summary.csv` | One row per diagram: completion state, number of stored blocks, smallest positive value. |
| `knots_gap_summary.csv`, `links_gap_summary.csv`, `alternating_gap_summary.csv` | The statistics shown in the figures. |
| `clean_knot_database.py`, `DATABASE-CLEANUP.json` | The September cleanup and its record. |
| `COMPLETION-VERIFICATION.json` | Comparison of the file at `ccfce9a` with the preceding snapshot: 37 knots and 3,689 rows added, every earlier row unchanged. |
| `frontier_diagnostics.json` | Modular ranks and floating-point residuals for eleven blocks at the low end of the gap distribution. |
| `check_claims.py` | Small deterministic checks of statements in the manuscript. |

To regenerate the figures and tables:

```bash
python Paper/Current/reproducibility/plot_submission_numerics.py
```

This needs Python 3.11 or later with NumPy, pandas and Matplotlib. It runs no
eigensolver.
