# Database verification and repair, November 2026

This folder holds the programs, batch files and results of the check of the
stored Betti numbers and spectral gaps carried out on 2026-11-22, and of the
comparison of the 14-crossing Betti numbers with khoca on 2026-11-23. The account
of what was found is in
[`NUMERICAL-REVIEW.md`](../../Paper/Current/reproducibility/NUMERICAL-REVIEW.md);
this page says what each file is and how to rerun it.

All programs open the databases read-only, except `apply_repairs.py`, which
writes a new copy and refuses to overwrite an existing file.

The repaired files now replace `knot_research.db`, `link_research.db` and
`twisted_unknot_research.db` in `databases/`. The files that were checked and
repaired are the ones at commit `9a06f279dc323ec97d996e21b98ee915020f32a6`,
called "the source files" below.

## Programs

| File | What it does | Needs |
|---|---|---|
| `compare_knotinfo.py` | Compares every stored knot Betti number with KnotInfo's integral Khovanov homology. | standard library, KnotInfo CSV |
| `compare_linkinfo.py` | Compares every stored link Betti table with LinkInfo's Khovanov polynomials, over orientations and mirrors. | standard library, LinkInfo CSV |
| `compare_khoca.py` | Compares every stored Betti number of a knot database with the rational Khovanov homology computed by the independent program khoca. Used for the 14-crossing databases, which KnotInfo does not cover. | khoca, py |
| `euler_check.py` | Checks that homology and chains have the same Euler characteristic for every diagram and quantum degree. | standard library |
| `exact_block_betti.py` | Exact rational Betti number of single bidegrees by integer rank. | python-flint |
| `recompute_rows.py` | Recomputes stored bidegrees: Betti number from exact modular ranks, gap read at the index the rank dictates, with a second independent gap computation. | NumPy, SciPy, `spectral_kh` |
| `apply_repairs.py` | Builds a repaired database copy from recomputed rows that pass every acceptance test, with a log of each old and new value. | standard library |
| `slim_database.py` | Writes a copy of a database without the `log_torsion` and `compute_time` columns, and optionally without named indexes, and verifies that nothing else changed. | standard library |
| `slurm/*.sbatch` | Slurm batch files for the recomputation, with site-specific options removed (8 CPUs and 32 GB per task). | |

The KnotInfo and LinkInfo tables are those shipped with the Python package
`database_knotinfo` (`knotinfo_data_complete.csv`, SHA-256 `ed4447ff…`, and
`linkinfo_data_complete.csv`, SHA-256 `308564aa…`; full digests are in the
comparison files).

The rank and Gram-matrix routines used by `recompute_rows.py` were written in
that script and afterwards moved unchanged to `src/spectral_kh/exact_rank.py`,
so that the database drivers use the same code.

## Results

On the source files:

| File | Content |
|---|---|
| `before_repair/knotinfo_comparison.json` | 12,155 of 12,965 knots agree; 1,774 rows on 810 knots disagree, all flagged `is_truncated`. Also the verdict on the 110 conflicting duplicate records. |
| `before_repair/linkinfo_comparison.json` | All 1,424 links agree. |
| `before_repair/euler_*.json` | Euler check. 804 knots fail in `knot_research.db`; no link or twisted unknot fails. |
| `euler_knot_research_14*.json` | Euler check of the two 14-crossing databases, which were not repaired. Only the partial records fail. |
| `khoca_knot_research_14*.json` | Comparison of the two 14-crossing databases with khoca 1.5, run on 2026-11-23. Every stored Betti number agrees: 275,896 rows (2,133 complete and 21 partial knots) and 277,468 rows (2,153 complete and 2 partial knots). |
| `exact_blocks_knots.json`, `exact_blocks_links.json` | Exact rational ranks for nine small frontier blocks, including `13a_4874` at (10, 31). |
| `recompute/knots_*.jsonl` | One line per recomputed row: 7,439 rows of 1,389 knots. The `knots_debug_*` files come from five extra workers that ran the same program with `--shards 5 --reverse` in a short-job queue; together with `knots_0` to `knots_39` they cover each knot once. |
| `recompute/links.jsonl` | 88 rows of 24 links. |
| `recompute/twisted_unknots.jsonl` | 141 rows of the diagrams `TU_10` to `TU_13`. |

On the repaired files now in `databases/`:

| File | Content |
|---|---|
| `databases/*.REPAIR.json` | Source and output checksums, counts of rows repaired and changed, and the post-conditions. |
| `databases/*.SLIM.json` | Checksums before and after the removal of the `log_torsion` and `compute_time` columns on 2026-11-23, for all five databases. For `knot_research.db` the index `unique_knot_bidegree` was removed as well. |
| `after_repair/knotinfo_comparison.json` | All 12,965 knots and 1,401,192 rows agree with KnotInfo. |
| `after_repair/linkinfo_comparison.json` | All 1,424 links agree with LinkInfo. |
| `after_repair/euler_*.json` | No diagram fails the Euler check. |

Each repaired database also carries the tables `row_repairs`,
`polynomial_repairs` and, for knots, `pd_code_fills`.

## Fields of a recomputed row

| Field | Meaning |
|---|---|
| `h`, `q` | Bidegree as labelled in the database. |
| `engine_h`, `engine_q` | The same bidegree in the current engine's labels (they differ only for twisted unknots). |
| `dimension`, `stored_dimension` | Chain-group dimension from the engine and from the database. |
| `d2_exactly_zero` | Whether the two differentials compose to zero in integer arithmetic. |
| `rank_out`, `rank_in`, `betti_mod_p` | Ranks modulo 65521 and the resulting Betti number; `null` for blocks too large for dense ranks. |
| `gap_from_out`, `gap_from_in`, `gap` | Smallest positive eigenvalue of each differential's Gram matrix at the rank-dictated index, and their minimum. |
| `zero_leak` | Largest eigenvalue at an index the rank says is zero; rounding error, to be compared with `gap`. |
| `laplacian_gap`, `laplacian_method` | The gap recomputed from the full Laplacian by a different route. |
| `stored_betti`, `stored_gap`, `stored_is_truncated` | Values in the database before the repair. |

## Rerunning

From the repository root, with the package installed (`pip install -e ".[checks]"`):

```bash
KI=$(python -c "import database_knotinfo, os; print(os.path.dirname(database_knotinfo.__file__))")/csv_data
OUT=$(mktemp -d)

# Independent comparisons of the current databases (minutes, standard library only)
python checks/database/compare_knotinfo.py --db databases/knot_research.db \
    --knotinfo $KI/knotinfo_data_complete.csv --out $OUT/knotinfo.json
python checks/database/compare_linkinfo.py --db databases/link_research.db \
    --linkinfo $KI/linkinfo_data_complete.csv --out $OUT/linkinfo.json
python checks/database/euler_check.py --db databases/knot_research.db --out $OUT/euler.json

# Rank routine against FLINT, and the engine regression tests
python checks/database/recompute_rows.py --self-test
python tests/test_exact_rank.py
```

To repeat the repair itself, first extract the source files:

```bash
SRC=9a06f279dc323ec97d996e21b98ee915020f32a6
mkdir $OUT/source
for n in knot_research link_research twisted_unknot_research; do
    git show $SRC:databases/$n.db > $OUT/source/$n.db
done

# Recompute the truncated rows (about 400 core-hours in total; see slurm/)
python checks/database/recompute_rows.py --db $OUT/source/knot_research.db \
    --shard 0 --shards 40 --out $OUT/knots_0.jsonl

# The same as a batch job, for all three databases
DB=$OUT/source sbatch --export=ALL checks/database/slurm/recompute_truncated.sbatch
DB=$OUT/source sbatch --export=ALL checks/database/slurm/recompute_twisted.sbatch

# Or build the repaired copy from the recorded results
python checks/database/apply_repairs.py --source $OUT/source/knot_research.db \
    --output $OUT/knot_research.db --report $OUT/knot_research.REPAIR.json \
    --recomputed "checks/database/recompute/knots_*.jsonl" \
    --knotinfo $KI/knotinfo_data_complete.csv --pd-csv data/knotinfo.csv
```

`apply_repairs.py` is deterministic apart from the timestamp and paths in its
report, so the last command reproduces the output SHA-256 recorded in
`databases/knot_research.REPAIR.json`. Removing the two columns then
reproduces the published file:

```bash
python checks/database/slim_database.py --source $OUT/knot_research.db \
    --output $OUT/knot_research.slim.db --report $OUT/knot_research.SLIM.json \
    --drop-index unique_knot_bidegree
sha256sum $OUT/knot_research.slim.db databases/knot_research.db
```
 Pass
`--cleanup Paper/Current/reproducibility/DATABASE-CLEANUP.json` to
`compare_knotinfo.py` on the source file to repeat the verdict on the 110
conflicting duplicate records.
