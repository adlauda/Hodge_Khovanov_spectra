#!/bin/bash
#SBATCH --job-name=kh-recheck
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-41
#SBATCH --output=recheck_%A_%a.log

# Recompute every bidegree flagged is_truncated=1 with exact ranks and
# rank-indexed gaps.  Tasks 0-39 are shards of the knot database; task 40 is
# the link database; task 41 is the twisted-unknot database (superseded by
# recompute_twisted.sbatch, which matches that database's bidegree labels).
# Databases are opened read-only.  Each task appends to its own JSON-lines
# file and resumes from it if resubmitted.
#
# Submit from the repository root.  DB must point to a folder holding the
# databases as they were before the repair (see ../README.md); the published
# files have no flagged rows left.  Add your scheduler's account and
# partition options, and load a Python with NumPy and SciPy.

set -euo pipefail
REPO=${REPO:-$SLURM_SUBMIT_DIR}
DB=${DB:?set DB to the folder with the pre-repair databases}
OUT=${OUT:-$REPO/recompute}
mkdir -p "$OUT"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH=$REPO/src
cd "$REPO"

T=$SLURM_ARRAY_TASK_ID
if [ "$T" -lt 40 ]; then
    python3 checks/database/recompute_rows.py --db "$DB/knot_research.db" \
        --shard "$T" --shards 40 --done-glob "$OUT/knots_*.jsonl" \
        --out "$OUT/knots_$T.jsonl"
elif [ "$T" -eq 40 ]; then
    python3 checks/database/recompute_rows.py --db "$DB/link_research.db" \
        --out "$OUT/links.jsonl"
else
    python3 checks/database/recompute_rows.py --db "$DB/twisted_unknot_research.db" \
        --rank-max 12000 --out "$OUT/twisted_unknots.jsonl"
fi
