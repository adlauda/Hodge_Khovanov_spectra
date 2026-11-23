#!/bin/bash
#SBATCH --job-name=kh-recheck
#SBATCH --account=lauda_1843
#SBATCH --partition=main
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-41
#SBATCH --output=/scratch1/%u/hks_review/logs/recheck_%A_%a.log

# Recompute every bidegree still flagged is_truncated=1 with exact ranks and
# rank-indexed gaps.  Tasks 0-39 are shards of the knot database; task 40 is
# the link database; task 41 is the twisted-unknot database.  Databases are
# opened read-only.  Each task appends to its own JSON-lines file and resumes
# from it if resubmitted.

module purge
module load usc

set -euo pipefail
REPO=/project2/lauda_1843/HKS_verify
OUT=/scratch1/$USER/hks_review/recompute
mkdir -p "$OUT"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH=$REPO/src
cd "$REPO"

T=$SLURM_ARRAY_TASK_ID
if [ "$T" -lt 40 ]; then
    python3 checks/database/recompute_rows.py --db databases/knot_research.db \
        --shard "$T" --shards 40 --done-glob "$OUT/knots_*.jsonl" \
        --out "$OUT/knots_$T.jsonl"
elif [ "$T" -eq 40 ]; then
    python3 checks/database/recompute_rows.py --db databases/link_research.db \
        --out "$OUT/links.jsonl"
else
    python3 checks/database/recompute_rows.py --db databases/twisted_unknot_research.db \
        --rank-max 12000 --out "$OUT/twisted_unknots.jsonl"
fi
