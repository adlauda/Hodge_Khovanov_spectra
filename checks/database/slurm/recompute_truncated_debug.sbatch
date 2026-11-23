#!/bin/bash
#SBATCH --job-name=kh-recheck-dbg
#SBATCH --account=lauda_1843
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch1/%u/hks_review/logs/recheck_debug_%A_%a.log

# One-hour passes over the truncated knot rows on the debug partition.
# Worker i of 5 takes every fifth diagram, starting from the END of the name
# list, and skips any diagram already recorded by any worker (debug or main).
# Resubmit until nothing is left; the debug QOS allows 5 jobs and 48 CPUs.

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

python3 checks/database/recompute_rows.py --db databases/knot_research.db \
    --shard "$SLURM_ARRAY_TASK_ID" --shards 5 --reverse \
    --done-glob "$OUT/knots_*.jsonl" --out "$OUT/knots_debug_$SLURM_ARRAY_TASK_ID.jsonl"
