#!/bin/bash
#SBATCH --job-name=kh-recheck-tu
#SBATCH --account=lauda_1075
#SBATCH --partition=main
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch1/%u/hks_review/logs/recheck_tu_%j.log

# Recompute the truncated twisted-unknot rows.  This database labels
# bidegrees with a different normalisation from the current engine, so
# --auto-shift matches labels by the unique constant offset that reproduces
# every stored chain dimension.  Blocks above 12000 use the stored Betti
# number and two matrix-free Lanczos runs instead of dense exact ranks.

module purge
module load usc

set -euo pipefail
REPO=/project2/lauda_1843/HKS_verify
OUT=/scratch1/$USER/hks_review/recompute
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH=$REPO/src
cd "$REPO"

python3 checks/database/recompute_rows.py --db databases/twisted_unknot_research.db \
    --auto-shift --rank-max 12000 --out "$OUT/twisted_unknots.jsonl"
