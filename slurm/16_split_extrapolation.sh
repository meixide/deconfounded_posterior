#!/bin/bash
#-----------------------------------------------------------------------
# Train on half the environments, predict into the unseen half: far targets
# (genuine extrapolation) against near ones.  See
# scripts/16_split_extrapolation.R for why LOEO cannot answer this.
#
#     sbatch --array=1-60%17 slurm/16_split_extrapolation.sh brfss gi_hd
#
# Tasks enumerate (replicate, role, rank) as reps x 2 x k, so the array range
# must equal --reps * 2 * --k (default 10 x 2 x 3 = 60).
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_split
#SBATCH --output=logs/split-%A_%a.out
#SBATCH --error=logs/split-%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=05:30:00

set -euo pipefail
module load cesga/system R/4.4.2
export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript "${BGI_ROOT}/scripts/16_split_extrapolation.R" \
        --task="${SLURM_ARRAY_TASK_ID:-1}" \
        --dataset="${1:-brfss}" --model="${2:-gi_hd_slope}" "${@:3}"
