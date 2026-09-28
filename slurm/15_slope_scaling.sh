#!/bin/bash
#-----------------------------------------------------------------------
# Why does the slope model fail on real data?  Varies E and n_e on one
# dataset with everything else held fixed.  See scripts/15_slope_scaling_diag.R
#
#     sbatch --array=1-12 slurm/15_slope_scaling.sh brfss
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_sscale
#SBATCH --output=logs/sscale-%A_%a.out
#SBATCH --error=logs/sscale-%A_%a.err
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

Rscript "${BGI_ROOT}/scripts/15_slope_scaling_diag.R" \
        --task="${SLURM_ARRAY_TASK_ID:-1}" --dataset="${1:-brfss}" "${@:2}"
