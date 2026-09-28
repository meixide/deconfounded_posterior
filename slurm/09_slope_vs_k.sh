#!/bin/bash
#-----------------------------------------------------------------------
# Slope versus covariance parameterisation, across degrees of environment
# covariance heterogeneity.  One array task = one (heterogeneity, replication)
# pair; each task fits both models to the same simulated data, so the
# comparison is paired.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-90 slurm/09_slope_vs_k.sh --n-rep=30
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_slopek
#SBATCH --output=logs/slopek-%A_%a.out
#SBATCH --error=logs/slopek-%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=03:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/slope_vs_k"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Host : $(hostname)"
echo "Task : ${SLURM_ARRAY_TASK_ID:-<none>}"
echo "Start: $(date)"

Rscript "${BGI_ROOT}/scripts/09_sim_slope_vs_k.R" \
    --chains="${SLURM_CPUS_PER_TASK:-4}" "$@"

echo "End  : $(date)"
