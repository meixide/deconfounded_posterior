#!/bin/bash
#-----------------------------------------------------------------------
# Per-parameter convergence diagnosis for the slope model on the case-study
# data.  See scripts/13_slope_convergence_diag.R for what it answers.
#
#     sbatch slurm/13_slope_diag.sh quiron "Balears, Illes"
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_sdiag
#SBATCH --output=logs/sdiag-%j.out
#SBATCH --error=logs/sdiag-%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=05:00:00

set -euo pipefail
module load cesga/system R/4.4.2
export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript "${BGI_ROOT}/scripts/13_slope_convergence_diag.R" \
        --dataset="${1:-quiron}" --target="${2:-Balears, Illes}" "${@:3}"
