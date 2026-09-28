#!/bin/bash
#-----------------------------------------------------------------------
# Environment-budget study: how many environments does reliable causal
# discovery need, beyond the formal minimum E = p + 1?
#
# One array task = one (E, n_e) cell, all replications of it.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-8 slurm/04_env_budget.sh --n-rep=20
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_envbudget
#SBATCH --output=logs/envbudget-%A_%a.out
#SBATCH --error=logs/envbudget-%A_%a.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=20:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/env_budget"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Host  : $(hostname)"
echo "Cell  : ${SLURM_ARRAY_TASK_ID:-<none>}"
echo "Start : $(date)"

Rscript "${BGI_ROOT}/scripts/04_sim_environment_budget.R" \
    --chains="${SLURM_CPUS_PER_TASK:-4}" "$@"

echo "End   : $(date)"
