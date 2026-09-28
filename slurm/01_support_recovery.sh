#!/bin/bash
#-----------------------------------------------------------------------
# Support-recovery / false-discovery simulation, as a SLURM job array.
#
# One array task = one (scenario, replication) pair = one BGI fit plus the
# baselines.  Tasks are independent and each writes its own CSV, so a failed
# or pre-empted task can be resubmitted on its own and completed tasks are
# skipped automatically.
#
# The array size must match the task table, which is
#     (number of scenarios) x (replications per scenario).
# Query it rather than hard-coding it:
#     Rscript -e 'source("R/setup.R"); bgi_setup("."); nrow(sim_task_table(20))'
# or use slurm/submit_support_recovery.sh, which does this for you.
#
# CESGA FinisTerrae III.  Submit from the new_code directory, e.g.
#     sbatch --array=1-80%40 slurm/01_support_recovery.sh --n-rep=20 --small
#
# Compile first, once, or every task in the array stops in its first seconds:
#
#     Rscript scripts/00_compile_models.R --force
#
# A task checks the compiled cache and refuses a stale one rather than
# rebuilding it, because eighty tasks writing the same cache would race.
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_support
#SBATCH --output=logs/support-%A_%a.out
#SBATCH --error=logs/support-%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/support_recovery"

# One MCMC chain per allocated core, and no hidden BLAS threading on top of
# it: nested parallelism oversubscribes the node and slows everything down.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Host       : $(hostname)"
echo "Array task : ${SLURM_ARRAY_TASK_ID:-<none>}"
echo "Cores      : ${SLURM_CPUS_PER_TASK:-1}"
echo "Started    : $(date)"

# SLURM_ARRAY_TASK_ID is read directly by the R script; extra flags such as
# --n-rep, --small, --chains or --cov-method are forwarded verbatim.
Rscript "${BGI_ROOT}/scripts/01_sim_support_recovery.R" \
    --chains="${SLURM_CPUS_PER_TASK:-4}" "$@"

echo "Finished   : $(date)"
