#!/bin/bash
#-----------------------------------------------------------------------
# Predictive coverage as the dimension grows (Table 1), as a SLURM job array.
#
# One array task = one (cell, replication) pair = one BGI fit plus the least
# squares baseline.  Tasks are independent and each writes its own CSV, so a
# failed or pre-empted task can be resubmitted on its own.
#
# The task table is (3 values of p) x (4 values of n) x (replications), so
# the default --n-rep=24 gives 288 rows: 12 * n_rep.
#
# A QoS normally caps how many jobs a user may have SUBMITTED, not merely
# running, and `--array=1-288%48` counts as 288 against that cap -- the %48
# throttles concurrency only.  Pass --chunk=N to put N consecutive rows in one
# array index, which divides the array size by N without changing the work:
#
#     rows   chunk   --array        one index is
#      288      1    1-288          1 replication
#      288      6    1-48           6 replications of one cell
#      288     24    1-12           one whole cell
#
# Rows are ordered cell-major, so a chunk never straddles two cells and its
# runtime is N times one replication of that cell.
#
# The p = 10, n = 2000 cell is by far the heaviest (E = 11 environments of
# 2000 observations each, so N = 22000); the p = 2 cells finish in seconds.
# The walltime below must cover N replications of that cell, so RAISE IT when
# raising --chunk.  Measure one first:
#     Rscript scripts/19_sim_dimension_sweep.R --task=265
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-48 slurm/19_dimension_sweep.sh --chunk=6
# then
#     Rscript scripts/20_aggregate_dimension_sweep.R
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_dimsweep
#SBATCH --output=logs/dimsweep-%A_%a.out
#SBATCH --error=logs/dimsweep-%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
# Sized for --chunk=6, i.e. six replications of the p = 10, n = 2000 cell, on
# a cluster core rather than a laptop one.  slurm/06_case_study.sh already
# uses six hours on this partition, so the ceiling is known to be accepted.
#SBATCH --time=06:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/dimension_sweep"

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
# --n-rep, --chains or --cov-method are forwarded verbatim.
Rscript "${BGI_ROOT}/scripts/19_sim_dimension_sweep.R" \
    --chains="${SLURM_CPUS_PER_TASK:-4}" "$@"

echo "Finished   : $(date)"
