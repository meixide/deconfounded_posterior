#!/bin/bash
#-----------------------------------------------------------------------
# Convergence gate for the case study.
#
# HANDOFF_REAL_DATA.md is explicit that no number from this model goes into
# the paper until Rhat has been read at full chain length.  The 200-iteration
# timing test established per-iteration cost and nothing else: 100 warmup
# iterations is far too few for Stan's mass-matrix adaptation on a model with
# ~700 parameters, so its Rhat of 2.77 is uninterpretable rather than alarming.
#
# Run this before 06_case_study.sh, and read the Rhat column.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch slurm/05_case_convergence.sh brfss
#     sbatch slurm/05_case_convergence.sh quiron
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_conv
#SBATCH --output=logs/conv-%j.out
#SBATCH --error=logs/conv-%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
# rstan parallelises across chains, not within; leaving BLAS threaded as well
# oversubscribes the allocation and slows everything down.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

DATASET="${1:-brfss}"
case "$DATASET" in
  brfss)  DATA="${BGI_ROOT}/../data/brfss/brfss2023_case.csv" ;;
  quiron) DATA="${BGI_ROOT}/../old_code/quiron/quiron_final.csv" ;;
  *) echo "Unknown dataset: $DATASET (expected brfss or quiron)" >&2; exit 2 ;;
esac

echo "Dataset : $DATASET"
echo "Data    : $DATA"
echo "Started : $(date)"

# Positional: $1 dataset, $2 iterations. Anything after those is forwarded to
# the R script, so --model=, --target= and friends can be set from sbatch.
Rscript "${BGI_ROOT}/tests/test_case_study_convergence.R" \
        --data="${DATA}" --dataset="${DATASET}" \
        --chains=4 --iter="${2:-2000}" "${@:3}"

echo "Finished: $(date)"
