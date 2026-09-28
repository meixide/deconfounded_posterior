#!/bin/bash
#-----------------------------------------------------------------------
# Signs of the causal coefficients in the BRFSS case study: one fit of the
# slope model on all 52 environments of the diabetic sample.  See
# scripts/18_case_signs.R for why a single fit suffices.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch slurm/18_case_signs.sh
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_signs
#SBATCH --output=logs/signs-%j.out
#SBATCH --error=logs/signs-%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
# A case-study fold took a median of about 20 minutes at 2000 iterations; this
# runs 4000 at adapt_delta 0.99, so allow the partition's 6 h ceiling.
#SBATCH --time=06:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Started : $(date)"

Rscript "${BGI_ROOT}/scripts/18_case_signs.R" \
        --data="${BGI_ROOT}/../data/brfss/brfss2023_case.csv" \
        --subset=diabetic --model=gi_hd_slope \
        --ncp=0 --eta-lkj=2 --adapt-delta=0.99 \
        --chains=4 --iter=4000 \
        --out="${BGI_ROOT}/results/case_signs"

echo "Finished: $(date)"
