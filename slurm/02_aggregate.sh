#!/bin/bash
#-----------------------------------------------------------------------
# Aggregate the per-task CSVs into the manuscript summary tables.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch slurm/02_aggregate.sh
# or let slurm/submit_support_recovery.sh chain it after the array.
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_aggregate
#SBATCH --output=logs/aggregate-%j.out
#SBATCH --error=logs/aggregate-%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:20:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"

Rscript "${BGI_ROOT}/scripts/02_aggregate_support_recovery.R" "$@"
