#!/bin/bash
#-----------------------------------------------------------------------
# Compile the Stan models once, before any simulation is submitted.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch slurm/00_compile.sh
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_compile
#SBATCH --output=logs/compile-%j.out
#SBATCH --error=logs/compile-%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:30:00

set -euo pipefail

module load cesga/system R/4.4.2

# SLURM starts the job in the submission directory; make the project root
# explicit so the scripts do not depend on where they were launched from.
export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"

echo "Host       : $(hostname)"
echo "Project    : ${BGI_ROOT}"
echo "Started    : $(date)"

Rscript "${BGI_ROOT}/scripts/00_compile_models.R" "$@"

echo "Finished   : $(date)"
