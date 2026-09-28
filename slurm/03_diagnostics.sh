#!/bin/bash
#-----------------------------------------------------------------------
# Validation and diagnostic runs.  These are long enough that they belong on
# a compute node rather than a login node.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch slurm/03_diagnostics.sh                 # all of them
#     sbatch slurm/03_diagnostics.sh coverage        # just one
#
# Valid selectors: recovery, fastref, fullcov, coverage
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_diag
#SBATCH --output=logs/diag-%j.out
#SBATCH --error=logs/diag-%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

WHICH="${1:-all}"

run_one () {
  echo
  echo "=============================================================="
  echo "$1  ($(date))"
  echo "=============================================================="
  shift
  Rscript "$@"
}

case "$WHICH" in
  recovery|all)
    run_one "Parameter recovery" \
      "${BGI_ROOT}/tests/test_recovery.R"
    ;;&
  fastref|all)
    run_one "Fast model vs reference implementation" \
      "${BGI_ROOT}/tests/test_fast_vs_reference.R"
    ;;&
  fullcov|all)
    run_one "Plug-in vs inferred covariances" \
      "${BGI_ROOT}/tests/test_plugin_vs_fullcov.R" --reps=20
    ;;&
  coverage|all)
    run_one "Sources of gamma under-coverage" \
      "${BGI_ROOT}/tests/test_gamma_coverage_sources.R" --reps=20
    ;;&
  *)
    ;;
esac

echo
echo "Finished: $(date)"
