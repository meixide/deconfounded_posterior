#!/bin/bash
#-----------------------------------------------------------------------
# Convenience launcher: works out the array size from the task table, submits
# the array, and chains the aggregation job so it runs once the array is done.
#
# Usage (from the new_code directory):
#     bash slurm/submit_support_recovery.sh [--n-rep=20] [--small]
#                                           [--max-concurrent=40]
#
# Everything after the recognised flags is forwarded to the R driver.
#-----------------------------------------------------------------------

set -euo pipefail

BGI_ROOT="${BGI_ROOT:-$PWD}"
export BGI_ROOT

N_REP=20
SMALL="FALSE"
MAX_CONCURRENT=40
PASSTHROUGH=()

for arg in "$@"; do
  case "$arg" in
    --n-rep=*)          N_REP="${arg#*=}"; PASSTHROUGH+=("$arg") ;;
    --small)            SMALL="TRUE";      PASSTHROUGH+=("$arg") ;;
    --max-concurrent=*) MAX_CONCURRENT="${arg#*=}" ;;
    *)                  PASSTHROUGH+=("$arg") ;;
  esac
done

if [ ! -f "${BGI_ROOT}/results/compiled/gi_hd.rds" ]; then
  echo "ERROR: models are not compiled." >&2
  echo "Run 'sbatch slurm/00_compile.sh' and wait for it to finish." >&2
  exit 1
fi

module load cesga/system R/4.4.2

# Ask the task table how many tasks there are, so the array bound can never
# drift out of step with the scenario grid.
N_TASKS=$(Rscript -e "
  suppressMessages({
    source('${BGI_ROOT}/R/setup.R')
    bgi_setup('${BGI_ROOT}', quiet = TRUE)
  })
  cat(nrow(sim_task_table(n_rep = ${N_REP}, small = ${SMALL})))
")

echo "Submitting ${N_TASKS} tasks (${N_REP} replications per scenario)."

ARRAY_JOB=$(sbatch --parsable \
  --array="1-${N_TASKS}%${MAX_CONCURRENT}" \
  "${BGI_ROOT}/slurm/01_support_recovery.sh" "${PASSTHROUGH[@]}")
echo "Array job      : ${ARRAY_JOB}"

AGG_JOB=$(sbatch --parsable \
  --dependency="afterany:${ARRAY_JOB}" \
  "${BGI_ROOT}/slurm/02_aggregate.sh")
echo "Aggregation job: ${AGG_JOB} (runs after the array finishes)"

echo
echo "Track progress with:  squeue -u \$USER"
echo "Results appear in  :  ${BGI_ROOT}/results/support_recovery/"
echo "Summaries in       :  ${BGI_ROOT}/results/summaries/"
