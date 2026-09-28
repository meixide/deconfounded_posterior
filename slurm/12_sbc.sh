#!/bin/bash
#-----------------------------------------------------------------------
# Simulation-based calibration for both parameterisations, with and without
# heterogeneous environment covariances.
#
# The grid is 2 models x 2 arms x 2 heterogeneity levels = 8 array tasks:
#
#   model       arm      hetero   what it answers
#   ---------------------------------------------------------------------
#   gi_hd       known    0.0      reference; this passed before (0.948)
#   gi_hd       plugin   0.0      reference; this failed before (0.784)
#   gi_hd       known    0.5      is the K-parameterisation itself sound when
#                                 the Sigma_e differ but are supplied?
#   gi_hd       plugin   0.5      how much worse does heterogeneity make the
#                                 plug-in problem?
#   gi_hd_slope known    0.0      implementation check, easy case
#   gi_hd_slope plugin   0.0      does the slope model lose anything to
#                                 plug-in when the Sigma_e are equal?
#   gi_hd_slope known    0.5      **the implementation check that matters**
#   gi_hd_slope plugin   0.5      **the headline: if this is uniform, the
#                                 reparameterisation is a fix, not a patch**
#
# Reading it: `known` failing means a bug in the Stan file, and nothing else
# in the table can be interpreted until it is fixed. `known` passing while
# `plugin` fails means the plug-in covariances are still doing damage.  Both
# passing for the slope model at hetero 0.5 is the result that closes §3.8.
#
# Cost: about 90s per fit with the two chains running in parallel on the two
# allocated cpus, so roughly 2.5-3h for 100 replications in one task.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-8 slurm/12_sbc.sh --reps=100 --iter=4000
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_sbc
#SBATCH --output=logs/sbc-%A_%a.out
#SBATCH --error=logs/sbc-%A_%a.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
# Home_FT2 is a shared volume that has run out of space under us before: on
# 2026-07-31 every task of an 8-task array died here within seconds with a bare
# "mkdir: write error", losing the whole array to a blip that cleared minutes
# later.  Retry a few times, and if it is still failing say so loudly rather
# than letting `set -e` kill the task with an unattributable one-liner.
for attempt in 1 2 3 4 5; do
    if mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/sbc" \
                "${BGI_ROOT}/results/summaries" 2>/dev/null &&
       touch "${BGI_ROOT}/results/sbc/.writable" 2>/dev/null; then
        rm -f "${BGI_ROOT}/results/sbc/.writable"
        break
    fi
    if [ "$attempt" = 5 ]; then
        echo "FATAL: ${BGI_ROOT}/results is not writable after 5 attempts." >&2
        echo "Shared-filesystem space or quota is the usual cause:" >&2
        df -h "${BGI_ROOT}" >&2 || true
        quota -s 2>/dev/null >&2 || true
        exit 1
    fi
    echo "WARN: results tree not writable (attempt ${attempt}); retrying." >&2
    sleep $(( attempt * 15 ))
done
# The two chains run in parallel, one per allocated cpu; keep BLAS single
# threaded so the allocation is not oversubscribed.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

MODELS=(gi_hd gi_hd gi_hd gi_hd gi_hd_slope gi_hd_slope gi_hd_slope gi_hd_slope)
ARMS=(known plugin known plugin known plugin known plugin)
HETERO=(0 0 0.5 0.5 0 0 0.5 0.5)

I=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
MODEL="${MODELS[$I]}"
ARM="${ARMS[$I]}"
H="${HETERO[$I]}"
TAG="${MODEL}_${ARM}_h${H}"

echo "Host : $(hostname)"
echo "Task : ${SLURM_ARRAY_TASK_ID:-<none>}  ${MODEL} / ${ARM} / hetero ${H}"
echo "Start: $(date)"

Rscript "${BGI_ROOT}/tests/test_sbc.R" \
    --model="${MODEL}" --arm="${ARM}" --hetero="${H}" --tag="${TAG}" "$@"

echo "End  : $(date)"
