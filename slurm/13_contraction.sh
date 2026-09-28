#!/bin/bash
#-----------------------------------------------------------------------
# Tier 1: posterior contraction in N (Theorems 1 and 2).
#
# One array task = one replication, running all seven sample sizes.  The
# sizes are grouped this way rather than given a task each because the
# design (w*, mu_e, Sigma_e, and hence c_in) is held fixed within a
# replication, so the seven fits of one replication belong together: the
# comparison the study makes is across N *within* a design, and a fit takes
# only a few seconds.
#
# 13_sim_contraction.R lays its grid out as expand.grid(rep_id, n_e), so the
# task ids belonging to replication r are r, r+R, r+2R, ... for R = n_rep.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-30 slurm/13_contraction.sh --n-rep=30
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_contr
#SBATCH --output=logs/contr-%A_%a.out
#SBATCH --error=logs/contr-%A_%a.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
# Home_FT2 is a shared volume that has run out of space under us before; see
# the note in 12_sbc.sh.  Retry rather than lose the array to a transient blip.
for attempt in 1 2 3 4 5; do
    if mkdir -p "${BGI_ROOT}/logs" "${BGI_ROOT}/results/contraction" 2>/dev/null &&
       touch "${BGI_ROOT}/results/contraction/.writable" 2>/dev/null; then
        rm -f "${BGI_ROOT}/results/contraction/.writable"
        break
    fi
    if [ "$attempt" = 5 ]; then
        echo "FATAL: ${BGI_ROOT}/results is not writable after 5 attempts." >&2
        df -h "${BGI_ROOT}" >&2 || true
        quota -s 2>/dev/null >&2 || true
        exit 1
    fi
    echo "WARN: results tree not writable (attempt ${attempt}); retrying." >&2
    sleep $(( attempt * 15 ))
done

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

N_SIZES=7
REP="${SLURM_ARRAY_TASK_ID:-1}"

# Pull --n-rep out of the passthrough arguments so the task ids can be built;
# it has to agree with what the R script uses to lay out the grid.
N_REP=30
for a in "$@"; do
    case "$a" in
        --n-rep=*) N_REP="${a#--n-rep=}" ;;
    esac
done

TASKS="${REP}"
for k in $(seq 1 $(( N_SIZES - 1 ))); do
    TASKS="${TASKS},$(( k * N_REP + REP ))"
done

echo "Host : $(hostname)"
echo "Rep  : ${REP} of ${N_REP}  -> tasks ${TASKS}"
echo "Start: $(date)"

Rscript "${BGI_ROOT}/scripts/13_sim_contraction.R" --task="${TASKS}" "$@"

echo "End  : $(date)"
