#!/bin/bash
#-----------------------------------------------------------------------
# Leave-one-environment-out case study, one held-out environment per array
# task.
#
# One task per fold rather than one job for all of them, for two reasons: the
# csic account shares an organisation-wide CPU pool that is often saturated, so
# small tasks start sooner than one large one; and a fold that fails its
# sampler diagnostics then fails alone instead of taking the run with it.
#
# `scripts/07_aggregate_case_study.R` combines the per-fold CSVs afterwards and
# is what produces the numbers for Section 3.2.
#
# CESGA FinisTerrae III.  Submit from the new_code directory:
#     sbatch --array=1-52 slurm/06_case_study.sh brfss              # slope model
#     sbatch --array=1-52 slurm/06_case_study.sh brfss gi_hd        # K-param
#     sbatch --array=1-51 slurm/06_case_study.sh quiron
#
# Check how many folds are eligible first:
#     Rscript scripts/05_case_study.R --data=... --summary-only
#-----------------------------------------------------------------------
#SBATCH --job-name=bgi_case
#SBATCH --output=logs/case-%A_%a.out
#SBATCH --error=logs/case-%A_%a.err
# `short` rather than `medium`: the medium QoS caps a user at 50 submitted
# jobs, and these arrays are 51-52 folds, so medium rejects them outright.
# short allows 100 jobs at up to 6 h, which is the binding constraint below.
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# The gate run peaked at 3.7 GB, so 4G x 4 is ample; asking for 32 GB on a
# saturated shared pool only delays scheduling.
#SBATCH --mem-per-cpu=4G
# A converged K-parameterisation fold took 255 s.  The slope model carries
# E x p = 624 extra slope deviations and HANDOFF §3.8 measured it up to 10x
# slower where the covariances genuinely differ, so the old 3 h would be
# cutting it fine.  6 h is the ceiling on the short partition, and a timed-out
# array task loses the whole fold.
#SBATCH --time=06:00:00

set -euo pipefail

module load cesga/system R/4.4.2

export BGI_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "${BGI_ROOT}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Positional arguments are the dataset, the model and the iteration count, in
# that order.  Anything starting with -- is not positional: it is forwarded to
# the R script untouched, so a run can be varied without editing this file.
# Separating by shape rather than by position is what keeps the two apart.
# Shifting past the dataset instead, as a first attempt did, put
# --v-prior-shape into MODEL and --v-prior-rate into ITER, and the job died
# loading a Stan model by that name.  EXTRA further down is set by the dataset
# entry itself and is a different thing.
POSITIONAL=()
USER_ARGS=()
for arg in ${@+"$@"}; do
  case "$arg" in
    --*) USER_ARGS+=("$arg") ;;
    *)   POSITIONAL+=("$arg") ;;
  esac
done
DATASET="${POSITIONAL[0]:-brfss}"
# DS is the loader's dataset name; DATASET is the label for the results
# directory. They differ when a run restricts an existing dataset rather than
# introducing a new one, and passing both as --dataset= would collide: the R
# flag parser takes the first match.
DS="$DATASET"
EXTRA=""
# Set when a dataset entry supplies its own --adapt-delta via EXTRA.
ADAPT=""
case "$DATASET" in
  brfss)  DATA="${BGI_ROOT}/../data/brfss/brfss2023_case.csv"
          MIN_N0=300 ;;
  # Adults 60+ with diagnosed diabetes. A tenth the size of the full survey,
  # so the selection rule of Section 2.2 can discriminate instead of declaring
  # every covariate a parent. min_n0 is 200 rather than 300 because the
  # subgroup is smaller; that keeps 48 of 52 states as eligible targets.
  # Adults with diagnosed diabetes, no age cut. The case-study sample: a
  # seventh the size of the full survey, so the selection rule discriminates,
  # while N*lambda_min = 33 matches the full survey's identification strength.
  # adapt_delta 0.99 because the smaller sample flattens the gamma/b_bar ridge
  # and 0.95 produced divergences on the age-restricted variant.
  brfss_diab) DATA="${BGI_ROOT}/../data/brfss/brfss2023_case.csv"
          MIN_N0=300; DS=brfss; EXTRA="--subset=diabetic"; ADAPT="--adapt-delta=0.99" ;;
  brfss_diab60) DATA="${BGI_ROOT}/../data/brfss/brfss2023_case.csv"
          MIN_N0=200; DS=brfss; EXTRA="--subset=older_diabetic" ;;
  quiron) DATA="${BGI_ROOT}/../old_code/quiron/quiron_final.csv"
          MIN_N0=300 ;;
  acs_tract) DATA="${BGI_ROOT}/../data/candidates/tract/acs_tract_2022.csv"
          MIN_N0=300 ;;
  acs_pums) DATA="${BGI_ROOT}/../data/candidates/acs/acs_slim.csv"
          MIN_N0=300 ;;
  # Communities and Crime has a median of 36 units per state, so the 300-row
  # floor used elsewhere would discard every fold. The whole dataset is 1,931
  # rows; 20 is the smallest target that still supports a Ledoit-Wolf Sigma_0
  # at p = 10.
  communities) DATA="${BGI_ROOT}/../data/candidates/cc/communities.data"
          MIN_N0=20 ;;
  *) echo "Unknown dataset: $DATASET (expected brfss or quiron)" >&2; exit 2 ;;
esac

TASK="${SLURM_ARRAY_TASK_ID:-1}"

# $3 overrides the iteration count. A longer run writes to its own results
# directory rather than overwriting the 2000-iteration one, so the two remain
# comparable -- the point of a longer run is to check whether the short one was
# adequate, which requires keeping both.
ITER="${POSITIONAL[2]:-2000}"
SUFFIX=""
[ "$ITER" != "2000" ] && SUFFIX="_iter${ITER}"

# Which Stan model. `gi_hd_slope` is the reparameterisation of HANDOFF §3.8 and
# is the default; pass `gi_hd` as $2 to reproduce the K-parameterisation. The
# results directory carries the model name so the two never overwrite one
# another — the whole point is to compare them on the same folds.
MODEL="${POSITIONAL[1]:-gi_hd_slope}"
OUT="${BGI_ROOT}/results/case_${DATASET}_${MODEL}${SUFFIX}/folds"

echo "Dataset   : $DATASET"
echo "Fold index: $TASK"
echo "Started   : $(date)"

# ncp = 0: the centred parameterisation.  slurm/05_case_convergence.sh measured
# it on this data at Rhat 1.0023 / ESS 2053 in 255 s, against Rhat 1.0224 /
# ESS 99 in 6963 s non-centred -- 27x slower *and* not mixed.  The simulations
# need ncp = 1; this data does not, and the default does not transfer.
#
# adapt_delta = 0.95 rather than the 0.9 default: the gate run threw a single
# divergence at 0.9 on an otherwise clean fold, and across 52 folds isolated
# divergences would either bias the summary or trip the exclusion rule.
Rscript "${BGI_ROOT}/scripts/05_case_study.R" \
        --data="${DATA}" --dataset="${DS}" ${EXTRA} \
        --fold-index="${TASK}" --min-n0="${MIN_N0}" \
        --model="${MODEL}" \
        --ncp=0 --eta-lkj=2 ${ADAPT:---adapt-delta=0.95} \
        --chains=4 --iter="${ITER}" --cores=4 \
        --max-target=2000 \
        ${USER_ARGS[@]+"${USER_ARGS[@]}"} \
        --out="${OUT}/fold_$(printf '%03d' "${TASK}")"

echo "Finished  : $(date)"
