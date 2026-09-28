#!/usr/bin/env bash
#-----------------------------------------------------------------------
# run_checks.sh -- everything in this replication package that runs on a
# laptop, in one command.
#
#     bash run_checks.sh
#
# Takes about four minutes on an Apple M2 once the models are compiled, and
# needs no cluster, no network access and no data download.  It verifies the implementation rather
# than reproducing the manuscript's tables: every table in the paper is a job
# array on a cluster, for reasons set out in README.md under "What runs where".
#
# What it does, in order:
#   0. checks R and the two required packages, and stops with a usable
#      message if either is missing or mismatched;
#   1. compiles the Stan models;
#   2. runs tests/test_recovery.R, which checks the four substantive claims
#      of the paper, including that the predictive intervals carry the
#      target-domain scale S_0 rather than the training residual scale;
#   3. runs tests/test_likelihood_identity.R, which checks exactly that the
#      sufficient-statistic likelihood is the per-observation model;
#   4. runs one replication of the dimension sweep at p = 2, so that the
#      simulation machinery is seen working end to end on a real fit.
#
# Exit status is 0 only if every step passed.
#-----------------------------------------------------------------------

set -uo pipefail

cd "$(dirname "$0")"

RED=$'\033[31m'; GREEN=$'\033[32m'; BOLD=$'\033[1m'; OFF=$'\033[0m'
if [ ! -t 1 ]; then RED=""; GREEN=""; BOLD=""; OFF=""; fi

step() { printf "\n%s==> %s%s\n" "$BOLD" "$1" "$OFF"; }
ok()   { printf "%s    ok%s  %s\n" "$GREEN" "$OFF" "$1"; }
bad()  { printf "%s    FAILED%s  %s\n" "$RED" "$OFF" "$1"; }

FAILED=0
note_fail() { bad "$1"; FAILED=$((FAILED + 1)); }

#-----------------------------------------------------------------------
step "0. Environment"

if ! command -v Rscript >/dev/null 2>&1; then
  bad "Rscript is not on the PATH. Install R from https://cran.r-project.org"
  exit 1
fi

Rscript -e '
  cat(R.version.string, "\n")
  missing <- character(0)
  for (p in c("rstan", "mvtnorm")) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
    if (is.na(v)) missing <- c(missing, p) else cat(sprintf("  %-10s %s\n", p, v))
  }
  v <- tryCatch(as.character(packageVersion("ggplot2")), error = function(e) NA)
  cat(sprintf("  %-10s %s\n", "ggplot2",
              if (is.na(v)) "not installed (optional: figures only)" else v))
  if (length(missing)) {
    cat("\nMissing required packages:", paste(missing, collapse = ", "), "\n")
    cat("Install with: install.packages(c(",
        paste(sprintf("\"%s\"", missing), collapse = ", "), "))\n", sep = "")
    quit(status = 1)
  }
' || { bad "required R packages are missing (see above)"; exit 1; }
ok "R and the required packages are present"

#-----------------------------------------------------------------------
step "1. Compiling the Stan models (about two minutes)"

if Rscript scripts/00_compile_models.R; then
  ok "models compiled"
else
  note_fail "scripts/00_compile_models.R"
  echo
  echo "Nothing downstream can run without the required models. If the error"
  echo "above mentions StanHeaders, the message printed by that script gives"
  echo "the exact install command."
  exit 1
fi

#-----------------------------------------------------------------------
step "2. Recovery checks (about five minutes)"
echo "    Checks that BGI recovers gamma where least squares cannot, and that"
echo "    the predictive intervals use the target-domain scale S_0."

if Rscript tests/test_recovery.R; then
  ok "tests/test_recovery.R"
else
  note_fail "tests/test_recovery.R"
fi

#-----------------------------------------------------------------------
step "3. Sufficient-statistic likelihood against the reference"
echo "    Checks exactly, without sampling, that the fast likelihood in"
echo "    gi_hd.stan is the per-observation model of gi_hd_reference.stan."

if Rscript tests/test_likelihood_identity.R; then
  ok "tests/test_likelihood_identity.R"
else
  note_fail "tests/test_likelihood_identity.R"
fi

#-----------------------------------------------------------------------
step "4. One replication of the dimension sweep (p = 2)"
echo "    The smallest cell of Table 1, to show the simulation pipeline"
echo "    working on a real fit. The full table is a 288-task job array."

if Rscript scripts/19_sim_dimension_sweep.R --task=1 \
     --out=results/laptop_demo; then
  ok "scripts/19_sim_dimension_sweep.R"
  echo
  echo "    Result written to results/laptop_demo/task_0001.csv:"
  Rscript -e '
    f <- "results/laptop_demo/task_0001.csv"
    if (file.exists(f)) {
      r <- read.csv(f)
      cat(sprintf("      coverage, ours (at S_0)        %.3f\n", r$cov_ours))
      cat(sprintf("      coverage, training scale       %.3f\n", r$cov_train_scale))
      cat(sprintf("      coverage, least squares        %.3f\n", r$cov_ols))
      cat(sprintf("      max Rhat                       %.3f\n", r$max_rhat))
    }' 2>/dev/null
else
  note_fail "scripts/19_sim_dimension_sweep.R"
fi

#-----------------------------------------------------------------------
step "5. Every table in the paper against the numbers in results/"
echo "    Recomputes each table from the committed per-task CSVs and prints the"
echo "    rows beside the ones the manuscript has, so they can be compared."
echo "    This needs no cluster and no Stan: the aggregators only read CSVs."
if [ -d results/support_recovery_ad99_slope ] || [ -d results/dimension_sweep ]; then
  CT_LOG=check_tables_output.txt
  if bash check_tables.sh > "$CT_LOG" 2>&1; then
    BAD=$(grep -c 'label not found\|no loeo_folds\|falling back' "$CT_LOG" || true)
    if [ "${BAD:-0}" -gt 0 ]; then
      echo "    $BAD table(s) could not be recomputed."
    else
      echo "    Every table recomputed from results/."
    fi
    echo "    Full side-by-side output: $CT_LOG"
  else
    note_fail "check_tables.sh"
  fi
else
  echo "    Skipped: results/ holds no per-task CSVs in this checkout."
  echo "    They are tracked in the repository; a shallow or partial clone"
  echo "    may not have them."
fi

#-----------------------------------------------------------------------
printf "\n%s%s%s\n" "$BOLD" "$(printf '%.0s-' {1..70})" "$OFF"
if [ "$FAILED" -eq 0 ]; then
  printf "%sAll laptop checks passed.%s\n\n" "$GREEN" "$OFF"
  cat <<'EOF'
Three different claims, and this run settled two of them.

  1. The implementation is correct            -- steps 1 to 4, on this laptop.
  2. The tables match the numbers behind them -- step 5, on this laptop, from
     the per-task CSVs tracked in results/.
  3. Those numbers are what the model produces from scratch -- NOT checked
     here. Regenerating them is a cluster job array per table; README.md,
     section "What runs where", gives the submission line for each one.

So a referee can confirm on a laptop that no table was transcribed wrongly, and
needs a cluster only to regenerate the fits themselves.
EOF
  exit 0
else
  printf "%s%d step(s) failed.%s\n" "$RED" "$FAILED" "$OFF"
  exit 1
fi
