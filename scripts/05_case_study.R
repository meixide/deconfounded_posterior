#!/usr/bin/env Rscript
## 05_case_study.R --------------------------------------------------------
##
## Leave-one-environment-out case study.
##
## Runs against either dataset of Section 3.2 — the proprietary quiron BMI data
## or the public BRFSS mirror — through the common `load_case_data()` in
## R/case_data.R.  Nothing below this line knows which.  See
## HANDOFF_REAL_DATA.md.
##
## This replaces the single-held-out-domain analysis of the submitted version.
## Referee 2's objection to that analysis was correct on two counts: coverage
## computed over individuals within one domain has a standard error far larger
## than it appears, because every interval shares one posterior; and a single
## held-out domain cannot validate a claim about generalising to new domains.
## Both are addressed by holding out each environment in turn and treating the
## environment as the unit of replication (see R/loeo.R).
##
## Usage:
##   Rscript scripts/05_case_study.R --data=PATH [--folds=all|N] [--max-n=NUM]
##
## Flags:
##   --data        path to the CSV
##   --dataset     "auto", "quiron" or "brfss"
##   --folds       "all", an integer N for the N largest environments, or a
##                 comma-separated list of environment names
##   --fold-index  run only the Nth largest eligible environment; for job
##                 arrays, pass SLURM_ARRAY_TASK_ID
##   --pa-numeric  BRFSS only; 0 expands the physical-activity score into
##                 dummies instead of treating it as numeric
##   --model       "gi_hd_slope" (default) or "gi_hd"
##   --sd-b-scale  half-normal scale for the spread of the per-environment
##                 slopes; slope model only
##   --min-n0      minimum size for an environment to serve as a target
##   --max-n       optionally subsample the *training* rows (see note below)
##   --max-target  cap on target rows scored per fold
##   --ncp         non-centring in [0, 1]; 0 for this data, see below
##   --v-prior-shape,--v-prior-rate  inverse-gamma prior on the free part of
##                 the residual variance; shape 0, the default, selects the
##                 improper Jeffreys prior.  See the note below.
##   --eta-lkj     LKJ concentration for the environment-mean correlations
##   --chains,--iter,--cores  passed to fit_bgi()
##   --summary-only  print the design summary and stop, fitting nothing
##   --out         output directory
##
## On `--ncp`: the simulations of Section 3.1 need `ncp = 1`, because their
## likelihood is weak and the centred geometry is then unsamplable.  On the
## case-study data the likelihood is about as informative as a likelihood gets
## and the ordering reverses — `tests/test_case_study_timing.R` measures the
## centred parameterisation at 13x faster per iteration, with no max-treedepth
## hits.  The default here is therefore 0, not the simulation default.
##
## A note on sample size.  The submitted script subsampled 2000 training and
## 2000 test rows from roughly 500,000.  That is no longer necessary: the
## likelihood in stan/gi_hd.stan depends on the data only through
## per-environment sufficient statistics, so the cost of a fit is O(E p^3)
## regardless of N.  Use all the data; `--max-n` exists only for smoke tests.

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)

data_path <- parse_flag(cli, "data",
                        file.path(dirname(root), "old_code", "quiron",
                                  "quiron_final.csv"))
dataset <- parse_flag(cli, "dataset", "auto")
## BRFSS only: treat the physical-activity category as a numeric score, as the
## submitted analysis treated `af`, or expand it into dummies.  Ignored by the
## quiron loader.
pa_numeric <- !identical(parse_flag(cli, "pa-numeric", "1"), "0")
## BRFSS only: restrict to a clinical subgroup. See `load_brfss_data()`.
subset_arg <- parse_flag(cli, "subset", "none")
folds_arg <- parse_flag(cli, "folds", "all")
fold_index <- parse_flag(cli, "fold-index", "")
min_n0 <- as.integer(parse_flag(cli, "min-n0", "300"))
max_n <- as.integer(parse_flag(cli, "max-n", "0"))
max_target <- as.integer(parse_flag(cli, "max-target", "2000"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
iter <- as.integer(parse_flag(cli, "iter", "2000"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
cov_method <- parse_flag(cli, "cov-method", "pooled")
## Which Stan model. `gi_hd_slope` is the reparameterisation of HANDOFF §3.8:
## each environment carries its own slope `b_e = Sigma_e^{-1} K`, so the noisy
## plug-in `Sigma_hat_e^{-1}` no longer enters the mean as a generated
## regressor. It beats `gi_hd` on gamma coverage, RMSE and false discoveries at
## every heterogeneity level tested (job 8729750), so it is the default here.
model_name <- parse_flag(cli, "model", "gi_hd_slope")
sd_b_scale <- as.numeric(parse_flag(cli, "sd-b-scale", "1"))
ncp <- as.numeric(parse_flag(cli, "ncp", "0"))
## Prior on the free part of the residual variance.  A shape of zero selects
## the Jeffreys prior p(v) propto 1/v, which is improper, and the posterior is
## then improper too whenever the target's Mahalanobis term exceeds every
## training one: the training conditional variances stay bounded away from zero
## as v -> 0, so the likelihood does not vanish there and the integral of dv/v
## diverges.  tests/test_posterior_propriety.R demonstrates it.  Any positive
## rate removes the problem, since exp(-rate/v) vanishes faster than v^{-a-1}
## grows; the shape controls how much the prior is felt where the data are.
v_prior_shape <- as.numeric(parse_flag(cli, "v-prior-shape", "3"))
v_prior_rate <- as.numeric(parse_flag(cli, "v-prior-rate", "2"))
eta_lkj <- as.numeric(parse_flag(cli, "eta-lkj", "2"))
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.95"))
seed <- as.integer(parse_flag(cli, "seed", "20260728"))
summary_only <- any(grepl("^--summary-only$", cli))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "case_study"))

## ---- Data ---------------------------------------------------------------
##
## `load_case_data()` lives in R/case_data.R and handles both the proprietary
## quiron file and the public BRFSS mirror.  It returns `x` (numeric design
## matrix, no intercept), `y`, `z` (environment labels), `covariate_names` and
## `confounders`; nothing below depends on which dataset produced them.

set.seed(seed)
dat <- load_case_data(data_path, dataset = dataset, pa_numeric = pa_numeric,
                      subset = subset_arg)

if (max_n > 0L && max_n < length(dat$y)) {
  message("Subsampling to ", max_n, " rows (smoke test only).")
  keep <- sort(sample(seq_along(dat$y), max_n))
  dat$x <- dat$x[keep, , drop = FALSE]
  dat$y <- dat$y[keep]
  dat$z <- dat$z[keep]
}

p <- ncol(dat$x)
env_sizes <- sort(table(dat$z), decreasing = TRUE)

message("\n--- Case study data ---")
message("observations      : ", length(dat$y))
message("covariates (p)    : ", p)
message("environments      : ", length(env_sizes))
message("environment sizes : min ", min(env_sizes), ", median ",
        stats::median(env_sizes), ", max ", max(env_sizes))
## The quantity the simulations identified as the one to watch: identification
## needs E >= p + 1, but reliable inference needs headroom above it.
message("between-env residual df, E - (p + 1) = ",
        length(env_sizes) - 1L - (p + 1L), " per fold")
if (length(env_sizes) - 1L - (p + 1L) < 5L) {
  warning("Few residual degrees of freedom in the between-environment ",
          "regression. scripts/04_sim_environment_budget.R quantifies what ",
          "this costs; interpret selections cautiously.", call. = FALSE)
}

## `--folds` takes "all", a count (the N largest environments), or an explicit
## comma-separated list of environment names.  `--fold-index` takes the Nth
## largest and exists so that a SLURM array can pass SLURM_ARRAY_TASK_ID
## straight through without knowing the environment names.
eligible <- names(env_sizes)[env_sizes >= min_n0]
folds <- if (nzchar(fold_index)) {
  i <- as.integer(fold_index)
  if (is.na(i) || i < 1L || i > length(eligible)) {
    stop("--fold-index=", fold_index, " is outside 1:", length(eligible),
         call. = FALSE)
  }
  eligible[i]
} else if (identical(folds_arg, "all")) {
  eligible
} else if (grepl("^[0-9]+$", folds_arg)) {
  utils::head(eligible, as.integer(folds_arg))
} else {
  requested <- trimws(strsplit(folds_arg, ",", fixed = TRUE)[[1]])
  unknown <- setdiff(requested, names(env_sizes))
  if (length(unknown) > 0L) {
    stop("Unknown environment(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  requested
}
message("eligible envs (n >= ", min_n0, "): ", length(eligible))
message("folds to evaluate : ", length(folds),
        if (length(folds) <= 4L) paste0(" (", paste(folds, collapse = ", "), ")")
        else "", "\n")

message("covariates:")
for (j in seq_len(p)) {
  col <- dat$x[, j]
  message(sprintf("  %-28s mean %8.3f  sd %7.3f", dat$covariate_names[j],
                  mean(col), stats::sd(col)))
}

if (summary_only) {
  message("\n--summary-only: stopping before the fits.")
  quit(save = "no", status = 0)
}

## ---- Fit ----------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

fold_results <- loeo_evaluate(
  dat$x, dat$y, dat$z,
  model = model,
  folds = folds,
  min_n0 = min_n0,
  max_target = max_target,
  cov_method = cov_method,
  ncp = ncp,
  v_prior_shape = v_prior_shape,
  v_prior_rate = v_prior_rate,
  eta_lkj = eta_lkj,
  sd_b_scale = sd_b_scale,
  adapt_delta = adapt_delta,
  chains = chains,
  iter = iter,
  cores = cores,
  seed = seed
)

bgi_write_csv(fold_results, file.path(out_dir, "loeo_folds.csv"))

## ---- Report -------------------------------------------------------------

## The per-fold CSV is already written, so a failure here loses no data. When a
## task runs a single fold and that fold misses the diagnostic thresholds,
## `loeo_summary()` correctly refuses to summarise it — but that must not mark
## the SLURM task FAILED, or a run whose data are intact looks like a run that
## crashed. The aggregate summary is produced by
## `scripts/07_aggregate_case_study.R` across all folds anyway.
summ <- tryCatch(loeo_summary(fold_results), error = function(e) {
  message("Per-fold summary skipped: ", conditionMessage(e))
  message("The fold CSV is written; aggregate with scripts/07_aggregate_case_study.R.")
  NULL
})
if (is.null(summ)) {
  quit(save = "no", status = 0)
}
stability <- loeo_selection_stability(fold_results, dat$covariate_names)
bgi_write_csv(summ, file.path(out_dir, "loeo_summary.csv"))
bgi_write_csv(stability, file.path(out_dir, "selection_stability.csv"))

cat("\n=== Leave-one-environment-out, ", nrow(fold_results),
    " held-out domains ===\n", sep = "")
cat("Standard errors are across environments, not across individuals:\n")
cat("within a domain every interval shares one posterior, so individual\n")
cat("coverage indicators are strongly dependent.\n\n")
print(summ, row.names = FALSE, digits = 3)

cat("\n=== Selection stability across held-out domains ===\n")
cat("Fraction of folds in which each covariate was selected.\n")
cat("A covariate reported as causal from a single split should be\n")
cat("selected in most folds; one that is not, was a split artefact.\n\n")
print(stability, row.names = FALSE, digits = 3)

cat("\n=== Predictive scale ===\n")
cat("S0 is the target-domain predictive sd, sigma_cond the training one.\n")
cat("They differ whenever Sigma_0 differs from the training covariances;\n")
cat("the submitted code used sigma_cond for both.\n\n")
scale_tab <- fold_results[, c("target_env", "sigma_cond_mean", "S0_mean",
                              "sigma_y_mean", "coverage")]
print(utils::head(scale_tab[order(-scale_tab$S0_mean), ], 10),
      row.names = FALSE, digits = 3)

cat("\nWritten to ", out_dir, "\n", sep = "")
