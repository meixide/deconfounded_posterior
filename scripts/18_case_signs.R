#!/usr/bin/env Rscript
## 18_case_signs.R ---------------------------------------------------------
##
## One fit of the slope model on the full BRFSS diabetic sample, all 52
## environments, no held-out state, to report the *sign* of every causal
## coefficient.
##
## Why this exists.  The leave-one-state-out study (05_case_study.R) stores only
## which covariates each fold selects, not their posterior summaries, so Table 4
## says that light drinking and heavy drinking are both selected but not in which
## direction each moves BMI.  The selected set is stable across folds, so one
## fit on all environments is enough to read the signs; it does not replace the
## leave-one-out evaluation.
##
## The alcohol categories are dummies against non-drinkers (R/case_data.R), so
## each coefficient is the effect of that band relative to drinking nothing.
##
## Usage (from new_code/):
##   Rscript scripts/18_case_signs.R --data=../data/brfss/brfss2023_case.csv
##
## Flags mirror 05_case_study.R; the defaults are the case-study configuration
## of Section 3.2 and Supplement C (subset diabetic, gi_hd_slope, ncp 0,
## adapt_delta 0.99, 4 chains of 4000 iterations).
##
## Output, in --out:
##   case_signs.csv     one row per covariate: posterior mean, 95% interval,
##                      P(gamma > 0), lfsr, declared sign at alpha
##   case_signs_draws.rds  gamma draws with covariate names, so nothing has to
##                      be refitted to compute a different summary

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
                        file.path(dirname(root), "data", "brfss",
                                  "brfss2023_case.csv"))
subset_arg <- parse_flag(cli, "subset", "diabetic")
model_name <- parse_flag(cli, "model", "gi_hd_slope")
cov_method <- parse_flag(cli, "cov-method", "pooled")
sd_b_scale <- as.numeric(parse_flag(cli, "sd-b-scale", "1"))
ncp <- as.numeric(parse_flag(cli, "ncp", "0"))
eta_lkj <- as.numeric(parse_flag(cli, "eta-lkj", "2"))
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.99"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
iter <- as.integer(parse_flag(cli, "iter", "4000"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
alpha <- as.numeric(parse_flag(cli, "alpha", "0.05"))
seed <- as.integer(parse_flag(cli, "seed", "20260728"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "case_signs"))

## ---- Data ---------------------------------------------------------------

set.seed(seed)
dat <- load_case_data(data_path, dataset = "brfss", pa_numeric = TRUE,
                      subset = subset_arg)
p <- ncol(dat$x)
message("observations : ", length(dat$y))
message("covariates   : ", p)
message("environments : ", length(unique(dat$z)))

## ---- Fit ----------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

## No target domain: every environment is used for training, and only the
## posterior of gamma is read.
fit <- fit_bgi(dat$x, dat$y, dat$z, x0 = NULL,
               model = model,
               cov_method = cov_method,
               sd_b_scale = sd_b_scale,
               ncp = ncp,
               eta_lkj = eta_lkj,
               adapt_delta = adapt_delta,
               chains = chains,
               iter = iter,
               cores = cores,
               seed = seed)

gamma <- fit$draws$gamma
colnames(gamma) <- dat$covariate_names
saveRDS(gamma, file.path(out_dir, "case_signs_draws.rds"))

## ---- Summaries ------------------------------------------------------------

## lfsr_j = min{P(gamma_j >= 0 | D), P(gamma_j <= 0 | D)}, equation (6).
p_pos <- colMeans(gamma > 0)
p_neg <- colMeans(gamma < 0)
lfsr <- pmin(p_pos, p_neg)
declared <- ifelse(lfsr < alpha, ifelse(p_pos > p_neg, "+", "-"), "abstain")

tab <- data.frame(
  covariate = dat$covariate_names,
  mean = colMeans(gamma),
  q025 = apply(gamma, 2, stats::quantile, probs = 0.025),
  q975 = apply(gamma, 2, stats::quantile, probs = 0.975),
  p_positive = p_pos,
  lfsr = lfsr,
  declared = declared,
  stringsAsFactors = FALSE
)
bgi_write_csv(tab, file.path(out_dir, "case_signs.csv"))

cat("\n=== Sampler diagnostics ===\n")
cat("max Rhat      : ", signif(fit$diagnostics$max_rhat, 4), "\n", sep = "")
cat("min bulk ESS  : ", round(fit$diagnostics$min_ess_bulk), "\n", sep = "")
cat("divergences   : ", fit$diagnostics$n_divergent, "\n", sep = "")

cat("\n=== Causal coefficients (BMI units per unit of covariate) ===\n")
cat("Alcohol bands are relative to non-drinkers. Declared sign at alpha = ",
    alpha, ".\n\n", sep = "")
print(tab, row.names = FALSE, digits = 3)

cat("\n=== Alcohol ===\n")
print(tab[grepl("^alcohol", tab$covariate), ], row.names = FALSE, digits = 3)

cat("\nWritten to ", out_dir, "\n", sep = "")
