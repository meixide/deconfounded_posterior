#!/usr/bin/env Rscript
## 08_case_baselines.R ----------------------------------------------------
##
## Leave-one-environment-out for the frequentist baselines alone.
##
## These need no MCMC — they are least squares and two-stage least squares —
## so the whole 52-fold sweep runs in seconds and does not have to wait behind
## the Bayesian fits.  It exists for two reasons.
##
## First, it is the cheap half of the comparison the revision checklist asks
## for: IV as a baseline, matching GI on `gamma` and losing on prediction.
## Running it separately means the baseline numbers are available while the
## sampler jobs are still queued.
##
## Second, it checks on real data the identity `tests/test_baseline_identities.R`
## proves in simulation: the frequentist GI estimate of `gamma` *is* 2SLS with
## the environment indicators as instruments.  If that identity holds fold by
## fold here, the "pooled OLS on X and the environment means" baseline Referee 2
## proposes and the multi-source IV comparison they ask for are the same
## estimator, and the paper can say so.
##
## Usage:
##   Rscript scripts/08_case_baselines.R --data=PATH [--dataset=auto]

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
dataset <- parse_flag(cli, "dataset", "auto")
## BRFSS only: treat the physical-activity category as a numeric score (as the
## submitted analysis treated `af`) or expand it into dummies.
pa_numeric <- !identical(parse_flag(cli, "pa-numeric", "1"), "0")
min_n0 <- as.integer(parse_flag(cli, "min-n0", "300"))
level <- as.numeric(parse_flag(cli, "level", "0.95"))
seed <- as.integer(parse_flag(cli, "seed", "20260728"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "case_baselines"))

set.seed(seed)
dat <- load_case_data(data_path, dataset = dataset, pa_numeric = pa_numeric)
p <- ncol(dat$x)
sizes <- sort(table(dat$z), decreasing = TRUE)
folds <- names(sizes)[sizes >= min_n0]

message("N ", length(dat$y), " | p ", p, " | E ", length(sizes),
        " | folds ", length(folds))

rows <- list()
gap <- numeric(0)
for (e0 in folds) {
  is0 <- dat$z == e0
  x_tr <- dat$x[!is0, , drop = FALSE]
  y_tr <- dat$y[!is0]
  z_tr <- dat$z[!is0]
  x0 <- dat$x[is0, , drop = FALSE]
  y0 <- dat$y[is0]

  ols <- fit_ols(x_tr, y_tr, x0, level = level)
  pgi <- fit_pooled_gi(x_tr, y_tr, z_tr, x0, level = level)
  iv <- fit_iv_2sls(x_tr, y_tr, z_tr, x0, level = level)

  ## The identity under test, on the scale of the coefficients themselves.
  gap <- c(gap, max(abs(pgi$gamma - iv$gamma)) / max(abs(iv$gamma)))

  row <- data.frame(target_env = e0, n_train = nrow(x_tr), n_target = nrow(x0),
                    stringsAsFactors = FALSE)
  for (nm in c("ols", "pooled_gi", "iv")) {
    f <- get(if (nm == "pooled_gi") "pgi" else nm)
    m <- interval_metrics(y0, f$pred_lower, f$pred_upper, level = level)
    row[[paste0(nm, "_coverage")]] <- m$coverage
    row[[paste0(nm, "_interval_score")]] <- m$interval_score
    row[[paste0(nm, "_rmse")]] <- sqrt(mean((f$pred_mean - y0)^2))
  }
  rows[[length(rows) + 1L]] <- row
}

res <- do.call(rbind, rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bgi_write_csv(res, file.path(out_dir, "baseline_folds.csv"))

cat("\n=== Frequentist baselines, leave-one-environment-out ===\n")
cat("Standard errors are across the ", nrow(res),
    " held-out environments.\n\n", sep = "")
summ <- do.call(rbind, lapply(c("ols", "pooled_gi", "iv"), function(m) {
  cov <- res[[paste0(m, "_coverage")]]
  data.frame(method = m,
             coverage = mean(cov),
             coverage_se = stats::sd(cov) / sqrt(length(cov)),
             interval_score = mean(res[[paste0(m, "_interval_score")]]),
             rmse = mean(res[[paste0(m, "_rmse")]]),
             stringsAsFactors = FALSE)
}))
print(summ[order(summ$interval_score), ], row.names = FALSE, digits = 4)
bgi_write_csv(summ, file.path(out_dir, "baseline_summary.csv"))

cat("\n=== Is frequentist GI the same estimator as 2SLS? ===\n")
cat(sprintf("max relative discrepancy in gamma over %d folds: %.3e\n",
            length(gap), max(gap)))
cat("Machine precision here means the two baselines Referee 2 asks to see\n")
cat("separately are algebraically one estimator, so reporting both adds\n")
cat("nothing beyond the check itself.\n")

cat("\nWritten to ", out_dir, "\n", sep = "")
