#!/usr/bin/env Rscript
## 17_aggregate_split.R ---------------------------------------------------
##
## Summarise the split-extrapolation runs of `scripts/16_split_extrapolation.R`.
##
## The question is whether the method earns its place against pooled OLS when
## the target is genuinely outside the training hull — which leave-one-
## environment-out cannot answer, because there the target sits inside the
## convex hull of the `E - 1` environments already seen and the pooled slope is
## near-optimal for it by construction.
##
## The estimand here is a *difference in differences*, and it has to be, because
## the far targets were selected by distance. Within replicate `r`, form the
## advantage of BGI over OLS separately on the far targets and on the near
## targets,
##
##     adv_r(role) = mean over targets of ( OLS metric - BGI metric )
##
## and take `adv_r(far) - adv_r(near)`. Reporting `adv_r(far)` alone would let
## the selection of far targets manufacture the result: far targets are harder
## for everything, so any method with wider intervals looks better there. The
## near targets are the control, drawn from the same training half in the same
## replicate, so the difference removes anything common to the split.
##
## Replicates are the unit of replication, exactly as environments are in
## `R/loeo.R`, so a paired t-interval over replicates is the right standard
## error and Referee 2's correlated-coverage objection is not reintroduced.
##
## Convergence is checked before anything is averaged. The slope model fails on
## a large fraction of leave-one-environment-out folds (see
## SLOPE_ON_REAL_DATA.md), so a pass rate is reported first and unconditionally.
##
## Usage:
##   Rscript scripts/17_aggregate_split.R --in=results/split_acs_tract_gi_hd_slope

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
# An aggregator reads CSVs; it neither fits nor simulates, so it does not
# need the fitting dependencies.  See `need_stan` in R/setup.R.
root <- bgi_bootstrap(need_stan = FALSE)

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
in_dir <- parse_flag(cli, "in", file.path(root, "results",
                                          "split_acs_tract_gi_hd_slope"))
rhat_max <- as.numeric(parse_flag(cli, "rhat-max", "1.01"))
ess_min <- as.numeric(parse_flag(cli, "ess-min", "400"))

files <- list.files(in_dir, "^task_.*\\.csv$", full.names = TRUE)
if (length(files) == 0L) {
  stop("No task_*.csv under ", in_dir, call. = FALSE)
}
d <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
d <- d[order(d$rep_id, d$role, d$rank), , drop = FALSE]

## ---- Convergence, before anything is averaged ---------------------------
d$ok <- !is.na(d$max_rhat) & d$max_rhat <= rhat_max & d$min_ess >= ess_min
cat(sprintf("\n=== %s ===\n", basename(in_dir)))
cat(sprintf("tasks %d | replicates %d | model %s | dataset %s\n",
            nrow(d), length(unique(d$rep_id)), d$model[1], d$dataset[1]))
cat(sprintf("convergence: %d of %d pass (Rhat <= %.3f, ESS >= %.0f)\n",
            sum(d$ok), nrow(d), rhat_max, ess_min))
cat(sprintf("  max Rhat %.3f | min ESS %.0f | total divergences %d\n",
            max(d$max_rhat, na.rm = TRUE), min(d$min_ess, na.rm = TRUE),
            sum(d$divergent)))
cat(sprintf("  pass rate by role: far %.2f, near %.2f\n",
            mean(d$ok[d$role == "far"]), mean(d$ok[d$role == "near"])))

if (mean(d$ok) < 0.5) {
  cat("\nFewer than half the fits converged. Everything below is reported for\n")
  cat("completeness but should not be quoted; see SLOPE_ON_REAL_DATA.md.\n")
}

## ---- The design, as drawn ------------------------------------------------
rep_tab <- unique(d[, c("rep_id", "n_train_env", "between_df", "lambda_min",
                        "lambda_cond", "n_lambda_min", "dist_min", "dist_max")])
cat("\n=== Identification of each training half (reported, not optimised) ===\n")
cat("The half is drawn at random; these are outputs of the draw.\n\n")
print(rep_tab, row.names = FALSE, digits = 3)

cat(sprintf("\ntarget distance: far %.2f (median), near %.3f (median), ratio %.0fx\n",
            stats::median(d$distance[d$role == "far"]),
            stats::median(d$distance[d$role == "near"]),
            stats::median(d$distance[d$role == "far"]) /
              stats::median(d$distance[d$role == "near"])))

## ---- Marginal performance ------------------------------------------------
use <- d[d$ok, , drop = FALSE]
if (nrow(use) == 0L) {
  stop("No converged task to summarise.", call. = FALSE)
}
cat("\n=== Performance by role, converged tasks only ===\n\n")
cat(sprintf("%-6s %6s %10s %10s %12s %12s %10s %10s\n", "role", "n",
            "BGI cov", "OLS cov", "BGI int.sc", "OLS int.sc", "BGI rmse",
            "OLS rmse"))
for (r in c("far", "near")) {
  u <- use[use$role == r, , drop = FALSE]
  if (nrow(u) == 0L) next
  cat(sprintf("%-6s %6d %10.4f %10.4f %12.4f %12.4f %10.4f %10.4f\n", r, nrow(u),
              mean(u$coverage), mean(u$ols_coverage),
              mean(u$interval_score), mean(u$ols_interval_score),
              mean(u$rmse), mean(u$ols_rmse)))
}

## ---- The difference in differences --------------------------------------
## Per replicate and role, the advantage of BGI over OLS; then far minus near.
adv <- function(u, metric) {
  mean(u[[paste0("ols_", metric)]] - u[[metric]])
}
reps <- sort(unique(use$rep_id))
rows <- list()
for (r in reps) {
  far <- use[use$rep_id == r & use$role == "far", , drop = FALSE]
  near <- use[use$rep_id == r & use$role == "near", , drop = FALSE]
  if (nrow(far) == 0L || nrow(near) == 0L) next
  rows[[length(rows) + 1L]] <- data.frame(
    rep_id = r, n_far = nrow(far), n_near = nrow(near),
    is_far = adv(far, "interval_score"), is_near = adv(near, "interval_score"),
    rmse_far = adv(far, "rmse"), rmse_near = adv(near, "rmse"),
    cov_far = mean(far$coverage), cov_near = mean(near$coverage),
    ols_cov_far = mean(far$ols_coverage), ols_cov_near = mean(near$ols_coverage),
    stringsAsFactors = FALSE)
}
if (length(rows) == 0L) {
  stop("No replicate has both a far and a near converged target.", call. = FALSE)
}
per_rep <- do.call(rbind, rows)

cat("\n=== Does the advantage grow with extrapolative distance? ===\n")
cat("Positive = BGI better. `far - near` is the estimand; `far` alone would\n")
cat("be confounded with the selection of far targets.\n\n")
report <- function(x, label) {
  if (length(x) < 2L || stats::sd(x) == 0) {
    cat(sprintf("  %-28s %+.4f  (n = %d, no interval)\n", label, mean(x),
                length(x)))
    return(invisible())
  }
  tt <- stats::t.test(x)
  cat(sprintf("  %-28s %+.4f  [%+.4f, %+.4f]  p = %.3f  (n = %d)\n",
              label, mean(x), tt$conf.int[1], tt$conf.int[2], tt$p.value,
              length(x)))
}
report(per_rep$is_far, "interval score, far")
report(per_rep$is_near, "interval score, near")
report(per_rep$is_far - per_rep$is_near, "DIFFERENCE (far - near)")
cat("\n")
report(per_rep$rmse_far, "RMSE, far")
report(per_rep$rmse_near, "RMSE, near")
report(per_rep$rmse_far - per_rep$rmse_near, "DIFFERENCE (far - near)")

cat("\n=== Coverage by role ===\n")
cat(sprintf("  far   : BGI %.4f  OLS %.4f\n", mean(per_rep$cov_far),
            mean(per_rep$ols_cov_far)))
cat(sprintf("  near  : BGI %.4f  OLS %.4f\n", mean(per_rep$cov_near),
            mean(per_rep$ols_cov_near)))

## ---- Continuous version --------------------------------------------------
## The dichotomy is a summary; the underlying relation is continuous, and a
## slope on log distance uses every target rather than only the extremes.
use$adv_is <- use$ols_interval_score - use$interval_score
use$adv_rmse <- use$ols_rmse - use$rmse
cat("\n=== Advantage against log distance, all converged targets ===\n")
for (nm in c("adv_is", "adv_rmse")) {
  m <- stats::lm(use[[nm]] ~ log(use$distance))
  ci <- stats::confint(m)[2, ]
  cat(sprintf("  %-10s slope %+.4f  [%+.4f, %+.4f]  p = %.3g  | Spearman %+.3f\n",
              nm, stats::coef(m)[2], ci[1], ci[2],
              summary(m)$coefficients[2, 4],
              stats::cor(log(use$distance), use[[nm]], method = "spearman")))
}
cat("\nA linear fit is easily dominated by one extreme target, so the Spearman\n")
cat("correlation is printed beside it; if the two disagree in sign, the slope\n")
cat("is a leverage artefact and should not be quoted.\n")

bgi_write_csv(per_rep, file.path(in_dir, "split_per_replicate.csv"))
bgi_write_csv(d, file.path(in_dir, "split_all_tasks.csv"))
cat("\nWritten to ", in_dir, "\n", sep = "")
