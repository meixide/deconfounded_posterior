#!/usr/bin/env Rscript
## 02_aggregate_support_recovery.R ----------------------------------------
##
## Collect the per-task CSVs written by 01_sim_support_recovery.R and produce
## the summary tables for the manuscript.
##
## Three tables are written:
##   support_recovery_by_method.csv   method x scenario, the selection metrics
##   calibration_by_scenario.csv      scenario, the coverage metrics
##   diagnostics_by_scenario.csv      fit quality, so that nothing is averaged
##                                    in silently
##
## Replications whose BGI fit failed its diagnostics (divergent transitions,
## or max Rhat above `--rhat-max`) are reported separately and, by default,
## excluded from the averages.  Reporting the exclusion count matters: a
## method that looks well calibrated only after discarding a third of its fits
## is not well calibrated.
##
## Usage:
##   Rscript scripts/02_aggregate_support_recovery.R [--in=DIR] [--out=DIR]
##                                                   [--rhat-max=1.01]
##                                                   [--keep-bad]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
## Reads the task CSVs and fits nothing, so it does not need Stan.
root <- bgi_bootstrap(need_stan = FALSE)

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}

cli <- commandArgs(trailingOnly = TRUE)
in_dir <- parse_flag(cli, "in", file.path(root, "results", "support_recovery"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "summaries"))
rhat_max <- as.numeric(parse_flag(cli, "rhat-max", "1.01"))
max_div_frac <- as.numeric(parse_flag(cli, "max-divergent-frac", "0.01"))
keep_bad <- "--keep-bad" %in% cli

files <- list.files(in_dir, pattern = "^task_[0-9]+\\.csv$", full.names = TRUE)
if (length(files) == 0L) {
  stop("No task CSVs found in ", in_dir, call. = FALSE)
}
message("Reading ", length(files), " task files from ", in_dir)

raw <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))

## Column renamed once the distinction between sign errors and false
## discoveries was made explicit: what the lfsr rules average is the posterior
## expected *sign-error* rate, not an FDR.  Accept the old name so that result
## directories produced before the rename still aggregate.
if ("posterior_efdr" %in% names(raw) && !"posterior_esr" %in% names(raw)) {
  names(raw)[names(raw) == "posterior_efdr"] <- "posterior_esr"
  message("Note: renamed legacy column 'posterior_efdr' to 'posterior_esr'.")
}
required <- c("tpr", "fpr", "fdp", "fwer", "sign_error", "posterior_esr")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols) > 0L) {
  stop("Task files are missing: ", paste(missing_cols, collapse = ", "),
       ". Re-run scripts/01_sim_support_recovery.R.", call. = FALSE)
}

## ---- Fit quality screening ---------------------------------------------

fits <- raw[!duplicated(raw$task_id), ]
## Screen on the divergence *rate*, as loeo_summary() already does for the case
## study with max_divergent_frac = 0.01.  Requiring zero divergences is a much
## harsher rule and was excluding fits that are fine by every other measure: in
## the weak-identifiability scenarios, replications with one to eleven
## divergences out of four thousand draws and Rhat below 1.01 were all being
## discarded, leaving three of ten.  A rate is what the recorded draw count is
## there for.
##
## Runs predating that column get it from a MANIFEST file in their own
## directory, read here.  The alternative -- passing --post-draws on the command
## line -- makes the published table depend on an argument the person
## reproducing it has to know, and check_tables.sh did not know it: it recomputed
## Table 2 under the zero-divergence gate, reported three of ten fits usable
## where the table used nine, and so failed to reproduce the very table it exists
## to verify.  A number the directory carries cannot be forgotten.
bgi_manifest_value <- function(dir, key) {
  path <- file.path(dir, "MANIFEST")
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  hit <- grep(paste0("^", key, "="), lines, value = TRUE)
  if (length(hit) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(sub(paste0("^", key, "="), "", hit[1])))
}
post_draws <- if ("bgi_post_draws" %in% names(fits)) {
  fits$bgi_post_draws
} else {
  from_flag <- as.numeric(parse_flag(cli, "post-draws", "NA"))
  from_manifest <- bgi_manifest_value(in_dir, "post_draws")
  rep(if (!is.na(from_flag)) from_flag else from_manifest, nrow(fits))
}
if (all(is.na(post_draws))) {
  message("No draw count recorded, no MANIFEST and --post-draws not given: ",
          "falling back to requiring zero divergences.  This is NOT the screen ",
          "the published table used.")
  div_ok <- fits$bgi_divergent == 0
} else {
  div_ok <- fits$bgi_divergent / post_draws <= max_div_frac
}
fits$ok <- div_ok & fits$bgi_max_rhat <= rhat_max
quality <- aggregate(
  cbind(n_fits = rep(1, nrow(fits)), n_ok = as.integer(fits$ok)) ~
    scenario_id + label + confounding + identifiability + n_e,
  data = fits, FUN = sum
)
quality$prop_ok <- quality$n_ok / quality$n_fits

timing <- aggregate(
  cbind(bgi_seconds, bgi_max_rhat, bgi_min_ess, cov_condition_max,
        cov_condition_target, cov_shrinkage_max, mu_min_sv, mu_condition) ~
    scenario_id,
  data = fits, FUN = mean
)
diagnostics <- merge(quality, timing, by = "scenario_id")

if (!keep_bad) {
  bad <- fits$task_id[!fits$ok]
  if (length(bad) > 0L) {
    message("Excluding ", length(bad), " of ", nrow(fits),
            " replications that failed the sampler diagnostics ",
            "(use --keep-bad to retain them).")
    raw <- raw[!raw$task_id %in% bad, ]
  }
}

## ---- Selection metrics, by method and scenario --------------------------

group <- c("scenario_id", "label", "confounding", "identifiability", "n_e",
           "p", "s0", "n_env", "method")
## `sign_error` and `fdp` measure different things and must always be read
## together: the first is what the lfsr rules control, the second is what they
## do not.  `posterior_esr` is the rule's own claim about the first, so
## posterior_esr versus sign_error is the calibration check.
metrics <- c("tpr", "fpr", "fdp", "fwer", "exact_recovery", "jaccard", "mcc",
             "n_selected", "sign_error", "posterior_esr")

by_method <- aggregate(
  raw[, metrics],
  by = raw[, group],
  FUN = function(v) mean(v, na.rm = TRUE)
)
counts <- aggregate(list(n_rep = raw$tpr), by = raw[, group], FUN = length)
by_method <- merge(by_method, counts, by = group)

## Monte Carlo standard error of the headline quantities, so that a reader can
## tell an effect from noise at this number of replications.
se_of <- function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
se_tab <- aggregate(raw[, c("tpr", "fdp", "fwer", "sign_error")],
                    by = raw[, group], FUN = se_of)
names(se_tab)[names(se_tab) %in% c("tpr", "fdp", "fwer", "sign_error")] <-
  c("tpr_se", "fdp_se", "fwer_se", "sign_error_se")
by_method <- merge(by_method, se_tab, by = group)
by_method <- by_method[order(by_method$scenario_id, by_method$method), ]

## ---- Calibration, by scenario ------------------------------------------

cal_vars <- c("bgi_gamma_coverage", "bgi_gamma_coverage_parents",
              "bgi_gamma_coverage_nulls", "bgi_gamma_rmse",
              "ols_gamma_rmse", "pooled_gi_gamma_rmse",
              "iv_gamma_rmse", "anchor_g8_gamma_rmse",
              "anchor_oracle_gamma_rmse", "anchor_oracle_gamma_anchor",
              "bgi_pred_coverage", "bgi_pred_width", "bgi_pred_rmse",
              "ols_pred_coverage", "ols_pred_width",
              "pooled_gi_pred_coverage", "pooled_gi_pred_width",
              "group_dro_gamma_rmse", "vrex_gamma_rmse", "wass_dro_gamma_rmse",
              unlist(lapply(
                c("bgi", "pooled_gi", "iv", "anchor_g8", "anchor_oracle",
                  "group_dro", "vrex", "wass_dro", "ols"),
                function(m) paste0(m, c("_pred_coverage", "_pred_width",
                                        "_pred_is", "_pred_rmse")))))
## Result directories produced before the distribution-shift baselines were
## added lack their columns; aggregate whatever is present rather than failing.
absent <- setdiff(cal_vars, names(raw))
if (length(absent) > 0L) {
  message("Note: these columns are absent and will be skipped (result ",
          "directory predates the distribution-shift baselines): ",
          paste(absent, collapse = ", "))
  cal_vars <- intersect(cal_vars, names(raw))
}
one_per_rep <- raw[!duplicated(raw$task_id), ]
calibration <- aggregate(
  one_per_rep[, cal_vars],
  by = one_per_rep[, c("scenario_id", "label", "confounding",
                       "identifiability", "n_e")],
  FUN = function(v) mean(v, na.rm = TRUE)
)
calibration <- calibration[order(calibration$scenario_id), ]

## ---- Write out ---------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
## Which Stan model produced these rows.  The two parameterisations are not
## the same model when the covariances differ across environments, and Table 3
## shows the choice matters for the calibration of gamma, so a summary that
## does not say which one it summarises invites the reader to assume.
if ("model" %in% names(raw)) {
  cat("\nStan model(s) behind these results: ",
      paste(sort(unique(raw$model)), collapse = ", "), "\n", sep = "")
  if (length(unique(raw$model)) > 1L) {
    cat("More than one: the summary mixes parameterisations.\n")
  }
}

bgi_write_csv(by_method, file.path(out_dir, "support_recovery_by_method.csv"))
bgi_write_csv(calibration, file.path(out_dir, "calibration_by_scenario.csv"))
bgi_write_csv(diagnostics, file.path(out_dir, "diagnostics_by_scenario.csv"))

## ---- Console report ----------------------------------------------------

cat("\n=== Support recovery (mean over replications) ===\n")
cat("sign_error: wrong direction among selected TRUE parents",
    "-- what the lfsr rules control\n")
cat("fdp       : exact zeros among all selected",
    "-- what they do NOT control\n\n")
show <- by_method[, c("label", "method", "n_rep", "tpr", "sign_error",
                      "fdp", "fpr", "fwer", "exact_recovery", "n_selected")]
show[, 4:10] <- round(show[, 4:10], 3)
print(show, row.names = FALSE)

## Is the rule's own claim about sign errors borne out?  posterior_esr is the
## posterior expected sign-error rate over the selected set; sign_error is the
## realised one.  A realised rate above the claim means the posterior is
## overconfident, which relabelling the estimand does not fix.
cat("\n=== Sign-error calibration (claimed vs realised) ===\n")
sg <- by_method[grepl("^bgi", by_method$method), ]
cal_sign <- data.frame(
  label = sg$label, method = sg$method,
  claimed = round(sg$posterior_esr, 4),
  realised = round(sg$sign_error, 4),
  mc_se = round(sg$sign_error_se, 4),
  overconfident = ifelse(is.na(sg$sign_error), NA,
                         sg$sign_error > sg$posterior_esr + sg$sign_error_se)
)
print(cal_sign, row.names = FALSE)

## Every method that produces an interval, on one footing.  Coverage cannot
## rank methods on its own -- a wide interval always covers -- so the interval
## score is the column to read.  `iv` and `pooled_gi` share `gamma` exactly and
## differ only by the K correction, so the gap between their rows measures what
## K contributes to prediction with the slope held fixed.
pred_methods <- c("bgi", "pooled_gi", "iv", "anchor_g8", "anchor_oracle",
                  "group_dro", "vrex", "wass_dro", "ols")
have <- pred_methods[
  paste0(pred_methods, "_pred_is") %in% names(calibration)]

if (length(have) > 0L) {
  cat("\n=== Target-domain prediction (nominal coverage 0.95;",
      "interval score lower is better) ===\n")
  for (lab in calibration$label) {
    row <- calibration[calibration$label == lab, ]
    tab <- data.frame(
      method = have,
      rmse = round(vapply(have, function(m)
        row[[paste0(m, "_pred_rmse")]], numeric(1)), 3),
      coverage = round(vapply(have, function(m)
        row[[paste0(m, "_pred_coverage")]], numeric(1)), 3),
      width = round(vapply(have, function(m)
        row[[paste0(m, "_pred_width")]], numeric(1)), 2),
      interval_score = round(vapply(have, function(m)
        row[[paste0(m, "_pred_is")]], numeric(1)), 2)
    )
    cat("\n", lab, "\n", sep = "")
    print(tab[order(tab$interval_score), ], row.names = FALSE)
  }
  if (all(c("iv_pred_is", "pooled_gi_pred_is") %in% names(calibration))) {
    cat("\nWhat the K correction buys, at identical gamma",
        "(iv and pooled_gi differ only by it):\n")
    k_gain <- data.frame(
      label = calibration$label,
      iv_is = round(calibration$iv_pred_is, 2),
      pooled_gi_is = round(calibration$pooled_gi_pred_is, 2),
      iv_rmse = round(calibration$iv_pred_rmse, 3),
      pooled_gi_rmse = round(calibration$pooled_gi_pred_rmse, 3)
    )
    print(k_gain, row.names = FALSE)
  }
} else {
  cat("\n=== Target-domain predictive coverage (nominal 0.95) ===\n")
  cov_cols <- intersect(
    c("bgi_pred_coverage", "ols_pred_coverage", "pooled_gi_pred_coverage"),
    names(calibration))
  cal_show <- calibration[, c("label", cov_cols)]
  cal_show[, -1] <- round(cal_show[, -1], 3)
  print(cal_show, row.names = FALSE)
}

cat("\n=== Estimation of gamma: RMSE by method ===\n")
rmse_cols <- intersect(
  c("bgi_gamma_rmse", "pooled_gi_gamma_rmse", "iv_gamma_rmse",
    "anchor_g8_gamma_rmse", "anchor_oracle_gamma_rmse", "ols_gamma_rmse"),
  names(calibration))
rmse_show <- calibration[, c("label", rmse_cols)]
rmse_show[, -1] <- round(rmse_show[, -1], 4)
print(rmse_show, row.names = FALSE)
if (all(c("pooled_gi_gamma_rmse", "iv_gamma_rmse") %in% names(calibration))) {
  cat("\n(pooled_gi and iv agree to",
      format(max(abs(calibration$pooled_gi_gamma_rmse -
                       calibration$iv_gamma_rmse), na.rm = TRUE), digits = 2),
      "-- they are the same estimator, see R/baselines.R)\n")
}
if ("anchor_oracle_gamma_anchor" %in% names(calibration)) {
  cat("\nOracle-tuned anchor strength by scenario",
      "(no single value is uniformly best,\nand the target labels it needs",
      "are unavailable in practice):\n")
  print(calibration[, c("label", "anchor_oracle_gamma_anchor")],
        row.names = FALSE)
}

cat("\n=== Credible/confidence interval coverage for gamma ===\n")
gam_show <- calibration[, c("label", "bgi_gamma_coverage")]
gam_show[, -1] <- round(gam_show[, -1], 3)
print(gam_show, row.names = FALSE)

## ---- Table 2 body, ready to paste ---------------------------------------
## Emitted rather than transcribed: the published table was once out of step
## with every surviving run, and a table a script writes cannot drift from the
## numbers behind it.

fmt2 <- function(x) sub("^0", "", formatC(round(x, 2), format = "f", digits = 2))
anchor_row <- parse_flag(cli, "anchor", "anchor_g8")

sel_rows <- list(
  c("bgi_sign",   "\\textbf{ours}, sign rule"),
  c("bgi_ci",     "\\textbf{ours}, credible interval"),
  c("bgi_rope",   "\\textbf{ours}, ROPE"),
  c("pooled_gi",  "frequentist GI"),
  c("ols",        "OLS"),
  c("ols_bh",     "OLS, Benjamini--Hochberg"),
  c(anchor_row,   "anchor regression"),
  c("group_dro",  "group DRO"),
  c("iv",         "instrumental variables"),
  c("icp",        "invariant causal prediction"))

labs <- c("conf0_strong_n200", "conf2_strong_n200",
          "conf0_weak_n200", "conf2_weak_n200")

cat("\n=== Table 2 body, paste into jcgs.tex ===\n\n")
cat("\\multicolumn{7}{c}{\\emph{support recovery}: true positive rate / false discovery proportion}\\\\\n\\hline\n")
for (r in sel_rows) {
  m <- r[1]
  cells <- vapply(labs, function(l) {
    row <- by_method[by_method$method == m & by_method$label == l, , drop = FALSE]
    if (nrow(row) == 0L) "---" else
      sprintf("%s / %s", formatC(row$tpr[1], format = "f", digits = 2), fmt2(row$fdp[1]))
  }, character(1))
  agg <- by_method[by_method$method == m, , drop = FALSE]
  cat(sprintf("%s & %s & %s & %s\\\\\n", r[2], paste(cells, collapse = " & "),
              fmt2(mean(agg$fwer)), fmt2(mean(agg$mcc))))
}

pred <- list(c("bgi", "\\textbf{ours}"), c("pooled_gi", "frequentist GI"),
             c("ols", "OLS"), c(anchor_row, "anchor regression"),
             c("group_dro", "group DRO"), c("vrex", "V-REx"),
             c("wass_dro", "Wasserstein DRO"), c("iv", "instrumental variables"))
cat("\\specialrule{1pt}{0pt}{0pt}\n")
cat("\\multicolumn{7}{c}{\\emph{prediction}: coverage / interval score, nominal $95\\%$}\\\\\n\\hline\n")
for (r in pred) {
  cc <- paste0(r[1], "_pred_coverage"); ii <- paste0(r[1], "_pred_is")
  if (!(cc %in% names(calibration))) cc <- paste0(cc, ".1")
  cells <- vapply(labs, function(l) {
    row <- calibration[calibration$label == l, , drop = FALSE]
    if (nrow(row) == 0L || !(ii %in% names(row))) "---" else
      sprintf("%s / %s", fmt2(row[[cc]][1]), formatC(row[[ii]][1], format = "f", digits = 1))
  }, character(1))
  cat(sprintf("%s & %s & --- & ---\\\\\n", r[2], paste(cells, collapse = " & ")))
}
cat("\n")

cat("\n=== Fit diagnostics ===\n")
diag_show <- diagnostics[, c("label", "n_fits", "prop_ok", "bgi_seconds",
                             "bgi_max_rhat", "bgi_min_ess", "mu_min_sv",
                             "cov_condition_max")]
diag_show[, -1] <- round(diag_show[, -1], 3)
print(diag_show, row.names = FALSE)

cat("\nSummaries written to ", out_dir, "\n", sep = "")
