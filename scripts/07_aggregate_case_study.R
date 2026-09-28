#!/usr/bin/env Rscript
## 07_aggregate_case_study.R ----------------------------------------------
##
## Combine the per-fold outputs of a `slurm/06_case_study.sh` job array into
## the three tables Section 3.2 reports.
##
## The point of the aggregation, and the reason it is a separate step rather
## than an average taken inside the fitting script, is Referee 2's first
## objection.  Coverage within one held-out domain is a mean over individuals
## whose intervals all come from a single posterior, so those indicators are
## strongly dependent and the binomial standard error over individuals is far
## too small.  The environment is the unit of replication, so the standard
## error comes from the spread *across* the per-fold coverages computed here.
##
## Usage:
##   Rscript scripts/07_aggregate_case_study.R --in=results/case_brfss/folds \
##                                             --out=results/case_brfss \
##                                             --data=../data/brfss/brfss2023_case.csv
##
## `--data` is optional and is used only to label the covariates in the
## stability table; without it they are reported positionally.  `--expected=N`
## warns when fewer than `N` folds were found, which is the failure mode that
## matters: 44 of 52 folds still average to something that looks fine.

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
in_dir <- parse_flag(cli, "in", file.path(root, "results", "case_brfss",
                                          "folds"))
out_dir <- parse_flag(cli, "out", dirname(in_dir))
expected <- as.integer(parse_flag(cli, "expected", "0"))

files <- list.files(in_dir, pattern = "^loeo_folds\\.csv$",
                    recursive = TRUE, full.names = TRUE)
if (length(files) == 0L) {
  stop("No loeo_folds.csv under ", in_dir, call. = FALSE)
}

folds <- do.call(rbind, lapply(files, utils::read.csv,
                               stringsAsFactors = FALSE))
folds <- folds[order(folds$target_env), , drop = FALSE]

message("Aggregated ", nrow(folds), " folds from ", length(files), " files.")
if (expected > 0L && nrow(folds) < expected) {
  ## A silently short array is the failure mode that matters here: 44 of 52
  ## folds still average to something that looks fine.
  warning("Expected ", expected, " folds but found ", nrow(folds),
          ". Check logs/case-*.err for tasks that died.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bgi_write_csv(folds, file.path(out_dir, "loeo_folds.csv"))

summ <- loeo_summary(folds)
bgi_write_csv(summ, file.path(out_dir, "loeo_summary.csv"))

## Covariate names, in decreasing order of preference: given explicitly, or
## recovered by re-reading the dataset header (cheap next to the fits, and it
## guarantees the labels match the columns the folds were actually fitted on).
covariate_names <- strsplit(
  parse_flag(cli, "covariates", ""), ",", fixed = TRUE)[[1]]
data_path <- parse_flag(cli, "data", "")
if (length(covariate_names) != folds$p[1] && nzchar(data_path)) {
  covariate_names <- tryCatch(
    load_case_data(data_path,
                   dataset = parse_flag(cli, "dataset", "auto"),
                   pa_numeric = !identical(
                     parse_flag(cli, "pa-numeric", "1"), "0"))$covariate_names,
    error = function(e) {
      message("Could not recover covariate names: ", conditionMessage(e))
      character(0)
    })
}
if (length(covariate_names) != folds$p[1]) {
  message("Falling back to positional covariate labels; pass --data=PATH ",
          "or --covariates=a,b,c for readable ones.")
  covariate_names <- NULL
}
stability <- loeo_selection_stability(folds, covariate_names)
bgi_write_csv(stability, file.path(out_dir, "selection_stability.csv"))

cat("\n=== Leave-one-environment-out, ", nrow(folds),
    " held-out domains ===\n", sep = "")
cat("Standard errors are across environments, not across individuals:\n")
cat("within a domain every interval shares one posterior, so individual\n")
cat("coverage indicators are strongly dependent.\n\n")
print(summ, row.names = FALSE, digits = 3)

cat("\n=== Sampler diagnostics ===\n")
cat(sprintf("max Rhat over folds    : %.4f\n",
            max(folds$max_rhat, na.rm = TRUE)))
cat(sprintf("min ESS over folds     : %.0f\n",
            min(folds$min_ess, na.rm = TRUE)))
cat(sprintf("median runtime (s)     : %.0f\n",
            stats::median(folds$runtime_sec)))
## Printed in full rather than summarised to a pass/fail count: the exclusion
## threshold is deliberately loose (see `loeo_summary()`), so the distribution
## is what tells a reader whether that was reasonable.
cat(sprintf("folds with divergences : %d of %d\n",
            sum(folds$divergent > 0), nrow(folds)))
if (!is.null(folds$post_draws)) {
  cat(sprintf("divergence rate        : max %.4f, total %d in %d draws\n",
              max(folds$divergent / folds$post_draws, na.rm = TRUE),
              sum(folds$divergent), sum(folds$post_draws)))
}
cat("divergences per fold   : ")
print(table(folds$divergent))

cat("\n=== Selection stability across held-out domains ===\n")
cat("Fraction of folds in which each covariate was selected.\n\n")
print(stability, row.names = FALSE, digits = 3)

cat("\n=== Predictive scale ===\n")
cat("S0 is the target-domain predictive sd, sigma_cond the training one.\n")
cat("The submitted code used sigma_cond for both, which is the error that\n")
cat("invalidated the reported 0.95.\n\n")
scale_tab <- folds[, c("target_env", "sigma_cond_mean", "S0_mean",
                       "sigma_y_mean", "coverage")]
print(utils::head(scale_tab[order(-scale_tab$S0_mean), ], 10),
      row.names = FALSE, digits = 3)

cat("\nWritten to ", out_dir, "\n", sep = "")
