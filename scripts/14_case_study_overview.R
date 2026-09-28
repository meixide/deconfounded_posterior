#!/usr/bin/env Rscript
## 14_case_study_overview.R -----------------------------------------------
##
## One table across every case-study run, for Section 3.2 and the response
## letter.
##
## Referee 2's objection to the submitted analysis was that coverage came from
## a single held-out domain, that the coverage events within it shared one
## posterior and so were strongly dependent, and that a single domain cannot
## distinguish a calibrated method from an uncalibrated one.  The reply is a
## design, not a number: hold out every environment in turn, make the
## environment the unit of replication, and take the standard error from the
## between-environment spread.  This script collects the result of that design
## over every dataset it has been run on, so the answer rests on several
## independent datasets rather than on one.
##
## It reports, per run:
##
##   convergence first, and unconditionally.  A coverage figure from chains
##   that did not mix is worse than no figure, so the diagnostics column is
##   printed next to the estimate rather than in an appendix.  Runs where a
##   large share of folds fail are shown with their coverage suppressed.
##
##   coverage with a between-environment standard error, and the minimum and
##   maximum over held-out domains.  The range is the part that answers the
##   referee: it shows directly what a single arbitrarily chosen domain could
##   have reported.
##
##   the same quantities for OLS, pooled GI and IV on identical folds.
##
## Usage:
##   Rscript scripts/14_case_study_overview.R [--min-pass=0.8] [--out=PATH]

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
## Below this pass rate a run is reported as not usable rather than summarised.
min_pass <- as.numeric(parse_flag(cli, "min-pass", "0.8"))
out_path <- parse_flag(cli, "out",
                       file.path(root, "results", "case_study_overview.csv"))

res_dir <- file.path(root, "results")
runs <- list.files(res_dir, "^case_", full.names = TRUE)
runs <- runs[dir.exists(file.path(runs, "folds"))]
if (length(runs) == 0L) {
  stop("No case-study runs with a folds/ directory under ", res_dir,
       call. = FALSE)
}

#' Split a results directory name into dataset and model.
#'
#' Directories are `case_<dataset>` for the default model and
#' `case_<dataset>_<model>` otherwise.
split_run <- function(nm) {
  body <- sub("^case_", "", nm)
  for (m in c("gi_hd_slope", "gi_hd_fullcov", "gi_hd")) {
    if (grepl(paste0("_", m, "$"), body)) {
      return(list(dataset = sub(paste0("_", m, "$"), "", body), model = m))
    }
  }
  list(dataset = body, model = "gi_hd")
}

rows <- list()
for (d in sort(runs)) {
  fs <- list.files(file.path(d, "folds"), "loeo_folds.csv",
                   recursive = TRUE, full.names = TRUE)
  if (length(fs) == 0L) next
  f <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
  who <- split_run(basename(d))

  pass <- !is.na(f$max_rhat) & f$max_rhat <= 1.01 &
    (is.null(f$min_ess) | f$min_ess >= 400)
  pass_rate <- mean(pass)
  usable <- pass_rate >= min_pass

  s <- if (usable) tryCatch(loeo_summary(f), error = function(e) NULL) else NULL
  pick <- function(m, col) {
    if (is.null(s) || !m %in% s$method) NA_real_ else s[[col]][s$method == m]
  }

  rows[[length(rows) + 1L]] <- data.frame(
    dataset = who$dataset,
    model = who$model,
    folds = nrow(f),
    p = f$p[1],
    E = f$n_train_env[1] + 1L,
    max_rhat = max(f$max_rhat, na.rm = TRUE),
    min_ess = min(f$min_ess, na.rm = TRUE),
    divergent = sum(f$divergent),
    pass_rate = pass_rate,
    usable = usable,
    bgi_coverage = pick("bgi", "coverage"),
    bgi_coverage_se = pick("bgi", "coverage_se"),
    bgi_coverage_min = pick("bgi", "coverage_min"),
    bgi_coverage_max = pick("bgi", "coverage_max"),
    bgi_interval_score = pick("bgi", "interval_score"),
    bgi_rmse = pick("bgi", "rmse"),
    ols_coverage = pick("ols", "coverage"),
    ols_interval_score = pick("ols", "interval_score"),
    ols_rmse = pick("ols", "rmse"),
    iv_coverage = pick("iv", "coverage"),
    iv_interval_score = pick("iv", "interval_score"),
    stringsAsFactors = FALSE
  )
}

tab <- do.call(rbind, rows)
tab <- tab[order(tab$dataset, tab$model), ]
bgi_write_csv(tab, out_path)

cat("\n=== Convergence, every run ===\n")
cat("A coverage estimate from chains that did not mix is worse than none, so\n")
cat("this is printed first and runs below the pass threshold are suppressed.\n\n")
print(tab[, c("dataset", "model", "folds", "p", "E", "max_rhat", "min_ess",
              "divergent", "pass_rate", "usable")],
      row.names = FALSE, digits = 4)

ok <- tab[tab$usable, , drop = FALSE]
if (nrow(ok) > 0L) {
  cat("\n=== Leave-one-environment-out coverage, usable runs ===\n")
  cat("Standard errors are across held-out environments, never across\n")
  cat("individuals: within a domain every interval shares one posterior.\n")
  cat("The min-max range is what a single arbitrarily chosen domain could\n")
  cat("have reported, which is the referee's point made quantitative.\n\n")
  print(ok[, c("dataset", "model", "folds", "bgi_coverage", "bgi_coverage_se",
               "bgi_coverage_min", "bgi_coverage_max", "ols_coverage",
               "iv_coverage")],
        row.names = FALSE, digits = 4)

  cat("\n=== Interval score and RMSE, usable runs ===\n")
  cat("Coverage alone is not evidence of a good interval: a wide enough\n")
  cat("interval always covers, and IV demonstrates exactly that below.\n\n")
  print(ok[, c("dataset", "model", "bgi_interval_score", "ols_interval_score",
               "iv_interval_score", "bgi_rmse", "ols_rmse")],
        row.names = FALSE, digits = 4)
}

bad <- tab[!tab$usable, , drop = FALSE]
if (nrow(bad) > 0L) {
  cat("\n=== Suppressed: fewer than ", round(100 * min_pass),
      "% of folds passed the diagnostics ===\n", sep = "")
  print(bad[, c("dataset", "model", "folds", "pass_rate", "max_rhat",
                "min_ess")],
        row.names = FALSE, digits = 4)
}

cat("\nWritten to ", out_path, "\n", sep = "")
