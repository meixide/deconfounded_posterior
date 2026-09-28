#!/usr/bin/env Rscript
## 20_aggregate_dimension_sweep.R -----------------------------------------
##
## Collect the per-task CSVs written by 19_sim_dimension_sweep.R and build
## Table 1 (`table:AEC`, Section 3.1.1).
##
## Two tables are written:
##   dimension_sweep.csv        cell-level averages, all columns
##   dimension_sweep_table1.txt the `n` x `p` grid in the "OLS/ours" layout
##                              the manuscript prints, ready to paste
##
## A third block is printed to stdout: coverage at `S_0` against coverage at
## the training residual scale, on the same fits.  That contrast is why the
## revised table differs from the submitted one, and it is the quantity
## Referee 2's minor comment (2) asks about.
##
## Replications whose fit failed its diagnostics (any divergent transition,
## or max Rhat above `--rhat-max`) are excluded from the averages and counted
## separately, on the same principle as 02_aggregate_support_recovery.R: a
## number that only looks right after silently dropping bad fits is not a
## number worth reporting.
##
## Usage:
##   Rscript scripts/20_aggregate_dimension_sweep.R [--in=DIR] [--out=DIR]

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
in_dir <- parse_flag(cli, "in", file.path(root, "results", "dimension_sweep"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "summaries"))
rhat_max <- as.numeric(parse_flag(cli, "rhat-max", "1.01"))

files <- list.files(in_dir, pattern = "^task_[0-9]+\\.csv$", full.names = TRUE)
if (length(files) == 0L) {
  stop("No task CSVs found in ", in_dir,
       ". Run scripts/19_sim_dimension_sweep.R first.")
}
## Task CSVs are written one per replication, possibly by different versions
## of 19_sim_dimension_sweep.R: a partial rerun after a column was added leaves
## a directory holding both schemas.  Aligning on the union of column names,
## rather than requiring every file to agree, means a rerun never invalidates
## the replications it did not touch.  A column a file predates reads as NA,
## which is the truth about it.
parts <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
all_names <- unique(unlist(lapply(parts, names)))
raw <- do.call(rbind, lapply(parts, function(d) {
  for (nm in setdiff(all_names, names(d))) d[[nm]] <- NA
  d[, all_names, drop = FALSE]
}))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Diagnostic gate ---------------------------------------------------

ok <- raw$divergences == 0 & raw$max_rhat <= rhat_max
dropped <- aggregate(cbind(n_total = rep(1, nrow(raw)),
                           n_dropped = as.integer(!ok)),
                     by = list(p = raw$p, n_e = raw$n_e), FUN = sum)
dat <- raw[ok, , drop = FALSE]
if (nrow(dat) == 0L) {
  stop("Every replication failed the diagnostic gate at rhat-max = ", rhat_max)
}

## ---- Cell averages -----------------------------------------------------

## Sampler settings are recorded per replication and averaged here only so
## that a cell using more than one setting is visible rather than silent.
## They are coordinates and step sizes, not model or data choices: the target
## posterior is the same whatever they are.
num_cols <- c("cov_ours", "score_ours", "width_ours", "rmse_ours",
              "cov_s0_scale", "cov_train_scale",
              "score_s0_scale", "score_train_scale",
              "cov_ols", "score_ols", "width_ols", "bgi_seconds")
settings_cols <- intersect(c("ncp", "adapt_delta"), names(dat))
num_cols <- c(num_cols, settings_cols)
## na.rm so that a cell mixing replications written before and after a column
## existed still reports the values that were recorded.
cells <- aggregate(dat[, num_cols],
                   by = list(p = dat$p, n_e = dat$n_e),
                   FUN = function(x) mean(x, na.rm = TRUE))
cells <- merge(cells, dropped, by = c("p", "n_e"))
cells <- cells[order(cells$p, cells$n_e), ]

bgi_write_csv(cells, file.path(out_dir, "dimension_sweep.csv"))

## ---- Table 1 in the manuscript's layout --------------------------------
## The manuscript prints one row per `n`, one column per `p`, each cell being
## "OLS/ours" to two decimals with the leading zero dropped.

fmt <- function(x) sub("^0", "", formatC(round(x, 2), format = "f", digits = 2))

p_vals <- sort(unique(cells$p))
n_vals <- sort(unique(cells$n_e))

lines <- c(
  sprintf("$n\\,\\backslash\\,p$ & %s \\\\ \\hline",
          paste(p_vals, collapse = " & ")))
for (n in n_vals) {
  entries <- vapply(p_vals, function(pp) {
    row <- cells[cells$p == pp & cells$n_e == n, , drop = FALSE]
    if (nrow(row) == 0L) "---" else
      sprintf("%s/%s", fmt(row$cov_ols), fmt(row$cov_ours))
  }, character(1))
  lines <- c(lines, sprintf("%-5d & %s \\\\", n, paste(entries, collapse = "   & ")))
}
writeLines(lines, file.path(out_dir, "dimension_sweep_table1.txt"))

cat("\nTable 1 body (OLS/ours), paste into jcgs.tex:\n\n")
cat(paste(lines, collapse = "\n"), "\n\n")

## ---- The correction, stated explicitly ---------------------------------

cat("Coverage at S_0 versus at the training residual scale (same fits):\n\n")
cmp <- cells[, c("p", "n_e", "cov_s0_scale", "cov_train_scale",
                 "score_s0_scale", "score_train_scale")]
cmp$cov_gain <- cmp$cov_s0_scale - cmp$cov_train_scale
print(format(cmp, digits = 3), row.names = FALSE)

cat(sprintf(
  "\nMean coverage: %.3f at S_0, %.3f at the training scale (%d cells).\n",
  mean(cmp$cov_s0_scale), mean(cmp$cov_train_scale), nrow(cmp)))

if (length(settings_cols) > 0L) {
  cat("\nSampler settings per cell (a non-integer ncp means the cell mixes",
      "\nreplications run under both parameterisations; NaN means none of the",
      "\nsurviving replications recorded it):\n\n")
  print(format(cells[, c("p", "n_e", settings_cols)], digits = 3),
        row.names = FALSE)
  n_unrecorded <- sum(!stats::complete.cases(dat[, settings_cols, drop = FALSE]))
  if (n_unrecorded > 0L) {
    cat(sprintf(
      "\n%d of %d surviving replications predate these columns.\n",
      n_unrecorded, nrow(dat)))
  }
}

if (any(dropped$n_dropped > 0L)) {
  cat("\nReplications excluded by the diagnostic gate:\n")
  bad <- dropped[dropped$n_dropped > 0L, ]
  print(bad, row.names = FALSE)
} else {
  cat("\nNo replication was excluded: every fit passed the diagnostic gate.\n")
}

message("\nWritten to ", out_dir)
