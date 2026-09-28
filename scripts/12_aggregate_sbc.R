#!/usr/bin/env Rscript
## 12_aggregate_sbc.R ------------------------------------------------------
##
## Collect the simulation-based calibration runs written by `slurm/12_sbc.sh`
## into one table, and draw the rank histograms.
##
## How to read the output.  SBC checks the one guarantee that is exact in
## finite samples: if theta* is drawn from the prior and the data from the
## likelihood, the rank of theta* among L posterior draws is uniform on
## {0, ..., L}.  So
##
##   uniformity_p > 0.05 and extreme_rank_ratio ~ 1   calibrated
##   extreme_rank_ratio > 1                           posterior too narrow
##   extreme_rank_ratio < 1                           posterior too wide
##
## `coverage_95` in this table is *prior-averaged* coverage, which is the
## quantity that must equal 0.95.  It is not the fixed-truth coverage reported
## by the simulation studies, and the two should not be compared: no theorem
## says fixed-truth coverage equals the nominal level, so a shortfall there is
## not by itself evidence of a defect.  This is the distinction that decides
## whether the residual gap in `results/slope_vs_k` is a bug or a property of
## shrinkage at that particular truth.
##
## Usage:
##   Rscript scripts/12_aggregate_sbc.R

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
# An aggregator reads CSVs; it neither fits nor simulates, so it does not
# need the fitting dependencies.  See `need_stan` in R/setup.R.
root <- bgi_bootstrap(need_stan = FALSE)

sum_dir <- file.path(root, "results", "summaries")
## `sbc_all.csv` is this script's own output, written into the same directory it
## reads.  Including it in the input made the second run of the script fail
## where the first had succeeded -- the combined file carries two columns the
## per-arm files do not, so rbind refused -- and the failure was silent to
## check_tables.sh, which discarded stderr and printed nothing under Table 3.
## A script must not read what it writes.
out_name <- "sbc_all.csv"
files <- setdiff(
  list.files(sum_dir, pattern = "^sbc_.*\\.csv$", full.names = TRUE),
  file.path(sum_dir, out_name))
if (length(files) == 0L) {
  stop("No per-arm sbc_*.csv in ", sum_dir, ". Has slurm/12_sbc.sh run?",
       call. = FALSE)
}

d <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
d <- d[order(d$model, d$hetero, d$arm), ]

## Monte Carlo error on the coverage, so a shortfall can be judged rather than
## eyeballed.  Coordinates within a replication share a posterior, so this is
## optimistic; treat it as a lower bound on the uncertainty.
d$coverage_se <- sqrt(d$coverage_95 * (1 - d$coverage_95) /
                        (d$usable_reps * d$p))
d$verdict <- ifelse(
  d$uniformity_p > 0.05 & abs(d$extreme_rank_ratio - 1) < 0.5, "calibrated",
  ifelse(d$extreme_rank_ratio > 1, "too narrow", "too wide"))

cat("\n=== Simulation-based calibration ===\n\n")
print(d[, c("model", "hetero", "arm", "usable_reps", "coverage_95",
            "coverage_se", "uniformity_p", "extreme_rank_ratio", "verdict")],
      row.names = FALSE, digits = 3)

bgi_write_csv(d, file.path(sum_dir, out_name))

## ---- Rank histograms ----------------------------------------------------
## Greyscale, base_size 11, as the editor requires: black and white printing
## and no text below 10pt.
rank_files <- list.files(file.path(root, "results", "sbc"),
                         pattern = "^ranks_.*\\.csv$", full.names = TRUE)
if (length(rank_files) > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
  rk <- do.call(rbind, lapply(rank_files, utils::read.csv,
                              stringsAsFactors = FALSE))
  rk$panel <- sprintf("%s | %s | h = %s", rk$model, rk$arm, rk$hetero)
  n_draws <- max(rk$n_draws)
  ## The 99% band for a uniform histogram, so a departure can be judged against
  ## sampling noise rather than by eye.
  per_panel <- stats::aggregate(rank ~ panel, rk, length)
  expected <- mean(per_panel$rank) / (n_draws + 1)
  band <- stats::qbinom(c(0.005, 0.995), round(mean(per_panel$rank)),
                        1 / (n_draws + 1))

  pl <- ggplot2::ggplot(rk, ggplot2::aes(x = rank)) +
    ggplot2::geom_histogram(bins = n_draws + 1, fill = "grey55",
                            colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = expected, linetype = "dashed",
                        linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = band, linetype = "dotted",
                        linewidth = 0.3) +
    ggplot2::facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    ggplot2::labs(x = "rank of the true value among posterior draws",
                  y = "count") +
    ggplot2::theme_bw(base_size = 11)

  out_pdf <- file.path(root, "results", "figures", "sbc_ranks.pdf")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_pdf, pl, width = 7, height = 1.9 * ceiling(
    length(unique(rk$panel)) / 2), device = grDevices::cairo_pdf)
  cat("\nRank histograms written to results/figures/sbc_ranks.pdf\n")
}

cat("\nWritten to results/summaries/sbc_all.csv\n")
