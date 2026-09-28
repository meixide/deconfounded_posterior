#!/usr/bin/env Rscript
## 10_environment_shift.R -------------------------------------------------
##
## How much do the environments actually differ?
##
## This is the diagnostic that decides how the case-study results should be
## read, and it was missing.  The correction term the method adds to a pooled
## regression is `K' Sigma_0^{-1} (X_0 - mu_0)`.  It can only do work to the
## extent that the target domain's covariate distribution differs from the
## training pool's.  Where environments are nearly exchangeable, `mu_0` sits on
## the pooled mean, the correction is near zero, and generative invariance
## necessarily reduces to ordinary least squares.  Finding no advantage in that
## regime is not evidence against the method; it is evidence that the dataset
## posed no problem for it to solve.
##
## Two numbers per dataset:
##
##   between/within  the ratio of between-environment variance in the covariate
##                   means to pooled within-environment variance, per covariate,
##                   on the standardised scale.  Small means the environments
##                   are near-replicates of each other.
##
##   Mahalanobis     the distance of each environment mean from the pooled mean
##                   in the metric of the within-environment covariance.
##                   Compare it with `p`, the expected squared distance of a
##                   single individual: an environment whose mean sits far
##                   closer to the centre than one person does is not a
##                   meaningfully shifted domain.
##
## Report this alongside the coverage table, so a reader can tell "the method
## did not help" from "there was nothing to help with".
##
## Usage:
##   Rscript scripts/10_environment_shift.R --data=PATH [--dataset=auto]

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
out_dir <- parse_flag(cli, "out", file.path(root, "results", "case_brfss"))

#' Between- versus within-environment variation in the covariate means.
#'
#' @param dat Output of `load_case_data()`.
#' @return A data frame with one row per covariate.
environment_shift <- function(dat) {
  xs <- scale(dat$x)
  z <- dat$z
  sw <- pooled_covariance(xs, z)
  mus <- t(vapply(split(seq_along(z), z),
                  function(i) colMeans(xs[i, , drop = FALSE]),
                  numeric(ncol(xs))))
  list(
    per_covariate = data.frame(
      covariate = dat$covariate_names,
      between_within = apply(mus, 2, stats::var) / diag(sw),
      stringsAsFactors = FALSE
    ),
    mahalanobis = stats::mahalanobis(mus, colMeans(mus), sw),
    p = ncol(dat$x)
  )
}

dat <- load_case_data(data_path, dataset = dataset)
sh <- environment_shift(dat)
tab <- sh$per_covariate[order(-sh$per_covariate$between_within), ]

cat(sprintf("\nN %d | p %d | E %d\n", length(dat$y), ncol(dat$x),
            length(unique(dat$z))))
cat("\n=== Between-environment variance as a fraction of within ===\n")
cat("Standardised scale, so the within-environment variance is about 1.\n\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  %-26s %.4f\n", tab$covariate[i], tab$between_within[i]))
}
cat(sprintf("\n  mean %.4f   max %.4f\n", mean(tab$between_within),
            max(tab$between_within)))

md <- sh$mahalanobis
cat("\n=== Mahalanobis distance of environment means from the pooled mean ===\n")
cat(sprintf("  median %.3f   max %.3f (%s)\n", stats::median(md), max(md),
            names(which.max(md))))
cat(sprintf("  expected squared distance of one individual: p = %d\n", sh$p))
cat(sprintf("  most extreme environment sits at %.1f%% of that\n",
            100 * max(md) / sh$p))

cat("\nRead the coverage table in this light. Where these numbers are small the\n")
cat("environments are near-replicates, the correction term has almost nothing\n")
cat("to correct, and the method is expected to coincide with pooled OLS.\n")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bgi_write_csv(tab, file.path(out_dir, "environment_shift.csv"))
bgi_write_csv(
  data.frame(environment = names(md), mahalanobis = as.numeric(md),
             stringsAsFactors = FALSE),
  file.path(out_dir, "environment_shift_mahalanobis.csv"))
cat("\nWritten to ", out_dir, "\n", sep = "")
