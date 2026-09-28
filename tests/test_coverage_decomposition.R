#!/usr/bin/env Rscript
## test_coverage_decomposition.R ------------------------------------------
##
## Why do the `gamma` credible intervals not attain 0.95 in simulation?
##
## First, a distinction that has to be kept straight, because the two things
## are routinely confused and only one of them is being measured here.
##
##   Bayesian (conditional) coverage.  P(gamma in C(D) | D) = 0.95, where the
##   probability is over the *posterior* given the observed data.  This is true
##   by construction: C(D) is built as a 0.95 posterior interval.  It is exact,
##   involves no asymptotics, and cannot fail.  Nothing in this script measures
##   it.
##
##   Frequentist coverage at a fixed truth.  P_{gamma*}(gamma* in C(D)), where
##   gamma* is held fixed and the probability is over repeated draws of D.
##   This is what a simulation study measures, it is what Section 3.1.2 of the
##   manuscript reports as "empirical coverage", and there is **no theorem
##   saying it must equal 0.95**.
##
## The two agree when averaged over the prior: if gamma ~ pi and D | gamma ~ p,
## then the joint probability P(gamma in C(D)) is exactly 0.95.  Our simulation
## does not draw gamma* from the prior — it fixes gamma*_j at +/-1 for parents
## and exactly 0 for nulls — so that guarantee does not apply.  At a fixed
## truth, agreement holds only asymptotically, under Bernstein-von Mises.
##
## So a value below 0.95 is not a contradiction.  The question is *why*, and
## there are only two mechanisms.  Writing g_hat for the posterior mean, s for
## the posterior sd, and S for the sampling sd of g_hat across replications:
##
##   bias           |E[g_hat] - gamma*| > 0
##   under-dispersion   s < S, i.e. the posterior is narrower than the actual
##                      sampling variability of its own centre
##
## This script estimates both, separately for true parents and true nulls,
## because shrinkage predicts they behave differently: a prior centred at zero
## pulls the parents (at +/-1) away from the truth and pulls the nulls (at 0)
## towards it.
##
## Usage:
##   Rscript tests/test_coverage_decomposition.R [--reps=30] [--confounding=2]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
n_reps <- as.integer(parse_flag(cli, "reps", "30"))
confounding <- as.numeric(parse_flag(cli, "confounding", "2"))
n_e <- as.integer(parse_flag(cli, "n-e", "200"))
n_env <- as.integer(parse_flag(cli, "n-env", "7"))
p <- as.integer(parse_flag(cli, "p", "6"))
s0 <- as.integer(parse_flag(cli, "s0", "3"))

model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))

post_mean <- matrix(NA_real_, n_reps, p)
post_sd <- matrix(NA_real_, n_reps, p)
truth_mat <- matrix(NA_real_, n_reps, p)
covered <- matrix(NA, n_reps, p)

for (r in seq_len(n_reps)) {
  bgi_set_seed(555000, r)
  dat <- simulate_gi_data(n_e = n_e, p = p, s0 = s0, n_env = n_env, q = 3,
                          confounding = confounding, n0 = 200)
  fit <- tryCatch(
    fit_bgi(dat$x, dat$y, dat$z, dat$x0, model = model,
            chains = 2, iter = 2000, seed = r),
    error = function(e) NULL)
  if (is.null(fit)) next

  g <- fit$draws$gamma
  post_mean[r, ] <- colMeans(g)
  post_sd[r, ] <- apply(g, 2, stats::sd)
  truth_mat[r, ] <- dat$truth$gamma
  lo <- apply(g, 2, stats::quantile, 0.025)
  hi <- apply(g, 2, stats::quantile, 0.975)
  covered[r, ] <- dat$truth$gamma >= lo & dat$truth$gamma <= hi
  if (r %% 5 == 0) message("  rep ", r, " of ", n_reps)
}

ok <- stats::complete.cases(post_mean)
post_mean <- post_mean[ok, , drop = FALSE]
post_sd <- post_sd[ok, , drop = FALSE]
truth_mat <- truth_mat[ok, , drop = FALSE]
covered <- covered[ok, , drop = FALSE]
n_ok <- nrow(post_mean)

## Pool coordinates by whether the truth is zero.  The identity of which
## coordinate is a parent changes between replications, so this is done
## element-wise rather than by column.
is_parent <- truth_mat != 0

summarise <- function(sel, label) {
  err <- (post_mean - truth_mat)[sel]
  s <- post_sd[sel]
  cov <- mean(covered[sel])
  ## Sampling sd of the posterior mean, computed after removing the truth so
  ## that parents with opposite signs can be pooled.
  S <- stats::sd(err)
  bias <- mean(err)
  cat(sprintf(
    "%-8s n=%4d | coverage %.3f | bias %+.4f | mean posterior sd %.4f | sampling sd of centre %.4f | ratio S/s %.2f\n",
    label, sum(sel), cov, bias, mean(s), S, S / mean(s)))
  ## Coverage that bias alone would produce if the posterior sd were honest.
  z <- bias / mean(s)
  cov_bias_only <- stats::pnorm(1.96 - z) - stats::pnorm(-1.96 - z)
  ## Coverage that under-dispersion alone would produce.
  cov_disp_only <- stats::pnorm(1.96 * mean(s) / S) -
    stats::pnorm(-1.96 * mean(s) / S)
  cat(sprintf("           predicted by bias alone %.3f | by under-dispersion alone %.3f\n",
              cov_bias_only, cov_disp_only))
  invisible(NULL)
}

cat(sprintf(
  "\n=== gamma coverage decomposition, %d usable replications ===\n", n_ok))
cat(sprintf("p = %d, s0 = %d, E = %d, n_e = %d, confounding = %g\n\n",
            p, s0, n_env, n_e, confounding))
cat("Frequentist coverage at a FIXED truth. Bayesian conditional coverage is\n")
cat("0.95 by construction and is not what this measures.\n\n")
summarise(is_parent, "parents")
summarise(!is_parent, "nulls")
summarise(matrix(TRUE, nrow(post_mean), p), "all")

cat("\nReading: S/s > 1 means the posterior is narrower than the sampling\n")
cat("variability of its own centre, i.e. genuinely over-confident. S/s ~ 1\n")
cat("with low coverage means the interval width is right and the centre is\n")
cat("displaced, i.e. bias -- which for a shrinkage prior is expected at a\n")
cat("fixed non-null truth and is not a defect of the implementation.\n")

out <- data.frame(
  group = c("parents", "nulls"),
  coverage = c(mean(covered[is_parent]), mean(covered[!is_parent])),
  bias = c(mean((post_mean - truth_mat)[is_parent]),
           mean((post_mean - truth_mat)[!is_parent])),
  mean_post_sd = c(mean(post_sd[is_parent]), mean(post_sd[!is_parent])),
  sampling_sd = c(stats::sd((post_mean - truth_mat)[is_parent]),
                  stats::sd((post_mean - truth_mat)[!is_parent])),
  stringsAsFactors = FALSE
)
out$ratio_S_over_s <- out$sampling_sd / out$mean_post_sd
bgi_write_csv(out, file.path(root, "results", "summaries",
                             "coverage_decomposition.csv"))
cat("\nWritten to results/summaries/coverage_decomposition.csv\n")
