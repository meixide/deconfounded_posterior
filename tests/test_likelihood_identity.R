#!/usr/bin/env Rscript
## test_likelihood_identity.R ---------------------------------------------
##
## `stan/gi_hd.stan` reduces the training likelihood to per-environment
## sufficient statistics, which makes one gradient evaluation independent of
## the sample size but also makes the code much less obviously correct.
## `stan/gi_hd_reference.stan` is the transparent per-observation
## implementation of the same model.
##
## This test checks the reduction *exactly*, without sampling anything.
##
## Why exactly rather than by comparing posteriors
## -----------------------------------------------
## The two models declare identical parameter blocks -- same names, types,
## order and dimensions -- so a point in the unconstrained space means the
## same thing to both.  If the reduction is correct then, for every such
## point,
##
##     log p_fast(theta | data) - log p_ref(theta | data) = constant,
##
## the constant being the part of the per-observation normalisation that the
## sufficient-statistic form drops.  It depends on the data and not on
## theta.  Evaluating both log-posteriors at a few dozen random points and
## checking that the difference has zero variance therefore settles the
## question outright, in seconds.
##
## Comparing fitted posteriors instead, as tests/test_fast_vs_reference.R
## does, answers the same question only up to Monte Carlo error, takes about
## forty-five minutes, and is the weaker evidence of the two: the reference
## implementation mixes roughly forty times more slowly, so at any practical
## number of iterations its posterior standard deviations are the noisier
## estimate and the comparison inherits that noise.  That test is still worth
## running before a release, and it checks the `generated quantities` block
## which this one does not reach.  This one is what belongs in a
## ten-minute check.
##
## Usage:  Rscript tests/test_likelihood_identity.R [--points=40]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

suppressMessages(library(rstan))

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
n_points <- as.integer(parse_flag(cli, "points", "40"))
## Zero variance up to double-precision accumulation over a few hundred
## terms.  The observed spread is ~1e-12 on data of this size; 1e-6 leaves
## several orders of magnitude of headroom while still failing loudly on any
## real discrepancy, which would show up at the scale of the log-likelihood.
tol <- as.numeric(parse_flag(cli, "tol", "1e-6"))

fast <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                       cache_dir = file.path(root, "results", "compiled"))
ref <- load_bgi_model(file.path(root, "stan", "gi_hd_reference.stan"),
                      cache_dir = file.path(root, "results", "compiled"))

## Same data construction as tests/test_fast_vs_reference.R, so that the two
## tests are demonstrably about the same comparison.
bgi_set_seed(4242, 1)
dat <- simulate_gi_data(n_e = 120, p = 3, s0 = 2, n_env = 4, q = 2,
                        confounding = 1.5, n0 = 100)
prep <- prepare_bgi_data(dat$x, dat$y, dat$z, dat$x0, cov_method = "pooled")
sd_fast <- prep$stan_data

xs <- sweep(sweep(dat$x, 2, prep$scaling$centre, "-"), 2,
            prep$scaling$scale, "/")
x0s <- sweep(sweep(dat$x0, 2, prep$scaling$centre, "-"), 2,
             prep$scaling$scale, "/")
sd_ref <- list(
  P = sd_fast$P, E = sd_fast$E, N = sd_fast$N,
  X = xs, Y = dat$y, Z = as.integer(as.factor(dat$z)),
  L_Sigma = sd_fast$L_Sigma, L_Sigma_bar = sd_fast$L_Sigma_bar,
  L_Sigma0 = sd_fast$L_Sigma0,
  hmu = sd_fast$hmu, eta_lkj = sd_fast$eta_lkj,
  sd_mu_scale = sd_fast$sd_mu_scale,
  a_tau = sd_fast$a_tau, b_tau = sd_fast$b_tau,
  ncp = sd_fast$ncp,
  N0 = sd_fast$N0, X0 = x0s, mu0 = sd_fast$mu0
)

## `log_prob` needs a stanfit to carry the model and data; two iterations is
## the cheapest way to obtain one and the draws are discarded.
skeleton <- function(model, data) {
  suppressWarnings(suppressMessages(
    rstan::sampling(model, data = data, chains = 1, iter = 2,
                    refresh = 0, seed = 1)))
}
f_fast <- skeleton(fast, sd_fast)
f_ref <- skeleton(ref, sd_ref)

n_upars <- rstan::get_num_upars(f_fast)
cat(sprintf("unconstrained dimension: fast %d, reference %d\n",
            n_upars, rstan::get_num_upars(f_ref)))
if (n_upars != rstan::get_num_upars(f_ref)) {
  cat("\nFAIL: the two models do not share a parameter space, so the",
      "reduction cannot be checked this way.\n")
  quit(status = 1)
}

set.seed(7)
diffs <- vapply(seq_len(n_points), function(i) {
  u <- stats::rnorm(n_upars, 0, 0.5)
  a <- rstan::log_prob(f_fast, u, adjust_transform = TRUE, gradient = FALSE)
  b <- rstan::log_prob(f_ref, u, adjust_transform = TRUE, gradient = FALSE)
  a - b
}, numeric(1))

finite <- is.finite(diffs)
if (sum(finite) < 2L) {
  cat("\nFAIL: fewer than two points gave a finite log-posterior.\n")
  quit(status = 1)
}
diffs <- diffs[finite]
spread <- max(diffs) - min(diffs)

cat(sprintf("points evaluated        : %d\n", length(diffs)))
cat(sprintf("lp_fast - lp_ref        : %.6f\n", mean(diffs)))
cat(sprintf("spread across points    : %.3e  (tolerance %.0e)\n", spread, tol))

if (spread < tol) {
  cat("\nPASS: the sufficient-statistic likelihood equals the",
      "per-observation one\n      up to an additive constant.\n")
  quit(status = 0)
}

cat("\nFAIL: the difference is not constant, so stan/gi_hd.stan is not the",
    "same\n      model as stan/gi_hd_reference.stan.\n")
quit(status = 1)
