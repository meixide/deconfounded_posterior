#!/usr/bin/env Rscript
## test_fast_vs_reference.R -----------------------------------------------
##
## `stan/gi_hd.stan` reduces the training likelihood to per-environment
## sufficient statistics, which makes one gradient evaluation independent of
## the sample size but also makes the code much less obviously correct.  This
## test fits the fast model and the transparent per-observation reference
## implementation to the same data, from the same seed, and checks that the
## posteriors agree to within Monte Carlo error.
##
## The comparison is on the standardised scale used inside Stan, so both fits
## receive byte-identical inputs.  `lp__` is deliberately not compared: the
## two log-posteriors differ by an additive constant that depends on the data
## but not on the parameters.
##
## This test takes about forty-five minutes, almost all of it the reference
## implementation.  tests/test_likelihood_identity.R settles the same
## question exactly and in seconds, by checking that the two log-posteriors
## differ by a constant across the parameter space, and is the one the
## ten-minute check in run_checks.sh runs.  What this test adds is that it
## exercises sampling and the generated quantities block, which the identity
## check does not reach; run it before a release, not on every checkout.
##
## Usage:  Rscript tests/test_fast_vs_reference.R [--iter=4000]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

fast <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                       cache_dir = file.path(root, "results", "compiled"))
ref <- load_bgi_model(file.path(root, "stan", "gi_hd_reference.stan"),
                      cache_dir = file.path(root, "results", "compiled"))

bgi_set_seed(4242, 1)
dat <- simulate_gi_data(n_e = 120, p = 3, s0 = 2, n_env = 4, q = 2,
                        confounding = 1.5, n0 = 100)

prep <- prepare_bgi_data(dat$x, dat$y, dat$z, dat$x0, cov_method = "pooled")
sd_fast <- prep$stan_data

## The reference model takes the raw standardised data rather than the
## sufficient statistics; everything else is shared verbatim.
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
  ## Required by the reference model's data block; it is part of the
  ## shared prior parameterisation and must match the fast model's value,
  ## or the two posteriors are not comparable in the first place.
  ncp = sd_fast$ncp,
  N0 = sd_fast$N0, X0 = x0s, mu0 = sd_fast$mu0
)

control <- list(adapt_delta = 0.95, max_treedepth = 12)
chains <- 4
iter <- 4000

t_fast <- system.time(
  fit_fast <- rstan::sampling(fast, data = sd_fast, chains = chains,
                              iter = iter, seed = 11, refresh = 0,
                              control = control)
)[["elapsed"]]
t_ref <- system.time(
  fit_ref <- rstan::sampling(ref, data = sd_ref, chains = chains,
                             iter = iter, seed = 11, refresh = 0,
                             control = control)
)[["elapsed"]]

pars <- c("alpha", "gamma", "K", "sigma_y", "sigma_cond", "S0")
s_fast <- rstan::summary(fit_fast, pars = pars)$summary
s_ref <- rstan::summary(fit_ref, pars = pars)$summary
common <- intersect(rownames(s_fast), rownames(s_ref))

## Difference in posterior means, expressed in units of the Monte Carlo
## standard error of that difference.  Values well below 4 mean the two
## samplers are exploring the same distribution.
mean_diff <- s_fast[common, "mean"] - s_ref[common, "mean"]
mcse <- sqrt(s_fast[common, "se_mean"]^2 + s_ref[common, "se_mean"]^2)
z <- mean_diff / mcse

## Posterior standard deviations, compared on the same footing as the means.
## A posterior sd is itself an estimate, with Monte Carlo standard error
## approximately sd / sqrt(2 (ESS - 1)).  Judging the two sds against a fixed
## relative tolerance ignores that, and so reports "disagreement" whenever the
## reference has not been run long enough -- which, since it mixes about forty
## times more slowly than the fast model, is the usual case at any practical
## number of iterations.
sd_ratio <- s_fast[common, "sd"] / s_ref[common, "sd"]
mcse_sd_fast <- s_fast[common, "sd"] / sqrt(2 * (s_fast[common, "n_eff"] - 1))
mcse_sd_ref <- s_ref[common, "sd"] / sqrt(2 * (s_ref[common, "n_eff"] - 1))
z_sd <- (s_fast[common, "sd"] - s_ref[common, "sd"]) /
  sqrt(mcse_sd_fast^2 + mcse_sd_ref^2)

report <- data.frame(
  parameter = common,
  mean_fast = round(s_fast[common, "mean"], 4),
  mean_ref = round(s_ref[common, "mean"], 4),
  z_mean_diff = round(z, 2),
  sd_ratio = round(sd_ratio, 3),
  z_sd_diff = round(z_sd, 2),
  ess_fast = round(s_fast[common, "n_eff"]),
  ess_ref = round(s_ref[common, "n_eff"]),
  row.names = NULL
)
print(report, row.names = FALSE)

cat(sprintf("\nRuntime: fast %.1f s, reference %.1f s (speed-up %.1fx)\n",
            t_fast, t_ref, t_ref / t_fast))

ok_mean <- max(abs(z)) < 4
ok_sd <- max(abs(z_sd)) < 4
min_ess_ref <- min(s_ref[common, "n_eff"])

cat(sprintf("max |z| on posterior means : %.2f  (tolerance 4)\n", max(abs(z))))
cat(sprintf("max |z| on posterior sds   : %.2f  (tolerance 4)\n", max(abs(z_sd))))
cat(sprintf("max relative sd discrepancy: %.3f  (reported, not a criterion)\n",
            max(abs(sd_ratio - 1))))
cat(sprintf("smallest reference ESS     : %.0f of %d draws\n",
            min_ess_ref, chains * iter / 2))

## An ESS this low means the reference's own posterior sds are too noisy for
## the comparison to have power.  Say so, rather than reporting the noise as
## a disagreement between the two models.
if (min_ess_ref < 100) {
  cat("\nThe reference implementation did not mix well enough at",
      sprintf("%d iterations\n", iter),
      "for this comparison to be informative. Raise --iter, or rely on\n",
      "tests/test_likelihood_identity.R, which settles the same question\n",
      "exactly and in seconds.\n")
}

if (ok_mean && ok_sd) {
  cat("\nPASS: the fast model reproduces the reference posterior.\n")
} else {
  cat("\nFAIL: the fast and reference posteriors disagree beyond Monte Carlo",
      "error.\n      Check tests/test_likelihood_identity.R first: if that",
      "passes, the two\n      models are the same and the disagreement here",
      "is a sampling problem.\n")
  quit(status = 1)
}
