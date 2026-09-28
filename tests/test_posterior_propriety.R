#!/usr/bin/env Rscript
## test_posterior_propriety.R ---------------------------------------------
##
## Is the posterior a probability distribution?
##
## gi_hd_slope.stan anchors the residual variance at the largest Mahalanobis
## term, sigma_y^2 = v_raw + max(max_e b_e' Sigma_e b_e, K_bar' Sigma_0^{-1}
## K_bar), and gives v_raw the Jeffreys prior p(v_raw) propto 1/v_raw, which is
## improper.  Whether the posterior is nonetheless proper turns on what the
## likelihood does as v_raw -> 0, and that depends on which term attains the
## maximum.
##
## If a training environment attains it, the conditional variance of that
## environment is exactly v_raw, so its contribution behaves like
## v_raw^{-n/2} exp(-SS / 2 v_raw), which vanishes faster than 1/v_raw
## diverges.  The posterior is proper.
##
## If the target attains it instead, every training conditional variance stays
## bounded away from zero as v_raw -> 0, the likelihood tends to a positive
## constant, and the remaining integral is int_0 dv/v.  The posterior is
## improper.
##
## The check is direct.  On the unconstrained scale v_raw = exp(u), the
## Jacobian log v exactly cancels the Jeffreys -log v, so log_prob(u) tends to
## the log-likelihood as u -> -infinity.  A limit that is finite means a
## density constant over an infinite range, hence not integrable.  Flat output
## in the first block below is the impropriety, not a numerical artefact.
##
## A proper inverse-gamma prior on v_raw, which fit_bgi() exposes as
## v_prior_shape and v_prior_rate, removes the problem; the default of zero
## selects Jeffreys.
##
## Usage:  Rscript tests/test_posterior_propriety.R

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()
suppressMessages(library(rstan))
m <- load_bgi_model(file.path(root, "stan", "gi_hd_slope.stan"),
                    cache_dir = file.path(root, "results", "compiled"))

## v_prior_shape is passed explicitly.  It used to be left to the default,
## which was zero; when the default became the proper inverse-gamma this test
## would have gone on printing a clean sweep and reporting nothing wrong, having
## quietly stopped exercising the prior it exists to indict.  A test of a default
## must not depend on that default.
make_fit <- function(shrink, shape = 0, rate = 0) {
  bgi_set_seed(11, 1)
  d <- simulate_gi_data(n_e=150, p=2, s0=2, n_env=3, q=2, confounding=1.5, n0=150)
  ctr <- colMeans(d$x0)
  d$x0 <- sweep(sweep(d$x0, 2, ctr, "-") * shrink, 2, ctr, "+")
  prep <- prepare_bgi_data(d$x, d$y, d$z, d$x0, cov_method="pooled",
                           v_prior_shape = shape, v_prior_rate = rate)
  suppressWarnings(suppressMessages(
    sampling(m, data=prep$stan_data, chains=1, iter=2, refresh=0, seed=1)))
}

find_v_idx <- function(f) {
  n <- get_num_upars(f); set.seed(1); u <- rnorm(n, 0, 0.3)
  v0 <- constrain_pars(f, u)$v_raw
  for (j in seq_len(n)) {
    u2 <- u; u2[j] <- u2[j] + 1
    if (abs(constrain_pars(f, u2)$v_raw - v0) > 1e-8) return(j)
  }
  NA
}

sweep_v <- function(f, j, tag) {
  n <- get_num_upars(f); set.seed(1); u <- rnorm(n, 0, 0.3)
  cat("\n", tag, "\n", sep="")
  cat(sprintf("  %-12s %-14s %s\n", "v_raw", "log_prob(u)", "S0 vs sigma_cond"))
  for (uv in c(0, -3, -6, -10, -15, -20, -30, -45, -60, -90)) {
    u[j] <- uv
    lp <- log_prob(f, u, adjust_transform=TRUE, gradient=FALSE)
    cp <- tryCatch(constrain_pars(f, u), error=function(e) NULL)
    cat(sprintf("  %-12.3e %-14.3f %s\n", exp(uv), lp,
        if (is.null(cp)) "(gq overflow)" else
          sprintf("%.5f / %.5f", cp$S0, cp$sigma_cond)))
  }
}

f1 <- make_fit(0.15); j1 <- find_v_idx(f1)
f2 <- make_fit(1.00); j2 <- find_v_idx(f2)
cat("v_raw is unconstrained coordinate", j1, "and", j2, "\n")
sweep_v(f1, j1,
        "JEFFREYS, TARGET CONCENTRATED (shrink 0.15): target term binds, so\n  log_prob should FLATTEN -- a density constant over an infinite range")
sweep_v(f2, j2,
        "JEFFREYS, TARGET AS TRAINING (shrink 1.00): a training env binds, so\n  log_prob should FALL away without limit")

## The third block is the remedy, and the reason the default changed.  Same data
## as the first, same binding term, the only difference being the prior the
## calibration study of tests/test_sbc.R validates.  log_prob must now fall away
## even though the target term still binds: the inverse-gamma vanishes at the
## origin like v^{-(shape+1)} exp(-rate/v), which no bounded likelihood can
## offset.
f3 <- make_fit(0.15, shape = 3, rate = 2); j3 <- find_v_idx(f3)
sweep_v(f3, j3,
        "INVERSE-GAMMA(3, 2), TARGET CONCENTRATED: same binding term, and\n  log_prob should FALL away -- this is the default the package now uses")
