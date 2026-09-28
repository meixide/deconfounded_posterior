#!/usr/bin/env Rscript
## check_slope_recovery.R -------------------------------------------------
##
## Recovery check for the varying-slope model under HETEROGENEOUS Sigma_e --
## the case the common-b model could not handle.
##
## Prints the true gamma alongside the posterior mean, the RMSE between them,
## and the sampler diagnostics (divergences, max Rhat, min bulk ESS) plus
## predictive coverage.
##
## Pass --force after editing stan/gi_hd_slope.stan; without it the cached
## binary under results/compiled/ is reused and source edits are NOT picked up.
##
## Usage:
##   BGI_ROOT=$PWD Rscript scripts/check_slope_recovery.R [--force]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

force <- "--force" %in% commandArgs(trailingOnly = TRUE)

stan_file <- file.path(root, "stan", "gi_hd_slope.stan")
cache_dir <- file.path(root, "results", "compiled")

compile_bgi_model(stan_file, cache_dir = cache_dir, force = force)
message("COMPILED")

m <- load_bgi_model(stan_file, cache_dir = cache_dir)

bgi_set_seed(9, 1)
d <- simulate_gi_data(n_e = 300, p = 4, s0 = 2, n_env = 9, q = 3,
                      confounding = 2, n0 = 300)
f <- fit_bgi(d$x, d$y, d$z, d$x0, model = m,
             chains = 2, iter = 1500, seed = 3, ncp = 0)

gamma_hat <- colMeans(f$draws$gamma)
rmse <- sqrt(mean((gamma_hat - d$truth$gamma)^2))

cat(sprintf("truth gamma : %s\n",
            paste(sprintf("%6.3f", d$truth$gamma), collapse = " ")))
cat(sprintf("slope model : %s | rmse %.4f\n",
            paste(sprintf("%6.3f", gamma_hat), collapse = " "), rmse))
cat(sprintf("div %d | Rhat %.3f | minESS %.0f | pred cov %.3f\n",
            f$diagnostics$n_divergent, f$diagnostics$max_rhat,
            f$diagnostics$min_ess_bulk,
            predictive_metrics(f$draws$y0_pred, d$y0)$coverage))
