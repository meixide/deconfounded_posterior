#!/usr/bin/env Rscript
## test_gamma_coverage_sources.R ------------------------------------------
##
## Where does the under-coverage of `gamma` come from?
##
## In the `conf2_strong` scenario of the support-recovery study, credible
## intervals for `gamma` cover at about 0.70 against a nominal 0.95, and the
## selection rules inherit the over-confidence as false discoveries.  This
## script isolates the candidate causes by varying one modelling choice at a
## time on the same replications:
##
##   base          the default: covariances plugged in, mu_e given the
##                 hierarchical prior N(hmu, Sigma_mu) of Section 2
##   flat_mu       sd_mu_scale raised so the mu_e prior is effectively flat.
##                 Isolates attenuation from shrinking the environment means
##                 towards the pooled mean: that shrinkage compresses exactly
##                 the between-environment spread that identifies
##                 Sigma^{-1} K, and gamma is the within-environment slope
##                 minus Sigma^{-1} K, so any attenuation of the latter lands
##                 on gamma as bias.
##   flat_ridge    a_tau raised so the ridge prior on (alpha, gamma, K) is
##                 effectively flat.  Isolates shrinkage of gamma itself.
##   fullcov       covariances inferred rather than plugged in.  Isolates the
##                 uncertainty discarded by conditioning on Sigma_hat.
##
## Coverage well below nominal indicates bias comparable to the posterior sd,
## so the variant that restores coverage identifies the source.
##
## Usage:
##   Rscript tests/test_gamma_coverage_sources.R [--reps=6]

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
n_reps <- as.integer(parse_flag(cli, "reps", "6"))
n_iter <- as.integer(parse_flag(cli, "iter", "1500"))

m_plugin <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                           cache_dir = file.path(root, "results", "compiled"))
m_full <- load_bgi_model(file.path(root, "stan", "gi_hd_fullcov.stan"),
                         cache_dir = file.path(root, "results", "compiled"))

variants <- list(
  base       = list(model = m_plugin, args = list()),
  flat_mu    = list(model = m_plugin, args = list(sd_mu_scale = 100)),
  flat_ridge = list(model = m_plugin, args = list(a_tau = 50, b_tau = 0.5)),
  fullcov    = list(model = m_full, args = list())
)

rows <- list()
for (r in seq_len(n_reps)) {
  bgi_set_seed(20260727, 100 + r)
  dat <- simulate_gi_data(n_e = 200, p = 6, s0 = 3, n_env = 7, q = 3,
                          confounding = 2, n0 = 500)
  ## Frequentist reference on the same data: it uses the sample environment
  ## means directly and applies no shrinkage to anything.
  pgi <- fit_pooled_gi(dat$x, dat$y, dat$z)
  pgi_cov <- mean(abs(pgi$gamma - dat$truth$gamma) <=
                    stats::qnorm(0.975) * pgi$se)

  for (nm in names(variants)) {
    v <- variants[[nm]]
    fit <- tryCatch(
      do.call(fit_bgi, c(list(x = dat$x, y = dat$y, z = dat$z, x0 = dat$x0,
                              model = v$model, chains = 2, iter = n_iter,
                              seed = r), v$args)),
      error = function(e) {
        message("  failed (", nm, ", rep ", r, "): ", conditionMessage(e))
        NULL
      })
    if (is.null(fit)) next

    pc <- parameter_coverage(fit$draws$gamma, dat$truth$gamma)
    bias <- mean(colMeans(fit$draws$gamma) - dat$truth$gamma)
    post_sd <- mean(apply(fit$draws$gamma, 2, stats::sd))
    sel <- select_parents(fit$draws$gamma, alpha = 0.05, rule = "sign")
    sm <- support_metrics(sel$selected, dat$truth$parents, 6)

    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, variant = nm,
      gamma_coverage = pc$gamma_coverage,
      gamma_rmse = pc$gamma_rmse,
      mean_bias = bias,
      mean_post_sd = post_sd,
      bias_over_sd = abs(bias) / post_sd,
      fdp = sm$fdp, tpr = sm$tpr,
      divergent = fit$diagnostics$n_divergent,
      pooled_gi_coverage = pgi_cov,
      stringsAsFactors = FALSE
    )
    message(sprintf("  rep %d %-11s coverage %.3f  |bias|/sd %.2f  fdp %.2f",
                    r, nm, pc$gamma_coverage, abs(bias) / post_sd, sm$fdp))
  }
}

res <- do.call(rbind, rows)
cat("\n=== Means by variant (nominal gamma coverage 0.95) ===\n")
summ <- aggregate(
  cbind(gamma_coverage, gamma_rmse, mean_post_sd, bias_over_sd, fdp, tpr,
        divergent, pooled_gi_coverage) ~ variant,
  data = res, FUN = mean)
print(summ, row.names = FALSE, digits = 3)

bgi_write_csv(res, file.path(root, "results", "summaries",
                             "gamma_coverage_sources.csv"))
cat("\nWritten to results/summaries/gamma_coverage_sources.csv\n")
