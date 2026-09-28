## Needs the optional model gi_hd_fullcov.stan, which is not built by default:
##     Rscript scripts/00_compile_models.R --with-optional
## It requires a newer StanHeaders than rstan 2.32.3 ships with.
#!/usr/bin/env Rscript
## test_plugin_vs_fullcov.R -----------------------------------------------
##
## Does inferring the covariance matrices fix the under-coverage of `gamma`?
##
## `stan/gi_hd.stan` conditions on the plug-in `Sigma_hat` as if it were
## known.  Because `gamma` is the within-environment slope minus
## `Sigma^{-1} K`, and `||Sigma^{-1} K||` is large under strong confounding,
## that discards a first-order source of uncertainty: in the support-recovery
## simulation, credible intervals for `gamma` covered at 0.70 against a
## nominal 0.95 in the `conf2_strong` scenario, and the selection rules
## inherited the over-confidence as false discoveries.
##
## `stan/gi_hd_fullcov.stan` gives `Sigma_1, ..., Sigma_E` an inverse-Wishart
## hierarchy and `Sigma_0` its own prior.  This script fits both to the same
## replications of the scenario where the problem is worst, and reports
## coverage, RMSE, false discovery proportion and runtime.
##
## Usage:
##   Rscript tests/test_plugin_vs_fullcov.R [--reps=5] [--iter=1500]

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
n_reps <- as.integer(parse_flag(cli, "reps", "5"))
n_iter <- as.integer(parse_flag(cli, "iter", "1500"))

models <- list(
  plugin = load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                          cache_dir = file.path(root, "results", "compiled")),
  fullcov = load_bgi_model(file.path(root, "stan", "gi_hd_fullcov.stan"),
                           cache_dir = file.path(root, "results", "compiled"))
)

rows <- list()
for (r in seq_len(n_reps)) {
  bgi_set_seed(20260727, 100 + r)
  dat <- simulate_gi_data(n_e = 200, p = 6, s0 = 3, n_env = 7, q = 3,
                          confounding = 2, n0 = 500)

  for (nm in names(models)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      fit_bgi(dat$x, dat$y, dat$z, dat$x0, model = models[[nm]],
              chains = 2, iter = n_iter, seed = r),
      error = function(e) {
        message("  fit failed (", nm, ", rep ", r, "): ",
                conditionMessage(e))
        NULL
      })
    if (is.null(fit)) next

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    pc <- parameter_coverage(fit$draws$gamma, dat$truth$gamma)
    sel <- select_parents(fit$draws$gamma, alpha = 0.05, rule = "sign")
    sm <- support_metrics(sel$selected, dat$truth$parents, 6)

    rows[[length(rows) + 1L]] <- data.frame(
      rep = r,
      model = nm,
      gamma_coverage = pc$gamma_coverage,
      gamma_rmse = pc$gamma_rmse,
      tpr = sm$tpr,
      fdp = sm$fdp,
      pred_coverage = predictive_metrics(fit$draws$y0_pred, dat$y0)$coverage,
      divergent = fit$diagnostics$n_divergent,
      max_rhat = fit$diagnostics$max_rhat,
      seconds = elapsed,
      stringsAsFactors = FALSE
    )
    message(sprintf("  rep %d %-8s  gamma coverage %.3f  fdp %.3f  %.0f s",
                    r, nm, pc$gamma_coverage, sm$fdp, elapsed))
  }
}

res <- do.call(rbind, rows)
cat("\n=== Per replication ===\n")
print(res, row.names = FALSE, digits = 3)

cat("\n=== Means by model (nominal gamma coverage 0.95) ===\n")
summ <- aggregate(
  cbind(gamma_coverage, gamma_rmse, tpr, fdp, pred_coverage, divergent,
        seconds) ~ model,
  data = res, FUN = mean)
print(summ, row.names = FALSE, digits = 3)

bgi_write_csv(res, file.path(root, "results", "summaries",
                             "plugin_vs_fullcov.csv"))
cat("\nWritten to results/summaries/plugin_vs_fullcov.csv\n")
