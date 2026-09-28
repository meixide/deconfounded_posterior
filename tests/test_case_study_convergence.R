#!/usr/bin/env Rscript
## test_case_study_convergence.R ------------------------------------------
##
## Does the model converge on the case-study data, and under which settings?
##
## `tests/test_case_study_timing.R` established that the centred
## parameterisation is about 13x faster per iteration on this design, which is
## the opposite of the simulation regime and is expected: non-centring helps
## when the likelihood is weak, and with 400,000 training rows it is anything
## but.  That test could not speak to *convergence*, because it ran only 100
## warmup iterations and Stan's mass-matrix adaptation needs roughly 150-200 to
## complete its windows.
##
## This runs full-length chains and reports Rhat and ESS honestly.  It also
## varies `eta_lkj`, because with 51 environment means informing a 12 x 12
## correlation matrix the `Sigma_mu` correlations are weakly identified; a
## large `eta_lkj` shrinks them towards diagonal and removes 66 parameters.
##
## Usage:
##   Rscript tests/test_case_study_convergence.R [--iter=2000] [--chains=2]

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
n_iter <- as.integer(parse_flag(cli, "iter", "2000"))
n_chains <- as.integer(parse_flag(cli, "chains", "2"))
data_path <- parse_flag(cli, "data",
                        file.path(dirname(root), "old_code", "quiron",
                                  "quiron_final.csv"))
dataset <- parse_flag(cli, "dataset", "auto")

dat <- load_case_data(data_path, dataset = dataset)

## Default to the largest environment, which is the hardest fold for the
## `generated quantities` block and the one the submitted analysis used.
target <- parse_flag(cli, "target",
                     names(sort(table(dat$z), decreasing = TRUE))[1])
if (!target %in% dat$z) {
  stop("No environment called ", target, ". Available: ",
       paste(utils::head(sort(unique(dat$z)), 10), collapse = ", "), " ...",
       call. = FALSE)
}
## Selectable rather than hard-coded: the paper now recommends the slope
## parameterisation, so the convergence gate has to be readable for that model
## and not only for the covariance one it was originally written against. The
## slope model also keeps Sigma_e^{-1} out of the conditional mean, which is
## the likelier cure for the treedepth saturation seen with `gi_hd`.
model_name <- parse_flag(cli, "model", "gi_hd_slope")
model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))
cat(sprintf("model: %s\n", model_name))

is_target <- dat$z == target
keep <- !is_target
i0 <- utils::head(which(is_target), 500)
x0 <- dat$x[i0, , drop = FALSE]
y0 <- dat$y[i0]

cat(sprintf("\nTarget %s | N_train %d | E_train %d | p %d\n",
            target, sum(keep), length(unique(dat$z[keep])), ncol(dat$x)))
cat(sprintf("%d chains x %d iterations\n\n", n_chains, n_iter))

settings <- list(
  list(tag = "centred, eta_lkj = 2",   ncp = 0, eta = 2),
  list(tag = "centred, eta_lkj = 30",  ncp = 0, eta = 30),
  list(tag = "non-centred, eta_lkj=2", ncp = 1, eta = 2)
)

for (cfg in settings) {
  t0 <- Sys.time()
  fit <- tryCatch(
    fit_bgi(dat$x[keep, , drop = FALSE], dat$y[keep], dat$z[keep], x0,
            model = model, ncp = cfg$ncp, eta_lkj = cfg$eta,
            chains = n_chains, iter = n_iter, seed = 1, cores = n_chains),
    error = function(e) {
      cat(sprintf("%-24s FAILED: %s\n", cfg$tag, conditionMessage(e))); NULL
    })
  if (is.null(fit)) next

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  pm <- predictive_metrics(fit$draws$y0_pred, y0)
  cat(sprintf(
    "%-24s | %6.1f s | div %4d | Rhat %.4f | minESS %6.0f | coverage %.3f | S0 %.2f vs sigma_cond %.2f\n",
    cfg$tag, elapsed, fit$diagnostics$n_divergent, fit$diagnostics$max_rhat,
    fit$diagnostics$min_ess_bulk, pm$coverage,
    mean(fit$draws$S0), mean(fit$draws$sigma_cond)))
  utils::flush.console()
}

cat("\nRead Rhat and ESS, not runtime: a fast fit that has not mixed is worth\n")
cat("nothing. Rhat should be below 1.01 and ESS in the hundreds before any\n")
cat("number from this model is quoted in the paper.\n")
