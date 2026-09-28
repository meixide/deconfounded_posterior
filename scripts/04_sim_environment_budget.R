#!/usr/bin/env Rscript
## 04_sim_environment_budget.R -------------------------------------------
##
## How many environments does reliable causal discovery actually need?
##
## Identification requires the environment means to span R^p, so E >= p + 1.
## The simulations at E = p + 1 exactly show credible intervals for gamma
## covering at about 0.75 against a nominal 0.95, with the *frequentist* GI
## estimator covering worse still (0.633).  Since the frequentist estimator
## has no prior, the shortfall is not a Bayesian artefact, and flattening
## either the mu_e prior or the ridge prior does not close it
## (tests/test_gamma_coverage_sources.R).
##
## The mechanism this script tests is regression dilution.  gamma is the
## within-environment slope minus Sigma^{-1} K, and Sigma^{-1} K is identified
## by regressing the E environment intercepts on the E environment means.  At
## E = p + 1 that regression is exactly saturated: p + 1 observations for
## p + 1 parameters, zero residual degrees of freedom.  The mu_hat_e are
## measured with error of order sqrt(Sigma / n_e), which attenuates the
## estimate of Sigma^{-1} K towards zero, and at the identifiability minimum
## there is no averaging to damp it.  The bias then lands on gamma.
##
## If that is the mechanism, coverage should improve in both directions:
##   * with E, as the between-environment regression gains residual degrees of
##     freedom;
##   * with n_e, as the measurement error in mu_hat_e shrinks.
##
## The output is the practical guidance Referee 1 asks for in point (3): what
## a practitioner should check before trusting a selection, beyond the formal
## spanning condition.
##
## Usage:
##   Rscript scripts/04_sim_environment_budget.R --cell=3 --n-rep=20
##   Rscript scripts/04_sim_environment_budget.R --cell=all --n-rep=2

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
n_rep <- as.integer(parse_flag(cli, "n-rep", "20"))
n_iter <- as.integer(parse_flag(cli, "iter", "2000"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
cores <- max(1L, min(cores, chains))
base_seed <- as.integer(parse_flag(cli, "seed", "20260728"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "env_budget"))

#' Grid of environment counts and per-environment sample sizes at p = 6.
env_budget_grid <- function(p = 6L) {
  g <- expand.grid(
    n_e = c(200L, 800L),
    n_env = c(p + 1L, p + 4L, 2L * p + 3L, 4L * p + 1L),
    stringsAsFactors = FALSE
  )
  g$p <- p
  g$cell_id <- seq_len(nrow(g))
  ## Total N is deliberately *not* held fixed: the question is what each of
  ## the two budgets buys, and confounding them would answer neither.
  g$label <- sprintf("E%02d_n%d", g$n_env, g$n_e)
  g
}

grid <- env_budget_grid()
slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
cell_arg <- if (nzchar(slurm_task)) slurm_task else parse_flag(cli, "cell", "1")
cells <- if (identical(cell_arg, "all")) {
  grid$cell_id
} else {
  as.integer(strsplit(cell_arg, ",", fixed = TRUE)[[1]])
}
stopifnot(all(cells %in% grid$cell_id))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))

for (cid in cells) {
  cell <- grid[grid$cell_id == cid, ]
  out_file <- file.path(out_dir, sprintf("cell_%02d.csv", cid))
  if (file.exists(out_file)) {
    message("Cell ", cid, " already done, skipping.")
    next
  }
  message(sprintf("[%s] cell %d: E = %d, n_e = %d, %d replications",
                  format(Sys.time(), "%H:%M:%S"), cid, cell$n_env, cell$n_e,
                  n_rep))

  rows <- list()
  for (r in seq_len(n_rep)) {
    bgi_set_seed(base_seed + 1000L * cid, r)
    dat <- simulate_gi_data(n_e = cell$n_e, p = cell$p, s0 = 3,
                            n_env = cell$n_env, q = 3, confounding = 2,
                            n0 = 500)

    fit <- tryCatch(
      fit_bgi(dat$x, dat$y, dat$z, dat$x0, model = model, chains = chains,
              iter = n_iter, seed = base_seed + r, cores = cores),
      error = function(e) {
        message("  rep ", r, " failed: ", conditionMessage(e)); NULL
      })
    if (is.null(fit)) next

    pc <- parameter_coverage(fit$draws$gamma, dat$truth$gamma)
    sel <- select_parents(fit$draws$gamma, alpha = 0.05, rule = "sign")
    sm <- support_metrics(sel$selected, dat$truth$parents, cell$p)

    ## Frequentist GI on the same data: no prior, no shrinkage, so it
    ## isolates whatever part of the shortfall is not Bayesian in origin.
    pgi <- fit_pooled_gi(dat$x, dat$y, dat$z)
    pgi_cov <- mean(abs(pgi$gamma - dat$truth$gamma) <=
                      stats::qnorm(0.975) * pgi$se)
    pgi_sm <- support_metrics(which(pgi$p_value < 0.05), dat$truth$parents,
                              cell$p)

    rows[[length(rows) + 1L]] <- data.frame(
      cell_id = cid, label = cell$label, n_env = cell$n_env, n_e = cell$n_e,
      p = cell$p, rep = r,
      ## Residual degrees of freedom in the between-environment regression
      ## that identifies Sigma^{-1} K: E - (p + 1).
      between_df = cell$n_env - (cell$p + 1L),
      bgi_gamma_coverage = pc$gamma_coverage,
      bgi_gamma_rmse = pc$gamma_rmse,
      bgi_post_sd = mean(apply(fit$draws$gamma, 2, stats::sd)),
      bgi_bias = mean(colMeans(fit$draws$gamma) - dat$truth$gamma),
      bgi_tpr = sm$tpr, bgi_fdp = sm$fdp, bgi_fwer = sm$fwer,
      pgi_gamma_coverage = pgi_cov,
      pgi_gamma_rmse = sqrt(mean((pgi$gamma - dat$truth$gamma)^2)),
      pgi_tpr = pgi_sm$tpr, pgi_fdp = pgi_sm$fdp,
      pred_coverage = predictive_metrics(fit$draws$y0_pred, dat$y0)$coverage,
      divergent = fit$diagnostics$n_divergent,
      max_rhat = fit$diagnostics$max_rhat,
      mu_min_sv = min(svd(scale(dat$env_means, scale = FALSE))$d),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) > 0L) {
    res <- do.call(rbind, rows)
    tmp <- paste0(out_file, ".tmp")
    utils::write.csv(res, tmp, row.names = FALSE)
    file.rename(tmp, out_file)
    message(sprintf("  BGI coverage %.3f, pooled GI coverage %.3f, BGI fdp %.3f",
                    mean(res$bgi_gamma_coverage),
                    mean(res$pgi_gamma_coverage), mean(res$bgi_fdp)))
  }
}

message("Results in ", out_dir)
