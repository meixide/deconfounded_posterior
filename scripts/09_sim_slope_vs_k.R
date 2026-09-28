#!/usr/bin/env Rscript
## 09_sim_slope_vs_k.R ----------------------------------------------------
##
## Does writing the model in slope coordinates fix the `gamma` under-coverage,
## and does it matter only when the environment covariances actually differ?
##
## Background.  `stan/gi_hd.stan` follows the companion paper: inner-product
## invariance makes K = Cov(eps_Y, X) common, so the slope on the centred
## covariates is Sigma_e^{-1} K and Sigma_e^{-1} sits inside the mean function.
## Sigma_e is a plug-in estimate, so Sigma_hat_e^{-1} is a noisy *generated
## regressor* whose error is not propagated; simulation-based calibration
## attributes the under-coverage of gamma entirely to this.
##
## `stan/gi_hd_slope.stan` gives each environment its own slope b_e, shrunk
## towards a common centre, so the likelihood contains no covariance at all.
## K_e = Sigma_e b_e is derived and K_bar is used for the target transfer,
## where Sigma_0^{-1} is genuinely needed and is confined to `generated
## quantities`.
##
## Prediction.  The two parameterisations are *identical* when the Sigma_e are
## equal, and should separate as they diverge.  Every simulation in this
## project before this one used a single common Sigma_e, so none of them could
## have detected the difference.
##
## A pilot at six replications gave gamma coverage 0.958 for both at
## heterogeneity 0, and 0.875 (K) against 0.958 (slope) at heterogeneity 0.5.
## The direction matches the theory, but with four coordinates sharing a
## posterior the effective replication count per cell is nearer six than
## twenty-four, so the standard error is about 0.12 and that gap settles
## nothing. This script runs the comparison at a sample size that can.
##
## One array task = one (heterogeneity, replication) pair.
##
## Usage:
##   Rscript scripts/09_sim_slope_vs_k.R --task=all --n-rep=3
##   sbatch --array=1-90 slurm/09_slope_vs_k.sh --n-rep=30

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
n_rep <- as.integer(parse_flag(cli, "n-rep", "30"))
p <- as.integer(parse_flag(cli, "p", "6"))
n_env <- as.integer(parse_flag(cli, "n-env", "12"))
n_e <- as.integer(parse_flag(cli, "n-e", "300"))
iter <- as.integer(parse_flag(cli, "iter", "2000"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
base_seed <- as.integer(parse_flag(cli, "seed", "31415"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "slope_vs_k"))

HETERO <- c(0, 0.25, 1)

tasks <- expand.grid(rep_id = seq_len(n_rep), hetero = HETERO)
tasks$task_id <- seq_len(nrow(tasks))

slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_arg <- if (nzchar(slurm_task)) slurm_task else parse_flag(cli, "task", "1")
task_ids <- if (identical(task_arg, "all")) tasks$task_id else
  as.integer(strsplit(task_arg, ",", fixed = TRUE)[[1]])
stopifnot(all(task_ids %in% tasks$task_id))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
models <- list(
  k_param = load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                           cache_dir = file.path(root, "results", "compiled")),
  slope = load_bgi_model(file.path(root, "stan", "gi_hd_slope.stan"),
                         cache_dir = file.path(root, "results", "compiled"))
)

for (tid in task_ids) {
  tk <- tasks[tasks$task_id == tid, ]
  out_file <- file.path(out_dir, sprintf("task_%04d.csv", tid))
  if (file.exists(out_file)) {
    message("Task ", tid, " already done."); next
  }

  bgi_set_seed(base_seed, tid)
  dat <- simulate_gi_data(n_e = n_e, p = p, s0 = ceiling(p / 2),
                          n_env = n_env, q = 3, confounding = 2, n0 = 400,
                          sigma_heterogeneity = tk$hetero)

  ## How much do the true slopes actually vary?  Recorded so the comparison can
  ## be read against the realised heterogeneity rather than the nominal knob.
  b_mat <- do.call(rbind, dat$truth$b_e)
  b_cv <- mean(apply(b_mat, 2, stats::sd) /
                 pmax(abs(colMeans(b_mat)), 1e-8))

  rows <- list()
  for (nm in names(models)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      fit_bgi(dat$x, dat$y, dat$z, dat$x0, model = models[[nm]],
              chains = chains, iter = iter, cores = cores,
              seed = base_seed + tid, ncp = 0),
      error = function(e) {
        message("  ", nm, " failed: ", conditionMessage(e)); NULL
      })
    if (is.null(fit)) next

    pc <- parameter_coverage(fit$draws$gamma, dat$truth$gamma)
    pm <- predictive_metrics(fit$draws$y0_pred, dat$y0)
    sel <- select_parents(fit$draws$gamma, alpha = 0.05, rule = "sign")
    sm <- support_metrics(sel$selected, dat$truth$parents, p)

    rows[[length(rows) + 1L]] <- data.frame(
      task_id = tid, rep_id = tk$rep_id, hetero = tk$hetero,
      b_cv = b_cv, model = nm, p = p, n_env = n_env, n_e = n_e,
      gamma_coverage = pc$gamma_coverage,
      gamma_coverage_parents = pc$gamma_coverage_parents,
      gamma_coverage_nulls = pc$gamma_coverage_nulls,
      gamma_rmse = pc$gamma_rmse,
      post_sd = mean(apply(fit$draws$gamma, 2, stats::sd)),
      pred_coverage = pm$coverage, pred_is = pm$interval_score,
      pred_rmse = pm$rmse,
      tpr = sm$tpr, fdp = sm$fdp,
      divergent = fit$diagnostics$n_divergent,
      max_rhat = fit$diagnostics$max_rhat,
      min_ess = fit$diagnostics$min_ess_bulk,
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) > 0L) {
    res <- do.call(rbind, rows)
    tmp <- paste0(out_file, ".tmp")
    utils::write.csv(res, tmp, row.names = FALSE)
    file.rename(tmp, out_file)
    message(sprintf("[%s] task %d h=%.2f | %s",
                    format(Sys.time(), "%H:%M:%S"), tid, tk$hetero,
                    paste(sprintf("%s cov %.3f", res$model,
                                  res$gamma_coverage), collapse = " | ")))
  }
}

message("Results in ", out_dir)
