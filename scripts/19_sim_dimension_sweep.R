#!/usr/bin/env Rscript
## 19_sim_dimension_sweep.R ----------------------------------------------
##
## Predictive coverage as the covariate dimension grows (Table 1,
## `table:AEC`, Section 3.1.1).
##
## Why this script exists
## ----------------------
## The dimension sweep reported in the submitted version was produced by a
## separate, earlier codebase, whose Stan model drew the predictive sample as
##
##     Y_pred[i] ~ normal(cond_mean, sigmay);        // training residual sd
##
## whereas the model in `stan/gi_hd.stan` draws
##
##     Y0_pred[i] ~ normal(f0[i], S0);               // target-domain scale
##
## with `S0^2 = sigma_y^2 - K' Sigma_0^{-1} K`.  Carrying the target-domain
## scale into the prediction is the contribution the revision claims, so a
## table of predictive coverage produced under the training scale cannot be
## used to support it.  This script reruns the sweep against the current
## implementation, on the same simulation infrastructure as every other
## experiment in the package.
##
## It also records what the *old* scale would have given on the very same
## fits, in the `cov_train_scale` column.  The two columns are paired -- same
## data, same posterior, same standard normal variates, only the predictive
## standard deviation differs -- so their difference isolates the effect of
## the correction rather than confounding it with Monte Carlo noise.  That
## comparison is the direct answer to Referee 2's minor comment (2), and it
## explains why the numbers in the revised Table 1 differ from the submitted
## ones.
##
## Design
## ------
## `p` in {2, 5, 10} crossed with `n_e` in {200, 500, 1000, 2000}, 24
## replications per cell, `E = p + 1` environments so that Assumption 1 holds
## at its minimum, `q = 3` hidden confounders loading on every covariate, and
## a target domain with both a shifted mean and an inflated dispersion.  All
## `p` covariates are parents (`s0 = p`), as in the mechanism the submitted
## version used.
##
## 3 * 4 * 24 = 288 tasks.  One task is one replication of one cell and is a
## single BGI fit, so the array is embarrassingly parallel.
##
## Usage:
##   Rscript scripts/19_sim_dimension_sweep.R --task=1
##   Rscript scripts/19_sim_dimension_sweep.R --task=all --n-rep=1   # smoke
##   sbatch --array=1-288 slurm/19_dimension_sweep.sh
##
## If the cluster caps submitted jobs per user, group replications so that the
## array fits under the cap; --chunk=6 gives 48 indices, --chunk=24 gives 12,
## one per cell:
##   sbatch --array=1-48 slurm/19_dimension_sweep.sh --chunk=6
##
## Then:
##   Rscript scripts/20_aggregate_dimension_sweep.R

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
n_rep <- as.integer(parse_flag(cli, "n-rep", "24"))
base_seed <- as.integer(parse_flag(cli, "seed", "20260923"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
iter <- as.integer(parse_flag(cli, "iter", "2000"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
cores <- max(1L, min(cores, chains))
## Replications per array task.  The task table has one row per
## (cell, replication) pair, but a cluster QoS usually caps how many jobs a
## user may have *submitted*, not merely running, and `--array=1-288%48`
## counts as 288 against that cap.  With `--chunk=N` one array index covers N
## consecutive rows, so the array becomes ceiling(rows / N) indices and the
## same work fits under the cap.  Rows are ordered cell-major, so a chunk
## holds replications of a single cell and its runtime is N times one
## replication of that cell.
chunk <- as.integer(parse_flag(cli, "chunk", "1"))
stopifnot(chunk >= 1L)
## Step size adaptation.  At p = 2 the design has only E = p + 1 = 3
## environments, which is Assumption 1 at its boundary and leaves the
## between-environment information about K as thin as it can be; the resulting
## funnel produces occasional divergences at the fit_bgi default of 0.9.  A
## divergence is a sampler-tuning problem here rather than a model failure,
## and discarding the replications that hit one selects on posterior geometry,
## hence on the data, which is exactly the bias the diagnostic gate exists to
## avoid.  Raising this is much cheaper than throwing fits away.
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.9"))
## Centred (0) versus non-centred (1) parameterisation of (alpha, gamma, K).
## This is a change of coordinates, not of model: gi_hd.stan scales the raw
## parameters so that the implied prior is the same for every ncp.  Which one
## samples well depends on how informative the likelihood is, as
## scripts/05_case_study.R records: non-centred when it is weak, centred when
## the data dominate, where it is far faster and stops hitting the treedepth
## ceiling.  The large-n, large-p cells of this sweep are in the second
## regime -- p = 10, n_e = 2000 is E = 11 environments of 2000 observations --
## so they want ncp = 0 even though the small cells want ncp = 1.
ncp <- as.numeric(parse_flag(cli, "ncp", "1"))
cov_method <- parse_flag(cli, "cov-method", "pooled")
model_name <- parse_flag(cli, "model", "gi_hd")
level <- as.numeric(parse_flag(cli, "level", "0.95"))
out_dir <- parse_flag(cli, "out",
                      file.path(root, "results", "dimension_sweep"))

#' The p-by-n grid of Table 1, expanded to one row per replication.
sweep_task_table <- function(n_rep = 24L,
                             p_grid = c(2L, 5L, 10L),
                             n_grid = c(200L, 500L, 1000L, 2000L)) {
  grid <- expand.grid(p = p_grid, n_e = n_grid,
                      KEEP.OUT.ATTRS = FALSE)
  ## E = p + 1 is the minimum number of environments for the augmented means
  ## to span R^{p+1}, which is Assumption 1 at its boundary.
  grid$n_env <- grid$p + 1L
  grid$q <- 3L
  ## Every covariate is a parent, as in the mechanism the submitted version
  ## used: this experiment is about predictive calibration, not selection.
  grid$s0 <- grid$p
  grid$n0 <- grid$n_e
  grid$confounding <- 1
  grid$gamma_signal <- 1
  grid$identifiability <- "strong"
  grid$cell_id <- seq_len(nrow(grid))
  grid$label <- sprintf("p%d_n%d", grid$p, grid$n_e)

  tasks <- grid[rep(seq_len(nrow(grid)), each = n_rep), , drop = FALSE]
  tasks$rep_id <- rep(seq_len(n_rep), times = nrow(grid))
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  tasks
}

slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_arg <- if (nzchar(slurm_task)) slurm_task else parse_flag(cli, "task", "1")

tasks <- sweep_task_table(n_rep = n_rep)
n_rows <- nrow(tasks)

task_ids <- if (identical(task_arg, "all")) {
  tasks$task_id
} else {
  requested <- as.integer(strsplit(task_arg, ",", fixed = TRUE)[[1]])
  stopifnot(!anyNA(requested), all(requested >= 1L))
  if (chunk > 1L) {
    n_chunks <- ceiling(n_rows / chunk)
    if (any(requested > n_chunks)) {
      stop("chunk index ", max(requested), " is past the end: with --chunk=",
           chunk, " and ", n_rows, " rows there are ", n_chunks,
           " array indices.")
    }
    unlist(lapply(requested, function(k) {
      seq.int((k - 1L) * chunk + 1L, min(k * chunk, n_rows))
    }))
  } else {
    requested
  }
}
stopifnot(all(task_ids %in% tasks$task_id))

if (chunk > 1L) {
  message(sprintf("chunk %d: this index covers rows %d-%d of %d",
                  chunk, min(task_ids), max(task_ids), n_rows))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

## ---- One replication ---------------------------------------------------

run_one_task <- function(task) {
  bgi_set_seed(base_seed, task$task_id)

  dat <- simulate_gi_data(
    n_e = task$n_e,
    p = task$p,
    s0 = task$s0,
    n_env = task$n_env,
    q = task$q,
    confounding = task$confounding,
    gamma_signal = task$gamma_signal,
    n0 = task$n0,
    identifiability = task$identifiability
  )

  t0 <- Sys.time()
  fit <- fit_bgi(
    x = dat$x, y = dat$y, z = dat$z, x0 = dat$x0,
    model = model,
    cov_method = cov_method,
    chains = chains,
    iter = iter,
    seed = base_seed + task$task_id,
    cores = cores,
    adapt_delta = adapt_delta,
    ncp = ncp
  )
  bgi_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## ---- Ours, at the scale the model actually uses ----------------------
  ours <- predictive_metrics(fit$draws$y0_pred, dat$y0, level = level)

  ## ---- Ours, rescored at the training residual scale -------------------
  ## Same fit, same conditional means, same standard normal variates; only
  ## the predictive sd changes.  This is what the submitted implementation
  ## reported, and the gap between the two coverages is the correction.
  n_draws <- nrow(fit$draws$f0)
  z_common <- matrix(stats::rnorm(n_draws * task$n0), nrow = n_draws)
  at_s0 <- predictive_metrics_at_scale(
    fit$draws$f0, fit$draws$S0, dat$y0, level = level, z = z_common)
  at_train <- predictive_metrics_at_scale(
    fit$draws$f0, fit$draws$sigma_cond, dat$y0, level = level, z = z_common)

  ## ---- Least squares baseline -----------------------------------------
  ## Large-sample normal prediction intervals, which is the comparison the
  ## submitted version reported alongside ours.
  ols <- fit_ols(dat$x, dat$y, dat$x0, level = level)
  ols_m <- interval_metrics(dat$y0, ols$pred_lower, ols$pred_upper,
                            level = level)

  diag <- fit$diagnostics

  data.frame(
    task_id = task$task_id,
    cell_id = task$cell_id,
    label = task$label,
    rep_id = task$rep_id,
    p = task$p,
    n_e = task$n_e,
    n_env = task$n_env,
    n0 = task$n0,
    cov_ours = ours$coverage,
    score_ours = ours$interval_score,
    width_ours = ours$mean_width,
    rmse_ours = ours$rmse,
    cov_s0_scale = at_s0$coverage,
    cov_train_scale = at_train$coverage,
    score_s0_scale = at_s0$interval_score,
    score_train_scale = at_train$interval_score,
    cov_ols = ols_m$coverage,
    score_ols = ols_m$interval_score,
    width_ols = ols_m$mean_width,
    divergences = diag$n_divergent,
    max_rhat = diag$max_rhat,
    min_ess_bulk = diag$min_ess_bulk,
    adapt_delta = adapt_delta,
    ncp = ncp,
    bgi_seconds = bgi_seconds,
    stringsAsFactors = FALSE
  )
}

## ---- Run ---------------------------------------------------------------

for (i in seq_along(task_ids)) {
  id <- task_ids[i]
  task <- tasks[tasks$task_id == id, , drop = FALSE]
  message(sprintf("[%s] task %d/%d  cell '%s'  rep %d",
                  format(Sys.time(), "%H:%M:%S"), i, length(task_ids),
                  task$label, task$rep_id))

  res <- run_one_task(task)
  utils::write.csv(
    res,
    file.path(out_dir, sprintf("task_%04d.csv", id)),
    row.names = FALSE)

  message(sprintf("  done in %.1f s (%d divergences, max Rhat %.3f)",
                  res$bgi_seconds, res$divergences, res$max_rhat))
}

message("Results written to ", out_dir)
