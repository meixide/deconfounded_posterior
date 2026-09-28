#!/usr/bin/env Rscript
## 01_sim_support_recovery.R ----------------------------------------------
##
## Causal support recovery and false discoveries under hidden confounding.
##
## This is the simulation requested by Referee 1, point (5).  The submitted
## manuscript evaluates predictive coverage only, and its data-generating
## process gives every covariate a nonzero slope, so a false discovery is
## impossible by construction and no statement about false discoveries can be
## read off the output.  Here `s0` of the `p` slopes are exactly zero, the
## hidden confounder loads on all `p` covariates, and the following are
## measured for each procedure:
##
##   * parent recovery  (TPR)
##   * per-null error   (FPR) -- what the coordinate-wise sign rule controls
##   * false discovery proportion (FDP) -- averaged over replications, the FDR
##   * familywise error (P(S_hat not a subset of pa(Y))) -- the guarantee that
##     invariant causal prediction targets, included so the comparison with
##     that literature is like for like
##   * sign errors among the selected coordinates
##
## alongside credible-interval coverage for `gamma` and predictive coverage in
## the shifted target domain, so that the selection results and the
## calibration results come from the same fits.
##
## Procedures compared
##   bgi_sign        BGI, local-false-sign-rate rule at alpha (Section 2.1)
##   bgi_ci          BGI, 0 outside the central (1 - alpha) credible interval
##   bgi_fdr         BGI, posterior expected sign-error control at q
##   bgi_rope        BGI, posterior mass outside (-delta, delta) exceeds
##                   1 - alpha; the rule that is well defined under a
##                   continuous prior (see R/selection.R)
##   ols             pooled OLS, t-test at alpha
##   ols_bh          pooled OLS, Benjamini-Hochberg at q
##   pooled_gi       frequentist GI as a pooled least squares fit, t-test
##   pooled_gi_bh    the same, Benjamini-Hochberg at q
##   iv              2SLS with the environment indicators as instruments
##   anchor_g2/8/32  anchor regression at three fixed regularisation strengths
##   anchor_oracle   anchor regression with gamma_anchor chosen to minimise
##                   target-domain RMSE.  Unattainable in practice -- the
##                   target is unlabelled -- and reported as an upper bound on
##                   what any tuning of anchor regression could achieve
##   group_dro       worst-case environment risk (no tuning parameter)
##   icp             invariant causal prediction at alpha
##
## V-REx and Wasserstein DRO are also fitted, but they are penalised point
## estimators with no valid Wald inference, so they appear in the predictive
## comparison only and not in the selection table.
##
## Note that `pooled_gi` and `iv` produce identical point estimates: the
## pooled-OLS-on-environment-means baseline *is* 2SLS with environment
## instruments (see R/baselines.R and tests/test_baseline_identities.R).  Both
## are reported because their standard errors differ, and the selection rules
## depend on those.
##
## Each SLURM array task runs one (scenario, replication) pair and writes one
## CSV, so tasks are independent, restartable and individually reproducible.
##
## Usage:
##   Rscript scripts/01_sim_support_recovery.R --task=17 --n-rep=20 [--small]
##   Rscript scripts/01_sim_support_recovery.R --task=all --n-rep=2 --small
##
## When SLURM_ARRAY_TASK_ID is set it overrides --task.

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

## ---- Arguments ---------------------------------------------------------

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}

cli <- commandArgs(trailingOnly = TRUE)
small <- "--small" %in% cli
n_rep <- as.integer(parse_flag(cli, "n-rep", "20"))
alpha <- as.numeric(parse_flag(cli, "alpha", "0.05"))
q_fdr <- as.numeric(parse_flag(cli, "q", "0.1"))
## Half-width of the region of practical equivalence for the `rope` rule,
## on the scale of gamma (the nonzero slopes have magnitude gamma_signal = 1).
rope_delta <- as.numeric(parse_flag(cli, "rope-delta", "0.15"))
base_seed <- as.integer(parse_flag(cli, "seed", "20260727"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
iter <- as.integer(parse_flag(cli, "iter", "2000"))
## Run the chains across the cores SLURM allocated to this task.  Chain seeds
## are seed + chain_id - 1 regardless of how many run at once, so this changes
## the wall clock and nothing else.
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
cores <- max(1L, min(cores, chains))
cov_method <- parse_flag(cli, "cov-method", "pooled")
## Either "gi_hd" (covariances plugged in, as in the manuscript) or
## "gi_hd_fullcov" (covariances inferred, so their estimation error
## propagates into the posterior for gamma).
model_name <- parse_flag(cli, "model", "gi_hd")
## Step size and parameterisation.  The weak-identifiability scenarios sit near
## the boundary of Assumption 1, where the funnel produces divergences at the
## fit_bgi default of 0.9 often enough that the diagnostic gate discards most
## of the replications: at 0.9 only one in ten survives in either weak cell,
## and a cell averaged over one replication is not an average.  The dimension
## sweep met the same thing and 0.99 largely fixed it.  Both are recorded in
## the output so that a table can always say what produced it.
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.9"))
ncp <- as.numeric(parse_flag(cli, "ncp", "1"))
out_dir <- parse_flag(cli, "out",
                      file.path(root, "results", "support_recovery"))

slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_arg <- if (nzchar(slurm_task)) slurm_task else parse_flag(cli, "task", "1")

tasks <- sim_task_table(n_rep = n_rep, small = small)
task_ids <- if (identical(task_arg, "all")) {
  tasks$task_id
} else {
  as.integer(strsplit(task_arg, ",", fixed = TRUE)[[1]])
}
stopifnot(all(task_ids %in% tasks$task_id))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

## ---- One replication ---------------------------------------------------

#' Run every procedure on one simulated data set and tabulate the results.
#'
#' @param task One row of `sim_task_table()`.
#' @return A data frame with one row per procedure, plus shared columns
#'   describing the scenario and the quality of the BGI fit.
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
  truth <- dat$truth
  p <- task$p

  ## ---- BGI ------------------------------------------------------------
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
  gamma_draws <- fit$draws$gamma
  gamma_hat_bgi <- colMeans(gamma_draws)

  ## ---- Baselines ------------------------------------------------------
  ols <- fit_ols(dat$x, dat$y, dat$x0)
  pgi <- fit_pooled_gi(dat$x, dat$y, dat$z, dat$x0)
  icp <- tryCatch(fit_icp(dat$x, dat$y, dat$z, alpha = alpha),
                  error = function(e) list(selected = integer(0)))

  ## Distribution-shift baselines.  IV needs E - 1 >= p; when the design does
  ## not identify it, record NA rather than silently substituting something
  ## else.
  iv <- tryCatch(fit_iv_2sls(dat$x, dat$y, dat$z, dat$x0),
                 error = function(e) NULL)

  ## Robustness / distribution-shift procedures.  group DRO has no tuning
  ## parameter; V-REx and Wasserstein DRO are penalised point estimators with
  ## no valid Wald inference, so they enter the prediction comparison only.
  gdro <- tryCatch(fit_group_dro(dat$x, dat$y, dat$z, dat$x0),
                   error = function(e) NULL)
  vrex <- tryCatch(fit_vrex(dat$x, dat$y, dat$z, dat$x0, lambda = 10),
                   error = function(e) NULL)
  wdro <- tryCatch(fit_dro_wasserstein(dat$x, dat$y, dat$x0, delta = 0.05),
                   error = function(e) NULL)

  anchor_grid <- c(2, 8, 32)
  anchors <- lapply(anchor_grid, function(g) {
    tryCatch(fit_anchor(dat$x, dat$y, dat$z, dat$x0, gamma_anchor = g),
             error = function(e) NULL)
  })
  names(anchors) <- paste0("anchor_g", anchor_grid)

  ## Oracle tuning: pick gamma_anchor by target-domain RMSE, which requires
  ## the target labels.  Not a usable procedure -- it is the ceiling that any
  ## practical tuning of anchor regression would have to reach.
  anchor_rmse <- vapply(anchors, function(a) {
    if (is.null(a)) NA_real_ else sqrt(mean((a$pred_mean - dat$y0)^2))
  }, numeric(1))
  anchor_oracle <- if (all(is.na(anchor_rmse))) NULL else
    anchors[[which.min(anchor_rmse)]]
  anchor_oracle_gamma <- if (is.null(anchor_oracle)) NA_real_ else
    anchor_grid[which.min(anchor_rmse)]

  ## ---- Selections -----------------------------------------------------
  sel <- list(
    bgi_sign     = select_parents(gamma_draws, alpha = alpha, rule = "sign"),
    bgi_ci       = select_parents(gamma_draws, alpha = alpha, rule = "ci"),
    bgi_fdr      = select_parents(gamma_draws, q = q_fdr, rule = "bayes_fdr"),
    bgi_rope     = select_parents(gamma_draws, alpha = alpha, rule = "rope",
                                  delta = rope_delta),
    ols          = list(selected = which(ols$p_value < alpha)),
    ols_bh       = list(selected = select_bh(ols$p_value, q = q_fdr)),
    pooled_gi    = list(selected = which(pgi$p_value < alpha)),
    pooled_gi_bh = list(selected = select_bh(pgi$p_value, q = q_fdr)),
    icp          = list(selected = icp$selected)
  )
  sel_extra <- list(
    iv            = iv,
    anchor_g2     = anchors[["anchor_g2"]],
    anchor_g8     = anchors[["anchor_g8"]],
    anchor_g32    = anchors[["anchor_g32"]],
    anchor_oracle = anchor_oracle,
    group_dro     = gdro
  )
  for (nm in names(sel_extra)) {
    f <- sel_extra[[nm]]
    sel[[nm]] <- list(selected = if (is.null(f)) integer(0) else
                        which(f$p_value < alpha))
  }
  point_est <- list(
    bgi_sign = gamma_hat_bgi, bgi_ci = gamma_hat_bgi, bgi_fdr = gamma_hat_bgi,
    bgi_rope = gamma_hat_bgi,
    ols = ols$gamma, ols_bh = ols$gamma,
    pooled_gi = pgi$gamma, pooled_gi_bh = pgi$gamma,
    icp = pgi$gamma
  )
  for (nm in names(sel_extra)) {
    f <- sel_extra[[nm]]
    point_est[[nm]] <- if (is.null(f)) rep(NA_real_, p) else f$gamma
  }

  rows <- do.call(rbind, lapply(names(sel), function(nm) {
    m <- support_metrics(sel[[nm]]$selected, truth$parents, p)
    m$method <- nm
    m$sign_error <- sign_error_rate(sel[[nm]]$selected, point_est[[nm]],
                                    truth$gamma)
    m$posterior_esr <- if (!is.null(sel[[nm]]$posterior_esr)) {
      sel[[nm]]$posterior_esr
    } else NA_real_
    m
  }))

  ## ---- Estimation and prediction quality ------------------------------
  bgi_par <- parameter_coverage(gamma_draws, truth$gamma)
  bgi_pred <- predictive_metrics(fit$draws$y0_pred, dat$y0)
  ## ---- Predictive comparison, every method on the same footing ---------
  ## Coverage alone cannot rank methods (a wide interval always covers), so
  ## the interval score is carried alongside it.  `iv` and `pooled_gi` share a
  ## gamma exactly and differ only by the K correction, which makes their
  ## predictive gap a controlled measurement of what K contributes.
  pred_table <- list(
    bgi = list(pred_mean = colMeans(fit$draws$y0_pred),
               pred_lower = apply(fit$draws$y0_pred, 2, stats::quantile, 0.025),
               pred_upper = apply(fit$draws$y0_pred, 2, stats::quantile, 0.975)),
    ols = ols, pooled_gi = pgi, iv = iv,
    anchor_g8 = anchors[["anchor_g8"]], anchor_oracle = anchor_oracle,
    group_dro = gdro, vrex = vrex, wass_dro = wdro
  )

  pred_cols <- list()
  for (nm in names(pred_table)) {
    f <- pred_table[[nm]]
    if (is.null(f) || is.null(f$pred_lower)) {
      m <- list(coverage = NA_real_, mean_width = NA_real_,
                interval_score = NA_real_)
      rmse <- NA_real_
    } else {
      m <- interval_metrics(dat$y0, f$pred_lower, f$pred_upper)
      rmse <- sqrt(mean((f$pred_mean - dat$y0)^2))
    }
    pred_cols[[paste0(nm, "_pred_coverage")]] <- m$coverage
    pred_cols[[paste0(nm, "_pred_width")]] <- m$mean_width
    pred_cols[[paste0(nm, "_pred_is")]] <- m$interval_score
    pred_cols[[paste0(nm, "_pred_rmse")]] <- rmse
  }
  pred_cols <- as.data.frame(pred_cols, stringsAsFactors = FALSE)

  ols_pred <- interval_metrics(dat$y0, ols$pred_lower, ols$pred_upper)
  pgi_pred <- interval_metrics(dat$y0, pgi$pred_lower, pgi$pred_upper)
  na_iv <- list(coverage = NA_real_, mean_width = NA_real_)
  iv_pred <- if (is.null(iv)) na_iv else
    interval_metrics(dat$y0, iv$pred_lower, iv$pred_upper)
  anc_pred <- if (is.null(anchor_oracle)) na_iv else
    interval_metrics(dat$y0, anchor_oracle$pred_lower, anchor_oracle$pred_upper)
  rmse_of <- function(f) if (is.null(f)) NA_real_ else
    sqrt(mean((f$gamma - truth$gamma)^2))

  shared <- data.frame(
    task_id = task$task_id,
    scenario_id = task$scenario_id,
    label = task$label,
    rep_id = task$rep_id,
    confounding = task$confounding,
    identifiability = task$identifiability,
    n_e = task$n_e,
    p = p,
    s0 = task$s0,
    n_env = task$n_env,
    model = model_name,
    adapt_delta = adapt_delta,
    ncp = ncp,
    alpha = alpha,
    q_fdr = q_fdr,
    ## Truth-side descriptors, useful when reading the results.
    k_norm = sqrt(sum(truth$k^2)),
    sigma_y_true = truth$sigma_y,
    ## How far the design is from the identification boundary: the smallest
    ## singular value of the centred environment-mean matrix.
    mu_min_sv = min(svd(scale(dat$env_means, scale = FALSE))$d),
    mu_condition = {
      d <- svd(scale(dat$env_means, scale = FALSE))$d
      max(d) / max(min(d), .Machine$double.eps)
    },
    ## BGI fit quality; rows with divergences or bad Rhat should be inspected
    ## rather than silently averaged in.
    bgi_seconds = bgi_seconds,
    bgi_divergent = fit$diagnostics$n_divergent,
    ## Recorded so the screen below can use a divergence *rate*: an
    ## absolute count says nothing without the number of draws it is out of.
    bgi_post_draws = fit$diagnostics$n_post_draws,
    bgi_max_rhat = fit$diagnostics$max_rhat,
    bgi_min_ess = fit$diagnostics$min_ess_bulk,
    cov_condition_max = max(fit$diagnostics$condition_train),
    cov_condition_target = fit$diagnostics$condition_target,
    cov_shrinkage_max = max(fit$diagnostics$shrinkage_train),
    ## Estimation and calibration, repeated on every row for convenience.
    bgi_gamma_coverage = bgi_par$gamma_coverage,
    bgi_gamma_coverage_parents = bgi_par$gamma_coverage_parents,
    bgi_gamma_coverage_nulls = bgi_par$gamma_coverage_nulls,
    bgi_gamma_rmse = bgi_par$gamma_rmse,
    ols_gamma_rmse = sqrt(mean((ols$gamma - truth$gamma)^2)),
    pooled_gi_gamma_rmse = sqrt(mean((pgi$gamma - truth$gamma)^2)),
    iv_gamma_rmse = rmse_of(iv),
    group_dro_gamma_rmse = rmse_of(gdro),
    vrex_gamma_rmse = rmse_of(vrex),
    wass_dro_gamma_rmse = rmse_of(wdro),
    anchor_g8_gamma_rmse = rmse_of(anchors[["anchor_g8"]]),
    anchor_oracle_gamma_rmse = rmse_of(anchor_oracle),
    anchor_oracle_gamma_anchor = anchor_oracle_gamma,
    iv_pred_coverage = iv_pred$coverage,
    iv_pred_width = iv_pred$mean_width,
    anchor_oracle_pred_coverage = anc_pred$coverage,
    anchor_oracle_pred_width = anc_pred$mean_width,
    bgi_pred_coverage = bgi_pred$coverage,
    bgi_pred_width = bgi_pred$mean_width,
    bgi_pred_rmse = bgi_pred$rmse,
    ols_pred_coverage = ols_pred$coverage,
    ols_pred_width = ols_pred$mean_width,
    pooled_gi_pred_coverage = pgi_pred$coverage,
    pooled_gi_pred_width = pgi_pred$mean_width,
    stringsAsFactors = FALSE
  )

  cbind(shared[rep(1, nrow(rows)), , drop = FALSE],
        pred_cols[rep(1, nrow(rows)), , drop = FALSE], rows)
}

## ---- Main loop ---------------------------------------------------------

for (tid in task_ids) {
  task <- tasks[tasks$task_id == tid, ]
  out_file <- file.path(out_dir, sprintf("task_%05d.csv", tid))

  if (file.exists(out_file)) {
    message("Task ", tid, " already done, skipping.")
    next
  }

  message(sprintf("[%s] task %d/%d  scenario '%s'  rep %d",
                  format(Sys.time(), "%H:%M:%S"), tid, nrow(tasks),
                  task$label, task$rep_id))

  res <- tryCatch(run_one_task(task), error = function(e) {
    message("  FAILED: ", conditionMessage(e))
    NULL
  })

  if (!is.null(res)) {
    ## Write atomically so that a killed job never leaves a truncated CSV
    ## that the aggregator would silently read.
    tmp <- paste0(out_file, ".tmp")
    utils::write.csv(res, tmp, row.names = FALSE)
    file.rename(tmp, out_file)
    message(sprintf("  done in %.1f s (%d divergences, max Rhat %.3f)",
                    res$bgi_seconds[1], res$bgi_divergent[1],
                    res$bgi_max_rhat[1]))
  }
}

message("Results written to ", out_dir)
