#!/usr/bin/env Rscript
## 16_split_extrapolation.R -----------------------------------------------
##
## Train on half the environments; predict into the half that was never seen,
## near targets and far targets alike.
##
## Why leave-one-environment-out cannot settle whether the method is worth
## using.  OLS is the population minimiser of `E(Y - X beta)^2` over whatever
## distribution it is fitted to.  Under LOEO with `E = 52` that distribution is
## a mixture of 51 environments, and the single held-out environment sits
## inside the convex hull of the 51 already seen.  The pooled slope is
## therefore close to optimal for the target *by construction*, and BGI has
## nothing left to recover.  Measured, that is exactly what happens: BGI and
## OLS tie on every dataset tried, and the tie survives at sixty times the
## domain shift (Communities and Crime), which no shift-based explanation
## accounts for.
##
## LOEO is still the right answer to Referee 2's objection, which is about
## *calibration*: coverage from a single held-out domain has a standard error
## far larger than it appears, and one domain cannot validate generalisation.
## Holding out every environment in turn, with the environment as the unit of
## replication, answers that. It simply answers a different question from "is
## this method worth using instead of pooled regression", which is what the
## Associate Editor asked.
##
## This script answers the second question. Each replicate splits the
## environments into a training half and a held-out half, ranks the held-out
## environments by how far their covariate mean sits from the training half,
## and evaluates the furthest and the nearest.  The far targets are genuine
## extrapolation: outside the convex hull of everything the model was fitted
## to.
##
## Two design choices that keep the comparison honest.
##
## Selecting targets by distance is itself a selection, so the *near* targets
## are evaluated too, from the same training half.  The claim under test is
## therefore "the advantage grows with extrapolative distance", a contrast
## within a replicate, rather than "BGI wins on far targets", which selection
## alone could manufacture.
##
## A single split is one draw.  Replicating over random halves restores the
## environment as the unit of replication, so this does not reintroduce the
## single-test-environment problem that the referee objected to.
##
## Usage, as an array (see slurm/16_split_extrapolation.sh):
##   Rscript scripts/16_split_extrapolation.R --task=1 --dataset=brfss
##
## Task numbering enumerates (replicate, target) pairs deterministically, so
## each array task reproduces the same split without coordination.

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
task <- as.integer(parse_flag(cli, "task", "1"))
dataset <- parse_flag(cli, "dataset", "brfss")
data_path <- parse_flag(cli, "data", switch(
  dataset,
  brfss = file.path(dirname(root), "data", "brfss", "brfss2023_case.csv"),
  quiron = file.path(dirname(root), "old_code", "quiron", "quiron_final.csv"),
  acs_tract = file.path(dirname(root), "data", "candidates", "tract",
                        "acs_tract_2022.csv"),
  acs_pums = file.path(dirname(root), "data", "candidates", "acs",
                       "acs_slim.csv"),
  communities = file.path(dirname(root), "data", "candidates", "cc",
                          "communities.data"),
  stop("Pass --data for dataset ", dataset)))
model_name <- parse_flag(cli, "model", "gi_hd_slope")
n_reps <- as.integer(parse_flag(cli, "reps", "10"))
k_each <- as.integer(parse_flag(cli, "k", "3"))
min_n0 <- as.integer(parse_flag(cli, "min-n0", "300"))
max_target <- as.integer(parse_flag(cli, "max-target", "2000"))
n_iter <- as.integer(parse_flag(cli, "iter", "2000"))
n_chains <- as.integer(parse_flag(cli, "chains", "4"))
ncp <- as.numeric(parse_flag(cli, "ncp", "0"))
sd_b_scale <- as.numeric(parse_flag(cli, "sd-b-scale", "1"))
v_prior <- parse_flag(cli, "v-prior", "off")
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.95"))
base_seed <- as.integer(parse_flag(cli, "seed", "20260804"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(n_chains))))
out_dir <- parse_flag(cli, "out",
                      file.path(root, "results",
                                paste0("split_", dataset, "_", model_name)))

dat <- load_case_data(data_path, dataset = dataset)
sizes <- table(dat$z)
eligible <- names(sizes)[sizes >= min_n0]
if (length(eligible) < 8L) {
  stop("Only ", length(eligible), " eligible environments; a half-split needs ",
       "more.", call. = FALSE)
}

#' The (replicate, role, rank) that a task index refers to.
#'
#' Enumerated rather than randomised per task so that every task can rebuild
#' the same split from the replicate index alone.
plan_of <- function(task, n_reps, k_each) {
  per_rep <- 2L * k_each
  rep_id <- ((task - 1L) %/% per_rep) + 1L
  within <- ((task - 1L) %% per_rep) + 1L
  role <- if (within <= k_each) "far" else "near"
  rank <- if (within <= k_each) within else within - k_each
  list(rep_id = rep_id, role = role, rank = rank)
}
n_tasks <- n_reps * 2L * k_each
if (task < 1L || task > n_tasks) {
  stop("--task must be in 1:", n_tasks, call. = FALSE)
}
plan <- plan_of(task, n_reps, k_each)

## ---- Rebuild this replicate's split -------------------------------------
set.seed(base_seed + plan$rep_id)
n_train_env <- floor(length(eligible) / 2)
train_env <- sort(sample(eligible, n_train_env))
hold_env <- setdiff(eligible, train_env)

## Distance of each held-out environment from the training half, in the metric
## of the pooled within-environment covariance of the training half.  This is
## the same quantity `scripts/10_environment_shift.R` reports, restricted to
## what the model will actually have seen.
in_train <- dat$z %in% train_env
xs_all <- scale(dat$x)
sw <- pooled_covariance(xs_all[in_train, , drop = FALSE], dat$z[in_train])
mu_train <- colMeans(xs_all[in_train, , drop = FALSE])
mu_hold <- t(vapply(hold_env, function(e)
  colMeans(xs_all[dat$z == e, , drop = FALSE]), numeric(ncol(xs_all))))
dist_hold <- stats::mahalanobis(mu_hold, mu_train, sw)
names(dist_hold) <- hold_env

## Identification strength of the training half, reported rather than optimised.
## The half is drawn at random; `lambda_min` of the averaged augmented
## second-moment matrix is the quantity Assumption 1 is about, and `N * lambda_min`
## its concentration-parameter analogue. Selecting halves to maximise it would
## buy little -- greedily dropping 32 of 52 environments gains 1.7x -- and would
## invite exactly the objection that the design was chosen to suit the method.
mu_tr_env <- t(vapply(train_env, function(e)
  colMeans(xs_all[dat$z == e, , drop = FALSE]), numeric(ncol(xs_all))))
S_train <- crossprod(cbind(1, mu_tr_env)) / length(train_env)
ev_train <- eigen(S_train, symmetric = TRUE)$values
lambda_min <- min(ev_train)
lambda_cond <- max(ev_train) / lambda_min

ordered_far <- names(sort(dist_hold, decreasing = TRUE))
target <- if (plan$role == "far") {
  ordered_far[plan$rank]
} else {
  rev(ordered_far)[plan$rank]
}

## ---- Fit ----------------------------------------------------------------
idx_tr <- which(in_train)
i0 <- which(dat$z == target)
set.seed(base_seed + task)
if (length(i0) > max_target) i0 <- sort(sample(i0, max_target))

model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

resid_var <- {
  xtr <- dat$x[idx_tr, , drop = FALSE]; ytr <- dat$y[idx_tr]
  ztr <- dat$z[idx_tr]
  rs <- unlist(lapply(split(seq_along(ztr), ztr), function(i) {
    if (length(i) <= ncol(xtr) + 2L) return(numeric(0))
    stats::residuals(stats::lm(ytr[i] ~ xtr[i, , drop = FALSE]))
  }), use.names = FALSE)
  stats::var(rs)
}
v_shape <- if (identical(v_prior, "auto")) 3 else 0
v_rate <- if (identical(v_prior, "auto")) 2 * resid_var else 0

cat(sprintf("\ntask %d | rep %d | %s rank %d | target %s\n",
            task, plan$rep_id, plan$role, plan$rank, target))
cat(sprintf("train envs %d | held-out envs %d | N_train %d | p %d\n",
            length(train_env), length(hold_env), length(idx_tr), ncol(dat$x)))
cat(sprintf("target distance %.3f (held-out range %.3f to %.3f) | n_target %d\n",
            dist_hold[[target]], min(dist_hold), max(dist_hold), length(i0)))
cat(sprintf("training half: lambda_min %.4g | cond %.0f | N*lambda_min %.1f | E-(p+1) %d\n\n",
            lambda_min, lambda_cond, length(idx_tr) * lambda_min,
            length(train_env) - (ncol(dat$x) + 1L)))

x_tr <- dat$x[idx_tr, , drop = FALSE]
y_tr <- dat$y[idx_tr]
z_tr <- dat$z[idx_tr]
x0 <- dat$x[i0, , drop = FALSE]
y0 <- dat$y[i0]

fit <- tryCatch(
  fit_bgi(x_tr, y_tr, z_tr, x0, model = model, ncp = ncp,
          sd_b_scale = sd_b_scale, v_prior_shape = v_shape,
          v_prior_rate = v_rate, chains = n_chains, iter = n_iter,
          cores = cores, adapt_delta = adapt_delta, seed = base_seed + task),
  error = function(e) { cat("fit failed: ", conditionMessage(e), "\n"); NULL })
if (is.null(fit)) quit(save = "no", status = 0)

pm <- predictive_metrics(fit$draws$y0_pred, y0)
row <- data.frame(
  task = task, rep_id = plan$rep_id, role = plan$role, rank = plan$rank,
  dataset = dataset, model = model_name, target_env = target,
  distance = dist_hold[[target]],
  n_train_env = length(train_env), n_train = length(idx_tr),
  n_target = length(i0), p = ncol(dat$x),
  between_df = length(train_env) - (ncol(dat$x) + 1L),
  lambda_min = lambda_min, lambda_cond = lambda_cond,
  n_lambda_min = length(idx_tr) * lambda_min,
  dist_min = min(dist_hold), dist_max = max(dist_hold),
  coverage = pm$coverage, interval_score = pm$interval_score, rmse = pm$rmse,
  S0_mean = mean(fit$draws$S0),
  sigma_cond_mean = mean(fit$draws$sigma_cond),
  divergent = fit$diagnostics$n_divergent,
  post_draws = fit$diagnostics$n_post_draws,
  max_rhat = fit$diagnostics$max_rhat,
  min_ess = fit$diagnostics$min_ess_bulk,
  n_chains = fit$diagnostics$n_chains,
  runtime_sec = fit$diagnostics$runtime_sec,
  stringsAsFactors = FALSE
)

for (nm in c("ols", "pooled_gi", "iv")) {
  f <- tryCatch(switch(nm,
                       ols = fit_ols(x_tr, y_tr, x0),
                       pooled_gi = fit_pooled_gi(x_tr, y_tr, z_tr, x0),
                       iv = fit_iv_2sls(x_tr, y_tr, z_tr, x0)),
                error = function(e) NULL)
  if (is.null(f) || is.null(f$pred_lower)) {
    row[[paste0(nm, "_coverage")]] <- NA_real_
    row[[paste0(nm, "_interval_score")]] <- NA_real_
    row[[paste0(nm, "_rmse")]] <- NA_real_
  } else {
    m <- interval_metrics(y0, f$pred_lower, f$pred_upper)
    row[[paste0(nm, "_coverage")]] <- m$coverage
    row[[paste0(nm, "_interval_score")]] <- m$interval_score
    row[[paste0(nm, "_rmse")]] <- sqrt(mean((f$pred_mean - y0)^2))
  }
}

print(t(row))
bgi_write_csv(row, file.path(out_dir, sprintf("task_%03d.csv", task)))
cat("\nwritten to ", file.path(out_dir, sprintf("task_%03d.csv", task)), "\n",
    sep = "")
