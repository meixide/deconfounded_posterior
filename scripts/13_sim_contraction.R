## 13_sim_contraction.R ---------------------------------------------------
##
## Tier 1: does the posterior actually contract at the rate Theorem 2 states?
##
## Every other simulation in this project measures something for which the
## paper proves nothing -- predictive coverage, false discovery proportions,
## repeated-sampling coverage of credible intervals at a fixed truth.  This
## one measures the quantity that appears on the left-hand side of the
## theorem, and nothing else.
##
## Theorem 2 (`thm:rate`) states that for r_N = M sqrt(log N / N),
##
##     Pi_N( ||w - w*|| > r_N | D_N )  <=  C N^{-c_in M^2 / 2},
##
## with w = (beta, K) = (alpha, gamma, K) in R^{2p+1} and
## c_in = lambda_min(S) / (2 sigma_y^2), S = (1/E) sum_e E_e[phi_e phi_e'],
## phi_e(X) = (1, X, X - mu_e)  (Supplement A, ll. 379-383).
##
## What is recorded per replication, all of it about the full vector w and
## not about gamma alone, because the theorem is about w:
##
##   * `tail_mass_M*`   Pi_N(||w - w*|| > M sqrt(log N / N)), the theorem's
##                      own left-hand side, for three values of M.  Theorem 2
##                      predicts this decays polynomially in N.
##   * `radius_q95`     the 95th posterior percentile of ||w - w*||.  This is
##                      the contraction radius; log-log against N it should
##                      have slope about -1/2, the sqrt(log N / N) rate being
##                      indistinguishable from N^{-1/2} over any feasible
##                      range of N.
##   * `post_mean_err`  ||E[w | D] - w*||, the point-estimate corollary of
##                      Theorem 1 stated at l. 574 of the manuscript.
##   * `mass_U`         Pi_N(||w - w*|| <= 0.5), a fixed neighbourhood.
##                      Theorem 1 says this tends to 1.
##   * `lambda_min_S`   the realised curvature constant, computed exactly from
##                      the generative (mu_e, Sigma_e).  Recorded so the rate
##                      can be read against the constant the theory says
##                      governs it rather than against a nominal knob.
##
## Design.  Everything the theorems hold fixed is held fixed: p and E are
## constant (Theorem 3 requires it and Theorem 2 is cleanest that way), the
## environment covariances are homogeneous, and identifiability is strong so
## that lambda_min(S) is comfortably away from zero and the asymptotic regime
## is reachable.  The only thing that moves is n_e, hence N = E * n_e.
##
## The K-parameterisation is used deliberately: w = (beta, K) is literally the
## parameter the theorems are about.  At `sigma_heterogeneity = 0` the slope
## and K parameterisations coincide, so this involves no commitment on the
## question 09_sim_slope_vs_k.R settles.
##
## One array task = one (n_e, replication) pair.
##
## Usage:
##   Rscript scripts/13_sim_contraction.R --task=1 --n-rep=1
##   sbatch --array=1-140 slurm/13_contraction.sh --n-rep=20

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
p <- as.integer(parse_flag(cli, "p", "4"))
n_env <- as.integer(parse_flag(cli, "n-env", "8"))
confounding <- as.numeric(parse_flag(cli, "confounding", "1"))
iter <- as.integer(parse_flag(cli, "iter", "2000"))
chains <- as.integer(parse_flag(cli, "chains", "4"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(chains))))
base_seed <- as.integer(parse_flag(cli, "seed", "27182"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "contraction"))

## Seven sizes spanning a factor of 64 in n_e, so that a -1/2 slope is
## separated from 0 or -1 by far more than the replication noise.
N_E <- c(25, 50, 100, 200, 400, 800, 1600)

tasks <- expand.grid(rep_id = seq_len(n_rep), n_e = N_E)
tasks$task_id <- seq_len(nrow(tasks))

## An explicit `--task=` wins over SLURM_ARRAY_TASK_ID.  The other scripts in
## this directory take the array id unconditionally because there one array
## task is one fit, but here one array task is one *replication* and therefore
## seven task ids, which the launcher passes as a comma-separated list.  With
## the usual precedence the array id silently overrides that list and every
## array task runs only its own id, so the sweep collapses onto the smallest
## sample size and the study measures nothing.
slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
cli_task <- parse_flag(cli, "task", "")
task_arg <- if (nzchar(cli_task)) cli_task else
  if (nzchar(slurm_task)) slurm_task else "1"
task_ids <- if (identical(task_arg, "all")) tasks$task_id else
  as.integer(strsplit(task_arg, ",", fixed = TRUE)[[1]])
stopifnot(all(task_ids %in% tasks$task_id))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))

#' The curvature constant of Theorem 2, computed exactly rather than estimated.
#'
#' `S_e = E_e[phi_e phi_e']` with `phi_e(X) = (1, X, X - mu_e)` and
#' `X ~ N(mu_e, Sigma_e)` has the closed form used below: the cross term
#' between `X` and `X - mu_e` is `Sigma_e`, and `E[(X - mu_e)] = 0` kills the
#' intercept block.  Averaging over environments gives `S`, and
#' `c_in = lambda_min(S) / (2 sigma_y^2)`.
curvature_constant <- function(mus, sigmas, sigma_y) {
  d <- 2 * ncol(mus) + 1
  S <- matrix(0, d, d)
  for (e in seq_len(nrow(mus))) {
    mu <- mus[e, ]
    sg <- sigmas[[e]]
    S <- S + rbind(
      c(1, mu, rep(0, length(mu))),
      cbind(mu, tcrossprod(mu) + sg, sg),
      cbind(rep(0, length(mu)), sg, sg))
  }
  S <- S / nrow(mus)
  lam <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  list(lambda_min = lam, c_in = lam / (2 * sigma_y^2))
}

for (tid in task_ids) {
  tk <- tasks[tasks$task_id == tid, ]
  out_file <- file.path(out_dir, sprintf("task_%04d.csv", tid))
  if (file.exists(out_file)) {
    message("Task ", tid, " already done."); next
  }

  ## Seed by replication, NOT by task.  `simulate_gi_data` draws the truth,
  ## the loadings, the environment means and the covariances before it draws a
  ## single observation, and none of those depend on `n_e`.  Seeding this way
  ## therefore gives every sample size within a replication *the same* design:
  ## same w*, same mu_e, same Sigma_e, hence the same c_in.  Only the number of
  ## observations changes, which is exactly the limit the theorems take.
  ##
  ## This is not a nicety.  Seeding by task instead lets lambda_min(S) vary by
  ## a factor of six across the sweep, and since the contraction rate depends
  ## on it, the N-effect the study exists to measure is swamped by design noise
  ## and the radius comes out non-monotone in N.
  bgi_set_seed(base_seed, tk$rep_id)
  dat <- simulate_gi_data(n_e = tk$n_e, p = p, s0 = ceiling(p / 2),
                          n_env = n_env, q = 3, confounding = confounding,
                          n0 = 200, sigma_heterogeneity = 0,
                          identifiability = "strong")

  N <- n_env * tk$n_e
  w_true <- c(dat$truth$alpha, dat$truth$gamma, dat$truth$k)

  ## `truth$sigma_e` is the per-environment covariance list, so this stays
  ## correct if the homogeneity above is ever relaxed.
  cc <- curvature_constant(
    mus = dat$env_means,
    sigmas = dat$truth$sigma_e,
    sigma_y = dat$truth$sigma_y)

  t0 <- Sys.time()
  fit <- tryCatch(
    fit_bgi(dat$x, dat$y, dat$z, dat$x0, model = model,
            chains = chains, iter = iter, cores = cores,
            seed = base_seed + tid, ncp = 0),
    error = function(e) {
      message("  fit failed: ", conditionMessage(e)); NULL
    })
  if (is.null(fit)) next

  ## Posterior draws of the full w = (alpha, gamma, K).
  w_draws <- cbind(fit$draws$alpha, fit$draws$gamma, fit$draws$K)
  dist <- sqrt(rowSums(sweep(w_draws, 2, w_true, "-")^2))

  r_N <- function(M) M * sqrt(log(N) / N)

  row <- data.frame(
    task_id = tid, rep_id = tk$rep_id, n_e = tk$n_e, N = N,
    p = p, n_env = n_env, confounding = confounding,
    lambda_min_S = cc$lambda_min, c_in = cc$c_in,
    sigma_y = dat$truth$sigma_y,
    r_N_M1 = r_N(1),
    tail_mass_M1 = mean(dist > r_N(1)),
    tail_mass_M2 = mean(dist > r_N(2)),
    tail_mass_M3 = mean(dist > r_N(3)),
    radius_q95 = stats::quantile(dist, 0.95, names = FALSE),
    ## The contraction radius expressed in units of the theorem's own rate.
    ## Theorem 2 is a statement that r_N = M sqrt(log N / N) suffices for some
    ## M, so if the rate is the right one this ratio stops growing with N;
    ## the raw tail masses below saturate at 1 whenever c_in is small, which
    ## it is here, and carry no information in that regime.
    M_star = stats::quantile(dist, 0.95, names = FALSE) / sqrt(log(N) / N),
    radius_q50 = stats::quantile(dist, 0.50, names = FALSE),
    post_mean_err = sqrt(sum((colMeans(w_draws) - w_true)^2)),
    mass_U = mean(dist <= 0.5),
    divergent = fit$diagnostics$n_divergent,
    max_rhat = fit$diagnostics$max_rhat,
    min_ess = fit$diagnostics$min_ess_bulk,
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    stringsAsFactors = FALSE)

  write.csv(row, out_file, row.names = FALSE)
  message(sprintf(
    "task %d | n_e %4d N %5d | radius95 %.4f | tail(M=1) %.3f | c_in %.4f",
    tid, tk$n_e, N, row$radius_q95, row$tail_mass_M1, row$c_in))
}
