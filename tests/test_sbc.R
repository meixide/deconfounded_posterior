#!/usr/bin/env Rscript
## test_sbc.R -------------------------------------------------------------
##
## Simulation-based calibration (Talts, Betancourt, Simpson, Vehtari and
## Gelman, 2018).
##
## ---- Why this and not the earlier coverage study ---------------------------
##
## The earlier simulations fix gamma* at +/-1 for parents and exactly 0 for
## nulls, then measure how often the credible interval contains it.  That is
## *frequentist* coverage at a fixed truth, and no theorem says it must equal
## 0.95: the Bayesian guarantee is either
##
##   (a) conditional on the data, P(gamma in C(D) | D) = 0.95, which is true by
##       construction and untestable by simulation; or
##   (b) averaged over the prior, P(gamma in C(D)) = 0.95 with the probability
##       taken jointly over gamma ~ pi and D | gamma, which is exact in finite
##       samples but requires gamma to be *drawn from the prior*.
##
## Fixing gamma* at +/-1 violates the condition for (b), so a shortfall there
## diagnoses nothing on its own.  This script instead respects (b): it draws
## every parameter from the model's own prior, generates data from the model's
## own likelihood, and checks the exact guarantee.
##
## The sharpest form of that check is rank uniformity.  If theta* ~ pi and
## D | theta* ~ p, and if the posterior sampler is correct, then the rank of
## theta* among L posterior draws is uniform on {0, ..., L}.  Deviations are
## diagnostic:
##
##   uniform      implementation correct
##   U-shaped     posterior too narrow -- over-confident
##   inverted-U   posterior too wide -- under-confident
##   sloped       posterior biased in one direction
##
## ---- The two arms -----------------------------------------------------------
##
## `known`     the true Sigma_e and the true prior mean for mu_e are handed to
##             the model, so the fitted model is *exactly* the generative one.
##             Rank uniformity here is a clean test of the likelihood, the
##             priors and the sampler. If this fails, there is a bug.
##
## `plugin`    Sigma_e and hmu are estimated from the data as in practice.  The
##             fitted model is then no longer the generative one, and the
##             guarantee no longer applies.  Comparing the two arms isolates
##             how much of the observed under-coverage is attributable to
##             conditioning on plug-in estimates as if they were known --
##             which is precisely Referee 1's point 4.
##
## `cut`       the nuisances are imputed from their first-stage posterior and
##             the second stage is run once per imputation, the draws pooled
##             (`fit_bgi_cut()`).  This targets the cut distribution of
##             Plummer (2015) and Jacob et al. (2017) in place of the naive
##             plug-in.  If the diagnosis is right, this arm should recover the
##             calibration that `plugin` loses.
##
## Two caveats, stated rather than hidden.  The paper's Jeffreys prior on the
## conditional variance is improper and cannot be drawn from, so both arms use
## a proper inverse-gamma in its place (`v_prior_shape`, `v_prior_rate`); this
## is a genuine, if small, departure from the analysis model.  And `hmu` is
## empirical Bayes by design, so even the `known` arm is only exact once it is
## supplied rather than estimated.
##
## Usage:
##   Rscript tests/test_sbc.R [--reps=100] [--arm=known|plugin|cut|both]
##                            [--model=gi_hd|gi_hd_slope] [--hetero=0]
##                            [--tag=<name>]
##
## `--hetero` controls whether the generative Sigma_e differ across
## environments.  It matters: with a common Sigma the slope and covariance
## parameterisations of the model are the same model, so running the comparison
## at `--hetero=0` alone answers nothing.

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
n_reps <- as.integer(parse_flag(cli, "reps", "100"))
arm_arg <- parse_flag(cli, "arm", "both")
p <- as.integer(parse_flag(cli, "p", "4"))
n_env <- as.integer(parse_flag(cli, "n-env", "8"))
n_e <- as.integer(parse_flag(cli, "n-e", "200"))
n_bins <- as.integer(parse_flag(cli, "bins", "20"))
## SBC needs the retained draws to be effectively independent: residual
## autocorrelation pulls ranks towards the middle and shows up as an
## extreme-rank ratio below 1, which is easy to misread as a well-behaved
## posterior that is merely conservative.
##
## The retained draws are spaced `N / n_bins` apart, and independence needs that
## spacing to exceed the autocorrelation time `N / ESS`; the two conditions
## coincide at `ESS = n_bins`.  The default below asks for five times that,
## which is margin without being punitive.  Dropping replications is not free --
## if low ESS is more likely for some prior draws than others, the surviving
## set is no longer a prior sample and SBC loses its guarantee -- so the drop
## counts are reported and `iter` is set high enough to keep them small.
iter <- as.integer(parse_flag(cli, "iter", "4000"))
## Run the two chains in parallel when the allocation allows it; SBC is a long
## sequence of small fits, so halving each one halves the whole run.
n_cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
min_ess <- as.numeric(parse_flag(cli, "min-ess", as.character(5 * n_bins)))
n_imp <- as.integer(parse_flag(cli, "n-imp", "20"))
## "gi_hd" is the K-parameterisation of the paper; "gi_hd_slope" writes the
## same sampling model in terms of the slope b = Sigma^{-1} K, which removes
## Sigma from the Y-likelihood entirely.
model_name <- parse_flag(cli, "model", "gi_hd")
slope_par <- identical(model_name, "gi_hd_slope")
## Heterogeneity of the environment covariances in the generative model.  0
## reproduces the original common-Sigma design; positive values draw each
## Sigma_e from an inverse-Wishart centred on a common matrix.  See the note on
## the `known` arm below.
hetero <- as.numeric(parse_flag(cli, "hetero", "0"))
tag <- parse_flag(cli, "tag", "")

## Prior hyperparameters, matching the model's defaults except for the
## conditional variance, which must be proper to be drawn from.
A_TAU <- 0.5
B_TAU <- 0.5
SD_MU_SCALE <- 2.5
ETA_LKJ <- 2
V_SHAPE <- 3
V_RATE <- 2
SD_B_SCALE <- 1     # half-normal scale for the spread of the b_e; the
                    #   `sd_b_scale` default in fit_bgi()

#' Draw a correlation matrix from LKJ(eta) by the onion method.
rlkj <- function(d, eta) {
  if (d == 1L) return(matrix(1, 1, 1))
  alpha <- eta + (d - 2) / 2
  r <- 2 * stats::rbeta(1, alpha, alpha) - 1
  L <- matrix(0, d, d)
  L[1, 1] <- 1
  L[2, 1] <- r
  L[2, 2] <- sqrt(1 - r^2)
  if (d > 2L) {
    for (m in 2:(d - 1)) {
      alpha <- alpha - 0.5
      y <- stats::rbeta(1, m / 2, alpha)
      u <- stats::rnorm(m)
      u <- u / sqrt(sum(u^2))
      L[m + 1, 1:m] <- sqrt(y) * u
      L[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  tcrossprod(L)
}

#' One draw from the model's joint prior, plus data from its likelihood.
#'
#' Everything here must match the Stan file exactly or SBC tests nothing.  Three
#' points where the two models genuinely differ:
#'
#'  * **What is primitive.**  `gi_hd.stan` puts the prior on `K` and derives the
#'    slope as `Sigma_e^{-1} K`, so the slope is common only when the `Sigma_e`
#'    are.  `gi_hd_slope.stan` puts the prior on the slopes themselves through
#'    `b_e ~ N(b_bar, diag(sd_b)^2)` and derives `K_e = Sigma_e b_e`.  Drawing a
#'    single common `b` for the slope model -- which is what an earlier version
#'    of this script did -- is the point `sd_b = 0`, a boundary of the prior
#'    with density zero.  SBC run that way is invalid, and any result from it
#'    must be discarded.
#'
#'  * **The conditional variance.**  The slope model anchors
#'    `sigma_y^2 = v_raw + max(max_e b_e' Sigma_e b_e, K_bar' Sigma_0^{-1}
#'    K_bar)` so that every environment's conditional variance stays positive,
#'    and the conditional variance in environment `e` is then
#'    `sigma_y^2 - b_e' Sigma_e b_e`.  With heterogeneous `b_e` these differ
#'    across environments, so the generative noise scale has to be recomputed
#'    per environment rather than being `v_raw` throughout.  `Sigma_0` is the
#'    pooled `Sigma_bar` here, because `fit_bgi()` falls back to it when no
#'    target sample is supplied.
#'
#'  * **Heterogeneous covariances.**  `hetero > 0` draws a separate `Sigma_e`
#'    per environment.  This is the configuration the two parameterisations
#'    actually disagree in; with a common `Sigma` they coincide and the
#'    comparison is empty.
sbc_draw <- function() {
  ## Covariances.  A common base matrix, optionally perturbed per environment.
  a <- matrix(stats::rnorm(p * p), p, p)
  sigma_common <- crossprod(a) / p + diag(p)
  sigma_e <- if (hetero <= 0) {
    rep(list(sigma_common), n_env)
  } else {
    nu <- p + 1 + 1 / hetero
    scale_mat <- (nu - p - 1) * sigma_common
    lapply(seq_len(n_env), function(e) {
      m <- solve(stats::rWishart(1, nu, solve(scale_mat))[, , 1])
      enforce_pd((m + t(m)) / 2, shrinkage = 0)
    })
  }
  sigma_bar <- Reduce(`+`, sigma_e) / n_env

  ## Hyperparameters, from the prior.
  v_raw <- 1 / stats::rgamma(1, shape = V_SHAPE, rate = V_RATE)
  tau2 <- stats::rgamma(1, A_TAU, 1) / stats::rgamma(1, B_TAU, 1)
  scale <- sqrt(v_raw) * sqrt(tau2)

  alpha <- stats::rnorm(1, 0, scale)
  gamma <- stats::rnorm(p, 0, scale)

  ## The prior sits on whichever object the model treats as primitive.
  if (slope_par) {
    b_bar <- stats::rnorm(p, 0, scale)
    sd_b <- abs(stats::rnorm(p, 0, SD_B_SCALE))
    b <- lapply(seq_len(n_env), function(e)
      b_bar + sd_b * stats::rnorm(p))
    k_e <- lapply(seq_len(n_env), function(e)
      as.vector(sigma_e[[e]] %*% b[[e]]))
    k <- Reduce(`+`, k_e) / n_env          # K_bar, the transferable object
  } else {
    k <- stats::rnorm(p, 0, scale)
    b <- lapply(sigma_e, function(s) as.vector(solve(s, k)))
    sd_b <- NULL
  }

  ## Conditional variances, computed the way the Stan file computes them.
  quad_e <- vapply(seq_len(n_env),
                   function(e) sum(b[[e]] * (sigma_e[[e]] %*% b[[e]])), 0)
  k0_quad <- sum(k * solve(sigma_bar, k))
  sigma_y_sq <- v_raw + max(max(quad_e), k0_quad)
  v_e <- sigma_y_sq - quad_e
  stopifnot(all(v_e > 0))

  ## Environment means, from their hierarchical prior, centred at a known hmu.
  hmu <- rep(0, p)
  sd_mu <- abs(stats::rnorm(p, 0, SD_MU_SCALE))
  sigma_mu <- diag(sd_mu, p) %*% rlkj(p, ETA_LKJ) %*% diag(sd_mu, p)
  mu <- mvtnorm::rmvnorm(n_env, mean = hmu, sigma = sigma_mu)

  ## Data, from the model's own likelihood.
  blocks <- lapply(seq_len(n_env), function(e) {
    xe <- mvtnorm::rmvnorm(n_e, mean = mu[e, ], sigma = sigma_e[[e]])
    centred <- sweep(xe, 2, mu[e, ], "-")
    mean_y <- alpha + as.vector(xe %*% gamma) +
      as.vector(centred %*% b[[e]])
    list(x = xe, y = stats::rnorm(n_e, mean_y, sqrt(v_e[e])))
  })
  x <- do.call(rbind, lapply(blocks, `[[`, "x"))
  y <- unlist(lapply(blocks, `[[`, "y"), use.names = FALSE)
  z <- rep(seq_len(n_env), each = n_e)

  list(x = x, y = y, z = z, gamma = gamma, alpha = alpha, k = k,
       sigma = sigma_e, hmu = hmu, v_raw = v_raw, sd_b = sd_b)
}

model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

run_arm <- function(arm) {
  ranks <- matrix(NA_integer_, n_reps, p)
  covered <- matrix(NA, n_reps, p)
  n_draws <- NA_integer_
  n_bad_rhat <- 0L
  n_bad_ess <- 0L

  for (r in seq_len(n_reps)) {
    bgi_set_seed(778000, r)
    d <- sbc_draw()

    args <- list(x = d$x, y = d$y, z = d$z, x0 = NULL, model = model,
                 chains = 2, iter = iter, cores = n_cores, seed = r, ncp = 0,
                 v_prior_shape = V_SHAPE, v_prior_rate = V_RATE,
                 sd_mu_scale = SD_MU_SCALE, eta_lkj = ETA_LKJ,
                 a_tau = A_TAU, b_tau = B_TAU, sd_b_scale = SD_B_SCALE)
    if (arm == "known") {
      args$sigma_known <- d$sigma
      args$hmu_known <- d$hmu
    }

    fit <- if (arm == "cut") {
      ## Fewer iterations per imputation: the pooled sample is the union over
      ## imputations, so the total number of draws is comparable.
      args$iter <- 600
      args$n_imp <- n_imp
      tryCatch(do.call(fit_bgi_cut, args), error = function(e) {
        message("  rep ", r, " failed: ", conditionMessage(e)); NULL
      })
    } else {
      tryCatch(do.call(fit_bgi, args), error = function(e) {
        message("  rep ", r, " failed: ", conditionMessage(e)); NULL
      })
    }
    if (is.null(fit)) next
    if (fit$diagnostics$max_rhat > 1.05) {
      message("  rep ", r, " dropped, Rhat ",
              round(fit$diagnostics$max_rhat, 3))
      n_bad_rhat <- n_bad_rhat + 1L
      next
    }
    ## Thinning to `n_bins` only yields independent draws if the effective
    ## sample size is comfortably above `n_bins`; otherwise the rank statistic
    ## is not uniform even for a correct posterior.
    if (fit$diagnostics$min_ess_bulk < min_ess) {
      message("  rep ", r, " dropped, bulk ESS ",
              round(fit$diagnostics$min_ess_bulk))
      n_bad_ess <- n_bad_ess + 1L
      next
    }

    g <- fit$draws$gamma
    ## Thin to `n_bins` effectively independent draws, as SBC requires: the
    ## rank statistic is uniform only for independent posterior samples.
    idx <- seq(1, nrow(g), length.out = n_bins)
    gt <- g[round(idx), , drop = FALSE]
    n_draws <- nrow(gt)
    ranks[r, ] <- colSums(sweep(gt, 2, d$gamma, "<"))
    lo <- apply(g, 2, stats::quantile, 0.025)
    hi <- apply(g, 2, stats::quantile, 0.975)
    covered[r, ] <- d$gamma >= lo & d$gamma <= hi
    if (r %% 10 == 0) message("  ", arm, ": rep ", r, " of ", n_reps)
  }

  ok <- stats::complete.cases(ranks)
  rr <- as.vector(ranks[ok, , drop = FALSE])
  cc <- as.vector(covered[ok, , drop = FALSE])

  ## Chi-square test of rank uniformity on {0, ..., n_draws}.
  tab <- table(factor(rr, levels = 0:n_draws))
  chi <- stats::chisq.test(tab)

  cat(sprintf(
    "\n%-7s | %3d usable reps, %d ranks | 95%% coverage %.3f (se %.3f) | uniformity chi2 p = %.3f\n",
    arm, sum(ok), length(rr), mean(cc),
    sqrt(mean(cc) * (1 - mean(cc)) / length(cc)), chi$p.value))
  cat(sprintf("        dropped: %d for Rhat > 1.05, %d for bulk ESS < %g\n",
              n_bad_rhat, n_bad_ess, min_ess))

  ## Shape of the departure, if any: compare the mass in the extreme bins with
  ## the mass in the middle.  Excess in the tails means the posterior is too
  ## narrow.
  edge <- mean(rr <= 0 | rr >= n_draws)
  edge_expected <- 2 / (n_draws + 1)
  cat(sprintf("        extreme-rank mass %.3f against %.3f expected (ratio %.2f)\n",
              edge, edge_expected, edge / edge_expected))
  list(arm = arm, reps = sum(ok), coverage = mean(cc),
       chisq_p = chi$p.value, edge_ratio = edge / edge_expected,
       ranks = rr, n_draws = n_draws,
       dropped_rhat = n_bad_rhat, dropped_ess = n_bad_ess)
}

arms <- if (arm_arg == "both") c("known", "plugin", "cut") else
  strsplit(arm_arg, ",", fixed = TRUE)[[1]]
cat(sprintf(
  "\nSBC [%s]: p = %d, E = %d, n_e = %d, Sigma heterogeneity %.2f, %d replications per arm\n",
  model_name, p, n_env, n_e, hetero, n_reps))
cat("Rank uniformity is exact in finite samples when the fitted model is the\n")
cat("generative one. Excess mass at the extreme ranks means the posterior is\n")
cat("too narrow.\n")

res <- lapply(arms, run_arm)

out <- do.call(rbind, lapply(res, function(x) data.frame(
  model = model_name, hetero = hetero, arm = x$arm,
  p = p, n_env = n_env, n_e = n_e,
  usable_reps = x$reps, coverage_95 = x$coverage,
  uniformity_p = x$chisq_p, extreme_rank_ratio = x$edge_ratio,
  dropped_rhat = x$dropped_rhat, dropped_ess = x$dropped_ess,
  stringsAsFactors = FALSE)))
## A distinct file per configuration, so parallel array tasks do not overwrite
## one another; `scripts/12_aggregate_sbc.R` collects them.
out_name <- if (nzchar(tag)) sprintf("sbc_%s.csv", tag) else "sbc.csv"
bgi_write_csv(out, file.path(root, "results", "summaries", out_name))

## The ranks themselves, so the histogram can be drawn and the shape of any
## departure judged rather than inferred from two summary numbers.
if (nzchar(tag)) {
  ranks_out <- do.call(rbind, lapply(res, function(x) data.frame(
    model = model_name, hetero = hetero, arm = x$arm,
    rank = x$ranks, n_draws = x$n_draws, stringsAsFactors = FALSE)))
  bgi_write_csv(ranks_out,
                file.path(root, "results", "sbc", sprintf("ranks_%s.csv", tag)))
}

cat("\n=== Summary ===\n")
print(out, row.names = FALSE, digits = 3)
cat("\nIf `known` is uniform (p > 0.05, ratio ~ 1) the implementation is\n")
cat("correct and the fixed-truth under-coverage reported elsewhere is a\n")
cat("property of shrinkage at that particular truth, not a defect.\n")
cat("If `plugin` departs while `known` does not, the cause is conditioning on\n")
cat("estimated covariances as if they were known.\n")
cat("If `cut` restores uniformity, multiple imputation over the first-stage\n")
cat("posterior is a sufficient fix and no change to the Stan model is needed.\n")
cat(sprintf("\nWritten to results/summaries/%s\n", out_name))
