## metrics.R --------------------------------------------------------------
##
## Summaries of support recovery and of predictive calibration.
##
## The support-recovery summaries deliberately separate three different error
## notions that the manuscript currently conflates:
##
##   fdp   the realised false discovery *proportion* of a single selected set;
##         averaging it over replications gives the FDR.
##   fpr   the proportion of true nulls that were selected; averaging it gives
##         the per-null error rate, which is what a coordinate-wise rule such
##         as `sign` controls.
##   fwer  the indicator that the selected set is not contained in pa(Y);
##         averaging it estimates 1 - P(S_hat subset of pa(Y)), the guarantee
##         that invariant causal prediction targets and the natural point of
##         comparison with it.
##
## A rule can look excellent on one and poor on another, so all three are
## reported.

#' Compare a selected set with the true parent set.
#'
#' @param selected Integer vector of selected covariate indices.
#' @param parents Integer vector of true parents.
#' @param p Total number of covariates.
#' @return A one-row data frame of recovery metrics.
support_metrics <- function(selected, parents, p) {
  selected <- as.integer(selected)
  parents <- as.integer(parents)
  nulls <- setdiff(seq_len(p), parents)

  tp <- length(intersect(selected, parents))
  fp <- length(setdiff(selected, parents))
  fn <- length(setdiff(parents, selected))
  tn <- length(nulls) - fp

  denom_mcc <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

  data.frame(
    n_selected = length(selected),
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    ## Power / parent recovery.
    tpr = if (length(parents) > 0) tp / length(parents) else NA_real_,
    ## Per-null error rate.
    fpr = if (length(nulls) > 0) fp / length(nulls) else NA_real_,
    ## Realised false discovery proportion.
    fdp = if (length(selected) > 0) fp / length(selected) else 0,
    ## Familywise error: did we claim any non-parent at all?
    fwer = as.integer(fp > 0),
    exact_recovery = as.integer(setequal(selected, parents)),
    jaccard = if (length(union(selected, parents)) > 0) {
      length(intersect(selected, parents)) / length(union(selected, parents))
    } else 1,
    mcc = if (denom_mcc > 0) (tp * tn - fp * fn) / denom_mcc else 0,
    stringsAsFactors = FALSE
  )
}

#' Sign-error rate among the selected true parents.
#'
#' The `sign` rule claims a sign for every coordinate it selects; this is the
#' realised frequency with which that claim is wrong, and is the quantity the
#' rule's alpha is meant to bound.
#'
#' Coordinates with `gamma_true == 0` are excluded.  Selecting one of them is
#' a false discovery, already counted by `fdp` and `fpr`, and scoring it here
#' as well would conflate two distinct errors: no sign is correct for an exact
#' zero, so including them would make the sign-error rate track the false
#' discovery rate rather than measure sign accuracy.
#'
#' @param selected Integer vector of selected indices.
#' @param gamma_hat Posterior mean (or point estimate) of `gamma`.
#' @param gamma_true True `gamma`.
#' @return The proportion of selected true parents given the wrong sign, or
#'   `NA` when no true parent was selected.
sign_error_rate <- function(selected, gamma_hat, gamma_true) {
  selected <- intersect(selected, which(gamma_true != 0))
  if (length(selected) == 0L) {
    return(NA_real_)
  }
  mean(sign(gamma_hat[selected]) != sign(gamma_true[selected]))
}

#' Coverage, width and interval score of a set of prediction intervals.
#'
#' Coverage alone cannot rank methods: an infinitely wide interval covers
#' perfectly.  The interval score (Gneiting and Raftery, 2007) is the proper
#' scoring rule for a central `1 - alpha` interval,
#'
#'   `IS = (u - l) + (2/alpha)(l - y) 1{y < l} + (2/alpha)(y - u) 1{y > u}`,
#'
#' which rewards narrow intervals and penalises misses.  Lower is better.  It
#' is the right headline number for a paper whose contribution is a predictive
#' *distribution* rather than a point prediction, and it is the comparison on
#' which the distribution-shift baselines can be judged at all, since most of
#' them supply only a point predictor plus a training residual scale.
#'
#' @param truth Observed target responses.
#' @param lower,upper Interval endpoints.
#' @param level Nominal level the interval was built at.
#' @return A one-row data frame.
interval_metrics <- function(truth, lower, upper, level = 0.95) {
  a <- 1 - level
  width <- upper - lower
  below <- pmax(lower - truth, 0)
  above <- pmax(truth - upper, 0)
  data.frame(
    coverage = mean(truth >= lower & truth <= upper),
    mean_width = mean(width),
    interval_score = mean(width + (2 / a) * (below + above)),
    stringsAsFactors = FALSE
  )
}

#' Coverage of marginal credible intervals for the entries of `gamma`.
#'
#' @param gamma_draws A `draws x p` matrix.
#' @param gamma_true True `gamma`.
#' @param level Nominal level.
parameter_coverage <- function(gamma_draws, gamma_true, level = 0.95) {
  a <- (1 - level) / 2
  lower <- apply(gamma_draws, 2, stats::quantile, probs = a)
  upper <- apply(gamma_draws, 2, stats::quantile, probs = 1 - a)
  inside <- gamma_true >= lower & gamma_true <= upper
  data.frame(
    gamma_coverage = mean(inside),
    gamma_coverage_parents = if (any(gamma_true != 0)) {
      mean(inside[gamma_true != 0])
    } else NA_real_,
    gamma_coverage_nulls = if (any(gamma_true == 0)) {
      mean(inside[gamma_true == 0])
    } else NA_real_,
    gamma_rmse = sqrt(mean((colMeans(gamma_draws) - gamma_true)^2)),
    stringsAsFactors = FALSE
  )
}

#' Predictive interval summaries from posterior predictive draws.
#'
#' @param y0_pred A `draws x n0` matrix of posterior predictive draws.
#' @param y0 Observed target responses.
#' @param level Nominal level.
predictive_metrics <- function(y0_pred, y0, level = 0.95) {
  a <- (1 - level) / 2
  lower <- apply(y0_pred, 2, stats::quantile, probs = a)
  upper <- apply(y0_pred, 2, stats::quantile, probs = 1 - a)
  m <- interval_metrics(y0, lower, upper)
  m$rmse <- sqrt(mean((colMeans(y0_pred) - y0)^2))
  m
}

#' Predictive metrics under a substituted predictive scale.
#'
#' Rebuilds the posterior predictive from the conditional mean draws `f0` with
#' the per-draw standard deviation `sigma_draws` in place of the one the model
#' actually used, and scores it.
#'
#' This exists to answer the question Referee 2's minor comment 2 poses
#' directly.  Having conceded that a pooled OLS regression on `X` and the
#' environment means reproduces the point predictions, the referee identifies
#' the added value of the Bayesian treatment as "(i) correct prediction
#' variance `sigma_Y^2 - K' Sigma_0^{-1} K` in the test domain rather than the
#' training residual variance".  That is a claim about calibration, and it is
#' testable: score the same fit twice, once at `S_0` and once at the training
#' scale the submitted code used, and compare.  Without this comparison the
#' paper can say the corrected formula is calibrated but not that it improves
#' on the formula it replaces.
#'
#' The same standard normal draws are reused across scales so that the
#' comparison is paired and the difference is not confounded with Monte Carlo
#' noise.
#'
#' @param f0 A `draws x n0` matrix of target-domain conditional means.
#' @param sigma_draws Length-`draws` vector of predictive standard deviations.
#' @param y0 Observed target responses.
#' @param level Nominal level.
#' @param z Optional `draws x n0` matrix of standard normal variates, so that
#'   several scales can be scored against identical randomness.
predictive_metrics_at_scale <- function(f0, sigma_draws, y0, level = 0.95,
                                        z = NULL) {
  f0 <- as.matrix(f0)
  if (is.null(z)) {
    z <- matrix(stats::rnorm(length(f0)), nrow(f0), ncol(f0))
  }
  y_pred <- f0 + sigma_draws * z
  predictive_metrics(y_pred, y0, level = level)
}
