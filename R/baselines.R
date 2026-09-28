## baselines.R ------------------------------------------------------------
##
## Competing procedures for recovering pa(Y) and for predicting in the target
## domain.  Four families:
##
##   naive
##     `fit_ols()`              pooled OLS.  Ignores confounding and ignores
##                              the environment structure entirely.
##   instrumental / GI
##     `fit_pooled_gi()`        frequentist GI as a pooled least squares fit of
##                              Y on X and the environment means.
##     `fit_iv_2sls()`          2SLS with the environment indicators as
##                              instruments.  Identical to the above in
##                              `gamma`; see the note there.
##   robustness / distribution shift
##     `fit_anchor()`           anchor regression (Rothenhausler et al., 2021).
##     `fit_group_dro()`        worst-case environment risk (Sagawa et al.,
##                              2020).  No tuning parameter.
##     `fit_vrex()`             risk extrapolation (Krueger et al., 2021).
##     `fit_dro_wasserstein()`  Wasserstein DRO least squares (Blanchet et al.,
##                              2019), the `sinha2020` family.
##   invariance-based selection
##     `fit_icp()`              invariant causal prediction (Peters et al.,
##                              2016).
##
## Two structural points these baselines make, both of which belong in the
## manuscript.
##
## First, `fit_ols()`, `fit_anchor()` and `fit_iv_2sls()` are one path, not
## three methods: anchor at `gamma_anchor = 1` is OLS and at
## `gamma_anchor -> infinity` is 2SLS, and 2SLS equals `fit_pooled_gi()`
## exactly.  GI's `gamma` sits at the 2SLS endpoint of that path.
##
## Second, and this is the substantive one: **none of these baselines
## estimates `K`.**  They differ only in how they estimate a slope, and every
## one of them predicts with that slope alone.  Since the target-domain
## conditional mean is
##
##     E[Y_0 | X_0] = alpha + gamma' X_0 + K' Sigma_0^{-1} (X_0 - mu_0),
##
## and `gamma*` is the causal slope rather than the target-domain population
## least-squares slope, predicting with any slope alone is suboptimal by
## construction under hidden confounding.  That is why the comparison has to
## be reported on predictive quantities -- coverage and interval score -- and
## not only on how well each method recovers `gamma`.
##
## Relatedly, `fit_vrex()` and `fit_dro_wasserstein()` return `p_value = NA`.
## They are penalised point estimators with no valid Wald inference, so they
## appear in the prediction comparison and not in the selection comparison.
## That most distribution-shift procedures supply a point predictor and
## nothing else is a finding worth stating, not a gap in the code.

#' Pooled ordinary least squares.
#'
#' @param x Training covariates, `N x p`.
#' @param y Training response.
#' @param x0 Target covariates; may be `NULL`.
#' @param level Nominal level for the prediction intervals.
#' @return A list with coefficient estimates, standard errors, two-sided
#'   p-values for `gamma_j = 0`, and target-domain prediction intervals.
fit_ols <- function(x, y, x0 = NULL, level = 0.95) {
  x <- as.matrix(x)
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  fit <- stats::lm(y ~ ., data = data.frame(y = y, x))
  co <- summary(fit)$coefficients

  slope_rows <- match(colnames(x), rownames(co))
  out <- list(
    alpha = unname(stats::coef(fit)[1]),
    gamma = unname(stats::coef(fit)[colnames(x)]),
    se = unname(co[slope_rows, "Std. Error"]),
    p_value = unname(co[slope_rows, "Pr(>|t|)"]),
    sigma = summary(fit)$sigma
  )

  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    colnames(x0) <- colnames(x)
    pred <- stats::predict(fit, newdata = as.data.frame(x0),
                           interval = "prediction", level = level)
    out$pred_mean <- unname(pred[, "fit"])
    out$pred_lower <- unname(pred[, "lwr"])
    out$pred_upper <- unname(pred[, "upr"])
  }
  out
}

#' Frequentist Generative Invariance as a pooled least squares fit.
#'
#' Regresses `Y` on `X` and on the environment-specific covariate means
#' `mu_hat_{Z_i}`.  Writing the GI conditional mean as
#' `alpha + (gamma + b)' X - b' mu_e` with `b = Sigma^{-1} K`, the coefficient
#' on `X` estimates `gamma + b` and the coefficient on the means estimates
#' `-b`, so `gamma_hat` is the sum of the two blocks.  Standard errors follow
#' from the corresponding linear contrast of the joint covariance matrix.
#'
#' ---- This estimator *is* two-stage least squares ---------------------------
#'
#' `gamma_hat` from this fit equals `fit_iv_2sls()` exactly, for every `E >=
#' p + 1`, over-identified cases included.  `tests/test_baseline_identities.R`
#' checks it numerically; the argument is Frisch-Waugh-Lovell:
#'
#'   Write `X = W + M` with `M = Pi_A X` the environment means and
#'   `W = X - M` the within-environment deviations.  `W` is orthogonal to the
#'   environment indicator space, hence to both `1` and `M`.  The regression of
#'   `Y` on `[1, X, M]` spans the same space as `[1, W, M]`, and since
#'   `X = W + M`, the coefficient on `W` is `c1` while the coefficient on `M`
#'   is `c1 + c2`.  Orthogonality of `W` to `[1, M]` makes that second
#'   coefficient the OLS fit of `Y` on `[1, M]` alone.  Two-stage least squares
#'   with the environment dummies as instruments has first-stage fitted values
#'   `Pi_A X = M`, so its second stage is precisely the OLS fit of `Y` on
#'   `[1, M]`.  Hence `c1 + c2 = gamma_hat_2SLS`.
#'
#' Three consequences for the manuscript:
#'
#'   * the baseline Referee 2 proposes ("pooled OLS on `X` and the environment
#'     means") and the IV comparison they ask to see extended beyond the
#'     single-source example are the *same* baseline;
#'   * the frequentist GI estimator of `gamma` is not a new estimator; it is
#'     IV with the environment as instrument.  Section 3.1.1 already observes
#'     this in one dimension and it holds in general.  **The coincidence is
#'     with the `gamma` part only, and `gamma` alone is the wrong thing to
#'     predict with.**  The target-domain conditional mean is
#'
#'         E[Y_0 | X_0] = alpha + gamma' X_0 + K' Sigma_0^{-1} (X_0 - mu_0),
#'
#'     and `gamma*` is the *causal* slope, not the population least-squares
#'     slope in the target domain.  Under hidden confounding those differ, and
#'     `K' Sigma_0^{-1} (X_0 - mu_0)` is exactly the correction between them.
#'     IV recovers `gamma*` and stops; predicting with it discards the
#'     correction.  That is why IV is causally right and predictively wrong,
#'     which is the point Section 3.1.1 makes for the single-source example and
#'     which `iv_pred_coverage` in the simulation output now makes for the
#'     multi-source one.  `K` is intrinsic to GI and has no IV counterpart;
#'   * combined with `fit_anchor()`, whose `gamma_anchor -> infinity` limit is
#'     also 2SLS, this places GI's `gamma` exactly at the endpoint of the
#'     anchor regression path.
#'
#' Caveat, reported here rather than hidden: the environment means are
#' estimated, and the standard errors below condition on them.  Measured
#' against the 2SLS standard errors for the identical point estimate, they run
#' 0.35 to 0.76 times as large, so they are materially anti-conservative.  This
#' is one reason the Bayesian treatment, which puts the `mu_e` in the model as
#' parameters, is not merely a reparameterisation of this fit.
#'
#' @param x Training covariates, `N x p`.
#' @param y Training response.
#' @param z Environment labels.
#' @param x0 Target covariates; may be `NULL`.
#' @param level Nominal level for the prediction intervals.
#' @return A list with `gamma`, `k_slope` (the estimate of `Sigma^{-1} K`),
#'   standard errors, p-values and target-domain prediction intervals.
fit_pooled_gi <- function(x, y, z, x0 = NULL, level = 0.95) {
  x <- as.matrix(x)
  p <- ncol(x)
  z <- as.integer(as.factor(z))

  ## Environment mean attached to each observation.
  env_means <- do.call(rbind, lapply(sort(unique(z)), function(e) {
    colMeans(x[z == e, , drop = FALSE])
  }))
  m <- env_means[z, , drop = FALSE]

  design <- cbind(x, m)
  colnames(design) <- c(paste0("x", seq_len(p)), paste0("m", seq_len(p)))
  fit <- stats::lm(y ~ ., data = data.frame(y = y, design))

  theta <- stats::coef(fit)
  v <- stats::vcov(fit)
  keep <- c(paste0("x", seq_len(p)), paste0("m", seq_len(p)))
  theta <- theta[keep]
  v <- v[keep, keep, drop = FALSE]

  ## gamma = (coef on X) + (coef on M), a linear contrast A theta.
  a <- cbind(diag(p), diag(p))
  gamma_hat <- as.vector(a %*% theta)
  gamma_var <- diag(a %*% v %*% t(a))
  se <- sqrt(pmax(gamma_var, 0))
  df <- stats::df.residual(fit)
  p_value <- 2 * stats::pt(abs(gamma_hat / se), df = df, lower.tail = FALSE)

  b_hat <- -as.vector(theta[paste0("m", seq_len(p))])

  out <- list(
    alpha = unname(stats::coef(fit)[1]),
    gamma = gamma_hat,
    se = se,
    p_value = p_value,
    k_slope = b_hat,
    sigma = summary(fit)$sigma
  )

  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    ## Target-domain conditional mean: alpha + gamma' X0 + b0' (X0 - mu0),
    ## with b0 = Sigma_0^{-1} K and K = Sigma_train b.  The pooled fit has no
    ## notion of Sigma_0, so the honest plug-in is
    ## b0 = Sigma_0^{-1} Sigma_train b.
    sigma_train <- stats::cov(x)
    sigma_0 <- stats::cov(x0)
    b0 <- tryCatch(
      as.vector(solve(sigma_0, sigma_train %*% b_hat)),
      error = function(e) b_hat
    )
    mu0 <- colMeans(x0)
    centred <- sweep(x0, 2, mu0, "-")
    out$pred_mean <- as.vector(out$alpha + x0 %*% gamma_hat + centred %*% b0)
    ## The pooled fit only has the training residual scale to offer.
    halfwidth <- stats::qnorm(1 - (1 - level) / 2) * out$sigma
    out$pred_lower <- out$pred_mean - halfwidth
    out$pred_upper <- out$pred_mean + halfwidth
  }
  out
}

#' Project onto the anchor space spanned by the environment indicators.
#'
#' For environment dummies `A`, the projection `Pi_A v` replaces each entry by
#' the mean of `v` within its environment.  Used by both `fit_iv_2sls()` and
#' `fit_anchor()`, which differ only in how much weight they put on this
#' component.
#'
#' @param v Numeric vector or matrix, observations in rows.
#' @param z Environment labels.
#' @return An object of the same shape as `v`.
anchor_projection <- function(v, z) {
  v <- as.matrix(v)
  out <- v
  for (e in unique(z)) {
    idx <- z == e
    out[idx, ] <- rep(colMeans(v[idx, , drop = FALSE]), each = sum(idx))
  }
  out
}

#' Two-stage least squares using the environment indicators as instruments.
#'
#' The multi-source counterpart of the IV comparison that the manuscript
#' currently shows only for the single-source example.  The environment label
#' is a valid instrument under the structural model of Section 2: it shifts
#' `X` but affects `Y` only through `X`.
#'
#' Identification needs at least as many instruments as covariates,
#' `E - 1 >= p`; at `E = p + 1` the system is exactly identified, which is the
#' weakest configuration in which IV is defined here and is also where it is
#' least stable.
#'
#' This is the `gamma -> infinity` limit of `fit_anchor()`.
#'
#' @param x Training covariates, `N x p`.
#' @param y Training response.
#' @param z Environment labels.
#' @param x0 Target covariates; may be `NULL`.
#' @param level Nominal level for the prediction intervals.
#' @return The same fields as `fit_ols()`.
fit_iv_2sls <- function(x, y, z, x0 = NULL, level = 0.95) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  z <- as.integer(as.factor(z))
  n_env <- length(unique(z))

  if (n_env - 1L < p) {
    stop(sprintf(
      "2SLS needs at least p instruments: E - 1 = %d < p = %d. The environment indicators cannot identify gamma here.",
      n_env - 1L, p), call. = FALSE)
  }

  ## First stage: the fitted values of X on the environment dummies are just
  ## the within-environment means.
  xhat <- anchor_projection(x, z)
  d_hat <- cbind(1, xhat)
  d_raw <- cbind(1, x)

  xtx <- crossprod(d_hat)
  coef <- tryCatch(as.vector(solve(xtx, crossprod(d_hat, y))),
                   error = function(e) {
                     ## Structurally this happens when the environment
                     ## means fail to span R^p.  Numerically it also happens
                     ## when the covariates differ wildly in scale, which is
                     ## easy to mistake for the former.
                     stop("2SLS normal equations are singular: the ",
                          "environment means are collinear, or the ",
                          "covariates differ too much in scale (largest sd / ",
                          "smallest sd = ",
                          signif(max(apply(x, 2, stats::sd)) /
                                   max(min(apply(x, 2, stats::sd)), 1e-300), 3),
                          ").", call. = FALSE)
                   })

  ## Residuals use the *original* covariates, as 2SLS requires.
  resid <- as.vector(y - d_raw %*% coef)
  df <- n - (p + 1L)
  s2 <- sum(resid^2) / df
  vcov <- s2 * solve(xtx)
  se <- sqrt(pmax(diag(vcov)[-1], 0))
  gamma_hat <- coef[-1]

  out <- list(
    alpha = coef[1],
    gamma = gamma_hat,
    se = se,
    p_value = 2 * stats::pt(abs(gamma_hat / se), df = df, lower.tail = FALSE),
    sigma = sqrt(s2)
  )

  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    out$pred_mean <- as.vector(out$alpha + x0 %*% gamma_hat)
    ## IV has no notion of the target covariance, so the only scale it can
    ## offer is the training residual one.  That is the point of the
    ## comparison, not an oversight.
    half <- stats::qnorm(1 - (1 - level) / 2) * out$sigma
    out$pred_lower <- out$pred_mean - half
    out$pred_upper <- out$pred_mean + half
  }
  out
}

#' Anchor regression (Rothenhausler, Meinshausen, Buhlmann and Peters, 2021).
#'
#' Minimises
#'   `|| (I - Pi_A)(Y - a - X b) ||^2 + gamma_anchor * || Pi_A (Y - a - X b) ||^2`
#' with `Pi_A` the projection onto the environment indicators.  Equivalently,
#' ordinary least squares after applying `T = (I - Pi_A) + sqrt(gamma) Pi_A` to
#' both sides.  `gamma_anchor = 1` recovers OLS and `gamma_anchor -> infinity`
#' recovers `fit_iv_2sls()`; intermediate values interpolate.
#'
#' This is the comparison the manuscript's Introduction sets up but never
#' runs.  Three claims made there are directly testable against it: that the
#' anchor regularisation "precludes asymptotically unbiased estimation of
#' causal parameters"; that GI departs from "predictive robustness as a
#' minimax optimization problem"; and that GI works "without requiring
#' hyperparameter tuning".  The third is the sharpest, because `gamma_anchor`
#' has no data-driven default when the target domain is unlabelled — which is
#' the setting the paper is about.
#'
#' @inheritParams fit_iv_2sls
#' @param gamma_anchor Regularisation strength; `1` is OLS.
#' @return The same fields as `fit_ols()`, plus `gamma_anchor`.
fit_anchor <- function(x, y, z, x0 = NULL, gamma_anchor = 8, level = 0.95) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(gamma_anchor > 0)

  w <- sqrt(gamma_anchor)
  pa_x <- anchor_projection(x, z)
  pa_y <- anchor_projection(matrix(y, ncol = 1), z)

  tx <- (x - pa_x) + w * pa_x
  ty <- as.vector((as.vector(y) - as.vector(pa_y)) + w * as.vector(pa_y))
  ## The intercept column lies in the anchor space, so T(1) = sqrt(gamma) * 1.
  t_design <- cbind(w, tx)

  qrd <- qr(t_design)
  if (qrd$rank < ncol(t_design)) {
    stop("Anchor design is rank deficient at gamma_anchor = ", gamma_anchor,
         call. = FALSE)
  }
  coef <- as.vector(qr.coef(qrd, ty))

  ## Report fit quality on the original scale, which is what predictions use.
  resid <- as.vector(y - cbind(1, x) %*% coef)
  df <- n - (p + 1L)
  s2 <- sum(resid^2) / df
  ## Sandwich-free standard errors on the transformed design; adequate for a
  ## baseline comparison and consistent with how anchor regression is usually
  ## reported.
  vcov <- s2 * chol2inv(qr.R(qrd))
  se <- sqrt(pmax(diag(vcov)[-1], 0))
  gamma_hat <- coef[-1]

  out <- list(
    alpha = coef[1],
    gamma = gamma_hat,
    se = se,
    p_value = 2 * stats::pt(abs(gamma_hat / se), df = df, lower.tail = FALSE),
    sigma = sqrt(s2),
    gamma_anchor = gamma_anchor
  )

  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    out$pred_mean <- as.vector(out$alpha + x0 %*% gamma_hat)
    half <- stats::qnorm(1 - (1 - level) / 2) * out$sigma
    out$pred_lower <- out$pred_mean - half
    out$pred_upper <- out$pred_mean + half
  }
  out
}

#' Design matrix and per-environment risks, shared by the robustness methods.
#'
#' @keywords internal
env_risks <- function(coef, d, y, z, envs) {
  r2 <- as.vector(y - d %*% coef)^2
  vapply(envs, function(e) mean(r2[z == e]), numeric(1))
}

#' Wrap a fitted coefficient vector in the standard baseline return shape.
#'
#' `p_value` is `NA` for the penalised and minimax estimators: they come with
#' no inferential statement, which is itself worth reporting.  A method that
#' returns only a point predictor cannot be compared on selection, only on
#' prediction.
#'
#' @keywords internal
robust_fit_output <- function(coef, x, y, x0, level, se = NULL) {
  p <- ncol(x)
  resid <- as.vector(y - cbind(1, x) %*% coef)
  sigma <- sqrt(sum(resid^2) / max(nrow(x) - p - 1L, 1L))
  gamma_hat <- coef[-1]

  out <- list(
    alpha = coef[1],
    gamma = gamma_hat,
    se = if (is.null(se)) rep(NA_real_, p) else se,
    p_value = if (is.null(se)) rep(NA_real_, p) else
      2 * stats::pnorm(abs(gamma_hat / se), lower.tail = FALSE),
    sigma = sigma
  )
  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    out$pred_mean <- as.vector(coef[1] + x0 %*% gamma_hat)
    half <- stats::qnorm(1 - (1 - level) / 2) * sigma
    out$pred_lower <- out$pred_mean - half
    out$pred_upper <- out$pred_mean + half
  }
  out
}

#' Group distributionally robust optimisation over environments.
#'
#' Minimises the worst-case environment risk, `min_beta max_e R_e(beta)`
#' (Sagawa, Koh, Hashimoto and Liang, 2020).  This is the concrete
#' instantiation of "predictive robustness as a minimax optimization problem"
#' that the manuscript's Introduction says GI departs from, so it belongs in
#' the comparison: the claim should be tested, not merely asserted.
#'
#' Unlike anchor regression, V-REx and Wasserstein DRO it has **no tuning
#' parameter**, which makes it the cleanest minimax comparator for a paper
#' whose selling point is that no tuning is required.
#'
#' Solved by exponentiated gradient ascent on the environment weights with a
#' weighted least squares step in between; the iterate attaining the smallest
#' worst-case risk is returned.
#'
#' @inheritParams fit_iv_2sls
#' @param n_iter Number of primal-dual iterations.
#' @param eta Step size for the weight update.
#' @return The same fields as `fit_ols()`; `se` and `p_value` come from the
#'   final weighted least squares fit and ignore the weight selection, so they
#'   are approximate.
fit_group_dro <- function(x, y, z, x0 = NULL, level = 0.95,
                          n_iter = 200L, eta = 0.5) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  z <- as.integer(as.factor(z))
  envs <- sort(unique(z))
  n_env <- length(envs)
  d <- cbind(1, x)

  w <- rep(1 / n_env, n_env)
  best_coef <- NULL
  best_val <- Inf
  best_ow <- NULL

  for (it in seq_len(n_iter)) {
    ## Observation weights: environment weight spread over its observations.
    ow <- numeric(n)
    for (j in seq_along(envs)) {
      idx <- z == envs[j]
      ow[idx] <- w[j] / sum(idx)
    }
    sw <- sqrt(ow)
    coef <- tryCatch(as.vector(qr.coef(qr(d * sw), y * sw)),
                     error = function(e) rep(NA_real_, p + 1L))
    if (anyNA(coef)) break

    risks <- env_risks(coef, d, y, z, envs)
    val <- max(risks)
    if (val < best_val) {
      best_val <- val
      best_coef <- coef
      best_ow <- ow
    }
    ## Scale-free update so eta does not depend on the units of Y.
    w <- w * exp(eta * risks / max(risks))
    w <- w / sum(w)
  }

  if (is.null(best_coef)) {
    stop("Group DRO failed to produce a usable iterate.", call. = FALSE)
  }

  ## Approximate weighted least squares standard errors at the final iterate.
  ## The optimisation weights sum to one; rescale them to mean one before
  ## forming the variance, so that D' W D is on the scale of D' D and the
  ## usual weighted least squares formula applies.  Without the rescaling the
  ## standard errors come out inflated by a factor of order n and the method
  ## never selects anything.
  se <- tryCatch({
    w_scaled <- best_ow * n / sum(best_ow)
    sw <- sqrt(w_scaled)
    xtx_inv <- chol2inv(qr.R(qr(d * sw)))
    resid <- as.vector(y - d %*% best_coef)
    s2 <- sum(w_scaled * resid^2) / max(n - p - 1L, 1L)
    sqrt(pmax(diag(s2 * xtx_inv)[-1], 0))
  }, error = function(e) NULL)

  robust_fit_output(best_coef, x, y, x0, level, se = se)
}

#' Risk extrapolation (V-REx).
#'
#' Minimises `sum_e R_e(beta) + lambda * Var_e(R_e(beta))` (Krueger et al.,
#' 2021): empirical risk plus a penalty on how much the risk varies across
#' environments.  This is the same family as the probabilistic-prediction
#' method of Henzi et al. cited in the Introduction, which penalises the
#' variance of a scoring rule across environments rather than of the squared
#' error.
#'
#' @inheritParams fit_iv_2sls
#' @param lambda Penalty on the across-environment variance of the risks.
#' @return The same fields as `fit_ols()`, with `p_value` all `NA`: a
#'   penalised point estimator carries no valid Wald inference.
fit_vrex <- function(x, y, z, x0 = NULL, level = 0.95, lambda = 10) {
  x <- as.matrix(x)
  z <- as.integer(as.factor(z))
  envs <- sort(unique(z))
  d <- cbind(1, x)

  objective <- function(coef) {
    risks <- env_risks(coef, d, y, z, envs)
    sum(risks) + lambda * stats::var(risks)
  }
  start <- as.vector(qr.coef(qr(d), y))
  start[is.na(start)] <- 0
  opt <- stats::optim(start, objective, method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-10))
  robust_fit_output(opt$par, x, y, x0, level)
}

#' Wasserstein distributionally robust linear regression.
#'
#' The Wasserstein-DRO least squares problem admits the closed-form reduction
#' `min_beta sqrt(MSE(beta)) + sqrt(delta) * ||beta||_2` (Blanchet, Kang and
#' Murthy, 2019), i.e. a square-root ridge.  This is the distributionally
#' robust optimisation family cited as `sinha2020` in the Introduction.
#'
#' Note it uses no environment structure at all: Wasserstein DRO hedges
#' against *any* nearby distribution rather than against the observed
#' heterogeneity, which is exactly the conservatism the manuscript argues
#' against. Including it makes that argument testable.
#'
#' @inheritParams fit_iv_2sls
#' @param delta Radius of the Wasserstein ambiguity set.
#' @return The same fields as `fit_ols()`, with `p_value` all `NA`.
fit_dro_wasserstein <- function(x, y, x0 = NULL, level = 0.95, delta = 0.05) {
  x <- as.matrix(x)
  d <- cbind(1, x)

  objective <- function(coef) {
    sqrt(mean(as.vector(y - d %*% coef)^2)) +
      sqrt(delta) * sqrt(sum(coef[-1]^2))
  }
  start <- as.vector(qr.coef(qr(d), y))
  start[is.na(start)] <- 0
  opt <- stats::optim(start, objective, method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-10))
  robust_fit_output(opt$par, x, y, x0, level)
}

#' Invariant causal prediction (subset-intersection version).
#'
#' For every subset `S` of the covariates, tests whether the regression of `Y`
#' on `X_S` is invariant across environments, by comparing the residuals
#' obtained in each environment with those obtained in the remaining
#' environments (a two-sample test of location and of scale, Bonferroni
#' corrected over the two tests and over environments).  The estimate is the
#' intersection of all non-rejected subsets, which yields
#' `P(S_hat subset of pa(Y)) >= 1 - alpha` when the invariance assumption
#' holds.
#'
#' Under hidden confounding the invariance assumption fails, and this
#' procedure is expected to reject every subset and return the empty set.
#' That behaviour is the point of including it.
#'
#' @param x Training covariates, `N x p`.
#' @param y Training response.
#' @param z Environment labels.
#' @param alpha Nominal level.
#' @param max_p Refuse to run when `p` exceeds this, since the subset search
#'   is exponential.
#' @return A list with the selected set and the number of accepted subsets.
fit_icp <- function(x, y, z, alpha = 0.05, max_p = 12) {
  x <- as.matrix(x)
  p <- ncol(x)
  if (p > max_p) {
    stop("fit_icp() enumerates 2^p subsets; p = ", p, " is too large.",
         call. = FALSE)
  }
  z <- as.integer(as.factor(z))
  envs <- sort(unique(z))
  colnames(x) <- paste0("x", seq_len(p))

  test_subset <- function(s) {
    ## Leave-one-environment-out residuals.  Column names are set explicitly:
    ## letting data.frame() derive them from the subsetting expression
    ## produces names that predict() then cannot match.
    nm <- colnames(x)[s]
    p_env <- vapply(envs, function(e) {
      inside <- z == e
      if (sum(inside) < 3L || sum(!inside) < (length(s) + 3L)) {
        return(1)
      }
      d_out <- data.frame(y = y[!inside])
      new_in <- data.frame(row.names = seq_len(sum(inside)))
      if (length(s) > 0L) {
        d_out[nm] <- x[!inside, s, drop = FALSE]
        new_in[nm] <- x[inside, s, drop = FALSE]
      }
      fit <- stats::lm(y ~ ., data = d_out)
      r_in <- y[inside] - as.vector(stats::predict(fit, newdata = new_in))
      r_out <- stats::residuals(fit)

      p_mean <- tryCatch(stats::t.test(r_in, r_out)$p.value,
                         error = function(e) 1)
      p_var <- tryCatch(stats::var.test(r_in, r_out)$p.value,
                        error = function(e) 1)
      min(1, 2 * min(p_mean, p_var))
    }, numeric(1))

    min(1, length(envs) * min(p_env))
  }

  accepted <- NULL
  n_accepted <- 0L
  for (size in 0:p) {
    combos <- utils::combn(p, size, simplify = FALSE)
    if (size == 0L) combos <- list(integer(0))
    for (s in combos) {
      if (test_subset(s) > alpha) {
        n_accepted <- n_accepted + 1L
        accepted <- if (is.null(accepted)) s else intersect(accepted, s)
        if (length(accepted) == 0L) {
          ## The intersection can only shrink, so we can stop early.
          return(list(selected = integer(0), n_accepted = NA_integer_,
                      stopped_early = TRUE))
        }
      }
    }
  }

  list(
    selected = if (is.null(accepted)) integer(0) else sort(as.integer(accepted)),
    n_accepted = n_accepted,
    stopped_early = FALSE
  )
}
