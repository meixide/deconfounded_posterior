## covariance.R -----------------------------------------------------------
##
## Plug-in estimators for the environment covariance matrices Sigma_e and for
## the target-domain Sigma_0.
##
## The BGI predictive distribution involves Sigma_0^{-1} explicitly and
## Sigma_e^{-1} through the training likelihood, so the whole procedure hinges
## on inverting covariance matrices.  With moderate p, or with environments
## that contribute few observations, the plain sample covariance is
## ill-conditioned and its inverse distorts both the predictive mean and the
## predictive variance.  These estimators let a user trade a little bias for a
## large reduction in that distortion, and let the simulations quantify the
## trade-off rather than assume it away.
##
## All estimators return a symmetric positive definite matrix and carry the
## realised shrinkage intensity in attribute "shrinkage".

#' Mean squared Frobenius deviation of the per-observation outer products from
#' their mean -- the quantity Ledoit-Wolf calls b-bar squared.
#'
#' Computed in closed form rather than by looping over observations.  Since
#'
#'   || x x' - S ||_F^2  =  (x'x)^2 - 2 x' S x + ||S||_F^2,
#'
#' the sum over `i` needs only row sums and one matrix product, which turns an
#' O(n) R-level loop with a p x p outer product per observation into two BLAS
#' calls.  On the case-study data (500,000 rows across 52 environments) the
#' loop version dominated the entire fit.
#'
#' @param xc Centred data matrix.
#' @param s_ml The 1/n scaled sample covariance of `xc`.
#' @keywords internal
mean_outer_product_variance <- function(xc, s_ml) {
  n <- nrow(xc)
  p <- ncol(xc)
  q <- rowSums(xc * xc)                       # x_i' x_i
  xsx <- rowSums((xc %*% s_ml) * xc)          # x_i' S x_i
  total <- sum(q * q) - 2 * sum(xsx) + n * sum(s_ml * s_ml)
  max(total, 0) / (n^2 * p)
}

#' Ledoit-Wolf linear shrinkage towards a scaled identity.
#'
#' Implements the well-conditioned estimator of Ledoit and Wolf (2004),
#' `(1 - rho) * S + rho * m * I`, with `m = tr(S)/p` and `rho` chosen to
#' minimise the expected squared Frobenius distance to the population
#' covariance.
#'
#' @param x Numeric matrix, observations in rows.
#' @return A `p x p` symmetric positive definite matrix.
shrink_ledoit_wolf <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  xc <- sweep(x, 2, colMeans(x), "-")

  s <- crossprod(xc) / n
  m <- sum(diag(s)) / p
  d2 <- sum((s - m * diag(p))^2) / p

  ## b2 estimates E||S - Sigma||_F^2 / p and is capped at d2, as in the
  ## original construction, so that rho stays in [0, 1].
  b2 <- min(mean_outer_product_variance(xc, s), d2)

  rho <- if (d2 > 0) b2 / d2 else 0
  out <- (1 - rho) * s + rho * m * diag(p)
  ## The sample covariance above uses the 1/n scaling; rescale to 1/(n-1) so
  ## that the estimator reduces to stats::cov() when rho = 0.
  out <- out * n / max(n - 1, 1)
  attr(out, "shrinkage") <- rho
  out
}

#' Oracle-approximating shrinkage (OAS).
#'
#' Chen, Wiesel, Eldar and Hero (2010).  Behaves better than Ledoit-Wolf when
#' the data really are Gaussian, which is the modelling assumption here.
#'
#' @inheritParams shrink_ledoit_wolf
shrink_oas <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  xc <- sweep(x, 2, colMeans(x), "-")

  s <- crossprod(xc) / n
  m <- sum(diag(s)) / p
  tr_s2 <- sum(s * s)
  tr2_s <- sum(diag(s))^2

  num <- (1 - 2 / p) * tr_s2 + tr2_s
  den <- (n + 1 - 2 / p) * (tr_s2 - tr2_s / p)
  rho <- if (den > 0) min(num / den, 1) else 1

  out <- (1 - rho) * s + rho * m * diag(p)
  out <- out * n / max(n - 1, 1)
  attr(out, "shrinkage") <- rho
  out
}

#' Ridge-regularised sample covariance.
#'
#' `S + lambda * I`, with `lambda` expressed as a fraction of the average
#' variance so that the estimator is scale-equivariant.
#'
#' @inheritParams shrink_ledoit_wolf
#' @param lambda Ridge intensity, relative to `tr(S)/p`.
shrink_ridge <- function(x, lambda = 0.01) {
  x <- as.matrix(x)
  p <- ncol(x)
  s <- stats::cov(x)
  out <- s + lambda * (sum(diag(s)) / p) * diag(p)
  attr(out, "shrinkage") <- lambda
  out
}

#' Ledoit-Wolf linear shrinkage towards an arbitrary target.
#'
#' Same construction as `shrink_ledoit_wolf()` but shrinking towards a
#' user-supplied matrix rather than a scaled identity.  The environment
#' covariances are shrunk towards the pooled within-environment covariance,
#' which is the target that matters here.
#'
#' Why this target.  The predictive variance in environment `e` is
#' `sigma_y^2 - K' Sigma_e^{-1} K`, so independent noise in each `Sigma_hat_e`
#' translates into noise in `E` separate variance constraints.  Under strong
#' confounding `K' Sigma_e^{-1} K` sits close to `sigma_y^2`, and the *most
#' extreme* of the `E` noisy estimates then drives the fit: the model is
#' forced to shrink `K` merely to keep one environment's conditional variance
#' positive, which biases `gamma` through `gamma = (within-environment slope)
#' - Sigma^{-1} K`.  Shrinking the `Sigma_hat_e` towards their pooled value
#' removes that artefact while still letting genuine heterogeneity through
#' when the environments are large enough to evidence it.
#'
#' @param x Numeric matrix, observations in rows.
#' @param target Symmetric positive definite shrinkage target.
#' @return A symmetric matrix with attribute "shrinkage".
shrink_to_target <- function(x, target) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  xc <- sweep(x, 2, colMeans(x), "-")

  ## Everything is computed on the unbiased 1/(n - 1) scale so that the
  ## dispersion estimate d2 and the returned matrix use the same convention
  ## as the target.
  s_ml <- crossprod(xc) / n
  s <- s_ml * n / max(n - 1, 1)
  d2 <- sum((s - target)^2) / p
  b2 <- min(mean_outer_product_variance(xc, s_ml), d2)

  rho <- if (d2 > 0) b2 / d2 else 1
  out <- (1 - rho) * s + rho * target
  attr(out, "shrinkage") <- rho
  out
}

#' Pooled within-environment covariance.
#'
#' @param x Numeric matrix, observations in rows.
#' @param z Environment labels.
#' @return The degrees-of-freedom weighted average of the within-environment
#'   sample covariances.
pooled_covariance <- function(x, z) {
  x <- as.matrix(x)
  envs <- sort(unique(z))
  acc <- matrix(0, ncol(x), ncol(x))
  df <- 0
  for (e in envs) {
    xe <- x[z == e, , drop = FALSE]
    ne <- nrow(xe)
    if (ne < 2L) next
    xc <- sweep(xe, 2, colMeans(xe), "-")
    acc <- acc + crossprod(xc)
    df <- df + ne - 1L
  }
  if (df <= 0) {
    stop("Cannot pool covariances: no environment has two observations.",
         call. = FALSE)
  }
  enforce_pd((acc / df + t(acc / df)) / 2, shrinkage = NA_real_)
}

#' Estimate a covariance matrix with the requested regularisation.
#'
#' @param x Numeric matrix, observations in rows.
#' @param method One of "sample", "ledoit_wolf", "oas", "ridge", "pooled".
#'   "pooled" requires `target` and shrinks towards it.
#' @param target Shrinkage target used when `method = "pooled"`.
#' @param lambda Ridge intensity, used only when `method = "ridge"`.
#' @param jitter Minimum eigenvalue enforced on the result, relative to the
#'   average variance.  Guards against exact singularity when `method =
#'   "sample"` and `n <= p`.
#' @return A symmetric positive definite matrix with attributes "shrinkage"
#'   and "condition" (the 2-norm condition number).
estimate_covariance <- function(x,
                                method = c("pooled", "ledoit_wolf", "sample",
                                           "oas", "ridge"),
                                target = NULL,
                                lambda = 0.01,
                                jitter = 1e-8) {
  method <- match.arg(method)
  x <- as.matrix(x)
  if (nrow(x) < 2L) {
    stop("estimate_covariance() needs at least two observations.", call. = FALSE)
  }
  if (method == "pooled" && is.null(target)) {
    stop('estimate_covariance(method = "pooled") needs a shrinkage target.',
         call. = FALSE)
  }

  out <- switch(
    method,
    sample      = {
      s <- stats::cov(x)
      attr(s, "shrinkage") <- 0
      s
    },
    ledoit_wolf = shrink_ledoit_wolf(x),
    oas         = shrink_oas(x),
    ridge       = shrink_ridge(x, lambda = lambda),
    pooled      = shrink_to_target(x, target)
  )

  out <- (out + t(out)) / 2
  out <- enforce_pd(out, jitter = jitter,
                    shrinkage = attr(out, "shrinkage"))
  out
}

#' Nudge a symmetric matrix to be numerically positive definite.
#'
#' @param s Symmetric matrix.
#' @param jitter Minimum eigenvalue, relative to the average variance.
#' @param shrinkage Value to carry through in the "shrinkage" attribute.
enforce_pd <- function(s, jitter = 1e-8, shrinkage = NA_real_) {
  p <- ncol(s)
  scale <- sum(diag(s)) / p
  ev <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
  floor_ev <- jitter * scale
  if (min(ev) < floor_ev) {
    s <- s + (floor_ev - min(ev)) * diag(p)
    ev <- ev + (floor_ev - min(ev))
  }
  attr(s, "shrinkage") <- shrinkage
  attr(s, "condition") <- max(ev) / min(ev)
  s
}

#' Lower Cholesky factor of the scatter matrix `S = sum_i (x_i - xbar)(...)'`.
#'
#' Returned as `C` with `S = C C'`, which is the form the full-covariance
#' model consumes.  A small ridge keeps `C` well defined when `n <= p`, where
#' the scatter matrix is singular; the model's inverse-Wishart prior is what
#' actually regularises in that regime, so the ridge only has to avoid a
#' numerical failure.
#'
#' @param x Numeric matrix, observations in rows.
#' @param jitter Ridge added to the diagonal, relative to the average variance.
scatter_chol <- function(x, jitter = 1e-6) {
  x <- as.matrix(x)
  p <- ncol(x)
  xc <- sweep(x, 2, colMeans(x), "-")
  s <- crossprod(xc)
  s <- (s + t(s)) / 2
  ridge <- jitter * max(sum(diag(s)) / p, .Machine$double.eps)
  t(chol(s + ridge * diag(p)))
}

#' Lower Cholesky factor, with a clear error if the matrix is not usable.
#'
#' @param s Symmetric positive definite matrix.
#' @param what Label used in the error message.
lower_chol <- function(s, what = "covariance matrix") {
  r <- tryCatch(chol(s), error = function(e) {
    stop(sprintf("%s is not positive definite: %s", what, conditionMessage(e)),
         call. = FALSE)
  })
  t(r)
}
