## simulate.R -------------------------------------------------------------
##
## Data-generating process for the simulation studies.
##
## The DGP follows the structural model of Section 2,
##
##   X  <- mu_e + Psi' H + delta,        delta ~ N(0, Sigma_V)
##   Y  <- alpha* + gamma*' X + phi' H + eps,        eps ~ N(0, sigma_eps^2)
##
## with H a vector of q hidden confounders shared by X and Y.  Two features
## matter for the support-recovery study and are absent from the DGP used in
## the submitted manuscript:
##
##  1. `gamma*` has exact zeros.  The submitted simulation drew every slope
##     from a fixed all-nonzero vector, so a false discovery was impossible by
##     construction and the false-discovery rate was not identifiable from the
##     output.  Here `s0` of the `p` slopes are nonzero and the remaining
##     `p - s0` are exactly zero, which is what makes `pa(Y)` a nontrivial
##     estimand.
##
##  2. The confounder loads on the null covariates too.  `Psi` is dense, so
##     every covariate -- including the ones with `gamma*_j = 0` -- is
##     marginally associated with `Y` through `H`.  A procedure that ignores
##     confounding is therefore expected to make false discoveries, which is
##     precisely the regime in which the selection rule of Section 2.1 needs
##     to be evaluated.
##
## `identifiability` controls the geometry of the environment means, which is
## what drives identification of (gamma*, K*).  Under "strong" the mu_e span
## R^p comfortably; under "weak" they are close to a common one-dimensional
## ray, so the design approaches the boundary of the identification condition
## discussed in Section 2.3.

#' Draw environment mean vectors with a controlled degree of collinearity.
#'
#' @param n_env Number of environments.
#' @param p Covariate dimension.
#' @param identifiability "strong" (means well spread over R^p) or "weak"
#'   (means nearly collinear).
#' @param weak_ratio Fraction of the mean variation that is off the common ray
#'   under "weak"; smaller values approach non-identifiability.
#' @return A `n_env x p` matrix, one mean vector per row.
draw_env_means <- function(n_env, p,
                           identifiability = c("strong", "weak"),
                           weak_ratio = 0.05) {
  identifiability <- match.arg(identifiability)
  base <- seq(-1, 1, length.out = p)

  if (identifiability == "strong") {
    ## Isotropic perturbations: the resulting means span R^p with a
    ## well-conditioned Gram matrix.
    shift <- matrix(stats::runif(n_env * p, -1.5, 1.5), nrow = n_env)
  } else {
    ## Nearly collinear: a shared direction with a small isotropic residual.
    direction <- stats::rnorm(p)
    direction <- direction / sqrt(sum(direction^2))
    magnitude <- stats::runif(n_env, -1.5, 1.5)
    shift <- outer(magnitude, direction) +
      weak_ratio * matrix(stats::rnorm(n_env * p), nrow = n_env)
  }

  sweep(shift, 2, base, "+")
}

#' Simulate one training/target data set from the hidden-confounding SCM.
#'
#' @param n_e Observations per training environment.
#' @param p Covariate dimension.
#' @param s0 Number of true parents, i.e. nonzero entries of `gamma*`.
#' @param n_env Number of training environments.  Identification of
#'   `(gamma*, K*)` needs the environment means to span R^p, hence
#'   `n_env >= p + 1`.
#' @param q Number of hidden confounders.
#' @param confounding Scalar multiplying the confounder loadings on `Y`.
#'   `confounding = 0` removes hidden confounding entirely.
#' @param gamma_signal Magnitude of the nonzero slopes.
#' @param n0 Target-domain sample size.
#' @param sigma_heterogeneity Degree to which the environment covariances
#'   differ. `0` makes every `Sigma_e` identical, which is what the earlier
#'   version of this function did and what makes the slope and covariance
#'   parameterisations of the model coincide. Positive values draw the
#'   idiosyncratic covariate covariance per environment from an inverse-Wishart
#'   centred on the common one, so `Sigma_e` genuinely varies while `Psi` and
#'   `phi` -- and therefore `K = Psi' phi` -- stay fixed. That is the
#'   configuration inner-product invariance actually describes, and the only
#'   one in which `b_e = Sigma_e^{-1} K` varies across environments.
#' @param identifiability Passed to `draw_env_means()`.
#' @param target_shift Extra dispersion of the target-domain mean relative to
#'   the training means.
#' @return A list with the training data (`x`, `y`, `z`), the target data
#'   (`x0`, `y0`), the truth (`gamma`, `alpha`, `parents`, `k`, `sigma_y`) and
#'   the environment means.
simulate_gi_data <- function(n_e = 200,
                             p = 6,
                             s0 = 3,
                             n_env = p + 1,
                             q = 3,
                             confounding = 1,
                             gamma_signal = 1,
                             n0 = 500,
                             sigma_heterogeneity = 0,
                             identifiability = c("strong", "weak"),
                             target_shift = 2) {
  identifiability <- match.arg(identifiability)
  stopifnot(s0 >= 0, s0 <= p, n_env >= 1, q >= 1)

  ## ---- Truth ----------------------------------------------------------
  parents <- if (s0 > 0) sort(sample.int(p, s0)) else integer(0)
  gamma_true <- numeric(p)
  if (s0 > 0) {
    ## Alternating signs so that the marginal association of a parent is not
    ## systematically aligned with the confounding bias.
    gamma_true[parents] <- gamma_signal * rep_len(c(1, -1), s0)
  }
  alpha_true <- 0.5

  ## Loadings of the hidden confounders.  Psi is dense, so H perturbs the null
  ## covariates as well as the parents.
  psi <- matrix(stats::rnorm(q * p), nrow = q, ncol = p)
  phi <- confounding * stats::rnorm(q)
  sigma_eps <- 1

  ## Idiosyncratic covariate noise: compound symmetry, as in Section 3.1.2.
  sigma_v <- 0.5 * matrix(1, p, p) + 0.5 * diag(p)

  ## ---- Population quantities implied by the truth ----------------------
  ## eps_Y = phi' H + eps, and X - mu_e = Psi' H + delta, so
  ##   K*      = Cov(eps_Y, X)  = phi' Psi          (a p-vector)
  ##   sigma_y^2 = Var(eps_Y)   = ||phi||^2 + sigma_eps^2
  ## and Sigma_e = Psi' Psi + Sigma_V for every environment.
  k_true <- as.vector(crossprod(phi, psi))
  sigma_y_true <- sqrt(sum(phi^2) + sigma_eps^2)
  sigma_x <- crossprod(psi) + sigma_v

  env_means <- draw_env_means(n_env, p, identifiability)

  ## ---- Training data ---------------------------------------------------
  draw_block <- function(n, mean_vec, cov_v) {
    h <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
    delta <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = cov_v)
    x <- sweep(h %*% psi + delta, 2, mean_vec, "+")
    y <- as.vector(alpha_true + x %*% gamma_true + h %*% phi +
                     stats::rnorm(n, sd = sigma_eps))
    list(x = x, y = y)
  }

  ## Environment-specific idiosyncratic covariance.  Psi and phi are held
  ## fixed so that K = Psi' phi remains common across environments, which is
  ## what inner-product invariance asserts; only Sigma_e moves.
  draw_sigma_v <- function() {
    if (sigma_heterogeneity <= 0) {
      return(sigma_v)
    }
    nu <- p + 1 + 1 / sigma_heterogeneity
    scale_mat <- (nu - p - 1) * sigma_v
    m <- solve(stats::rWishart(1, nu, solve(scale_mat))[, , 1])
    enforce_pd((m + t(m)) / 2)
  }
  sigma_v_env <- lapply(seq_len(n_env), function(e) draw_sigma_v())

  blocks <- lapply(seq_len(n_env), function(e) {
    draw_block(n_e, env_means[e, ], sigma_v_env[[e]])
  })
  x <- do.call(rbind, lapply(blocks, `[[`, "x"))
  y <- unlist(lapply(blocks, `[[`, "y"), use.names = FALSE)
  z <- rep(seq_len(n_env), each = n_e)

  ## ---- Target data -----------------------------------------------------
  ## The target domain shifts both the mean and the dispersion of X, so that
  ## Sigma_0 differs from every Sigma_e.  This is what makes the distinction
  ## between the training conditional variance and S_0^2 = sigma_y^2 -
  ## K*' Sigma_0^{-1} K* visible in the predictive intervals.
  mu0_true <- seq(-1, 1, length.out = p) +
    target_shift * stats::runif(p, -1, 1)
  sigma_v0 <- sigma_v + 0.5 * diag(p)
  target <- draw_block(n0, mu0_true, sigma_v0)

  list(
    x = x,
    y = y,
    z = z,
    x0 = target$x,
    y0 = target$y,
    truth = list(
      alpha = alpha_true,
      gamma = gamma_true,
      parents = parents,
      k = k_true,
      sigma_y = sigma_y_true,
      sigma_x = sigma_x,
      sigma_x0 = crossprod(psi) + sigma_v0,
      psi = psi,
      phi = phi,
      ## Per-environment covariances and the induced slopes, so a test can
      ## check directly how much b_e = Sigma_e^{-1} K actually varies.
      sigma_e = lapply(sigma_v_env, function(sv) crossprod(psi) + sv),
      b_e = lapply(sigma_v_env, function(sv)
        as.vector(solve(crossprod(psi) + sv, k_true)))
    ),
    env_means = env_means,
    mu0 = mu0_true,
    settings = list(
      n_e = n_e, p = p, s0 = s0, n_env = n_env, q = q,
      confounding = confounding, gamma_signal = gamma_signal, n0 = n0,
      sigma_heterogeneity = sigma_heterogeneity,
      identifiability = identifiability
    )
  )
}
