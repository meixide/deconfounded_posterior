## fit_bgi.R --------------------------------------------------------------
##
## User-facing wrapper around the Stan models in ../stan/.
##
## Responsibilities:
##   * standardise the covariates, fit on the standardised scale, and map the
##     posterior back to the original scale;
##   * supply the covariance matrices Sigma_1, ..., Sigma_E and Sigma_0,
##     either as regularised plug-ins (stan/gi_hd.stan) or as quantities the
##     model infers (stan/gi_hd_fullcov.stan);
##   * reduce the training data to per-environment sufficient statistics;
##   * return posterior draws for (alpha, gamma, K, sigma_y, S_0) and for the
##     target-domain predictive distribution.
##
## Standardisation matters here for two reasons.  The ridge prior of Section 2
## puts a common scale tau * sigma_y on every entry of gamma and K, which is
## only sensible once the covariates share a scale; and centring/scaling
## improves the conditioning of the Sigma_e that the model has to invert.
## Because the map is linear the back-transformation is exact:
##
##   X~ = (X - c) / d   =>   gamma~ = d * gamma,  K~ = K / d,
##                           alpha~ = alpha + gamma' c
##
## (element-wise in d), so posterior draws on the original scale are recovered
## without any approximation.

#' Compile a Stan model once and cache the result on disk.
#'
#' Compiling inside a parallel worker is both wasteful and racy: every worker
#' pays the compilation cost, and concurrent workers can collide on rstan's
#' auto-write cache.  Batch jobs should call this once up front (see
#' `scripts/00_compile_models.R`) and the workers should read the cache.
#'
#' @param stan_file Path to the `.stan` file.
#' @param cache_dir Directory holding the compiled `.rds`.
#' @param force Recompile even if an up-to-date cache exists.
#' @return A `stanmodel` object.
compile_bgi_model <- function(stan_file,
                              cache_dir = file.path(dirname(dirname(
                                normalizePath(stan_file))), "results",
                                "compiled"),
                              force = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(cache_dir,
                     paste0(tools::file_path_sans_ext(basename(stan_file)),
                            ".rds"))

  if (!force && file.exists(cache) &&
      file.mtime(cache) > file.mtime(stan_file)) {
    model <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(model)) {
      return(model)
    }
  }

  message("Compiling ", basename(stan_file), " ...")
  model <- rstan::stan_model(file = stan_file,
                             model_name = tools::file_path_sans_ext(
                               basename(stan_file)))
  saveRDS(model, cache)
  model
}

#' Load a previously compiled model, failing loudly if the cache is missing.
#'
#' @inheritParams compile_bgi_model
load_bgi_model <- function(stan_file,
                           cache_dir = file.path(dirname(dirname(
                             normalizePath(stan_file))), "results",
                             "compiled")) {
  cache <- file.path(cache_dir,
                     paste0(tools::file_path_sans_ext(basename(stan_file)),
                            ".rds"))
  if (!file.exists(cache)) {
    stop("No compiled model at ", cache,
         ". Run scripts/00_compile_models.R first.", call. = FALSE)
  }
  ## A stale cache is silent and produces results from the wrong model, which
  ## is the worst possible failure mode in a reproducibility package.
  if (file.exists(stan_file) && file.mtime(stan_file) > file.mtime(cache)) {
    stop(basename(stan_file), " is newer than its compiled cache. ",
         "Run: Rscript scripts/00_compile_models.R --force", call. = FALSE)
  }
  readRDS(cache)
}

#' Assemble the Stan data list for `gi_hd.stan`.
#'
#' @param x Training covariates, `N x p`.
#' @param y Training response, length `N`.
#' @param z Integer environment labels in `1:E`.
#' @param x0 Target covariates, `n0 x p`; may be `NULL`.
#' @param cov_method Covariance regulariser, see `estimate_covariance()`.
#' @param cov_lambda Ridge intensity when `cov_method = "ridge"`.
#' @param standardize Centre and scale the covariates before fitting.
#' @param eta_lkj,sd_mu_scale,sd_sigma_scale,a_tau,b_tau Prior hyperparameters.
#' @param ncp Degree of non-centring for `(alpha, gamma, K)`, in `[0, 1]`.
#'   The implied prior is the same for every value; only the sampler geometry
#'   changes.  `1` (the default) is what makes the near-collinear regime of
#'   Section 2.3 samplable at all; `0` is faster when the environment means
#'   are well separated and the likelihood is strongly informative.
#' @param v_prior_shape,v_prior_rate Inverse-gamma prior on the free part of
#'   the conditional variance.  The default `3` and `2` is the prior the
#'   calibration study of `tests/test_sbc.R` validates, with mean and standard
#'   deviation both one on a standardised response.
#'
#'   `0` selects the Jeffreys prior `p(v) propto 1/v` used in the first
#'   submission, and it is retained only to reproduce those fits.  Under it the
#'   posterior is improper whenever the target Mahalanobis term can exceed every
#'   training one, which is an open region of the parameter space for any target
#'   covariance that is not uniformly larger than the training covariances.
#'   `tests/test_posterior_propriety.R` exhibits the divergence; the sampler
#'   gives no sign of it, since chains can mix perfectly well inside a
#'   non-integrable density.  Because the inverse-gamma vanishes at the origin
#'   like `v^{-(shape+1)} exp(-rate/v)`, a positive shape removes the problem
#'   whatever the anchoring does.
#' @param sigma_known Optional known covariance: either a single `p x p`
#'   matrix used for every environment, or a list of `n_env` matrices, one per
#'   environment.  Used in place of the plug-in estimate for every environment
#'   and for the target.  Lets a calibration study separate implementation
#'   error from plug-in error, and the list form lets it do so under
#'   heterogeneous covariances.
#' @param hmu_known Optional known prior mean for the `mu_e`, in place of the
#'   pooled sample mean.
#' @param compute_log_lik Return pointwise log-likelihoods (needed by `loo`).
#' @return A list with elements `stan_data`, `scaling` and `diagnostics`.
prepare_bgi_data <- function(x, y, z, x0 = NULL,
                             cov_method = "pooled",
                             cov_lambda = 0.01,
                             standardize = TRUE,
                             eta_lkj = 2,
                             sd_mu_scale = 2.5,
                             sd_sigma_scale = 2.5,
                             sd_b_scale = 1,
                             a_tau = 0.5,
                             b_tau = 0.5,
                             ncp = 1,
                             v_prior_shape = 3,
                             v_prior_rate = 2,
                             sigma_known = NULL,
                             hmu_known = NULL,
                             compute_log_lik = FALSE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  z <- as.integer(as.factor(z))
  p <- ncol(x)
  n <- nrow(x)

  stopifnot(length(y) == n, length(z) == n)
  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    stopifnot(ncol(x0) == p)
  }

  ## ---- Standardisation -------------------------------------------------
  ## The centre and scale come from the pooled training covariates only; the
  ## target covariates are mapped with the same transformation, so nothing
  ## about the target distribution leaks into the training fit.
  if (standardize) {
    centre <- colMeans(x)
    scale_x <- apply(x, 2, stats::sd)
    scale_x[scale_x <= 0 | !is.finite(scale_x)] <- 1
  } else {
    centre <- rep(0, p)
    scale_x <- rep(1, p)
  }
  xs <- sweep(sweep(x, 2, centre, "-"), 2, scale_x, "/")
  x0s <- if (is.null(x0)) {
    NULL
  } else {
    sweep(sweep(x0, 2, centre, "-"), 2, scale_x, "/")
  }

  ## ---- Per-environment covariances and sufficient statistics -----------
  envs <- sort(unique(z))
  n_env <- length(envs)
  n_e <- integer(n_env)
  xbar <- matrix(0, n_env, p)
  xx <- array(0, dim = c(n_env, p, p))
  xy <- matrix(0, n_env, p)
  ysum <- numeric(n_env)
  yy <- numeric(n_env)
  l_sigma <- array(0, dim = c(n_env, p, p))
  chol_scatter <- array(0, dim = c(n_env, p, p))
  cond_e <- numeric(n_env)
  shrink_e <- numeric(n_env)

  ## The pooled within-environment covariance is both the shrinkage target for
  ## the Sigma_e and the reference against which the model parameterises the
  ## conditional variance, so it is computed first.
  ## `sigma_known` may be a single matrix, used for every environment, or a
  ## list of `n_env` matrices, one per environment.  The list form is what
  ## `tests/test_sbc.R` needs to run its `known` arm under heterogeneous
  ## covariances: with a single matrix the generative model is forced to have a
  ## common Sigma_e, which is the one configuration in which the slope and
  ## covariance parameterisations of the model coincide, and so the one
  ## configuration in which the comparison cannot say anything.
  sigma_known_list <- NULL
  if (!is.null(sigma_known)) {
    sigma_known_list <- if (is.list(sigma_known)) sigma_known else
      rep(list(sigma_known), n_env)
    if (length(sigma_known_list) != n_env) {
      stop(sprintf(
        "sigma_known has %d matrices for %d environments.",
        length(sigma_known_list), n_env), call. = FALSE)
    }
  }

  sigma_bar <- if (!is.null(sigma_known_list)) {
    ## The pooled reference is the average of the known covariances, which is
    ## what `pooled_covariance()` estimates when the environments are balanced.
    m <- Reduce(`+`, sigma_known_list) / n_env
    enforce_pd(m / tcrossprod(scale_x), shrinkage = 0)
  } else {
    pooled_covariance(xs, z)
  }
  l_sigma_bar <- lower_chol(sigma_bar, "pooled within-environment Sigma")

  for (j in seq_along(envs)) {
    idx <- which(z == envs[j])
    if (length(idx) <= p) {
      warning(sprintf(
        "Environment %d has %d observations for p = %d; the plug-in covariance is rank deficient and is being regularised.",
        envs[j], length(idx), p), call. = FALSE)
    }
    xe <- xs[idx, , drop = FALSE]
    ye <- y[idx]

    sigma_e <- if (!is.null(sigma_known_list)) {
      ## Known covariance: bypass estimation entirely.  Used by
      ## tests/test_sbc.R to separate implementation error from the effect of
      ## plugging in an estimate.  On the standardised scale the supplied
      ## matrix must be transformed the same way the data were.
      m <- sigma_known_list[[j]] / tcrossprod(scale_x)
      attr(m, "shrinkage") <- 0
      enforce_pd((m + t(m)) / 2)
    } else {
      estimate_covariance(xe, method = cov_method,
                          target = sigma_bar,
                          lambda = cov_lambda)
    }
    cond_e[j] <- attr(sigma_e, "condition")
    shrink_e[j] <- attr(sigma_e, "shrinkage")
    l_sigma[j, , ] <- lower_chol(sigma_e,
                                 sprintf("Sigma for environment %d", envs[j]))

    n_e[j] <- length(idx)
    xbar[j, ] <- colMeans(xe)
    cp <- crossprod(xe)
    xx[j, , ] <- (cp + t(cp)) / 2      # exact symmetry for quad_form_sym()
    xy[j, ] <- as.vector(crossprod(xe, ye))
    ysum[j] <- sum(ye)
    yy[j] <- sum(ye^2)

    ## Scatter matrix, as a lower Cholesky factor C with S = C C'.  Only the
    ## full-covariance model uses it, but it is cheap to always supply.
    chol_scatter[j, , ] <- scatter_chol(xe)
  }

  ## ---- Target-domain covariance ---------------------------------------
  ## Sigma_0 must not be shrunk towards the training covariance: the whole
  ## point is that the target domain has shifted, and pooling it back towards
  ## the training environments would erase exactly the difference the
  ## predictive variance is supposed to reflect.  Ledoit-Wolf towards a scaled
  ## identity is used instead whenever the caller asked for pooling.
  cov_method0 <- if (identical(cov_method, "pooled")) "ledoit_wolf" else
    cov_method
  if (is.null(x0s)) {
    ## With no target sample the natural reference is the pooled training
    ## covariance; the predictive block is then switched off (N0 = 0).
    sigma0 <- sigma_bar
    mu0 <- colMeans(xs)
    x0s <- matrix(0, 0, p)
  } else if (!is.null(sigma_known_list)) {
    ## No known target covariance is supplied separately, so the pooled
    ## reference stands in for it, exactly as in the no-target branch above.
    sigma0 <- sigma_bar
    mu0 <- colMeans(x0s)
  } else {
    sigma0 <- estimate_covariance(x0s, method = cov_method0,
                                  lambda = cov_lambda)
    mu0 <- colMeans(x0s)
  }
  cond0 <- attr(sigma0, "condition")
  l_sigma0 <- lower_chol(sigma0, "Sigma for the target domain")

  stan_data <- list(
    P = p,
    E = n_env,
    n_e = n_e,
    N = n,
    xbar = xbar,
    xx = xx,
    xy = xy,
    ## rstan passes a length-one R vector as a Stan scalar, so a model
    ## declaring vector[P] or vector[E] fails to initialise whenever P or E is
    ## one.  as.array() keeps the dimension and is a no-op otherwise.  Without
    ## it the package cannot fit a single covariate at all, which is exactly
    ## the case the one-dimensional illustration needs.
    ysum = as.array(ysum),
    yy = as.array(yy),
    L_Sigma = l_sigma,
    L_Sigma_bar = l_sigma_bar,
    L_Sigma0 = l_sigma0,
    hmu = as.array(if (is.null(hmu_known)) colMeans(xs) else
      (hmu_known - centre) / scale_x),  # pooled mean of X, as in Section 2
    eta_lkj = eta_lkj,
    sd_mu_scale = sd_mu_scale,
    a_tau = a_tau,
    b_tau = b_tau,
    ncp = ncp,
    v_prior_shape = v_prior_shape,
    v_prior_rate = v_prior_rate,
    N0 = nrow(x0s),
    X0 = x0s,
    mu0 = as.array(mu0),
    compute_log_lik = as.integer(compute_log_lik),
    X = if (compute_log_lik) xs else matrix(0, 0, p),
    Y = if (compute_log_lik) y else numeric(0),
    Z = if (compute_log_lik) z else integer(0),

    ## Extra inputs consumed only by stan/gi_hd_fullcov.stan.  rstan ignores
    ## list elements that the model does not declare, so both models can be
    ## fed the same list.
    chol_scatter = chol_scatter,
    xbar0 = mu0,
    chol_scatter0 = if (nrow(x0s) > 1L) {
      scatter_chol(x0s)
    } else {
      diag(sqrt(diag(sigma0)), p)
    },
    n0 = nrow(x0s),
    sd_sigma_scale = sd_sigma_scale,
    sd_b_scale = sd_b_scale
  )

  list(
    stan_data = stan_data,
    scaling = list(centre = centre, scale = scale_x, standardize = standardize),
    diagnostics = list(
      condition_train = cond_e,
      condition_target = cond0,
      condition_pooled = attr(sigma_bar, "condition"),
      shrinkage_train = shrink_e,
      shrinkage_target = attr(sigma0, "shrinkage"),
      n_per_env = n_e
    )
  )
}

#' Fit the Bayesian Generative Invariance model.
#'
#' @inheritParams prepare_bgi_data
#' @param model A compiled `stanmodel`; obtained from `load_bgi_model()`.
#'   Either `gi_hd.stan` (covariances plugged in) or `gi_hd_fullcov.stan`
#'   (covariances inferred).  The same data list serves both.
#' @param chains,iter,warmup,seed,cores Passed to `rstan::sampling()`.
#' @param adapt_delta,max_treedepth NUTS control parameters.
#' @param ... Further arguments to `rstan::sampling()`.
#' @return An object of class `bgi_fit`: a list with posterior draws on the
#'   original covariate scale (`draws`), the raw `stanfit`, sampler
#'   diagnostics and the covariance diagnostics from `prepare_bgi_data()`.
fit_bgi <- function(x, y, z, x0 = NULL,
                    model,
                    cov_method = "pooled",
                    cov_lambda = 0.01,
                    standardize = TRUE,
                    eta_lkj = 2,
                    sd_mu_scale = 2.5,
                    sd_sigma_scale = 2.5,
                    sd_b_scale = 1,
                    a_tau = 0.5,
                    b_tau = 0.5,
                    ncp = 1,
                    v_prior_shape = 3,
                    v_prior_rate = 2,
                    sigma_known = NULL,
                    hmu_known = NULL,
                    compute_log_lik = FALSE,
                    chains = 4,
                    iter = 2000,
                    warmup = floor(iter / 2),
                    seed = 1,
                    cores = 1,
                    adapt_delta = 0.9,
                    max_treedepth = 12,
                    ...) {
  prep <- prepare_bgi_data(x, y, z, x0,
                           cov_method = cov_method,
                           cov_lambda = cov_lambda,
                           standardize = standardize,
                           eta_lkj = eta_lkj,
                           sd_mu_scale = sd_mu_scale,
                           sd_sigma_scale = sd_sigma_scale,
                           sd_b_scale = sd_b_scale,
                           a_tau = a_tau,
                           ncp = ncp,
                           v_prior_shape = v_prior_shape,
                           v_prior_rate = v_prior_rate,
                           sigma_known = sigma_known,
                           hmu_known = hmu_known,
                           b_tau = b_tau,
                           compute_log_lik = compute_log_lik)

  fit <- rstan::sampling(
    model,
    data = prep$stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    cores = cores,
    refresh = 0,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    ...
  )

  draws <- extract_bgi_draws(fit, prep$scaling)

  structure(
    list(
      draws = draws,
      stanfit = fit,
      scaling = prep$scaling,
      diagnostics = c(prep$diagnostics, sampler_diagnostics(fit)),
      stan_data = prep$stan_data
    ),
    class = "bgi_fit"
  )
}

#' Map posterior draws back to the original covariate scale.
#'
#' @param fit A `stanfit` from `gi_hd.stan` or `gi_hd_reference.stan`.
#' @param scaling The `scaling` element returned by `prepare_bgi_data()`.
#' @return A list of posterior draw matrices on the original scale.
extract_bgi_draws <- function(fit, scaling) {
  pars <- c("alpha", "gamma", "K", "sigma_y", "sigma_cond", "S0", "tau")
  available <- intersect(pars, fit@model_pars)
  post <- rstan::extract(fit, pars = available, permuted = TRUE)

  d <- scaling$scale
  centre <- scaling$centre

  ## gamma~ = d * gamma  =>  gamma = gamma~ / d, element-wise per draw.
  gamma <- sweep(post$gamma, 2, d, "/")
  ## K~ = K / d  =>  K = d * K~.
  k <- sweep(post$K, 2, d, "*")
  ## alpha~ = alpha + gamma' centre  =>  alpha = alpha~ - gamma' centre.
  alpha <- as.vector(post$alpha) - as.vector(gamma %*% centre)

  out <- list(
    alpha = alpha,
    gamma = gamma,
    K = k,
    sigma_y = as.vector(post$sigma_y),
    sigma_cond = as.vector(post$sigma_cond),
    S0 = as.vector(post$S0),
    tau = as.vector(post$tau)
  )

  ## The predictive draws are already on the response scale, which the
  ## covariate standardisation leaves untouched.
  if ("Y0_pred" %in% fit@model_pars) {
    y0 <- rstan::extract(fit, pars = "Y0_pred", permuted = TRUE)$Y0_pred
    if (length(y0) > 0) {
      out$y0_pred <- y0
    }
  }
  if ("f0" %in% fit@model_pars) {
    f0 <- rstan::extract(fit, pars = "f0", permuted = TRUE)$f0
    if (length(f0) > 0) {
      out$f0 <- f0
    }
  }
  out
}

#' Collect the sampler diagnostics that should always be reported.
#'
#' @param fit A `stanfit`.
sampler_diagnostics <- function(fit) {
  summ <- rstan::summary(fit)$summary
  keep <- grep("^(alpha|gamma|K|sigma_y|sigma_cond|S0|tau)(\\[|$)",
               rownames(summ))
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

  list(
    ## Recorded because Rhat is undefined for a single chain, and downstream
    ## has to tell "not assessable" apart from "did not mix".
    n_chains = length(sp),
    ## Post-warmup draws, so that a divergence *rate* can be formed downstream.
    n_post_draws = sum(vapply(sp, nrow, numeric(1))),
    n_divergent = sum(vapply(sp, function(s) sum(s[, "divergent__"]),
                             numeric(1))),
    max_treedepth_hit = sum(vapply(
      sp, function(s) sum(s[, "treedepth__"] >= max(s[, "treedepth__"])),
      numeric(1))),
    max_rhat = suppressWarnings(max(summ[keep, "Rhat"], na.rm = TRUE)),
    min_ess_bulk = suppressWarnings(min(summ[keep, "n_eff"], na.rm = TRUE)),
    runtime_sec = sum(rstan::get_elapsed_time(fit))
  )
}

#' Cut posterior for the plug-in nuisances, by multiple imputation.
#'
#' The model conditions on `Sigma_e` and on the prior mean `hmu` as if they
#' were known, having estimated both from the same data.  Simulation-based
#' calibration (`tests/test_sbc.R`) shows what that costs: with the true values
#' supplied the posterior is calibrated (0.948 against nominal 0.95, ranks
#' uniform), and with the plug-ins conditioned on it is not (0.784, rank
#' histogram carrying 2.5 times the expected mass in its extreme bins).  The
#' posterior is not biased; it is too narrow, because it does not carry the
#' first-stage uncertainty.
#'
#' This targets the **cut distribution** (Plummer 2015; Jacob, Murray, Holmes
#' and Robert 2017) in place of the naive plug-in,
#'
#'   p_cut(theta | D)  =  int p(theta | D, Sigma, hmu) p(Sigma, hmu | X) d(...)
#'
#' estimated by multiple imputation: draw `n_imp` values of the nuisances from
#' their first-stage posterior, run the second stage once per draw with that
#' value held fixed, and pool the posterior draws.  The pooled sample is a
#' Monte Carlo approximation to the mixture, so no asymptotics are involved.
#'
#' Two properties make this the right correction here rather than a sandwich
#' adjustment.  It is exact in finite samples, which matters because the
#' variance of `Sigma^{-1}` at `n_e = 200` and `p` up to 12 is far from its
#' asymptotic value; and `Sigma` stays *fixed within each run*, so the
#' per-environment sufficient-statistic likelihood is untouched and the cost
#' per imputation is the cost of the existing fit.  The imputations are
#' independent and embarrassingly parallel.
#'
#' What it does not do: the cut posterior is not the full Bayes posterior.
#' Feedback from `Y` to `Sigma` is deliberately severed.  That is a modelling
#' choice, and it is the appropriate one when the first stage is credible on
#' its own — here `Sigma_e` is informed by `n_e` direct observations of `X`,
#' so the information lost by cutting is small.
#'
#' @param x,y,z,x0 As in `fit_bgi()`.
#' @param model Compiled `stanmodel`.
#' @param n_imp Number of imputations. 25 is usually ample; the Monte Carlo
#'   error of the pooled posterior falls as `1/sqrt(n_imp)`.
#' @param sigma_common Draw a single covariance shared by all environments
#'   (`TRUE`, appropriate when the environments genuinely share one) or one per
#'   environment (`FALSE`).
#' @param impute_hmu Also impute the prior mean `hmu`.
#' @param ... Passed to `fit_bgi()`.
#' @return A `bgi_fit`-like list whose `draws` pool all imputations, with an
#'   extra `n_imp` element and per-imputation diagnostics.
fit_bgi_cut <- function(x, y, z, x0 = NULL, model, n_imp = 25L,
                        sigma_common = TRUE, impute_hmu = TRUE, ...) {
  x <- as.matrix(x)
  z <- as.integer(as.factor(z))
  n <- nrow(x)
  p <- ncol(x)
  envs <- sort(unique(z))
  n_env <- length(envs)

  ## ---- First-stage posterior for the nuisances -------------------------
  ## Under the standard non-informative prior p(Sigma) proportional to
  ## |Sigma|^{-(p+1)/2}, the posterior given the covariate data is
  ## inverse-Wishart with the within-group scatter as scale.  Drawing from it
  ## is what propagates the first-stage uncertainty.
  scatter <- matrix(0, p, p)
  for (e in envs) {
    xe <- x[z == e, , drop = FALSE]
    xc <- sweep(xe, 2, colMeans(xe), "-")
    scatter <- scatter + crossprod(xc)
  }
  df_w <- n - n_env
  if (df_w <= p + 1L) {
    stop("Too few degrees of freedom (", df_w, ") to impute a ", p, "x", p,
         " covariance.", call. = FALSE)
  }
  scatter <- (scatter + t(scatter)) / 2
  scatter_inv <- solve(scatter)

  ## The prior mean for the mu_e is estimated from the E environment means, so
  ## its uncertainty is governed by the between-environment spread and by E,
  ## not by the total sample size.
  env_means <- do.call(rbind, lapply(envs, function(e) {
    colMeans(x[z == e, , drop = FALSE])
  }))
  hmu_hat <- colMeans(env_means)
  hmu_var <- if (n_env > 1L) stats::cov(env_means) / n_env else
    diag(0, p)

  draw_nuisance <- function() {
    ## W ~ Wishart(df, scatter^{-1}) implies W^{-1} ~ InvWishart(df, scatter)
    ## with E[W^{-1}] = scatter / (df - p - 1), which is the posterior mean of
    ## Sigma under the non-informative prior. No further scaling.
    sigma <- solve(stats::rWishart(1, df_w, scatter_inv)[, , 1])
    sigma <- (sigma + t(sigma)) / 2
    hmu <- if (impute_hmu && n_env > 1L) {
      as.vector(mvtnorm::rmvnorm(1, hmu_hat, hmu_var))
    } else {
      hmu_hat
    }
    list(sigma = sigma, hmu = hmu)
  }

  ## ---- Second stage, once per imputation -------------------------------
  pooled <- NULL
  diags <- list()
  n_ok <- 0L
  for (m in seq_len(n_imp)) {
    nu <- draw_nuisance()
    fit <- tryCatch(
      fit_bgi(x, y, z, x0, model = model,
              sigma_known = nu$sigma, hmu_known = nu$hmu, ...),
      error = function(e) {
        message("  imputation ", m, " failed: ", conditionMessage(e)); NULL
      })
    if (is.null(fit)) next
    n_ok <- n_ok + 1L
    diags[[length(diags) + 1L]] <- fit$diagnostics

    if (is.null(pooled)) {
      pooled <- fit$draws
    } else {
      for (nm in names(pooled)) {
        pooled[[nm]] <- if (is.matrix(pooled[[nm]])) {
          rbind(pooled[[nm]], fit$draws[[nm]])
        } else {
          c(pooled[[nm]], fit$draws[[nm]])
        }
      }
    }
  }
  if (n_ok == 0L) {
    stop("Every imputation failed.", call. = FALSE)
  }

  structure(
    list(
      draws = pooled,
      n_imp = n_ok,
      diagnostics = list(
        n_divergent = sum(vapply(diags, `[[`, numeric(1), "n_divergent")),
        max_rhat = max(vapply(diags, `[[`, numeric(1), "max_rhat")),
        min_ess_bulk = min(vapply(diags, `[[`, numeric(1), "min_ess_bulk")),
        runtime_sec = sum(vapply(diags, `[[`, numeric(1), "runtime_sec"))
      )
    ),
    class = c("bgi_cut_fit", "bgi_fit")
  )
}
