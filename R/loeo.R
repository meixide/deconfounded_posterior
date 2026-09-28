## loeo.R -----------------------------------------------------------------
##
## Leave-one-environment-out evaluation.
##
## This exists because of Referee 2's criticism of Section 3.2:
##
##   "the reported 'remarkable empirical coverage of 0.95' seems to be based on
##    a single test environment (Madrid). Coverage is computed across
##    individuals within Madrid, but all prediction intervals share the same
##    posterior [...] making the coverage events highly correlated across
##    individuals. More fundamentally, this measures within-environment
##    marginal coverage for one specific held-out domain."
##
## The criticism is correct on both counts, and they need different fixes.
##
## 1. Correlated coverage within an environment.  Every interval in a held-out
##    domain is built from one posterior draw set, so the individual coverage
##    indicators are not independent and their mean has a standard error far
##    larger than sqrt(p(1-p)/n_0) would suggest.  Quoting a single figure to
##    two decimals from one domain overstates the precision considerably.
##    The fix is to treat the *environment* as the unit of replication:
##    compute within-environment coverage, then summarise across environments
##    and take the Monte Carlo error from the between-environment spread.
##
## 2. One held-out domain cannot validate generalisation.  A well-calibrated
##    and a badly-calibrated method can both land on 0.95 in any particular
##    domain.  The fix is to hold out every environment in turn.
##
## `loeo_evaluate()` does both.  It also records, per fold, the diagnostics
## that the simulations showed matter: the residual degrees of freedom
## `E - (p + 1)` in the between-environment regression, the conditioning of the
## plugged-in covariances, and the sampler diagnostics.  A fold that failed its
## diagnostics must not be silently averaged in.

#' Leave-one-environment-out evaluation of the BGI predictive distribution.
#'
#' @param x Covariates, `N x p`.
#' @param y Response, length `N`.
#' @param z Environment labels, length `N`.
#' @param model Compiled `stanmodel` from `load_bgi_model()`.
#' @param folds Environments to hold out in turn.  Defaults to every
#'   environment with at least `min_n0` observations.
#' @param min_n0 Minimum size for an environment to be used as a target.
#' @param min_train_env Minimum number of training environments required; a
#'   fold with fewer is skipped rather than fitted, since identification needs
#'   the environment means to span R^p.
#' @param level Nominal level for the predictive intervals.
#' @param max_target Optionally subsample the held-out environment to this many
#'   rows, to bound the cost of the `generated quantities` block.  The fit
#'   itself does not depend on it.
#' @param baselines Also fit OLS, pooled GI and IV on each fold.
#' @param ... Passed to `fit_bgi()`.
#' @return A data frame with one row per fold.
loeo_evaluate <- function(x, y, z,
                          model,
                          folds = NULL,
                          min_n0 = 200L,
                          min_train_env = NULL,
                          level = 0.95,
                          max_target = 2000L,
                          baselines = TRUE,
                          ...) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  z <- as.character(z)
  p <- ncol(x)
  envs <- sort(unique(z))
  sizes <- table(z)

  if (is.null(min_train_env)) {
    ## Identification needs the environment means to span R^p, so E >= p + 1.
    ## The simulations show that the *minimum* is not enough for reliable
    ## inference, so this is a floor, not a recommendation.
    min_train_env <- p + 1L
  }
  if (is.null(folds)) {
    folds <- names(sizes)[sizes >= min_n0]
  }
  folds <- intersect(folds, envs)
  if (length(folds) == 0L) {
    stop("No environment has at least ", min_n0, " observations.", call. = FALSE)
  }

  rows <- list()
  for (e0 in folds) {
    is_target <- z == e0
    train_env <- setdiff(envs, e0)
    if (length(train_env) < min_train_env) {
      message("Skipping fold ", e0, ": only ", length(train_env),
              " training environments for p = ", p, ".")
      next
    }

    x_tr <- x[!is_target, , drop = FALSE]
    y_tr <- y[!is_target]
    z_tr <- z[!is_target]

    idx0 <- which(is_target)
    if (length(idx0) > max_target) {
      idx0 <- sort(sample(idx0, max_target))
    }
    x0 <- x[idx0, , drop = FALSE]
    y0 <- y[idx0]

    message(sprintf("[%s] fold %-28s  n_train %7d  E_train %3d  n_target %5d",
                    format(Sys.time(), "%H:%M:%S"), e0, nrow(x_tr),
                    length(train_env), length(idx0)))

    fit <- tryCatch(
      fit_bgi(x_tr, y_tr, z_tr, x0, model = model, ...),
      error = function(err) {
        message("  fit failed: ", conditionMessage(err)); NULL
      })
    if (is.null(fit)) next

    pm <- predictive_metrics(fit$draws$y0_pred, y0, level = level)

    ## Score the same fit at the predictive scales it did *not* use, so that
    ## the corrected S_0 can be compared with the training-residual scale the
    ## submitted code used (`normal_rng(cond_mean, sigmay)`).  Paired: one set
    ## of standard normal draws serves every scale.
    naive <- list(sigma_y = NA_real_, sigma_cond = NA_real_)
    naive_is <- list(sigma_y = NA_real_, sigma_cond = NA_real_)
    if (!is.null(fit$draws$f0)) {
      f0 <- fit$draws$f0
      zmat <- matrix(stats::rnorm(length(f0)), nrow(f0), ncol(f0))
      for (nm in c("sigma_y", "sigma_cond")) {
        m <- predictive_metrics_at_scale(f0, fit$draws[[nm]], y0,
                                         level = level, z = zmat)
        naive[[nm]] <- m$coverage
        naive_is[[nm]] <- m$interval_score
      }
    }

    gamma_draws <- fit$draws$gamma
    sel <- select_parents(gamma_draws, alpha = 1 - level, rule = "sign")

    row <- data.frame(
      target_env = e0,
      n_train = nrow(x_tr),
      n_target = length(idx0),
      n_train_env = length(train_env),
      p = p,
      ## The quantity the simulations identified as the one to monitor.
      between_df = length(train_env) - (p + 1L),
      coverage = pm$coverage,
      ## The same posterior scored at the two scales the model did not use.
      coverage_at_sigma_y = naive$sigma_y,
      coverage_at_sigma_cond = naive$sigma_cond,
      interval_score_at_sigma_y = naive_is$sigma_y,
      interval_score_at_sigma_cond = naive_is$sigma_cond,
      mean_width = pm$mean_width,
      interval_score = pm$interval_score,
      rmse = pm$rmse,
      n_selected = length(sel$selected),
      selected = paste(sel$selected, collapse = "|"),
      S0_mean = mean(fit$draws$S0),
      sigma_cond_mean = mean(fit$draws$sigma_cond),
      sigma_y_mean = mean(fit$draws$sigma_y),
      cov_condition_max = max(fit$diagnostics$condition_train),
      cov_condition_target = fit$diagnostics$condition_target,
      divergent = fit$diagnostics$n_divergent,
      post_draws = fit$diagnostics$n_post_draws,
      max_rhat = fit$diagnostics$max_rhat,
      min_ess = fit$diagnostics$min_ess_bulk,
      n_chains = fit$diagnostics$n_chains,
      runtime_sec = fit$diagnostics$runtime_sec,
      stringsAsFactors = FALSE
    )

    if (baselines) {
      ols <- tryCatch(fit_ols(x_tr, y_tr, x0), error = function(e) NULL)
      pgi <- tryCatch(fit_pooled_gi(x_tr, y_tr, z_tr, x0),
                      error = function(e) NULL)
      iv <- tryCatch(fit_iv_2sls(x_tr, y_tr, z_tr, x0),
                     error = function(e) NULL)
      for (nm in c("ols", "pooled_gi", "iv")) {
        f <- get(if (nm == "pooled_gi") "pgi" else nm)
        if (is.null(f) || is.null(f$pred_lower)) {
          row[[paste0(nm, "_coverage")]] <- NA_real_
          row[[paste0(nm, "_interval_score")]] <- NA_real_
          row[[paste0(nm, "_rmse")]] <- NA_real_
        } else {
          m <- interval_metrics(y0, f$pred_lower, f$pred_upper, level = level)
          row[[paste0(nm, "_coverage")]] <- m$coverage
          row[[paste0(nm, "_interval_score")]] <- m$interval_score
          row[[paste0(nm, "_rmse")]] <- sqrt(mean((f$pred_mean - y0)^2))
        }
      }
    }

    rows[[length(rows) + 1L]] <- row
    message(sprintf("  coverage %.3f  width %.2f  interval score %.2f  (%d divergences)",
                    pm$coverage, pm$mean_width, pm$interval_score,
                    fit$diagnostics$n_divergent))
  }

  if (length(rows) == 0L) {
    stop("Every fold failed.", call. = FALSE)
  }
  do.call(rbind, rows)
}

#' Summarise leave-one-environment-out folds, treating the environment as the
#' unit of replication.
#'
#' The standard error is taken from the spread *across* environments, which is
#' the only defensible choice: within an environment the coverage indicators
#' share a posterior and are strongly dependent, so a binomial standard error
#' computed over individuals is far too small.
#'
#' @param folds Output of `loeo_evaluate()`.
#' @param drop_bad Exclude folds whose sampler diagnostics failed.
#' @param max_divergent_frac Largest tolerated fraction of post-warmup
#'   transitions that diverged.
#'
#'   Deliberately loose, because the alternative is worse.  Over the 52 BRFSS
#'   folds, divergences are sporadic (0 to 8 per 4,000 draws) and carry no
#'   signal: the count correlates with nothing measured — target covariance
#'   conditioning `r = -0.25`, `Rhat` `r = 0.15`, coverage `r = 0.003` — and
#'   mean coverage is the same with them and without (0.950 against 0.954).
#'   The worst fold on this criterion still had `Rhat = 1.0027` and
#'   `ESS = 2055`.
#'
#'   Excluding those folds would therefore discard sound fits on a diagnostic
#'   unrelated to the estimate, and would risk selecting the summary on domain
#'   shift itself — precisely the dimension the study is about.  The threshold
#'   is set where a divergence rate really does indicate a posterior that
#'   cannot be trusted, and everything below it is reported rather than
#'   filtered.  Every fold's count is in `loeo_folds.csv` either way.
#' @param min_ess_floor Smallest tolerated bulk ESS.  Nothing in the BRFSS run
#'   comes near it (the minimum over 52 folds is ~1,100); it guards against a
#'   future run where the chains stick.
#' @param rhat_max Threshold for `max_rhat`.
#' @return A one-row-per-method data frame.
loeo_summary <- function(folds, drop_bad = TRUE, rhat_max = 1.01,
                         max_divergent_frac = 0.01,
                         min_ess_floor = 400) {
  n_all <- nrow(folds)
  if (drop_bad) {
    ## Rhat is undefined for a single chain, so an NA there means "not
    ## assessable", not "failed".  With two or more chains an NA means the
    ## chains did not mix and the fold must go.  Distinguishing the two matters
    ## because a smoke test run at one chain would otherwise drop every fold.
    rhat_ok <- if (is.null(folds$n_chains)) {
      is.na(folds$max_rhat) | folds$max_rhat <= rhat_max
    } else {
      ifelse(folds$n_chains <= 1L & is.na(folds$max_rhat),
             TRUE, !is.na(folds$max_rhat) & folds$max_rhat <= rhat_max)
    }
    ## Divergences are judged as a rate, not as a count: `post_draws` is what
    ## the fold actually sampled, so the threshold means the same thing however
    ## long the chains were run.
    post <- if (is.null(folds$post_draws)) rep(NA_real_, nrow(folds)) else
      folds$post_draws
    div_frac <- ifelse(is.na(post) | post <= 0, NA_real_,
                       folds$divergent / post)
    div_ok <- ifelse(is.na(div_frac), folds$divergent %in% 0,
                     div_frac <= max_divergent_frac)

    ess_ok <- if (is.null(folds$min_ess)) TRUE else
      is.na(folds$min_ess) | folds$min_ess >= min_ess_floor

    ## `%in% TRUE` rather than the bare logical: an NA in the subscript
    ## silently yields a row of NAs instead of dropping the row.
    ok <- (div_ok & rhat_ok & ess_ok) %in% TRUE
    if (any(!ok)) {
      message("Excluding ", sum(!ok), " of ", n_all,
              " folds that failed the sampler diagnostics.")
    }
    if (any(ok & folds$divergent > 0)) {
      ## Tolerated is not the same as absent, and the paper should say so.
      message("Retained ", sum(ok & folds$divergent > 0),
              " fold(s) with a divergence rate below ", max_divergent_frac,
              "; max retained rate ",
              signif(max(div_frac[ok], na.rm = TRUE), 3), ".")
    }
    folds <- folds[ok, , drop = FALSE]
  }
  if (nrow(folds) == 0L) {
    stop("No fold passed the sampler diagnostics (", n_all, " attempted). ",
         "Read the per-fold CSV before summarising: a run where every fold ",
         "failed is a modelling problem, not a summary to be averaged.",
         call. = FALSE)
  }

  methods <- c("bgi", "ols", "pooled_gi", "iv")
  get_col <- function(m, what) {
    nm <- if (m == "bgi") what else paste0(m, "_", what)
    if (nm %in% names(folds)) folds[[nm]] else rep(NA_real_, nrow(folds))
  }

  out <- do.call(rbind, lapply(methods, function(m) {
    cov <- get_col(m, "coverage")
    if (all(is.na(cov))) {
      return(NULL)
    }
    data.frame(
      method = m,
      n_folds = sum(!is.na(cov)),
      coverage = mean(cov, na.rm = TRUE),
      ## Across-environment standard error: environments are the replicates.
      coverage_se = stats::sd(cov, na.rm = TRUE) / sqrt(sum(!is.na(cov))),
      coverage_min = min(cov, na.rm = TRUE),
      coverage_max = max(cov, na.rm = TRUE),
      interval_score = mean(get_col(m, "interval_score"), na.rm = TRUE),
      rmse = mean(get_col(m, "rmse"), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out)) {
    stop("No method produced a coverage column; nothing to summarise.",
         call. = FALSE)
  }
  out[order(out$interval_score), ]
}

#' Stability of the selected covariate set across held-out environments.
#'
#' The manuscript reports one selected set from one split. Whether that set is
#' stable when the held-out domain changes is a much stronger claim, and it is
#' the one a reader will care about. Reports, for each covariate, the fraction
#' of folds in which it was selected.
#'
#' @param folds Output of `loeo_evaluate()`.
#' @param covariate_names Optional names for the columns of `x`.
loeo_selection_stability <- function(folds, covariate_names = NULL) {
  p <- folds$p[1]
  picks <- lapply(strsplit(folds$selected, "|", fixed = TRUE), as.integer)
  freq <- vapply(seq_len(p), function(j) {
    mean(vapply(picks, function(s) j %in% s, logical(1)))
  }, numeric(1))
  data.frame(
    covariate = if (is.null(covariate_names)) paste0("x", seq_len(p)) else
      covariate_names,
    selection_frequency = freq,
    stringsAsFactors = FALSE
  )[order(-freq), ]
}
