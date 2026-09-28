#!/usr/bin/env Rscript
## test_recovery.R --------------------------------------------------------
##
## End-to-end checks that the pipeline recovers what it claims to recover.
## Each check is cheap enough to run on a login node and is meant to be run
## after any change to the model or to the wrappers.
##
##   1. Standardisation round trip.  The map used by extract_bgi_draws() must
##      leave the fitted conditional mean function unchanged.  Note that the
##      two *posteriors* are not identical: standardising rescales the ridge
##      prior, which is the reason for standardising in the first place.  What
##      has to hold exactly is the algebra of the back-transformation, and
##      what has to hold approximately is that both fits recover the truth.
##   2. No confounding.  With confounding = 0 the truth is K = 0, and both BGI
##      and OLS should recover gamma.  If BGI fails here, the problem is in
##      the model rather than in the confounding correction.
##   3. Strong confounding.  BGI should recover gamma where OLS cannot, and
##      the true parents should be selected while the true nulls are not.
##   4. Predictive variance.  The target-domain predictive sd S_0 must differ
##      from the training conditional sd, and posterior predictive intervals
##      should cover at close to the nominal rate.  This is the property the
##      submitted implementation did not have: it sampled the predictive
##      draws with the training residual sd.
##
## Usage:  Rscript tests/test_recovery.R

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))

failures <- character(0)
check <- function(name, condition, detail = "") {
  status <- if (isTRUE(condition)) "PASS" else "FAIL"
  cat(sprintf("  [%s] %s%s\n", status, name,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!isTRUE(condition)) failures <<- c(failures, name)
  invisible(condition)
}

## ---- 1. Standardisation round trip -------------------------------------

cat("\n1. Standardisation round trip\n")

## (a) Exact algebra.  With X~ = (X - c) / D the transformation is
##     gamma~ = D gamma,  K~ = D^{-1} K,  alpha~ = alpha + gamma' c,
##     Sigma~ = D^{-1} Sigma D^{-1},  mu~ = D^{-1} (mu - c),
## and the conditional mean alpha + gamma' x + K' Sigma_0^{-1} (x - mu_0) must
## be invariant.  Checked on arbitrary values, so it tests the map itself
## rather than one fit.
set.seed(2024)
p_t <- 5
d_scale <- runif(p_t, 0.4, 3)
c_shift <- rnorm(p_t)
alpha_t <- 0.7
gamma_t <- rnorm(p_t)
k_t <- rnorm(p_t)
sigma0_t <- crossprod(matrix(rnorm(p_t * p_t), p_t)) + diag(p_t)
mu0_t <- rnorm(p_t)
x_t <- rnorm(p_t)

d_mat <- diag(d_scale)
d_inv <- diag(1 / d_scale)
mean_original <- alpha_t + sum(gamma_t * x_t) +
  as.numeric(crossprod(k_t, solve(sigma0_t, x_t - mu0_t)))
mean_standardised <- (alpha_t + sum(gamma_t * c_shift)) +
  sum((d_mat %*% gamma_t) * ((x_t - c_shift) / d_scale)) +
  as.numeric(crossprod(d_inv %*% k_t,
                       solve(d_inv %*% sigma0_t %*% d_inv,
                             (x_t - c_shift) / d_scale - (mu0_t - c_shift) / d_scale)))
check("conditional mean is invariant under standardisation",
      abs(mean_original - mean_standardised) < 1e-8,
      sprintf("difference = %.2e", abs(mean_original - mean_standardised)))

## (b) Both fits recover the truth.  The posteriors themselves differ, because
## standardising rescales the ridge prior; only the recovered truth is
## comparable.
bgi_set_seed(101, 1)
d1 <- simulate_gi_data(n_e = 200, p = 4, s0 = 2, n_env = 5, q = 3,
                       confounding = 1, n0 = 200)
f_std <- fit_bgi(d1$x, d1$y, d1$z, d1$x0, model = model, standardize = TRUE,
                 chains = 2, iter = 2000, seed = 3)
f_raw <- fit_bgi(d1$x, d1$y, d1$z, d1$x0, model = model, standardize = FALSE,
                 chains = 2, iter = 2000, seed = 3)
err_std <- sqrt(mean((colMeans(f_std$draws$gamma) - d1$truth$gamma)^2))
err_raw <- sqrt(mean((colMeans(f_raw$draws$gamma) - d1$truth$gamma)^2))
err_ols <- sqrt(mean((fit_ols(d1$x, d1$y)$gamma - d1$truth$gamma)^2))
## E = p + 1 here, the minimum for identifiability, so the between-environment
## information about K is thin and neither parameterisation is sharp.  What
## the check establishes is that neither is *broken*: both beat OLS, and they
## agree with each other to within a factor of two.
check("both parameterisations behave",
      max(err_std, err_raw) < err_ols &&
        max(err_std, err_raw) / min(err_std, err_raw) < 2,
      sprintf("rmse standardised %.4f, raw %.4f, OLS %.4f",
              err_std, err_raw, err_ols))

## ---- 2. No confounding -------------------------------------------------

cat("\n2. No hidden confounding (K = 0)\n")
bgi_set_seed(102, 1)
d2 <- simulate_gi_data(n_e = 400, p = 4, s0 = 2, n_env = 5, q = 3,
                       confounding = 0, n0 = 300)
f2 <- fit_bgi(d2$x, d2$y, d2$z, d2$x0, model = model, chains = 4, iter = 2000,
              seed = 4)
err2 <- sqrt(mean((colMeans(f2$draws$gamma) - d2$truth$gamma)^2))
ols2 <- fit_ols(d2$x, d2$y, d2$x0)
err2_ols <- sqrt(mean((ols2$gamma - d2$truth$gamma)^2))
check("BGI recovers gamma", err2 < 0.1,
      sprintf("rmse = %.4f (OLS %.4f, which is the correct model here)",
              err2, err2_ols))
## K* = 0 exactly in this scenario.  The check is whether zero is credible,
## not whether the posterior mean is numerically small: with E = p + 1
## environments the between-environment information about K is thin, so a
## posterior mean of a few tenths is entirely consistent with K* = 0.
k_z <- abs(colMeans(f2$draws$K)) / apply(f2$draws$K, 2, stats::sd)
check("zero is credible for every entry of K", max(k_z) < 3,
      sprintf("max |posterior mean| / posterior sd = %.2f (max |K| = %.3f)",
              max(k_z), max(abs(colMeans(f2$draws$K)))))

## ---- 3. Strong confounding ---------------------------------------------

cat("\n3. Strong hidden confounding\n")
bgi_set_seed(103, 1)
d3 <- simulate_gi_data(n_e = 400, p = 4, s0 = 2, n_env = 9, q = 3,
                       confounding = 2, n0 = 300)
f3 <- fit_bgi(d3$x, d3$y, d3$z, d3$x0, model = model, chains = 4, iter = 2000,
              seed = 5)
ols3 <- fit_ols(d3$x, d3$y, d3$x0)
err3 <- sqrt(mean((colMeans(f3$draws$gamma) - d3$truth$gamma)^2))
err3_ols <- sqrt(mean((ols3$gamma - d3$truth$gamma)^2))
check("BGI beats OLS on gamma", err3 < err3_ols,
      sprintf("BGI rmse %.4f vs OLS %.4f", err3, err3_ols))

sel3 <- select_parents(f3$draws$gamma, alpha = 0.05, rule = "sign")
m3 <- support_metrics(sel3$selected, d3$truth$parents, 4)
sel3_ols <- which(ols3$p_value < 0.05)
m3_ols <- support_metrics(sel3_ols, d3$truth$parents, 4)
check("BGI selects all true parents", m3$tpr == 1,
      sprintf("selected {%s}, parents {%s}",
              paste(sel3$selected, collapse = ","),
              paste(d3$truth$parents, collapse = ",")))
check("BGI makes no false discovery", m3$fp == 0,
      sprintf("BGI fp = %d, OLS fp = %d", m3$fp, m3_ols$fp))

## ---- 4. Predictive variance --------------------------------------------

cat("\n4. Target-domain predictive variance\n")
s0_mean <- mean(f3$draws$S0)
cond_mean <- mean(f3$draws$sigma_cond)
sy_mean <- mean(f3$draws$sigma_y)
check("S0 differs from the training conditional sd",
      abs(s0_mean - cond_mean) / cond_mean > 0.02,
      sprintf("S0 = %.3f, sigma_cond = %.3f, sigma_y = %.3f (true %.3f)",
              s0_mean, cond_mean, sy_mean, d3$truth$sigma_y))
check("S0 lies between sigma_cond and sigma_y",
      s0_mean >= cond_mean - 1e-6 && s0_mean <= sy_mean + 1e-6,
      "the target domain is more dispersed, so less of eps_Y is explained")

pm3 <- predictive_metrics(f3$draws$y0_pred, d3$y0)
om3 <- interval_metrics(d3$y0, ols3$pred_lower, ols3$pred_upper)
check("predictive coverage is near nominal",
      abs(pm3$coverage - 0.95) < 0.08,
      sprintf("BGI %.3f (width %.2f) vs OLS %.3f (width %.2f)",
              pm3$coverage, pm3$mean_width, om3$coverage, om3$mean_width))

## ---- Report ------------------------------------------------------------

cat("\n")
if (length(failures) == 0L) {
  cat("All recovery checks passed.\n")
} else {
  cat("FAILED checks: ", paste(failures, collapse = "; "), "\n", sep = "")
  quit(status = 1)
}
