/* -----------------------------------------------------------------------
 * BGI with environment-specific confounding corrections, in slope coordinates.
 *
 * ---- What changes and why -------------------------------------------------
 *
 * `gi_hd.stan` follows the companion paper exactly: inner-product invariance
 * makes K = Cov(eps_Y, X) common across environments, so the slope on the
 * centred covariates is Sigma_e^{-1} K and the conditional mean is
 *
 *     alpha + gamma' X + K' Sigma_e^{-1} (X - mu_e).
 *
 * That forces Sigma_e^{-1} into the mean function.  Since Sigma_e is a plug-in
 * estimate, Sigma_hat_e^{-1} is a *generated regressor*: noisy, and with its
 * noise not propagated.  Simulation-based calibration quantifies the damage --
 * 0.948 coverage when the true Sigma is supplied, 0.784 when it is estimated,
 * with 2.5 times the expected mass in the extreme rank bins.
 *
 * This model writes the conditional mean in terms of the slope
 *
 *     b_e := Sigma_e^{-1} K_e,        so that     K_e = Sigma_e b_e,
 *
 * giving
 *
 *     alpha + gamma' X + b_e' (X - mu_e)
 *       =  [alpha - b_e' mu_e]  +  (gamma + b_e)' X,
 *
 * which contains **no covariance at all**.  `b_e` is a regression slope and is
 * the quantity the likelihood actually identifies; K_e is derived.
 *
 * ---- The assumption this trades ------------------------------------------
 *
 * This is not a free reparameterisation and should not be presented as one.
 * A single common `b` would assume the Sigma_e are equal, which in general
 * they are not.  Instead the b_e are given a hierarchical prior,
 *
 *     b_e ~ N(b_bar, diag(sd_b)^2),
 *
 * so every environment carries its own confounding correction, shrunk towards
 * a common centre.  In exchange, inner-product invariance is relaxed: instead
 * of K being exactly common, the *slope* of eps_Y on the centred covariates is
 * common in distribution.  The two coincide when the Sigma_e coincide and
 * differ otherwise.  Both say "the confounding mechanism is stable across
 * environments"; this one says it in regression coordinates.
 *
 * This is the relaxation the manuscript's Discussion already proposes -- "each
 * environment could be assigned a distinct K^z_*, all shrunk through a prior".
 * It does introduce a new assumption, though, and an earlier version of this
 * comment denied it: once the K_e differ, which K belongs to an unlabelled
 * target is no longer settled by the model, and predicting with K_bar asserts
 * K_0 = K_bar at every draw.  Nothing in the hierarchical prior implies that,
 * and the target responses that would test it are never observed.  Section 2.2
 * of the manuscript states it as an assumption and says the transport results
 * proved for an exactly common K_* do not carry over.
 *
 * Identification of gamma does not depend on the relaxation.  Within
 * environment e the data identify the slope s_e = gamma + b_e and the
 * intercept a_e = alpha - b_e' mu_e; eliminating b_e gives
 *
 *     a_e + s_e' mu_e  =  alpha + gamma' mu_e,
 *
 * one equation per environment, solvable for (alpha, gamma) exactly when the
 * augmented environment means span R^{p+1} -- the same condition as before
 * (Assumption 1 of the manuscript), unchanged.
 *
 * ---- Where Sigma^{-1} survives, and why it must -----------------------------
 *
 * In the *likelihood* Sigma now enters only through quadratic forms
 * `b' Sigma b`, never through its inverse.  That closes the
 * generated-regressor channel and makes the positivity constraint on the
 * conditional variances far better conditioned.
 *
 * It does **not** disappear from the target-domain prediction, and an earlier
 * version of this model that removed it there was wrong.  Writing the target
 * correction as `b_bar' (X_0 - mu_0)` transfers the training *slope* to the
 * target, which is only correct if Sigma_0 equals the training covariances --
 * precisely what fails under covariate shift.  Measured, that error dropped
 * predictive coverage to 0.607.
 *
 * The transferable object is the covariance K, not the slope.  So the model
 * derives K_e = Sigma_e b_e per environment, averages to K_bar, and predicts
 * with `K_bar' Sigma_0^{-1} (X_0 - mu_0)`.  The Sigma^{-1} in the manuscript's
 * Eq. (2) is not an arbitrary modelling choice: it is what makes the correction
 * transfer across domains.
 *
 * It does NOT follow that the inverse stays out of the likelihood, and an
 * earlier version of this comment claimed that it did.  `max_quad` below is the
 * larger of the training quadratic forms and `b0_quad = K_bar' Sigma_0^{-1}
 * K_bar`, and `sigma_y_sq = v_raw + max_quad` is what every training
 * conditional variance `v_e` is built from.  So on the set where the target
 * term dominates, Sigma_0^{-1} enters the likelihood and the posterior for
 * gamma depends on the target covariance.  The fit is transductive.  That is
 * defensible where X_0 is observed, which is the setting of the paper, and
 * Section 2.2 of the manuscript says so; what is not defensible is calling it
 * confined to `generated quantities`.
 * --------------------------------------------------------------------- */

functions {
  /**
   * Quadratic form v' Sigma v from the lower Cholesky factor L of Sigma.
   */
  real quad_form_chol(vector v, matrix L) {
    return dot_self(L' * v);
  }

  /**
   * Mahalanobis norm v' Sigma^{-1} v.  Used only for the mu_e likelihood,
   * where Sigma_e is a fixed weight and no parameter is divided by it.
   */
  real quad_form_inv_chol(vector v, matrix L) {
    return dot_self(mdivide_left_tri_low(L, v));
  }
}

data {
  int<lower=1> P;                       // number of covariates
  int<lower=1> E;                       // number of training environments
  array[E] int<lower=1> n_e;            // observations per environment
  int<lower=1> N;                       // sum(n_e)

  // Per-environment sufficient statistics of the training data.
  array[E] vector[P] xbar;
  array[E] matrix[P, P] xx;             // sum_i X_ei X_ei'  (symmetric)
  array[E] vector[P] xy;
  vector[E] ysum;
  vector[E] yy;

  // Plugged-in covariances, as lower Cholesky factors.  Used only in
  // quadratic forms and in the mu_e likelihood -- never inverted in the mean.
  array[E] matrix[P, P] L_Sigma;
  matrix[P, P] L_Sigma_bar;
  matrix[P, P] L_Sigma0;

  // Hyperparameters.
  vector[P] hmu;
  real<lower=0> eta_lkj;
  real<lower=0> sd_mu_scale;
  real<lower=0> sd_b_scale;             // half-normal scale for the spread of
                                        //   the environment-specific slopes
  real<lower=0> a_tau;
  real<lower=0> b_tau;
  real<lower=0, upper=1> ncp;

  real<lower=0> v_prior_shape;          // 0 selects the Jeffreys prior
  real<lower=0> v_prior_rate;

  // Target domain.
  int<lower=0> N0;
  matrix[N0, P] X0;
  vector[P] mu0;

  int<lower=0, upper=1> compute_log_lik;
  matrix[compute_log_lik ? N : 0, P] X;
  vector[compute_log_lik ? N : 0] Y;
  array[compute_log_lik ? N : 0] int<lower=1, upper=E> Z;
}

parameters {
  real alpha_raw;
  vector[P] gamma_raw;
  vector[P] b_bar_raw;                  // population-level slope

  matrix[E, P] b_dev;                   // environment deviations, non-centred
  vector<lower=0>[P] sd_b;

  real<lower=0> v_raw;                  // free part of the error variance

  array[E] vector[P] mu;
  cholesky_factor_corr[P] L_R;
  vector<lower=0>[P] sd_mu;

  real<lower=0> tau2_num;
  real<lower=0> tau2_den;
}

transformed parameters {
  real<lower=0> tau2 = tau2_num / tau2_den;
  real<lower=0> tau = sqrt(tau2);
  real prior_scale = sqrt(v_raw) * tau;
  real ncp_mult = pow(prior_scale, ncp);

  real alpha = alpha_raw * ncp_mult;
  vector[P] gamma = gamma_raw * ncp_mult;
  vector[P] b_bar = b_bar_raw * ncp_mult;

  array[E] vector[P] b;                 // environment-specific slopes
  vector[P] K_bar;                      // transferable confounding covariance
  real b0_quad;                         // K_bar' Sigma_0^{-1} K_bar
  real max_quad;
  real<lower=0> sigma_y_sq;
  real<lower=0> sigma_y;
  real<lower=0> S0;                     // target-domain predictive sd
  real<lower=0> sigma_cond;

  {
    vector[P] k_acc = rep_vector(0, P);
    real acc = 0;
    for (e in 1:E) {
      b[e] = b_bar + sd_b .* to_vector(b_dev[e]);
      // K_e = Sigma_e b_e is the covariance Cov(eps_Y, X) in environment e.
      // It, not the slope, is the object inner-product invariance holds fixed
      // and therefore the object that transfers to a domain with a different
      // covariance.
      k_acc += multiply_lower_tri_self_transpose(L_Sigma[e]) * b[e];
      acc = fmax(acc, quad_form_chol(b[e], L_Sigma[e]));
    }
    K_bar = k_acc / E;
    b0_quad = quad_form_inv_chol(K_bar, L_Sigma0);
    max_quad = fmax(acc, b0_quad);
  }

  // sigma_y^2 = Var(eps_Y).  In environment e the conditional variance is
  // sigma_y^2 - K_e' Sigma_e^{-1} K_e = sigma_y^2 - b_e' Sigma_e b_e, and in
  // the target it is sigma_y^2 - K_bar' Sigma_0^{-1} K_bar.  Anchoring at the
  // largest keeps every one positive.
  sigma_y_sq = v_raw + max_quad;
  sigma_y = sqrt(sigma_y_sq);
  S0 = sqrt(sigma_y_sq - b0_quad);
  sigma_cond = sqrt(v_raw);
}

model {
  // ---- Priors -----------------------------------------------------------
  if (v_prior_shape > 0) {
    v_raw ~ inv_gamma(v_prior_shape, v_prior_rate);
  } else {
    target += -log(v_raw);
  }

  tau2_num ~ gamma(a_tau, 1);
  tau2_den ~ gamma(b_tau, 1);

  {
    real raw_sd = pow(prior_scale, 1 - ncp);
    alpha_raw ~ normal(0, raw_sd);
    gamma_raw ~ normal(0, raw_sd);
    b_bar_raw ~ normal(0, raw_sd);
  }

  to_vector(b_dev) ~ std_normal();
  sd_b ~ normal(0, sd_b_scale);

  L_R ~ lkj_corr_cholesky(eta_lkj);
  sd_mu ~ normal(0, sd_mu_scale);
  mu ~ multi_normal_cholesky(hmu, diag_pre_multiply(sd_mu, L_R));

  // ---- Likelihood -------------------------------------------------------
  for (e in 1:E) {
    real ne = n_e[e];
    vector[P] g_e = gamma + b[e];               // slope on the raw covariates
    real a_e = alpha - dot_product(b[e], mu[e]);
    real v_e = sigma_y_sq - quad_form_chol(b[e], L_Sigma[e]);
    real ss;

    // X_ei ~ N(mu_e, Sigma_e), up to a constant free of parameters.
    target += -0.5 * ne * quad_form_inv_chol(xbar[e] - mu[e], L_Sigma[e]);

    ss = yy[e]
         - 2 * a_e * ysum[e]
         - 2 * dot_product(g_e, xy[e])
         + ne * square(a_e)
         + 2 * ne * a_e * dot_product(g_e, xbar[e])
         + quad_form_sym(xx[e], g_e);
    target += -0.5 * ne * (log(2 * pi()) + log(v_e)) - 0.5 * ss / v_e;
  }
}

generated quantities {
  vector[N0] f0;
  vector[N0] Y0_pred;
  // The target correction is K_bar' Sigma_0^{-1} (X_0 - mu_0).  Sigma_0^{-1}
  // is unavoidable here and should not be removed: it is precisely what makes
  // the correction transfer to a domain whose covariance differs from the
  // training ones.  This is not its only appearance, though -- see the header:
  // `b0_quad` carries it into `sigma_y_sq`, and so into the likelihood,
  // whenever the target term is the larger one in `max_quad`.
  vector[P] K = K_bar;
  matrix[P, P] Sigma_mu = multiply_lower_tri_self_transpose(
                            diag_pre_multiply(sd_mu, L_R));
  vector[compute_log_lik ? N : 0] log_lik;

  for (i in 1:N0) {
    f0[i] = alpha + dot_product(gamma, X0[i])
            + dot_product(K_bar,
                          mdivide_right_tri_low(
                            mdivide_left_tri_low(L_Sigma0, X0[i]' - mu0)',
                            L_Sigma0)');
    Y0_pred[i] = normal_rng(f0[i], S0);
  }

  if (compute_log_lik) {
    array[E] real a_all;
    array[E] real sd_all;
    for (e in 1:E) {
      a_all[e] = alpha - dot_product(b[e], mu[e]);
      sd_all[e] = sqrt(sigma_y_sq - quad_form_chol(b[e], L_Sigma[e]));
    }
    for (i in 1:N) {
      int e = Z[i];
      log_lik[i] =
        normal_lpdf(Y[i] | a_all[e] + dot_product(gamma + b[e], X[i]),
                    sd_all[e])
        + multi_normal_cholesky_lpdf(X[i]' | mu[e], L_Sigma[e]);
    }
  }
}
