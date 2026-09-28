/* -----------------------------------------------------------------------
 * Reference implementation of the BGI model.
 *
 * Same posterior as `gi_hd.stan`, written observation by observation with no
 * sufficient-statistic algebra and no precomputation.  It is deliberately
 * slow and deliberately transparent: its only purpose is to validate the fast
 * model.  `tests/test_fast_vs_reference.R` fits both to the same data and
 * checks that the posteriors agree.
 *
 * The two log-posteriors differ by an additive constant that depends on the
 * data but not on the parameters (the fast model drops the term
 * -0.5 * tr(Sigma_e^{-1} S_e) from the covariate likelihood), so `lp__` is
 * shifted while the posterior distribution is identical.
 * --------------------------------------------------------------------- */

functions {
  vector chol_solve(matrix L, vector v) {
    return mdivide_right_tri_low(mdivide_left_tri_low(L, v)', L)';
  }
  real quad_form_inv_chol(vector v, matrix L) {
    return dot_self(mdivide_left_tri_low(L, v));
  }
  real max_variance_excess(vector K, array[] matrix L, matrix L_bar) {
    real base = quad_form_inv_chol(K, L_bar);
    real out = 0;
    for (m in 1:size(L)) {
      out = fmax(out, quad_form_inv_chol(K, L[m]) - base);
    }
    return out;
  }
}

data {
  int<lower=1> P;
  int<lower=1> E;
  int<lower=1> N;
  matrix[N, P] X;
  vector[N] Y;
  array[N] int<lower=1, upper=E> Z;

  array[E] matrix[P, P] L_Sigma;
  matrix[P, P] L_Sigma_bar;
  matrix[P, P] L_Sigma0;

  vector[P] hmu;
  real<lower=0> eta_lkj;
  real<lower=0> sd_mu_scale;
  real<lower=0> a_tau;
  real<lower=0> b_tau;
  real<lower=0, upper=1> ncp;           // degree of non-centring; see gi_hd.stan

  // The prior on the free part of the conditional variance, mirroring
  // gi_hd.stan.  This model exists to validate that model's likelihood, which
  // it can only do if the two differ by a constant; carrying a different prior
  // makes the difference vary with v_raw instead, and the identity test then
  // fails for a reason that has nothing to do with the likelihood.  That is
  // what happened when the package default changed from Jeffreys to a proper
  // inverse-gamma and this file was left behind.
  real<lower=0> v_prior_shape;          // 0 selects the Jeffreys prior
  real<lower=0> v_prior_rate;

  int<lower=0> N0;
  matrix[N0, P] X0;
  vector[P] mu0;
}

transformed data {
  array[E + 1] matrix[P, P] L_all;
  for (e in 1:E) {
    L_all[e] = L_Sigma[e];
  }
  L_all[E + 1] = L_Sigma0;
}

parameters {
  real alpha_raw;
  vector[P] gamma_raw;
  vector[P] K_raw;
  real<lower=0> v_raw;
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
  vector[P] K = K_raw * ncp_mult;
  real base = quad_form_inv_chol(K, L_Sigma_bar);
  real<lower=0> v = v_raw + max_variance_excess(K, L_all, L_Sigma_bar);
  real<lower=0> sigma_cond = sqrt(v);
  real<lower=0> sigma_y_sq = v + base;
  real<lower=0> sigma_y = sqrt(sigma_y_sq);
  real<lower=0> S0 = sqrt(sigma_y_sq - quad_form_inv_chol(K, L_Sigma0));
}

model {
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
    K_raw ~ normal(0, raw_sd);
  }
  L_R ~ lkj_corr_cholesky(eta_lkj);
  sd_mu ~ normal(0, sd_mu_scale);
  mu ~ multi_normal_cholesky(hmu, diag_pre_multiply(sd_mu, L_R));

  for (i in 1:N) {
    int e = Z[i];
    vector[P] b_e = chol_solve(L_Sigma[e], K);
    real v_e = v - (dot_product(K, b_e) - base);
    X[i]' ~ multi_normal_cholesky(mu[e], L_Sigma[e]);
    Y[i] ~ normal(alpha + dot_product(gamma, X[i])
                  + dot_product(b_e, X[i]' - mu[e]), sqrt(v_e));
  }
}

generated quantities {
  vector[N0] f0;
  vector[N0] Y0_pred;
  for (i in 1:N0) {
    f0[i] = alpha + dot_product(gamma, X0[i])
            + dot_product(K, chol_solve(L_Sigma0, X0[i]' - mu0));
    Y0_pred[i] = normal_rng(f0[i], S0);
  }
}
