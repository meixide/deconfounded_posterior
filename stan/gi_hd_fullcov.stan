/* -----------------------------------------------------------------------
 * Bayesian Generative Invariance with the covariance matrices treated as
 * unknown parameters rather than plugged in.
 *
 * Same sampling model as `gi_hd.stan`:
 *
 *   X_ei             ~ N_p(mu_e, Sigma_e)
 *   Y_ei | X_ei      ~ N(alpha + gamma' X_ei + K' Sigma_e^{-1} (X_ei - mu_e),
 *                        sigma_y^2 - K' Sigma_e^{-1} K)
 *
 * but with Sigma_1, ..., Sigma_E and Sigma_0 given priors and inferred, so
 * that their estimation error propagates into the posterior for gamma.
 *
 * ---- Why this model exists ------------------------------------------------
 *
 * The plug-in version conditions on Sigma_hat as if it were known.  Since
 * gamma is the within-environment slope minus Sigma^{-1} K, any error in
 * Sigma_hat lands on gamma, and under strong confounding ||Sigma^{-1} K|| is
 * large, so that error is a first-order contribution to the uncertainty in
 * gamma -- not the asymptotically negligible one the plug-in argument
 * assumes at moderate n_e.  Empirically the plug-in model's credible
 * intervals for gamma cover well below their nominal level in that regime,
 * and the selection rules of Section 2.1 inherit the over-confidence as false
 * discoveries.
 *
 * ---- Priors ---------------------------------------------------------------
 *
 *   Sigma_e | nu, Sigma_bar ~ InvWishart(nu, (nu - p - 1) Sigma_bar)
 *
 * so that E[Sigma_e] = Sigma_bar: the environment covariances are shrunk
 * towards a common covariance whose value, and whose strength nu, are learned
 * from the data.  This is the model-based version of the ad hoc shrinkage
 * used by the plug-in variant.
 *
 * Sigma_bar and Sigma_0 each get an LKJ correlation prior and a half-normal
 * scale prior.  Sigma_0 is *not* tied to Sigma_bar: the target domain has
 * shifted, and pooling it back would erase the difference that the predictive
 * variance S_0^2 = sigma_y^2 - K' Sigma_0^{-1} K exists to express.  Sigma_0
 * is informed by the target covariate sample, which is observed.
 *
 * ---- Computation ----------------------------------------------------------
 *
 * The likelihood still uses only per-environment sufficient statistics.  The
 * covariate part needs the scatter matrices, supplied as Cholesky factors
 * C_e with S_e = C_e C_e', since
 *
 *   sum_i (x_i - mu)' Sigma^{-1} (x_i - mu)
 *       = sum of squares of L^{-1} C  +  n (xbar - mu)' Sigma^{-1} (xbar - mu).
 *
 * Cost is O(E p^3) per gradient evaluation, independent of N, as in
 * `gi_hd.stan`.
 * --------------------------------------------------------------------- */

functions {
  vector chol_solve(matrix L, vector v) {
    return mdivide_right_tri_low(mdivide_left_tri_low(L, v)', L)';
  }
  real quad_form_inv_chol(vector v, matrix L) {
    return dot_self(mdivide_left_tri_low(L, v));
  }
}

data {
  int<lower=1> P;
  int<lower=1> E;
  array[E] int<lower=1> n_e;
  int<lower=1> N;

  // Sufficient statistics of the training data.
  array[E] vector[P] xbar;
  array[E] matrix[P, P] chol_scatter;   // C_e with S_e = C_e C_e'
  array[E] matrix[P, P] xx;             // sum_i X_ei X_ei'
  array[E] vector[P] xy;
  vector[E] ysum;
  vector[E] yy;

  // Target domain.
  int<lower=0> N0;
  matrix[N0, P] X0;
  vector[P] xbar0;
  matrix[P, P] chol_scatter0;           // C_0 with S_0 = C_0 C_0'
  int<lower=0> n0;                      // target sample size (may exceed N0)

  // Hyperparameters.
  vector[P] hmu;
  real<lower=0> eta_lkj;
  real<lower=0> sd_mu_scale;
  real<lower=0> sd_sigma_scale;         // half-normal scale for sd(Sigma_bar),
                                        //   sd(Sigma_0)
  real<lower=0> a_tau;
  real<lower=0> b_tau;
  real<lower=0, upper=1> ncp;           // degree of non-centring; see gi_hd.stan
}

parameters {
  // Non-centred with respect to the common ridge scale; see gi_hd.stan.
  real alpha_raw;
  vector[P] gamma_raw;
  vector[P] K_raw;

  // Environment covariances and their common centre.
  array[E] cholesky_factor_cov[P] L_Sigma;
  cholesky_factor_corr[P] L_R_bar;
  vector<lower=0>[P] sd_bar;
  real<lower=0> nu_excess;              // nu = P + 1 + nu_excess

  // Target covariance.
  cholesky_factor_corr[P] L_R0;
  vector<lower=0>[P] sd0;

  // Conditional error variance, referenced to Sigma_bar.  The lower bound
  // keeps every v_e and S_0^2 positive; because Sigma_bar is now a parameter
  // the bound moves with it, which is exactly the dependence the plug-in
  // model was missing.
  real<lower=0> v_raw;

  array[E] vector[P] mu;
  cholesky_factor_corr[P] L_R;
  vector<lower=0>[P] sd_mu;

  real<lower=0> tau2_num;
  real<lower=0> tau2_den;
}

transformed parameters {
  real<lower=P + 1> nu = P + 1 + nu_excess;
  matrix[P, P] L_Sigma_bar = diag_pre_multiply(sd_bar, L_R_bar);
  matrix[P, P] L_Sigma0 = diag_pre_multiply(sd0, L_R0);

  real<lower=0> tau2 = tau2_num / tau2_den;
  real<lower=0> tau = sqrt(tau2);
  real prior_scale = sqrt(v_raw) * tau;
  real ncp_mult = pow(prior_scale, ncp);

  real alpha = alpha_raw * ncp_mult;
  vector[P] gamma = gamma_raw * ncp_mult;
  vector[P] K = K_raw * ncp_mult;

  real base = quad_form_inv_chol(K, L_Sigma_bar);
  real max_excess;
  real v;
  real<lower=0> sigma_y_sq;
  real<lower=0> sigma_y;
  real<lower=0> S0;
  real<lower=0> sigma_cond;

  {
    real acc = quad_form_inv_chol(K, L_Sigma0) - base;
    for (e in 1:E) {
      acc = fmax(acc, quad_form_inv_chol(K, L_Sigma[e]) - base);
    }
    max_excess = fmax(acc, 0);
  }
  v = max_excess + v_raw;
  sigma_cond = sqrt(v);
  sigma_y_sq = v + base;
  sigma_y = sqrt(sigma_y_sq);
  S0 = sqrt(sigma_y_sq - quad_form_inv_chol(K, L_Sigma0));
}

model {
  // ---- Priors -----------------------------------------------------------
  // Jeffreys prior on the free part of the residual variance.
  target += -log(v_raw);

  tau2_num ~ gamma(a_tau, 1);
  tau2_den ~ gamma(b_tau, 1);

  {
    real raw_sd = pow(prior_scale, 1 - ncp);
    alpha_raw ~ normal(0, raw_sd);
    gamma_raw ~ normal(0, raw_sd);
    K_raw ~ normal(0, raw_sd);
  }

  // Covariance hierarchy.
  L_R_bar ~ lkj_corr_cholesky(eta_lkj);
  sd_bar ~ normal(0, sd_sigma_scale);
  L_R0 ~ lkj_corr_cholesky(eta_lkj);
  sd0 ~ normal(0, sd_sigma_scale);
  nu_excess ~ gamma(2, 0.1);            // weakly informative, mean 20

  // Sigma_e ~ InvWishart(nu, (nu - P - 1) Sigma_bar), so E[Sigma_e] = Sigma_bar.
  for (e in 1:E) {
    L_Sigma[e] ~ inv_wishart_cholesky(nu, sqrt(nu - P - 1) * L_Sigma_bar);
  }

  // Environment means.
  L_R ~ lkj_corr_cholesky(eta_lkj);
  sd_mu ~ normal(0, sd_mu_scale);
  mu ~ multi_normal_cholesky(hmu, diag_pre_multiply(sd_mu, L_R));

  // ---- Likelihood -------------------------------------------------------
  for (e in 1:E) {
    real ne = n_e[e];
    vector[P] b_e = chol_solve(L_Sigma[e], K);
    vector[P] g_e = gamma + b_e;
    real a_e = alpha - dot_product(b_e, mu[e]);
    real v_e = v - (dot_product(K, b_e) - base);
    real ss;

    // X_ei ~ N(mu_e, Sigma_e), through the scatter matrix.
    target += -ne * sum(log(diagonal(L_Sigma[e])))
              - 0.5 * dot_self(to_vector(
                        mdivide_left_tri_low(L_Sigma[e], chol_scatter[e])))
              - 0.5 * ne * quad_form_inv_chol(xbar[e] - mu[e], L_Sigma[e]);

    // Y_ei ~ N(a_e + g_e' X_ei, v_e), through sufficient statistics.
    ss = yy[e]
         - 2 * a_e * ysum[e]
         - 2 * dot_product(g_e, xy[e])
         + ne * square(a_e)
         + 2 * ne * a_e * dot_product(g_e, xbar[e])
         + quad_form_sym(xx[e], g_e);
    target += -0.5 * ne * (log(2 * pi()) + log(v_e)) - 0.5 * ss / v_e;
  }

  // Target covariates: X_0i ~ N(mu_0, Sigma_0).  mu_0 is integrated out under
  // a flat prior, which contributes the usual -0.5 log|Sigma_0 / n_0| term;
  // dropping it changes the posterior for Sigma_0 only at order 1 / n_0.
  if (n0 > 0) {
    target += -n0 * sum(log(diagonal(L_Sigma0)))
              - 0.5 * dot_self(to_vector(
                        mdivide_left_tri_low(L_Sigma0, chol_scatter0)));
  }
}

generated quantities {
  vector[N0] f0;
  vector[N0] Y0_pred;
  matrix[P, P] Sigma_bar = multiply_lower_tri_self_transpose(L_Sigma_bar);
  matrix[P, P] Sigma0 = multiply_lower_tri_self_transpose(L_Sigma0);

  for (i in 1:N0) {
    f0[i] = alpha + dot_product(gamma, X0[i])
            + dot_product(K, chol_solve(L_Sigma0, X0[i]' - xbar0));
    Y0_pred[i] = normal_rng(f0[i], S0);
  }
}
