/* -----------------------------------------------------------------------
 * Bayesian Generative Invariance (BGI) -- multivariate, multi-environment.
 *
 * Sampling model (training environment e = 1, ..., E):
 *
 *   X_ei             ~ N_p(mu_e, Sigma_e)
 *   Y_ei | X_ei      ~ N(alpha + gamma' X_ei + K' Sigma_e^{-1} (X_ei - mu_e),
 *                        sigma_y^2 - K' Sigma_e^{-1} K)
 *
 * Target (test) environment e = 0:
 *
 *   Y_0 | X_0        ~ N(alpha + gamma' X_0 + K' Sigma_0^{-1} (X_0 - mu_0),
 *                        sigma_y^2 - K' Sigma_0^{-1} K)
 *
 * K = Cov(eps_Y, X) is a *covariance*, not a slope: the slope on the centred
 * covariates is Sigma_e^{-1} K.  The training likelihood and the target
 * predictive formula therefore use the same object, and sigma_y^2 =
 * Var(eps_Y) is the marginal (not the conditional) error variance.
 *
 * ---- Parameterisation of the variance -------------------------------------
 *
 * The free scale parameter is `v`, the conditional variance measured against
 * the pooled within-environment covariance Sigma_bar:
 *
 *   sigma_y^2 = v + K' Sigma_bar^{-1} K
 *   v_e       = sigma_y^2 - K' Sigma_e^{-1} K = v - delta_e
 *   S_0^2     = sigma_y^2 - K' Sigma_0^{-1} K = v - delta_0
 *   delta_j  := K' Sigma_j^{-1} K - K' Sigma_bar^{-1} K
 *
 * This is an exact reparameterisation of (sigma_y^2, K), but a much better
 * conditioned one.  Parameterising by sigma_y^2 directly requires
 * sigma_y^2 > max_j K' Sigma_j^{-1} K, a bound on the *level* of the
 * Mahalanobis norms.  Under strong confounding that level sits just below
 * sigma_y^2, so the largest of the E + 1 noisy plug-in estimates binds the
 * constraint and the fit shrinks K simply to keep one environment's
 * conditional variance positive -- which biases gamma, since gamma is the
 * within-environment slope minus Sigma^{-1} K.  Parameterising by v turns the
 * bound into one on the *differences* delta_j, which are of the order of the
 * covariance estimation error rather than of sigma_y^2, and which shrinkage
 * of the Sigma_hat_e towards Sigma_bar drives towards zero (see
 * R/covariance.R).
 *
 * The causal Mahalanobis condition K' Sigma_j^{-1} K < sigma_y^2 is thus
 * imposed exactly for every environment including the target.  It is not an
 * extra assumption on the data-generating process: it holds automatically at
 * the true K whenever (X, Y) is jointly non-degenerate Gaussian.
 *
 * ---- Computational notes --------------------------------------------------
 *
 * The training log-likelihood depends on the data only through
 * per-environment sufficient statistics, so one gradient evaluation costs
 * O(E p^3) instead of O(N p^2).  Every Sigma is plugged in and supplied as a
 * Cholesky factor, so all factorisations happen once, in R.
 * --------------------------------------------------------------------- */

functions {
  /**
   * Solve Sigma x = v given the lower Cholesky factor L of Sigma.
   */
  vector chol_solve(matrix L, vector v) {
    return mdivide_right_tri_low(mdivide_left_tri_low(L, v)', L)';
  }

  /**
   * Mahalanobis norm v' Sigma^{-1} v from the lower Cholesky factor L.
   */
  real quad_form_inv_chol(vector v, matrix L) {
    return dot_self(mdivide_left_tri_low(L, v));
  }

  /**
   * Lower bound for v: the largest excess Mahalanobis norm of K over the
   * pooled reference, floored at zero so that v itself stays positive.
   */
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
  int<lower=1> P;                       // number of covariates
  int<lower=1> E;                       // number of training environments
  array[E] int<lower=1> n_e;            // observations per environment
  int<lower=1> N;                       // sum(n_e)

  // Per-environment sufficient statistics of the training data.
  array[E] vector[P] xbar;              // mean of X within environment e
  array[E] matrix[P, P] xx;             // sum_i X_ei X_ei'  (symmetric)
  array[E] vector[P] xy;                // sum_i X_ei Y_ei
  vector[E] ysum;                       // sum_i Y_ei
  vector[E] yy;                         // sum_i Y_ei^2

  // Plugged-in covariances, supplied as lower Cholesky factors.
  array[E] matrix[P, P] L_Sigma;        // chol(Sigma_e), e = 1, ..., E
  matrix[P, P] L_Sigma_bar;             // chol(pooled within-environment Sigma)
  matrix[P, P] L_Sigma0;                // chol(Sigma_0), target domain

  // Hyperparameters.
  vector[P] hmu;                        // prior mean for mu_e (pooled mean of X)
  real<lower=0> eta_lkj;                // LKJ shape for corr(Sigma_mu);
                                        //   large values -> diagonal Sigma_mu
  real<lower=0> sd_mu_scale;            // half-normal scale for sd(Sigma_mu)
  real<lower=0> a_tau;                  // Beta-prime shapes for tau^2
  real<lower=0> b_tau;                  //   (a = b = 1/2 -> half-Cauchy on tau)

  // Prior on the conditional variance.  Zero selects the Jeffreys prior
  // p(v_raw) propto 1/v_raw used throughout the paper, which is improper.
  // A positive shape selects inv_gamma(v_prior_shape, v_prior_rate) instead,
  // which is proper and therefore samplable -- required by simulation-based
  // calibration, where parameters must be drawn from the prior.  It has no
  // effect on any analysis that leaves it at zero.
  real<lower=0> v_prior_shape;
  real<lower=0> v_prior_rate;

  // Degree of non-centring for (alpha, gamma, K): 1 fully non-centred,
  // 0 fully centred, values in between partially so.  Non-centring is what
  // makes the near-collinear regime of Section 2.3 samplable at all, but it
  // is slower when the likelihood is strongly informative, which is the
  // usual centred/non-centred trade-off.  See README.md.
  real<lower=0, upper=1> ncp;

  // Target domain.
  int<lower=0> N0;                      // number of target covariate vectors
  matrix[N0, P] X0;                     // target covariates
  vector[P] mu0;                        // target mean of X

  // Optional pointwise log-likelihood (for loo / WAIC).
  int<lower=0, upper=1> compute_log_lik;
  matrix[compute_log_lik ? N : 0, P] X;
  vector[compute_log_lik ? N : 0] Y;
  array[compute_log_lik ? N : 0] int<lower=1, upper=E> Z;
}

transformed data {
  array[E + 1] matrix[P, P] L_all;      // all factors constraining v
  matrix[N0, P] X0_centred_prec;        // rows: (Sigma_0^{-1} (X_0i - mu_0))'

  for (e in 1:E) {
    L_all[e] = L_Sigma[e];
  }
  L_all[E + 1] = L_Sigma0;

  for (i in 1:N0) {
    X0_centred_prec[i] = chol_solve(L_Sigma0, X0[i]' - mu0)';
  }
}

parameters {
  // (alpha, gamma, K) are non-centred with respect to their common prior
  // scale.  Centred, the half-Cauchy tau creates a funnel that the sampler
  // cannot traverse when the environment means are near-collinear -- exactly
  // the weak-identifiability regime of Section 2.3, where in the centred
  // parameterisation every fit failed its diagnostics.
  real alpha_raw;
  vector[P] gamma_raw;
  vector[P] K_raw;

  // Free part of the conditional error variance.  Declared before K so that
  // the prior scale below is well defined; the Mahalanobis excess is added
  // back in `transformed parameters`, which is what keeps every v_e and
  // S_0^2 strictly positive.
  real<lower=0> v_raw;

  array[E] vector[P] mu;                // environment-specific covariate means

  cholesky_factor_corr[P] L_R;          // correlation of the mu_e prior
  vector<lower=0>[P] sd_mu;             // scales of the mu_e prior

  real<lower=0> tau2_num;               // tau^2 = tau2_num / tau2_den is
  real<lower=0> tau2_den;               //   Beta-prime(a_tau, b_tau)
}

transformed parameters {
  real<lower=0> tau2 = tau2_num / tau2_den;
  real<lower=0> tau = sqrt(tau2);
  real prior_scale = sqrt(v_raw) * tau;  // ridge scale, as in Section 2
  // Partial non-centring: the raw parameters carry prior sd
  // prior_scale^(1 - ncp) and are multiplied by prior_scale^ncp, so the
  // implied prior on (alpha, gamma, K) is N(0, prior_scale^2) for every ncp.
  real ncp_mult = pow(prior_scale, ncp);

  real alpha = alpha_raw * ncp_mult;
  vector[P] gamma = gamma_raw * ncp_mult;
  vector[P] K = K_raw * ncp_mult;

  real base = quad_form_inv_chol(K, L_Sigma_bar);
  real<lower=0> v = v_raw + max_variance_excess(K, L_all, L_Sigma_bar);
  real<lower=0> sigma_cond = sqrt(v);   // training-domain conditional sd
  real<lower=0> sigma_y_sq = v + base;
  real<lower=0> sigma_y = sqrt(sigma_y_sq);
  // Target-domain predictive sd: S_0 = sqrt(sigma_y^2 - K' Sigma_0^{-1} K).
  real<lower=0> S0 = sqrt(sigma_y_sq - quad_form_inv_chol(K, L_Sigma0));
}

model {
  // ---- Priors -----------------------------------------------------------
  if (v_prior_shape > 0) {
    v_raw ~ inv_gamma(v_prior_shape, v_prior_rate);
  } else {
    // Jeffreys prior for the residual scale: p(v_raw) propto 1/v_raw.
    target += -log(v_raw);
  }

  tau2_num ~ gamma(a_tau, 1);
  tau2_den ~ gamma(b_tau, 1);

  // Ridge-type priors.  Together with the definitions above these imply
  // alpha, gamma, K ~ N(0, v_raw tau^2) -- the ridge prior of Section 2 with
  // the residual scale -- for any value of ncp.
  {
    real raw_sd = pow(prior_scale, 1 - ncp);
    alpha_raw ~ normal(0, raw_sd);
    gamma_raw ~ normal(0, raw_sd);
    K_raw ~ normal(0, raw_sd);
  }

  L_R ~ lkj_corr_cholesky(eta_lkj);
  sd_mu ~ normal(0, sd_mu_scale);       // half-normal: sd_mu is declared > 0
  mu ~ multi_normal_cholesky(hmu, diag_pre_multiply(sd_mu, L_R));

  // ---- Likelihood -------------------------------------------------------
  {
    for (e in 1:E) {
      real ne = n_e[e];
      vector[P] b_e = chol_solve(L_Sigma[e], K);    // slope on centred X
      vector[P] g_e = gamma + b_e;                  // slope on raw X
      real a_e = alpha - dot_product(b_e, mu[e]);   // effective intercept
      real v_e = v - (dot_product(K, b_e) - base);  // conditional variance
      real ss;

      // X_ei ~ N(mu_e, Sigma_e).  Up to a constant free of parameters this is
      // equivalent to mu_e | xbar_e ~ N(xbar_e, Sigma_e / n_e).
      target += -0.5 * ne * quad_form_inv_chol(xbar[e] - mu[e], L_Sigma[e]);

      // Y_ei ~ N(a_e + g_e' X_ei, v_e), through sufficient statistics:
      //   SS = sum y^2 - 2 a sum y - 2 g' sum(x y)
      //        + n a^2 + 2 n a g' xbar + g' (sum x x') g
      ss = yy[e]
           - 2 * a_e * ysum[e]
           - 2 * dot_product(g_e, xy[e])
           + ne * square(a_e)
           + 2 * ne * a_e * dot_product(g_e, xbar[e])
           + quad_form_sym(xx[e], g_e);
      target += -0.5 * ne * (log(2 * pi()) + log(v_e)) - 0.5 * ss / v_e;
    }
  }
}

generated quantities {
  vector[N0] f0;                        // target-domain conditional mean
  vector[N0] Y0_pred;                   // posterior predictive draws for Y_0
  matrix[P, P] Sigma_mu = multiply_lower_tri_self_transpose(
                            diag_pre_multiply(sd_mu, L_R));
  vector[compute_log_lik ? N : 0] log_lik;

  f0 = alpha + X0 * gamma + X0_centred_prec * K;
  for (i in 1:N0) {
    Y0_pred[i] = normal_rng(f0[i], S0);
  }

  if (compute_log_lik) {
    array[E] vector[P] b_all;
    array[E] real a_all;
    array[E] real sd_all;
    for (e in 1:E) {
      b_all[e] = chol_solve(L_Sigma[e], K);
      a_all[e] = alpha - dot_product(b_all[e], mu[e]);
      sd_all[e] = sqrt(v - (dot_product(K, b_all[e]) - base));
    }
    for (i in 1:N) {
      int e = Z[i];
      log_lik[i] =
        normal_lpdf(Y[i] | a_all[e] + dot_product(gamma + b_all[e], X[i]),
                    sd_all[e])
        + multi_normal_cholesky_lpdf(X[i]' | mu[e], L_Sigma[e]);
    }
  }
}
