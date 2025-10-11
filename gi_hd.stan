data {
  int<lower=0> N;                // number of observations
  array[N] int<lower=1> Z;       // instrument 
  int<lower=0> NZ;               // number of environments (randomization)
  int<lower=0> P;                // number of covariates

  matrix[N, P] X;                // N x P matrix of covariates
  vector[N] Y;                   // outcome
  
  // Prior parameters
  vector[P] hmu;             // prior mean for each environment-covariate combination

  // New data points for prediction
  int<lower=0> N0;               // number of new X0 points
  matrix[N0, P] X0;              // new X values for prediction
  
  // External parameters
  matrix[P, P] ivar_X0;  // variance matrix of X0
  vector[P] mu0;                 // mean of X0 
    matrix[P, P] var_X[NZ];  // Array of N matrices
    matrix[P,P] avg_var_X;
 
}

parameters {
  vector[P] beta;                // coefficient vector
  real betaintercept;
  vector[P] k;                   // interaction coefficient vector
  real log_sigmay;          // outcome variance
  real<lower=0> uppertau2;
  real<lower=0> lowertau2;
  matrix[NZ, P] mu;              // environment-specific means for each covariate
  corr_matrix[P] R;      // Correlation matrix
  vector<lower=0>[P] sds; // Standard deviations
}

transformed parameters {
  cov_matrix[P] Sigma;
  Sigma = quad_form_diag(R, sds); // Convert to covariance matrix
  real<lower=0> sigmay = exp(log_sigmay);
  real<lower=0> tau2= uppertau2 / lowertau2;
 real<lower=0> tau=sqrt(tau2);
}

model {
  // Priors
  uppertau2 ~ gamma(0.5, 1);  
  lowertau2 ~ gamma(0.5, 1);  
  
  for (z in 1:NZ) {
    mu[z,] ~ multi_normal(hmu,Sigma);  // Normal prior for each environment-covariate mean
  }
  beta ~ normal(0,sigmay*tau);
  betaintercept ~ normal(0,sigmay*tau);
  k ~ normal(0,sigmay*tau);
  

  // Likelihood
  for (i in 1:N) {
    
        X[i,] ~ multi_normal(mu[Z[i],], var_X[Z[i]]); 

    // Outcome model with multivariate interaction
    Y[i] ~ normal(betaintercept +
      dot_product(beta, X[i,]) +   
      dot_product(k, X[i,] - mu[Z[i],]), 
      sigmay
    );
  }
}

generated quantities {
  // Predictive distribution for Y given X0
  vector[N0] Y_pred;
  
 
  
  // Matrix multiplication for variance ratio
  matrix[P, P] var_ratio = avg_var_X * ivar_X0;
  
  for (i in 1:N0) {
    // Conditional mean for Y given X0[i]
    real cond_mean = betaintercept + dot_product(beta, X0[i,]) + 
                     dot_product(k, var_ratio * (X0[i,]' - mu0));
    
    // Sample Y_pred[i] from the normal distribution
    Y_pred[i] = normal_rng(cond_mean, sigmay);
  }
}

