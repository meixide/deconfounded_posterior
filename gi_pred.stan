data {
  int<lower=0> N;                // number of observations
  array[N] int<lower=1> Z;       // instrument 
  int<lower=0> NZ;               // number of environments (randomization)

  vector[N] X;
  vector[N] Y;                   // outcome
  
  // Prior parameters
  real hmu;  
  real<lower=0> hsigma;          // prior std for mu

   // New data points for prediction
  int<lower=0> N0;               // number of new X0 points
  vector[N0] X0;                 // new X values for prediction
  
  // External parameters
  real var_X;                    // variance of X (given externally)
  real var_X0;                   // variance of X0 (given externally)
  real mu0;                      // mean of X0 (given externally)
  real<lower=0> sigma; // group-specific standard deviations

}

parameters {
  real beta;                   
  real k;
  real<lower=0> sigmay;
  array[NZ] real mu; 
}

transformed parameters {
  real K = var_X*k;
}

model {
  // Priors
  mu ~ normal(hmu, 1);


  // Likelihood
  for (i in 1:N) {
    X[i] ~ normal(mu[Z[i]], sigma); 
    Y[i] ~ normal(beta * X[i] + k * (X[i] - mu[Z[i]]), sigmay);
  }
}
generated quantities {
  // Predictive distribution for Y given X0
  vector[N0] Y_pred;
  for (i in 1:N0) {
    // Conditional mean for Y given X0[i]
    real cond_mean = beta * X0[i] + k * (var_X / var_X0) * (X0[i] - mu0);
    // Sample Y_pred[i] from the normal distribution
    Y_pred[i] = normal_rng(cond_mean, sigmay);
  }
}

