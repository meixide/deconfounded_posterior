// Stan model for Bayesian Generative Invariance
data {
  int<lower=0> N;                // number of observations
  array[N] int<lower=1> Z;   // instrument 
  int<lower=0> NZ;   // number of environments (randomization)

  vector[N] X;
  vector[N] Y;                   // outcome
  
  // Prior parameters
  real hmu;  
  real<lower=0> hsigma;          // covariance matrix for multivariate normal
  
}

parameters {

  real beta;                   
  real k;
  array[NZ] real mu; 

}

model {
  mu ~ normal(hmu, hsigma);


  // Likelihood
  for (i in 1:N) {
    X[i] ~ normal(mu[Z[i]],1); 
    Y[i] ~ normal(beta * X[i] + k * (X[i] - mu[Z[i]]), 1);
  }
  
}

