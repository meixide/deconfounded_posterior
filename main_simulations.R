#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Required arguments: n p b CH num_simulations\n  n: number of observations\n  p: number of covariates\n  b: number of cores for parallel\n num_simulations: number of simulations")
}

p <- as.numeric(args[2])  # number of covariates

n <- floor(as.numeric(args[1])/(p+1))  # number of observations per environment
b <- as.numeric(args[3])  # number of cores for parallel processing
num_simulations <- as.numeric(args[4])  # number of simulations

folder_name <- sprintf("np_%d_%d", as.numeric(args[1]), p)


# Create the folder if it doesn't exist
if (!dir.exists(folder_name)) {
  dir.create(folder_name, recursive = TRUE, showWarnings = TRUE)
  if (dir.exists(folder_name)) {
    message(sprintf("Folder '%s' created successfully.", folder_name))
  } else {
    stop(sprintf("Failed to create folder '%s'.", folder_name))
  }
} else {
  message(sprintf("Folder '%s' already exists.", folder_name))
}

print(num_simulations)
print('Available cores:')
library(parallel)
library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)

print("Available cores:")
print(parallel::detectCores())

list_to_matrix <- function(input_list) {
  if (length(unique(sapply(input_list, length))) > 1) {
    stop("All vectors in the list must have the same length")
  }
  matrix_output <- do.call(rbind, input_list)
  return(matrix_output)
}

# Function to run one simulation
run_simulation <- function(sim_num) {
  q = 3  # number of hidden confounders
  num_environments = p + 1  # number of different environments
  
  # Fixed true coefficient vector
  true_beta <- c(-0.24, 0.87, -0.91, -0.66, 0.06, -0.13, -0.1, -0.72, -0.54, -0.71)[1:p]
  
  # Create correlation matrix for E
  Sigma_E <- matrix(0.5, nrow = p, ncol = p)
  diag(Sigma_E) <- 1
  base_mean <- seq(-1, 1, length.out = p)

  # Generate environment-specific mean vectors
  env_means <- lapply(1:num_environments, function(env) {
    base_mean + 1 * runif(p, -1, 1)
  })
  
  e <- rnorm(n)
  Psi <- matrix(rnorm(p * q), nrow = p, ncol = q)
  phi <- rnorm(q)
  
  # Container for environment-specific results
  sigmaw = list()
  
  # Generate initial environment
  E <- mvtnorm::rmvnorm(n, mean = env_means[[1]], sigma = Sigma_E)
  H <- matrix(rnorm(n * q), nrow = n, ncol = q)
  X <- E + H %*% t(Psi)
  Y <- X %*% true_beta + H %*% phi + e
  sigmaw[[1]] = var(X)
  
  # Generate additional environments
  for (env in 2:num_environments) {   
    e <- rnorm(n)
    H <- matrix(rnorm(n * q), nrow = n, ncol = q)
    E <- mvtnorm::rmvnorm(n, mean = env_means[[env]], sigma = Sigma_E)
    Xe <- E + H %*% t(Psi)
    X <- rbind(X, Xe)
    Ye <- Xe %*% true_beta + H %*% phi + e
    Y <- c(Y, Ye)
    sigmaw[[env]] = var(Xe)
  }
  
  Z = rep(1:num_environments, each=n)
  
  # Generate test data
  e <- rnorm(n)
  H <- matrix(rnorm(n * q), nrow = n, ncol = q)
  E <- mvtnorm::rmvnorm(n, 
                        mean = base_mean + 3 * runif(p, -1, 1), 
                        sigma = Sigma_E + 0.5*diag(p))
  X0 <- E + H %*% t(Psi)
  Y0 <- X0 %*% true_beta + H %*% phi + e
  
  options(mc.cores = 1)
  rstan_options(threads_per_chain = 1)  # No within-chain threading
  # Prepare Stan data
  stan_data <- list(
    N = n*num_environments,
    Z = Z,
    NZ = num_environments,
    P = p,
    X = X,
    Y = Y,
    hmu = Reduce('+',env_means)/num_environments,
    N0 = n,
    X0 = X0,
    ivar_X0 = solve(var(X0)),
    mu0 = apply(X0, 2, mean),
    var_X = sigmaw,
    avg_var_X = Reduce('+',sigmaw)/num_environments
  )
  message(sprintf("Started job %d", sim_num))
  
  # Fit Stan model
  fit <- stan(
    file = "gi_hd.stan",
    data = stan_data,
    chains = 4,
    iter = 1000,
    warmup = 500
  )
  message(sprintf("Finished job %d", sim_num))
  
  # Extract predictions
  posterior_samples <- extract(fit)
  Y_pred_samples <- posterior_samples$Y_pred
  Y_pred_lower <- apply(Y_pred_samples, 2, function(x) quantile(x, 0.025))
  Y_pred_upper <- apply(Y_pred_samples, 2, function(x) quantile(x, 0.975))
  
  # Calculate coverage for Stan model
  covours <- mean(Y0 >= Y_pred_lower & Y0 <= Y_pred_upper)
  
  # Calculate coverage for OLS
  ols <- lm(Y~., data=data.frame(Y,X))
  X0df <- as.data.frame(X0)
  colnames(X0df) <- colnames(data.frame(X))
  predictiols <- predict(ols, newdata = X0df, interval = "prediction")
  covols <- mean(Y0 >= predictiols[,"lwr"] & Y0 <= predictiols[,"upr"])
  
  # Save results to file
  filename <- sprintf("%s/simulationnn_%d.txt", folder_name, sim_num)
  write.table(
    data.frame(covours = covours, covols = covols),
    file = filename,
    row.names = FALSE,
    sep = "\t"
  )
  
  return(c(covours, covols))
}


results <- mclapply(1:num_simulations, run_simulation, mc.cores = b)
results_matrix <- do.call(rbind, results)
colnames(results_matrix) <- c("covours", "covols")
print(colMeans(results_matrix))