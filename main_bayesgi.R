n=1000
set.seed(12025)
# Create a logarithmic grid of lambdas from 0.01 to 0.1
lambda_min <- 0.1
lambda_max <- 1
n_points <- 9  # Number of points in the grid

post_beta=matrix(0,nrow=n_points,ncol=4000)
post_k=matrix(0,nrow=n_points,ncol=4000)

# Generate the logarithmic grid
log_lambda_grid <- seq(log10(lambda_min), log10(lambda_max), length.out = n_points)
lambdas <- 10^log_lambda_grid
# Simulate exogenous random variables
H <- rnorm(n, mean = 0, sd = 1)         # H ~ N(0, 1)
epsilon1 <- rnorm(n, mean = 0, sd = 1)  # eps1 ~ N(0, 1)
epsilon2 <- rnorm(n, mean = 0, sd = 1)  # eps2 ~ N(0, 1)
epsilon=H + 2 * epsilon2
for (i in 1:n_points) {

# Compute X and Y
X <- H + 2 *(lambdas[i] + epsilon1)
Y <- 2 * X + epsilon

ols=lm(Y~-1 + X)
ols
cov(X,epsilon)
# Prepare data for Stan
stan_data <- list(
  N = n,
  Z = rep(1,n),
  NZ = 1,
  X = X,
  Y = Y,
  hmu = mean(X),
  hsigma = sd(X)
)
library(rstan)
# Fit the Stan model
fit <- stan(
  file = "zero_prior.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)
library(rstan)
options(mc.cores = 2)



# Print the summary of the fit
print(fit)

#pairs(fit)


# Extract samples from the posterior distribution of beta
posterior_samples <- extract(fit)
beta_samples <- posterior_samples$beta
k_samples <- posterior_samples$k

post_beta[i,]=beta_samples
post_k[i,]=k_samples

#hist(beta_samples)
#abline(v=coef(ols))
}


post=post_beta
# Define colors and transparency
colors <- c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), rgb(0, 1, 0, 0.5), rgb(1, 0, 1, 0.5),rgb(1, 1, 1, 0.5))
colors=c(colors,colors,colors)
# Create an empty plot

png(
  filename = "post_beta.png",
  width = 5, height = 4, units = "in", res = 300, # 6x6 inches, 300 dpi
  pointsize = 10                                # optional: controls text size
)
hist(post[1, ], breaks = 30, col = colors[1], xlim = c(0,2.2), ylim = c(0, 12), main = "Histograms of posterior samples", xlab = "beta", ylab = "Frequency", border = NA,freq=FALSE)

# Add histograms for each row
for (i in 1:nrow(post)) {
  hist(post[i, ], breaks = 30, col = colors[i], add = TRUE, border = NA,freq=FALSE)
}

# Add a legend
legend("topleft", legend = paste("lambda", round(lambdas,2)), fill = colors, cex = 0.8)
abline(v=mean(Y)/mean(X))

dev.off()

