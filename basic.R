set.seed('2025')

n_source <- 500
n0 <- 200

sdh <- 0.5

# Generate synthetic data
H <- rnorm(n_source, mean = 0, sd = sdh)  # Continuous confounder H
X <- 0.5 * H + rnorm(n_source, sd = 0.1) + 2  # Continuous predictor X, influenced by H
Y <- 1 * X - 2 * H + rnorm(n_source, sd = 0.01)  # Outcome Y, influenced negatively by X and positively by H

ols <- lm(Y ~ -1 + X)

summary(ols)
mean(X)

H0 <- rnorm(n0, mean = 0, sd = sdh)  # Continuous confounder H
X0 <- 0.5 * H0 + rnorm(n0, sd = 0.1) + 5  # Continuous predictor X, influenced by H
Y0 <- 1 * X0 - 2 * H0 + rnorm(n0, sd = 0.01)  # Outcome Y, influenced negatively by X and positively by H

# Prepare data for Stan
stan_data <- list(
  N = n_source,
  Z = rep(1, n_source),
  NZ = 1,
  X = X,
  Y = Y,
  hmu = mean(X),
  hsigma = sd(X),
  N0 = n0,
  X0 = X0,
  var_X = var(X),
  var_X0 = var(X0),
  mu0 = mean(X0),
  sigma = var(X)
)

library(rstan)
options(mc.cores = 2)

# Fit the Stan model
fit <- stan(
  file = "gi_pred.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed=2025
)

# Print the summary of the fit
print(fit)

# Figure 3a

pairs(fit, pars = c('beta', 'K'))

# Set up output to PNG file
png(
  filename = "stan.png",
  width = 5, height = 5, units = "in", res = 300, # 6x6 inches, 300 dpi
  pointsize = 10                                # optional: controls text size
)

# Plot
pairs(fit, pars = c("beta", "K"))

# Close device (important!)
dev.off()



# Figure 3b

posterior_samples <- rstan::extract(fit)

# Extract Y_pred (assuming it's a matrix with dimensions [iterations, N0])
Y_pred_samples <- posterior_samples$Y_pred

# Compute posterior means and credible intervals
Y_pred_means <- colMeans(Y_pred_samples)
Y_pred_lower <- apply(Y_pred_samples, 2, function(x) quantile(x, 0.025))
Y_pred_upper <- apply(Y_pred_samples, 2, function(x) quantile(x, 0.975))

sum(Y_pred_upper > Y0 & Y_pred_lower < Y0) / n0

# Load necessary library
library(ggplot2)

# Combine into a data frame
data <- data.frame(
  X0 = X0,
  Y_pred_mean = Y_pred_means,
  Y_pred_lower = Y_pred_lower,
  Y_pred_upper = Y_pred_upper,
  Y0 = Y0
)

# ggplot with shaded ribbon
p <- ggplot(data, aes(x = X0)) +
  geom_ribbon(aes(ymin = Y_pred_lower, ymax = Y_pred_upper), fill = "grey", alpha = 0.5) +
  geom_point(aes(y = Y_pred_mean), color = "black", size = 2) +
  geom_point(aes(y = Y0), color = "chartreuse3", size = 2, shape = 17) +
  geom_abline(intercept = 0, slope = coef(ols), color = "red", linetype = "dashed", linewidth = 1) +
  geom_abline(intercept = 0, slope = mean(Y) / mean(X), color = "orange", linetype = "dashed", linewidth = 1) +
  labs(
    title = "95% predictive posterior intervals",
    x = "X",
    y = "Y"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line()
  )



# Set up output to PNG file
png(
  filename = "posterior_predictive_plot.png",
  width = 5, height = 5, units = "in", res = 300, # 6x6 inches, 300 dpi
  pointsize = 10                                # optional: controls text size
)

p

dev.off()


# Figure 1

olsint=lm(Y ~ X)
data <- data.frame(X, Y)

betagi=mean(Y)/mean(X)
kgi=mean((Y-betagi*X)*(X-mean(X)))


# Set up output to PNG file
png(
  filename = "firstf.png",
  width = 5, height = 4, units = "in", res = 300, # 6x6 inches, 300 dpi
  pointsize = 10                                # optional: controls text size
)

plot(X0, predict(olsint, newdata = data.frame(X = X0)), xlim = c(1, 6), ylim = c(-10, 8), pch = 20,
     ylab = 'Y', xlab = 'X')

points(X, Y)
abline(a=0,b=1, lty = 2)
points(X0, Y0, pch = 4)
points(X0,betagi*X0 + (kgi/var(X0))*(X0 - mean(X0)),col='red',pch=20)
rug(X0)

dev.off()

