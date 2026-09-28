#!/usr/bin/env Rscript
## 21_one_dim_illustration.R ----------------------------------------------
##
## The one-dimensional illustration of Supplement E, computed with the model
## the rest of the paper uses.
##
## Why this script exists
## ----------------------
## The figures in Supplement E were produced by an earlier, single-source
## script against a Stan file (`gi_pred.stan`) that is not part of this
## package.  That script wrote the conditional mean as
##
##     gamma * X + (K / var(X)) * (X - mu),
##
## passing `var_X` and `var_X0` in as data, so `K` there was the covariance
## and `K / var(X)` the slope.  The supplement, however, printed the model as
## `gamma * X + K * (X - mu)`, which makes `K` a slope and is not what was
## fitted: with that data-generating process the slope is
## `-0.25 / 0.0725 = -3.45`, not the `-0.25` the figure shows.  The number was
## right and the equation was not.
##
## `stan/gi_hd.stan` takes the same convention, explicitly: it forms
## `b_e = Sigma_e^{-1} K` and calls it the slope on centred X.  Recomputing
## the illustration with it therefore settles the discrepancy in favour of the
## figure, and puts the supplement on the same computational core as every
## other result.
##
## Two deliberate departures from the original script
## --------------------------------------------------
## First, two training environments rather than one.  The old script dropped
## the intercept so that a single environment would identify the rest;
## `gi_hd.stan` carries an intercept, and with it Assumption 1 needs the
## augmented means to span R^2, hence E >= 2.  Keeping the intercept and
## adding an environment is the honest way round, and costs the illustration
## nothing: its point is that instrumental variables recover the causal slope
## and still predict badly once the domain moves, which needs one covariate
## and a shifted target, not a single source.
##
## Second, the target domain is more dispersed than the training ones.  In the
## original script `X0` had the same idiosyncratic variance as `X`, so
## `S_0` and the training conditional scale coincided and the figure could not
## show the correction this paper is about.  Here `Var(X_0) > Var(X_e)`, as in
## every other experiment in the package, which separates them.
##
## Usage:
##   Rscript scripts/21_one_dim_illustration.R [--out=DIR] [--seed=N]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(parse_flag(cli, "seed", "2025"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "one_dim"))
## Two environments at the boundary of Assumption 1 give the same funnel the
## p = 2 cells of the dimension sweep meet, so the step size is adapted the
## same way there.
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.99"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Data-generating process -------------------------------------------
## H is the hidden confounder, loading on both X and Y, and it is the only
## source of dependence between X and the error of Y.

set.seed(seed)

sd_h <- 0.5
sd_x_train <- 0.1
sd_x_target <- 0.3          # target domain is the more dispersed one
sd_y <- 0.01
gamma_true <- 1
env_shift <- c(2, 3.5)      # two training environments
target_shift <- 5
n_e <- 500
n0 <- 200

draw <- function(n, shift, sd_x) {
  h <- stats::rnorm(n, 0, sd_h)
  x <- 0.5 * h + stats::rnorm(n, 0, sd_x) + shift
  y <- gamma_true * x - 2 * h + stats::rnorm(n, 0, sd_y)
  list(x = x, y = y)
}

train <- lapply(env_shift, function(s) draw(n_e, s, sd_x_train))
x <- matrix(unlist(lapply(train, `[[`, "x")), ncol = 1)
y <- unlist(lapply(train, `[[`, "y"))
z <- rep(seq_along(env_shift), each = n_e)

target <- draw(n0, target_shift, sd_x_target)
x0 <- matrix(target$x, ncol = 1)
y0 <- target$y

## ---- What the truth implies --------------------------------------------
## eps_Y = -2H + noise, so Cov(eps_Y, X) = -2 * 0.5 * Var(H) and the slope in
## an environment is that covariance divided by the variance of X there.

var_h <- sd_h^2
k_true <- -2 * 0.5 * var_h
var_x_train <- 0.25 * var_h + sd_x_train^2
var_x_target <- 0.25 * var_h + sd_x_target^2
sigma_y_true <- sqrt(4 * var_h + sd_y^2)
s0_true <- sqrt(sigma_y_true^2 - k_true^2 / var_x_target)
cond_true <- sqrt(sigma_y_true^2 - k_true^2 / var_x_train)

cat("Truth implied by the design\n")
cat(sprintf("  gamma                     %8.4f\n", gamma_true))
cat(sprintf("  K (a covariance)          %8.4f\n", k_true))
cat(sprintf("  slope in training, K/var  %8.4f\n", k_true / var_x_train))
cat(sprintf("  slope in target,   K/var  %8.4f\n", k_true / var_x_target))
cat(sprintf("  sigma_Y                   %8.4f\n", sigma_y_true))
cat(sprintf("  training conditional sd   %8.4f\n", cond_true))
cat(sprintf("  S_0 (target scale)        %8.4f\n", s0_true))

## ---- Fit ----------------------------------------------------------------

model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))
fit <- fit_bgi(x = x, y = y, z = z, x0 = x0, model = model,
               cov_method = "pooled", chains = 4, iter = 2000,
               seed = seed, cores = 4, adapt_delta = adapt_delta)

post <- function(v) c(mean = mean(v), lo = stats::quantile(v, 0.025, names = FALSE),
                      hi = stats::quantile(v, 0.975, names = FALSE))
g <- post(as.vector(fit$draws$gamma))
k <- post(as.vector(fit$draws$K))
s0 <- post(fit$draws$S0)
sc <- post(fit$draws$sigma_cond)

cat("\nPosterior, mean [95% interval]\n")
cat(sprintf("  gamma      %8.4f  [%7.4f, %7.4f]   truth %7.4f\n",
            g[1], g[2], g[3], gamma_true))
cat(sprintf("  K          %8.4f  [%7.4f, %7.4f]   truth %7.4f\n",
            k[1], k[2], k[3], k_true))
cat(sprintf("  S_0        %8.4f  [%7.4f, %7.4f]   truth %7.4f\n",
            s0[1], s0[2], s0[3], s0_true))
cat(sprintf("  sigma_cond %8.4f  [%7.4f, %7.4f]   truth %7.4f\n",
            sc[1], sc[2], sc[3], cond_true))
cat(sprintf("\ndivergences %d | max Rhat %.4f | min bulk ESS %.0f\n",
            fit$diagnostics$n_divergent, fit$diagnostics$max_rhat,
            fit$diagnostics$min_ess_bulk))

## ---- Baselines ----------------------------------------------------------

ols <- fit_ols(x, y, x0)
iv <- tryCatch(fit_iv_2sls(x, y, z, x0), error = function(e) NULL)
cat(sprintf("\nOLS slope %.4f | IV slope %s | truth %.4f\n",
            ols$gamma,
            if (is.null(iv)) "not identified" else sprintf("%.4f", iv$gamma),
            gamma_true))

## ---- Coverage -----------------------------------------------------------

pm <- predictive_metrics(fit$draws$y0_pred, y0)
cat(sprintf("predictive coverage %.3f at interval score %.3f\n",
            pm$coverage, pm$interval_score))

saveRDS(list(fit = fit$draws, diagnostics = fit$diagnostics,
             ols = ols, iv = iv, x = x, y = y, z = z, x0 = x0, y0 = y0,
             truth = list(gamma = gamma_true, k = k_true, s0 = s0_true,
                          cond = cond_true)),
        file.path(out_dir, "one_dim.rds"))
## ---- Figures ------------------------------------------------------------
## Greyscale, every series named in a legend, text at 10pt or larger: the
## journal prints in black and white.

gam <- as.vector(fit$draws$gamma)
kk <- as.vector(fit$draws$K)

png(file.path(out_dir, "stan.png"), width = 5, height = 5, units = "in",
    res = 300, pointsize = 11)
op <- par(no.readonly = TRUE)
layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),
       widths = c(3, 1), heights = c(1, 3))
par(mar = c(4.2, 4.2, 0.4, 0.4))
plot(gam, kk, pch = 20, cex = 0.35, col = grey(0.25),
     xlab = expression(gamma), ylab = expression(K))
abline(v = gamma_true, h = k_true, lty = 2)
par(mar = c(0.3, 4.2, 0.6, 0.4))
hg <- hist(gam, breaks = 40, plot = FALSE)
barplot(hg$counts, space = 0, col = grey(0.75), border = grey(0.35), axes = FALSE)
par(mar = c(4.2, 0.3, 0.4, 0.6))
hk <- hist(kk, breaks = 40, plot = FALSE)
barplot(hk$counts, space = 0, col = grey(0.75), border = grey(0.35),
        axes = FALSE, horiz = TRUE)
par(op); dev.off()

png(file.path(out_dir, "posterior_predictive_plot.png"), width = 6.2,
    height = 4.4, units = "in", res = 300, pointsize = 11)
par(mar = c(4.4, 4.4, 0.6, 0.6))
lo <- apply(fit$draws$y0_pred, 2, stats::quantile, 0.025)
hi <- apply(fit$draws$y0_pred, 2, stats::quantile, 0.975)
mid <- colMeans(fit$draws$y0_pred)
o <- order(x0[, 1])
plot(x0[, 1], y0, type = "n", xlab = "X", ylab = "Y",
     ylim = range(lo, hi, y0))
polygon(c(x0[o, 1], rev(x0[o, 1])), c(lo[o], rev(hi[o])),
        col = grey(0.85), border = NA)
points(x0[, 1], mid, pch = 20, cex = 0.7)
points(x0[, 1], y0, pch = 2, cex = 0.7)
abline(a = ols$alpha, b = ols$gamma, lty = 2, lwd = 1.6)
if (!is.null(iv)) abline(a = iv$alpha, b = iv$gamma, lty = 3, lwd = 1.6)
## Bottom left is the only region the band and the responses leave clear.
legend("bottomleft", bty = "n", cex = 0.9, inset = c(0, 0.01),
       legend = c("predictive mean", "unseen responses",
                  "95% predictive band", "least squares", "instrumental variables"),
       pch = c(20, 2, 15, NA, NA), lty = c(NA, NA, NA, 2, 3),
       lwd = c(NA, NA, NA, 1.6, 1.6),
       col = c("black", "black", grey(0.85), "black", "black"))
dev.off()
## Third figure: what the two procedures predict in the target domain, against
## the causal truth.  Same fit, same conventions as the caption in Supplement E.
png(file.path(out_dir, "firsf.png"), width = 5.6, height = 4.4, units = "in",
    res = 300, pointsize = 11)
par(mar = c(4.4, 4.4, 0.6, 0.6))
f0 <- colMeans(fit$draws$f0)
plot(x[, 1], y, pch = 1, cex = 0.5, col = grey(0.45),
     xlim = range(x[, 1], x0[, 1]), ylim = range(y, y0, ols$pred_mean),
     xlab = "X", ylab = "Y")
points(x0[, 1], y0, pch = 4, cex = 0.7)
points(x0[, 1], ols$pred_mean, pch = 2, cex = 0.7)
points(x0[, 1], f0, pch = 20, cex = 0.7)
abline(a = 0, b = gamma_true, lty = 2)
rug(x0[, 1])
legend("topleft", bty = "n", cex = 0.9,
       legend = c("training sample", "test responses",
                  "least squares prediction", "our prediction", "causal truth"),
       pch = c(1, 4, 2, 20, NA), lty = c(NA, NA, NA, NA, 2),
       col = c(grey(0.45), "black", "black", "black", "black"))
dev.off()
message("Figures written to ", out_dir)
message("\nWritten to ", out_dir)
