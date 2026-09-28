#!/usr/bin/env Rscript
## 13_slope_convergence_diag.R --------------------------------------------
##
## Why does `gi_hd_slope` fail to mix on the case-study data, and what does it
## contaminate?
##
## The leave-one-environment-out run of the slope model reports a maximum
## `Rhat` of 42.4 on quiron with a minimum bulk ESS of 2.002 across four
## chains.  An ESS of two is not slow mixing; it is four chains sitting in four
## places and never moving between them.  Covariance conditioning does not
## explain it — the correlation between `log Rhat` and `log` condition number
## over the 51 folds is 0.016, and the best and worst folds have the same
## conditioning.
##
## The suspected mechanism is the variance anchor.  `gi_hd_slope.stan` sets
##
##     sigma_y^2 = v_raw + max_e { b_e' Sigma_e b_e }   (and the target term)
##
## through `fmax`.  That is non-differentiable wherever the arg-max switches
## between environments.  At the `E = 12` of `scripts/09_sim_slope_vs_k.R` the
## boundary is rarely contested; at `E = 52` with `E x p = 624` free slope
## deviations it is contested constantly, and each arg-max regime is a separate
## mode.
##
## What matters for the paper is *which* parameters this reaches.  If only the
## anchored scale parameters (`sigma_y`, `S0`) fail while `gamma` mixes, the
## damage is confined to the predictive width and the causal claims survive.
## If `gamma` fails, nothing from this model can be reported on this data.
## This script answers that by refitting one fold and printing `Rhat` and ESS
## per parameter block, plus the per-chain posterior means that reveal whether
## the chains are in genuinely different places.
##
## Usage:
##   Rscript scripts/13_slope_convergence_diag.R --dataset=quiron \
##           --target="Balears, Illes" [--model=gi_hd_slope] [--iter=2000]

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
dataset <- parse_flag(cli, "dataset", "quiron")
data_path <- parse_flag(cli, "data", switch(
  dataset,
  quiron = file.path(dirname(root), "old_code", "quiron", "quiron_final.csv"),
  brfss = file.path(dirname(root), "data", "brfss", "brfss2023_case.csv"),
  stop("Pass --data for dataset ", dataset)))
model_name <- parse_flag(cli, "model", "gi_hd_slope")
n_iter <- as.integer(parse_flag(cli, "iter", "2000"))
n_chains <- as.integer(parse_flag(cli, "chains", "4"))
ncp <- as.numeric(parse_flag(cli, "ncp", "0"))
sd_b_scale <- as.numeric(parse_flag(cli, "sd-b-scale", "1"))
adapt_delta <- as.numeric(parse_flag(cli, "adapt-delta", "0.95"))
## `auto` replaces the improper Jeffreys prior on the residual variance with a
## proper inverse-gamma centred on the pooled within-environment residual
## variance.  The suspicion is that `v_raw` is only weakly identified under the
## slope parameterisation -- it is the difference `v_e - b0_quad + q_e` of two
## large, freely-growing quantities -- so the `1/v_raw` spike at zero dominates
## and drags the predictive scale to zero.  Under `gi_hd` the same prior is
## harmless because `v_raw` there simply *is* the residual variance.
v_prior <- parse_flag(cli, "v-prior", "off")
seed <- as.integer(parse_flag(cli, "seed", "1"))
max_target <- as.integer(parse_flag(cli, "max-target", "2000"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(n_chains))))

dat <- load_case_data(data_path, dataset = dataset)
target <- parse_flag(cli, "target",
                     names(sort(table(dat$z), decreasing = TRUE))[1])
if (!target %in% dat$z) {
  stop("No environment called ", target, call. = FALSE)
}

model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

is0 <- dat$z == target
i0 <- which(is0)
if (length(i0) > max_target) i0 <- sort(sample(i0, max_target))

cat(sprintf("\n%s | model %s | target %s | N_train %d | E_train %d | p %d\n",
            dataset, model_name, target, sum(!is0),
            length(unique(dat$z[!is0])), ncol(dat$x)))
cat(sprintf("%d chains x %d iterations | ncp %.0f | sd_b_scale %.2f | adapt_delta %.2f\n\n",
            n_chains, n_iter, ncp, sd_b_scale, adapt_delta))

## Pooled within-environment residual variance: the scale the prior is
## centred on.  From per-environment OLS fits, so it makes no reference to the
## model under test.
resid_var <- {
  xtr <- dat$x[!is0, , drop = FALSE]; ytr <- dat$y[!is0]; ztr <- dat$z[!is0]
  rs <- unlist(lapply(split(seq_along(ztr), ztr), function(i) {
    if (length(i) <= ncol(xtr) + 2L) return(numeric(0))
    stats::residuals(stats::lm(ytr[i] ~ xtr[i, , drop = FALSE]))
  }), use.names = FALSE)
  stats::var(rs)
}
## inverse-gamma(3, 2 s^2) has mean s^2 and finite variance.
v_shape <- if (identical(v_prior, "auto")) 3 else 0
v_rate <- if (identical(v_prior, "auto")) 2 * resid_var else 0
cat(sprintf("pooled within-env residual sd %.3f | v-prior %s%s | seed %d\n\n",
            sqrt(resid_var), v_prior,
            if (identical(v_prior, "auto"))
              sprintf(" -> inv-gamma(%.0f, %.3f)", v_shape, v_rate) else "",
            seed))

fit <- fit_bgi(dat$x[!is0, , drop = FALSE], dat$y[!is0], dat$z[!is0],
               dat$x[i0, , drop = FALSE],
               model = model, ncp = ncp, sd_b_scale = sd_b_scale,
               v_prior_shape = v_shape, v_prior_rate = v_rate,
               chains = n_chains, iter = n_iter, cores = cores,
               adapt_delta = adapt_delta, seed = seed)

summ <- rstan::summary(fit$stanfit)$summary
blocks <- c("alpha", "gamma", "K", "b_bar", "sd_b", "sigma_y", "sigma_cond",
            "S0", "tau", "v_raw", "mu")
cat("=== Rhat and ESS by parameter block ===\n")
cat(sprintf("%-12s %6s %10s %10s %10s\n", "block", "n", "max Rhat",
            "min ESS", "median ESS"))
for (b in blocks) {
  idx <- grep(paste0("^", b, "(\\[|$)"), rownames(summ))
  if (length(idx) == 0L) next
  cat(sprintf("%-12s %6d %10.3f %10.0f %10.0f\n", b, length(idx),
              max(summ[idx, "Rhat"], na.rm = TRUE),
              min(summ[idx, "n_eff"], na.rm = TRUE),
              stats::median(summ[idx, "n_eff"], na.rm = TRUE)))
}

## Per-chain means for the worst offenders: if the chains disagree on a
## parameter's location rather than merely exploring it slowly, they are in
## different modes and more iterations will not help.
cat("\n=== per-chain posterior means, worst-mixing scalars ===\n")
scalars <- rownames(summ)[grep("^(sigma_y|sigma_cond|S0|tau|alpha|v_raw)$",
                               rownames(summ))]
worst_g <- rownames(summ)[grep("^gamma\\[", rownames(summ))]
worst_g <- worst_g[order(-summ[worst_g, "Rhat"])][seq_len(min(3, length(worst_g)))]
arr <- as.array(fit$stanfit, pars = c(scalars, worst_g))
for (nm in dimnames(arr)[[3]]) {
  cm <- colMeans(arr[, , nm, drop = FALSE][, , 1])
  cat(sprintf("  %-14s Rhat %7.3f | chain means %s\n", nm,
              summ[nm, "Rhat"],
              paste(sprintf("%.3f", cm), collapse = "  ")))
}

sp <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
cat(sprintf("\ndivergences %d | max treedepth hits %d | runtime %.0f s\n",
            fit$diagnostics$n_divergent,
            sum(vapply(sp, function(s) sum(s[, "treedepth__"] >= 10), numeric(1))),
            fit$diagnostics$runtime_sec))

cat("\nRead the block table first. If `gamma` mixes and only the anchored\n")
cat("scale parameters do not, the causal conclusions survive and the damage\n")
cat("is confined to the predictive width. If `gamma` does not mix, nothing\n")
cat("from this model can be quoted on this data.\n")
