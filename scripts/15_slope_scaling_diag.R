#!/usr/bin/env Rscript
## 15_slope_scaling_diag.R ------------------------------------------------
##
## What breaks `gi_hd_slope` on real data?
##
## The slope model was validated at `E = 12`, `p = 6`, `n_e = 300`
## (`scripts/09_sim_slope_vs_k.R`, 90 tasks, zero divergences, max `Rhat`
## 1.019). On all five real datasets it fails, with a minimum bulk ESS of about
## 2 across four chains on four of them — four chains sitting in four places.
##
## The obvious explanation, more covariates, is refuted by the runs themselves:
## BRFSS has the most covariates at `p = 13` and behaves best (max `Rhat`
## 1.07), while ACS PUMS has the fewest at `p = 7` and behaves worst (77.9).
## The ordering is close to reversed, so `p` is not the driver.
##
## Real data differs from that simulation on two axes at once, and the runs
## cannot separate them because no dataset varies one while holding the other:
##
##   E     12 in simulation; 34 to 52 in the real runs.  Suspected mechanism:
##         `sigma_y^2 = v_raw + max_e { b_e' Sigma_e b_e }` is built with
##         `fmax`, which is non-differentiable wherever the arg-max moves
##         between environments.  Each arg-max regime is effectively a separate
##         mode, and the number of candidates for the maximum grows with `E`.
##
##   n_e   300 in simulation; 1,164 to 172,962 in the real runs.  A very
##         informative likelihood pins each `b_e` tightly, so the `E` terms
##         competing for the maximum are themselves sharp, which would make the
##         mode structure worse rather than better.
##
## This varies the two on one dataset, holding the covariates, the fold and
## every prior fixed, so the cause can be attributed rather than guessed. If
## the `E = 12`, `n_e = 300` cell reproduces the simulation's clean behaviour
## and the failure appears as either axis grows, the mechanism is located.
##
## Usage, as a 12-task array (see slurm/15_slope_scaling.sh):
##   Rscript scripts/15_slope_scaling_diag.R --task=1 [--dataset=brfss]

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
task <- as.integer(parse_flag(cli, "task", "1"))
dataset <- parse_flag(cli, "dataset", "brfss")
data_path <- parse_flag(cli, "data", switch(
  dataset,
  brfss = file.path(dirname(root), "data", "brfss", "brfss2023_case.csv"),
  quiron = file.path(dirname(root), "old_code", "quiron", "quiron_final.csv"),
  acs_tract = file.path(dirname(root), "data", "candidates", "tract",
                        "acs_tract_2022.csv"),
  stop("Pass --data for dataset ", dataset)))
model_name <- parse_flag(cli, "model", "gi_hd_slope")
n_iter <- as.integer(parse_flag(cli, "iter", "1000"))
n_chains <- as.integer(parse_flag(cli, "chains", "4"))
cores <- as.integer(parse_flag(
  cli, "cores", Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(n_chains))))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "slope_scaling"))

## The grid.  `n_e = Inf` means every row of the environment.  `sd_b_scale`
## varies only at full `E`, to check whether a tighter prior on the spread of
## the slopes removes the problem without changing the geometry.
grid <- expand.grid(E_sub = c(12L, 26L, 52L), n_sub = c(300L, 3000L, NA),
                    sd_b = 1, stringsAsFactors = FALSE)
grid <- rbind(grid,
              data.frame(E_sub = 52L, n_sub = NA, sd_b = c(0.25, 0.05)),
              data.frame(E_sub = 12L, n_sub = NA, sd_b = 0.25))
if (task < 1L || task > nrow(grid)) {
  stop("--task must be in 1:", nrow(grid), call. = FALSE)
}
cfg <- grid[task, ]

set.seed(20260804 + task)
dat <- load_case_data(data_path, dataset = dataset)
sizes <- sort(table(dat$z), decreasing = TRUE)

## Hold out the largest environment throughout, so the fold is identical in
## every cell and only the training design changes.
target <- names(sizes)[1]
train_pool <- setdiff(names(sizes), target)
E_use <- min(cfg$E_sub, length(train_pool))
keep_env <- utils::head(train_pool, E_use)

idx <- which(dat$z %in% keep_env)
if (!is.na(cfg$n_sub)) {
  idx <- unlist(lapply(split(idx, dat$z[idx]), function(i) {
    if (length(i) > cfg$n_sub) sort(sample(i, cfg$n_sub)) else i
  }), use.names = FALSE)
}
i0 <- which(dat$z == target)
if (length(i0) > 1000L) i0 <- sort(sample(i0, 1000L))

model <- load_bgi_model(file.path(root, "stan", paste0(model_name, ".stan")),
                        cache_dir = file.path(root, "results", "compiled"))

cat(sprintf("\ntask %d | %s | %s | E %d | n_e %s | sd_b_scale %.2f\n",
            task, dataset, model_name, E_use,
            if (is.na(cfg$n_sub)) "all" else format(cfg$n_sub), cfg$sd_b))
cat(sprintf("N_train %d | p %d | median n_e %.0f\n",
            length(idx), ncol(dat$x), stats::median(table(dat$z[idx]))))

t0 <- Sys.time()
fit <- tryCatch(
  fit_bgi(dat$x[idx, , drop = FALSE], dat$y[idx], dat$z[idx],
          dat$x[i0, , drop = FALSE], model = model,
          ncp = 0, sd_b_scale = cfg$sd_b,
          chains = n_chains, iter = n_iter, cores = cores,
          adapt_delta = 0.95, seed = 1),
  error = function(e) {
    cat("FAILED: ", conditionMessage(e), "\n"); NULL
  })
if (is.null(fit)) quit(save = "no", status = 0)

summ <- rstan::summary(fit$stanfit)$summary
block <- function(b) {
  i <- grep(paste0("^", b, "(\\[|$)"), rownames(summ))
  if (length(i) == 0L) return(c(NA_real_, NA_real_))
  c(max(summ[i, "Rhat"], na.rm = TRUE), min(summ[i, "n_eff"], na.rm = TRUE))
}
g <- block("gamma"); sy <- block("sigma_y"); sd_b <- block("sd_b")
b_bar <- block("b_bar"); s0 <- block("S0")
pm <- predictive_metrics(fit$draws$y0_pred, dat$y[i0])

row <- data.frame(
  task = task, dataset = dataset, model = model_name,
  E = E_use, n_sub = if (is.na(cfg$n_sub)) Inf else cfg$n_sub,
  sd_b_scale = cfg$sd_b, N_train = length(idx), p = ncol(dat$x),
  rhat_gamma = g[1], ess_gamma = g[2],
  rhat_b_bar = b_bar[1], ess_b_bar = b_bar[2],
  rhat_sd_b = sd_b[1], ess_sd_b = sd_b[2],
  rhat_sigma_y = sy[1], ess_sigma_y = sy[2],
  rhat_S0 = s0[1], ess_S0 = s0[2],
  max_rhat = fit$diagnostics$max_rhat,
  min_ess = fit$diagnostics$min_ess_bulk,
  divergent = fit$diagnostics$n_divergent,
  coverage = pm$coverage,
  runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
  stringsAsFactors = FALSE
)
print(t(row))
bgi_write_csv(row, file.path(out_dir, sprintf("cell_%02d.csv", task)))
cat("\nwritten to ", file.path(out_dir, sprintf("cell_%02d.csv", task)), "\n",
    sep = "")
