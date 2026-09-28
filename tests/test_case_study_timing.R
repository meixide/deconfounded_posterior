#!/usr/bin/env Rscript
## test_case_study_timing.R -----------------------------------------------
##
## How long does one leave-one-environment-out fold actually take, and which
## parameterisation should the case study use?
##
## The case-study design (E = 51, p = 12, N > 400,000) is an order of
## magnitude larger than anything the simulations exercised, and a first smoke
## test did not finish a single fold in 25 minutes.  Rather than guess, this
## measures where the time goes.
##
## The hypothesis under test is that non-centring is the culprit.  Non-centring
## helps when the likelihood is weak and hurts when it is strong, and with
## 400,000 training rows the likelihood here is about as informative as it
## gets.  If that is right, `ncp = 0` should be dramatically faster and the
## tree depth should drop.
##
## Usage:
##   Rscript tests/test_case_study_timing.R [--data=PATH] [--iter=200]
##                                          [--envs=0] [--target=Madrid]
##
## `--envs=N` keeps only the N largest training environments, which isolates
## the cost of E from the cost of N.

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
data_path <- parse_flag(cli, "data",
                        file.path(dirname(root), "old_code", "quiron",
                                  "quiron_final.csv"))
n_iter <- as.integer(parse_flag(cli, "iter", "200"))
n_envs <- as.integer(parse_flag(cli, "envs", "0"))
target <- parse_flag(cli, "target", "Madrid")

## Reuse the loader from the case-study driver without executing the driver.
src <- readLines(file.path(root, "scripts", "05_case_study.R"))
start <- grep("^load_case_data <- function", src)
ends <- grep("^\\}$", src)
eval(parse(text = paste(src[start:min(ends[ends > start])], collapse = "\n")))

dat <- load_case_data(data_path)
model <- load_bgi_model(file.path(root, "stan", "gi_hd.stan"),
                        cache_dir = file.path(root, "results", "compiled"))

is_target <- dat$z == target
if (!any(is_target)) {
  stop("Target environment '", target, "' not present.", call. = FALSE)
}
keep <- !is_target

if (n_envs > 0L) {
  sizes <- sort(table(dat$z[keep]), decreasing = TRUE)
  chosen <- names(sizes)[seq_len(min(n_envs, length(sizes)))]
  keep <- keep & dat$z %in% chosen
}

x_tr <- dat$x[keep, , drop = FALSE]
y_tr <- dat$y[keep]
z_tr <- dat$z[keep]
x0 <- dat$x[head(which(is_target), 300), , drop = FALSE]

cat(sprintf("\nTarget %s | N_train %d | E_train %d | p %d | iter %d, 1 chain\n\n",
            target, nrow(x_tr), length(unique(z_tr)), ncol(x_tr), n_iter))

## Time the data preparation separately: it loops over environments and
## touches every row, so it is a plausible suspect independent of the sampler.
t_prep <- system.time(
  prep <- prepare_bgi_data(x_tr, y_tr, z_tr, x0, cov_method = "pooled")
)[["elapsed"]]
cat(sprintf("prepare_bgi_data: %.1f s\n\n", t_prep))

results <- list()
for (ncp in c(1, 0)) {
  t0 <- Sys.time()
  fit <- tryCatch(
    fit_bgi(x_tr, y_tr, z_tr, x0, model = model, ncp = ncp,
            chains = 1, iter = n_iter, warmup = floor(n_iter / 2), seed = 1),
    error = function(e) {
      cat("  ncp = ", ncp, " failed: ", conditionMessage(e), "\n", sep = "")
      NULL
    })
  if (is.null(fit)) next
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sp <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)[[1]]

  cat(sprintf(
    "ncp = %d | %7.1f s | mean treedepth %4.1f | max-depth hits %3d/%d | div %3d | Rhat %.3f | min ESS %5.0f\n",
    ncp, elapsed, mean(sp[, "treedepth__"]),
    sum(sp[, "treedepth__"] >= 10), nrow(sp),
    fit$diagnostics$n_divergent, fit$diagnostics$max_rhat,
    fit$diagnostics$min_ess_bulk))
  results[[as.character(ncp)]] <- list(seconds = elapsed, fit = fit)
}

if (length(results) == 2L) {
  speedup <- results[["1"]]$seconds / results[["0"]]$seconds
  cat(sprintf("\nCentred is %.1fx %s than non-centred.\n", max(speedup, 1 / speedup),
              if (speedup > 1) "faster" else "slower"))
  g1 <- colMeans(results[["1"]]$fit$draws$gamma)
  g0 <- colMeans(results[["0"]]$fit$draws$gamma)
  cat(sprintf("max |gamma difference| between parameterisations: %.4f\n",
              max(abs(g1 - g0))))
  cat("(they target the same posterior, so this should be Monte Carlo noise)\n")
}

best <- min(vapply(results, function(r) r$seconds, numeric(1)))
cat(sprintf("\nExtrapolated: one chain of 2000 iterations ~ %.1f min;\n", 
            best * (2000 / n_iter) / 60))
cat(sprintf("4 chains in parallel on 4 cores ~ the same; 51 folds ~ %.1f h\n",
            best * (2000 / n_iter) / 60 * 51 / 60))
cat("\nNOTE: warmup here is n_iter/2. Stan's adaptation needs roughly 150-200\n")
cat("warmup iterations to complete its windows, so short runs cannot diagnose\n")
cat("convergence -- only per-iteration cost. Judge Rhat from a full-length run.\n")
