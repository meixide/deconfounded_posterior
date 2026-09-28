#!/usr/bin/env Rscript
## test_baseline_identities.R ---------------------------------------------
##
## The distribution-shift baselines are not independent methods; they sit on a
## single path, and the identities between them are worth locking in because
## they are the cleanest way to say where GI sits in the literature.
##
##   anchor(gamma_anchor = 1)         == pooled OLS
##   anchor(gamma_anchor -> infinity) == 2SLS with environment instruments
##   pooled GI                        == 2SLS with environment instruments
##
## The third is the substantive one.  It says that the baseline Referee 2
## proposes ("pooled OLS on X and the environment means") and the IV comparison
## they ask to see beyond the single-source example are the same estimator, and
## that the frequentist GI estimate of `gamma` is IV with the environment as
## instrument.  Proof by Frisch-Waugh-Lovell is in R/baselines.R.
##
## These are exact algebraic identities, so the tolerances are numerical, not
## statistical.
##
## Usage:  Rscript tests/test_baseline_identities.R

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "..", "scripts", "_bootstrap.R"))
})
root <- bgi_bootstrap()

failures <- character(0)
check <- function(name, ok, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!isTRUE(ok)) failures <<- c(failures, name)
}

## Exactly identified (E = p + 1) and over-identified (E > p + 1), plus a case
## with unequal environment sizes, since that is where a projection written
## carelessly would break.
configs <- list(
  list(p = 4, n_env = 5,  n_e = 250, tag = "exactly identified"),
  list(p = 4, n_env = 12, n_e = 150, tag = "over-identified"),
  list(p = 6, n_env = 7,  n_e = 200, tag = "exactly identified, p = 6"),
  list(p = 3, n_env = 9,  n_e = 120, tag = "over-identified, p = 3")
)

for (cfg in configs) {
  cat("\n", cfg$tag, sprintf(" (p = %d, E = %d)\n", cfg$p, cfg$n_env), sep = "")
  bgi_set_seed(31337, cfg$n_env * 10 + cfg$p)
  dat <- simulate_gi_data(n_e = cfg$n_e, p = cfg$p, s0 = 2,
                          n_env = cfg$n_env, q = 3, confounding = 1.5,
                          n0 = 100)

  ols <- fit_ols(dat$x, dat$y)
  iv <- fit_iv_2sls(dat$x, dat$y, dat$z)
  pgi <- fit_pooled_gi(dat$x, dat$y, dat$z)
  a_one <- fit_anchor(dat$x, dat$y, dat$z, gamma_anchor = 1)
  a_inf <- fit_anchor(dat$x, dat$y, dat$z, gamma_anchor = 1e10)

  scale_g <- max(1, max(abs(ols$gamma)))

  check("anchor(1) == OLS",
        max(abs(a_one$gamma - ols$gamma)) / scale_g < 1e-8,
        sprintf("max rel diff %.2e", max(abs(a_one$gamma - ols$gamma)) / scale_g))

  check("anchor(inf) == 2SLS",
        max(abs(a_inf$gamma - iv$gamma)) / scale_g < 1e-5,
        sprintf("max rel diff %.2e", max(abs(a_inf$gamma - iv$gamma)) / scale_g))

  check("pooled GI == 2SLS",
        max(abs(pgi$gamma - iv$gamma)) / scale_g < 1e-8,
        sprintf("max rel diff %.2e", max(abs(pgi$gamma - iv$gamma)) / scale_g))

  ## The point estimates coincide but the standard errors do not: pooled GI
  ## conditions on the environment means as if they were fixed.
  se_ratio <- max(pgi$se / iv$se)
  check("pooled GI standard errors are anti-conservative",
        se_ratio < 1,
        sprintf("max se ratio pooled_gi / 2SLS = %.3f", se_ratio))
}

## Unequal environment sizes: the projection must weight by environment
## membership, not assume a balanced design.
cat("\nunequal environment sizes\n")
bgi_set_seed(31337, 999)
dat <- simulate_gi_data(n_e = 200, p = 4, s0 = 2, n_env = 6, q = 3,
                        confounding = 1.5, n0 = 100)
keep <- unlist(lapply(sort(unique(dat$z)), function(e) {
  idx <- which(dat$z == e)
  idx[seq_len(round(length(idx) * stats::runif(1, 0.4, 1)))]
}))
iv_u <- fit_iv_2sls(dat$x[keep, ], dat$y[keep], dat$z[keep])
pgi_u <- fit_pooled_gi(dat$x[keep, ], dat$y[keep], dat$z[keep])
check("pooled GI == 2SLS with unbalanced environments",
      max(abs(pgi_u$gamma - iv_u$gamma)) / max(1, max(abs(iv_u$gamma))) < 1e-8,
      sprintf("sizes %s", paste(table(dat$z[keep]), collapse = "/")))

## 2SLS must refuse to run rather than return nonsense when under-identified.
cat("\nunder-identification is refused\n")
bgi_set_seed(31337, 4242)
dat_u <- simulate_gi_data(n_e = 200, p = 6, s0 = 2, n_env = 3, q = 3,
                          confounding = 1, n0 = 50)
check("2SLS errors when E - 1 < p",
      inherits(tryCatch(fit_iv_2sls(dat_u$x, dat_u$y, dat_u$z),
                        error = function(e) e), "error"),
      "E = 3, p = 6")

cat("\n")
if (length(failures) == 0L) {
  cat("All baseline identities hold.\n")
} else {
  cat("FAILED: ", paste(unique(failures), collapse = "; "), "\n", sep = "")
  quit(status = 1)
}
