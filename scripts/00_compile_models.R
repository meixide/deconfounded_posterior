#!/usr/bin/env Rscript
## 00_compile_models.R ----------------------------------------------------
##
## Compile the Stan models once and cache them under results/compiled/.
##
## This must run to completion before any simulation is launched.  Compiling
## inside a parallel worker makes every worker pay the cost, and concurrent
## workers race on rstan's on-disk cache, which is a common cause of
## otherwise inexplicable failures in array jobs.
##
## Two models are needed for every result in the manuscript and two are used
## only by the robustness checks in tests/.  A failure in the second group
## must not stop the first from being built, because that would leave a fresh
## checkout with no usable model at all: the loop therefore records failures
## and reports them at the end rather than aborting on the first one.
##
## Only the models the manuscript depends on are built by default.  There is one
## optional model, `gi_hd_fullcov.stan`, which two standalone tests use and no
## table does; it needs a newer StanHeaders than rstan 2.32.3 ships with, so on
## a stock installation it fails.  Reporting that failure to someone who only
## wants to reproduce the paper is noise that reads like a broken package, so it
## is built only when asked for.
##
## Usage:
##   Rscript scripts/00_compile_models.R [--force] [--with-optional]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

cli <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% cli
## Required-only is the default; --required-only is still accepted so that any
## existing invocation keeps working.
required_only <- !("--with-optional" %in% cli)

## ---- Toolchain check ---------------------------------------------------
## rstan and StanHeaders are versioned together, and an installation that
## pairs a recent rstan with a stale StanHeaders compiles most models happily
## and then fails on whichever one uses a newer Stan function.  Saying so here
## is much cheaper than letting the reviewer read a `stanc` type error.

rstan_v <- packageVersion("rstan")
headers_v <- packageVersion("StanHeaders")
message(sprintf("rstan %s | StanHeaders %s", rstan_v, headers_v))

toolchain_warning <- NULL
if (package_version(headers_v)$major != package_version(rstan_v)$major ||
    package_version(headers_v)$minor != package_version(rstan_v)$minor) {
  toolchain_warning <- sprintf(paste0(
    "rstan %s is paired with StanHeaders %s.  These are released together ",
    "and are expected to match on major.minor.  A stale StanHeaders is the ",
    "usual reason a model using a recent Stan function fails to compile.\n",
    "  Fix with:  install.packages(c(\"StanHeaders\", \"rstan\"), ",
    "repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\")))"),
    rstan_v, headers_v)
  message("\nNOTE: ", toolchain_warning, "\n")
}

## ---- Models ------------------------------------------------------------
## `required` marks the models without which nothing in the manuscript can be
## reproduced.  `gi_hd_fullcov` gives the environment covariances an
## inverse-Wishart prior instead of plugging them in; it needs Stan >= 2.30
## for `inv_wishart_cholesky` and is used only by tests/test_plugin_vs_fullcov.R
## and tests/test_gamma_coverage_sources.R.

models <- list(
  list(file = "gi_hd.stan",           required = TRUE,
       note = "main model, every table"),
  list(file = "gi_hd_slope.stan",     required = TRUE,
       note = "slope parameterisation, Sections 2.4 and 3.2"),
  list(file = "gi_hd_reference.stan", required = TRUE,
       note = "reference implementation, tests/test_fast_vs_reference.R"),
  list(file = "gi_hd_fullcov.stan",   required = FALSE,
       note = "covariances inferred; needs Stan >= 2.30")
)

if (required_only) {
  models <- Filter(function(m) m$required, models)
}

failures <- list()

for (m in models) {
  t0 <- Sys.time()
  ok <- tryCatch({
    compile_bgi_model(file.path(root, "stan", m$file),
                      cache_dir = file.path(root, "results", "compiled"),
                      force = force)
    TRUE
  }, error = function(e) {
    failures[[m$file]] <<- conditionMessage(e)
    FALSE
  })

  if (ok) {
    message(sprintf("  %-24s ready in %.1f s", m$file,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  } else {
    message(sprintf("  %-24s FAILED (%s)", m$file,
                    if (m$required) "required" else "optional"))
  }
}

## ---- Report ------------------------------------------------------------
## Progress above goes through message(), hence to stderr, which is R's
## convention.  A batch scheduler usually sends stderr and stdout to separate
## files, so the summary is repeated on stdout: otherwise the .out file of a
## compile job is empty and says nothing about whether it worked.

say <- function(...) cat(..., "\n", sep = "")

say("")
for (m in models) {
  say(sprintf("  %-24s %s", m$file,
              if (is.null(failures[[m$file]])) "ok" else
                if (m$required) "FAILED (required)" else "failed (optional)"))
}
say("")
say("Compiled models are in ", file.path(root, "results", "compiled"))

if (length(failures) == 0L) {
  say("RESULT: every model compiled.")
}

if (length(failures) > 0L) {
  required_failed <- vapply(models, function(m) {
    m$required && !is.null(failures[[m$file]])
  }, logical(1))

  message("\n", strrep("-", 70))
  for (f in names(failures)) {
    message("\n", f, " did not compile:\n",
            paste0("  ", strsplit(failures[[f]], "\n")[[1]], collapse = "\n"))
  }
  if (!is.null(toolchain_warning)) {
    message("\n", toolchain_warning)
  }
  message(strrep("-", 70))

  if (any(required_failed)) {
    say("RESULT: a required model failed to compile.")
    stop("A required model failed to compile; nothing downstream will run.")
  }
  say("RESULT: every required model compiled; one optional model did not.")
  message(
    "\nEvery required model compiled.  The failure above is in an optional\n",
    "model, requested with --with-optional and used only by\n",
    "tests/test_plugin_vs_fullcov.R and tests/test_gamma_coverage_sources.R;\n",
    "every table in the manuscript can still be reproduced.")
}
