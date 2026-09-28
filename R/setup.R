## setup.R ----------------------------------------------------------------
##
## Single entry point for every script: resolves paths, loads the library,
## checks dependencies and fixes global options.  Source this and nothing
## else.
##
## Dependencies are deliberately minimal.  Only `rstan` and `mvtnorm` are
## required; `ggplot2` is used by the plotting script alone and its absence is
## reported rather than fatal.  Everything else is base R.  Code follows
## tidyverse naming and layout conventions (snake_case, one idea per
## function, explicit argument names) without taking a dependency on the
## tidyverse itself, so that the whole pipeline runs in a clean R
## installation on a cluster node with no network access.

BGI_REQUIRED_PACKAGES <- c("rstan", "mvtnorm")
BGI_SUGGESTED_PACKAGES <- c("ggplot2")

#' Locate the project root from any working directory.
#'
#' Looks for the directory that contains both `R/` and `stan/`, starting from
#' the `BGI_ROOT` environment variable, then the current directory, then its
#' ancestors.
bgi_root <- function() {
  from_env <- Sys.getenv("BGI_ROOT", unset = "")
  candidates <- c(if (nzchar(from_env)) from_env, ".", "..", "../..")
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R")) &&
        dir.exists(file.path(cand, "stan"))) {
      return(normalizePath(cand))
    }
  }
  stop("Cannot locate the project root. Set BGI_ROOT to the new_code ",
       "directory, or run from inside it.", call. = FALSE)
}

#' Check that the required packages are installed, with an actionable message.
bgi_check_dependencies <- function(required = BGI_REQUIRED_PACKAGES,
                                   suggested = BGI_SUGGESTED_PACKAGES) {
  missing <- required[!vapply(required, requireNamespace, logical(1),
                              quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c(",
         paste0('"', missing, '"', collapse = ", "), "))", call. = FALSE)
  }
  absent <- suggested[!vapply(suggested, requireNamespace, logical(1),
                              quietly = TRUE)]
  if (length(absent) > 0L) {
    message("Optional packages not installed (plotting will be skipped): ",
            paste(absent, collapse = ", "))
  }
  invisible(TRUE)
}

#' Load the project library and set global options.
#'
#' @param root Project root; defaults to `bgi_root()`.
#' @param quiet Suppress the startup banner.
#' @return The project root, invisibly.
#' @param need_stan Whether this script will fit a model or simulate data.
#'   The aggregation scripts do neither: they read the CSVs the fitting scripts
#'   wrote.  Requiring the fitting dependencies of them would stop a reader
#'   from looking at results on a machine without them, for no reason -- and it
#'   did, on a cluster login shell with no R module loaded, where the
#'   aggregator refused to print a table it had every number for.  `mvtnorm`
#'   goes with `rstan` here because its only uses are drawing simulated
#'   covariates and one sampling step inside the fit.
bgi_setup <- function(root = bgi_root(), quiet = FALSE, need_stan = TRUE) {
  bgi_check_dependencies(
    required = if (need_stan) BGI_REQUIRED_PACKAGES else
      setdiff(BGI_REQUIRED_PACKAGES, c("rstan", "mvtnorm")))

  for (f in c("covariance.R", "simulate.R", "fit_bgi.R", "selection.R",
              "baselines.R", "metrics.R", "scenarios.R", "loeo.R",
              "case_data.R")) {
    source(file.path(root, "R", f), local = FALSE)
  }

  if (need_stan) {
    suppressPackageStartupMessages(library(rstan))
    ## Each replication is one process; within-process parallelism is handled
    ## explicitly by the caller, so the defaults are set to be conservative.
    rstan::rstan_options(auto_write = FALSE)
    options(mc.cores = 1L)
  }

  if (!quiet) {
    message("BGI project root: ", root)
    message(if (need_stan)
              paste0("rstan ", utils::packageVersion("rstan"), " | R ",
                     getRversion())
            else paste0("R ", getRversion(), " (no Stan needed here)"))
  }
  invisible(root)
}

#' Deterministic per-replication seed.
#'
#' Uses L'Ecuyer-CMRG streams so that replications are reproducible whether
#' they run sequentially, in forked workers or in independent SLURM array
#' tasks.  The stream is a pure function of `(base_seed, index)`, so a
#' replication can be re-run in isolation and will reproduce exactly.
#'
#' @param base_seed Integer seed for the whole study.
#' @param index Replication index (1-based).
bgi_set_seed <- function(base_seed, index) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  seed <- .Random.seed
  for (i in seq_len(index)) {
    seed <- parallel::nextRNGStream(seed)
  }
  assign(".Random.seed", seed, envir = globalenv())
  invisible(seed)
}

#' Write a data frame as CSV, creating the directory if needed.
#'
#' @param x Data frame.
#' @param path Destination path.
bgi_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}
