## _bootstrap.R -----------------------------------------------------------
##
## Locate the project root and load the library.  Every script in scripts/
## starts with:
##
##   source(file.path(dirname(sys.frame(1)$ofile %||% "."), "_bootstrap.R"))
##
## which is fragile, so instead each script inlines the two lines below.  This
## file exists so the logic is written once and can be sourced interactively.

bgi_bootstrap <- function(need_stan = TRUE) {
  from_env <- Sys.getenv("BGI_ROOT", unset = "")
  root <- if (nzchar(from_env)) {
    normalizePath(from_env)
  } else {
    args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(args) > 0L) {
      normalizePath(file.path(dirname(sub("^--file=", "", args[1])), ".."))
    } else {
      normalizePath(".")
    }
  }
  if (!dir.exists(file.path(root, "R"))) {
    stop("Project root ", root, " does not contain R/. ",
         "Set BGI_ROOT to the new_code directory.", call. = FALSE)
  }
  source(file.path(root, "R", "setup.R"))
  get("bgi_setup")(root, need_stan = need_stan)
  root
}
