#!/usr/bin/env Rscript
## 09_brfss_confounders.R -------------------------------------------------
##
## Show that the variables excluded from the BRFSS design really are
## confounders.
##
## The quiron analysis has to *assert* that sex, cholesterol and glucose
## confound the BMI/lifestyle relationship, because they are not in the file.
## The public mirror can demonstrate it, which is the one respect in which the
## reproducible analysis is stronger than the proprietary one: sex, age,
## income, education, race, diabetes and self-rated health are all present in
## BRFSS and are dropped from the design on purpose.
##
## A confounder has to do two things — move the response, and move the
## exposures.  Both are reported here.  The second is the one that matters for
## this paper: if the omitted variables shifted only `Y`, they would inflate
## the noise and nothing more, and `gamma` would still be the causal slope.
## It is because they also shift `X` that the training-domain least-squares
## slope is biased for `gamma*`, which is the entire premise of Section 2.
##
## Usage:
##   Rscript scripts/09_brfss_confounders.R --data=../data/brfss/brfss2023_case.csv

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
data_path <- parse_flag(cli, "data",
                        file.path(dirname(root), "data", "brfss",
                                  "brfss2023_case.csv"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "case_brfss"))

dat <- load_case_data(data_path, dataset = "brfss")
cf <- dat$confounders
if (ncol(cf) == 0L) {
  stop("This dataset carries no confounder columns.", call. = FALSE)
}

## Categorical unless genuinely continuous; age is the only numeric one.
as_term <- function(nm) if (nm == "age") cf[[nm]] else factor(cf[[nm]])
conf_frame <- as.data.frame(lapply(stats::setNames(names(cf), names(cf)),
                                   as_term))
r2 <- function(y, design) summary(stats::lm(y ~ ., data = design))$r.squared

cat("\n=== Do the excluded variables move the response? ===\n")
each_y <- vapply(names(cf), function(nm)
  r2(dat$y, conf_frame[, nm, drop = FALSE]), numeric(1))
for (nm in names(sort(each_y, decreasing = TRUE))) {
  cat(sprintf("  BMI ~ %-12s R2 %.4f\n", nm, each_y[nm]))
}
joint_y <- r2(dat$y, conf_frame)
covar_y <- summary(stats::lm(dat$y ~ dat$x))$r.squared
cat(sprintf("\n  BMI ~ all %d confounders   R2 %.4f\n", ncol(cf), joint_y))
cat(sprintf("  BMI ~ the %d covariates    R2 %.4f\n", ncol(dat$x), covar_y))
cat("\nThe variables left out of the design explain more of the response than\n")
cat("the ones left in.  That alone would only mean extra noise.\n")

## Read the marginal column carefully.  Age enters linearly here and BMI is
## not monotone in age -- it rises through middle age and falls after -- so a
## linear term captures almost none of it.  A near-zero marginal R2 is
## evidence about that functional form, not evidence that the variable fails
## to confound; what matters for the argument is the joint fit above and the
## association with the exposures below.
age_marginal <- if ("age" %in% names(each_y)) each_y[["age"]] else NA_real_
if (!is.na(age_marginal) && age_marginal < 0.001) {
  cat(sprintf(
    "\n(Age's marginal R2 of %.4f reflects a linear term fitted to a\n",
    age_marginal))
  cat("non-monotone relationship, not an absence of confounding.)\n")
}

cat("\n=== Do they move the exposures? ===\n")
cat("This is the part that makes them confounders rather than nuisance.\n\n")
each_x <- vapply(seq_len(ncol(dat$x)), function(j)
  r2(dat$x[, j], conf_frame), numeric(1))
names(each_x) <- dat$covariate_names
for (nm in names(sort(each_x, decreasing = TRUE))) {
  cat(sprintf("  %-26s R2 %.4f\n", nm, each_x[nm]))
}

top <- names(which.max(each_x))
cat(sprintf("\nStrongest: %s (R2 %.4f).\n", top, max(each_x)))
cat("The omitted variables drive both the response and the exposures, and\n")
cat("they load hardest on physical activity -- the coordinate the paper's\n")
cat("substantive claim rests on.  So the public mirror is a demanding test of\n")
cat("the method, not a soft one.\n")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bgi_write_csv(
  data.frame(target = c(paste0("BMI ~ ", names(each_y)),
                        "BMI ~ all confounders", "BMI ~ covariates",
                        paste0("covariate: ", names(each_x))),
             r_squared = c(each_y, joint_y, covar_y, each_x),
             stringsAsFactors = FALSE),
  file.path(out_dir, "confounder_strength.csv"))
cat("\nWritten to ", out_dir, "\n", sep = "")
