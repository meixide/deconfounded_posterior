## 14_aggregate_contraction.R --------------------------------------------
##
## Aggregate the Tier 1 contraction study of 13_sim_contraction.R.
##
## The comparison is *within* a replication: 13_sim_contraction.R holds the
## design (w*, mu_e, Sigma_e, and therefore c_in) fixed across the seven
## sample sizes of a replication and varies only n_e.  So the log-log slope is
## fitted per replication and then summarised across replications, rather than
## fitted once to the pooled cloud, which would mix the across-design variation
## in c_in back in.
##
## What Theorem 2 predicts.  With r_N = M sqrt(log N / N),
##
##     d log r_N / d log N = -1/2 + 1/(2 log N),
##
## which over N = 200 to 12800 runs from about -0.41 to -0.45.  A contraction
## radius tracking the rate should give a slope in that neighbourhood; a slope
## at or below it means the posterior contracts at least as fast as the rate.
##
## Usage:
##   Rscript scripts/14_aggregate_contraction.R

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
# An aggregator reads CSVs; it neither fits nor simulates, so it does not
# need the fitting dependencies.  See `need_stan` in R/setup.R.
root <- bgi_bootstrap(need_stan = FALSE)

in_dir <- file.path(root, "results", "contraction")
files <- list.files(in_dir, pattern = "^task_\\d+\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0L)
d <- do.call(rbind, lapply(files, read.csv))
d <- d[order(d$rep_id, d$N), ]

cat(sprintf("\n%d rows, %d replications, sizes %s\n",
            nrow(d), length(unique(d$rep_id)),
            paste(sort(unique(d$n_e)), collapse = ", ")))

## The design must be constant within a replication; if it is not, the seeding
## fix in 13_sim_contraction.R has regressed and the slopes below are not
## interpretable.
cin_spread <- tapply(d$c_in, d$rep_id, function(v) diff(range(v)))
cat(sprintf("max within-replication spread of c_in: %.3g %s\n",
            max(cin_spread),
            if (max(cin_spread) < 1e-8) "(constant, as intended)" else
              "*** NOT CONSTANT -- check the seeding ***"))

## Convergence filter.  A fit whose chains did not mix has a posterior spread
## that means nothing, and because the radius enters on a log scale a single
## failure moves the group mean a long way: in the first run of this study one
## fit at N = 12800 came back with R-hat 23.6 and a radius of 16.5 against a
## typical 1.0, which on its own reversed the direction of the last point.
## Drop non-converged fits and say how many, rather than reporting a rate
## computed partly from sampler failure.
rhat_max <- 1.05
bad <- is.na(d$max_rhat) | d$max_rhat > rhat_max
if (any(bad)) {
  cat(sprintf("\ndropping %d of %d fits with R-hat > %.2f or NA:\n",
              sum(bad), nrow(d), rhat_max))
  print(d[bad, c("rep_id", "N", "max_rhat", "radius_q95")],
        digits = 3, row.names = FALSE)
}
d <- d[!bad, ]

cat("\n-- by sample size, averaged over converged replications --\n")
agg <- aggregate(
  cbind(radius_q95, M_star, post_mean_err, mass_U, tail_mass_M1,
        max_rhat, min_ess, seconds) ~ N,
  data = d, FUN = mean)
print(agg, digits = 3, row.names = FALSE)

## Per-replication log-log slope.
slopes <- sapply(split(d, d$rep_id), function(s)
  unname(coef(lm(log(radius_q95) ~ log(N), data = s))[2]))
tt <- t.test(slopes)
cat(sprintf(
  "\nlog-log slope of radius_q95 on N, per replication:\n  mean %.3f  [%.3f, %.3f]  sd %.3f  (n = %d)\n",
  mean(slopes), tt$conf.int[1], tt$conf.int[2], sd(slopes), length(slopes)))
Nr <- range(d$N)
cat(sprintf("  Theorem 2 predicts %.3f to %.3f over N = %d to %d\n",
            -0.5 + 1 / (2 * log(Nr[2])), -0.5 + 1 / (2 * log(Nr[1])),
            Nr[1], Nr[2]))

## Is M_star bounded?  A rate that is too slow shows up as M_star growing.
ms <- sapply(split(d, d$rep_id), function(s)
  unname(coef(lm(log(M_star) ~ log(N), data = s))[2]))
tm <- t.test(ms)
cat(sprintf(
  "\nlog-log slope of M_star on N: mean %.3f [%.3f, %.3f]\n  %s\n",
  mean(ms), tm$conf.int[1], tm$conf.int[2],
  if (tm$conf.int[2] < 0)
    "M_star declines: contraction is at least as fast as the stated rate."
  else if (tm$conf.int[1] > 0)
    "*** M_star grows: contraction is SLOWER than the stated rate. ***"
  else "M_star is flat: contraction matches the stated rate."))

out <- file.path(root, "results", "summaries", "contraction.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(agg, out, row.names = FALSE)
cat(sprintf("\nWritten to %s\n", out))
