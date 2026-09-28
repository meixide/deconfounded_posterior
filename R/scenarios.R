## scenarios.R ------------------------------------------------------------
##
## The scenario grid for the support-recovery study, and the mapping from a
## flat SLURM array index to a (scenario, replication) pair.
##
## Keeping the grid in one function means the driver, the aggregator and the
## SLURM submission script all agree on how many tasks exist and what each
## one does, without any of them hard-coding a count.

#' Scenario grid for the support-recovery study.
#'
#' The grid crosses three factors that the reviewers' concern turns on:
#'
#'   confounding      strength of the hidden confounder's effect on Y.  At 0
#'                    there is no confounding and every method should behave;
#'                    the interesting regime is where OLS is badly biased.
#'   identifiability  geometry of the environment means.  "strong" satisfies
#'                    the spanning condition comfortably; "weak" sits near the
#'                    boundary discussed in Section 2.3, where the prior takes
#'                    over.  Support recovery under "weak" is what tells us
#'                    whether near-violations of identifiability produce false
#'                    discoveries or merely a loss of power.
#'   n_e              observations per environment, which controls both the
#'                    posterior contraction rate and the conditioning of the
#'                    plugged-in Sigma_e.
#'
#' @param small Reduce the grid to a quick smoke-test configuration.
#' @return A data frame with one row per scenario.
sim_scenarios <- function(small = FALSE) {
  if (small) {
    grid <- expand.grid(
      confounding = c(0, 2),
      identifiability = c("strong", "weak"),
      n_e = 200L,
      p = 6L,
      s0 = 3L,
      stringsAsFactors = FALSE
    )
  } else {
    grid <- expand.grid(
      confounding = c(0, 1, 2),
      identifiability = c("strong", "weak"),
      n_e = c(100L, 400L),
      p = 6L,
      s0 = 3L,
      stringsAsFactors = FALSE
    )
  }
  ## E = p + 1 is the minimum for the environment means to span R^p.  Running
  ## at the minimum makes the between-environment information about K as thin
  ## as it can be while still being identified, which is the regime the
  ## identifiability discussion is about.
  grid$n_env <- grid$p + 1L
  grid$q <- 3L
  grid$n0 <- 500L
  grid$gamma_signal <- 1
  grid$scenario_id <- seq_len(nrow(grid))
  grid$label <- sprintf("conf%s_%s_n%d",
                        format(grid$confounding, trim = TRUE),
                        grid$identifiability, grid$n_e)
  grid[, c("scenario_id", "label", "confounding", "identifiability",
           "n_e", "n0", "p", "s0", "n_env", "q", "gamma_signal")]
}

#' Full task list: every (scenario, replication) pair.
#'
#' @param n_rep Replications per scenario.
#' @param small Passed to `sim_scenarios()`.
#' @return A data frame with one row per SLURM array task.
sim_task_table <- function(n_rep = 20L, small = FALSE) {
  scen <- sim_scenarios(small = small)
  tasks <- scen[rep(seq_len(nrow(scen)), each = n_rep), , drop = FALSE]
  tasks$rep_id <- rep(seq_len(n_rep), times = nrow(scen))
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  tasks
}
