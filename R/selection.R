## selection.R ------------------------------------------------------------
##
## Decision rules for causal parenthood, and the error criteria they control.
##
## Inferential target.  Under the structural model of Section 2 the estimand is
##
##     pa(Y) = { j in 1:p : gamma*_j != 0 },
##
## the set of covariates that appear with a nonzero coefficient in the
## structural assignment for Y.  This is a statement about gamma*, and it is
## identified -- jointly with K* -- as soon as the environment means span R^p.
## It is *not* a statement about conditional independence or about invariance
## of Y | X_S across environments, which is what the invariance-based
## literature targets and which hidden confounding destroys.
##
## ---- The `sign` rule as a Bayes rule ---------------------------------------
##
## The rule of Section 2.1 is not a test and alpha is not a level.  Both facts
## follow from writing the decision problem out, which is worth doing because
## the derivation fixes what alpha means.
##
## The decision for each coordinate has *three* actions, not two: declare the
## effect positive, declare it negative, or decline to call it.  Two actions
## would not do -- under 0-1 sign loss with only {+, -} the Bayes rule is
## always "declare the sign with the larger posterior probability", which never
## abstains, selects every coordinate, and produces no threshold at all.
##
## With action space {+, -, empty} and losses
##
##     L(+, gamma)     = 1{gamma < 0}
##     L(-, gamma)     = 1{gamma > 0}
##     L(empty, gamma) = lambda                (constant, 0 < lambda < 1/2)
##
## the posterior expected losses are P(gamma < 0 | D), P(gamma > 0 | D) and
## lambda.  Declaring the majority sign costs exactly
##
##     lfsr_j = min{ P(gamma_j > 0 | D), P(gamma_j < 0 | D) },
##
## so the Bayes rule is: declare sgn(gamma_j) when lfsr_j < lambda, otherwise
## abstain.  That is the rule implemented below, with alpha = lambda.
##
## The consequence worth stating in the manuscript: **alpha is a loss ratio,
## not an error rate.**  It is the cost of declining to call a coordinate
## relative to the cost of getting its sign wrong.  Setting alpha = 0.05
## asserts that a sign error is twenty times as costly as an abstention.  That
## is a preference, declared openly, and it is exactly why the quantity being
## thresholded is not a p-value.  (lambda >= 1/2 would never abstain, since
## lfsr <= 1/2 always.)
##
## What the rules control.  None of the rules below is a frequentist test.
## Each thresholds a posterior probability:
##
##   `sign`  selects j when the local false sign rate
##              lfsr_j = min{ P(gamma_j > 0 | D), P(gamma_j < 0 | D) }
##           falls below alpha.  Bayes rule for the loss above; it bounds the
##           *posterior* probability of a sign error at alpha for each selected
##           coordinate.  This is the rule stated in Section 2.1.
##
##   `ci`    selects j when 0 lies outside the central (1 - alpha) credible
##           interval for gamma_j.  For a central interval this is exactly
##           `sign` at level alpha / 2, so the two rules are the same family
##           read at different levels; reporting both makes the relationship
##           explicit rather than leaving it implicit.
##
##   `bayes_fdr` orders the coordinates by lfsr and selects the largest
##           initial segment whose average lfsr does not exceed q.  Since the
##           lfsr are posterior probabilities of a *sign* error, their average
##           over the selected set is the posterior expected proportion of
##           sign errors (the direct posterior probability approach of Newton
##           et al., 2004).  Neither `sign` nor `ci` bounds any proportion --
##           they act coordinate by coordinate.
##
##   `rope`  selects j when the posterior mass outside (-delta, delta) exceeds
##           1 - alpha.  See below for why this rule exists.
##
## ---- What the lfsr rules do *not* control ---------------------------------
##
## An important limitation, and one the simulations make concrete.  Under a
## continuous prior on gamma, P(gamma_j = 0 | D) = 0 for every j, so no rule
## built from the posterior of gamma alone can bound the probability of
## selecting a coordinate whose true value is *exactly* zero.  The lfsr bounds
## the probability of attributing the wrong *sign*, which is a different
## quantity: for gamma_j = 0 exactly, no sign is correct, and the lfsr of a
## true null is not required to be near 1/2.
##
## The consequence is visible in `scripts/01_sim_support_recovery.R`: the
## `bayes_fdr` rule holds its posterior expected sign-error rate at or below q
## by construction while the realised proportion of exactly-zero coefficients
## among the selected set runs far higher.  Both numbers are correct; they
## measure different things.
##
## Three coherent ways out, and the manuscript should pick one explicitly:
##
##   * Keep the sign target and say so.  Declare directional confidence to be
##     the estimand, drop the claim to be testing membership of pa(Y), and
##     report the sign-error rate as the headline quantity.  This is defensible
##     -- see `sign_error` in the simulation output, which is 0.000 in three of
##     four scenarios -- but note that it is a genuine change of estimand, not
##     a rephrasing.  In a structural causal model a non-parent has
##     gamma*_j = 0 *exactly*, by construction of the assignment, so the usual
##     "no effect is ever exactly zero" argument is weaker here than in, say,
##     genomics: the zeros are structural rather than idealised.  Whichever
##     estimand is chosen, the false-parent rate still has to be reported,
##     because it is the first thing a reader will want to know.
##   * Restate the target as practical significance.  If the estimand is
##     "gamma_j is not negligible" rather than "gamma_j is exactly zero", use
##     `rope` with a delta set on the scale of the application.  Well defined
##     under a continuous prior.
##   * Put an atom at zero.  A spike-and-slab or similar selection prior makes
##     P(gamma_j = 0 | D) positive and posterior inclusion probabilities
##     meaningful.  Section 2.1 declines this for stated reasons; the cost of
##     declining it is that FDR over exact nulls is not available.
##
## How the sign target actually behaves, from `02_aggregate_support_recovery.R`
## ("Sign-error calibration"), which reports the posterior expected sign-error
## rate over the selected set against the realised one:
##
##   scenario           claimed   realised   usable fits
##   conf0_strong       0.0032     0.0000     10 / 10
##   conf2_strong       0.0048     0.0000     10 / 10
##   conf0_weak         0.0065     0.0000      6 / 10
##   conf2_weak         0.0180     0.3333      3 / 10
##
## In every scenario with enough usable fits to judge, the `sign` rule is
## *conservative* on its own criterion: it never got a direction wrong.  That
## is the empirical case for making direction the estimand.
##
## A caution that survives the reframing.  Nothing here establishes that the
## sign target is calibrated where identification is weak: `conf2_weak` retains
## only 3 of 10 fits, and its realised rate carries a Monte Carlo standard
## error as large as the estimate, so it supports no conclusion either way.
## The point is that reframing the estimand cannot *by itself* buy calibration
## -- an overconfident posterior understates the lfsr just as it understates
## anything else -- so the weak-identification regime needs settling on its own
## terms.  See `scripts/04_sim_environment_budget.R`.
##
## Frequentist behaviour of every rule here is measured, not assumed.

#' Local false sign rate of each posterior coordinate.
#'
#' @param draws A `draws x p` matrix of posterior samples.
#' @return A numeric vector of length `p`.
local_false_sign_rate <- function(draws) {
  draws <- as.matrix(draws)
  p_pos <- colMeans(draws > 0)
  p_neg <- colMeans(draws < 0)
  pmin(p_pos, p_neg)
}

#' Select causal parents from posterior draws of `gamma`.
#'
#' @param gamma_draws A `draws x p` matrix of posterior samples of `gamma`.
#' @param alpha Level for the `sign`, `ci` and `rope` rules.
#' @param q Target posterior expected sign-error rate for `bayes_fdr`.
#' @param rule One of "sign", "ci", "bayes_fdr", "rope".
#' @param delta Half-width of the region of practical equivalence, used only
#'   by `rule = "rope"`.  Interpreted on the scale of `gamma`, so it should be
#'   set from what counts as a negligible effect in the application.
#' @return A list with the selected index set, the per-coordinate lfsr, and
#'   the posterior expected sign-error rate over the selection.
select_parents <- function(gamma_draws, alpha = 0.05, q = 0.1,
                           rule = c("sign", "ci", "bayes_fdr", "rope"),
                           delta = 0.1) {
  rule <- match.arg(rule)
  gamma_draws <- as.matrix(gamma_draws)
  p <- ncol(gamma_draws)
  lfsr <- local_false_sign_rate(gamma_draws)

  selected <- switch(
    rule,
    sign = which(lfsr < alpha),
    rope = {
      ## Posterior probability that gamma_j is practically zero.  Unlike the
      ## lfsr this is a quantity a continuous posterior can bound, because it
      ## refers to an interval rather than to a point.
      p_rope <- colMeans(abs(gamma_draws) <= delta)
      which(p_rope < alpha)
    },
    ci = {
      lower <- apply(gamma_draws, 2, stats::quantile, probs = alpha / 2)
      upper <- apply(gamma_draws, 2, stats::quantile, probs = 1 - alpha / 2)
      which(lower > 0 | upper < 0)
    },
    bayes_fdr = {
      ord <- order(lfsr)
      running <- cumsum(lfsr[ord]) / seq_len(p)
      k <- which(running <= q)
      if (length(k) == 0L) integer(0) else sort(ord[seq_len(max(k))])
    }
  )

  list(
    selected = as.integer(selected),
    lfsr = lfsr,
    ## Posterior expected proportion of *sign* errors over the selection.
    ## Not the posterior probability of having selected an exact zero -- see
    ## the note at the top of this file.
    posterior_esr = if (length(selected) == 0L) 0 else mean(lfsr[selected])
  )
}

#' Benjamini-Hochberg step-up procedure.
#'
#' Used to give the p-value based baselines a false-discovery-rate rule
#' comparable to `bayes_fdr`, so that the comparison is not between a
#' per-coordinate rule and a set-level rule.
#'
#' @param p_values Numeric vector of p-values.
#' @param q Target FDR.
#' @return Integer vector of selected indices.
select_bh <- function(p_values, q = 0.1) {
  m <- length(p_values)
  ord <- order(p_values)
  thresh <- q * seq_len(m) / m
  below <- which(p_values[ord] <= thresh)
  if (length(below) == 0L) {
    return(integer(0))
  }
  sort(ord[seq_len(max(below))])
}
