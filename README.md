# Predictive posteriors under hidden confounding — replication code

Code accompanying the revision of JCGS-25-897, *Predictive posteriors under
hidden confounding*.

Everything here runs from a clean R installation with two packages
(`rstan`, `mvtnorm`); `ggplot2` is optional and used only for figures. There
are no network calls at run time and no absolute paths.

---

## Quick start

Everything in this package that runs on a laptop, in one command:

```bash
cd new_code
bash run_checks.sh
```

About four minutes on an Apple M2 (*measured*, end to end, models already
compiled; add two minutes for the first compile). No cluster, no network
access, no data download. It checks R and the two required packages, compiles the Stan
models, runs the recovery and reference tests, fits one replication of
the smallest cell of Table 1 so that the simulation pipeline is seen working
end to end, and finally recomputes **every table in the paper** from the
per-task CSVs tracked in `results/`, printing the recomputed rows beside the
ones the manuscript prints. It exits non-zero if any step fails.

Three claims are worth keeping apart, because the cost of checking them differs
by four orders of magnitude:

| Claim | Where it is checked | Cost |
|---|---|---|
| The implementation is correct | `run_checks.sh`, steps 1–4 | minutes, laptop |
| No table was transcribed wrongly | `run_checks.sh` step 5, or `bash check_tables.sh` alone | seconds, laptop |
| The numbers are what the model produces from scratch | one job array per table | hours, cluster |

So a referee can confirm on a laptop that every figure in every table is the
figure the recorded fits produce, and needs a cluster only to regenerate those
fits. The numeric results are in the repository for exactly this reason; the
section after next says which directory belongs to which table.

---

## What runs where

The distinction matters for a reader who wants to check the code rather than
re-run the study, so it is stated exactly rather than left to be discovered.

Timings marked *measured* were taken on an Apple M2 laptop, 8 cores, R 4.3.2
with rstan 2.32.3, running one task at a time. Cluster walltimes are the
`#SBATCH --time` ceilings in `slurm/`, which are set for the heaviest task in
each array and are not typical runtimes.

### Runs on a laptop

| What | Command | Cost |
|---|---|---|
| All of the below, in order | `bash run_checks.sh` | **4 min** (*measured*) |
| Compile the Stan models | `Rscript scripts/00_compile_models.R` | 30–60 s per model (*measured*) |
| The four substantive claims, including that predictions carry `S_0` and not the training residual scale | `Rscript tests/test_recovery.R` | ~5 min |
| Sufficient-statistic likelihood equals the per-observation model, checked exactly | `Rscript tests/test_likelihood_identity.R` | ~2 s (*measured*) |
| One replication of the smallest cell of Table 1 | `Rscript scripts/19_sim_dimension_sweep.R --task=1` | 6–12 s (*measured*) |
| **Every table in the paper against `results/`**, rows printed side by side | `bash check_tables.sh` | ~40 s (*measured*) |
| Case-study design summary, no fitting | `Rscript scripts/05_case_study.R --data=... --summary-only` | seconds, but needs the BRFSS download |

### Two checks of the same thing, and why

`stan/gi_hd.stan` reduces the training likelihood to per-environment
sufficient statistics. That is what makes a gradient evaluation independent
of the sample size, and it is also the part of the code least obviously
correct, so it is checked twice.

`tests/test_likelihood_identity.R` is the one `run_checks.sh` runs. The two
models declare identical parameter blocks, so a point in the unconstrained
space means the same thing to both; if the reduction is right, the two
log-posteriors differ at every such point by the same constant. Checking that
across forty random points takes about two seconds and settles the question
outright. On the data used there the difference is constant to 1e-11.

`tests/test_fast_vs_reference.R` fits both models and compares the
posteriors. It additionally exercises sampling and the `generated quantities`
block, which the identity check does not reach, but it takes about
forty-five minutes and answers the question only up to Monte Carlo error: the
reference implementation mixes some forty times more slowly, so its posterior
standard deviations are the noisier estimate and the comparison inherits that
noise. Run it before a release; do not put it in a reviewer's first four
minutes.

### Needs a cluster

Every table in the manuscript is here. The reason is arithmetic, not
convenience: one replication of the support-recovery study is a BGI fit plus
nine baselines and took **54 minutes** (*measured*), and Table 2 is forty of
them, so reproducing that one table serially on a laptop is about 36 hours.

| Result | Submit | Array | Walltime ceiling |
|---|---|---|---|
| **Table 1** — predictive coverage as `p` grows (`table:AEC`) | `sbatch --array=1-48 slurm/19_dimension_sweep.sh --chunk=6` | 48 | see note |
| **Table 2** — support recovery and predictive intervals (`tab:support`) | `sbatch --array=1-40%40 slurm/01_support_recovery.sh --n-rep=10 --small` | 40 | 2 h |
| **Table 3** — simulation-based calibration (`tab:sbc`) | `sbatch slurm/12_sbc.sh` | — | 12 h |
| **Table 4** — BRFSS case study (`tab:brfss`) | `sbatch --array=1-52 slurm/06_case_study.sh brfss` | 52 | 6 h |
| Supplement, contraction rate (`tab:contractionfull`) | `sbatch slurm/13_contraction.sh` | — | 4 h |
| Supplement, environment budget | `sbatch --array=... slurm/04_env_budget.sh` | — | 20 h |
| Supplement, split extrapolation | `sbatch --array=... slurm/16_split_extrapolation.sh` | — | 5.5 h |

A QoS normally caps how many jobs a user may have *submitted*, not merely
running, so a 288-index array can be refused outright even with `%48`
throttling concurrency. On CESGA the `short` QoS allows 100 submitted and 50
running per user, which is why Table 1 is submitted with `--chunk=6`: that
puts six consecutive rows in one index, giving 48 indices that fit under both
caps and all run at once. Raise `--chunk` to shrink the array further, and
raise the job's `--time` with it, since one index then covers proportionally
more replications of the same cell.

Each array task writes its own CSV and skips work already done, so a
pre-empted or failed task can be resubmitted alone. Aggregate afterwards:

```bash
Rscript scripts/20_aggregate_dimension_sweep.R     # Table 1
Rscript scripts/02_aggregate_support_recovery.R    # Table 2
Rscript scripts/12_aggregate_sbc.R                 # Table 3
Rscript scripts/07_aggregate_case_study.R --in=results/case_brfss/folds
Rscript scripts/14_aggregate_contraction.R
```

The full submission recipe, including the module loads and the convergence
gate that should be read before launching the case-study array, is in
[Running on CESGA FinisTerrae III](#running-on-cesga-finisterrae-iii).

### If you only have a laptop and want to see a table form

Reduce the replication count rather than the sampler settings, and read the
result as a smoke test rather than as the manuscript's numbers:

```bash
# Two cells of Table 1, one replication each: minutes, not hours.
Rscript scripts/19_sim_dimension_sweep.R --task=1,25
Rscript scripts/20_aggregate_dimension_sweep.R
```

Cutting `--iter` instead will trip the diagnostic gate in the aggregators,
which excludes replications with divergent transitions or `Rhat > 1.01` and
tells you how many it dropped. That is deliberate: a number that only looks
right after silently discarding bad fits is not worth reporting.

---

## Layout

```
.
├── run_checks.sh               Every laptop-runnable check, one command.
├── check_tables.sh             Recomputes each table in the paper from
│                               results/ and prints the rows beside the
│                               manuscript's. Also lists figures the prose
│                               quotes that appear in no table.
├── stan/
│   ├── gi_hd_slope.stan        The model the paper recommends, and the one
│                               behind Table 2, Table 4 and the case study.
│                               Carries a slope b_e per environment under a
│                               hierarchical prior and derives K_e = Sigma_e b_e,
│                               so no covariance is inverted in the likelihood.
│                               Its header records which assumption this costs.
│   ├── gi_hd.stan              The covariance parameterisation: a common K,
│                               with b_e = Sigma_e^{-1} K derived, as the
│                               manuscript's Eq. (2) writes it. Behind Table 1,
│                               the contraction table and the figures.
│   ├── gi_hd_fullcov.stan      Same model with Sigma_e and Sigma_0 inferred
│                               rather than plugged in.
│   ├── gi_hd_reference.stan    Transparent per-observation model, used only
│                               to validate the fast one.
│   └── unit_variance.stan      One-parameter model, used only by
│                               scripts/22_identifiability_figure.R.
├── R/
│   ├── setup.R                 Paths, dependency checks, RNG streams.
│   ├── simulate.R              Data-generating process.
│   ├── covariance.R            Plug-in covariance estimators.
│   ├── fit_bgi.R               Model compilation, data prep, sampling.
│   ├── selection.R             Causal-parenthood decision rules.
│   ├── baselines.R             OLS, pooled GI, invariant causal prediction.
│   ├── metrics.R               Recovery and calibration summaries.
│   ├── scenarios.R             Scenario grid and SLURM task table.
│   ├── loeo.R                  Leave-one-environment-out evaluation.
│   └── case_data.R             Loading for both case-study datasets.
├── scripts/
│   ├── 00_compile_models.R
│   ├── 01_sim_support_recovery.R      Support recovery / false discoveries.
│   ├── 02_aggregate_support_recovery.R
│   ├── 03_plots.R                     Greyscale, >=10pt figures.
│   ├── 04_sim_environment_budget.R    How many environments are enough?
│   ├── 05_case_study.R                Leave-one-environment-out case study.
│   ├── 06_prepare_brfss.py            BRFSS XPORT -> CSV (no dependencies).
│   ├── 07_aggregate_case_study.R      Combine per-fold array outputs.
│   ├── 19_sim_dimension_sweep.R       Predictive coverage vs p (Table 1).
│   ├── 20_aggregate_dimension_sweep.R Builds Table 1; also reports what
│                                      the training-scale predictive sd
│                                      would have given on the same fits.
│   ├── 21_one_dim_illustration.R       stan.png, posterior_predictive_plot.png
│                                      and firsf.png: the single-predictor
│                                      figures of Section 2 and Supplement E.
│   └── 22_identifiability_figure.R     post_beta.png, the progressive-violation
│                                      figure of the supplement. Standalone,
│                                      and older than the rest of the package:
│                                      it fits stan/unit_variance.stan and goes
│                                      through neither R/setup.R nor the
│                                      sufficient-statistic likelihood.
├── slurm/
│   ├── 00_compile.sh
│   ├── 01_support_recovery.sh
│   ├── 02_aggregate.sh
│   ├── 03_diagnostics.sh
│   ├── 04_env_budget.sh
│   ├── 05_case_convergence.sh
│   ├── 06_case_study.sh
│   ├── 19_dimension_sweep.sh
│   └── submit_support_recovery.sh
├── tests/
│   ├── test_recovery.R                Parameter recovery, predictive variance.
│   ├── test_likelihood_identity.R     Fast likelihood == reference, exactly
│   │                                  (log_prob differs by a constant).
│   ├── test_fast_vs_reference.R       Same question via sampling; ~45 min,
│   │                                  run before a release, not per checkout.
│   ├── test_plugin_vs_fullcov.R       Does inferring Sigma help? (No.)
│   ├── test_gamma_coverage_sources.R  Which modelling choice causes the bias?
│   ├── test_baseline_identities.R     anchor/IV/pooled-GI identities
│   ├── test_case_study_timing.R       Per-iteration cost, ncp = 0 vs 1.
│   └── test_case_study_convergence.R  Rhat at full length. Run before quoting.
└── results/                    Created at run time.
```

Every script resolves the project root itself, or honours `BGI_ROOT`, so
nothing depends on the working directory.

---

## The model

For training environments `e = 1, ..., E`:

```
X_ei         ~ N_p(mu_e, Sigma_e)
Y_ei | X_ei  ~ N(alpha + gamma' X_ei + K' Sigma_e^{-1} (X_ei - mu_e),
                 sigma_y^2 - K' Sigma_e^{-1} K)
```

and in the target domain `e = 0`:

```
Y_0 | X_0    ~ N(alpha + gamma' X_0 + K' Sigma_0^{-1} (X_0 - mu_0),
                 sigma_y^2 - K' Sigma_0^{-1} K)
```

`K = Cov(eps_Y, X)` is a **covariance**, not a slope; the slope on the centred
covariates is `Sigma_e^{-1} K`. Writing it this way makes the training
likelihood and the target predictive formula use the same object, and makes
`sigma_y^2 = Var(eps_Y)` the marginal error variance throughout.

### Three quantities that are easy to confuse

| Symbol | Meaning | In the code |
|---|---|---|
| `sigma_y` | marginal sd of `eps_Y` | `sigma_y` |
| `sigma_cond` | conditional sd in the training domain | `sigma_cond` |
| `S0` | conditional sd in the **target** domain | `S0` |

Posterior predictive draws for `Y_0` use `S0`. This is the point of the
method and it is worth stating plainly: sampling them with the training
residual sd instead gives intervals of the wrong width whenever `Sigma_0`
differs from `Sigma_e`.

### Parameterisation of the variance

The free scale parameter is `v`, the conditional variance referenced to the
pooled within-environment covariance `Sigma_bar`:

```
sigma_y^2 = v + K' Sigma_bar^{-1} K
v_e       = v - (K' Sigma_e^{-1} K - K' Sigma_bar^{-1} K)
S_0^2     = v - (K' Sigma_0^{-1} K - K' Sigma_bar^{-1} K)
```

This is an exact reparameterisation of `(sigma_y^2, K)`, and it is the one
that works. Parameterising by `sigma_y^2` directly requires
`sigma_y^2 > max_j K' Sigma_j^{-1} K`, a constraint on the *level* of the
Mahalanobis norms. Under strong confounding that level sits just below
`sigma_y^2`, so the largest of the `E + 1` noisy plug-in estimates binds the
constraint, and the fit shrinks `K` purely to keep one environment's
conditional variance positive. Because `gamma` is the within-environment slope
minus `Sigma^{-1} K`, that shrinkage lands directly on `gamma` as bias. With
`v` as the free parameter the bound applies to *differences* of Mahalanobis
norms, which are of the order of the covariance estimation error.

In one representative configuration (`p = 4`, `E = 15`, `n_e = 1000`, strong
confounding) this reparameterisation together with pooled-target shrinkage
takes the RMSE of `gamma` from 0.24 to 0.057, against 0.025 for the
frequentist pooled estimator on the same data.

The **causal Mahalanobis condition** `K' Sigma_j^{-1} K < sigma_y^2` is
enforced exactly by this bound. It is not an extra assumption on the
data-generating process — it holds automatically at the true `K` whenever
`(X, Y)` is jointly non-degenerate Gaussian — but it must hold numerically at
every point the sampler visits, which is what the bound guarantees.

### Covariance estimation

`R/covariance.R` provides `sample`, `ledoit_wolf`, `oas`, `ridge` and
`pooled`. The default for the environment covariances is `pooled`:
Ledoit-Wolf shrinkage of each `Sigma_hat_e` towards the pooled
within-environment covariance, with a data-driven intensity.

`Sigma_0` is never shrunk towards the training covariance — the target domain
has shifted, and pooling it back would erase the very difference the
predictive variance is meant to express. It uses Ledoit-Wolf towards a scaled
identity instead.

Condition numbers and realised shrinkage intensities are returned in
`fit$diagnostics` and recorded in every simulation output row, so the
sensitivity of the results to covariance estimation can be read off directly
rather than assumed away.

### Plugging in versus inferring the covariances

`stan/gi_hd.stan` conditions on `Sigma_hat` as if it were known, which is what
the manuscript does. `stan/gi_hd_fullcov.stan` instead gives
`Sigma_1, ..., Sigma_E` an inverse-Wishart hierarchy centred on a common
`Sigma_bar`, and `Sigma_0` its own LKJ prior, so their estimation error
propagates into the posterior for `gamma`.

The motivating worry is that because
`gamma = (within-environment slope) - Sigma^{-1} K`, error in `Sigma_hat`
lands directly on `gamma`, and under strong confounding `||Sigma^{-1} K||` is
large enough that this could be a first-order contribution rather than the
asymptotically negligible one the plug-in argument assumes.

**Measured, it is not.** Over 20 matched replications at `p = 6`, `E = 7`,
`n_e = 200` with strong confounding (`tests/test_plugin_vs_fullcov.R`):

| model | `gamma` coverage | RMSE | FDP | seconds |
|---|---|---|---|---|
| plug-in | **0.750** | 0.261 | 0.200 | 107 |
| full covariance | 0.692 | 0.318 | 0.277 | 288 |

Inferring the covariances was slightly *worse* on every metric at 2.7x the
runtime. The model is kept because being able to report that removing the
plug-in assumption changes nothing is a stronger answer to Referee 1's point
(4) than arguing the error is asymptotically negligible — but it is not a fix,
and the manuscript should not present it as one.

Select between them with `--model=gi_hd` (default) or `--model=gi_hd_fullcov`.
Both consume the same data list; the full-covariance model additionally reads
the scatter matrices, which `prepare_bgi_data()` always supplies.

Both keep the sufficient-statistic likelihood, so both cost `O(E p^3)` per
gradient evaluation.

---

## What the selection rules control

The estimand is

```
pa(Y) = { j : gamma*_j != 0 }
```

— the covariates entering the structural assignment for `Y` with a nonzero
coefficient. It is a statement about `gamma*`, identified jointly with `K*`
once the environment means span `R^p`. It is *not* a statement about
conditional independence, nor about invariance of `Y | X_S` across
environments, which is what the invariance literature targets and what hidden
confounding destroys.

None of the rules in `R/selection.R` is a frequentist test. Each thresholds a
posterior probability:

| Rule | Selects `j` when | Bounds |
|---|---|---|
| `sign` | `lfsr_j < alpha` | posterior probability of a **sign** error, per coordinate |
| `ci` | `0` outside the central `1 - alpha` interval | same family; equals `sign` at `alpha / 2` |
| `bayes_fdr` | largest initial segment of sorted `lfsr` with mean `<= q` | posterior expected **sign-error** proportion |
| `rope` | posterior mass outside `(-delta, delta)` exceeds `1 - alpha` | posterior probability that `gamma_j` is negligible |

where `lfsr_j = min{P(gamma_j > 0 | D), P(gamma_j < 0 | D)}`.

### `alpha` is a loss ratio, not a level

The `sign` rule is the Bayes rule for a **three-action** decision — declare
positive, declare negative, or decline to call — under

```
L(+, gamma) = 1{gamma < 0},  L(-, gamma) = 1{gamma > 0},  L(abstain) = lambda
```

Declaring the majority sign costs exactly `lfsr_j`, so the Bayes rule is
"declare `sgn(gamma_j)` when `lfsr_j < lambda`, else abstain", with
`alpha = lambda`. Two actions would not produce a threshold at all: under 0-1
sign loss with only `{+, -}` the Bayes rule never abstains and selects
everything.

So `alpha` is the cost of declining to call a coordinate relative to the cost
of getting its sign wrong. Setting `alpha = 0.05` asserts a sign error is
twenty times as costly as an abstention — a declared preference, which is
precisely why the thresholded quantity is not a p-value.

**A limitation worth stating plainly.** Under a continuous prior on `gamma`,
`P(gamma_j = 0 | D) = 0` for every `j`, so no rule built from the posterior of
`gamma` alone can bound the probability of selecting a coordinate whose true
value is *exactly* zero. The `lfsr` bounds the probability of attributing the
wrong **sign**, which is a different quantity — for `gamma_j = 0` exactly, no
sign is correct.

Note this bites harder in a structural causal model than in the settings the
lfsr literature usually addresses: a non-parent has `gamma*_j = 0` *exactly*,
by construction of the structural assignment, so the zeros are structural
rather than an idealisation of "small".

The simulation makes this concrete, and the aggregator now prints both side by
side so they cannot be confused. In `conf2_strong`, the `sign` rule attains

| quantity | value |
|---|---|
| realised sign error among selected true parents | **0.000** |
| exact zeros among all selected (FDP) | **0.190** |

The method is essentially perfect at what the rule controls and mediocre at
what it does not. `bayes_fdr` likewise holds its posterior expected sign-error
rate at or below `q = 0.1` by construction while its FDP over exact zeros runs
around 0.43. Both numbers are right; they measure different things, and
reporting only one of them is what makes the claims read as overstated.

Two coherent resolutions, and the manuscript should choose one explicitly:
restate the target as "not negligible" and use `rope`, or introduce an atom at
zero (spike-and-slab) so that inclusion probabilities exist. Section 2.1
declines the second; the cost of declining it is that FDR over exact nulls is
not available from this posterior.

Frequentist behaviour of every rule is measured, not assumed, in
`scripts/01_sim_support_recovery.R`.

---

## The support-recovery simulation

`scripts/01_sim_support_recovery.R` addresses Referee 1's point (5).

The data-generating process gives `s0` of the `p` slopes a nonzero value and
sets the rest to exactly zero, and the hidden confounder loads on **all** `p`
covariates. False discoveries are therefore possible and attributable, which
they are not under a DGP in which every slope is nonzero.

Reported per replication: TPR, FPR, realised FDP, familywise error
`1{S_hat not a subset of pa(Y)}`, exact recovery, Jaccard, MCC, sign errors
among selected coordinates, credible-interval coverage for `gamma`, and
predictive coverage in the shifted target domain — from the same fits, so
selection and calibration cannot be reported from different runs.

Procedures compared:

| Name | Description |
|---|---|
| `bgi_sign`, `bgi_ci`, `bgi_fdr`, `bgi_rope` | this paper, four decision rules |
| `ols`, `ols_bh` | pooled OLS, t-test and Benjamini-Hochberg |
| `pooled_gi`, `pooled_gi_bh` | frequentist GI as a pooled least squares fit |
| `iv` | 2SLS with the environment indicators as instruments |
| `anchor_g2`, `anchor_g8`, `anchor_g32` | anchor regression at three fixed strengths |
| `anchor_oracle` | anchor regression with `gamma_anchor` tuned on target labels |
| `icp` | invariant causal prediction |

### The baselines lie on one path

`tests/test_baseline_identities.R` verifies to machine precision, for both
exactly-identified and over-identified designs:

```
anchor(gamma_anchor = 1)          ==  pooled OLS
anchor(gamma_anchor -> infinity)  ==  2SLS with environment instruments
pooled GI                         ==  2SLS with environment instruments
```

The last is an exact algebraic identity (Frisch-Waugh-Lovell; proof in
`R/baselines.R`). So the "pooled OLS on `X` and the environment means"
baseline Referee 2 proposes and the multi-source IV comparison they ask for
are the **same estimator**, and the frequentist GI estimate of `gamma` is IV
with the environment as instrument — which Section 3.1.1 of the manuscript
already observes in one dimension.

The practical reading: anchor regression interpolates from OLS to 2SLS, and
GI's `gamma` sits at the 2SLS endpoint.

**The coincidence is with `gamma` only, and `gamma` alone is the wrong thing
to predict with.** The target-domain conditional mean is

```
E[Y_0 | X_0] = alpha + gamma' X_0 + K' Sigma_0^{-1} (X_0 - mu_0)
```

`gamma*` is the *causal* slope, not the population least-squares slope in the
target domain; under hidden confounding those differ, and the `K` term is
exactly the correction between them. IV recovers `gamma*` and stops, so
predicting with it discards the correction — causally right, predictively
wrong. `K` is intrinsic to GI and has no counterpart in any of these
baselines, which is why the comparison must be reported on predictive
quantities and not only on `gamma`.

`anchor_oracle` picks `gamma_anchor` by target-domain RMSE using the target
*labels*. That is not a usable procedure — the target is unlabelled, which is
the whole setting — and it is reported as the ceiling any practical tuning of
anchor regression would have to reach. It is the fair test of the claim that
GI needs no hyperparameter tuning.

The scenario grid crosses confounding strength, identifiability geometry
(`strong` versus near-collinear environment means, Section 2.3) and `n_e`.

---

## The case study (Section 3.2)

Referee 2 made three separable criticisms of the submitted Section 3.2, all
correct: coverage was computed over individuals within one held-out domain
whose intervals all share a single posterior, so those indicators are strongly
dependent and the quoted precision is an illusion; one held-out domain cannot
validate a claim about generalising to new domains; and the data are
proprietary, so nothing can be replicated.

`R/loeo.R` and `scripts/05_case_study.R` answer the first two by holding out
every environment in turn and taking the standard error from the spread
*across* environments. The third is answered by running the same analysis on a
public dataset.

### Two datasets, one pipeline

| | quiron | BRFSS 2023 |
|---|---|---|
| response | BMI | BMI (`_BMI5`) |
| environment | Spanish province of work | US state |
| `N` | 507,123 | 295,773 |
| `p` | 12 | 13 |
| `E` | 52 | 52 |
| `E_train - (p + 1)` | 38 | 37 |
| public? | no | yes |

`load_case_data()` in `R/case_data.R` handles both and returns the same four
elements, so nothing downstream knows which dataset it has.

```bash
# BRFSS: download LLCP2023XPT.zip from
#   https://www.cdc.gov/brfss/annual_data/annual_2023.html
# then, with no R or Python dependencies at all:
python3 scripts/06_prepare_brfss.py --xpt=../data/brfss/LLCP2023.XPT \
                                    --out=../data/brfss/brfss2023_case.csv

Rscript scripts/05_case_study.R --data=../data/brfss/brfss2023_case.csv \
                                --summary-only          # design summary only

sbatch slurm/05_case_convergence.sh brfss                # gate: read Rhat
sbatch --array=1-52 slurm/06_case_study.sh brfss         # one fold per task
Rscript scripts/07_aggregate_case_study.R --in=results/case_brfss/folds
```

`scripts/06_prepare_brfss.py` is a SAS XPORT v5 reader written against the
standard library alone, because the cluster's R has neither `foreign` nor
`haven` and a replication package should not acquire a dependency for one
file-format step. It checks the row count it produces against the count implied
by the file geometry, which is what catches a truncated download.

### How faithful the mirror is

BRFSS rotates its content, and no recent year carries both a graded
physical-activity measure and the sleep questions. 2023 is used because
physical activity is the variable the manuscript's selection claim is about;
2022 and 2024 collapse it to a single binary. Sleep quality and duration
therefore have no counterpart, and e-cigarette and smokeless-tobacco use stand
in as the two additional lifestyle factors.

| quiron | BRFSS 2023 | columns |
|---|---|---|
| `consumo_alcohol`, 6 levels | drinks/week `_DRNKWK2` cut at 0, 1, 3, 7, 14 | 5 |
| `fumador`, 3 levels | `_SMOKER3`, 4 levels, reference never | 3 |
| `af`, numeric 1-5 | `_PACAT3` reversed to a numeric 1-4 score | 1 |
| `calidad_sueño` | none; `_PASTRNG` muscle strengthening instead | 1 |
| `duracion_sueño` | none; `USENOW3` and `_CURECI2` instead | 3 |

`load_brfss_data(pa_numeric = FALSE)` expands the activity score into dummies,
which is the sensitivity check for treating an ordinal scale as continuous —
the same open question the quiron `af` raises.

### The mirror is better than the original in one respect

Sex, age, income, education, race, diabetes and self-rated health are *in* the
BRFSS file. They are excluded from the design on purpose, exactly as sex,
cholesterol and glucose are absent from the quiron data — but here the
confounding story can be checked rather than asserted:

```
BMI ~ the 7 excluded confounders     R^2 = 0.106
BMI ~ the 13 lifestyle covariates    R^2 = 0.047
```

and they predict the covariates too, most strongly physical activity
(`R^2 = 0.127`) — the variable the paper reports as causal. So the omitted
variables drive both the response and the exposures, and they load hardest on
the one coordinate the substantive claim depends on. `load_case_data()` returns
them in `confounders` for this check; they never enter `x`.

Two features of the public data to state in the paper rather than leave
implicit. Kentucky and Pennsylvania did not field enough of the 2023 survey to
appear, so `E = 52` is 48 states plus DC, Guam, Puerto Rico and the Virgin
Islands. And the complete-case filter is mildly state-dependent — retention
runs from 0.611 (Oklahoma) to 0.767 (Montana), sd 0.032 — which shifts the
`mu_e` slightly for reasons that are not substantive.

---

## Running on CESGA FinisTerrae III

```bash
cd new_code

sbatch slurm/00_compile.sh                 # once; wait for it to finish
bash  slurm/submit_support_recovery.sh --n-rep=20 --max-concurrent=40
```

`submit_support_recovery.sh` queries the task table for the array size,
submits the array, and chains the aggregation job with an `afterany`
dependency.

Design notes:

- **One array task = one replication.** Tasks are independent, each writes its
  own CSV atomically, and completed tasks are skipped on resubmission, so a
  pre-empted or failed task can simply be resubmitted.
- **Models are compiled once, up front.** Compiling inside workers makes every
  worker pay the cost and lets concurrent workers race on rstan's on-disk
  cache.
- **One MCMC chain per allocated core**, with `OMP_NUM_THREADS=1` and friends,
  so nested BLAS threading does not oversubscribe the node.
- **Seeds are pure functions of `(base_seed, task_id)`** via L'Ecuyer-CMRG
  streams, so a replication reproduces identically whether it runs in the
  array, serially, or on its own.

Partition `short` (6 h limit) is enough at the default grid size; use
`--partition=medium` for larger runs.

Resubmit only what failed:

```bash
sbatch --array=17,42,103 slurm/01_support_recovery.sh --n-rep=20
```

---

## Changes from the code used for the submitted version

Recorded here because several are substantive rather than cosmetic.

**Correctness**

1. Posterior predictive draws for `Y_0` used the training residual sd
   `sigma_y` rather than `S_0 = sqrt(sigma_y^2 - K' Sigma_0^{-1} K)`. They now
   use `S_0`.
2. `K` was a free slope on `(X - mu_e)` in the likelihood but was treated as a
   covariance in the prediction block, reconciled by multiplying through by
   `avg_var_X %*% ivar_X0`. That is only valid when every `Sigma_e` is equal.
   `K` is now a covariance throughout, with slope `Sigma_e^{-1} K` per
   environment.
3. `sds`, the scale of the hierarchical prior on `mu_e`, was declared but
   given no prior, i.e. an improper flat prior on `(0, inf)`. It now has a
   half-normal prior.
4. The causal Mahalanobis condition was not enforced, so conditional variances
   could go negative during sampling. It is now imposed exactly.
5. `hmu` was set to the average of the *true* environment means
   (`main_claude.R:114`), which do not exist outside a simulation. It is now
   `colMeans(x)`, the pooled sample mean, as the paper specifies. The
   numerical effect is negligible — the two differ by `O(N^{-1/2})`, and
   `hmu` is only the prior mean for `mu_e`, which the likelihood dominates —
   but the simulation script could not be pointed at real data as written,
   and the real-data script (`main_quiron.R:69`) already did this correctly,
   so the two disagreed.
6. The environment covariances were plugged in without regularisation, which
   under strong confounding biases `gamma` through the variance constraint
   (see above). Shrinkage towards the pooled covariance is now the default.

**Reproducibility**

7. No seed was set anywhere in the simulation scripts. Seeds are now
   deterministic functions of `(base_seed, task_id)`.
8. `#SBATCH --job-name="$1_$2"` in `run_simulation.sh` never expanded — SBATCH
   directives are comments. Job naming is now handled at submission.
9. The Stan file mixed current `array[N] int` syntax with the `matrix[P,P]
   var_X[NZ]` form removed in Stan 2.33. Current syntax throughout.
10. The model was recompiled inside every `mclapply` worker. It is compiled
    once and cached.

**Performance**

11. The likelihood looped over all `N` observations with a `multi_normal` call
    per row. It now uses per-environment sufficient statistics, so one
    gradient evaluation costs `O(E p^3)` instead of `O(N p^2)` — the cost no
    longer grows with sample size.
12. `corr_matrix` + `lkj_corr` replaced by `cholesky_factor_corr` +
    `lkj_corr_cholesky`; `Sigma_0^{-1}(X_0 - mu_0)` precomputed in
    `transformed data` since it is free of parameters.

`stan/gi_hd_reference.stan` is a direct per-observation implementation of the
same posterior. `tests/test_fast_vs_reference.R` fits both to identical data
and checks that the posterior means agree to within Monte Carlo error and the
posterior sds to within 10%, which is what makes item 11 safe to rely on.

---

## Open issue: `gamma` under-coverage under strong confounding

Recorded here rather than buried, because it bears directly on what the
manuscript can claim.

In the `conf2_strong` scenario (`p = 6`, `s0 = 3`, `E = 7`, `n_e = 200`,
strong hidden confounding), 95% credible intervals for `gamma` cover at about
0.70, and the `sign` selection rule reaches a false discovery proportion of
0.19 against the frequentist pooled GI baseline's 0.11. The residual gap is
bias, not variance: the posterior sd is smaller than the RMSE.

Ruled out, each over 20 replications:

| candidate | result |
|---|---|
| sampler geometry | real but insufficient — on the same grid and seeds, non-centring took coverage 0.479 to 0.700, and rescued the weak-identifiability scenarios where previously *every* fit failed diagnostics |
| plug-in covariance uncertainty | not the cause — inferring `Sigma_e`, `Sigma_0` gave 0.692 against 0.750, at 2.7x runtime (`tests/test_plugin_vs_fullcov.R`) |
| `mu_e` prior shrinkage | not the cause — flattening it gives 0.775 against 0.750, inside Monte Carlo error |
| ridge shrinkage of `gamma` | not the cause — flattening it gives 0.775, and costs RMSE and divergences |

The decisive observation: **the frequentist GI estimator, with no prior and no
shrinkage anywhere, under-covers worse (0.633) than the Bayesian one (0.750)**.
The shortfall is not a Bayesian artefact and no prior will remove it.

The mechanism is regression dilution in the design. `Sigma^{-1} K` is
identified by regressing the `E` environment intercepts on the `E` environment
means; at `E = p + 1` that regression is exactly saturated — `p + 1`
observations, `p + 1` parameters, zero residual degrees of freedom. The
`mu_hat_e` carry error of order `sqrt(Sigma / n_e)`, which attenuates
`Sigma^{-1} K` towards zero with nothing to average it away, and since
`gamma` is the within-environment slope minus it, the attenuation lands on
`gamma`.

Whether BGI's treatment of the `mu_e` as parameters recovers part of this is
**not established**: two runs at the same `E` and `n_e` gave BGI 0.750 vs
pooled GI 0.633, and BGI 0.717 vs pooled GI 0.733. The ordering reverses and
the standard error is near 0.08. Both under-cover; that is the supportable
claim.

`scripts/04_sim_environment_budget.R` tests the prediction that coverage
recovers with `E` and with `n_e`:

```bash
sbatch --array=1-8 slurm/04_env_budget.sh --n-rep=20
```

The practical guidance this yields is what Referee 1 asks for in point (3):
`E >= p + 1` is necessary for identification but not sufficient for reliable
selection, and the quantity to monitor is the residual degrees of freedom
`E - (p + 1)`, not merely whether the environment means span `R^p`.

## Reproducing individual results

```bash
# One specific replication, in isolation.
Rscript scripts/01_sim_support_recovery.R --task=17 --n-rep=20

# Sensitivity to the covariance estimator.
Rscript scripts/01_sim_support_recovery.R --task=all --n-rep=5 --small \
        --cov-method=sample --out=results/sens_sample

# Aggregate any results directory.
Rscript scripts/02_aggregate_support_recovery.R --in=results/sens_sample
```

Flags accepted by `01_sim_support_recovery.R`: `--task`, `--n-rep`, `--small`,
`--alpha`, `--q`, `--seed`, `--chains`, `--iter`, `--cov-method`, `--out`.

## Session information used for the reported runs

```
R 4.4.2, rstan 2.32.6 (Stan 2.32.2), mvtnorm 1.3.3
CESGA FinisTerrae III, module cesga/system R/4.4.2
```
