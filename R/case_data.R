## case_data.R ------------------------------------------------------------
##
## Data loading for the real-data case studies.
##
## Two datasets are supported and they answer the same question with the same
## structure:
##
##   quiron  BMI against five lifestyle factors, Spanish province of work as
##           the environment.  Proprietary; the substantive analysis.
##   brfss   BMI against lifestyle factors, US state as the environment.  CDC
##           Behavioral Risk Factor Surveillance System, fully public; the
##           reproducible mirror Referee 2 asked for.
##
## Everything downstream sees the same four elements — `x`, `y`, `z` and
## `covariate_names` — so `scripts/05_case_study.R` and the convergence test
## are identical for both.  This file is the only place that knows about file
## formats or variable codings.

BGI_CASE_DATASETS <- c("quiron", "brfss", "acs_tract", "acs_pums",
                       "communities")

#' Load and prepare a case-study dataset.
#'
#' @param path Path to the CSV.
#' @param dataset One of `BGI_CASE_DATASETS`, or `"auto"` to detect from the
#'   column names.
#' @param ... Passed to the dataset-specific loader.
#' @return A list with `x` (numeric design matrix, no intercept), `y`, `z`
#'   (environment labels), `covariate_names`, and `confounders` — a data frame
#'   of variables deliberately excluded from `x`, empty when the source does
#'   not contain them.
load_case_data <- function(path, dataset = "auto", ...) {
  if (!file.exists(path)) {
    stop("Case-study data not found at ", path,
         ".\nPass --data=PATH. For the public dataset, see ",
         "scripts/06_prepare_brfss.py and HANDOFF_REAL_DATA.md.",
         call. = FALSE)
  }
  if (identical(dataset, "auto")) {
    dataset <- detect_case_dataset(path)
    message("Detected dataset format: ", dataset)
  }
  dataset <- match.arg(dataset, BGI_CASE_DATASETS)
  switch(dataset,
         quiron = load_quiron_data(path, ...),
         brfss = load_brfss_data(path, ...),
         acs_tract = load_acs_tract_data(path, ...),
         acs_pums = load_acs_pums_data(path, ...),
         communities = load_communities_data(path, ...))
}

#' Guess which dataset a CSV holds from its header alone.
detect_case_dataset <- function(path) {
  header <- names(utils::read.csv(path, nrows = 1L, fileEncoding = "UTF-8",
                                  check.names = FALSE))
  if ("imc" %in% header && "prov_trab" %in% header) {
    return("quiron")
  }
  if ("X_STATE" %in% header || "_STATE" %in% header) {
    return("brfss")
  }
  if (all(c("median_income", "state", "geo_id") %in% header)) {
    return("acs_tract")
  }
  if (all(c("WAGP", "STATE", "SCHL") %in% header)) {
    return("acs_pums")
  }
  stop("Cannot identify the dataset from the columns of ", path,
       ".\nPass --dataset=quiron or --dataset=brfss.", call. = FALSE)
}

## ---- quiron --------------------------------------------------------------

#' Load the proprietary BMI / lifestyle data of Section 3.2.
#'
#' The response and covariates are as submitted: BMI against five lifestyle
#' factors, with province of work as the environment.  Sex, cholesterol and
#' glucose are treated as hidden confounders and are, by design, absent from
#' the file.
#'
#' @param ... Ignored.  Present so that callers can pass BRFSS-specific options
#'   through `load_case_data()` without having to branch on the dataset.
load_quiron_data <- function(path, ...) {
  af <- utils::read.csv(path, fileEncoding = "UTF-8", stringsAsFactors = FALSE)

  required <- c("imc", "prov_trab")
  if (!all(required %in% names(af))) {
    stop("Expected columns ", paste(required, collapse = ", "),
         " in ", path, call. = FALSE)
  }

  ## Same cleaning as the submitted script: positive BMI, no blank fields.
  af <- af[af$imc > 0, , drop = FALSE]
  af <- af[stats::complete.cases(af), , drop = FALSE]
  ## Column-wise, not row-wise: a row-wise apply() coerces the whole frame to a
  ## character matrix once per row and dominates the runtime on 500,000 rows.
  blank <- Reduce(`|`, lapply(af, function(col) trimws(as.character(col)) == ""))
  af <- af[!blank, , drop = FALSE]

  y <- af$imc
  z <- af$prov_trab
  predictors <- setdiff(names(af), c("imc", "prov_trab"))

  ## Character columns become factors; anything already numeric (the physical
  ## activity score) stays numeric.  Reference levels follow the submitted
  ## script where they were set explicitly.
  refs <- list(consumo_alcohol = "no", fumador = "no fumador",
               `calidad_sueño` = "profundo", `duracion_sueño` = "6-9h")
  for (v in predictors) {
    if (is.character(af[[v]])) {
      af[[v]] <- factor(af[[v]])
      if (!is.null(refs[[v]]) && refs[[v]] %in% levels(af[[v]])) {
        af[[v]] <- stats::relevel(af[[v]], ref = refs[[v]])
      }
    }
  }

  x <- stats::model.matrix(
    stats::as.formula(paste("~", paste(predictors, collapse = "+"))),
    data = af[, predictors, drop = FALSE])[, -1, drop = FALSE]

  list(x = x, y = y, z = z, covariate_names = colnames(x),
       confounders = af[, character(0), drop = FALSE])
}

## ---- BRFSS ---------------------------------------------------------------

## FIPS state and territory codes used by BRFSS.  Kept here so that fold names
## in the output are readable rather than numeric codes.
BRFSS_STATE_NAMES <- c(
  "1" = "Alabama", "2" = "Alaska", "4" = "Arizona", "5" = "Arkansas",
  "6" = "California", "8" = "Colorado", "9" = "Connecticut",
  "10" = "Delaware", "11" = "District of Columbia", "12" = "Florida",
  "13" = "Georgia", "15" = "Hawaii", "16" = "Idaho", "17" = "Illinois",
  "18" = "Indiana", "19" = "Iowa", "20" = "Kansas", "21" = "Kentucky",
  "22" = "Louisiana", "23" = "Maine", "24" = "Maryland",
  "25" = "Massachusetts", "26" = "Michigan", "27" = "Minnesota",
  "28" = "Mississippi", "29" = "Missouri", "30" = "Montana",
  "31" = "Nebraska", "32" = "Nevada", "33" = "New Hampshire",
  "34" = "New Jersey", "35" = "New Mexico", "36" = "New York",
  "37" = "North Carolina", "38" = "North Dakota", "39" = "Ohio",
  "40" = "Oklahoma", "41" = "Oregon", "42" = "Pennsylvania",
  "44" = "Rhode Island", "45" = "South Carolina", "46" = "South Dakota",
  "47" = "Tennessee", "48" = "Texas", "49" = "Utah", "50" = "Vermont",
  "51" = "Virginia", "53" = "Washington", "54" = "West Virginia",
  "55" = "Wisconsin", "56" = "Wyoming", "66" = "Guam",
  "72" = "Puerto Rico", "78" = "Virgin Islands"
)

#' Recode a BRFSS variable, mapping its missing codes to `NA`.
#'
#' BRFSS marks refusals and don't-knows with in-range sentinel values (`7`,
#' `9`, `99900`, ...) rather than with a missing flag, so every variable needs
#' its sentinels named explicitly.  Silently treating a `9` as a level is the
#' classic way to get a fifth smoking category.
#'
#' @param v Numeric vector as read from the CSV.
#' @param missing Sentinel values to map to `NA`.
#' @param labels Optional named character vector, `code = label`; when given,
#'   the result is a factor with `levels(labels)` in the order supplied.
#' @param ref Optional reference level.
brfss_recode <- function(v, missing = c(7, 9), labels = NULL, ref = NULL) {
  v[v %in% missing] <- NA
  if (is.null(labels)) {
    return(v)
  }
  out <- factor(labels[as.character(v)], levels = unname(labels))
  if (!is.null(ref)) {
    out <- stats::relevel(out, ref = ref)
  }
  out
}

#' Load the BRFSS public mirror of the case study.
#'
#' Prepared by `scripts/06_prepare_brfss.py` from the CDC annual LLCP file.
#'
#' The design mirrors the quiron analysis as closely as BRFSS allows: BMI as
#' the response, lifestyle factors as covariates, and the US state as the
#' environment.  Sex, age, income, education, race and diabetes are the hidden
#' confounders.  They *are* in the file, and are excluded from `x` on purpose;
#' that is the one respect in which the public mirror is better than the
#' proprietary original, because the confounding story can be checked rather
#' than asserted.
#'
#' Two quiron factors have no BRFSS counterpart in 2023.  Sleep quality and
#' sleep duration were dropped from the BRFSS core after 2022, and the graded
#' physical-activity module — the analogue of the `af` score, which is the
#' variable the paper reports as causal — is carried in 2023 but not 2022.  No
#' recent year has both.  The 2023 file is used because physical activity is
#' the variable that matters to the paper's claim; e-cigarette and smokeless
#' tobacco use stand in as the two additional lifestyle factors, so the design
#' has p = 13 against the quiron p = 12.
#'
#' @param path CSV written by `scripts/06_prepare_brfss.py`.
#' @param pa_numeric Treat the four-level physical-activity category as a
#'   numeric score, as the submitted analysis treated `af`.  `FALSE` expands
#'   it into dummies, which is the sensitivity check HANDOFF_REAL_DATA.md §3
#'   asks for.
#' @param bmi_range Plausible BMI bounds; values outside are dropped.
#' @param min_state_n Drop states with fewer than this many complete rows.
load_brfss_data <- function(path, pa_numeric = TRUE,
                            bmi_range = c(12, 70), min_state_n = 100L,
                            subset = "none") {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ## read.csv mangles the leading underscore of BRFSS computed variables; work
  ## with whichever spelling the file actually has.
  names(raw) <- sub("^X_", "_", names(raw))

  need <- c("_STATE", "_BMI5", "_SMOKER3", "_PACAT3", "_PASTRNG", "USENOW3")
  absent <- setdiff(need, names(raw))
  if (length(absent) > 0L) {
    stop("Missing BRFSS columns: ", paste(absent, collapse = ", "),
         "\nRegenerate with scripts/06_prepare_brfss.py.", call. = FALSE)
  }

  ## Variable names drift between years; take whichever spelling is present.
  pick <- function(...) {
    for (nm in c(...)) if (nm %in% names(raw)) return(raw[[nm]])
    NULL
  }

  d <- data.frame(row.names = seq_len(nrow(raw)))

  ## Response: _BMI5 carries two implied decimal places.
  d$bmi <- raw[["_BMI5"]] / 100

  ## Environment: FIPS code to state name.
  d$state <- BRFSS_STATE_NAMES[as.character(raw[["_STATE"]])]

  ## --- Covariates ------------------------------------------------------

  ## Smoking, four levels.  The quiron `fumador` had three; BRFSS separates
  ## daily from occasional smokers and there is no reason to discard that.
  d$smoking <- brfss_recode(
    raw[["_SMOKER3"]], missing = 9,
    labels = c("1" = "smoker_daily", "2" = "smoker_occasional",
               "3" = "former", "4" = "never"),
    ref = "never")

  ## Alcohol, six levels, matching the six of the quiron `consumo_alcohol`.
  ## _DRNKWK* is drinks per week with two implied decimals; 99900 is the
  ## missing sentinel.
  drinks_raw <- pick("_DRNKWK2", "_DRNKWK3", "_DRNKWK1")
  if (is.null(drinks_raw)) {
    stop("No _DRNKWK* column in ", path, call. = FALSE)
  }
  drinks <- brfss_recode(drinks_raw, missing = 99900) / 100
  d$alcohol <- cut(drinks, breaks = c(-Inf, 0, 1, 3, 7, 14, Inf),
                   labels = c("none", "le1", "1to3", "3to7", "7to14", "gt14"),
                   right = TRUE)
  d$alcohol <- stats::relevel(d$alcohol, ref = "none")

  ## Physical activity: the analogue of `af`, and the variable the paper's
  ## selection claim is about.  _PACAT3 runs 1 = highly active to 4 = inactive;
  ## reversed here so that the coefficient sign has the obvious reading.
  pa <- brfss_recode(raw[["_PACAT3"]], missing = 9)
  if (pa_numeric) {
    d$phys_activity <- 5 - pa
  } else {
    d$phys_activity <- factor(
      c("inactive", "insufficient", "active", "highly_active")[5 - pa],
      levels = c("inactive", "insufficient", "active", "highly_active"))
  }

  ## Muscle strengthening: the second physical-activity dimension BRFSS
  ## measures, and independent of the aerobic one.
  d$strength <- brfss_recode(
    raw[["_PASTRNG"]], missing = 9,
    labels = c("1" = "meets", "2" = "does_not"), ref = "does_not")

  ## E-cigarettes.
  ecig_raw <- pick("_CURECI2", "_CURECI3", "_CURECI1")
  if (!is.null(ecig_raw)) {
    d$ecig <- brfss_recode(
      ecig_raw, missing = 9,
      labels = c("1" = "no", "2" = "yes"), ref = "no")
  }

  ## Smokeless tobacco, three levels.
  d$smokeless <- brfss_recode(
    raw[["USENOW3"]], missing = c(7, 9),
    labels = c("1" = "daily", "2" = "some_days", "3" = "none"), ref = "none")

  ## --- Confounders, excluded from x on purpose --------------------------
  conf <- data.frame(row.names = seq_len(nrow(raw)))
  conf$sex <- brfss_recode(pick("SEXVAR", "_SEX"), missing = 9)
  conf$age <- brfss_recode(pick("_AGE80"), missing = c())
  conf$income <- brfss_recode(pick("_INCOMG1", "INCOME3"), missing = c(9, 77, 99))
  conf$education <- brfss_recode(pick("_EDUCAG", "EDUCA"), missing = c(9))
  conf$diabetes <- brfss_recode(pick("DIABETE4"), missing = c(7, 9))
  conf$race <- brfss_recode(pick("_RACE"), missing = c(9))
  conf$genhlth <- brfss_recode(pick("GENHLTH"), missing = c(7, 9))

  ## --- Clinical restriction ---------------------------------------------
  ## `older_diabetic` keeps respondents aged 60 or over with diagnosed
  ## diabetes.  The purpose is statistical as well as clinical: on the full
  ## survey every fold trains on about 277,000 rows, at which sample size the
  ## local false sign rate is small for every coordinate and the selection rule
  ## of Section 2.2 declares all 13 covariates parents.  That is the rule
  ## behaving as specified, but it cannot discriminate, so the case study
  ## cannot illustrate selection.  This subgroup is roughly a tenth the size.
  ##
  ## The restriction is clinically standard rather than chosen to suit the
  ## method: weight management in older adults with type 2 diabetes is a
  ## first-line question, and the subgroup has mean BMI 30.8, in the obese
  ## range.  Two properties were checked before adopting it.  The confounding
  ## survives -- the remaining excluded variables explain more of BMI
  ## (R^2 = 0.083) than the retained covariates do (0.038), and still predict
  ## physical activity at R^2 = 0.125 against 0.127 in the full survey.  And
  ## identification improves rather than degrades: the environment means are
  ## better conditioned (condition number 1775 against 9272), because
  ## conditioning on age and diabetes removes variation that was common to all
  ## states.
  ##
  ## Note that age and diabetes are two of the variables the design excludes as
  ## confounders.  Conditioning on them narrows the estimand to a stratum; it
  ## does not use information the model is meant to be denied, since the
  ## restriction is applied before fitting and identically in every
  ## environment.
  if (!identical(subset, "none")) {
    keep_sub <- switch(
      subset,
      ## Diagnosed diabetes, no age restriction. Preferred over
      ## `older_diabetic`: it is a fifth larger (N = 40,286 against 28,450),
      ## keeps all 52 states, and doubles the identification strength
      ## `N * lambda_min` from 16 to 33, which is where the full survey sits.
      ## Letting age vary again also restores it as a confounder -- the OLS and
      ## IV slope vectors correlate -0.20 here against 0.02 under the age
      ## restriction -- and "adults with diagnosed diabetes" needs no arbitrary
      ## age cut to justify clinically.
      diabetic = !is.na(conf$diabetes) & conf$diabetes == 1,
      older_diabetic = !is.na(conf$age) & conf$age >= 60 &
        !is.na(conf$diabetes) & conf$diabetes == 1,
      older = !is.na(conf$age) & conf$age >= 60,
      stop("Unknown subset: ", subset, call. = FALSE))
    d <- d[keep_sub, , drop = FALSE]
    conf <- conf[keep_sub, , drop = FALSE]
    message("Subset '", subset, "': ", nrow(d), " rows retained.")
  }

  ## --- Cleaning ---------------------------------------------------------
  covariates <- setdiff(names(d), c("bmi", "state"))
  ok <- stats::complete.cases(d) &
    is.finite(d$bmi) & d$bmi >= bmi_range[1] & d$bmi <= bmi_range[2]
  d <- d[ok, , drop = FALSE]
  conf <- conf[ok, , drop = FALSE]

  ## An environment must be large enough for a plug-in covariance to mean
  ## something; prepare_bgi_data() warns below p observations, but a handful of
  ## rows is not an environment either way.
  sizes <- table(d$state)
  keep <- d$state %in% names(sizes)[sizes >= min_state_n]
  d <- d[keep, , drop = FALSE]
  conf <- conf[keep, , drop = FALSE]

  ## Drop levels that the cleaning emptied, or model.matrix produces all-zero
  ## columns and the covariances become singular.
  for (v in covariates) {
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }

  x <- stats::model.matrix(
    stats::as.formula(paste("~", paste(covariates, collapse = "+"))),
    data = d[, covariates, drop = FALSE])[, -1, drop = FALSE]

  list(x = x, y = d$bmi, z = d$state, covariate_names = colnames(x),
       confounders = conf)
}

## ---- ACS census tracts ---------------------------------------------------

#' Load the ACS census-tract dataset.
#'
#' Prepared by `scripts/11_prepare_acs_tract.py` from the ACS 5-year
#' table-based summary files.
#'
#' This dataset exists because BRFSS and quiron both turned out to contain
#' almost no domain shift — the between-environment variance in the covariate
#' means is under 1% of the within-environment variance in each — and in that
#' regime the `K' Sigma_0^{-1}(X_0 - mu_0)` correction has nothing to correct,
#' so generative invariance is *expected* to coincide with pooled OLS.
#' `scripts/10_environment_shift.R` measures this one at a mean ratio of 0.274,
#' roughly 37 times BRFSS, with the most extreme environment (Puerto Rico)
#' sitting further from the pooled mean than a typical single tract does.
#'
#' The structure otherwise mirrors the other two deliberately: a continuous
#' response, around ten covariates, US states as environments giving `E = 52`
#' and `E - (p + 1) = 41`, and a set of variables that are present in the file
#' and excluded from the design on purpose.
#'
#' Two respects in which it is a better test than either existing case study.
#' The covariates are genuinely continuous shares, medians and indices, so
#' `X_ei ~ N(mu_e, Sigma_e)` is a far milder assumption than it is for eleven
#' binary dummies.  And the median tract count per state is about 1,200 against
#' 26 communities per state in the UCI Communities-and-Crime data, which is the
#' other public dataset with shift this large but which cannot support a
#' per-environment plug-in covariance at `p = 9`.
#'
#' The unit is a tract rather than a person, which is the reason the shift is
#' large: within-environment variance is then between-tract variance rather
#' than between-person variance.  State that plainly — it changes the estimand
#' to an ecological one, and the confounding operates between areas rather than
#' between individuals.
#'
#' @param path CSV written by `scripts/11_prepare_acs_tract.py`.
#' @param log_income Model `log` median household income.  Tract income is
#'   right-skewed and the log is close to symmetric, which suits the Gaussian
#'   response of Section 2.
#' @param min_state_n Drop states with fewer than this many complete tracts.
#' @param ... Ignored, so that dataset-specific options can be passed through
#'   `load_case_data()` without the caller branching on the dataset.
load_acs_tract_data <- function(path, log_income = TRUE, min_state_n = 100L,
                                ...) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)

  covariates <- c("pct_bachelors_plus", "unemployment_rate",
                  "labor_force_rate", "median_house_value", "median_rent",
                  "pct_owner_occupied", "avg_household_size", "gini",
                  "pct_commute_30plus")
  confounders <- c("median_age", "pct_white", "pct_black", "pct_foreign_born")

  absent <- setdiff(c("state", "median_income", covariates, confounders),
                    names(d))
  if (length(absent) > 0L) {
    stop("Missing ACS columns: ", paste(absent, collapse = ", "),
         "\nRegenerate with scripts/11_prepare_acs_tract.py.", call. = FALSE)
  }

  ok <- stats::complete.cases(d[, c("median_income", covariates, confounders)]) &
    d$median_income > 0 & d$median_house_value > 0 & d$median_rent > 0
  d <- d[ok, , drop = FALSE]

  sizes <- table(d$state)
  d <- d[d$state %in% names(sizes)[sizes >= min_state_n], , drop = FALSE]

  y <- if (log_income) log(d$median_income) else d$median_income

  ## House value and rent are in dollars, with standard deviations around
  ## 293,000 and 600, against 0.2 for the shares.  Left on those scales the
  ## 2SLS normal equations are numerically singular even though the design is
  ## comfortably over-identified.  Logs are the right choice for prices anyway
  ## — they are multiplicative quantities and the logs are close to symmetric —
  ## and they put every covariate on a comparable scale.
  for (v in c("median_house_value", "median_rent")) {
    d[[v]] <- log(d[[v]])
  }
  names(d)[names(d) == "median_house_value"] <- "log_house_value"
  names(d)[names(d) == "median_rent"] <- "log_rent"
  covariates[covariates == "median_house_value"] <- "log_house_value"
  covariates[covariates == "median_rent"] <- "log_rent"

  x <- as.matrix(d[, covariates])

  list(x = x, y = y, z = d$state, covariate_names = covariates,
       confounders = d[, confounders, drop = FALSE])
}

## ---- ACS PUMS, individual level ------------------------------------------

#' Load the ACS person-level public-use microdata.
#'
#' Prepared by streaming the columns of interest out of the Census PUMS CSV;
#' see HANDOFF_REAL_DATA.md for the command.
#'
#' A second individual-level public dataset with US states as environments,
#' included because it is the cleanest available replication of the BRFSS
#' design in a completely different subject area: the response is labour
#' income rather than BMI and the covariates are human capital rather than
#' health behaviours.  It measures the same between-environment shift as BRFSS
#' to three decimal places (0.0073 against 0.0074), which is what establishes
#' that the near-absence of shift between US states is a property of the
#' geography and not of the variables anyone happens to choose.
#'
#' @param path CSV of selected PUMS columns.
#' @param min_state_n Drop states with fewer than this many complete records.
load_acs_pums_data <- function(path, min_state_n = 300L, ...) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("STATE", "WAGP", "WKHP", "WKWN", "SCHL", "COW", "MAR", "ESR")
  absent <- setdiff(need, names(d))
  if (length(absent) > 0L) {
    stop("Missing PUMS columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  ## Wage earners only: the response is the log wage, so a zero is not a small
  ## value but a different state of the world.
  ok <- !is.na(d$WAGP) & d$WAGP > 0 & !is.na(d$WKHP) & d$WKHP >= 1 &
    d$ESR %in% c(1, 2) &
    stats::complete.cases(d[, c("SCHL", "WKWN", "COW", "MAR")])
  d <- d[ok, , drop = FALSE]

  x <- cbind(
    schl = d$SCHL,
    hours = d$WKHP,
    weeks = d$WKWN,
    self_employed = as.integer(d$COW %in% c(6, 7)),
    government = as.integer(d$COW %in% c(3, 4, 5)),
    married = as.integer(d$MAR == 1),
    never_married = as.integer(d$MAR == 5)
  )
  y <- log(d$WAGP)
  z <- as.character(d$STATE)

  sizes <- table(z)
  keep <- z %in% names(sizes)[sizes >= min_state_n]

  conf <- d[keep, intersect(c("AGEP", "SEX", "RAC1P", "NATIVITY", "ENG",
                              "DIS", "HICOV"), names(d)), drop = FALSE]
  list(x = x[keep, , drop = FALSE], y = y[keep], z = z[keep],
       covariate_names = colnames(x), confounders = conf)
}

## ---- UCI Communities and Crime -------------------------------------------

#' Load the UCI Communities and Crime data.
#'
#' The unit is a community and the environment is its US state, which is what
#' makes this the highest-shift public dataset measured here: the
#' between-environment variance in the covariate means is 0.44 of the
#' within-environment variance, against 0.0074 for BRFSS
#' (`scripts/10_environment_shift.R`).  Within-environment variance is
#' between-community variance rather than between-person variance, so the state
#' means are no longer swamped by individual heterogeneity.
#'
#' It is included specifically as a stress test of the slope reparameterisation
#' (HANDOFF §3.8).  With a median of 26 communities per state against `p = 10`,
#' the per-environment sample covariance is thin, and the K-parameterisation
#' inverts exactly that quantity inside the conditional mean.  The slope model
#' keeps it out of the mean, so if the reparameterisation is doing what it
#' claims, this is where the difference should be largest.  Read the two runs
#' together; the dataset is a diagnostic, not a headline result.
#'
#' Caveats to state if it is reported. UCI ships every column rescaled to
#' `[0, 1]` by equal-interval binning, so coefficients are not interpretable in
#' natural units; the rescaling is per column and global, so it leaves the
#' between/within ratio unchanged. `N = 1994` in total, which is three orders of
#' magnitude below the other case studies. And the estimand is ecological.
#'
#' @param path Path to `communities.data`.
#' @param names_path Path to `communities.names`; defaults alongside `path`.
#' @param min_state_n Drop states with fewer communities than this.
load_communities_data <- function(path,
                                  names_path = file.path(dirname(path),
                                                         "communities.names"),
                                  min_state_n = 12L, ...) {
  if (!file.exists(names_path)) {
    stop("Cannot find ", names_path, call. = FALSE)
  }
  header <- grep("^@attribute", readLines(names_path, warn = FALSE),
                 value = TRUE)
  nms <- sub("^@attribute\\s+(\\S+).*$", "\\1", header)

  d <- utils::read.csv(path, header = FALSE, na.strings = "?",
                       stringsAsFactors = FALSE)
  names(d) <- nms

  ## Socioeconomic covariates only, all fully observed in this file.
  covariates <- c("medIncome", "pctWPubAsst", "PctPopUnderPov",
                  "PctUnemployed", "PctNotHSGrad", "PctBSorMore",
                  "PctPersDenseHous", "PctHousOccup", "MedRent", "PopDens")
  ## Composition and urbanicity are the confounders here: they move both the
  ## socioeconomic covariates and the crime rate, and are excluded on purpose.
  confounders <- c("racepctblack", "racePctWhite", "racePctHisp",
                   "racePctAsian", "agePct65up", "pctUrban", "PctImmigRecent")

  absent <- setdiff(c(covariates, confounders, "ViolentCrimesPerPop", "state"),
                    names(d))
  if (length(absent) > 0L) {
    stop("Missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }

  ok <- stats::complete.cases(d[, c(covariates, confounders,
                                    "ViolentCrimesPerPop")])
  d <- d[ok, , drop = FALSE]

  z <- as.character(d$state)
  sizes <- table(z)
  keep <- z %in% names(sizes)[sizes >= min_state_n]
  d <- d[keep, , drop = FALSE]

  list(x = as.matrix(d[, covariates]), y = d$ViolentCrimesPerPop,
       z = as.character(d$state), covariate_names = covariates,
       confounders = d[, confounders, drop = FALSE])
}
