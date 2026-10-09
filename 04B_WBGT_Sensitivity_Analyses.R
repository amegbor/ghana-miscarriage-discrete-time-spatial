# GitHub copy: portable paths; restricted data are not distributed.
# ==============================================================================
# 04B_WBGT_Sensitivity_Analyses.R
# Ghana miscarriage analysis: missing and corrected WBGT sensitivity analyses
#
# Designed for interactive RStudio on NYU Torch HPC.
# Input: rebuilt Script 02 outputs under <project root>/Results/Models
# Output: <project root>/Results; see README for path configuration
#
# This script does NOT rerun the completed primary HI analyses. It adds:
#   1. Population-averaged complementary-log-log GEE using WBGT
#   2. Pregnancy-level binary logistic spline model using early-pregnancy WBGT
#   3. Primary and sensitivity discrete-time WBGT mixed models:
#        - current-month WBGT, df = 4
#        - one-month-lagged WBGT
#        - prior cumulative-mean WBGT
#        - exclusion of losses before gestational month 2
#        - strict residence restriction
#        - logit link
#        - unweighted analysis
#        - linear raw-unit WBGT
#        - spline df = 3, 5, and 6
#   4. Fixed-effect, model-fit, omnibus-Wald, and local-contrast exports
#
# Heat index and WBGT are alternative thermal metrics and are never entered in
# the same model. Clean cooking fuel and current marital status are excluded.
# Region random effects are excluded. Mixed models retain mother and cluster
# random intercepts.
# ==============================================================================

# ---- 1. Packages and directories ---------------------------------------------

required_packages <- c(
  "dplyr", "purrr", "readr", "tibble", "stringr", "tidyr",
  "splines", "glmmTMB", "geepack", "emmeans", "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(splines)
  library(glmmTMB)
  library(geepack)
  library(emmeans)
  library(openxlsx)
})

# Run from the repository root, or set GHANA_MISCARRIAGE_ROOT explicitly.
# Private inputs and outputs can live entirely outside the Git checkout.
project_dir <- normalizePath(
  Sys.getenv("GHANA_MISCARRIAGE_ROOT", unset = "."),
  winslash = "/", mustWork = TRUE
)

results_dir <- file.path(project_dir, "Results")
models_dir <- file.path(results_dir, "Models")
tables_dir <- file.path(results_dir, "Tables")
checkpoint_dir <- file.path(
  models_dir,
  "WBGT_sensitivity_checkpoints_no_cooking_v2"
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Use a new checkpoint version whenever the formula or variable coding changes.
# This prevents older HI- or cooking-fuel-adjusted models from being reloaded.


# ---- 2. Load the analysis-ready discrete-time data ---------------------------

# Require the rebuilt primary complete-case dataset from corrected Script 02.
primary_file <- file.path(models_dir, "Ghana_DT_AnalysisReady_PrimaryCC.rds")
if (!file.exists(primary_file)) {
  stop("Run corrected Script 02 first. Missing dataset: ", primary_file)
}

message("Loading: ", primary_file)
DT <- readRDS(primary_file)
if (!isTRUE(attr(DT, "ghana_cooking_fuel_excluded"))) {
  stop("Run corrected Script 02 first; this dataset has not been rebuilt without cooking fuel.")
}
DT <- as.data.frame(DT)

# Support both WBT and WBGT naming while preserving WBT as the analysis name.
if (!"WBT" %in% names(DT) && "WBGT" %in% names(DT)) {
  DT$WBT <- DT$WBGT
}
if (!"WBT_lag1" %in% names(DT) && "WBGT_lag1" %in% names(DT)) {
  DT$WBT_lag1 <- DT$WBGT_lag1
}
if (!"WBT_prior_cummean" %in% names(DT) &&
    "WBGT_prior_cummean" %in% names(DT)) {
  DT$WBT_prior_cummean <- DT$WBGT_prior_cummean
}

required_core_variables <- c(
  "miscarriage_event", "preg_uid", "mother_uid", "CLUSTID",
  "gest_month", "WBT", "CO", "PM", "EVI"
)

missing_core_variables <- setdiff(required_core_variables, names(DT))
if (length(missing_core_variables) > 0L) {
  stop(
    "The analysis dataset is missing: ",
    paste(missing_core_variables, collapse = ", ")
  )
}

# Select a valid weight variable and scale it to mean 1 for numerical stability.
weight_candidates <- c("analysis_weight", "v005", "sampling_weight")
weight_variable <- weight_candidates[weight_candidates %in% names(DT)][1]

if (is.na(weight_variable)) {
  stop("No analysis weight was found.")
}

DT$analysis_weight <- as.numeric(DT[[weight_variable]])
DT$analysis_weight <- DT$analysis_weight /
  mean(DT$analysis_weight, na.rm = TRUE)

# Ensure binary outcome and grouping variables are coded appropriately.
DT$miscarriage_event <- as.integer(DT$miscarriage_event)
DT$mother_uid <- factor(DT$mother_uid)
DT$CLUSTID <- factor(DT$CLUSTID)

if (!"gest_month_f" %in% names(DT)) {
  DT$gest_month_f <- factor(DT$gest_month)
}

if (!"survey_year_f" %in% names(DT) && "survey_year" %in% names(DT)) {
  DT$survey_year_f <- factor(DT$survey_year)
}

# Remove unused levels before model construction.
DT <- droplevels(DT)


# ---- 3. Reconstruct lagged and cumulative WBGT if necessary -----------------

DT <- DT %>%
  arrange(preg_uid, gest_month)

if (!"WBT_lag1" %in% names(DT)) {
  DT <- DT %>%
    group_by(preg_uid) %>%
    mutate(WBT_lag1 = lag(WBT, 1L)) %>%
    ungroup()
}

if (!"WBT_prior_cummean" %in% names(DT)) {
  DT <- DT %>%
    group_by(preg_uid) %>%
    mutate(
      .wbt_running_sum = cumsum(replace(WBT, is.na(WBT), 0)),
      .wbt_running_n = cumsum(!is.na(WBT)),
      .wbt_running_mean = if_else(
        .wbt_running_n > 0,
        .wbt_running_sum / .wbt_running_n,
        NA_real_
      ),
      WBT_prior_cummean = lag(.wbt_running_mean, 1L)
    ) %>%
    ungroup() %>%
    select(-.wbt_running_sum, -.wbt_running_n, -.wbt_running_mean)
}

# A separate current-and-prior cumulative mean is also retained for descriptive
# checks, but the primary cumulative sensitivity uses prior months only to avoid
# including the concurrent exposure twice.
if (!"WBT_cummean" %in% names(DT)) {
  DT <- DT %>%
    group_by(preg_uid) %>%
    mutate(
      .wbt_running_sum = cumsum(replace(WBT, is.na(WBT), 0)),
      .wbt_running_n = cumsum(!is.na(WBT)),
      WBT_cummean = if_else(
        .wbt_running_n > 0,
        .wbt_running_sum / .wbt_running_n,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    select(-.wbt_running_sum, -.wbt_running_n)
}


# ---- 4. Final adjustment terms -----------------------------------------------

candidate_factor_terms <- c(
  "parity_cat",
  "previous_miscarriage_f",
  "interpregnancy_interval_f",
  "education_f",
  "wealth_f",
  "no_weekly_media_f",
  "health_insurance_f",
  "sanitation_f",
  "water_improved_f",
  "internet_ever_f"
)

available_factor_terms <- intersect(candidate_factor_terms, names(DT))

# Explicit exclusions required by the final manuscript specification.
available_factor_terms <- setdiff(
  available_factor_terms,
  c(
    "clean_cooking_fuel_f", "ph_cook_clean", "cook_fuel_type",
    "marital_status_f", "v501", "q703", "residence_f"
  )
)

for (variable in available_factor_terms) {
  DT[[variable]] <- droplevels(factor(DT[[variable]]))
}

# strata_f already contains survey year; use only one representation.
if ("strata_f" %in% names(DT) &&
    nlevels(droplevels(factor(DT$strata_f))) > 1L) {
  design_terms <- c("gest_month_f", "strata_f")
} else {
  design_terms <- c("gest_month_f", "survey_year_f")
}

# Residence is not entered separately when strata are included because the
# survey strata already encode geographic/residence design information.
if (!"strata_f" %in% design_terms && "residence_f" %in% names(DT)) {
  available_factor_terms <- c(available_factor_terms, "residence_f")
}

age_term <- if ("age_at_preg_start" %in% names(DT)) {
  "splines::ns(age_at_preg_start, df = 4)"
} else {
  character(0)
}

coexposure_spline_terms <- c(
  "splines::ns(CO, df = 4)",
  "splines::ns(PM, df = 4)",
  "splines::ns(EVI, df = 4)"
)

fixed_adjustment_terms <- c(
  design_terms,
  age_term,
  available_factor_terms,
  coexposure_spline_terms
)

random_terms <- c(
  "(1 | mother_uid)",
  "(1 | CLUSTID)"
)

formula_text <- function(outcome, terms) {
  stats::as.formula(
    paste(outcome, "~", paste(unique(terms), collapse = " + "))
  )
}


# ---- 5. General model utilities ----------------------------------------------

variables_in_formula <- function(formula, weight = NULL) {
  unique(c(all.vars(formula), weight))
}

prepare_model_data <- function(data, formula, weight = NULL) {
  variables <- variables_in_formula(formula, weight)
  variables <- intersect(variables, names(data))

  keep <- complete.cases(data[, variables, drop = FALSE])
  output <- droplevels(data[keep, , drop = FALSE])

  if (!is.null(weight)) {
    output <- output[
      is.finite(output[[weight]]) & output[[weight]] > 0,
      ,
      drop = FALSE
    ]
  }

  output
}

attach_analysis_counts <- function(model, data) {
  attr(model, "analysis_n") <- nrow(data)
  attr(model, "analysis_pregnancies") <- dplyr::n_distinct(data$preg_uid)
  attr(model, "analysis_events") <- sum(data$miscarriage_event == 1L)
  model
}

fit_or_load_glmmTMB <- function(
    formula,
    data,
    model_name,
    link = "cloglog",
    weight = "analysis_weight") {

  checkpoint_file <- file.path(
    checkpoint_dir,
    paste0(model_name, ".rds")
  )

  if (file.exists(checkpoint_file)) {
    saved_model <- readRDS(checkpoint_file)
    if (inherits(saved_model, "glmmTMB")) {
      message("Loading checkpoint: ", model_name)
      return(saved_model)
    }
    warning("Invalid checkpoint ignored: ", checkpoint_file)
  }

  model_data <- prepare_model_data(data, formula, weight)

  message(
    model_name, ": ", format(nrow(model_data), big.mark = ","),
    " records; ",
    format(n_distinct(model_data$preg_uid), big.mark = ","),
    " pregnancies; ",
    sum(model_data$miscarriage_event == 1L), " events."
  )

  model <- glmmTMB::glmmTMB(
    formula = formula,
    data = model_data,
    family = binomial(link = link),
    weights = if (is.null(weight)) NULL else model_data[[weight]],
    control = glmmTMBControl(
      optimizer = nlminb,
      optCtrl = list(iter.max = 20000, eval.max = 30000),
      rank_check = "adjust"
    )
  )

  model <- attach_analysis_counts(model, model_data)
  attr(model, "analysis_link") <- link
  attr(model, "analysis_weight_variable") <- ifelse(
    is.null(weight), "None", weight
  )

  saveRDS(model, checkpoint_file, compress = FALSE)
  message("Saved checkpoint: ", checkpoint_file)
  model
}

safe_logLik <- function(model) {
  tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_)
}

safe_information_criterion <- function(model, function_name) {
  tryCatch(
    as.numeric(get(function_name, mode = "function")(model)),
    error = function(e) NA_real_
  )
}

extract_model_fit <- function(model, model_name) {
  tibble(
    model = model_name,
    records = attr(model, "analysis_n") %||% NA_integer_,
    pregnancies = attr(model, "analysis_pregnancies") %||% NA_integer_,
    events = attr(model, "analysis_events") %||% NA_integer_,
    link = attr(model, "analysis_link") %||% NA_character_,
    weight_variable = attr(model, "analysis_weight_variable") %||%
      NA_character_,
    log_likelihood = safe_logLik(model),
    AIC = safe_information_criterion(model, "AIC"),
    BIC = safe_information_criterion(model, "BIC"),
    convergence_code = model$fit$convergence %||% NA_integer_,
    optimizer_message = model$fit$message %||% NA_character_,
    positive_definite_Hessian = model$sdr$pdHess %||% NA
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

extract_fixed_ratios <- function(model, model_name, effect_measure = "HR") {
  coefficient_matrix <- summary(model)$coefficients$cond

  tibble(
    model = model_name,
    effect_measure = effect_measure,
    term = rownames(coefficient_matrix),
    estimate_log_scale = coefficient_matrix[, "Estimate"],
    standard_error = coefficient_matrix[, "Std. Error"],
    estimate = exp(estimate_log_scale),
    lower_95 = exp(estimate_log_scale - 1.96 * standard_error),
    upper_95 = exp(estimate_log_scale + 1.96 * standard_error),
    p_value = coefficient_matrix[, "Pr(>|z|)"],
    estimate_95_CI = sprintf(
      "%.2f (%.2f–%.2f)", estimate, lower_95, upper_95
    ),
    p_value_formatted = if_else(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  )
}

extract_spline_wald <- function(model, model_name, exposure) {
  beta <- glmmTMB::fixef(model)$cond
  V <- as.matrix(vcov(model)$cond)

  pattern <- paste0(
    "^splines::ns\\(", exposure, ", df = [0-9]+\\)"
  )

  selected_terms <- grep(pattern, names(beta), value = TRUE)

  # Exclude interaction terms should this function be reused later.
  selected_terms <- selected_terms[!grepl(":", selected_terms)]

  if (length(selected_terms) == 0L) {
    return(tibble(
      model = model_name,
      exposure = exposure,
      test_df = NA_integer_,
      wald_chi_square = NA_real_,
      p_value = NA_real_
    ))
  }

  beta_selected <- beta[selected_terms]
  V_selected <- V[selected_terms, selected_terms, drop = FALSE]

  statistic <- tryCatch(
    as.numeric(t(beta_selected) %*% solve(V_selected, beta_selected)),
    error = function(e) NA_real_
  )

  tibble(
    model = model_name,
    exposure = exposure,
    test_df = length(selected_terms),
    wald_chi_square = statistic,
    p_value = pchisq(statistic, df = length(selected_terms), lower.tail = FALSE)
  )
}

# emmeans obtains the correct spline basis from the fitted model and propagates
# the full covariance of the two predictions. This avoids the earlier
# '.model_weight not found' prediction error.
contrast_original_unit <- function(
    model,
    data,
    exposure,
    model_name,
    effect_measure,
    increment = 1) {

  from_value <- stats::median(data[[exposure]], na.rm = TRUE)
  to_value <- from_value + increment

  at_values <- setNames(list(c(from_value, to_value)), exposure)

  reference_grid <- emmeans::emmeans(
    model,
    specs = exposure,
    at = at_values,
    type = "link"
  )

  comparison <- emmeans::contrast(
    reference_grid,
    method = list("increase" = c(-1, 1)),
    adjust = "none"
  )

  result <- as.data.frame(summary(comparison, infer = c(TRUE, TRUE)))

  tibble(
    model = model_name,
    exposure = exposure,
    contrast = paste0("+", increment, " original unit"),
    from_value = from_value,
    to_value = to_value,
    effect_measure = effect_measure,
    estimate_log_scale = result$estimate[1],
    standard_error = result$SE[1],
    estimate = exp(result$estimate[1]),
    lower_95 = exp(result$lower.CL[1]),
    upper_95 = exp(result$upper.CL[1]),
    p_value = result$p.value[1],
    estimate_95_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      exp(result$estimate[1]),
      exp(result$lower.CL[1]),
      exp(result$upper.CL[1])
    )
  )
}


# ---- 6. WBGT discrete-time sensitivity models -------------------------------

early_loss_ids <- DT %>%
  filter(miscarriage_event == 1L, gest_month - 1L < 2L) %>%
  distinct(preg_uid) %>%
  pull(preg_uid)

DT_exclude_early_loss <- DT %>%
  filter(!preg_uid %in% early_loss_ids)

if ("residence_eligible_strict" %in% names(DT)) {
  DT_strict_residence <- DT %>%
    filter(residence_eligible_strict %in% c(TRUE, 1L))
} else {
  warning(
    "residence_eligible_strict is absent; Strict_Residence_WBGT will use DT."
  )
  DT_strict_residence <- DT
}

scenario_definitions <- list(
  Primary_nonlinear_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 4)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Lag1_WBGT = list(
    exposure = "WBT_lag1",
    formula = formula_text(
      "miscarriage_event",
      c(
        "splines::ns(WBT_lag1, df = 4)",
        fixed_adjustment_terms,
        random_terms
      )
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  PriorCumulative_WBGT = list(
    exposure = "WBT_prior_cummean",
    formula = formula_text(
      "miscarriage_event",
      c(
        "splines::ns(WBT_prior_cummean, df = 4)",
        fixed_adjustment_terms,
        random_terms
      )
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Exclude_Loss_Before_2_Months_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 4)", fixed_adjustment_terms, random_terms)
    ),
    data = DT_exclude_early_loss,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Strict_Residence_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 4)", fixed_adjustment_terms, random_terms)
    ),
    data = DT_strict_residence,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Logit_Link_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 4)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "logit",
    weight = "analysis_weight",
    effect_measure = "Interval OR"
  ),
  Unweighted_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 4)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = NULL,
    effect_measure = "HR"
  ),
  Linear_Raw_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("WBT", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Spline_df3_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 3)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Spline_df5_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 5)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  ),
  Spline_df6_WBGT = list(
    exposure = "WBT",
    formula = formula_text(
      "miscarriage_event",
      c("splines::ns(WBT, df = 6)", fixed_adjustment_terms, random_terms)
    ),
    data = DT,
    link = "cloglog",
    weight = "analysis_weight",
    effect_measure = "HR"
  )
)

wbgt_sensitivity_models <- list()
wbgt_model_fit <- list()
wbgt_fixed_effects <- list()
wbgt_contrasts <- list()
wbgt_omnibus_tests <- list()

for (scenario_name in names(scenario_definitions)) {
  specification <- scenario_definitions[[scenario_name]]

  fitted_model <- fit_or_load_glmmTMB(
    formula = specification$formula,
    data = specification$data,
    model_name = scenario_name,
    link = specification$link,
    weight = specification$weight
  )

  wbgt_sensitivity_models[[scenario_name]] <- fitted_model
  wbgt_model_fit[[scenario_name]] <- extract_model_fit(
    fitted_model,
    scenario_name
  )
  wbgt_fixed_effects[[scenario_name]] <- extract_fixed_ratios(
    fitted_model,
    scenario_name,
    specification$effect_measure
  )

  wbgt_contrasts[[scenario_name]] <- tryCatch(
    contrast_original_unit(
      model = fitted_model,
      data = specification$data,
      exposure = specification$exposure,
      model_name = scenario_name,
      effect_measure = specification$effect_measure,
      increment = 1
    ),
    error = function(e) {
      warning("Contrast failed for ", scenario_name, ": ", conditionMessage(e))
      tibble(
        model = scenario_name,
        exposure = specification$exposure,
        contrast = "+1 original unit",
        from_value = NA_real_, to_value = NA_real_,
        effect_measure = specification$effect_measure,
        estimate_log_scale = NA_real_, standard_error = NA_real_,
        estimate = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_,
        p_value = NA_real_, estimate_95_CI = NA_character_
      )
    }
  )

  if (grepl("splines::ns", paste(deparse(specification$formula), collapse = ""))) {
    wbgt_omnibus_tests[[scenario_name]] <- extract_spline_wald(
      fitted_model,
      scenario_name,
      specification$exposure
    )
  }
}

wbgt_model_fit_df <- bind_rows(wbgt_model_fit)
wbgt_fixed_effects_df <- bind_rows(wbgt_fixed_effects)
wbgt_contrasts_df <- bind_rows(wbgt_contrasts)
wbgt_omnibus_tests_df <- bind_rows(wbgt_omnibus_tests)


# ---- 7. Population-averaged WBGT GEE ----------------------------------------

gee_formula <- formula_text(
  "miscarriage_event",
  c(
    "WBT",
    fixed_adjustment_terms
  )
)

gee_data <- prepare_model_data(DT, gee_formula, "analysis_weight") %>%
  arrange(mother_uid, preg_uid, gest_month)

message(
  "GEE_WBGT: ", format(nrow(gee_data), big.mark = ","), " records; ",
  format(n_distinct(gee_data$preg_uid), big.mark = ","), " pregnancies; ",
  sum(gee_data$miscarriage_event == 1L), " events."
)

gee_wbgt_model <- geepack::geeglm(
  formula = gee_formula,
  id = mother_uid,
  data = gee_data,
  family = binomial(link = "cloglog"),
  weights = analysis_weight,
  corstr = "exchangeable",
  std.err = "san.se",
  scale.fix = TRUE
)

gee_coefficient_matrix <- summary(gee_wbgt_model)$coefficients

gee_wbgt_fixed_effects <- tibble(
  model = "Population_averaged_GEE_WBGT",
  effect_measure = "HR",
  term = rownames(gee_coefficient_matrix),
  estimate_log_HR = gee_coefficient_matrix[, "Estimate"],
  robust_SE = gee_coefficient_matrix[, "Std.err"],
  HR = exp(estimate_log_HR),
  lower_95 = exp(estimate_log_HR - 1.96 * robust_SE),
  upper_95 = exp(estimate_log_HR + 1.96 * robust_SE),
  p_value = gee_coefficient_matrix[, "Pr(>|W|)"],
  estimate_95_CI = sprintf("%.2f (%.2f–%.2f)", HR, lower_95, upper_95),
  p_value_formatted = if_else(
    p_value < 0.001, "<0.001", sprintf("%.3f", p_value)
  )
)

gee_wbgt_summary <- gee_wbgt_fixed_effects %>%
  filter(term == "WBT")

saveRDS(
  gee_wbgt_model,
  file.path(models_dir, "Population_averaged_GEE_WBGT.rds"),
  compress = FALSE
)


# ---- 8. Pregnancy-level binary WBGT model -----------------------------------

first_nonmissing <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) NA else observed[1]
}

binary_required_variables <- unique(c(
  "preg_uid", "mother_uid", "CLUSTID", "miscarriage_event",
  "analysis_weight", "gest_month", "WBT", "CO", "PM", "EVI",
  design_terms[design_terms != "gest_month_f"],
  "age_at_preg_start", available_factor_terms
))
binary_required_variables <- intersect(binary_required_variables, names(DT))

if (!"outcome_type" %in% names(DT)) {
  stop("Binary sensitivity requires outcome_type to identify live birth and miscarriage.")
}
binary_data <- DT %>%
  filter(outcome_type %in% c("live_birth", "miscarriage")) %>%
  arrange(preg_uid, gest_month) %>%
  group_by(preg_uid) %>%
  summarise(
    miscarriage = max(miscarriage_event, na.rm = TRUE),
    mother_uid = first_nonmissing(mother_uid),
    CLUSTID = first_nonmissing(CLUSTID),
    analysis_weight = first_nonmissing(analysis_weight),
    WBT_early = mean(WBT[gest_month <= 3], na.rm = TRUE),
    CO_early = mean(CO[gest_month <= 3], na.rm = TRUE),
    PM_early = mean(PM[gest_month <= 3], na.rm = TRUE),
    EVI_early = mean(EVI[gest_month <= 3], na.rm = TRUE),
    across(
      all_of(setdiff(
        binary_required_variables,
        c(
          "preg_uid", "mother_uid", "CLUSTID", "miscarriage_event",
          "analysis_weight", "gest_month", "WBT", "CO", "PM", "EVI"
        )
      )),
      first_nonmissing
    ),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(WBT_early, CO_early, PM_early, EVI_early),
           ~ ifelse(is.nan(.x), NA_real_, .x)),
    mother_uid = factor(mother_uid),
    CLUSTID = factor(CLUSTID)
  ) %>%
  droplevels()

binary_design_terms <- intersect(
  c("survey_year_f", "strata_f"),
  names(binary_data)
)

binary_factor_terms <- intersect(available_factor_terms, names(binary_data))
binary_age_term <- if ("age_at_preg_start" %in% names(binary_data)) {
  "splines::ns(age_at_preg_start, df = 4)"
} else {
  character(0)
}

binary_formula <- formula_text(
  "miscarriage",
  c(
    "splines::ns(WBT_early, df = 4)",
    "splines::ns(CO_early, df = 4)",
    "splines::ns(PM_early, df = 4)",
    "splines::ns(EVI_early, df = 4)",
    binary_design_terms,
    binary_age_term,
    binary_factor_terms,
    "(1 | mother_uid)",
    "(1 | CLUSTID)"
  )
)

binary_checkpoint <- file.path(
  checkpoint_dir,
  "Pregnancy_level_binary_logit_WBGT.rds"
)

if (file.exists(binary_checkpoint) &&
    inherits(readRDS(binary_checkpoint), "glmmTMB")) {
  binary_wbgt_model <- readRDS(binary_checkpoint)
  message("Loading checkpoint: Pregnancy_level_binary_logit_WBGT")
} else {
  binary_model_data <- prepare_model_data(
    binary_data,
    binary_formula,
    "analysis_weight"
  )

  message(
    "Pregnancy_level_binary_logit_WBGT: ",
    format(nrow(binary_model_data), big.mark = ","), " pregnancies; ",
    sum(binary_model_data$miscarriage == 1L), " miscarriages."
  )

  binary_wbgt_model <- glmmTMB(
    formula = binary_formula,
    data = binary_model_data,
    family = binomial(link = "logit"),
    weights = analysis_weight,
    control = glmmTMBControl(
      optimizer = nlminb,
      optCtrl = list(iter.max = 20000, eval.max = 30000),
      rank_check = "adjust"
    )
  )

  binary_wbgt_model <- attach_analysis_counts(
    binary_wbgt_model,
    binary_model_data %>%
      mutate(
        preg_uid = row_number(),
        miscarriage_event = miscarriage
      )
  )
  attr(binary_wbgt_model, "analysis_link") <- "logit"
  attr(binary_wbgt_model, "analysis_weight_variable") <- "analysis_weight"

  saveRDS(binary_wbgt_model, binary_checkpoint, compress = FALSE)
}

binary_wbgt_fixed_effects <- extract_fixed_ratios(
  binary_wbgt_model,
  "Pregnancy_level_binary_logit_WBGT",
  "OR"
)

binary_wbgt_omnibus <- bind_rows(
  extract_spline_wald(
    binary_wbgt_model,
    "Pregnancy_level_binary_logit_WBGT",
    "WBT_early"
  ),
  extract_spline_wald(
    binary_wbgt_model,
    "Pregnancy_level_binary_logit_WBGT",
    "CO_early"
  ),
  extract_spline_wald(
    binary_wbgt_model,
    "Pregnancy_level_binary_logit_WBGT",
    "PM_early"
  ),
  extract_spline_wald(
    binary_wbgt_model,
    "Pregnancy_level_binary_logit_WBGT",
    "EVI_early"
  )
)

binary_wbgt_contrast <- contrast_original_unit(
  model = binary_wbgt_model,
  data = binary_data,
  exposure = "WBT_early",
  model_name = "Pregnancy_level_binary_logit_WBGT",
  effect_measure = "OR",
  increment = 1
)


# ---- 9. Export all results ---------------------------------------------------

write_csv(
  wbgt_model_fit_df,
  file.path(tables_dir, "WBGT_sensitivity_model_fit.csv")
)
write_csv(
  wbgt_fixed_effects_df,
  file.path(tables_dir, "WBGT_sensitivity_fixed_effects.csv")
)
write_csv(
  wbgt_contrasts_df,
  file.path(tables_dir, "WBGT_sensitivity_local_contrasts.csv")
)
write_csv(
  wbgt_omnibus_tests_df,
  file.path(tables_dir, "WBGT_sensitivity_spline_omnibus_tests.csv")
)
write_csv(
  gee_wbgt_fixed_effects,
  file.path(tables_dir, "Sensitivity_GEE_WBGT_fixed_effects.csv")
)
write_csv(
  gee_wbgt_summary,
  file.path(tables_dir, "Sensitivity_GEE_WBGT_exposure_result.csv")
)
write_csv(
  binary_wbgt_fixed_effects,
  file.path(tables_dir, "Sensitivity_binary_WBGT_fixed_effects_OR.csv")
)
write_csv(
  binary_wbgt_omnibus,
  file.path(tables_dir, "Sensitivity_binary_WBGT_omnibus_tests.csv")
)
write_csv(
  binary_wbgt_contrast,
  file.path(tables_dir, "Sensitivity_binary_WBGT_local_contrast.csv")
)

saveRDS(
  wbgt_sensitivity_models,
  file.path(models_dir, "WBGT_sensitivity_models_complete.rds"),
  compress = FALSE
)

saveRDS(
  binary_wbgt_model,
  file.path(models_dir, "Pregnancy_level_binary_logit_WBGT.rds"),
  compress = FALSE
)

# A single workbook for supplementary-table preparation.
workbook_file <- file.path(
  tables_dir,
  "WBGT_sensitivity_results_for_supplement.xlsx"
)

openxlsx::write.xlsx(
  list(
    DTS_local_contrasts = wbgt_contrasts_df,
    DTS_omnibus_tests = wbgt_omnibus_tests_df,
    DTS_model_fit = wbgt_model_fit_df,
    DTS_fixed_effects = wbgt_fixed_effects_df,
    GEE_exposure_result = gee_wbgt_summary,
    GEE_fixed_effects = gee_wbgt_fixed_effects,
    Binary_local_contrast = binary_wbgt_contrast,
    Binary_omnibus_tests = binary_wbgt_omnibus,
    Binary_fixed_effects = binary_wbgt_fixed_effects
  ),
  file = workbook_file,
  overwrite = TRUE
)


# ---- 10. Essential diagnostics ----------------------------------------------

diagnostic_summary <- imap_dfr(
  wbgt_sensitivity_models,
  function(model, model_name) {
    gradient <- tryCatch(
      max(abs(model$sdr$gradient.fixed), na.rm = TRUE),
      error = function(e) NA_real_
    )

    tibble(
      model = model_name,
      convergence_code = model$fit$convergence %||% NA_integer_,
      optimizer_message = model$fit$message %||% NA_character_,
      positive_definite_Hessian = model$sdr$pdHess %||% NA,
      maximum_absolute_gradient = gradient
    )
  }
)

binary_diagnostic_summary <- tibble(
  model = "Pregnancy_level_binary_logit_WBGT",
  convergence_code = binary_wbgt_model$fit$convergence %||% NA_integer_,
  optimizer_message = binary_wbgt_model$fit$message %||% NA_character_,
  positive_definite_Hessian = binary_wbgt_model$sdr$pdHess %||% NA,
  maximum_absolute_gradient = tryCatch(
    max(abs(binary_wbgt_model$sdr$gradient.fixed), na.rm = TRUE),
    error = function(e) NA_real_
  )
)

diagnostic_summary <- bind_rows(
  diagnostic_summary,
  binary_diagnostic_summary
)

write_csv(
  diagnostic_summary,
  file.path(tables_dir, "WBGT_sensitivity_model_diagnostics.csv")
)

print(diagnostic_summary, n = Inf, width = Inf)
print(wbgt_contrasts_df, n = Inf, width = Inf)
print(wbgt_omnibus_tests_df, n = Inf, width = Inf)
print(gee_wbgt_summary, width = Inf)
print(binary_wbgt_contrast, width = Inf)
print(binary_wbgt_omnibus, n = Inf, width = Inf)

message("WBGT sensitivity analyses completed.")
message("Results workbook: ", workbook_file)
