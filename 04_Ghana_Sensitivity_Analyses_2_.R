# GitHub copy: portable paths; restricted data are not distributed.
# ============================================================================
# 04_Ghana_Sensitivity_Analyses.R
# Sensitivity, robustness, and effect-modification analyses
#
# Final implementation aligned with Script 02:
#   - sparse glmmTMB mother and DHS-cluster random intercepts
#   - natural cubic splines for nonlinear exposure-response functions
#   - complementary log-log link for discrete-time hazard models
#   - raw exposure units; linear EVI estimates are per 0.1-unit increase
#   - no region random intercept, current marital status, or cooking fuel
# ============================================================================

rm(list = ls())
gc(reset = TRUE)


# ---- 0. Packages ------------------------------------------------------------

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "splines",
  "glmmTMB", "geepack", "ggplot2", "patchwork", "MASS"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install before running Script 04: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(splines)
  library(glmmTMB)
  library(geepack)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260928)


# ---- 1. Open OnDemand paths and settings -----------------------------------

# Run from the repository root, or set GHANA_MISCARRIAGE_ROOT explicitly.
# Private inputs and outputs can live entirely outside the Git checkout.
project_dir <- normalizePath(
  Sys.getenv("GHANA_MISCARRIAGE_ROOT", unset = "."),
  winslash = "/", mustWork = TRUE
)

results_dir <- file.path(project_dir, "Results")
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")
models_dir <- file.path(results_dir, "Models")
checkpoint_dir <- file.path(models_dir, "Sensitivity_glmmTMB_checkpoints_no_cooking_v2")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

input_all <- file.path(
  models_dir, "Ghana_DT_AnalysisReady_AllRecords.rds"
)
input_primary <- file.path(
  models_dir, "Ghana_DT_AnalysisReady_PrimaryCC.rds"
)
input_primary_models <- file.path(
  models_dir, "Primary_discrete_time_spline_GLMM_models.rds"
)

thread_text <- Sys.getenv("GHANA_N_THREADS", unset = "12")
N_THREADS <- suppressWarnings(as.integer(thread_text))
if (is.na(N_THREADS) || N_THREADS < 1L) N_THREADS <- 12L

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.finite(detected_cores) && detected_cores > 0L) {
  N_THREADS <- min(N_THREADS, as.integer(detected_cores))
}

refit_text <- tolower(Sys.getenv("GHANA_REFIT_SENSITIVITY", unset = "false"))
REFIT_MODELS <- refit_text %in% c("true", "1", "yes", "y")

options(mc.cores = N_THREADS)

message(
  "Project root: ", project_dir,
  "\nExecution: Open OnDemand interactive RStudio",
  "\nglmmTMB threads: ", N_THREADS,
  "\nRefit checkpoints: ", REFIT_MODELS
)


# ---- 2. Helpers -------------------------------------------------------------

formula_text <- function(outcome, terms) {
  terms <- unique(terms[!is.na(terms) & nzchar(terms)])
  as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

model_variables <- function(formula) unique(all.vars(formula))

ns_term <- function(variable, df = 4L) {
  paste0("splines::ns(", variable, ", df = ", as.integer(df), ")")
}

null_default <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  names(which.max(table(x)))[1]
}

save_high_resolution <- function(plot, stem, width = 11, height = 7) {
  ggsave(
    file.path(figures_dir, paste0(stem, ".png")),
    plot = plot, width = width, height = height,
    units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figures_dir, paste0(stem, ".tiff")),
    plot = plot, width = width, height = height,
    units = "in", dpi = 600, compression = "lzw", bg = "white"
  )
}

fit_or_load_glmmTMB <- function(
    formula,
    data,
    model_name,
    link = "cloglog",
    weight_variable = "analysis_weight"
) {
  outcome <- all.vars(formula)[1]
  needed <- unique(c(model_variables(formula), weight_variable))
  missing_needed <- setdiff(needed, names(data))

  if (length(missing_needed) > 0L) {
    stop(
      model_name, " is missing: ",
      paste(missing_needed, collapse = ", ")
    )
  }

  checkpoint_file <- file.path(checkpoint_dir, paste0(model_name, ".rds"))

  if (file.exists(checkpoint_file) && !REFIT_MODELS) {
    message("Loading checkpoint: ", model_name)
    model <- readRDS(checkpoint_file)

    if (!inherits(model, "glmmTMB")) {
      stop("Invalid checkpoint: ", checkpoint_file)
    }

    # up2date() has returned diagnostic character values in some glmmTMB
    # versions. Never overwrite a valid fitted model unless the returned
    # object is itself a glmmTMB model.
    updated_model <- try(
      glmmTMB::up2date(model),
      silent = TRUE
    )
    if (inherits(updated_model, "glmmTMB")) {
      model <- updated_model
    }

    return(model)
  }

  dat <- data[
    complete.cases(data[, needed, drop = FALSE]),
    ,
    drop = FALSE
  ]
  dat <- droplevels(dat)
  dat$.model_weight <- as.numeric(dat[[weight_variable]])

  if (nrow(dat) == 0L) stop(model_name, " has no complete records.")
  if (!all(dat[[outcome]] %in% c(0, 1))) {
    stop(model_name, " outcome is not binary 0/1.")
  }
  if (sum(dat[[outcome]], na.rm = TRUE) == 0L) {
    stop(model_name, " has no events.")
  }
  if (
    any(!is.finite(dat$.model_weight)) ||
    any(dat$.model_weight <= 0)
  ) {
    stop(model_name, " contains invalid model weights.")
  }

  message(
    model_name, ": ", format(nrow(dat), big.mark = ","), " records; ",
    format(n_distinct(dat$preg_uid), big.mark = ","), " pregnancies; ",
    sum(dat[[outcome]], na.rm = TRUE), " events."
  )

  model <- glmmTMB::glmmTMB(
    formula = formula,
    family = binomial(link = link),
    data = dat,
    weights = .model_weight,
    na.action = na.exclude,
    control = glmmTMB::glmmTMBControl(
      parallel = N_THREADS,
      optCtrl = list(iter.max = 1000L, eval.max = 1000L)
    )
  )

  attr(model, "model_name") <- model_name
  attr(model, "analysis_n") <- nrow(dat)
  attr(model, "analysis_pregnancies") <- n_distinct(dat$preg_uid)
  attr(model, "analysis_events") <- sum(dat[[outcome]], na.rm = TRUE)
  attr(model, "analysis_link") <- link
  attr(model, "weight_variable") <- weight_variable

  if (!is.null(model$fit$convergence) && model$fit$convergence != 0L) {
    warning(model_name, " convergence code: ", model$fit$convergence)
  }
  if (!is.null(model$sdr$pdHess) && !isTRUE(model$sdr$pdHess)) {
    warning(model_name, " has a non-positive-definite Hessian.")
  }

  saveRDS(model, checkpoint_file, compress = FALSE)
  message("Saved checkpoint: ", checkpoint_file)
  model
}

extract_model_fit <- function(model, model_name) {
  ll <- suppressWarnings(
    tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_)
  )

  analysis_n <- attr(model, "analysis_n")
  if (is.null(analysis_n) || length(analysis_n) == 0L) {
    analysis_n <- if (!is.null(model$frame)) {
      nrow(model$frame)
    } else {
      NA_integer_
    }
  }

  analysis_pregnancies <- attr(model, "analysis_pregnancies")
  if (is.null(analysis_pregnancies) || length(analysis_pregnancies) == 0L) {
    analysis_pregnancies <- NA_integer_
  }

  analysis_events <- attr(model, "analysis_events")
  if (is.null(analysis_events) || length(analysis_events) == 0L) {
    analysis_events <- NA_integer_
  }

  analysis_link <- attr(model, "analysis_link")
  if (is.null(analysis_link) || length(analysis_link) == 0L) {
    analysis_link <- model$modelInfo$family$link
  }

  weight_variable <- attr(model, "weight_variable")
  if (is.null(weight_variable) || length(weight_variable) == 0L) {
    weight_variable <- NA_character_
  }

  tibble(
    model = model_name,
    records = as.integer(analysis_n),
    pregnancies = as.integer(analysis_pregnancies),
    events = as.integer(analysis_events),
    link = as.character(analysis_link),
    weight_variable = as.character(weight_variable),
    log_likelihood = ll,
    AIC = if (is.finite(ll)) AIC(model) else NA_real_,
    BIC = if (is.finite(ll)) BIC(model) else NA_real_,
    convergence_code = model$fit$convergence,
    optimizer_message = model$fit$message,
    positive_definite_Hessian = model$sdr$pdHess
  )
}

extract_fixed_ratios <- function(
    model,
    model_name,
    effect_measure = "HR"
) {
  tab <- as.data.frame(summary(model)$coefficients$cond)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  estimate <- tab[["Estimate"]]
  standard_error <- tab[["Std. Error"]]
  p_value <- tab[["Pr(>|z|)"]]

  tibble(
    model = model_name,
    effect_measure = effect_measure,
    term = tab$term,
    estimate_log_scale = estimate,
    standard_error = standard_error,
    ratio = exp(estimate),
    lower_95 = exp(estimate - 1.96 * standard_error),
    upper_95 = exp(estimate + 1.96 * standard_error),
    p_value = p_value,
    p_value_formatted = format_p(p_value)
  )
}

# Local raw-unit contrast using the full covariance between two fixed-effect
# predictions. For cloglog models, exponentiating the link-scale contrast gives
# the conditional monthly hazard ratio.
contrast_raw_unit <- function(
    model,
    data,
    exposure,
    model_name,
    effect_measure = "HR",
    unit_increase = 1,
    modifier = NULL,
    category = NULL
) {
  # Reconstruct the exact complete-case analysis sample used by the fitted
  # model. This is essential for lagged/prior-cumulative specifications, where
  # early gestational months can be absent and their factor levels were dropped
  # before fitting.
  weight_variable <- attr(model, "weight_variable")
  if (is.null(weight_variable) || length(weight_variable) == 0L) {
    weight_variable <- character(0)
  }

  needed <- unique(c(
    model_variables(formula(model)),
    weight_variable
  ))
  missing_needed <- setdiff(needed, names(data))
  if (length(missing_needed) > 0L) {
    stop(
      model_name, " contrast data are missing: ",
      paste(missing_needed, collapse = ", ")
    )
  }

  analysis_data <- data[
    complete.cases(data[, needed, drop = FALSE]),
    ,
    drop = FALSE
  ]
  analysis_data <- droplevels(analysis_data)

  if (nrow(analysis_data) == 0L) {
    stop(model_name, " has no complete records for prediction.")
  }

  observed <- as.numeric(analysis_data[[exposure]])
  observed <- observed[is.finite(observed)]
  if (length(observed) == 0L) {
    stop("No finite values for ", exposure, ".")
  }

  from_value <- median(observed)
  to_value <- from_value + unit_increase

  if (to_value > max(observed)) {
    to_value <- from_value
    from_value <- from_value - unit_increase
  }
  if (from_value < min(observed)) {
    stop("Range of ", exposure, " is smaller than the contrast.")
  }

  reference_row <- analysis_data[1, , drop = FALSE]

  for (variable in names(reference_row)) {
    x <- analysis_data[[variable]]
    if (is.factor(x)) {
      reference_row[[variable]] <- factor(
        mode_value(x), levels = levels(x)
      )
    } else if (is.numeric(x) && variable != "miscarriage_event") {
      reference_row[[variable]] <- median(x, na.rm = TRUE)
    }
  }

  if (!is.null(modifier)) {
    reference_row[[modifier]] <- factor(
      category, levels = levels(analysis_data[[modifier]])
    )
  }
  prediction_data <- reference_row[rep(1, 2), , drop = FALSE]
  prediction_data[[exposure]] <- c(from_value, to_value)
  # The fitted call evaluates weights = .model_weight during prediction.
  prediction_data$.model_weight <- 1

  prediction <- predict(
    model,
    newdata = prediction_data,
    type = "link",
    re.form = NA,
    allow.new.levels = TRUE,
    se.fit = TRUE,
    cov.fit = TRUE
  )

  fitted_link <- as.numeric(prediction$fit)
  prediction_covariance <- as.matrix(prediction$cov.fit)
  estimate_log_scale <- fitted_link[2] - fitted_link[1]
  contrast_variance <-
    prediction_covariance[2, 2] +
    prediction_covariance[1, 1] -
    2 * prediction_covariance[2, 1]
  standard_error <- sqrt(max(0, contrast_variance))
  z_value <- estimate_log_scale / standard_error
  p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)

  tibble(
    model = model_name,
    exposure = exposure,
    contrast = paste0("+", unit_increase, " original unit"),
    from_value = from_value,
    to_value = to_value,
    effect_measure = effect_measure,
    estimate_log_scale = estimate_log_scale,
    standard_error = standard_error,
    estimate = exp(estimate_log_scale),
    lower_95 = exp(estimate_log_scale - 1.96 * standard_error),
    upper_95 = exp(estimate_log_scale + 1.96 * standard_error),
    p_value = p_value
  )
}

joint_wald_test <- function(model, pattern, model_name, modifier) {
  beta_all <- glmmTMB::fixef(model)$cond
  covariance_all <- as.matrix(vcov(model)$cond)
  index <- grep(pattern, names(beta_all))

  if (length(index) == 0L) {
    return(tibble(
      model = model_name,
      modifier = modifier,
      degrees_freedom = 0L,
      Wald_chi_square = NA_real_,
      p_value = NA_real_
    ))
  }

  beta <- beta_all[index]
  covariance <- covariance_all[index, index, drop = FALSE]
  inverse_covariance <- tryCatch(
    solve(covariance),
    error = function(e) MASS::ginv(covariance)
  )
  statistic <- as.numeric(t(beta) %*% inverse_covariance %*% beta)

  tibble(
    model = model_name,
    modifier = modifier,
    degrees_freedom = length(index),
    Wald_chi_square = statistic,
    p_value = pchisq(
      statistic, df = length(index), lower.tail = FALSE
    )
  )
}


# ---- 3. Load finalized Script 02 outputs -----------------------------------

required_inputs <- c(input_all, input_primary, input_primary_models)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
  stop(
    "Missing Script 02 output(s): ",
    paste(missing_inputs, collapse = ", ")
  )
}

DT <- readRDS(input_all)
DT_primary <- readRDS(input_primary)
primary_models <- readRDS(input_primary_models)

# Formula changes cannot recover pregnancies excluded by the old fuel filter.
# Require both datasets to have been rebuilt by the corrected Script 02.
if (!isTRUE(attr(DT, "ghana_cooking_fuel_excluded")) ||
    !isTRUE(attr(DT_primary, "ghana_cooking_fuel_excluded"))) {
  stop(
    "Rerun the cooking-fuel-corrected Script 02 before Script 04. ",
    "The existing analysis-ready datasets may retain the old fuel restriction."
  )
}
if (any(vapply(primary_models, function(model) {
  any(grepl("cook|fuel", all.vars(formula(model)), ignore.case = TRUE))
}, logical(1)))) {
  stop("Primary spline models still contain cooking fuel; rerun corrected Script 02.")
}

required_primary_models <- c(
  "Spline_M3_multi_HI", "Spline_M3_multi_WBT"
)
missing_primary_models <- setdiff(
  required_primary_models, names(primary_models)
)

if (length(missing_primary_models) > 0L) {
  stop(
    "Primary model file is missing: ",
    paste(missing_primary_models, collapse = ", ")
  )
}

DT <- DT %>%
  mutate(EVI_per_0_1 = EVI / 0.1, unit_weight = 1)

DT_primary <- DT_primary %>%
  mutate(EVI_per_0_1 = EVI / 0.1, unit_weight = 1)

required_core <- c(
  "miscarriage_event", "preg_uid", "mother_uid", "CLUSTID",
  "analysis_weight", "gest_month_f", "survey_year_f",
  "HI", "WBT", "CO", "PM", "EVI"
)
missing_core <- setdiff(required_core, names(DT_primary))

if (length(missing_core) > 0L) {
  stop(
    "Primary dataset is missing: ",
    paste(missing_core, collapse = ", ")
  )
}

# Current marital status and cooking fuel are intentionally excluded.
factor_candidates <- c(
  "education_f", "wealth_f", "no_weekly_media_f",
  "health_insurance_f", "sanitation_f", "water_improved_f",
  "internet_ever_f", "residence_f"
)

available_factor_terms <- intersect(
  factor_candidates, names(DT_primary)
)
available_factor_terms <- available_factor_terms[
  vapply(
    DT_primary[available_factor_terms],
    function(x) nlevels(droplevels(factor(x))) > 1L,
    logical(1)
  )
]

# strata_f was constructed in Script 02 as survey year x source stratum.
if (
  "strata_f" %in% names(DT_primary) &&
  nlevels(droplevels(DT_primary$strata_f)) > 1L
) {
  design_terms <- c("gest_month_f", "strata_f")
  available_factor_terms <- setdiff(
    available_factor_terms, "residence_f"
  )
} else {
  design_terms <- c("gest_month_f", "survey_year_f")
}

reproductive_terms <- c(
  ns_term("age_at_preg_start", 4L),
  "parity_cat",
  "previous_miscarriage_f",
  "interpregnancy_interval_f"
)

random_terms <- c("(1 | mother_uid)", "(1 | CLUSTID)")

adjustment_fixed_terms <- unique(c(
  design_terms, reproductive_terms, available_factor_terms
))
adjustment_terms <- unique(c(
  adjustment_fixed_terms, random_terms
))

# Residence-interaction models use survey year and residence explicitly.
residence_fixed_terms <- unique(c(
  "gest_month_f",
  "survey_year_f",
  reproductive_terms,
  setdiff(available_factor_terms, "residence_f"),
  if ("residence_f" %in% names(DT_primary)) "residence_f" else NULL
))
residence_adjustment_terms <- unique(c(
  residence_fixed_terms, random_terms
))

write_csv(
  tibble(
    component = c(
      "Design terms", "Reproductive terms", "Other adjustments",
      "Random effects", "Excluded covariate", "Linear EVI scale"
    ),
    specification = c(
      paste(design_terms, collapse = " + "),
      paste(reproductive_terms, collapse = " + "),
      paste(available_factor_terms, collapse = " + "),
      paste(random_terms, collapse = " + "),
      "Current marital status; cooking fuel",
      "Per 0.1 EVI unit"
    )
  ),
  file.path(tables_dir, "Sensitivity_model_specification.csv")
)


# ---- 4. Formula builder -----------------------------------------------------

build_environmental_formula <- function(
    thermal = "HI",
    co = "CO",
    pm = "PM",
    evi = "EVI",
    spline_df = 4L,
    nonlinear = TRUE,
    covariate_terms = adjustment_terms,
    outcome = "miscarriage_event"
) {
  exposures <- c(thermal, co, pm, evi)

  environmental_terms <- if (nonlinear) {
    vapply(
      exposures, ns_term, character(1), df = spline_df
    )
  } else {
    exposures
  }

  formula_text(
    outcome,
    c(environmental_terms, covariate_terms)
  )
}


# ---- 5. Primary-model local contrasts --------------------------------------

sensitivity_contrasts <- list()

sensitivity_contrasts$Primary_HI <- contrast_raw_unit(
  primary_models$Spline_M3_multi_HI,
  DT_primary,
  "HI",
  "Primary nonlinear HI",
  "HR",
  1
)

sensitivity_contrasts$Primary_WBGT <- contrast_raw_unit(
  primary_models$Spline_M3_multi_WBT,
  DT_primary,
  "WBT",
  "Primary nonlinear WBGT",
  "HR",
  1
)


# ---- 6. Principal discrete-time sensitivity models -------------------------

scenario_definitions <- list()

lag_variables <- c(
  "HI_lag1", "WBT_lag1", "CO_lag1", "PM_lag1", "EVI_lag1"
)

if (all(lag_variables %in% names(DT))) {
  scenario_definitions$Lag1_HI <- list(
    data = DT,
    formula = build_environmental_formula(
      "HI_lag1", "CO_lag1", "PM_lag1", "EVI_lag1"
    ),
    exposure = "HI_lag1",
    link = "cloglog",
    weight = "analysis_weight"
  )

  scenario_definitions$Lag1_WBGT <- list(
    data = DT,
    formula = build_environmental_formula(
      "WBT_lag1", "CO_lag1", "PM_lag1", "EVI_lag1"
    ),
    exposure = "WBT_lag1",
    link = "cloglog",
    weight = "analysis_weight"
  )
}

prior_variables <- c(
  "HI_prior_cummean", "CO_prior_cummean",
  "PM_prior_cummean", "EVI_prior_cummean"
)

if (all(prior_variables %in% names(DT))) {
  scenario_definitions$PriorCumulative_HI <- list(
    data = DT,
    formula = build_environmental_formula(
      "HI_prior_cummean", "CO_prior_cummean",
      "PM_prior_cummean", "EVI_prior_cummean"
    ),
    exposure = "HI_prior_cummean",
    link = "cloglog",
    weight = "analysis_weight"
  )
}

if (all(c("miscarriage_event", "gest_month") %in% names(DT))) {
  excluded_early_ids <- DT %>%
    filter(
      miscarriage_event == 1L,
      gest_month - 1L < 2L
    ) %>%
    distinct(preg_uid) %>%
    pull(preg_uid)

  scenario_definitions$Exclude_Loss_Before_2_Months <- list(
    data = DT_primary %>%
      filter(!preg_uid %in% excluded_early_ids),
    formula = build_environmental_formula(),
    exposure = "HI",
    link = "cloglog",
    weight = "analysis_weight"
  )
}

if ("residence_eligible_strict" %in% names(DT_primary)) {
  scenario_definitions$Strict_Residence <- list(
    data = DT_primary %>%
      filter(residence_eligible_strict %in% TRUE),
    formula = build_environmental_formula(),
    exposure = "HI",
    link = "cloglog",
    weight = "analysis_weight"
  )
}

scenario_definitions$Logit_Link <- list(
  data = DT_primary,
  formula = build_environmental_formula(),
  exposure = "HI",
  link = "logit",
  weight = "analysis_weight"
)

scenario_definitions$Unweighted <- list(
  data = DT_primary,
  formula = build_environmental_formula(),
  exposure = "HI",
  link = "cloglog",
  weight = "unit_weight"
)

# Raw linear terms; EVI is per 0.1-unit increase.
scenario_definitions$Linear_Raw_HI <- list(
  data = DT_primary,
  formula = build_environmental_formula(
    "HI", "CO", "PM", "EVI_per_0_1",
    nonlinear = FALSE
  ),
  exposure = "HI",
  link = "cloglog",
  weight = "analysis_weight"
)

# Primary nonlinear models use df=4.
for (df_value in c(3L, 5L, 6L)) {
  scenario_name <- paste0("Spline_df", df_value, "_HI")

  scenario_definitions[[scenario_name]] <- list(
    data = DT_primary,
    formula = build_environmental_formula(
      spline_df = df_value
    ),
    exposure = "HI",
    link = "cloglog",
    weight = "analysis_weight"
  )
}

sensitivity_models <- list()
sensitivity_model_fit <- list()
sensitivity_fixed_effects <- list()

for (scenario_name in names(scenario_definitions)) {
  specification <- scenario_definitions[[scenario_name]]

  fitted_model <- fit_or_load_glmmTMB(
    specification$formula,
    specification$data,
    scenario_name,
    specification$link,
    specification$weight
  )

  sensitivity_models[[scenario_name]] <- fitted_model
  sensitivity_model_fit[[scenario_name]] <- extract_model_fit(
    fitted_model, scenario_name
  )

  effect_measure <- if (
    identical(specification$link, "logit")
  ) "Interval OR" else "HR"

  sensitivity_fixed_effects[[scenario_name]] <- extract_fixed_ratios(
    fitted_model, scenario_name, effect_measure
  )

  sensitivity_contrasts[[scenario_name]] <- contrast_raw_unit(
    fitted_model,
    specification$data,
    specification$exposure,
    scenario_name,
    effect_measure,
    1
  )
}

write_csv(
  bind_rows(sensitivity_model_fit),
  file.path(tables_dir, "Sensitivity_model_fit.csv")
)
write_csv(
  bind_rows(sensitivity_fixed_effects),
  file.path(tables_dir, "Sensitivity_fixed_effects.csv")
)
write_csv(
  bind_rows(sensitivity_contrasts),
  file.path(tables_dir, "Sensitivity_exposure_contrasts.csv")
)
saveRDS(
  sensitivity_models,
  file.path(models_dir, "Sensitivity_discrete_time_GLMM_models.rds"),
  compress = FALSE
)


# ---- 7. Population-averaged GEE --------------------------------------------

gee_terms <- unique(c(
  "HI", "CO", "PM", "EVI_per_0_1",
  design_terms,
  ns_term("age_at_preg_start", 4L),
  "parity_cat", "previous_miscarriage_f",
  "interpregnancy_interval_f",
  available_factor_terms
))

gee_formula <- formula_text("miscarriage_event", gee_terms)
gee_needed <- unique(c(
  model_variables(gee_formula),
  "mother_uid", "preg_uid", "gest_month", "analysis_weight"
))

gee_data <- DT_primary[
  complete.cases(DT_primary[, gee_needed, drop = FALSE]),
  ,
  drop = FALSE
] %>%
  droplevels() %>%
  arrange(mother_uid, preg_uid, gest_month)

gee_model <- geepack::geeglm(
  formula = gee_formula,
  family = binomial(link = "cloglog"),
  data = gee_data,
  id = mother_uid,
  weights = analysis_weight,
  corstr = "exchangeable",
  std.err = "san.se"
)

gee_tab <- as.data.frame(summary(gee_model)$coefficients)
gee_tab$term <- rownames(gee_tab)
rownames(gee_tab) <- NULL
gee_p <- gee_tab[["Pr(>|W|)"]]

gee_results <- tibble(
  model = "Population-averaged GEE",
  term = gee_tab$term,
  estimate_log_HR = gee_tab$Estimate,
  robust_SE = gee_tab$Std.err,
  HR = exp(gee_tab$Estimate),
  lower_95 = exp(gee_tab$Estimate - 1.96 * gee_tab$Std.err),
  upper_95 = exp(gee_tab$Estimate + 1.96 * gee_tab$Std.err),
  p_value = gee_p,
  p_value_formatted = format_p(gee_p)
)

write_csv(
  gee_results,
  file.path(tables_dir, "Sensitivity_GEE_HRs.csv")
)
saveRDS(
  gee_model,
  file.path(models_dir, "Sensitivity_GEE_model.rds"),
  compress = FALSE
)


# ---- 8. Secondary pregnancy-level binary analysis --------------------------

if ("outcome_type" %in% names(DT)) {
  binary_factor_terms <- intersect(
    available_factor_terms, names(DT)
  )

  binary_data <- DT_primary %>%
    filter(outcome_type %in% c("live_birth", "miscarriage")) %>%
    arrange(preg_uid, gest_month) %>%
    group_by(preg_uid) %>%
    summarise(
      miscarriage = as.integer(first(outcome_type) == "miscarriage"),
      mother_uid = first(mother_uid),
      CLUSTID = first(CLUSTID),
      survey_year_f = first(survey_year_f),
      analysis_weight = first(analysis_weight),
      age_at_preg_start = first(age_at_preg_start),
      parity_cat = first(parity_cat),
      previous_miscarriage_f = first(previous_miscarriage_f),
      interpregnancy_interval_f = first(interpregnancy_interval_f),
      across(all_of(binary_factor_terms), first),
      HI_early = mean(HI[gest_month <= 3], na.rm = TRUE),
      CO_early = mean(CO[gest_month <= 3], na.rm = TRUE),
      PM_early = mean(PM[gest_month <= 3], na.rm = TRUE),
      EVI_early = mean(EVI[gest_month <= 3], na.rm = TRUE),
      .groups = "drop"
    )

  if ("strata_f" %in% names(DT)) {
    pregnancy_strata <- DT %>%
      distinct(preg_uid, strata_f)
    binary_data <- binary_data %>%
      left_join(pregnancy_strata, by = "preg_uid")
  }

  binary_data <- binary_data %>%
    mutate(
      EVI_early_per_0_1 = EVI_early / 0.1,
      mother_uid = factor(mother_uid),
      CLUSTID = factor(CLUSTID)
    )

  if (
    "strata_f" %in% names(binary_data) &&
    nlevels(droplevels(factor(binary_data$strata_f))) > 1L
  ) {
    binary_design_terms <- "strata_f"
  } else {
    binary_design_terms <- "survey_year_f"
  }

  binary_formula <- formula_text(
    "miscarriage",
    unique(c(
      ns_term("HI_early", 4L),
      ns_term("CO_early", 4L),
      ns_term("PM_early", 4L),
      ns_term("EVI_early", 4L),
      binary_design_terms,
      ns_term("age_at_preg_start", 4L),
      "parity_cat", "previous_miscarriage_f",
      "interpregnancy_interval_f",
      binary_factor_terms,
      "(1 | mother_uid)", "(1 | CLUSTID)"
    ))
  )

  binary_model <- fit_or_load_glmmTMB(
    binary_formula,
    binary_data,
    "Pregnancy_level_binary_logit",
    "logit",
    "analysis_weight"
  )

  binary_fixed <- extract_fixed_ratios(
    binary_model, "Pregnancy_level_binary_logit", "OR"
  )
  binary_contrast <- contrast_raw_unit(
    binary_model,
    binary_data,
    "HI_early",
    "Pregnancy_level_binary_logit",
    "OR",
    1
  )

  write_csv(
    binary_fixed,
    file.path(tables_dir, "Sensitivity_binary_fixed_effect_ORs.csv")
  )
  write_csv(
    binary_contrast,
    file.path(tables_dir, "Sensitivity_binary_HI_contrast_OR.csv")
  )
  saveRDS(
    binary_model,
    file.path(models_dir, "Sensitivity_binary_logit_model.rds"),
    compress = FALSE
  )
}


# ---- 9. Socioeconomic effect modification ---------------------------------

# HI and WBGT are evaluated in separate models, never entered together.
for (thermal in c("HI", "WBT")) {
  thermal_label <- if (thermal == "WBT") "WBGT" else "HI"
  base_linear_environmental_terms <- c(
    thermal, "CO", "PM", "EVI_per_0_1"
  )

  interaction_models <- list()
  interaction_tests <- list()
  interaction_coefficients <- list()

  interaction_modifiers <- intersect(
    c("wealth_f", "health_insurance_f"),
    names(DT_primary)
  )

  for (modifier in interaction_modifiers) {
    model_name <- paste0(thermal_label, "_by_", modifier)

    interaction_formula <- formula_text(
      "miscarriage_event",
      c(
        base_linear_environmental_terms,
        adjustment_terms,
        paste0(thermal, ":", modifier)
      )
    )

    fitted_model <- fit_or_load_glmmTMB(
      interaction_formula,
      DT_primary,
      model_name,
      "cloglog",
      "analysis_weight"
    )

    interaction_models[[model_name]] <- fitted_model
    interaction_tests[[model_name]] <- joint_wald_test(
      fitted_model,
      paste0("(^", thermal, ":", modifier, "|", modifier, ".*:", thermal, "$)"),
      model_name,
      modifier
    )
    interaction_coefficients[[model_name]] <- extract_fixed_ratios(
      fitted_model, model_name, "HR"
    ) %>%
      filter(grepl(paste0("(^", thermal, "$|", thermal, ":|:", thermal, "$)"), term))
  }

  # Residence modification replaces strata_f with survey year plus residence.
  if (
    "residence_f" %in% names(DT_primary) &&
    nlevels(droplevels(factor(DT_primary$residence_f))) > 1L
  ) {
    residence_model_name <- paste0(thermal_label, "_by_residence_f")

    residence_formula <- formula_text(
      "miscarriage_event",
      c(
        base_linear_environmental_terms,
        residence_adjustment_terms,
        paste0(thermal, ":residence_f")
      )
    )

    residence_model <- fit_or_load_glmmTMB(
      residence_formula,
      DT_primary,
      residence_model_name,
      "cloglog",
      "analysis_weight"
    )

    interaction_models[[residence_model_name]] <- residence_model
    interaction_tests[[residence_model_name]] <- joint_wald_test(
      residence_model,
      paste0("(^", thermal, ":residence_f|residence_f.*:", thermal, "$)"),
      residence_model_name,
      "residence_f"
    )
    interaction_coefficients[[residence_model_name]] <- extract_fixed_ratios(
      residence_model, residence_model_name, "HR"
    ) %>%
      filter(grepl(paste0("(^", thermal, "$|", thermal, ":|:", thermal, "$)"), term))
  }

  category_contrasts <- list()
  for (model_name in names(interaction_models)) {
    modifier <- sub(paste0("^", thermal_label, "_by_"), "", model_name)
    for (category in levels(droplevels(factor(DT_primary[[modifier]])))) {
      category_contrasts[[paste(model_name, category)]] <- contrast_raw_unit(
        interaction_models[[model_name]], DT_primary, thermal, model_name,
        "HR", 1, modifier = modifier, category = category
      ) %>% mutate(modifier = modifier, category = category)
    }
  }
  if (length(category_contrasts) > 0L) {
    write_csv(bind_rows(category_contrasts), file.path(
      tables_dir, paste0("Effect_modification_", thermal_label, "_category_specific_HRs.csv")
    ))
  }

  if (length(interaction_models) > 0L) {
    write_csv(
      bind_rows(interaction_tests),
      file.path(tables_dir, paste0("Effect_modification_", thermal_label, "_joint_Wald_tests.csv"))
    )
    write_csv(
      bind_rows(interaction_coefficients),
      file.path(tables_dir, paste0("Effect_modification_", thermal_label, "_coefficients.csv"))
    )
    saveRDS(
      interaction_models,
      file.path(models_dir, paste0("Effect_modification_", thermal_label, "_models.rds")),
      compress = FALSE
    )
  }

}

# ---- 10. EVI-response heterogeneity by residence ----------------------------

if (
  "residence_f" %in% names(DT_primary) &&
  nlevels(droplevels(factor(DT_primary$residence_f))) > 1L
) {
  evi_residence_formula <- formula_text(
    "miscarriage_event",
    c(
      ns_term("HI", 4L),
      ns_term("CO", 4L),
      ns_term("PM", 4L),
      paste0(ns_term("EVI", 4L), " * residence_f"),
      residence_adjustment_terms
    )
  )

  evi_residence_model <- fit_or_load_glmmTMB(
    evi_residence_formula,
    DT_primary,
    "EVI_by_residence",
    "cloglog",
    "analysis_weight"
  )

  evi_residence_test <- joint_wald_test(
    evi_residence_model,
    "ns\\(EVI.*:residence_f|residence_f.*:ns\\(EVI",
    "EVI_by_residence",
    "residence_f"
  )

  write_csv(
    extract_fixed_ratios(
      evi_residence_model, "EVI_by_residence", "HR"
    ),
    file.path(
      tables_dir,
      "Effect_modification_EVI_by_residence_coefficients.csv"
    )
  )
  write_csv(
    evi_residence_test,
    file.path(
      tables_dir,
      "Effect_modification_EVI_by_residence_Wald_test.csv"
    )
  )
  saveRDS(
    evi_residence_model,
    file.path(
      models_dir,
      "Effect_modification_EVI_by_residence_model.rds"
    ),
    compress = FALSE
  )
}


# ---- 11. Differential under-ascertainment scenarios ------------------------

reporting_required <- c("wealth_f", "internet_ever_f")

if (all(reporting_required %in% names(DT_primary))) {
  wealth_text <- tolower(as.character(DT_primary$wealth_f))
  internet_text <- tolower(as.character(DT_primary$internet_ever_f))

  high_wealth <- grepl("richest|highest|^5$", wealth_text)
  internet_access <- grepl(
    "yes|ever|daily|weekly|^1$", internet_text
  )

  if (!any(high_wealth, na.rm = TRUE)) {
    wealth_levels <- levels(droplevels(factor(DT_primary$wealth_f)))
    high_wealth <- (
      as.character(DT_primary$wealth_f) == tail(wealth_levels, 1)
    )
  }

  if (!any(internet_access, na.rm = TRUE)) {
    internet_levels <- levels(
      droplevels(factor(DT_primary$internet_ever_f))
    )
    internet_access <- (
      as.character(DT_primary$internet_ever_f) ==
        tail(internet_levels, 1)
    )
  }

  DT_reporting <- DT_primary %>%
    mutate(
      high_reporting_access = high_wealth | internet_access
    )

  reporting_group_counts <- DT_reporting %>%
    distinct(preg_uid, high_reporting_access) %>%
    count(high_reporting_access, name = "pregnancies")

  write_csv(
    reporting_group_counts,
    file.path(tables_dir, "Sensitivity_reporting_group_counts.csv")
  )

  reporting_scenarios <- tibble(
    scenario = paste0("LowAccess_Se", c(50, 60, 70, 80)),
    sensitivity_low_access = c(0.50, 0.60, 0.70, 0.80),
    sensitivity_high_access = 0.95,
    assumed_specificity = 1.00
  )

  for (thermal in c("HI", "WBT")) {
    thermal_label <- if (thermal == "WBT") "WBGT" else "HI"
  reporting_formula <- formula_text(
    "miscarriage_event",
    c(thermal, "CO", "PM", "EVI_per_0_1", adjustment_terms)
  )

  reporting_models <- list()
  reporting_results <- list()

  for (row_index in seq_len(nrow(reporting_scenarios))) {
    scenario <- reporting_scenarios[row_index, ]
    scenario_name <- as.character(scenario$scenario[[1]])
    sensitivity_low <- scenario$sensitivity_low_access[[1]]
    sensitivity_high <- scenario$sensitivity_high_access[[1]]

    scenario_data <- DT_reporting %>%
      mutate(
        assumed_sensitivity = if_else(
          high_reporting_access,
          sensitivity_high,
          sensitivity_low
        ),
        # Event-record inverse-sensitivity weighting. This is an
        # assumption-based analysis, not outcome validation.
        ascertainment_factor = if_else(
          miscarriage_event == 1,
          1 / assumed_sensitivity,
          1
        ),
        bias_analysis_weight =
          analysis_weight * ascertainment_factor
      )

    fitted_model <- fit_or_load_glmmTMB(
      reporting_formula,
      scenario_data,
      paste0("Reporting_", thermal_label, "_", scenario_name),
      "cloglog",
      "bias_analysis_weight"
    )

    reporting_models[[scenario_name]] <- fitted_model
    reporting_results[[scenario_name]] <- extract_fixed_ratios(
      fitted_model, scenario_name, "HR"
    ) %>%
      filter(term == thermal) %>%
      mutate(
        sensitivity_low_access = sensitivity_low,
        sensitivity_high_access = sensitivity_high,
        assumed_specificity =
          scenario$assumed_specificity[[1]],
        interpretation =
          "Assumption-based event weighting; not outcome validation"
      )
  }

  write_csv(
    reporting_scenarios,
    file.path(tables_dir, paste0("Sensitivity_reporting_", thermal_label, "_assumptions.csv"))
  )
  write_csv(
    bind_rows(reporting_results),
    file.path(
      tables_dir,
      paste0("Sensitivity_reporting_underascertainment_", thermal_label, "_results.csv")
    )
  )
  saveRDS(
    reporting_models,
    file.path(
      models_dir,
      paste0("Sensitivity_reporting_underascertainment_", thermal_label, "_models.rds")
    ),
    compress = FALSE
  )
  }
}


# ---- 12. Compact comparison figure -----------------------------------------

contrast_plot_data <- bind_rows(sensitivity_contrasts) %>%
  filter(
    is.finite(estimate),
    is.finite(lower_95),
    is.finite(upper_95),
    estimate > 0,
    lower_95 > 0,
    upper_95 > 0
  ) %>%
  mutate(
    result = sprintf(
      "%.2f (%.2f–%.2f)", estimate, lower_95, upper_95
    ),
    model = factor(model, levels = rev(unique(model)))
  )

comparison_plot <- ggplot(
  contrast_plot_data,
  aes(
    x = estimate,
    y = model,
    xmin = lower_95,
    xmax = upper_95,
    shape = effect_measure,
    color = effect_measure
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_errorbarh(height = 0, linewidth = 0.75) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  scale_color_manual(
    values = c("HR" = "#B22222", "Interval OR" = "#0072B2")
  ) +
  scale_shape_manual(
    values = c("HR" = 16, "Interval OR" = 17)
  ) +
  labs(
    x = paste0(
      "Ratio for a one-unit increase in the thermal exposure ",
      "(log scale)"
    ),
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Sensitivity of Thermal-Exposure Associations",
    subtitle = paste0(
      "Nonlinear estimates are local contrasts from the exposure median; ",
      "the logit model is an interval odds ratio"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

save_high_resolution(
  comparison_plot,
  "Sensitivity_Exposure_Contrast_Forest",
  width = 10,
  height = 7.5
)


# ---- 13. Reproducibility record --------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "SessionInfo_04.txt")
)

message(
  "Script 04 completed.",
  "\nTables: ", tables_dir,
  "\nFigures: ", figures_dir,
  "\nCheckpoints: ", checkpoint_dir,
  "\nUnder-ascertainment results are assumption-based scenarios only."
)
