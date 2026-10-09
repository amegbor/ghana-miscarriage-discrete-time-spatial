# GitHub copy: portable paths; restricted data are not distributed.
# ============================================================================
# 02_Ghana_DiscreteTime_Analysis.R
# Primary discrete-time survival analysis of reported miscarriage in Ghana
#
# Input created by:
#   01_Ghana_Pregnancy_DataPrep_DiscreteTime_EnvLinkage.R
#
# Main estimand:
#   Conditional monthly miscarriage hazard from the estimated conception month
#   through six completed months (seven calendar-month risk intervals).
#   Sparse glmmTMB complementary-log-log mixed models include mother and
#   cluster random intercepts. Exponentiated linear exposure
#   coefficients are reported as discrete-time hazard ratios (HRs). Nonlinear
#   natural-spline effects are HR contrasts relative to the exposure median.
# ============================================================================

rm(list = ls())
gc(reset = TRUE)

# ---- 0. Packages ------------------------------------------------------------

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "haven",
  "glmmTMB", "ggplot2", "patchwork", "survey", "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(haven)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(survey)
})

options(survey.lonely.psu = "adjust")
set.seed(20260921)


# ---- 1. Open OnDemand RStudio paths and user settings ----------------------

# Run from the repository root, or set GHANA_MISCARRIAGE_ROOT explicitly.
# Private inputs and outputs can live entirely outside the Git checkout.
project_dir <- normalizePath(
  Sys.getenv("GHANA_MISCARRIAGE_ROOT", unset = "."),
  winslash = "/", mustWork = TRUE
)
data_dir <- Sys.getenv(
  "GHANA_DATA_DIR", unset = file.path(project_dir, "data")
)
results_dir <- file.path(project_dir, "Results")
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")
models_dir <- file.path(results_dir, "Models")
# Use a dedicated checkpoint namespace for the final analytic specification.
checkpoint_dir <- file.path(models_dir, "glmmTMB_checkpoints_final_v1")

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "Ghana_Discrete_Time_Monthly.rds")

# Natural cubic spline degrees of freedom for nonlinear exposure effects.
NS_DF_PRIMARY <- 4L

# Checkpoint behavior. By default, valid completed glmmTMB models are loaded
# from disk. Set GHANA_REFIT_MODELS=true before sourcing the script to refit
# every model from scratch.
REFIT_MODELS <- tolower(Sys.getenv("GHANA_REFIT_MODELS", unset = "false")) %in%
  c("true", "1", "yes")

# Respect the CPU allocation supplied by Slurm/OOD. With a 16-core allocation,
# the default is 12 glmmTMB threads, leaving capacity for RStudio and system
# processes. GHANA_N_THREADS can override this value but is still capped at the
# allocated and physically detected core counts.
allocated_text <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
ALLOCATED_CORES <- suppressWarnings(as.integer(allocated_text))
default_threads <- if (
  is.finite(ALLOCATED_CORES) && ALLOCATED_CORES > 0L
) {
  min(12L, ALLOCATED_CORES)
} else {
  4L
}

thread_text <- Sys.getenv(
  "GHANA_N_THREADS",
  unset = as.character(default_threads)
)
N_THREADS <- suppressWarnings(as.integer(thread_text))
if (is.na(N_THREADS) || N_THREADS < 1L) N_THREADS <- default_threads
if (is.finite(ALLOCATED_CORES) && ALLOCATED_CORES > 0L) {
  N_THREADS <- min(N_THREADS, ALLOCATED_CORES)
}
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.finite(detected_cores) && detected_cores > 0L) {
  N_THREADS <- min(N_THREADS, as.integer(detected_cores))
}
options(mc.cores = N_THREADS)

message(
  "Ghana miscarriage project root: ", project_dir,
  "\nExecution mode: Open OnDemand RStudio/interactive R",
  "\nglmmTMB threads: ", N_THREADS,
  " (must not exceed the cores requested in OOD)",
  "\nRefit completed checkpoints: ", REFIT_MODELS
)


# ---- 2. Helper functions ----------------------------------------------------

first_existing <- function(data, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% names(data)]

  if (length(hit) == 0) {
    if (required) {
      stop(
        "Could not find ", label %||% "a required variable", ". Tried: ",
        paste(candidates, collapse = ", ")
      )
    }
    return(NULL)
  }

  hit[1]
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_numeric_clean <- function(x) {
  if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)
  suppressWarnings(as.numeric(x))
}

as_factor_clean <- function(x) {
  if (inherits(x, "haven_labelled")) {
    x <- haven::as_factor(x, levels = "labels")
  } else {
    x <- factor(x)
  }
  droplevels(x)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(which.max(table(x)))[1]
}

formula_text <- function(outcome, terms) {
  as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

model_variables <- function(formula) {
  unique(all.vars(formula))
}

fit_or_load_glmmTMB <- function(formula, data, model_name) {
  needed <- unique(c(model_variables(formula), "analysis_weight"))
  missing_needed <- setdiff(needed, names(data))

  if (length(missing_needed) > 0) {
    stop(model_name, " is missing variables: ", paste(missing_needed, collapse = ", "))
  }

  checkpoint_file <- file.path(checkpoint_dir, paste0(model_name, ".rds"))

  if (file.exists(checkpoint_file) && !REFIT_MODELS) {
    message("Loading completed checkpoint: ", model_name)
    model <- readRDS(checkpoint_file)
    updated_model <- try(glmmTMB::up2date(model), silent = TRUE)
    if (inherits(updated_model, "glmmTMB")) model <- updated_model

    if (!inherits(model, "glmmTMB")) {
      stop("Checkpoint is not a glmmTMB model: ", checkpoint_file)
    }

    return(model)
  }

  dat <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  dat <- droplevels(dat)

  if (sum(dat$miscarriage_event, na.rm = TRUE) == 0) {
    stop(model_name, " has no miscarriage events after complete-case restriction.")
  }

  message(
    model_name, ": ", format(nrow(dat), big.mark = ","), " person-months; ",
    format(dplyr::n_distinct(dat$preg_uid), big.mark = ","), " pregnancies; ",
    sum(dat$miscarriage_event), " events."
  )

  model <- glmmTMB::glmmTMB(
    formula = formula,
    family = binomial(link = "cloglog"),
    data = dat,
    weights = analysis_weight,
    na.action = na.exclude,
    # glmmTMB constructs the random-effects design matrices sparsely by
    # default. Keeping the comparatively small fixed-effects matrix dense
    # avoids Matrix::sparse.model.matrix failures with ns() and wide factors.
    control = glmmTMB::glmmTMBControl(
      parallel = N_THREADS,
      optCtrl = list(iter.max = 1000L, eval.max = 1000L)
    )
  )

  attr(model, "analysis_n") <- nrow(dat)
  attr(model, "analysis_pregnancies") <- dplyr::n_distinct(dat$preg_uid)
  attr(model, "analysis_events") <- sum(dat$miscarriage_event)
  attr(model, "model_name") <- model_name

  if (!is.null(model$fit$convergence) && model$fit$convergence != 0L) {
    warning(model_name, " returned convergence code ", model$fit$convergence)
  }
  if (!is.null(model$sdr$pdHess) && !isTRUE(model$sdr$pdHess)) {
    warning(model_name, " has a non-positive-definite Hessian.")
  }

  saveRDS(model, checkpoint_file, compress = FALSE)
  message("Saved checkpoint: ", checkpoint_file)
  model
}

extract_fixed_effects <- function(model, model_name) {
  tab <- as.data.frame(summary(model)$coefficients$cond)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  tibble(
    model = model_name,
    term = tab$term,
    estimate = tab$Estimate,
    std_error = tab$`Std. Error`,
    z_value = tab$`z value`,
    exp_coefficient = exp(tab$Estimate),
    lower_95 = exp(tab$Estimate - 1.96 * tab$`Std. Error`),
    upper_95 = exp(tab$Estimate + 1.96 * tab$`Std. Error`),
    p_value = tab$`Pr(>|z|)`,
    interpretation = if_else(
      grepl("splines::ns\\(", tab$term),
      "Spline-basis coefficient; do not interpret as a standalone HR",
      "Exponentiated fixed-effect coefficient"
    )
  )
}

extract_spline_wald <- function(model, model_name, exposure) {
  beta <- glmmTMB::fixef(model)$cond
  covariance <- as.matrix(vcov(model)$cond)
  prefix <- paste0("splines::ns(", exposure, ", df = ")
  coefficient_names <- names(beta)[startsWith(names(beta), prefix)]

  if (length(coefficient_names) == 0L) {
    stop("No spline coefficients found for ", exposure, " in ", model_name)
  }

  beta_spline <- beta[coefficient_names]
  covariance_spline <- covariance[
    coefficient_names, coefficient_names, drop = FALSE
  ]
  degrees_freedom <- qr(covariance_spline)$rank
  wald_statistic <- tryCatch(
    as.numeric(crossprod(beta_spline, solve(covariance_spline, beta_spline))),
    error = function(e) NA_real_
  )

  tibble(
    model = model_name,
    exposure = exposure,
    spline_df = length(coefficient_names),
    test_df = degrees_freedom,
    wald_chi_square = wald_statistic,
    p_value = if_else(
      is.finite(wald_statistic),
      pchisq(wald_statistic, df = degrees_freedom, lower.tail = FALSE),
      NA_real_
    ),
    test = "Joint Wald test of all natural-spline basis coefficients"
  )
}

extract_model_fit <- function(fitted_model, model_name) {
  tibble(
    model = model_name,
    person_months = attr(fitted_model, "analysis_n"),
    pregnancies = attr(fitted_model, "analysis_pregnancies"),
    events = attr(fitted_model, "analysis_events"),
    log_likelihood = as.numeric(logLik(fitted_model)),
    AIC = AIC(fitted_model),
    BIC = BIC(fitted_model),
    convergence_code = fitted_model$fit$convergence,
    positive_definite_Hessian = fitted_model$sdr$pdHess
  )
}

# Correct nonlinear HR contrasts. Unlike exp(centered smooth), this explicitly
# compares each exposure value with the sample median while holding other terms
# fixed. The full model covariance matrix supplies the confidence interval.
predict_hr_curve <- function(model, data, exposure, model_name, n = 250L) {
  observed <- data[[exposure]]
  limits <- quantile(observed, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  x_values <- seq(limits[1], limits[2], length.out = n)
  reference_x <- median(observed, na.rm = TRUE)

  # The analytic data are complete for every model variable, although source
  # columns retained from DHS may still contain missing values.
  reference_row <- data[1, , drop = FALSE]

  for (variable in names(reference_row)) {
    x <- data[[variable]]

    if (is.factor(x)) {
      reference_row[[variable]] <- factor(
        mode_value(x), levels = levels(x)
      )
    } else if (is.numeric(x) && variable != "miscarriage_event") {
      reference_row[[variable]] <- median(x, na.rm = TRUE)
    }
  }

  newdata <- reference_row[rep(1, n), , drop = FALSE]
  newdata[[exposure]] <- x_values

  reference_data <- reference_row
  reference_data[[exposure]] <- reference_x

  prediction_data <- bind_rows(newdata, reference_data)
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
  reference_index <- n + 1L
  log_hr <- fitted_link[seq_len(n)] - fitted_link[reference_index]
  contrast_variance <-
    diag(prediction_covariance)[seq_len(n)] +
    prediction_covariance[reference_index, reference_index] -
    2 * prediction_covariance[seq_len(n), reference_index]
  se <- sqrt(pmax(0, contrast_variance))

  tibble(
    model = model_name,
    exposure = exposure,
    exposure_value = x_values,
    reference_value = reference_x,
    HR = exp(log_hr),
    lower_95 = exp(log_hr - 1.96 * se),
    upper_95 = exp(log_hr + 1.96 * se)
  )
}

make_curve_panel <- function(curve, data, exposure, x_label, title) {
  percentile_data <- tibble(
    probability = factor(c("25th", "75th", "90th"), levels = c("25th", "75th", "90th")),
    value = as.numeric(quantile(data[[exposure]], c(0.25, 0.75, 0.90), na.rm = TRUE))
  )

  rug_data <- tibble(value = data[[exposure]]) %>%
    filter(is.finite(value))

  if (nrow(rug_data) > 4000L) {
    rug_data <- slice_sample(rug_data, n = 4000L)
  }

  ggplot(curve, aes(x = exposure_value, y = HR)) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F4A6A6", alpha = 0.50) +
    geom_line(color = "#B22222", linewidth = 0.9) +
    geom_hline(yintercept = 1, linetype = 3, color = "grey35") +
    geom_vline(
      data = percentile_data,
      aes(xintercept = value, color = probability),
      linetype = 2,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_rug(data = rug_data, aes(x = value), inherit.aes = FALSE, alpha = 0.12, sides = "b") +
    scale_color_manual(values = c("25th" = "#1B9E77", "75th" = "#D95F02", "90th" = "#7570B3")) +
    labs(x = x_label, y = "Hazard ratio (95% CI)", title = title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

save_high_resolution <- function(plot, stem, width = 11, height = 8.5) {
  ggsave(
    file.path(figures_dir, paste0(stem, ".png")),
    plot = plot, width = width, height = height, units = "in", dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(figures_dir, paste0(stem, ".tiff")),
    plot = plot, width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
}

environmental_vif <- function(data, variables, specification) {
  map_dfr(variables, function(outcome) {
    predictors <- setdiff(variables, outcome)
    fit <- lm(
      reformulate(predictors, response = outcome),
      data = data[, c(outcome, predictors), drop = FALSE]
    )
    r_squared <- summary(fit)$r.squared
    tibble(
      specification = specification,
      exposure = outcome,
      R_squared_against_other_exposures = r_squared,
      VIF = 1 / (1 - r_squared)
    )
  })
}


# ---- 3. Read and harmonize the revised person-month data -------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, "\nRun script 01 first.")
}

DT_raw <- readRDS(input_file)

required_core <- c(
  "preg_uid", "mother_uid", "CLUSTID", "gest_month",
  "miscarriage_event", "survey_year", "age_at_preg_start",
  "prior_pregnancies", "parity_before_pregnancy", "previous_miscarriage",
  "interpregnancy_interval_months"
)

missing_core <- setdiff(required_core, names(DT_raw))
if (length(missing_core) > 0) {
  stop(
    "The revised person-month dataset is missing: ",
    paste(missing_core, collapse = ", ")
  )
}

# Permit either the revised short names or the original manuscript names.
exposure_candidates <- list(
  HI = c("HI", "HI_T1"),
  WBT = c("WBT", "WBT_T1"),
  CO = c("CO", "CO_T1"),
  PM = c("PM", "PM_T1"),
  EVI = c("EVI", "EVI_T1")
)

# Preserve every environmental exposure in its original measurement unit.
# No IQR or z-score transformation is used in this Ghana analysis.

# Prespecified covariates for adjustment and complete-case eligibility.
factor_candidates <- list(
  education_f = c("v106", "education"),
  wealth_f = c("v190", "wealth"),
  no_weekly_media_f = c("rc_media_none", "no_weekly_media"),
  health_insurance_f = c("rc_hins_any", "health_insurance"),
  sanitation_f = c("ph_sani_improve", "sanitation"),
  water_improved_f = c("ph_wtr_improve", "water_improved"),
  internet_ever_f = c("rc_intr_ever", "internet_ever"),
  residence_f = c("v025", "URBAN_RURA", "place_residence")
)

DT <- DT_raw
variable_map <- list()

for (new_name in names(exposure_candidates)) {
  source_name <- first_existing(
    DT, exposure_candidates[[new_name]], required = TRUE, label = new_name
  )
  DT[[new_name]] <- as_numeric_clean(DT[[source_name]])
  variable_map[[length(variable_map) + 1L]] <- tibble(
    analysis_variable = new_name, source_variable = source_name
  )
}

for (new_name in names(factor_candidates)) {
  source_name <- first_existing(
    DT, factor_candidates[[new_name]], required = FALSE, label = new_name
  )

  if (!is.null(source_name)) {
    x <- DT[[source_name]]
    if (new_name == "sanitation_f") {
      x_num <- as_numeric_clean(x)
      x_num[x_num == 16] <- 1
      x <- factor(
        x_num,
        levels = c(1, 2, 3),
        labels = c("Improved", "Unimproved", "Open defecation")
      )
    }
    DT[[new_name]] <- as_factor_clean(x)
    variable_map[[length(variable_map) + 1L]] <- tibble(
      analysis_variable = new_name, source_variable = source_name
    )
  }
}

variable_map <- bind_rows(variable_map)
write_csv(variable_map, file.path(tables_dir, "Variable_mapping.csv"))

# DHS sample weights. v005 is normally stored with six implied decimals.
weight_source <- first_existing(
  DT, c("v005", "sample_weight", "sampling_weight", "weight"),
  required = TRUE, label = "DHS sampling weight"
)

DT$raw_sample_weight <- as_numeric_clean(DT[[weight_source]])
if (median(DT$raw_sample_weight, na.rm = TRUE) > 100) {
  DT$raw_sample_weight <- DT$raw_sample_weight / 1e6
}

DT <- DT %>%
  mutate(
    miscarriage_event = as.integer(miscarriage_event),
    gest_month_f = factor(gest_month, levels = 1:7),
    survey_year_f = factor(survey_year),
    mother_uid = factor(mother_uid),
    CLUSTID = factor(CLUSTID),
    parity_cat = cut(
      parity_before_pregnancy,
      breaks = c(-Inf, 0, 1, 2, Inf),
      labels = c("0", "1", "2", "3+"),
      right = TRUE
    ),
    previous_miscarriage_f = factor(
      previous_miscarriage,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    interpregnancy_interval_f = factor(
      case_when(
        prior_pregnancies == 0 ~ "No previous pregnancy",
        is.na(interpregnancy_interval_months) ~ "Missing/unknown",
        interpregnancy_interval_months < 6 ~ "<6 months",
        interpregnancy_interval_months < 12 ~ "6-11 months",
        interpregnancy_interval_months < 24 ~ "12-23 months",
        interpregnancy_interval_months < 60 ~ "24-59 months",
        TRUE ~ "60+ months"
      ),
      levels = c(
        "No previous pregnancy", "<6 months", "6-11 months",
        "12-23 months", "24-59 months", "60+ months", "Missing/unknown"
      )
    )
  )

# Normalize weights at the pregnancy level within each survey, then carry the
# same normalized weight to every monthly risk interval. Normalizing directly
# on the expanded person-month file would allow longer pregnancies to influence
# the normalization constant more heavily.
weight_consistency <- DT %>%
  group_by(preg_uid) %>%
  summarise(n_weights = n_distinct(raw_sample_weight), .groups = "drop")

if (any(weight_consistency$n_weights > 1)) {
  stop("At least one pregnancy has more than one DHS sampling weight.")
}

pregnancy_weights <- DT %>%
  distinct(preg_uid, survey_year_f, raw_sample_weight) %>%
  group_by(survey_year_f) %>%
  mutate(analysis_weight = raw_sample_weight / mean(raw_sample_weight, na.rm = TRUE)) %>%
  ungroup() %>%
  select(preg_uid, analysis_weight)

DT <- DT %>%
  select(-any_of("analysis_weight")) %>%
  left_join(pregnancy_weights, by = "preg_uid")

# Add the DHS stratum as a fixed design term when available.
strata_source <- first_existing(
  DT, c("v022", "v023", "strata", "sample_stratum"), required = FALSE
)

if (!is.null(strata_source)) {
  DT$strata_f <- interaction(
    DT$survey_year_f,
    as_factor_clean(DT[[strata_source]]),
    drop = TRUE,
    lex.order = TRUE
  )
}


# ---- 4. Essential person-month quality checks -------------------------------

pregnancy_event_check <- DT %>%
  group_by(preg_uid) %>%
  summarise(events = sum(miscarriage_event, na.rm = TRUE), .groups = "drop")

if (any(pregnancy_event_check$events > 1)) {
  stop("At least one pregnancy contributes more than one miscarriage event.")
}

if (all(c("interval_date", "outcome_date") %in% names(DT))) {
  event_date_check <- DT %>%
    filter(miscarriage_event == 1) %>%
    summarise(all_match = all(interval_date == outcome_date, na.rm = TRUE)) %>%
    pull(all_match)

  if (!isTRUE(event_date_check)) {
    stop("At least one miscarriage event interval does not match its outcome month.")
  }
}

flow_table <- tibble(
  stage = c(
    "Person-month records from script 01",
    "Unique pregnancies from script 01",
    "Unique mothers",
    "Observed miscarriage events"
  ),
  n = c(
    nrow(DT),
    n_distinct(DT$preg_uid),
    n_distinct(DT$mother_uid),
    sum(DT$miscarriage_event, na.rm = TRUE)
  )
)
write_csv(flow_table, file.path(tables_dir, "Cohort_flow_initial.csv"))

if ("duration_for_calendar" %in% names(DT)) {
  duration_distribution <- DT %>%
    filter(miscarriage_event == 1) %>%
    count(duration_for_calendar, name = "miscarriages") %>%
    mutate(percent = 100 * miscarriages / sum(miscarriages))

  write_csv(
    duration_distribution,
    file.path(tables_dir, "Miscarriage_gestational_duration_distribution.csv")
  )
}


# ---- 5. Model terms and complete-case analytic dataset ----------------------

available_factor_terms <- intersect(names(factor_candidates), names(DT))
available_factor_terms <- available_factor_terms[
  vapply(DT[available_factor_terms], function(x) nlevels(droplevels(x)) > 1, logical(1))
]

reproductive_terms <- c(
  "splines::ns(age_at_preg_start, df = 4)",
  "parity_cat",
  "previous_miscarriage_f",
  "interpregnancy_interval_f"
)

# strata_f is explicitly constructed as survey_year x source stratum. It
# therefore already contains survey-wave information and must not be entered
# together with survey_year_f, which creates an exactly redundant column.
if ("strata_f" %in% names(DT) && nlevels(DT$strata_f) > 1) {
  design_terms <- c("gest_month_f", "strata_f")
} else {
  design_terms <- c("gest_month_f", "survey_year_f")
}

# Survey strata already encode the urban/rural residence component. Including
# residence_f with strata_f creates a redundant fixed-effect column.
if ("strata_f" %in% design_terms) {
  available_factor_terms <- setdiff(
    available_factor_terms,
    "residence_f"
  )
}

adjustment_fixed_terms <- unique(c(
  design_terms,
  reproductive_terms,
  available_factor_terms
))
wave_only_fixed_terms <- design_terms

# Use one common analytic sample for all primary models. This prevents changes
# in model coefficients from being driven only by different missing-data sets.
complete_case_variables <- unique(c(
  "miscarriage_event", "preg_uid", "mother_uid", "CLUSTID", "analysis_weight",
  "gest_month_f", "survey_year_f", "age_at_preg_start", "parity_cat",
  "previous_miscarriage_f", "interpregnancy_interval_f",
  names(exposure_candidates), available_factor_terms,
  if ("strata_f" %in% design_terms) "strata_f" else NULL
))

DT_primary <- DT[complete.cases(DT[, complete_case_variables, drop = FALSE]), , drop = FALSE]
DT_primary <- droplevels(DT_primary)

flow_table <- bind_rows(
  flow_table,
  tibble(
    stage = c(
      "Primary complete-case person-month records",
      "Primary complete-case pregnancies",
      "Primary complete-case miscarriage events"
    ),
    n = c(
      nrow(DT_primary),
      n_distinct(DT_primary$preg_uid),
      sum(DT_primary$miscarriage_event)
    )
  )
)
write_csv(flow_table, file.path(tables_dir, "Cohort_flow_primary.csv"))

# Mark analysis-ready datasets for downstream compatibility checks.
attr(DT, "ghana_analysis_ready_final") <- TRUE
attr(DT_primary, "ghana_analysis_ready_final") <- TRUE
saveRDS(DT, file.path(models_dir, "Ghana_DT_AnalysisReady_AllRecords.rds"))
saveRDS(DT_primary, file.path(models_dir, "Ghana_DT_AnalysisReady_PrimaryCC.rds"))


# ---- 6. Exposure correlation and collinearity diagnostics ------------------

correlation_exposures <- names(exposure_candidates)
correlation_data <- DT_primary[, correlation_exposures, drop = FALSE]

correlation_spearman <- cor(
  correlation_data,
  method = "spearman",
  use = "pairwise.complete.obs"
)

write_csv(
  as.data.frame(correlation_spearman) %>% rownames_to_column("exposure"),
  file.path(tables_dir, "Exposure_correlation_Spearman.csv")
)

# Pairwise two-sided p-values and pairwise sample sizes. exact = FALSE is used
# because the environmental series contain ties and the analytic dataset is
# large. These tests are descriptive: they do not account for repeated
# person-months, spatial dependence, or the complex survey design.
correlation_p_values <- matrix(
  NA_real_,
  nrow = length(correlation_exposures),
  ncol = length(correlation_exposures),
  dimnames = list(correlation_exposures, correlation_exposures)
)
correlation_pair_n <- matrix(
  NA_integer_,
  nrow = length(correlation_exposures),
  ncol = length(correlation_exposures),
  dimnames = list(correlation_exposures, correlation_exposures)
)

for (row_index in seq_along(correlation_exposures)) {
  for (column_index in seq_len(row_index)) {
    row_variable <- correlation_exposures[row_index]
    column_variable <- correlation_exposures[column_index]
    complete_pair <- is.finite(correlation_data[[row_variable]]) &
      is.finite(correlation_data[[column_variable]])
    pair_n <- sum(complete_pair)

    correlation_pair_n[row_index, column_index] <- pair_n
    correlation_pair_n[column_index, row_index] <- pair_n

    if (row_index == column_index) {
      next
    }

    if (
      pair_n >= 3L &&
        n_distinct(correlation_data[[row_variable]][complete_pair]) > 1L &&
        n_distinct(correlation_data[[column_variable]][complete_pair]) > 1L
    ) {
      pair_test <- suppressWarnings(
        cor.test(
          correlation_data[[row_variable]][complete_pair],
          correlation_data[[column_variable]][complete_pair],
          method = "spearman",
          alternative = "two.sided",
          exact = FALSE
        )
      )
      correlation_p_values[row_index, column_index] <- pair_test$p.value
      correlation_p_values[column_index, row_index] <- pair_test$p.value
    }
  }
}

correlation_long <- as.data.frame(
  as.table(correlation_spearman),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rename(
    exposure_row = Var1,
    exposure_column = Var2,
    spearman_rho = Freq
  ) %>%
  left_join(
    as.data.frame(as.table(correlation_p_values), stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      rename(
        exposure_row = Var1,
        exposure_column = Var2,
        p_value = Freq
      ),
    by = c("exposure_row", "exposure_column")
  ) %>%
  left_join(
    as.data.frame(as.table(correlation_pair_n), stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      rename(
        exposure_row = Var1,
        exposure_column = Var2,
        pairwise_n = Freq
      ),
    by = c("exposure_row", "exposure_column")
  ) %>%
  mutate(
    row_index = match(exposure_row, correlation_exposures),
    column_index = match(exposure_column, correlation_exposures),
    significance = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    p_value_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", formatC(p_value, format = "f", digits = 3))
    ),
    cell_label = case_when(
      exposure_row == exposure_column ~
        paste0("rho = ", formatC(spearman_rho, format = "f", digits = 2)),
      TRUE ~ paste0(
        "rho = ", formatC(spearman_rho, format = "f", digits = 2),
        significance, "\n", p_value_label
      )
    ),
    text_color = if_else(abs(spearman_rho) >= 0.55, "white", "black")
  )

write_csv(
  correlation_long %>%
    select(
      exposure_row, exposure_column, spearman_rho,
      p_value, significance, pairwise_n
    ),
  file.path(tables_dir, "Exposure_correlation_Spearman_with_pvalues.csv")
)

exposure_display_labels <- c(
  HI = "Heat index",
  WBT = "WBGT",
  CO = "CO",
  PM = "PM2.5",
  EVI = "EVI"
)

correlation_plot <- correlation_long %>%
  filter(row_index >= column_index) %>%
  ggplot(aes(x = exposure_column, y = exposure_row, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = cell_label, color = text_color),
    lineheight = 0.95,
    size = 3.6,
    fontface = "bold"
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(
    limits = correlation_exposures,
    labels = exposure_display_labels,
    position = "top",
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(correlation_exposures),
    labels = exposure_display_labels,
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = "Spearman Correlations Among Environmental Exposures",
    subtitle = paste0(
      "Cells report rho and two-sided p-values ",
      "(* p < 0.05; ** p < 0.01; *** p < 0.001)"
    ),
    caption = paste0(
      "Pairwise-complete observations from the primary analytic sample. ",
      "P-values are descriptive and are not adjusted for repeated person-months, ",
      "spatial dependence, survey clustering, or multiple comparisons."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0, size = 8),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figures_dir, "Exposure_correlation_Spearman.tiff"),
  plot = correlation_plot,
  device = "tiff",
  width = 8.5,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

vif_table <- bind_rows(
  environmental_vif(
    DT_primary, c("HI", "WBT", "CO", "PM", "EVI"),
    "All exposures together - diagnostic only"
  ),
  environmental_vif(
    DT_primary, c("HI", "CO", "PM", "EVI"),
    "Primary HI multi-exposure family"
  ),
  environmental_vif(
    DT_primary, c("WBT", "CO", "PM", "EVI"),
    "Primary WBGT multi-exposure family"
  )
)

write_csv(vif_table, file.path(tables_dir, "Exposure_VIF_diagnostics.csv"))


# ---- 7. Linear exposure-effect multilevel discrete-time models --------------

# These are linear exposure-effect specifications of the same discrete-time
# cloglog model, not Gaussian linear regressions. Exposures remain in their
# indicated measurement units throughout. EVI is divided by 0.1 so that its
# linear coefficient is directly interpretable per 0.1-unit increase.
#
# Unadjusted:
#   raw exposure + gestational month + survey wave/design + mother and cluster
#   random intercepts.
#
# Adjusted:
#   unadjusted structure + reproductive and socioeconomic covariates.
#
# CLUSTID is preferred over raw v021 because it includes survey year and avoids
# cluster-number collisions between the 2017 and 2022 surveys.

linear_exposure_dictionary <- tibble(
  exposure = c("HI", "WBT", "CO", "PM", "EVI"),
  exposure_label = c(
    "Heat index",
    "WBGT",
    "Carbon monoxide",
    "PM2.5",
    "Enhanced Vegetation Index"
  ),
  linear_term = c("HI", "WBT", "CO", "PM", "EVI_per_0_1"),
  original_unit = c(
    "1 degree C",
    "1 degree C",
    "1 ppbv",
    "1 microgram/m^3",
    "0.1 EVI unit"
  )
)

# Create the reporting-scale EVI variable before validating model terms.
# One unit of EVI_per_0_1 corresponds to 0.1 unit of the original EVI.
DT_primary$EVI_per_0_1 <- DT_primary$EVI / 0.1

missing_linear_terms <- setdiff(
  linear_exposure_dictionary$linear_term,
  names(DT_primary)
)

if (length(missing_linear_terms) > 0) {
  stop(
    "The following raw exposure variables are missing: ",
    paste(missing_linear_terms, collapse = ", ")
  )
}


# ---- 7.1 Create the linear-model analysis data ------------------------------

DT_linear <- DT_primary

DT_linear$mother_uid <- droplevels(factor(DT_linear$mother_uid))
DT_linear$CLUSTID <- droplevels(factor(DT_linear$CLUSTID))
DT_linear <- droplevels(DT_linear)

# Scaling by 0.1 changes only the reporting increment: a one-unit increase in
# EVI_per_0_1 equals a 0.1-unit increase in the original EVI variable.
DT_linear$EVI_per_0_1 <- DT_linear$EVI / 0.1

# ---- 7.2 Model specifications ----------------------------------------------

linear_random_terms <- c(
  "(1 | mother_uid)",
  "(1 | CLUSTID)"
)

# "Unadjusted" means no individual or household confounders. The baseline
# hazard, survey design terms, and multilevel structure remain fundamental
# components of the discrete-time model.
linear_unadjusted_terms <- unique(c(
  design_terms,
  linear_random_terms
))

linear_adjusted_terms <- unique(c(
  adjustment_fixed_terms,
  linear_random_terms
))


# ---- 7.3 Fit one raw-scale linear model for each exposure ------------------

linear_unadjusted_models <- list()
linear_adjusted_models <- list()

for (row_index in seq_len(nrow(linear_exposure_dictionary))) {
  exposure_id <- linear_exposure_dictionary$exposure[row_index]
  exposure_term <- linear_exposure_dictionary$linear_term[row_index]

  model_suffix <- if (exposure_id == "EVI") {
    "EVI_per_0_1"
  } else {
    exposure_id
  }

  linear_unadjusted_models[[exposure_id]] <- fit_or_load_glmmTMB(
    formula = formula_text(
      "miscarriage_event",
      c(exposure_term, linear_unadjusted_terms)
    ),
    data = DT_linear,
    model_name = paste0("Linear_unadjusted_", model_suffix)
  )

  linear_adjusted_models[[exposure_id]] <- fit_or_load_glmmTMB(
    formula = formula_text(
      "miscarriage_event",
      c(exposure_term, linear_adjusted_terms)
    ),
    data = DT_linear,
    model_name = paste0("Linear_adjusted_", model_suffix)
  )
}


# ---- 7.5 Extract raw-unit exposure HRs -------------------------------------

extract_linear_exposure_hr <- function(
    model,
    model_type,
    exposure_id,
    exposure_label,
    exposure_term,
    original_unit
) {
  parameter_table <- as.data.frame(summary(model)$coefficients$cond)
  parameter_table$term <- rownames(parameter_table)
  rownames(parameter_table) <- NULL

  exposure_row <- parameter_table[
    parameter_table$term == exposure_term,
    ,
    drop = FALSE
  ]

  if (nrow(exposure_row) != 1L) {
    stop(
      "Could not uniquely identify term ", exposure_term,
      " in the ", model_type, " model."
    )
  }

  estimate <- exposure_row$Estimate[1]
  standard_error <- exposure_row[["Std. Error"]][1]
  p_value <- exposure_row[["Pr(>|z|)"]][1]

  tibble(
    exposure = exposure_id,
    exposure_label = exposure_label,
    original_unit = original_unit,
    model_type = model_type,
    exposure_term = exposure_term,
    estimate_log_HR = estimate,
    standard_error = standard_error,
    HR_per_unit = exp(estimate),
    lower_95 = exp(estimate - 1.96 * standard_error),
    upper_95 = exp(estimate + 1.96 * standard_error),
    p_value = p_value,
    person_months = attr(model, "analysis_n"),
    pregnancies = attr(model, "analysis_pregnancies"),
    miscarriage_events = attr(model, "analysis_events")
  )
}

linear_results <- map_dfr(
  seq_len(nrow(linear_exposure_dictionary)),
  function(row_index) {
    exposure_id <- linear_exposure_dictionary$exposure[row_index]
    exposure_label <- linear_exposure_dictionary$exposure_label[row_index]
    exposure_term <- linear_exposure_dictionary$linear_term[row_index]
    original_unit <- linear_exposure_dictionary$original_unit[row_index]

    bind_rows(
      extract_linear_exposure_hr(
        model = linear_unadjusted_models[[exposure_id]],
        model_type = "Unadjusted",
        exposure_id = exposure_id,
        exposure_label = exposure_label,
        exposure_term = exposure_term,
        original_unit = original_unit
      ),
      extract_linear_exposure_hr(
        model = linear_adjusted_models[[exposure_id]],
        model_type = "Adjusted",
        exposure_id = exposure_id,
        exposure_label = exposure_label,
        exposure_term = exposure_term,
        original_unit = original_unit
      )
    )
  }
) %>%
  mutate(
    p_value_formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ formatC(p_value, format = "f", digits = 3)
    ),
    result = paste0(
      formatC(HR_per_unit, format = "f", digits = 2),
      " (",
      formatC(lower_95, format = "f", digits = 2),
      "-",
      formatC(upper_95, format = "f", digits = 2),
      ")"
    ),
    exposure_order = match(exposure, linear_exposure_dictionary$exposure),
    model_order = match(model_type, c("Unadjusted", "Adjusted"))
  ) %>%
  arrange(exposure_order, model_order)

if (nrow(linear_results) != 2L * nrow(linear_exposure_dictionary)) {
  stop(
    "Linear-result extraction returned ", nrow(linear_results),
    " rows; expected ", 2L * nrow(linear_exposure_dictionary), "."
  )
}


# ---- 7.6 Model storage and fit information ---------------------------------

all_linear_models <- c(
  setNames(
    linear_unadjusted_models,
    paste0("Linear_unadjusted_", names(linear_unadjusted_models))
  ),
  setNames(
    linear_adjusted_models,
    paste0("Linear_adjusted_", names(linear_adjusted_models))
  )
)

linear_model_fit <- imap_dfr(
  all_linear_models,
  function(fitted_model, model_name) {
    tibble(
      model = model_name,
      person_months = attr(fitted_model, "analysis_n"),
      pregnancies = attr(fitted_model, "analysis_pregnancies"),
      miscarriage_events = attr(fitted_model, "analysis_events"),
      log_likelihood = as.numeric(logLik(fitted_model)),
      AIC = AIC(fitted_model),
      BIC = BIC(fitted_model),
      convergence_code = fitted_model$fit$convergence,
      positive_definite_Hessian = fitted_model$sdr$pdHess
    )
  }
)

saveRDS(
  all_linear_models,
  file.path(models_dir, "Linear_exposure_discrete_time_models.rds"),
  compress = FALSE
)


# ---- 7.7 Export results to an Excel workbook --------------------------------

linear_model_specification <- tibble(
  component = c(
    "Outcome",
    "Link function",
    "Exposure scale",
    "Unadjusted model",
    "Adjusted model",
    "Mother level",
    "Cluster level",
    "Sampling weights",
    "Interpretation"
  ),
  specification = c(
    "Monthly miscarriage event indicator",
    "Complementary log-log glmmTMB mixed model",
    paste0(
      "Raw units with no IQR or z-score transformation; EVI is divided by ",
      "0.1 and reported per 0.1-unit increase"
    ),
    paste(linear_unadjusted_terms, collapse = " + "),
    paste(linear_adjusted_terms, collapse = " + "),
    "Sparse random intercept: (1 | mother_uid)",
    "Sparse random intercept: (1 | CLUSTID)",
    "Normalized DHS/GMHS analysis_weight",
    "Exponentiated coefficient is the monthly hazard ratio per indicated increment"
  )
)

linear_excel_results <- linear_results %>%
  select(
    exposure,
    exposure_label,
    original_unit,
    model_type,
    estimate_log_HR,
    standard_error,
    HR_per_unit,
    lower_95,
    upper_95,
    p_value,
    p_value_formatted,
    result,
    person_months,
    pregnancies,
    miscarriage_events
  )

linear_excel_file <- file.path(
  tables_dir,
  "Linear_exposure_discrete_time_results.xlsx"
)

openxlsx::write.xlsx(
  x = list(
    Linear_HRs = linear_excel_results,
    Exposure_scales = linear_exposure_dictionary,
    Model_fit = linear_model_fit,
    Specification = linear_model_specification
  ),
  file = linear_excel_file,
  overwrite = TRUE
)


# ---- 7.8 Forest plot of raw-unit linear exposure HRs -----------------------

linear_plot_labels <- paste0(
  linear_exposure_dictionary$exposure_label,
  " (",
  linear_exposure_dictionary$original_unit,
  ")"
)

linear_forest_data <- linear_results %>%
  mutate(
    model_type = factor(model_type, levels = c("Unadjusted", "Adjusted")),
    exposure_plot_label = factor(
      paste0(exposure_label, " (", original_unit, ")"),
      levels = rev(linear_plot_labels)
    )
  ) %>%
  filter(is.finite(HR_per_unit), HR_per_unit > 0)

# Do not remove a point estimate merely because its confidence limits are not
# finite. Such models remain visible as points and are flagged in Model_fit;
# confidence intervals are drawn only where both limits are valid.
linear_forest_ci_data <- linear_forest_data %>%
  filter(
    is.finite(lower_95),
    is.finite(upper_95),
    lower_95 > 0,
    upper_95 > 0
  )

linear_forest_plot <- ggplot(
  linear_forest_data,
  aes(
    x = HR_per_unit,
    y = exposure_plot_label,
    color = model_type,
    shape = model_type,
    group = model_type
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.7,
    color = "grey35"
  ) +
  geom_errorbarh(
    data = linear_forest_ci_data,
    aes(xmin = lower_95, xmax = upper_95),
    position = position_dodge(width = 0.60),
    height = 0,
    linewidth = 0.75
  ) +
  geom_point(
    position = position_dodge(width = 0.60),
    size = 2.8
  ) +
  scale_x_log10(
    breaks = scales::breaks_log(n = 7),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  scale_color_manual(
    values = c("Unadjusted" = "#0072B2", "Adjusted" = "#D55E00")
  ) +
  scale_shape_manual(values = c("Unadjusted" = 16, "Adjusted" = 17)) +
  labs(
    x = "Monthly hazard ratio per indicated increment (log scale; 95% CI)",
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Linear Associations Between Environmental Exposures and Miscarriage",
    subtitle = paste0(
      "EVI per 0.1-unit increase; other exposures in their indicated units"
    ),
    caption = paste0(
      "Both specifications include gestational-month and survey-design terms ",
      "plus mother and cluster random intercepts."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0, size = 8),
    legend.position = "bottom"
  )

save_high_resolution(
  linear_forest_plot,
  "Figure_Linear_Exposure_HRs_Forest",
  width = 9,
  height = 6.5
)

message(
  "Linear raw-unit exposure models completed.\n",
  "Excel results: ", linear_excel_file, "\n",
  "Forest plot: ",
  file.path(figures_dir, "Figure_Linear_Exposure_HRs_Forest.tiff")
)


# Free the ten linear model objects after their combined RDS, Excel workbook,
# and forest plot have been saved. Individual checkpoints remain available.
rm(linear_unadjusted_models, linear_adjusted_models, all_linear_models)
invisible(gc(full = TRUE))


# ---- 8. Fast primary nonlinear spline mixed models -------------------------

# Fastest substantive nonlinear option:
#   1. Fit natural-spline models only for the prespecified thermal metrics HI
#      and WBGT, each with design-only and fully adjusted covariates.
#   2. Retain two adjusted multi-exposure models in which CO, PM, and EVI also
#      enter as natural cubic splines.
#   3. Use the same sparse mother and cluster random intercepts as the
#      linear models. HI and WBGT are not entered together because they are
#      alternative, highly correlated measures of thermal stress.
primary_thermal_exposures <- c("HI", "WBT")

missing_primary_thermal <- setdiff(primary_thermal_exposures, names(DT_linear))
if (length(missing_primary_thermal) > 0) {
  stop(
    "Missing primary thermal exposure(s): ",
    paste(missing_primary_thermal, collapse = ", ")
  )
}

natural_spline_term <- function(exposure) {
  paste0("splines::ns(", exposure, ", df = ", NS_DF_PRIMARY, ")")
}

nonlinear_wave_terms <- unique(c(
  wave_only_fixed_terms,
  linear_random_terms
))

nonlinear_adjusted_terms <- unique(c(
  adjustment_fixed_terms,
  linear_random_terms
))

single_wave_models <- list()
single_adjusted_models <- list()

for (exposure in primary_thermal_exposures) {
  spline_term <- natural_spline_term(exposure)

  single_wave_models[[exposure]] <- fit_or_load_glmmTMB(
    formula_text(
      "miscarriage_event",
      c(spline_term, nonlinear_wave_terms)
    ),
    DT_linear,
    paste0("Spline_M1_wave_", exposure)
  )

  single_adjusted_models[[exposure]] <- fit_or_load_glmmTMB(
    formula_text(
      "miscarriage_event",
      c(spline_term, nonlinear_adjusted_terms)
    ),
    DT_linear,
    paste0("Spline_M2_adjusted_", exposure)
  )
}

full_hi_terms <- c(
  natural_spline_term("HI"),
  natural_spline_term("CO"),
  natural_spline_term("PM"),
  natural_spline_term("EVI"),
  nonlinear_adjusted_terms
)

full_wbt_terms <- c(
  natural_spline_term("WBT"),
  natural_spline_term("CO"),
  natural_spline_term("PM"),
  natural_spline_term("EVI"),
  nonlinear_adjusted_terms
)

full_hi_model <- fit_or_load_glmmTMB(
  formula_text("miscarriage_event", full_hi_terms),
  DT_linear,
  "Spline_M3_multi_HI"
)

full_wbt_model <- fit_or_load_glmmTMB(
  formula_text("miscarriage_event", full_wbt_terms),
  DT_linear,
  "Spline_M3_multi_WBT"
)

all_models <- c(
  setNames(
    single_wave_models,
    paste0("Spline_M1_wave_", names(single_wave_models))
  ),
  setNames(
    single_adjusted_models,
    paste0("Spline_M2_adjusted_", names(single_adjusted_models))
  ),
  list(
    Spline_M3_multi_HI = full_hi_model,
    Spline_M3_multi_WBT = full_wbt_model
  )
)

saveRDS(
  all_models,
  file.path(models_dir, "Primary_discrete_time_spline_GLMM_models.rds"),
  compress = FALSE
)

message(
  "Fast nonlinear stage completed: ", length(all_models),
  " checkpointed spline mixed models."
)


# ---- 9. Export nonlinear model results and diagnostics ----------------------

fixed_effect_results <- imap_dfr(all_models, extract_fixed_effects)
fit_results <- imap_dfr(all_models, extract_model_fit)

model_spline_exposures <- list(
  Spline_M1_wave_HI = "HI",
  Spline_M1_wave_WBT = "WBT",
  Spline_M2_adjusted_HI = "HI",
  Spline_M2_adjusted_WBT = "WBT",
  Spline_M3_multi_HI = c("HI", "CO", "PM", "EVI"),
  Spline_M3_multi_WBT = c("WBT", "CO", "PM", "EVI")
)

spline_omnibus_results <- imap_dfr(
  model_spline_exposures,
  function(exposures, model_name) {
    map_dfr(
      exposures,
      ~ extract_spline_wald(all_models[[model_name]], model_name, .x)
    )
  }
)

random_effect_results <- imap_dfr(
  all_models,
  function(model, model_name) {
    imap_dfr(
      VarCorr(model)$cond,
      function(covariance_matrix, grouping_variable) {
        covariance_matrix <- as.matrix(covariance_matrix)
        tibble(
          model = model_name,
          grouping_variable = grouping_variable,
          random_effect_term = rownames(covariance_matrix),
          variance = diag(covariance_matrix),
          standard_deviation = sqrt(diag(covariance_matrix))
        )
      }
    )
  }
)

write_csv(
  fixed_effect_results,
  file.path(tables_dir, "Primary_spline_GLMM_fixed_effects.csv")
)
write_csv(
  spline_omnibus_results,
  file.path(tables_dir, "Primary_spline_omnibus_tests.csv")
)
write_csv(
  random_effect_results,
  file.path(tables_dir, "Primary_spline_GLMM_random_effects.csv")
)
write_csv(fit_results, file.path(tables_dir, "Primary_spline_GLMM_model_fit.csv"))


# ---- 10. Weighted pregnancy-level descriptive summaries --------------------

pregnancy_summary <- DT_primary %>%
  arrange(preg_uid, gest_month) %>%
  group_by(preg_uid) %>%
  summarise(
    miscarriage = max(miscarriage_event),
    CLUSTID = first(CLUSTID),
    survey_year_f = first(survey_year_f),
    analysis_weight = first(analysis_weight),
    across(all_of(available_factor_terms), first),
    across(
      all_of(names(exposure_candidates)),
      ~ mean(.x[gest_month <= 3], na.rm = TRUE),
      .names = "{.col}_early_pregnancy"
    ),
    .groups = "drop"
  )

if ("strata_f" %in% names(DT_primary)) {
  pregnancy_summary <- DT_primary %>%
    distinct(preg_uid, strata_f) %>%
    left_join(pregnancy_summary, by = "preg_uid")
}

if ("strata_f" %in% names(pregnancy_summary)) {
  descriptive_design <- svydesign(
    ids = ~CLUSTID,
    strata = ~strata_f,
    weights = ~analysis_weight,
    data = pregnancy_summary,
    nest = TRUE
  )
} else {
  descriptive_design <- svydesign(
    ids = ~CLUSTID,
    weights = ~analysis_weight,
    data = pregnancy_summary,
    nest = TRUE
  )
}

miscarriage_prevalence <- svymean(~miscarriage, descriptive_design, na.rm = TRUE)
weighted_outcome <- tibble(
  estimate = as.numeric(coef(miscarriage_prevalence)),
  standard_error = as.numeric(SE(miscarriage_prevalence)),
  lower_95 = estimate - 1.96 * standard_error,
  upper_95 = estimate + 1.96 * standard_error
)
write_csv(weighted_outcome, file.path(tables_dir, "Weighted_miscarriage_proportion.csv"))

weighted_exposure_summary <- map_dfr(names(exposure_candidates), function(exposure) {
  variable <- paste0(exposure, "_early_pregnancy")
  mean_result <- svymean(as.formula(paste0("~", variable)), descriptive_design, na.rm = TRUE)
  quantile_result <- svyquantile(
    as.formula(paste0("~", variable)),
    descriptive_design,
    quantiles = c(0.25, 0.50, 0.75),
    ci = FALSE,
    na.rm = TRUE
  )

  tibble(
    exposure = exposure,
    weighted_mean = as.numeric(coef(mean_result)),
    weighted_SE = as.numeric(SE(mean_result)),
    weighted_p25 = as.numeric(quantile_result)[1],
    weighted_median = as.numeric(quantile_result)[2],
    weighted_p75 = as.numeric(quantile_result)[3]
  )
})
write_csv(weighted_exposure_summary, file.path(tables_dir, "Weighted_exposure_summary.csv"))


# ---- 11. Publication-quality exposure-response figures ---------------------

exposure_labels <- c(
  HI = "Heat Index (degrees C)",
  WBT = "Wet-bulb globe temperature (degrees C)",
  CO = "Carbon monoxide (ppbv)",
  PM = expression(PM[2.5] ~ (mu * g/m^3)),
  EVI = "Enhanced Vegetation Index"
)

wave_curves <- imap(single_wave_models, function(model, exposure) {
  predict_hr_curve(
    model, DT_linear, exposure, paste0("Spline_M1_wave_", exposure)
  )
})

wave_panels <- imap(wave_curves, function(curve, exposure) {
  make_curve_panel(
    curve, DT_linear, exposure, exposure_labels[[exposure]],
    paste0(exposure, ": wave-adjusted")
  )
})

figure_2 <- wrap_plots(wave_panels, ncol = 2) +
  plot_annotation(
    title = "Monthly environmental exposure and reported miscarriage hazard",
    subtitle = "Survey-wave adjusted discrete-time models; curves restricted to the 1st-99th percentiles"
  )
save_high_resolution(figure_2, "Figure_2_Wave_Adjusted_Curves")

adjusted_curves <- imap(single_adjusted_models, function(model, exposure) {
  predict_hr_curve(
    model, DT_linear, exposure, paste0("Spline_M2_adjusted_", exposure)
  )
})

adjusted_panels <- imap(adjusted_curves, function(curve, exposure) {
  make_curve_panel(
    curve, DT_linear, exposure, exposure_labels[[exposure]],
    paste0(exposure, ": confounder-adjusted")
  )
})

figure_3 <- wrap_plots(adjusted_panels, ncol = 2) +
  plot_annotation(
    title = "Adjusted monthly environmental exposure-response relationships",
    subtitle = "Discrete-time models with reproductive, socioeconomic, survey-design, mother, and cluster terms"
  )
save_high_resolution(figure_3, "Figure_3_Adjusted_Curves")

make_full_model_figure <- function(model, exposures, model_name, stem) {
  panels <- map(exposures, function(exposure) {
    curve <- predict_hr_curve(model, DT_linear, exposure, model_name)
    make_curve_panel(
      curve, DT_linear, exposure, exposure_labels[[exposure]], exposure
    )
  })

  combined <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = paste0("Mutually adjusted exposure-response relationships: ", model_name),
      subtitle = "Hazard ratios are contrasts with the exposure-specific median; rugs show observed support"
    )

  save_high_resolution(combined, stem)
}

make_full_model_figure(
  full_hi_model,
  c("HI", "CO", "PM", "EVI"),
  "Spline GLMM HI model family",
  "Figure_4A_MultiExposure_HI_Curves"
)

make_full_model_figure(
  full_wbt_model,
  c("WBT", "CO", "PM", "EVI"),
  "Spline GLMM WBGT model family",
  "Figure_4B_MultiExposure_WBGT_Curves"
)


# ---- 12. Reproducibility record --------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(results_dir, "SessionInfo_02.txt"))

message(
  "Script 02 completed. Results are in: ", results_dir,
  "\nPrimary models: ",
  file.path(models_dir, "Primary_discrete_time_spline_GLMM_models.rds"),
  "\nIndividual checkpoints: ", checkpoint_dir
)
