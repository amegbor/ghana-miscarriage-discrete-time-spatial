# GitHub copy: portable paths; restricted data are not distributed.
# =============================================================================
# 05_Ghana_Primary_Model_Diagnostics_Extraction.R
#
# Purpose
# -------
# Recover the fitted primary linear glmmTMB models from saved checkpoint RDS
# files, document which checkpoint version was selected, verify that each model
# uses the final manuscript specification, extract fit/convergence diagnostics,
# and combine the results with the already-exported primary nonlinear model-fit
# table.
#
# This script DOES NOT fit or refit any statistical models.
#
# Final manuscript specification checked here
# ------------------------------------------------
#   * Mother-level random intercept included
#   * DHS-cluster random intercept included
#   * Region random intercept excluded
#   * Cooking fuel excluded
#   * EVI expressed per 0.1-unit increase in the linear EVI models
#
# Selection rule
# --------------
# When several checkpoint versions exist for the same model, the script uses
# only the corrected Script 02 checkpoint directory; legacy folders are
# supplied in Section 2.
# =============================================================================


# ---- 0. Clean session and load packages ------------------------------------

rm(list = ls())
gc(reset = TRUE)

required_packages <- c(
  "dplyr",
  "purrr",
  "tibble",
  "readr",
  "stringr",
  "glmmTMB"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
  library(tibble)
  library(readr)
  library(stringr)
  library(glmmTMB)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}


# ---- 1. Project paths -------------------------------------------------------

# Run from the repository root, or set GHANA_MISCARRIAGE_ROOT explicitly.
# Private inputs and outputs can live entirely outside the Git checkout.
project_dir <- normalizePath(
  Sys.getenv("GHANA_MISCARRIAGE_ROOT", unset = "."),
  winslash = "/", mustWork = TRUE
)

results_dir <- file.path(project_dir, "Results")
models_dir <- file.path(results_dir, "Models")
tables_dir <- file.path(results_dir, "Tables")
diagnostics_dir <- file.path(results_dir, "Diagnostics")

dir.create(
  tables_dir,
  recursive = FALSE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(models_dir)) {
  stop("Models directory not found: ", models_dir)
}

nonlinear_fit_file <- file.path(
  tables_dir,
  "Primary_spline_GLMM_model_fit.csv"
)

message(
  "Project directory: ", project_dir,
  "\nModels directory: ", models_dir,
  "\nTables directory: ", tables_dir
)


# Overrides are intentionally disabled: this audit uses only the current run.
checkpoint_overrides <- character(0)

# ---- 3. Inventory all primary linear checkpoints ---------------------------

expected_linear_models <- c(
  "Linear_unadjusted_HI",
  "Linear_adjusted_HI",
  "Linear_unadjusted_WBT",
  "Linear_adjusted_WBT",
  "Linear_unadjusted_CO",
  "Linear_adjusted_CO",
  "Linear_unadjusted_PM",
  "Linear_adjusted_PM",
  "Linear_unadjusted_EVI_per_0_1",
  "Linear_adjusted_EVI_per_0_1"
)

current_checkpoint_dir <- file.path(models_dir, "glmmTMB_checkpoints_no_cooking_v1")
if (!dir.exists(current_checkpoint_dir)) {
  stop("Run corrected Script 02 first; current checkpoint directory is missing.")
}
linear_checkpoint_files <- list.files(
  current_checkpoint_dir,
  pattern = paste0(
    "^Linear_(unadjusted|adjusted)_",
    "(HI|WBT|CO|PM|EVI_per_0_1)\\.rds$"
  ),
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = FALSE
)

if (length(linear_checkpoint_files) == 0L) {
  stop(
    "No primary linear checkpoint files were found under: ",
    models_dir
  )
}

safe_read_checkpoint <- function(checkpoint_file) {
  tryCatch(
    readRDS(checkpoint_file),
    error = function(e) structure(
      list(error_message = conditionMessage(e)),
      class = "checkpoint_read_error"
    )
  )
}

checkpoint_objects_for_validation <- map(
  linear_checkpoint_files,
  safe_read_checkpoint
)

linear_checkpoint_inventory <- tibble(
  inventory_index = seq_along(linear_checkpoint_files),
  model = tools::file_path_sans_ext(
    basename(linear_checkpoint_files)
  ),
  checkpoint_version = basename(
    dirname(linear_checkpoint_files)
  ),
  checkpoint_file = linear_checkpoint_files,
  modified = file.info(linear_checkpoint_files)$mtime,
  file_size_bytes = file.info(linear_checkpoint_files)$size,
  object_class = map_chr(
    checkpoint_objects_for_validation,
    ~ paste(class(.x), collapse = ", ")
  ),
  valid_glmmTMB = map_lgl(
    checkpoint_objects_for_validation,
    ~ inherits(.x, "glmmTMB")
  ),
  read_error = map_chr(
    checkpoint_objects_for_validation,
    function(object) {
      if (inherits(object, "checkpoint_read_error")) {
        object$error_message
      } else {
        NA_character_
      }
    }
  )
) %>%
  arrange(
    model,
    desc(modified)
  )

write_csv(
  linear_checkpoint_inventory,
  file.path(
    diagnostics_dir,
    "Primary_linear_checkpoint_complete_inventory.csv"
  )
)

message(
  "Identified ",
  nrow(linear_checkpoint_inventory),
  " checkpoint files for ",
  n_distinct(linear_checkpoint_inventory$model),
  " unique primary linear models; ",
  sum(linear_checkpoint_inventory$valid_glmmTMB),
  " files contain valid glmmTMB objects."
)

missing_expected_models <- setdiff(
  expected_linear_models,
  unique(linear_checkpoint_inventory$model)
)

if (length(missing_expected_models) > 0L) {
  stop(
    "The following expected models have no checkpoint: ",
    paste(missing_expected_models, collapse = ", ")
  )
}

models_without_valid_checkpoints <- expected_linear_models[
  !expected_linear_models %in%
    linear_checkpoint_inventory$model[
      linear_checkpoint_inventory$valid_glmmTMB
    ]
]

if (length(models_without_valid_checkpoints) > 0L) {
  stop(
    "No valid glmmTMB checkpoint was found for: ",
    paste(models_without_valid_checkpoints, collapse = ", "),
    ". Inspect Primary_linear_checkpoint_complete_inventory.csv."
  )
}


# ---- 4. Select one checkpoint for each model -------------------------------

selected_linear_files <- linear_checkpoint_inventory %>%
  filter(valid_glmmTMB) %>%
  group_by(model) %>%
  slice_max(
    order_by = modified,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

if (length(checkpoint_overrides) > 0L) {
  invalid_override_names <- setdiff(
    names(checkpoint_overrides),
    expected_linear_models
  )

  if (length(invalid_override_names) > 0L) {
    stop(
      "Unrecognised checkpoint override names: ",
      paste(invalid_override_names, collapse = ", ")
    )
  }

  missing_override_files <- checkpoint_overrides[
    !file.exists(checkpoint_overrides)
  ]

  if (length(missing_override_files) > 0L) {
    stop(
      "The following override files do not exist: ",
      paste(missing_override_files, collapse = ", ")
    )
  }

  for (model_name in names(checkpoint_overrides)) {
    override_file <- unname(
      checkpoint_overrides[[model_name]]
    )

    selected_linear_files <- selected_linear_files %>%
      filter(model != model_name) %>%
      bind_rows(
        tibble(
          inventory_index = NA_integer_,
          model = model_name,
          checkpoint_version = basename(dirname(override_file)),
          checkpoint_file = override_file,
          modified = file.info(override_file)$mtime,
          file_size_bytes = file.info(override_file)$size,
          object_class = paste(
            class(safe_read_checkpoint(override_file)),
            collapse = ", "
          ),
          valid_glmmTMB = inherits(
            safe_read_checkpoint(override_file),
            "glmmTMB"
          ),
          read_error = NA_character_
        )
      )
  }
}

if (any(!selected_linear_files$valid_glmmTMB)) {
  stop(
    "A manually overridden checkpoint is not a valid glmmTMB object."
  )
}

selected_linear_files <- selected_linear_files %>%
  mutate(
    model_order = match(model, expected_linear_models),
    selection_method = if_else(
      model %in% names(checkpoint_overrides),
      "Manual override",
      "Current corrected Script 02 checkpoint"
    )
  ) %>%
  arrange(model_order)

if (
  nrow(selected_linear_files) != length(expected_linear_models) ||
  anyDuplicated(selected_linear_files$model)
) {
  stop(
    "Checkpoint selection did not produce exactly one file for each model."
  )
}

write_csv(
  selected_linear_files,
  file.path(
    diagnostics_dir,
    "Primary_linear_checkpoint_selected_files.csv"
  )
)

print(
  selected_linear_files %>%
    select(
      model,
      checkpoint_version,
      modified,
      selection_method,
      checkpoint_file
    ),
  n = Inf,
  width = Inf
)


# ---- 5. Load and validate selected model objects ----------------------------

# Reuse the exact object that passed validation in Section 3. This prevents a
# stale workspace object or a file changed between inventory and extraction
# from creating a disagreement between validation and model loading. Manual
# overrides, if supplied, are read directly because they are not members of the
# original inventory list.
linear_models <- setNames(
  map2(
    selected_linear_files$inventory_index,
    selected_linear_files$checkpoint_file,
    function(inventory_index, checkpoint_file) {
      if (!is.na(inventory_index)) {
        checkpoint_objects_for_validation[[inventory_index]]
      } else {
        safe_read_checkpoint(checkpoint_file)
      }
    }
  ),
  selected_linear_files$model
)

model_class_check <- imap_dfr(
  linear_models,
  function(fitted_model, model_name) {
    tibble(
      model = model_name,
      object_class = paste(class(fitted_model), collapse = ", "),
      valid_glmmTMB = inherits(fitted_model, "glmmTMB")
    )
  }
)

if (any(!model_class_check$valid_glmmTMB)) {
  print(model_class_check, n = Inf)
  stop("One or more selected checkpoints are not valid glmmTMB models.")
}


# ---- 6. Audit final model formulas -----------------------------------------

linear_formula_audit <- imap_dfr(
  linear_models,
  function(fitted_model, model_name) {
    formula_text <- paste(
      deparse(formula(fitted_model)),
      collapse = " "
    )

    selected_row <- selected_linear_files %>%
      filter(model == model_name)

    is_evi_model <- grepl("_EVI_per_0_1$", model_name)

    tibble(
      model = model_name,
      checkpoint_version = selected_row$checkpoint_version,
      checkpoint_file = selected_row$checkpoint_file,
      valid_glmmTMB = inherits(fitted_model, "glmmTMB"),
      includes_mother_random_effect = grepl(
        "mother_uid",
        formula_text,
        fixed = TRUE
      ),
      includes_cluster_random_effect = grepl(
        "CLUSTID",
        formula_text,
        fixed = TRUE
      ),
      includes_region_random_effect = grepl(
        "region_uid",
        formula_text,
        fixed = TRUE
      ),
      includes_cooking_fuel = grepl(
        "clean_cooking_fuel|cook_fuel",
        formula_text,
        ignore.case = TRUE
      ),
      EVI_scaled_per_0_1 = if (is_evi_model) {
        grepl("EVI_per_0_1", formula_text, fixed = TRUE)
      } else {
        NA
      },
      final_specification_pass =
        includes_mother_random_effect &
        includes_cluster_random_effect &
        !includes_region_random_effect &
        !includes_cooking_fuel &
        if_else(
          is_evi_model,
          EVI_scaled_per_0_1,
          TRUE
        ),
      formula = formula_text
    )
  }
)

write_csv(
  linear_formula_audit,
  file.path(
    diagnostics_dir,
    "Primary_linear_model_formula_audit.csv"
  )
)

print(
  linear_formula_audit,
  n = Inf,
  width = Inf
)

if (any(!linear_formula_audit$final_specification_pass)) {
  warning(
    paste0(
      "At least one selected model failed the final-specification audit. ",
      "Inspect Primary_linear_model_formula_audit.csv before using the ",
      "diagnostic results in the manuscript."
    )
  )
}


# ---- 7. Safe diagnostic helper functions -----------------------------------

safe_numeric <- function(expression) {
  tryCatch(
    as.numeric(expression),
    error = function(e) NA_real_
  )
}

safe_integer <- function(expression) {
  tryCatch(
    as.integer(expression),
    error = function(e) NA_integer_
  )
}

safe_log_likelihood <- function(model) {
  tryCatch(
    as.numeric(logLik(model)),
    error = function(e) NA_real_
  )
}

safe_AIC <- function(model) {
  tryCatch(
    as.numeric(AIC(model)),
    error = function(e) NA_real_
  )
}

safe_BIC <- function(model) {
  tryCatch(
    as.numeric(BIC(model)),
    error = function(e) NA_real_
  )
}

safe_nobs <- function(model) {
  tryCatch(
    as.integer(stats::nobs(model)),
    error = function(e) NA_integer_
  )
}

safe_maximum_gradient <- function(model) {
  tryCatch(
    {
      gradient <- model$sdr$gradient.fixed

      if (is.null(gradient) || length(gradient) == 0L) {
        return(NA_real_)
      }

      if (all(!is.finite(gradient))) {
        return(NA_real_)
      }

      max(abs(gradient), na.rm = TRUE)
    },
    error = function(e) NA_real_
  )
}

safe_covariance_matrix <- function(model) {
  tryCatch(
    as.matrix(vcov(model)$cond),
    error = function(e) NULL
  )
}

safe_optimizer_message <- function(model) {
  tryCatch(
    {
      output <- model$fit$message
      if (is.null(output) || length(output) == 0L) {
        NA_character_
      } else {
        paste(as.character(output), collapse = "; ")
      }
    },
    error = function(e) NA_character_
  )
}


# ---- 8. Extract linear model fit and convergence diagnostics ----------------

# This is a screening threshold rather than a universal pass/fail rule. Models
# exceeding it should be reviewed with glmmTMB::diagnose() and their parameter
# estimates inspected for separation or boundary variance estimates.
GRADIENT_REVIEW_THRESHOLD <- 0.01

linear_model_diagnostics <- imap_dfr(
  linear_models,
  function(fitted_model, model_name) {
    covariance_matrix <- safe_covariance_matrix(fitted_model)

    formula_audit_row <- linear_formula_audit %>%
      filter(model == model_name)

    selected_row <- selected_linear_files %>%
      filter(model == model_name)

    analysis_n <- attr(fitted_model, "analysis_n") %||%
      safe_nobs(fitted_model)
    analysis_pregnancies <- attr(
      fitted_model,
      "analysis_pregnancies"
    ) %||% NA_integer_
    analysis_events <- attr(
      fitted_model,
      "analysis_events"
    ) %||% NA_integer_

    maximum_gradient <- safe_maximum_gradient(fitted_model)

    convergence_code <- tryCatch(
      as.integer(fitted_model$fit$convergence),
      error = function(e) NA_integer_
    )

    positive_definite_hessian <- tryCatch(
      as.logical(fitted_model$sdr$pdHess),
      error = function(e) NA
    )

    covariance_contains_NA <- if (is.null(covariance_matrix)) {
      NA
    } else {
      anyNA(covariance_matrix)
    }

    covariance_contains_nonfinite <- if (is.null(covariance_matrix)) {
      NA
    } else {
      any(!is.finite(covariance_matrix))
    }

    tibble(
      model = model_name,
      model_family = "Primary linear",
      adjustment = if_else(
        grepl("^Linear_unadjusted_", model_name),
        "Unadjusted",
        "Adjusted"
      ),
      exposure = sub(
        "^Linear_(unadjusted|adjusted)_",
        "",
        model_name
      ),
      person_months = as.integer(analysis_n),
      pregnancies = as.integer(analysis_pregnancies),
      miscarriage_events = as.integer(analysis_events),
      log_likelihood = safe_log_likelihood(fitted_model),
      AIC = safe_AIC(fitted_model),
      BIC = safe_BIC(fitted_model),
      convergence_code = convergence_code,
      optimizer_message = safe_optimizer_message(fitted_model),
      positive_definite_Hessian = positive_definite_hessian,
      maximum_absolute_gradient = maximum_gradient,
      gradient_requires_review = if_else(
        is.na(maximum_gradient),
        NA,
        maximum_gradient > GRADIENT_REVIEW_THRESHOLD
      ),
      covariance_contains_NA = covariance_contains_NA,
      covariance_contains_nonfinite = covariance_contains_nonfinite,
      final_specification_pass =
        formula_audit_row$final_specification_pass,
      checkpoint_version = selected_row$checkpoint_version,
      checkpoint_modified = selected_row$modified,
      checkpoint_file = selected_row$checkpoint_file
    )
  }
)

linear_model_diagnostics <- linear_model_diagnostics %>%
  mutate(
    convergence_diagnostic_pass =
      convergence_code == 0L &
      positive_definite_Hessian %in% TRUE &
      covariance_contains_NA %in% FALSE &
      covariance_contains_nonfinite %in% FALSE,
    overall_review_required =
      !final_specification_pass |
      !convergence_diagnostic_pass |
      gradient_requires_review %in% TRUE
  )

write_csv(
  linear_model_diagnostics,
  file.path(
    tables_dir,
    "Primary_linear_GLMM_model_fit_and_diagnostics.csv"
  )
)

print(
  linear_model_diagnostics,
  n = Inf,
  width = Inf
)


# ---- 9. Create a concise manuscript-table source for linear models ----------

linear_table_source <- linear_model_diagnostics %>%
  transmute(
    panel = "Panel A. Primary linear models",
    model,
    exposure,
    adjustment,
    person_months,
    pregnancies,
    miscarriage_events,
    log_likelihood,
    AIC,
    BIC,
    convergence_code,
    positive_definite_Hessian,
    maximum_absolute_gradient,
    covariance_valid =
      covariance_contains_NA %in% FALSE &
      covariance_contains_nonfinite %in% FALSE,
    final_specification_pass,
    overall_review_required
  )

write_csv(
  linear_table_source,
  file.path(
    tables_dir,
    "Supplementary_Table_Primary_Linear_Diagnostics_Source.csv"
  )
)


# ---- 10. Import and standardise nonlinear model-fit diagnostics -------------

nonlinear_models_file <- file.path(models_dir, "Primary_discrete_time_spline_GLMM_models.rds")
if (!file.exists(nonlinear_models_file)) {
  stop("Run corrected Script 02 first; nonlinear models are missing.")
}
primary_data_file <- file.path(models_dir, "Ghana_DT_AnalysisReady_PrimaryCC.rds")
if (!file.exists(primary_data_file) ||
    !isTRUE(attr(readRDS(primary_data_file), "ghana_cooking_fuel_excluded"))) {
  stop("Run corrected Script 02 first; the primary dataset is stale or missing.")
}
nonlinear_models <- readRDS(nonlinear_models_file)
nonlinear_formula_audit <- imap_dfr(nonlinear_models, function(model, model_name) {
  text <- paste(deparse(formula(model)), collapse = " ")
  tibble(
    model = model_name,
    formula = text,
    final_specification_pass =
      grepl("mother_uid", text, fixed = TRUE) &
      grepl("CLUSTID", text, fixed = TRUE) &
      !grepl("region_uid|cook|fuel", text, ignore.case = TRUE)
  )
})
write_csv(nonlinear_formula_audit, file.path(diagnostics_dir, "Primary_nonlinear_formula_audit.csv"))

if (file.exists(nonlinear_fit_file)) {
  nonlinear_model_fit <- read_csv(
    nonlinear_fit_file,
    show_col_types = FALSE
  )

  required_nonlinear_columns <- c(
    "model",
    "person_months",
    "pregnancies",
    "events",
    "log_likelihood",
    "AIC",
    "BIC",
    "convergence_code",
    "positive_definite_Hessian"
  )

  missing_nonlinear_columns <- setdiff(
    required_nonlinear_columns,
    names(nonlinear_model_fit)
  )

  if (length(missing_nonlinear_columns) > 0L) {
    stop(
      "Primary_spline_GLMM_model_fit.csv is missing: ",
      paste(missing_nonlinear_columns, collapse = ", ")
    )
  }

  nonlinear_table_source <- nonlinear_model_fit %>%
    left_join(nonlinear_formula_audit %>% select(model, final_specification_pass), by = "model") %>%
    mutate(
      panel = "Panel B. Primary nonlinear models",
      exposure = case_when(
        grepl("WBT", model, fixed = TRUE) ~ "WBGT",
        grepl("HI", model, fixed = TRUE) ~ "Heat Index",
        TRUE ~ NA_character_
      ),
      adjustment = case_when(
        grepl("M1_wave", model, fixed = TRUE) ~ "Survey-wave adjusted",
        grepl("M2_adjusted", model, fixed = TRUE) ~ "Confounder adjusted",
        grepl("M3_multi", model, fixed = TRUE) ~
          "Mutually adjusted multi-exposure",
        TRUE ~ NA_character_
      ),
      miscarriage_events = events,
      maximum_absolute_gradient = NA_real_,
      covariance_valid = NA,
      overall_review_required =
        is.na(final_specification_pass) | !final_specification_pass |
        convergence_code != 0L |
        !positive_definite_Hessian
    ) %>%
    select(
      panel,
      model,
      exposure,
      adjustment,
      person_months,
      pregnancies,
      miscarriage_events,
      log_likelihood,
      AIC,
      BIC,
      convergence_code,
      positive_definite_Hessian,
      maximum_absolute_gradient,
      covariance_valid,
      final_specification_pass,
      overall_review_required
    )

  write_csv(
    nonlinear_table_source,
    file.path(
      tables_dir,
      "Supplementary_Table_Primary_Nonlinear_Diagnostics_Source.csv"
    )
  )

  primary_model_diagnostics_combined <- bind_rows(
    linear_table_source,
    nonlinear_table_source
  )

  write_csv(
    primary_model_diagnostics_combined,
    file.path(
      tables_dir,
      "Supplementary_Table_Primary_Model_Diagnostics_Combined.csv"
    )
  )

  message(
    "Combined primary linear and nonlinear diagnostic source table saved."
  )
} else {
  warning(
    "Nonlinear fit file was not found: ",
    nonlinear_fit_file,
    ". Linear diagnostic outputs were still created."
  )
}


# ---- 11. Produce an audit summary ------------------------------------------

audit_summary <- tibble(
  item = c(
    "Checkpoint files inventoried",
    "Unique primary linear models found",
    "Models selected",
    "Valid glmmTMB objects",
    "Models passing final formula specification",
    "Models with convergence code 0",
    "Models with positive-definite Hessian",
    "Models with valid covariance matrix",
    paste0(
      "Models with maximum absolute gradient > ",
      GRADIENT_REVIEW_THRESHOLD
    ),
    "Models requiring review"
  ),
  n = c(
    nrow(linear_checkpoint_inventory),
    n_distinct(linear_checkpoint_inventory$model),
    nrow(selected_linear_files),
    sum(model_class_check$valid_glmmTMB),
    sum(linear_formula_audit$final_specification_pass),
    sum(linear_model_diagnostics$convergence_code == 0L, na.rm = TRUE),
    sum(
      linear_model_diagnostics$positive_definite_Hessian %in% TRUE,
      na.rm = TRUE
    ),
    sum(
      linear_model_diagnostics$covariance_contains_NA %in% FALSE &
        linear_model_diagnostics$covariance_contains_nonfinite %in% FALSE,
      na.rm = TRUE
    ),
    sum(
      linear_model_diagnostics$gradient_requires_review %in% TRUE,
      na.rm = TRUE
    ),
    sum(
      linear_model_diagnostics$overall_review_required %in% TRUE,
      na.rm = TRUE
    )
  )
)

write_csv(
  audit_summary,
  file.path(
    diagnostics_dir,
    "Primary_linear_model_diagnostics_audit_summary.csv"
  )
)

print(audit_summary, n = Inf)


# ---- 12. Save reproducibility information ----------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    diagnostics_dir,
    "SessionInfo_Primary_Model_Diagnostics_Extraction.txt"
  )
)

message(
  "Diagnostic extraction completed without refitting models.",
  "\nMain linear diagnostic file: ",
  file.path(
    tables_dir,
    "Primary_linear_GLMM_model_fit_and_diagnostics.csv"
  ),
  "\nFormula audit: ",
  file.path(
    diagnostics_dir,
    "Primary_linear_model_formula_audit.csv"
  ),
  "\nSelected checkpoint audit: ",
  file.path(
    diagnostics_dir,
    "Primary_linear_checkpoint_selected_files.csv"
  )
)
