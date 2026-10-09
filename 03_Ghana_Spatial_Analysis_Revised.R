# GitHub copy: portable paths; restricted data are not distributed.
# =============================================================================
# 03_Ghana_Spatial_Analysis_Revised.R
# Joint Bayesian discrete-time survival and INLA-SPDE spatial analysis
#
# Outputs
#   * Refit and checkpoint HI and WBGT spatial models
#   * Main WBGT posterior-mean fields with Ghana regional boundaries
#   * Supplementary HI posterior-mean fields with regional boundaries
#   * Fully adjusted HI-WBGT comparison and uncertainty figures
#   * Separate, clearly labelled HI and WBGT Leaflet HTML maps
#
# The Leaflet maps use explicit OpenStreetMap raster tiles. OpenStreetMap does
# not require an API key. The saved HTML files require internet access when
# opened because map tiles are retrieved at viewing time.
# =============================================================================

rm(list = ls())
gc(reset = TRUE)


# ---- 0. Packages ------------------------------------------------------------

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "sf", "ggplot2",
  "patchwork", "INLA", "fmesher", "Matrix", "geodata", "leaflet",
  "htmlwidgets", "htmltools", "viridisLite", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Script 03: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(INLA)
  library(fmesher)
  library(Matrix)
  library(leaflet)
  library(htmlwidgets)
})

set.seed(20261001)


# ---- 1. Paths and analysis settings ----------------------------------------

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
spatial_dir <- file.path(results_dir, "Spatial")
checkpoint_dir <- file.path(models_dir, "INLA_SPDE_checkpoints_final_v1")
gadm_cache_dir <- file.path(project_dir, "GADM")

walk(
  c(
    results_dir, tables_dir, figures_dir, models_dir, spatial_dir,
    checkpoint_dir, gadm_cache_dir
  ),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

# Require the rebuilt primary complete-case dataset from corrected Script 02.
input_file <- file.path(models_dir, "Ghana_DT_AnalysisReady_PrimaryCC.rds")
if (!file.exists(input_file)) {
  stop("Run corrected Script 02 first. Missing dataset: ", input_file)
}

# This revision intentionally refits all spatial models. Set the environment
# variable GHANA_REFIT_SPATIAL_MODELS=false only when checkpoints should be used.
REFIT_MODELS <- tolower(
  Sys.getenv("GHANA_REFIT_SPATIAL_MODELS", unset = "true")
) %in% c("true", "1", "yes")

thread_text <- Sys.getenv("GHANA_N_THREADS", unset = "4")
N_THREADS <- suppressWarnings(as.integer(thread_text))
if (is.na(N_THREADS) || N_THREADS < 1L) N_THREADS <- 4L

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.finite(detected_cores) && detected_cores > 0L) {
  N_THREADS <- min(N_THREADS, as.integer(detected_cores))
}

options(mc.cores = N_THREADS)
INLA::inla.setOption(num.threads = N_THREADS)

PROJECTED_CRS <- 32630
MESH_MAX_EDGE_INNER <- 25000
MESH_MAX_EDGE_OUTER <- 75000
MESH_CUTOFF <- 5000
MESH_OFFSET_INNER <- 25000
MESH_OFFSET_OUTER <- 100000
SPATIAL_RANGE0 <- 50000
SPATIAL_RANGE_PROB <- 0.50
SPATIAL_SD0 <- 1
SPATIAL_SD_PROB <- 0.01
IID_SD0 <- 1
IID_SD_PROB <- 0.01
MAP_CELL_SIZE <- 5000

posterior_sample_text <- Sys.getenv(
  "GHANA_POSTERIOR_SAMPLES",
  unset = "500"
)
N_POSTERIOR_SAMPLES <- suppressWarnings(as.integer(posterior_sample_text))
if (is.na(N_POSTERIOR_SAMPLES) || N_POSTERIOR_SAMPLES < 2L) {
  N_POSTERIOR_SAMPLES <- 500L
}

message(
  "Input: ", input_file,
  "\nRefit models: ", REFIT_MODELS,
  "\nINLA threads: ", N_THREADS,
  "\nPosterior samples per uncertainty map: ", N_POSTERIOR_SAMPLES
)


# ---- 2. General helper functions -------------------------------------------

formula_from_terms <- function(fixed_terms) {
  fixed_string <- paste(c("0 + Intercept", fixed_terms), collapse = " + ")

  as.formula(
    paste0(
      "miscarriage_event ~ ", fixed_string,
      " + f(mother_idx, model = 'iid', hyper = iid_hyper)",
      " + f(cluster_idx, model = 'iid', hyper = iid_hyper)",
      " + f(spatial_field, model = spde_model)"
    )
  )
}

fit_spatial_model <- function(
    model_name,
    fixed_terms,
    data,
    A_observation,
    spde_model,
    force_refit = REFIT_MODELS) {

  checkpoint_file <- file.path(
    checkpoint_dir,
    paste0(model_name, ".rds")
  )

  if (!force_refit && file.exists(checkpoint_file)) {
    message("Loading checkpoint: ", checkpoint_file)
    return(readRDS(checkpoint_file))
  }

  message(
    "Fitting ", model_name, ": ",
    format(nrow(data), big.mark = ","), " person-months; ",
    format(n_distinct(data$preg_uid), big.mark = ","), " pregnancies; ",
    sum(data$miscarriage_event, na.rm = TRUE), " events."
  )

  fixed_variable_names <- unique(
    c(fixed_terms, "analysis_weight", "mother_idx", "cluster_idx")
  )

  fixed_effects <- as.data.frame(
    data[, fixed_variable_names, drop = FALSE]
  )
  fixed_effects$Intercept <- 1

  spatial_index <- INLA::inla.spde.make.index(
    name = "spatial_field",
    n.spde = spde_model$n.spde
  )

  model_stack <- INLA::inla.stack(
    data = list(miscarriage_event = data$miscarriage_event),
    A = list(1, A_observation),
    effects = list(fixed_effects, spatial_index),
    tag = model_name
  )

  stack_data <- INLA::inla.stack.data(model_stack)

  fit <- INLA::inla(
    formula = formula_from_terms(fixed_terms),
    family = "binomial",
    Ntrials = 1,
    data = stack_data,
    weights = stack_data$analysis_weight,
    control.family = list(link = "cloglog"),
    control.predictor = list(
      A = INLA::inla.stack.A(model_stack),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      config = TRUE,
      return.marginals.predictor = FALSE
    ),
    control.inla = list(
      strategy = "adaptive",
      int.strategy = "ccd"
    ),
    verbose = TRUE
  )

  attr(fit, "model_name") <- model_name
  attr(fit, "fixed_terms") <- fixed_terms
  saveRDS(fit, checkpoint_file, compress = FALSE)
  message("Saved checkpoint: ", checkpoint_file)
  fit
}

extract_fixed_hr <- function(fit, model_name) {
  out <- as.data.frame(fit$summary.fixed)
  out$term <- rownames(out)

  tibble(
    model = model_name,
    term = out$term,
    posterior_mean_log_hr = out$mean,
    posterior_sd = out$sd,
    HR = exp(out$mean),
    lower_95 = exp(out$`0.025quant`),
    median_HR = exp(out$`0.5quant`),
    upper_95 = exp(out$`0.975quant`)
  )
}

extract_fit_criteria <- function(fit, model_name) {
  tibble(
    model = model_name,
    DIC = fit$dic$dic,
    effective_parameters_DIC = fit$dic$p.eff,
    WAIC = fit$waic$waic,
    effective_parameters_WAIC = fit$waic$p.eff,
    mean_negative_log_CPO = if (!is.null(fit$cpo$cpo)) {
      mean(-log(fit$cpo$cpo), na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

extract_hyperparameters <- function(fit, model_name) {
  out <- as.data.frame(fit$summary.hyperpar)
  out$parameter <- rownames(out)
  out$model <- model_name
  as_tibble(out) %>% relocate(model, parameter)
}

project_field_mean <- function(fit, A_prediction) {
  as.numeric(
    A_prediction %*% fit$summary.random$spatial_field$mean
  )
}

sample_projected_field <- function(fit, A_prediction, n_samples) {
  samples <- INLA::inla.posterior.sample(
    n = n_samples,
    result = fit
  )

  latent_names <- rownames(samples[[1]]$latent)
  field_rows <- grep("^spatial_field", latent_names)
  expected_n <- nrow(fit$summary.random$spatial_field)

  if (length(field_rows) != expected_n) {
    stop(
      "Spatial field could not be identified in posterior samples. Found ",
      length(field_rows), "; expected ", expected_n, "."
    )
  }

  field_draws <- vapply(
    samples,
    function(draw) as.numeric(draw$latent[field_rows, 1]),
    numeric(length(field_rows))
  )

  projected_draws <- as.matrix(A_prediction %*% field_draws)

  tibble(
    posterior_mean = rowMeans(projected_draws),
    posterior_sd = apply(projected_draws, 1, sd),
    probability_positive = rowMeans(projected_draws > 0),
    lower_95 = apply(projected_draws, 1, quantile, probs = 0.025),
    upper_95 = apply(projected_draws, 1, quantile, probs = 0.975)
  )
}

save_high_resolution <- function(plot, stem, width = 11, height = 8.5) {
  ggsave(
    file.path(figures_dir, paste0(stem, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  ggsave(
    file.path(figures_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}


# ---- 3. Read GADM Ghana boundaries -----------------------------------------

ghana_regions_file <- file.path(
  gadm_cache_dir,
  "GADM_4.1_Ghana_ADM1.gpkg"
)

if (file.exists(ghana_regions_file)) {
  ghana_regions <- st_read(ghana_regions_file, quiet = TRUE)
} else {
  message("Downloading Ghana level-1 boundaries from GADM...")

  ghana_regions <- geodata::gadm(
    country = "GHA",
    level = 1,
    path = gadm_cache_dir,
    version = "4.1",
    resolution = 1
  ) %>%
    st_as_sf() %>%
    st_make_valid()

  st_write(
    ghana_regions,
    ghana_regions_file,
    layer = "Ghana_ADM1",
    delete_dsn = TRUE,
    quiet = TRUE
  )
}

region_name_candidates <- c(
  "NAME_1", "VARNAME_1", "NL_NAME_1", "ADM1_NAME", "ADM1NAME",
  "REGION_NAME", "REGION", "Region", "region", "NAME"
)

region_name_variable <- region_name_candidates[
  region_name_candidates %in% names(ghana_regions)
][1]

if (is.na(region_name_variable)) {
  stop(
    "Could not identify the GADM regional-name field. Available fields: ",
    paste(names(ghana_regions), collapse = ", ")
  )
}

ghana_regions_utm <- ghana_regions %>%
  st_make_valid() %>%
  st_transform(PROJECTED_CRS) %>%
  mutate(region_label = as.character(.data[[region_name_variable]])) %>%
  select(region_label) %>%
  group_by(region_label) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_make_valid()

ghana_boundary_utm <- ghana_regions_utm %>%
  summarise(do_union = TRUE) %>%
  st_make_valid()

if (nrow(ghana_regions_utm) != 16L) {
  warning(
    "Expected 16 Ghana regions but obtained ",
    nrow(ghana_regions_utm), "."
  )
}

region_label_points <- suppressWarnings(
  st_point_on_surface(ghana_regions_utm)
)
region_label_xy <- st_coordinates(region_label_points)

region_labels <- region_label_points %>%
  mutate(
    label_x = region_label_xy[, 1],
    label_y = region_label_xy[, 2]
  )


# ---- 4. Prepare spatial analysis data --------------------------------------

DT <- readRDS(input_file)
if (!isTRUE(attr(DT, "ghana_analysis_ready_final"))) {
  stop("Run corrected Script 02 first; this dataset lacks the final analysis-ready compatibility marker.")
}
DT <- as.data.frame(DT)

required_coordinates <- c("LONGNUM", "LATNUM")
missing_coordinates <- setdiff(required_coordinates, names(DT))

if (length(missing_coordinates) > 0L) {
  stop(
    "Spatial coordinates are missing: ",
    paste(missing_coordinates, collapse = ", ")
  )
}

DT_spatial <- DT %>%
  filter(
    is.finite(LONGNUM),
    is.finite(LATNUM),
    !(LONGNUM == 0 & LATNUM == 0)
  ) %>%
  mutate(
    age_10y = (
      age_at_preg_start - median(age_at_preg_start, na.rm = TRUE)
    ) / 10,
    age_10y_sq = age_10y^2,
    EVI_per_0_1 = EVI / 0.1
  ) %>%
  droplevels()

# Current marital status is intentionally excluded from all
# final manuscript models. Residence is not added separately when strata_f
# already encodes the survey design/residence stratification.

cluster_coordinate_check <- DT_spatial %>%
  distinct(CLUSTID, LONGNUM, LATNUM) %>%
  count(CLUSTID, name = "n_coordinate_pairs")

if (any(cluster_coordinate_check$n_coordinate_pairs != 1L)) {
  inconsistent_coordinates <- DT_spatial %>%
    semi_join(
      filter(cluster_coordinate_check, n_coordinate_pairs != 1L),
      by = "CLUSTID"
    ) %>%
    distinct(CLUSTID, LONGNUM, LATNUM)

  write_csv(
    inconsistent_coordinates,
    file.path(tables_dir, "Clusters_with_inconsistent_coordinates.csv")
  )

  stop(
    "One or more DHS clusters have multiple coordinate pairs. See the QC file."
  )
}

cluster_coordinates <- DT_spatial %>%
  distinct(CLUSTID, LONGNUM, LATNUM)

cluster_points <- st_as_sf(
  cluster_coordinates,
  coords = c("LONGNUM", "LATNUM"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(PROJECTED_CRS)

inside_country <- lengths(
  st_covered_by(cluster_points, ghana_boundary_utm)
) > 0

if (any(!inside_country)) {
  write_csv(
    st_drop_geometry(cluster_points[!inside_country, ]),
    file.path(tables_dir, "Clusters_outside_Ghana_boundary.csv")
  )
}

cluster_points <- cluster_points[inside_country, ]
valid_clusters <- as.character(cluster_points$CLUSTID)

DT_spatial <- DT_spatial %>%
  filter(as.character(CLUSTID) %in% valid_clusters) %>%
  droplevels()

xy_lookup <- cluster_points %>%
  mutate(
    x_utm = st_coordinates(.)[, 1],
    y_utm = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(CLUSTID, x_utm, y_utm)

n_before_join <- nrow(DT_spatial)

DT_spatial <- DT_spatial %>%
  left_join(xy_lookup, by = "CLUSTID") %>%
  mutate(
    mother_idx = as.integer(factor(mother_uid)),
    cluster_idx = as.integer(factor(CLUSTID)),
    analysis_weight = analysis_weight / mean(analysis_weight, na.rm = TRUE)
  )

if (nrow(DT_spatial) != n_before_join) {
  stop("The coordinate join changed the person-month row count.")
}

spatial_qc <- tibble(
  item = c(
    "Analysis person-months",
    "Pregnancies",
    "Mothers",
    "DHS clusters",
    "Miscarriage events"
  ),
  n = c(
    nrow(DT_spatial),
    n_distinct(DT_spatial$preg_uid),
    n_distinct(DT_spatial$mother_uid),
    n_distinct(DT_spatial$CLUSTID),
    sum(DT_spatial$miscarriage_event, na.rm = TRUE)
  )
)

write_csv(
  spatial_qc,
  file.path(tables_dir, "Spatial_sample_QC.csv")
)

observation_coordinates <- as.matrix(
  DT_spatial[, c("x_utm", "y_utm")]
)
mesh_coordinates <- st_coordinates(cluster_points)


# ---- 5. Construct SPDE mesh -------------------------------------------------

mesh <- fmesher::fm_mesh_2d_inla(
  loc = mesh_coordinates,
  boundary = ghana_boundary_utm,
  max.edge = c(MESH_MAX_EDGE_INNER, MESH_MAX_EDGE_OUTER),
  cutoff = MESH_CUTOFF,
  offset = c(MESH_OFFSET_INNER, MESH_OFFSET_OUTER),
  crs = fmesher::fm_crs(ghana_boundary_utm)
)

spde_model <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  constr = TRUE,
  prior.range = c(SPATIAL_RANGE0, SPATIAL_RANGE_PROB),
  prior.sigma = c(SPATIAL_SD0, SPATIAL_SD_PROB)
)

iid_hyper <- list(
  prec = list(
    prior = "pc.prec",
    param = c(IID_SD0, IID_SD_PROB)
  )
)

A_observation <- INLA::inla.spde.make.A(
  mesh = mesh,
  loc = observation_coordinates
)

mesh_summary <- tibble(
  parameter = c(
    "mesh_vertices", "mesh_triangles", "inner_max_edge_m",
    "outer_max_edge_m", "cutoff_m", "inner_offset_m", "outer_offset_m",
    "projected_crs_epsg", "spatial_prior_range_m",
    "spatial_prior_P_range_below", "spatial_prior_sd",
    "spatial_prior_P_sd_above"
  ),
  value = c(
    mesh$n, nrow(mesh$graph$tv), MESH_MAX_EDGE_INNER,
    MESH_MAX_EDGE_OUTER, MESH_CUTOFF, MESH_OFFSET_INNER,
    MESH_OFFSET_OUTER, PROJECTED_CRS, SPATIAL_RANGE0,
    SPATIAL_RANGE_PROB, SPATIAL_SD0, SPATIAL_SD_PROB
  )
)

write_csv(
  mesh_summary,
  file.path(tables_dir, "SPDE_mesh_and_prior_specification.csv")
)

tiff(
  file.path(figures_dir, "Supplementary_Figure_S1_SPDE_Mesh.tiff"),
  width = 10,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw"
)
plot(mesh, asp = 1, main = "SPDE mesh and Ghana DHS cluster locations")
plot(st_geometry(ghana_boundary_utm), add = TRUE, border = "black", lwd = 1.2)
points(
  mesh_coordinates,
  pch = 16,
  cex = 0.35,
  col = adjustcolor("#B22222", 0.65)
)
dev.off()


# ---- 6. Define and refit spatial models ------------------------------------

if (
  "strata_f" %in% names(DT_spatial) &&
  n_distinct(DT_spatial$strata_f, na.rm = TRUE) > 1L
) {
  design_terms <- c("gest_month_f", "strata_f")
} else {
  design_terms <- c("gest_month_f", "survey_year_f")
}

social_candidates <- c(
  "age_10y", "age_10y_sq", "parity_cat", "previous_miscarriage_f",
  "interpregnancy_interval_f", "education_f", "wealth_f",
  "no_weekly_media_f", "health_insurance_f", "sanitation_f",
  "water_improved_f", "internet_ever_f"
)

social_terms <- intersect(
  social_candidates,
  names(DT_spatial)
)

environmental_hi_terms <- c("HI", "CO", "PM", "EVI_per_0_1")
environmental_wbgt_terms <- c("WBT", "CO", "PM", "EVI_per_0_1")

required_model_variables <- unique(
  c(
    "miscarriage_event", "analysis_weight", "mother_uid", "CLUSTID",
    "preg_uid", design_terms, environmental_hi_terms,
    environmental_wbgt_terms, social_terms
  )
)

missing_model_variables <- setdiff(
  required_model_variables,
  names(DT_spatial)
)

if (length(missing_model_variables) > 0L) {
  stop(
    "Spatial model variables are missing: ",
    paste(missing_model_variables, collapse = ", ")
  )
}

model_missingness <- vapply(
  DT_spatial[, required_model_variables, drop = FALSE],
  function(x) sum(is.na(x)),
  integer(1)
)

if (any(model_missingness > 0L)) {
  missingness_table <- tibble(
    variable = names(model_missingness),
    missing_n = as.integer(model_missingness)
  ) %>%
    filter(missing_n > 0L)

  write_csv(
    missingness_table,
    file.path(tables_dir, "Spatial_model_variable_missingness.csv")
  )

  stop("Missing data remain in spatial model variables. See the QC table.")
}

model_term_sets <- list(
  A_baseline = design_terms,
  B_environmental_HI = c(design_terms, environmental_hi_terms),
  B_environmental_WBGT = c(design_terms, environmental_wbgt_terms),
  C_social = c(design_terms, social_terms),
  D_full_HI = c(design_terms, environmental_hi_terms, social_terms),
  D_full_WBGT = c(design_terms, environmental_wbgt_terms, social_terms)
)

spatial_models <- imap(
  model_term_sets,
  ~ fit_spatial_model(
    model_name = .y,
    fixed_terms = .x,
    data = DT_spatial,
    A_observation = A_observation,
    spde_model = spde_model,
    force_refit = REFIT_MODELS
  )
)

saveRDS(
  spatial_models,
  file.path(models_dir, "Joint_SPDE_discrete_time_models_HI_WBGT_refit.rds"),
  compress = FALSE
)

saveRDS(
  mesh,
  file.path(models_dir, "Ghana_SPDE_mesh_refit.rds"),
  compress = FALSE
)

fixed_results <- imap_dfr(spatial_models, extract_fixed_hr)
fit_results <- imap_dfr(spatial_models, extract_fit_criteria)
hyperparameter_results <- imap_dfr(
  spatial_models,
  extract_hyperparameters
)

write_csv(
  fixed_results,
  file.path(tables_dir, "SPDE_fixed_effect_HRs_HI_WBGT_refit.csv")
)
write_csv(
  fit_results,
  file.path(tables_dir, "SPDE_model_comparison_HI_WBGT_refit.csv")
)
write_csv(
  hyperparameter_results,
  file.path(tables_dir, "SPDE_hyperparameters_HI_WBGT_refit.csv")
)


# ---- 7. Prediction grid and projected posterior means ----------------------

grid_polygons <- st_make_grid(
  ghana_boundary_utm,
  cellsize = MAP_CELL_SIZE,
  what = "polygons",
  square = TRUE
)

grid_sf <- st_sf(
  grid_id = seq_along(grid_polygons),
  geometry = grid_polygons
)

grid_sf <- suppressWarnings(
  st_intersection(grid_sf, ghana_boundary_utm)
) %>%
  st_collection_extract("POLYGON")

grid_sf <- grid_sf[!st_is_empty(grid_sf), ]
grid_prediction_points <- suppressWarnings(st_point_on_surface(grid_sf))
grid_coordinates <- st_coordinates(grid_prediction_points)

A_prediction <- INLA::inla.spde.make.A(
  mesh = mesh,
  loc = grid_coordinates
)

project_model_set <- function(models) {
  output <- grid_sf

  for (model_name in names(models)) {
    output[[model_name]] <- project_field_mean(
      models[[model_name]],
      A_prediction
    )
  }

  output
}

hi_model_names <- c(
  "A_baseline", "B_environmental_HI", "C_social", "D_full_HI"
)
wbgt_model_names <- c(
  "A_baseline", "B_environmental_WBGT", "C_social", "D_full_WBGT"
)

hi_models <- spatial_models[hi_model_names]
wbgt_models <- spatial_models[wbgt_model_names]

grid_means_hi <- project_model_set(hi_models)
grid_means_wbgt <- project_model_set(wbgt_models)

hi_limits <- range(
  st_drop_geometry(grid_means_hi)[hi_model_names],
  na.rm = TRUE
)
wbgt_limits <- range(
  st_drop_geometry(grid_means_wbgt)[wbgt_model_names],
  na.rm = TRUE
)


# ---- 8. Static maps with Ghana regional boundaries -------------------------

make_regional_map_panel <- function(
    grid_data,
    variable,
    panel_title,
    limits,
    show_region_labels = TRUE) {

  panel <- ggplot() +
    geom_sf(
      data = grid_data,
      aes(fill = .data[[variable]]),
      color = NA
    ) +
    geom_sf(
      data = ghana_regions_utm,
      fill = NA,
      color = "grey20",
      linewidth = 0.35
    ) +
    geom_sf(
      data = ghana_boundary_utm,
      fill = NA,
      color = "black",
      linewidth = 0.7
    ) +
    scale_fill_viridis_c(
      option = "C",
      limits = limits,
      oob = scales::squish,
      name = "Residual\nlog-hazard",
      na.value = "transparent"
    ) +
    coord_sf(datum = NA) +
    labs(title = panel_title) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  if (show_region_labels) {
    panel <- panel +
      geom_text(
        data = st_drop_geometry(region_labels),
        aes(x = label_x, y = label_y, label = region_label),
        size = 2.05,
        fontface = "bold",
        check_overlap = TRUE
      )
  }

  panel
}

panel_titles_hi <- c(
  A_baseline = "A. Baseline",
  B_environmental_HI = "B. Environmental: Heat Index",
  C_social = "C. Reproductive and socioeconomic",
  D_full_HI = "D. Fully adjusted: Heat Index"
)

panel_titles_wbgt <- c(
  A_baseline = "A. Baseline",
  B_environmental_WBGT = "B. Environmental: WBGT",
  C_social = "C. Reproductive and socioeconomic",
  D_full_WBGT = "D. Fully adjusted: WBGT"
)

hi_panels <- map(
  hi_model_names,
  ~ make_regional_map_panel(
    grid_means_hi,
    .x,
    panel_titles_hi[[.x]],
    hi_limits
  )
)

wbgt_panels <- map(
  wbgt_model_names,
  ~ make_regional_map_panel(
    grid_means_wbgt,
    .x,
    panel_titles_wbgt[[.x]],
    wbgt_limits
  )
)

figure_hi_regions <- wrap_plots(
  hi_panels,
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(
    title = "Residual spatial variation in miscarriage hazard",
    subtitle = "Heat Index spatial-model sequence"
  ) &
  theme(legend.position = "right")

figure_wbgt_regions <- wrap_plots(
  wbgt_panels,
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(
    title = "Residual spatial variation in miscarriage hazard",
    subtitle = "WBGT spatial-model sequence"
  ) &
  theme(legend.position = "right")

# WBGT is the selected main spatial specification.
save_high_resolution(
  figure_wbgt_regions,
  "Figure_6_WBGT_SPDE_Posterior_Mean_Fields_Regions",
  width = 11,
  height = 8.5
)

# HI is retained as the corresponding supplementary analysis.
save_high_resolution(
  figure_hi_regions,
  "Supplementary_Figure_S2_HI_SPDE_Posterior_Mean_Fields_Regions",
  width = 11,
  height = 8.5
)

# Direct comparison of the two fully adjusted thermal specifications. Use a
# common colour scale only for this comparison.
thermal_limits <- range(
  c(
    grid_means_hi$D_full_HI,
    grid_means_wbgt$D_full_WBGT
  ),
  na.rm = TRUE
)

figure_thermal_comparison <- wrap_plots(
  list(
    make_regional_map_panel(
      grid_means_hi,
      "D_full_HI",
      "A. Fully adjusted Heat Index model",
      thermal_limits
    ),
    make_regional_map_panel(
      grid_means_wbgt,
      "D_full_WBGT",
      "B. Fully adjusted WBGT model",
      thermal_limits
    )
  ),
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(
    title = "Fully adjusted residual spatial fields by thermal metric"
  ) &
  theme(legend.position = "right")

save_high_resolution(
  figure_thermal_comparison,
  "Supplementary_Figure_HI_WBGT_Fully_Adjusted_Comparison_Regions",
  width = 12,
  height = 6.5
)


# ---- 9. Posterior uncertainty maps for HI and WBGT -------------------------

make_uncertainty_figure <- function(
    fit,
    metric_label,
    output_stem) {

  uncertainty <- sample_projected_field(
    fit,
    A_prediction,
    N_POSTERIOR_SAMPLES
  )

  uncertainty_grid <- grid_sf
  for (variable in names(uncertainty)) {
    uncertainty_grid[[variable]] <- uncertainty[[variable]]
  }

  mean_limits <- range(
    uncertainty_grid$posterior_mean,
    na.rm = TRUE
  )
  sd_limits <- range(
    uncertainty_grid$posterior_sd,
    na.rm = TRUE
  )

  panels <- list(
    make_regional_map_panel(
      uncertainty_grid,
      "posterior_mean",
      "A. Posterior mean spatial effect",
      mean_limits,
      show_region_labels = FALSE
    ),
    make_regional_map_panel(
      uncertainty_grid,
      "posterior_sd",
      "B. Posterior uncertainty",
      sd_limits,
      show_region_labels = FALSE
    ),
    make_regional_map_panel(
      uncertainty_grid,
      "probability_positive",
      "C. Probability spatial effect > 0",
      c(0, 1),
      show_region_labels = FALSE
    )
  )

  figure <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title = paste0(
        "Residual spatial pattern and uncertainty: ",
        metric_label,
        " model"
      ),
      subtitle = paste0(
        "Joint posterior summaries based on ",
        N_POSTERIOR_SAMPLES,
        " latent-field samples"
      )
    )

  save_high_resolution(
    figure,
    output_stem,
    width = 14,
    height = 5.5
  )

  list(
    figure = figure,
    grid = uncertainty_grid
  )
}

wbgt_uncertainty_output <- make_uncertainty_figure(
  spatial_models$D_full_WBGT,
  "WBGT",
  "Figure_7_WBGT_SPDE_FullModel_Uncertainty_Regions"
)

hi_uncertainty_output <- make_uncertainty_figure(
  spatial_models$D_full_HI,
  "Heat Index",
  "Supplementary_Figure_S3_HI_SPDE_FullModel_Uncertainty_Regions"
)


# ---- 10. API-free Leaflet map function -------------------------------------

create_leaflet_spatial_map <- function(
    grid_data,
    model_names,
    model_labels,
    metric_label,
    limits,
    output_filename) {

  grid_ll <- grid_data %>%
    st_make_valid() %>%
    st_transform(4326)

  regions_ll <- ghana_regions_utm %>%
    st_make_valid() %>%
    st_transform(4326)

  boundary_ll <- ghana_boundary_utm %>%
    st_make_valid() %>%
    st_transform(4326)

  palette <- leaflet::colorNumeric(
    palette = viridisLite::plasma(256),
    domain = limits,
    na.color = "transparent"
  )

  map <- leaflet::leaflet(
    options = leaflet::leafletOptions(
      preferCanvas = TRUE,
      minZoom = 5
    )
  ) %>%
    # Explicit OpenStreetMap tiles: no API key is required. Do not replace
    # this with Carto provider tiles, which caused the earlier API-key message.
    leaflet::addTiles(
      urlTemplate = "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
      attribution = paste0(
        "&copy; <a href='https://www.openstreetmap.org/copyright'>",
        "OpenStreetMap contributors</a>"
      ),
      options = leaflet::tileOptions(
        maxZoom = 18,
        crossOrigin = TRUE
      ),
      group = "OpenStreetMap"
    )

  for (model_name in model_names) {
    layer_label <- unname(model_labels[[model_name]])
    layer_data <- grid_ll
    layer_data$map_value <- layer_data[[model_name]]
    layer_data$popup_text <- paste0(
      "<strong>", layer_label, "</strong><br>",
      "Residual log-hazard: ",
      formatC(layer_data$map_value, digits = 3, format = "f")
    )

    map <- map %>%
      leaflet::addPolygons(
        data = layer_data,
        group = layer_label,
        stroke = FALSE,
        smoothFactor = 0,
        fillColor = ~ palette(map_value),
        fillOpacity = 0.72,
        popup = ~ popup_text,
        highlightOptions = leaflet::highlightOptions(
          weight = 1.5,
          color = "#FFFFFF",
          fillOpacity = 0.88,
          bringToFront = TRUE
        )
      )
  }

  map <- map %>%
    leaflet::addPolygons(
      data = regions_ll,
      group = "Ghana regional boundaries",
      fill = FALSE,
      color = "#303030",
      opacity = 0.95,
      weight = 1.2,
      label = ~ region_label
    ) %>%
    leaflet::addPolygons(
      data = boundary_ll,
      group = "Ghana national boundary",
      fill = FALSE,
      color = "#000000",
      opacity = 1,
      weight = 2
    ) %>%
    leaflet::addLegend(
      position = "bottomright",
      pal = palette,
      values = limits,
      title = "Residual log-hazard",
      opacity = 0.85,
      labFormat = leaflet::labelFormat(digits = 2)
    ) %>%
    leaflet::addLayersControl(
      baseGroups = "OpenStreetMap",
      overlayGroups = c(
        unname(model_labels),
        "Ghana regional boundaries",
        "Ghana national boundary"
      ),
      options = leaflet::layersControlOptions(collapsed = FALSE)
    ) %>%
    leaflet::hideGroup(
      unname(model_labels[model_names[-length(model_names)]])
    ) %>%
    leaflet::addControl(
      html = paste0(
        "<div style='background:rgba(255,255,255,0.94);padding:9px 12px;",
        "border-radius:4px;box-shadow:0 1px 5px rgba(0,0,0,0.35);",
        "font-family:Arial,sans-serif;'>",
        "<strong>Residual spatial variation in miscarriage hazard</strong>",
        "<br><span style='font-size:12px;'>",
        metric_label,
        " spatial model</span></div>"
      ),
      position = "topleft"
    )

  bbox <- st_bbox(boundary_ll)

  map <- map %>%
    leaflet::fitBounds(
      lng1 = unname(bbox[["xmin"]]),
      lat1 = unname(bbox[["ymin"]]),
      lng2 = unname(bbox[["xmax"]]),
      lat2 = unname(bbox[["ymax"]])
    )

  output_path <- file.path(spatial_dir, output_filename)

  htmlwidgets::saveWidget(
    widget = map,
    file = output_path,
    selfcontained = TRUE,
    title = paste0(
      "Residual spatial variation in miscarriage hazard: ",
      metric_label
    )
  )

  message("Saved ", metric_label, " interactive map: ", output_path)
  map
}

hi_leaflet_labels <- c(
  A_baseline = "A. Baseline",
  B_environmental_HI = "B. Environmental: Heat Index",
  C_social = "C. Reproductive and socioeconomic",
  D_full_HI = "D. Fully adjusted: Heat Index"
)

wbgt_leaflet_labels <- c(
  A_baseline = "A. Baseline",
  B_environmental_WBGT = "B. Environmental: WBGT",
  C_social = "C. Reproductive and socioeconomic",
  D_full_WBGT = "D. Fully adjusted: WBGT"
)

hi_leaflet <- create_leaflet_spatial_map(
  grid_data = grid_means_hi,
  model_names = hi_model_names,
  model_labels = hi_leaflet_labels,
  metric_label = "Heat Index",
  limits = hi_limits,
  output_filename = "Interactive_HI_SPDE_Posterior_Mean_Fields.html"
)

wbgt_leaflet <- create_leaflet_spatial_map(
  grid_data = grid_means_wbgt,
  model_names = wbgt_model_names,
  model_labels = wbgt_leaflet_labels,
  metric_label = "WBGT",
  limits = wbgt_limits,
  output_filename = "Interactive_WBGT_SPDE_Posterior_Mean_Fields.html"
)


# ---- 11. Save map-ready spatial outputs ------------------------------------

spatial_gpkg <- file.path(
  spatial_dir,
  "SPDE_spatial_field_predictions_HI_WBGT_refit.gpkg"
)

st_write(
  grid_means_hi,
  spatial_gpkg,
  layer = "posterior_mean_HI_models",
  delete_layer = TRUE,
  quiet = TRUE
)

st_write(
  grid_means_wbgt,
  spatial_gpkg,
  layer = "posterior_mean_WBGT_models",
  delete_layer = TRUE,
  quiet = TRUE
)

st_write(
  hi_uncertainty_output$grid,
  spatial_gpkg,
  layer = "full_HI_uncertainty",
  delete_layer = TRUE,
  quiet = TRUE
)

st_write(
  wbgt_uncertainty_output$grid,
  spatial_gpkg,
  layer = "full_WBGT_uncertainty",
  delete_layer = TRUE,
  quiet = TRUE
)


# ---- 12. Reproducibility record --------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "SessionInfo_03_Spatial_Refit.txt")
)

message(
  "Script 03 completed successfully.",
  "\nStatic figures: ", figures_dir,
  "\nInteractive HTML maps and map-ready outputs: ", spatial_dir,
  "\nModel checkpoints: ", checkpoint_dir
)

