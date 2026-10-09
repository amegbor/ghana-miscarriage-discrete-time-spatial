# Environmental and Socioeconomic Determinants of Miscarriage in Ghana

## Discrete-Time Survival and Spatial Analysis

Analysis code for pooled pregnancy histories from the 2017 Ghana Maternal Health Survey and the 2022 Ghana Demographic and Health Survey. The workflow fits weighted multilevel discrete-time survival models, nonlinear exposure-response models, sensitivity analyses, and Bayesian INLA-SPDE spatial models.

This repository contains code only. Survey microdata, GPS files, derived pregnancy records, environmental extracts, fitted models, and generated outputs are not distributed. Obtain survey data and GPS access through the DHS Program and follow the applicable data-use agreements. Environmental inputs must be obtained from their original providers.

## Scope and status

These are portable copies of the corrected analysis scripts. Cooking fuel is excluded from adjustment and complete-case eligibility. Models retain mother- and cluster-level random intercepts; region random intercepts are excluded. HI and WBGT are evaluated separately. Linear heat effects use original units (per 1 degree Celsius), not IQR scaling; linear EVI effects use a 0.1-unit increment.

The cooking-fuel correction requires rebuilding the primary sample and refitting downstream analyses. Previously reported sample sizes and estimates must be updated after the corrected runs. These public copies have received static checks but have not been executed end to end. Package versions and session information from the actual analysis environment should be recorded for a reproducibility release; no dependency lockfile is supplied here.

The obsolete `02B_Ghana_Refit_Noncompliant_Linear_Models.R` is intentionally omitted. Corrected Script 02 fits the primary linear models directly.

## Private input

The starting input is `Ghana_Discrete_Time_Monthly.rds`, produced by the upstream pregnancy reconstruction and environmental linkage workflow. It is not included, and these analysis scripts do not reconstruct it from raw DHS files. The previously saved `Ghana_DT_AnalysisReady_PrimaryCC.rds` alone is insufficient to rebuild pregnancies excluded by the old cooking-fuel filter.

By default, run from the repository root and place the starting input at:

```text
data/Ghana_Discrete_Time_Monthly.rds
```

The input must contain the harmonized monthly risk set, pregnancy and mother identifiers, cluster identifiers, outcome, weights, exposure series, timing variables, and adjustment variables expected by Script 02. Spatial analysis additionally requires linked cluster coordinates. Each script performs its own input checks.

Script 02 writes the rebuilt analysis-ready datasets and fitted models under `Results/Models`. Later scripts require these outputs and reject datasets lacking the corrected-sample marker. Do not substitute older files named `Ghana_DT_AnalysisReady_PrimaryCC_CookingCorrected.rds`.

## Paths

No user-specific HPC paths are embedded in these scripts. All scripts use the same project-root setting:

| Variable | Default | Purpose |
| --- | --- | --- |
| `GHANA_MISCARRIAGE_ROOT` | Current working directory | Project root; outputs go to `Results` beneath it |
| `GHANA_DATA_DIR` | `<project root>/data` | Directory containing the starting monthly RDS input; used by Script 02 |
| `GHANA_N_THREADS` | Script-specific default | Requested model thread count; set to match available resources |

The project root must already exist. For a private working directory outside the repository, set both paths explicitly. The following R settings can be used before sourcing scripts in RStudio:

```r
Sys.setenv(
  GHANA_MISCARRIAGE_ROOT = "/path/to/private/project",
  GHANA_DATA_DIR = "/path/to/private/inputs",
  GHANA_N_THREADS = "4"
)
```

On a Unix shell, use `export` before the commands below. The old `GHANA_MISCARRIAGE_DIR` setting is not used by these copies.

## Run order

Run sequentially, allowing each script to finish successfully before continuing:

```bash
Rscript 02_Ghana_DiscreteTime_Analysis.R
Rscript 04_Ghana_Sensitivity_Analyses_2_.R
Rscript 04B_WBGT_Sensitivity_Analyses.R
Rscript 03_Ghana_Spatial_Analysis_Revised.R
Rscript 05_Ghana_Primary_Model_Diagnostics_Extraction.R
```

| Script | Purpose |
| --- | --- |
| `02_Ghana_DiscreteTime_Analysis.R` | Rebuild the common analytic sample; fit linear, spline, and multi-exposure primary models; produce descriptive tables and figures |
| `04_Ghana_Sensitivity_Analyses_2_.R` | Main HI sensitivities, HI/WBGT effect modification and category-specific HRs, EVI-by-residence analysis, and HI/WBGT assumption-based reporting-bias scenarios |
| `04B_WBGT_Sensitivity_Analyses.R` | WBGT sensitivity, population-averaged GEE, and pregnancy-level binary models |
| `03_Ghana_Spatial_Analysis_Revised.R` | INLA-SPDE spatial models, residual spatial fields, uncertainty, and map exports |
| `05_Ghana_Primary_Model_Diagnostics_Extraction.R` | Extract primary diagnostics from the current corrected checkpoints and audit linear and nonlinear formulas |

With the risk-set indexing used by the upstream preparation workflow, `gest_month = 1` represents zero completed months. Both sensitivity scripts exclude entire pregnancies with a recorded miscarriage at `gest_month - 1 < 2`. Binary-model exposures average available observations within risk intervals 1-3; losses before interval 3 contribute only their available intervals. Binary analyses include live birth and reported miscarriage.

## Dependencies

The scripts check required packages before running. The combined dependencies are:

```text
dplyr, tidyr, purrr, tibble, readr, haven, stringr, splines,
glmmTMB, geepack, emmeans, MASS, survey, openxlsx,
ggplot2, patchwork, sf, INLA, fmesher, Matrix, geodata,
leaflet, htmlwidgets, htmltools, viridisLite, scales
```

Use the official INLA installation instructions for INLA. Spatial packages may require system libraries. Script 03 downloads Ghana administrative boundaries when its local cache is absent; interactive basemaps require internet access. The boundary cache and generated spatial outputs are ignored by Git.

## Checkpoints and outputs

Corrected runs use new checkpoint namespaces, so older v6/v7 fits are not selected. Script 05 only inventories `glmmTMB_checkpoints_no_cooking_v1` and recognizes `EVI_per_0_1` filenames. Do not manually repopulate the corrected directories with legacy models.

To force refitting after another data or model change, set `GHANA_REFIT_MODELS=true` for Script 02 and `GHANA_REFIT_SENSITIVITY=true` for Script 04. Script 03 refits by default; keep `GHANA_REFIT_SPATIAL_MODELS=true`. Script 04B resumes from its corrected checkpoint directory; use a fresh checkpoint-directory name in that script after further data or specification changes.

Generated results retain the established `Results/Tables`, `Results/Figures`, `Results/Models`, `Results/Spatial`, and `Results/Diagnostics` locations. Keep outputs private by default; this repository does not automatically publish tables or maps.

## Before uploading

The supplied `.gitignore` excludes data directories, RDS files (including `extract_WBGT.rds` and cluster-month files), survey files, shapefile components, spatial databases, rasters, archives, tabular extracts, checkpoints, and outputs. Review `git status --short` before committing. Ignore rules do not remove files already tracked by Git or erase earlier commits. Do not add restricted inputs with `git add -f`.

No software license is assigned by this package. Add a license selected by the repository owner before treating the code as an openly licensed release.
