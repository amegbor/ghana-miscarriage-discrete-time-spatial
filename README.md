# Environmental and Socioeconomic Determinants of Miscarriage in Ghana

R code for the analysis of monthly heat, air pollution and vegetation exposures and miscarriage, using pooled pregnancy histories from the 2017 Ghana Maternal Health Survey and the 2022 Ghana Demographic and Health Survey. The analysis uses weighted multilevel discrete-time survival models, natural-spline exposure-response models, sensitivity analyses and Bayesian INLA-SPDE spatial models.

## Data

This repository contains code only. No survey data, GPS files, environmental data or model outputs are included.

- Survey data and GPS coordinates are available from the [DHS Program](https://dhsprogram.com) and are subject to its terms of use.
- Environmental data are available from their original providers, as cited in the paper.

The scripts start from `Ghana_Discrete_Time_Monthly.rds`, a monthly person-pregnancy file produced by an earlier data-preparation step that is not included here. Place it at:

```text
data/Ghana_Discrete_Time_Monthly.rds
```

## Setup

By default, the scripts read from `data/` and write results to `Results/` in the working directory. To use other locations, set these before running:

```r
Sys.setenv(
  GHANA_MISCARRIAGE_ROOT = "/path/to/project",
  GHANA_DATA_DIR = "/path/to/data",
  GHANA_N_THREADS = "4"
)
```

## Run order

```bash
Rscript 02_Ghana_DiscreteTime_Analysis.R
Rscript 04_Ghana_Sensitivity_Analyses_2_.R
Rscript 04B_WBGT_Sensitivity_Analyses.R
Rscript 03_Ghana_Spatial_Analysis_Revised.R
Rscript 05_Ghana_Primary_Model_Diagnostics_Extraction.R
```

| Script | What it does |
| --- | --- |
| `02` | Builds the analytic sample; fits the linear, spline and multi-exposure models; produces descriptive tables and figures |
| `04` | Heat Index sensitivity analyses, effect modification, EVI by residence and reporting-bias scenarios |
| `04B` | WBGT sensitivity analyses, GEE and pregnancy-level models |
| `03` | Spatial models and maps |
| `05` | Model diagnostics |

Each script needs the outputs of the scripts before it.

## Requirements

R packages:

```text
dplyr, tidyr, purrr, tibble, readr, haven, stringr, splines,
glmmTMB, geepack, emmeans, MASS, survey, openxlsx,
ggplot2, patchwork, sf, INLA, fmesher, Matrix, geodata,
leaflet, htmlwidgets, htmltools, viridisLite, scales
```

Install INLA from [r-inla.org](https://www.r-inla.org/download-install). Script 03 downloads Ghana's administrative boundaries on first run, so it needs internet access.



## Code was tidied for public release with AI assistance and checked by the authors

## License

MIT. The license covers the code only, not the survey or environmental data.
