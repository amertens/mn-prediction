# Micronutrient Deficiency Prediction Pipeline

Sub-national prediction of micronutrient deficiency prevalence using survey biomarker data and remotely sensed environmental covariates. The pipeline combines individual-level SuperLearner models with area-level ecological models to generate Admin-2 prevalence estimates, including for areas with no survey coverage.

## Overview

```
Survey biomarker data  ──┐
                         ├──> SuperLearner ──> Admin-1/Admin-2 prevalence
GEE raster covariates  ──┘                         │
                                                   ▼
                                    Area-level model (glmnet)
                                            │
                                            ▼
                              Predictions for unsurveyed areas
```

**Countries:** Gambia (active), with Ghana, Sierra Leone, and Malawi in development.

**Outcomes:** Vitamin A deficiency (children & women), iron deficiency anaemia (children & women).

## Quick Start

### Prerequisites

```r
# Core packages
install.packages(c("targets", "here", "dplyr", "tidyr", "readr", "ggplot2",
                    "sf", "geodata", "terra", "exactextractr", "glmnet",
                    "srvyr", "survey", "pROC", "ggrepel", "scales", "viridis",
                    "patchwork", "haven", "labelled", "caret", "recipes",
                    "future.apply", "data.table"))

# tlverse (SuperLearner)
# See https://tlverse.org for installation instructions
install.packages(c("sl3", "origami"), repos = "https://tlverse.r-universe.dev")

# Additional
install.packages(c("ck37r", "washb"), repos = "https://cloud.r-project.org")
```

### Run the pipeline

```r
# From the project root:
targets::tar_make()
```

That's it. The pipeline defaults to **fast mode**, which runs in ~10-30 minutes and is sufficient for development and debugging.

## Pipeline Modes

The pipeline has two modes that control the SuperLearner stack complexity and bootstrap replication count:

### Fast mode (default)

For development, debugging, and checking relative performance across outcomes. Results are directionally correct but confidence intervals are rough.

```r
# Fast mode is the default — just run:
targets::tar_make()
```

| Parameter | Value |
|-----------|-------|
| SL stack | 3 learners (mean + glmnet + ranger) |
| CV folds | 5 |
| Admin-1 bootstrap | 10 reps |
| Admin-2 SL bootstrap | 5 reps |
| Area-level bootstrap | 20 reps |
| Approximate runtime | 10-30 minutes |

### Full mode

For publication-quality results with stable confidence intervals and the complete learner library.

```r
# Set the environment variable before running:
Sys.setenv(PIPELINE_MODE = "full")
targets::tar_make()
```

Or from the command line:

```bash
PIPELINE_MODE=full Rscript -e 'targets::tar_make()'
```

| Parameter | Value |
|-----------|-------|
| SL stack | 5+ learners with screener pipelines (RF, lasso, correlation) |
| CV folds | 5 |
| Admin-1 bootstrap | 200 reps |
| Admin-2 SL bootstrap | 50 reps |
| Area-level bootstrap | 500 reps |
| Approximate runtime | 2-4 hours |

### Switching between modes

When you switch modes, you need to invalidate the cached parameters so downstream targets re-execute with the new settings:

```r
Sys.setenv(PIPELINE_MODE = "full")
targets::tar_invalidate(names = "pipeline_params")
targets::tar_make()
```

## Pipeline Structure

The pipeline is managed by [{targets}](https://docs.ropensci.org/targets/), which tracks dependencies between steps and only re-runs what has changed. Visualize the dependency graph with:

```r
targets::tar_visnetwork()
```

### Targets (per country x outcome)

| Target | Description | Expensive? |
|--------|-------------|-----------|
| `merged_{country}` | Load merged survey + covariate dataset | No |
| `outcome_data_{country}_{outcome}` | Filter to population, select predictors, remove leakage | No |
| `sl_fit_{country}_{outcome}` | Fit continuous + binary SuperLearner models | **Yes** |
| `cv_perf_{country}_{outcome}` | Extract cross-validated performance metrics | No |
| `admin1_prev_{country}_{outcome}` | Aggregate predictions to Admin-1 prevalence | No |
| `bootstrap_ci_{country}_{outcome}` | Cluster bootstrap for Admin-1/national CIs | **Yes** |
| `ablation_{country}_{outcome}` | Domain ablation (leave-one-domain-out importance) | **Yes** |
| `svy_admin2_{country}_{outcome}` | Survey-weighted Admin-2 prevalence (design-based) | No |
| `admin2_sl_{country}_{outcome}` | Aggregate individual SL predictions to Admin-2 | No |
| `admin2_error_{country}_{outcome}` | Error metrics: SL vs survey at Admin-2 | No |
| `area_model_{country}_{outcome}` | Area-level elastic net: GEE rasters to Admin-2 prevalence | No |
| `plot_*_{country}_{outcome}` | Maps, scatter plots, forest plots | No |

### Summary targets

| Target | Description |
|--------|-------------|
| `cv_perf` | Combined CV performance table across all outcomes |
| `admin2_error_all` | Combined Admin-2 error table |
| `ablation_all` | Combined domain ablation results |
| `save_tables` | Write combined CSV files to `results/tables/` |

## Project Layout

```
mn-prediction/
├── _targets.R              # Pipeline definition (start here)
├── R/                      # Pipeline function modules
│   ├── config.R            #   Country configs, outcome definitions, mode params
│   ├── data_prep.R         #   Load data, build outcome-specific datasets
│   ├── sl_fitting.R        #   SuperLearner fitting (fast/full stacks)
│   ├── admin1_analysis.R   #   Admin-1 aggregation and CV performance
│   ├── bootstrap.R         #   Cluster bootstrap uncertainty
│   ├── domain_ablation.R   #   Leave-one-domain-out feature importance
│   ├── admin2_analysis.R   #   Admin-2 analysis + area-level GEE model
│   └── plotting.R          #   Maps, scatter plots, forest plots
├── src/                    # Legacy source code and country-specific scripts
│   ├── analysis/           #   Core analysis scripts (01-05) + sl_helpers.R
│   ├── Gambia/             #   Gambia data cleaning, merging, DHS aggregation
│   ├── Ghana/              #   Ghana analysis scripts
│   └── ...
├── scripts/                # Standalone scripts and runners
│   ├── 06_admin2_predictions_map.R   # Admin-2 pipeline (standalone version)
│   └── run_full_pipeline.R           # Legacy sequential runner
├── data/                   # Input data (gitignored)
│   ├── IPD/                #   Individual participant data by country
│   └── {Country}_GEE_rasters/  #   GEE-derived raster covariates (.tif)
├── results/                # Pipeline outputs
│   ├── models/             #   Saved SL model objects (.rds)
│   ├── tables/             #   CSV tables (CV performance, error, ablation)
│   ├── figures/            #   PNG figures (maps, scatter plots, forest plots)
│   └── admin2/             #   Admin-2 specific outputs
├── docs/                   # Reports and presentations
│   ├── stakeholder_walkthrough_proxy_vad.qmd  # Stakeholder-facing walkthrough
│   └── helpers/            #   Data preparation for QMD documents
└── _targets/               # Targets cache (gitignored)
```

## Useful `targets` Commands

```r
# Run the pipeline (only outdated targets):
targets::tar_make()

# See what's outdated:
targets::tar_outdated()

# Visualize dependency graph:
targets::tar_visnetwork()

# Read a cached result:
targets::tar_read(cv_perf)
targets::tar_read(area_model_gambia_child_vitA)

# Force a specific target to re-run:
targets::tar_invalidate(names = "sl_fit_gambia_child_vitA")
targets::tar_make()

# See errors from the last run:
targets::tar_meta(fields = any_of("error"), complete_only = TRUE)

# Clean the entire cache (re-runs everything):
targets::tar_destroy()
```

## Adding a New Country

1. Add a new entry to `get_country_configs()` in `R/config.R`:

```r
Ghana = list(
  country     = "Ghana",
  gadm_code   = "GHA",
  data_path   = here::here("data", "IPD", "Ghana", "Ghana_merged_dataset.rds"),
  raster_dir  = here::here("data", "Ghana_GEE_rasters"),
  cluster_id  = "gw_cnum",
  admin1_col  = "Admin1",
  admin2_col  = "Admin2",
  # ... (see Gambia entry for full template)
  outcomes = list(
    child_vitA = list(tag = "child_vitA", ...)
  )
)
```

2. Run `targets::tar_make()` -- new targets are auto-generated.

## Adding a New Outcome

1. Add to the `outcomes` list in the relevant country config in `R/config.R`.
2. Run `targets::tar_make()`.

## Key Methodological Notes

- **Ecological prediction model:** Area-level predictors map to area-level prevalence. This is explicitly *not* individual risk prediction or causal inference.
- **Survey design:** Stratified cluster sampling with `srvyr::as_survey_design()`. PSU = cluster number, weights = survey weights.
- **SuperLearner:** Cluster-blocked cross-validation via `origami::make_folds(cluster_ids = ...)`. NNLS metalearner.
- **Area-level model:** Elastic net (`glmnet`, alpha = 0.5) trained on survey-weighted Admin-2 prevalence ~ GEE raster zonal means. Enables prediction to Admin-2 areas with no survey data.
- **Bootstrap:** Cluster-level resampling (resample PSUs with replacement, refit SL, predict on original data).
