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

**Countries:** Gambia, Ghana, Sierra Leone, and Malawi (all active).

**Outcomes:** Vitamin A deficiency (children & women), iron deficiency anaemia (children & women); several countries also have folate, vitamin B12, and zinc.

## Quick Start

### Prerequisites

```r
# Core packages
install.packages(c("targets", "here", "dplyr", "tidyr", "readr", "ggplot2",
                    "sf", "geodata", "terra", "exactextractr", "glmnet",
                    "srvyr", "survey", "pROC", "ggrepel", "scales", "viridis",
                    "patchwork", "haven", "labelled", "caret", "recipes",
                    "future.apply", "data.table"))

# SuperLearner stack (mlr3) — the active modelling backend
install.packages(c("mlr3", "mlr3learners", "mlr3extralearners",
                    "ranger", "xgboost", "dbarts", "kernlab", "origami"))
# mlr3superlearner is GitHub-only:
# remotes::install_github("nt-williams/mlr3superlearner")

# Preprocessing / screening helpers
install.packages("ck37r", repos = "https://cloud.r-project.org")
install.packages("washb", repos = "https://cloud.r-project.org")

# Note: the legacy sl3 path (R/sl_fitting.R, R/bootstrap.R) is retained for
# reference only and is NOT part of the active `targets` graph.
```

### Run the pipeline

```r
# From the project root:
targets::tar_make()
```

That's it. The pipeline defaults to **fast mode**, which runs in ~10-30 minutes and is sufficient for development and debugging.

## Pipeline Modes

The pipeline has two modes that control the SuperLearner stack complexity:

### Fast mode (default)

For development, debugging, and checking relative performance across outcomes. Results are directionally correct but confidence intervals are rough.

```r
# Fast mode is the default — just run:
targets::tar_make()
```

| Parameter | Value |
|-----------|-------|
| SL stack | 5 learners (mean + lasso + elastic net + ranger + xgboost) |
| CV folds | 5 |
| Uncertainty | Split + locally-adaptive conformal intervals |
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
| SL stack | ~16 learners (glmnet ×3, ranger ×3, xgboost ×2, BART ×3, Gaussian process, + class-weighted variants for rare outcomes) |
| CV folds | 5 |
| Uncertainty | Split + locally-adaptive conformal intervals |
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
| `cv_perf_{country}_{outcome}` | Extract out-of-fold (cross-validated) performance metrics | No |
| `admin1_prev_{country}_{outcome}` | Aggregate predictions to Admin-1 prevalence | No |
| `conformal_ci_{country}_{outcome}` | Conformal prediction intervals (Admin-1/national) | No |
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
│   ├── config.R              #   Country configs, outcome definitions, mode params
│   ├── data_prep.R           #   Load data, build outcome-specific datasets
│   ├── mlr3_fitting.R        #   SuperLearner fitting (active, mlr3) + out-of-fold preds
│   ├── sl_fitting.R          #   Legacy sl3 fitting (reference only)
│   ├── conformal.R           #   Conformal prediction intervals (active uncertainty)
│   ├── bootstrap.R           #   Legacy cluster bootstrap (reference only)
│   ├── admin1_analysis.R     #   Admin-1 aggregation and CV performance
│   ├── admin2_analysis.R     #   Admin-2 analysis + area-level GEE model
│   ├── area_level_comparison.R #  Area-level SL + cross-country LOCO (GEE-only)
│   ├── transportability_area.R #  Universal area-level transportability model
│   ├── domain_ablation.R     #   Leave-one-domain-out feature importance
│   └── plotting.R            #   Maps, scatter plots, forest plots
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
- **SuperLearner:** `mlr3superlearner` (`R/mlr3_fitting.R`) with cluster-blocked cross-validation (cluster IDs passed via `group=`, so a PSU's observations stay together across folds). Discrete SuperLearner — CV selects the single best learner. The legacy `sl3` path (`R/sl_fitting.R`) is retained for reference only.
- **Out-of-fold predictions:** Reported CV performance and the residuals feeding the conformal intervals use genuine out-of-fold predictions (`res$yhat_full`), recomputed by resampling the fitted learners (`mlr3_oof_predictions()`), since `mlr3superlearner` does not retain its internal CV predictions. `res$yhat_insample` keeps the resubstitution predictions for reference and for permutation-importance baselines.
- **Uncertainty:** Conformal prediction intervals (`R/conformal.R`) — split (constant-width) and locally-adaptive variants — built from out-of-fold residuals. (The earlier cluster-bootstrap path in `R/bootstrap.R` is retained for reference only.)
- **Individual-level area aggregation:** Out-of-fold individual predictions are survey-weighted and aggregated to Admin-1/Admin-2 prevalence.
- **Area-level model:** Elastic net (`glmnet`, alpha = 0.5) trained on survey-weighted Admin-2 prevalence ~ GEE raster zonal means. Enables prediction to Admin-2 areas with no survey data.
- **Cross-country transportability (`R/transportability_area.R`):** A universal, parsimonious within-country-centered elastic net on harmonized GEE + IHME + Malaria-Atlas + food-security proxies, validated by leave-one-country-out CV (`area_transport_*` targets).
