# Proxy Models for Micronutrient Deficiency — Tutorial

A simulation-based, fully self-contained tutorial teaching biostatistical
researchers how to build proxy models for micronutrient deficiency using
cluster-sampled survey microdata, geospatial covariates, and Super Learner
ensembles.

## What this tutorial covers

| Section | Topic |
|---------|-------|
| A | Motivation — why proxy models, three prediction targets |
| B | Simulated data generating process (spatial field, cluster survey) |
| C | Preprocessing: MCAR/MAR missingness, imputation, standardisation |
| D | Super Learner — theory, cluster-blocked CV, continuous + binary models |
| E | Four Admin1 prevalence estimation approaches + comparison |
| F | Uncertainty: cluster bootstrap, delta method, GAM-Monte Carlo |
| G | Spatial autocorrelation, wall-to-wall grid maps, Moran's I |
| H | Loss functions and rare outcomes (PR-AUC, focal loss, Brier) |
| I | Post-hoc threshold tuning with nested CV |
| J | Variable importance and grouped domain contribution |
| K | Alternative models: spatial GAM, Fay-Herriot, geostatistical |
| L | Practical checklist and common failure modes |

## File structure

```
mn_proxy_tutorial/
├── mn_proxy_tutorial.qmd          # Main tutorial (render this)
├── README.md                      # This file
└── R/
    ├── utils.R                    # Metrics, calibration, loss functions
    ├── simulate_data.R            # DGP: survey + covariate grid
    ├── preprocess.R               # Missingness, imputation, scaling
    ├── fit_superlearner.R         # Cluster CV, SL wrappers, SL fitting
    ├── mapping_and_aggregation.R  # 4 prevalence approaches, grid predict
    ├── uncertainty.R              # Bootstrap, delta method, GAM-MC
    └── spatial_methods.R         # Moran's I, spatial CV, Fay-Herriot
```

## Requirements

### R version
R ≥ 4.2.0 recommended.

### Required CRAN packages

Install all at once:

```r
install.packages(c(
  # Core modeling
  "SuperLearner",  # Super Learner ensemble
  "glmnet",        # Elastic net base learner
  "ranger",        # Random forest base learner
  "mgcv",          # GAM with spatial smooths
  "MASS",          # mvrnorm for spatial simulation

  # Spatial
  "sf",            # Spatial features (polygons, points)
  "spdep",         # Moran's I, spatial weights

  # Survey
  "srvyr",         # Tidy survey design (optional — used in 06_ script)
  "survey",        # Survey design and estimation
  "sae",           # Fay-Herriot small area estimation

  # Performance metrics
  "pROC",          # ROC / AUC
  "PRROC",         # Precision-Recall AUC

  # Data manipulation
  "dplyr",
  "tidyr",
  "purrr",

  # Visualisation
  "ggplot2",
  "patchwork",
  "scales",
  "viridis",

  # Output
  "knitr",
  "kableExtra",

  # Optional — improves GLMM in section K
  "lme4"
))
```

### Optional: xgboost

For gradient boosting as an additional base learner:

```r
install.packages("xgboost")
# Then add "SL.xgboost" to SL_LIBRARY in the QMD setup chunk
```

### Quarto

Install Quarto from <https://quarto.org/docs/get-started/>.

Verify installation:

```bash
quarto --version   # should print 1.3 or higher
```

## How to render the tutorial

### From the terminal

```bash
cd mn_proxy_tutorial
quarto render mn_proxy_tutorial.qmd
```

The output file `mn_proxy_tutorial.html` will appear in the same directory.
Open it in any browser.

### From R / RStudio

```r
setwd("mn_proxy_tutorial")   # working directory MUST be here
quarto::quarto_render("mn_proxy_tutorial.qmd")
```

Or use the **Render** button in RStudio with the .qmd open.

### Expected runtime

| Section | Approximate time |
|---------|-----------------|
| Data simulation | < 1 s |
| Preprocessing | < 1 s |
| SL continuous model (n=1440, 5 learners, 5 folds) | 3–8 min |
| SL binary model | 3–8 min |
| Fast bootstrap (B=100) | 30–60 s |
| GAM uncertainty (n_sim=300) | 10–30 s |
| Grid prediction | 30–60 s |
| Everything else | < 1 s each |
| **Total** | **~10–20 min** |

Caching is enabled (`cache: true` in the YAML header). After the first render,
subsequent renders complete in < 1 minute.  Delete the `_cache/` directory to
force a full rerun.

## Reproducibility

All stochastic operations use fixed seeds:

| Constant | Value |
|----------|-------|
| `SEED` | 42 |
| Simulation seed | 42 |
| Grid simulation seed | 99 |
| Bootstrap seed | 42 |
| Fold seed | 42 |

Set `SEED <- 42` at the top of the setup chunk to reproduce all results exactly.

## renv (optional)

To snapshot the exact package versions used:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

This creates `renv.lock`. To restore on another machine:

```r
renv::restore()
```

## Adapting to real data

To swap in real micronutrient survey data:

1. Replace the `simulate_mn_survey()` call in Section B with your own
   `readRDS()` / `read.csv()` call.
2. Ensure your data has columns matching those described in Section B
   (or update `get_predictor_cols()` / `make_domain_groups()` in `preprocess.R`).
3. Replace `simulate_covariate_grid()` in Section G with real raster extracts
   (use `terra::extract()` or `exactextractr::exact_extract()`).
4. Replace `simulate_admin1_polygons()` with real boundary shapefiles via
   `sf::st_read()` or `geodata::gadm()`.

## Key design decisions

- **CRAN-only packages** — no GitHub-only dependencies required.
- **SuperLearner** (CRAN) is the primary ensemble framework.  Section D
  notes that `sl3` (GitHub/tlverse) offers a more modern API if available.
- **Cluster-blocked CV** is the default throughout — individual-level
  randomisation is explicitly warned against.
- **No external data** — everything is simulated and deterministic.
