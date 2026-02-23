# Spatial ensemble learning for prevalence prediction and mapping (simulated example)

This folder contains a fully reproducible, end-to-end R workflow illustrating:

1. Simulation of cluster-level survey prevalence data with residual spatial autocorrelation
2. Spatial block cross-validation
3. Four modeling strategies for incorporating spatial information into an ensemble workflow
4. Cross-validated predictive accuracy metrics and calibration diagnostics
5. Wall-to-wall mapping with uncertainty using a spatial second-stage model

## Files

- `00_setup.R`  
  Installs (if needed) and loads packages, sets a reproducible seed, and defines paths.

- `01_simulate_data.R`  
  Simulates:
  - cluster survey data: `(Y, N)`, lon/lat, non-spatial covariates  
  - wall-to-wall raster covariates over the domain  
  - a latent spatial Gaussian process surface to induce residual spatial correlation

- `02_fit_models_cv.R`  
  Fits and compares:
  - (A) Non-spatial Super Learner (stacking of non-spatial learners)
  - (B) Super Learner with coordinates included in base learners
  - (C) Two-stage: non-spatial Super Learner + spatial meta-learner (GAM "gp" smooth)
  - (D) Single spatial model with linear mean + spatial random field (GAM "gp" smooth)

  Uses spatial block cross-validation to generate out-of-fold predictions and compute:
  - Brier score (MSE on probability scale)
  - Cross-validated R² (weighted)
  - Calibration intercept/slope and an expected calibration error (ECE)

- `03_mapping_uncertainty.R`  
  Demonstrates mapping and uncertainty for Method (C):
  - refit stage-1 ensemble on full data
  - refit stage-2 spatial meta-learner on full data
  - predict on a grid with wall-to-wall raster covariates
  - produce a posterior-approximate mean and standard deviation map via simulation

## How to run

From an R session with working directory set to this folder:

```r
source("00_setup.R")
source("01_simulate_data.R")
source("02_fit_models_cv.R")
source("03_mapping_uncertainty.R")
```

Outputs (saved in `outputs/`) include:
- `cv_metrics.csv`
- `cv_predictions.csv`
- `map_mean.tif`, `map_sd.tif`
- `map_mean.png`, `map_sd.png`
