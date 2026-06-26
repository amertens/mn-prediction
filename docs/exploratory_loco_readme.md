# LOCO transportability sandbox

A scratch workspace for **out-of-the-box** approaches to admin-2 micronutrient
deficiency prediction across countries. Outside the `targets` pipeline so we
can iterate fast without invalidating production targets.

## Goal

Find a **parsimonious** model that **transports across countries** — i.e.,
beats the current best LOCO Pearson r (~0.5 for iron, ~0.3 for vitamin A)
without needing the 17-hour SuperLearner.

## What's here

- `00_setup.R` — load data, common helpers, the LOCO evaluator used by every
  hypothesis script
- `01_baseline.R` — re-run current best methods (penalized, mixed, two_stage,
  forward, gam, betareg) to confirm the baseline
- `02_transferability_filter.R` — H1: keep only features whose
  univariate correlation with the outcome is *stable across training
  countries*
- `03_synthetic_control_knn.R` — H2: kNN in covariate space, no
  regression at all
- `04_spatial_gam_coords.R` — H3: thin-plate spline on lat/lon
  directly, ignore most covariates
- `05_outcome_shared_features.R` — H4: features that predict ALL 4
  outcomes well — borrow strength across outcomes
- `06_domain_adversarial.R` — H5: train features to be predictive of
  outcome but NOT predictive of country
- `findings_log.md` — running log of what worked and didn't

## How to run

Each script is standalone. Source `00_setup.R` first (the others source it).

```r
source("sandbox/00_setup.R")
source("sandbox/01_baseline.R")   # baseline numbers
source("sandbox/02_transferability_filter.R")
# ... etc
```

Results land in `sandbox/results/` as CSVs.

## Ground rules

- No DHS-derived features (`gw_*`, `dhs_*`) in LOCO evaluation — held-out
  countries lack DHS by construction. All experiments use GEE-only
  predictors.
- Headline metric is **Pearson r between predicted and observed admin-2
  prevalence** under leave-one-country-out. MAE / bias reported alongside.
- Compute should stay under ~5 minutes per script; if a method takes
  longer, it's not parsimonious enough for the parsimony test.
