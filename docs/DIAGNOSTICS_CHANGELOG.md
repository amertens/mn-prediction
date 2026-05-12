# Post-Hoc Diagnostics Upgrade

**Date:** 2026-03-25
**Scope:** Lightweight diagnostics and reporting from existing predictions only.
**No SL refitting required.**

## What was added

### New file: `R/diagnostics.R`

Computes from existing cross-validated predictions (`sl_fit$res$Y` and `sl_fit$res$yhat_full`):

- **PR-AUC** (precision-recall area under curve) via ROCR — more informative than ROC-AUC for rare outcomes
- **Brier score** and **Brier skill score** (vs. prevalence-only baseline)
- **Calibration intercept and slope** via logistic recalibration
- **Binned calibration tables** (decile-level predicted vs observed)
- **PR curve data** for plotting
- **ROC curve data** for plotting
- **Prevalence and event counts** per model
- Continuous model diagnostics: RMSE, MAE, R2, prevalence recovery (observed vs predicted at cutoff)

### New targets in `_targets.R`

- `diagnostics_{country}_{outcome}` — per-model diagnostics (cheap, reads from cached `sl_fit`)
- `diagnostics_all` — combined summary (binary metrics, continuous metrics, calibration tables, PR/ROC curves)
- Diagnostics tables saved to `results/tables/diagnostics_binary.csv` and `diagnostics_continuous.csv`
- Added `ROCR` to pipeline packages

### Updated: `docs/results_report.qmd`

New section "Model Diagnostics: Beyond AUC" (#sec-diagnostics) with:

- Plain-language explanation of discrimination vs calibration vs prevalence recovery
- Extended diagnostics summary table (n, events, prevalence, ROC-AUC, PR-AUC, Brier, Brier skill, calibration slope/intercept)
- Faceted reliability (calibration) plots
- Faceted precision-recall curves with no-skill baselines
- Continuous model diagnostics with prevalence recovery
- Concrete strongest/weakest model examples with interpretation

### Updated: `docs/mn_prediction_slides.qmd`

New slide "Evaluating Models for Rare Outcomes" with:

- Three-dimension framework (discrimination, calibration, prevalence recovery)
- Table of key metrics and what they measure
- Explanation of why standard AUC misleads for rare outcomes
- Speaker notes on area-level vs individual-level interpretation

## What was NOT changed

- No SL models were refit
- No learner library changes
- No new ablation analyses
- No permutation importance
- All diagnostics computed from existing cached `sl_fit$res` predictions
