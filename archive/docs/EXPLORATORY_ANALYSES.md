# Exploratory analyses: findings index (2026-06)

The project ran two scratch workspaces outside the `targets` pipeline to iterate
quickly on transportability and preprocessing. They have been retired: their
**findings writeups are preserved in `docs/`** (linked below), their **scripts are
archived** under `archive/sandbox/` and `archive/sandbox_fe/` (and remain in git
history), and their large regenerable caches and per-experiment result CSVs were
removed. Nothing learned was lost; what follows is the index of where each finding
now lives and what was promoted into the production pipeline.

## `sandbox/` — LOCO transportability experiments
Goal: find a parsimonious model that transports across countries (better
leave-one-country-out, or LOCO, Pearson r) without the 17-hour SuperLearner.
Full writeups:

- [`exploratory_loco_findings.md`](exploratory_loco_findings.md) — the running
  findings log (method-by-method LOCO results, bootstrap CIs, promotions).
- [`exploratory_loco_readme.md`](exploratory_loco_readme.md) — the hypothesis
  map (H1–H6) and ground rules.

Headline findings:
- The strongest LOCO model is the simplest: a **thin-plate spline on admin-2
  polygon centroids** (mean r ~ 0.285), beating all 152 earth-observation
  covariates. Adding the **top-5 univariate-LOCO SoilGrids predictors** (soil
  Al / Ca / Mg, linear) lifts it to mean r ~ 0.33 — the best parsimonious model
  found. Adding more features degrades it (small-n degrees-of-freedom).
- Soil micronutrient content is the closest covariate to a causal pathway (soil
  nutrient → crop nutrient → intake → biomarker); geography carries the rest.
- Bootstrap CIs show the wins concentrate in iron and vitamin A for Gambia and
  Ghana; Malawi and low-prevalence vitamin A do not transport. Report point
  estimates with CIs, not clean scalars.
- **Promoted to the pipeline:** `fit_predict_invariance_filter` (method 15) and
  `fit_predict_combined_filter` (method 16) in `R/benchmark_models.R`, plus
  `bootstrap_loco_ci()` extensions. The spatial-spline + soil recipe is the
  manuscript's parsimonious result; the all-methods comparison figure/data are
  now at `results/figures/all_methods_child_iron.png` and
  `results/tables/all_methods_child_iron.csv`.

## `sandbox_fe/` — feature engineering & small-area estimation
Goal: learn what preprocessing the data structure actually rewards, evaluated
under honest leakage-free CV. Full writeups:

- [`exploratory_fe_findings.md`](exploratory_fe_findings.md) — preprocessing and
  feature-engineering investigation.
- [`exploratory_fe_sae_findings.md`](exploratory_fe_sae_findings.md) — the
  Fay-Herriot small-area-estimation prototype.

Headline findings:
- **Effective sample size is the number of areas (14–87), not individuals.** The
  predictors are almost all admin-2-level, so ~100% of their variance is
  between-area. This is the root cause of weak and unstable transfer, and it
  motivates modelling at the area level.
- **Rank/quantile normalization beats z-scoring** (+0.01–0.02 AUC, free);
  aggressive dimensionality reduction (~27 features match the full ~1000 within
  country); median imputation is sufficient; unsupervised PCA loadings are **not**
  transportable (components flip sign across countries).
- **Select predictor domains by outcome biology:** vitamin A transfers on
  environment (GEE); iron on malaria (MAP) + modelled health burden (IHME).
  Pooling all domains dilutes the signal.
- **SAE prototype:** admin-2 is below the surveys' design resolution (1–3 clusters
  per area), so direct area estimation is non-viable and a model-based approach
  (Fay-Herriot / multilevel) is mandatory; there is a large in-sample-vs-leave-
  area-out optimism gap (report leave-area-out CV); FH borrows strength correctly
  and tightens intervals 19–33%, though honest area intervals remain wide.
- These findings shaped the corrected (P1–P8) methods now in the pipeline
  (in-fold preprocessing, cluster/spatial-block CV, partial-pooling SAE,
  sampling-error-aware metrics) and are recorded in the project memory.

## Old pipeline artifacts removed
- **`_targets_fast_backup/`** (729 MB) — a stale backup of the fast-mode targets
  store, superseded by `_targets_full/` (the analysis of record) and fully
  regenerable; removed.

## Flagged, intentionally NOT removed
- **`src/`** — the legacy per-country data-cleaning and merge scripts. These are
  the upstream data-preparation provenance that regenerates the pipeline's input
  tables, and `src/0-functions.R` is still sourced by `_targets.R`. There are
  clear within-`src` versioned duplicates (`*_V2.R`) worth a dedicated, reviewed
  pruning pass, but the tree must not be deleted wholesale.
- **`national_prediction/`** — a separate sub-pipeline (own `_targets.R`) whose
  `clean_brinda_data.rds` is used by `docs/brinda_vmnis_loco_validation_plan.qmd`
  and `scripts/pull_vmnis_validation.R`; clarify active/superseded status before
  any removal.
