# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

Sub-national prediction of micronutrient deficiency prevalence (Vitamin A, iron
anaemia, and in some countries folate/B12/zinc) using DHS/MICS survey biomarker
data and remotely sensed (GEE) environmental covariates, producing Admin-2
prevalence estimates for areas with and without survey coverage. Countries:
Gambia, Ghana, Sierra Leone, Malawi. This is an **ecological prediction model**
(area-level predictors → area-level prevalence) — explicitly not individual
risk prediction or causal inference.

Two parallel analyses share the same pipeline:
- **Primary — Admin-2 area-level small-area estimation (SAE).** Models fit
  directly at Admin-2 level on survey-weighted prevalences + GEE raster
  aggregates: area-level SuperLearner, Fay-Herriot, SL→BYM2, benchmarked
  against standard SAE methods. These are the headline results and the
  dashboard default (`R/benchmark_models.R`, the "PRIMARY ANALYSIS" section of
  `_targets.R`).
- **Sensitivity — individual-level SuperLearner.** A person-level ensemble
  aggregated up to Admin-2. Retained for comparison but not the headline
  estimator (shrinkage of a noisy individual classifier biases area means).
  Code lives in `R/sensitivity/` and `sensitivity/`; see `sensitivity/README.md`.

## Commands

```r
# Run the pipeline (fast mode, default — 10-30 min sequential):
targets::tar_make()

# Parallel (see _targets.R header for worker counts):
targets::tar_make_future(workers = 4)

# Full mode (publication-quality, 2-4 hrs) — set before invalidating/running:
Sys.setenv(PIPELINE_MODE = "full")
targets::tar_invalidate(names = "pipeline_params")   # required when switching modes
targets::tar_make()
# or from the shell:
PIPELINE_MODE=full Rscript -e 'targets::tar_make()'

# Run/inspect a single target (closest thing to "run a single test" here):
targets::tar_make(names = "sl_fit_gambia_child_vitA")
targets::tar_read(cv_perf)
targets::tar_read(area_model_gambia_child_vitA)
targets::tar_invalidate(names = "sl_fit_gambia_child_vitA")   # force a re-run

targets::tar_visnetwork()                                     # dependency graph
targets::tar_outdated()                                       # what would rerun
targets::tar_meta(fields = any_of("error"), complete_only = TRUE)  # last-run errors
targets::tar_destroy()                                        # nuke the whole cache
```

There is no test suite (no `testthat`) — correctness is checked by running
targets and inspecting `results/tables/`, `tar_meta()` errors, and (for the
dashboard) a manual `shiny::runApp("dashboard")` / `testServer` smoke check.

## Architecture

The whole analysis is one [`{targets}`](https://docs.ropensci.org/targets/)
pipeline (`_targets.R`, ~1500 lines) that dependency-tracks and caches every
step; only targets whose inputs changed re-run. Functions are split across
`R/*.R` and auto-sourced recursively via `tar_source("R/")` — adding a
function to any file under `R/` (including `R/sensitivity/` and
`R/corrected/`) makes it available to the pipeline with no explicit import.

**Targets are generated programmatically per country × outcome** (not written
out by hand): `_targets.R` loops over `get_country_configs()` and each
country's `outcomes` list (both defined in `R/config.R`) and builds target
names like `sl_fit_{country}_{outcome}`, `cv_perf_{country}_{outcome}`,
`admin2_error_{country}_{outcome}` via `tar_target_raw()` + `substitute()`.
To add a country or outcome, edit `R/config.R` and re-run `tar_make()` — new
targets appear automatically. See the README's "Adding a New Country" /
"Adding a New Outcome" sections for the config shape.

`_targets.R` is organized into clearly delimited sections (search for the
`# ── ... ──` banners): static/shared targets, per-country targets, the
area-level model (all country × outcome combos), area-level comparison /
LOCO targets, the **PRIMARY ANALYSIS** SAE benchmark suite
(`R/benchmark_models.R`), cross-country transportability (LOCO CV), Côte
d'Ivoire out-of-sample prediction, summary/rollup targets, and finally the
**CORRECTED METHODS (P1–P8)** section — a from-scratch reimplementation of
the fitting/evaluation logic (folded in from a formerly-separate
`_targets_corrected.R`, now archived) that fixes specific methodological bugs
(honest CV, no leakage, proper centering, calibration, decision-value
framing). Its functions live in `R/corrected/` and it emits
`corrected_methods_comparison`, written to `results/tables/corrected/*.csv`
and `dashboard/data/methods_comparison.rds` for side-by-side comparison
against the production numbers. Treat `R/corrected/*` as the "getting it
right, audited" layer, not a redundant duplicate — when fixing a
methodological bug, check whether the fix belongs in the corrected layer, the
main pipeline, or both.

**Key methodological pieces** (each with a dedicated `R/` module):
- Survey design: stratified cluster sampling via `srvyr::as_survey_design()`
  (PSU = cluster number, weights = survey weights).
- SuperLearner fitting: `mlr3superlearner` (`R/mlr3_fitting.R`) with
  cluster-blocked CV (PSU ids passed via `group=` so a cluster's observations
  never split across folds); discrete SuperLearner (CV picks the single best
  learner). The legacy `sl3`-based path (`R/sl_fitting.R`) is reference-only
  and not part of the active DAG.
- Out-of-fold predictions: CV performance and conformal-interval residuals use
  genuine out-of-fold predictions (`res$yhat_full`), recomputed via
  `mlr3_oof_predictions()` because `mlr3superlearner` doesn't retain internal
  CV predictions. `res$yhat_insample` is resubstitution-only, kept for
  permutation-importance baselines.
- Uncertainty: conformal prediction intervals (`R/conformal.R`) — split
  (constant-width) and locally-adaptive variants — built from out-of-fold
  residuals; no refitting needed. The older cluster-bootstrap path
  (`R/bootstrap.R`) is reference-only.
- Area-level model: elastic net (`glmnet`, alpha = 0.5) on survey-weighted
  Admin-2 prevalence ~ GEE raster zonal means — this is what lets the
  pipeline predict into Admin-2 areas with no survey data.
- Cross-country transportability (`R/transportability_area.R`): a universal,
  parsimonious within-country-centered elastic net over harmonized
  GEE + IHME + Malaria-Atlas + food-security proxies, validated by
  leave-one-country-out CV.

**Downstream consumers.** `results/tables/`, `results/figures/`, and
`results/models/` are written by summary targets and read by both the Quarto
docs (`docs/*.qmd`) and the dashboard's `data-raw/` builders. The
`dashboard/` app (Shiny, bslib + leaflet, modular under `dashboard/R/mod_*.R`)
is a separate consumer, not part of the `targets` DAG — its `data-raw/`
scripts read pipeline outputs and pre-build `.rds` bundles into
`dashboard/data/`, which `app.R`/`global.R` load at Shiny startup. When
pipeline output shapes change, the corresponding `dashboard/data-raw/*`
builder usually needs re-running before the dashboard reflects it. The
individual-level SL sensitivity results surface in the dashboard as a
labelled "sensitivity" layer (`mod_map_explorer.R`) alongside the SAE
default — see `sensitivity/README.md` for the full list of shared artifacts.

**Legacy/reference code**, not part of the active pipeline unless noted
above: `src/` (country-specific data cleaning/merging scripts — some are
still the provenance for `data/IPD/*`, so check before deleting), `archive/`
(retired files, moved not deleted — see `archive/ARCHIVE_MANIFEST.md` to
restore anything), `national_prediction/` (a separate sub-pipeline with its
own `_targets.R`), `mn_proxy_tutorial/` (self-contained teaching artifact),
`simplified subset/` (a scoped-down build for a specific use case).

## Working conventions

- Prefer editing/adding functions in `R/` over inlining logic in `_targets.R`
  — everything under `R/` is auto-sourced, so new functions need no wiring.
- When changing fitting/evaluation logic, check whether `R/corrected/` has a
  parallel implementation that also needs the fix (or is the *source* of
  truth the main pipeline should be reconciled against).
- Mode switches (`PIPELINE_MODE`) require `tar_invalidate(names =
  "pipeline_params")` or the pipeline will silently keep using cached
  fast-mode (or full-mode) results.
- Files with `_tmp_`, `_run_` prefixes or living in `archive/`/`sandbox*` are
  throwaway/retired — don't build on them without checking
  `archive/ARCHIVE_MANIFEST.md` first.
