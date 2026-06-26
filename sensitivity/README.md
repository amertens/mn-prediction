# Sensitivity analysis — individual-level SuperLearner

The **primary analysis** for this project is the Admin-2 **area-level small-area
estimation (SAE)**: models fit directly at the Admin-2 level on survey-weighted
prevalences and earth-observation aggregates (area-level SuperLearner,
Fay–Herriot, SL→BYM2), benchmarked against standard SAE methods. Those are the
headline results, the manuscript's primary tables, and the dashboard default.

This folder collects the **individual-level SuperLearner** analysis, which is
now a **sensitivity analysis**: a person-level ensemble that predicts each
individual's deficiency status from proxy covariates and is then aggregated up
to Admin-2. It is retained because it shares the same inputs as the primary
analysis and provides a useful contrast — but aggregating person-level
predictions introduces systematic bias (shrinkage of a noisy individual-level
classifier propagates to the area means), so it is **not** the headline
estimator.

See the call rationale (2026-06-15 BMGF bi-weekly): SAE methods tend to match or
beat the SuperLearner for area-level prediction and transportability, and the
recommended direction is to combine them (SL predictions feeding an SAE model,
e.g. SL→BYM2) rather than to lead with the individual-level SL.

## Where the code lives

| Component | Location | Notes |
|---|---|---|
| SL fitting engine | [`R/sensitivity/sl_fitting.R`](../R/sensitivity/sl_fitting.R), [`R/sensitivity/mlr3_fitting.R`](../R/sensitivity/mlr3_fitting.R) | Still auto-sourced by `tar_source("R/")` (recursive), so the pipeline is unaffected by the move. |
| Pipeline targets | `_targets.R` → "DYNAMIC TARGET FACTORY" | `sl_fit_*`, conformal CIs, ablations, diagnostics, SHAP, national estimates, and `aggregate_admin2_sl`. Marked as the sensitivity analysis in-file. |
| Standalone scripts | [`09_cluster_cv_sensitivity.R`](09_cluster_cv_sensitivity.R) | Cluster-CV sensitivity for the individual-level SL. |
| Dashboard layer | `dashboard/R/mod_map_explorer.R` | "Individual SuperLearner — sensitivity (surveyed districts)" map layer; the SAE layers are the default. |

## Results produced by this analysis

The individual-level SL still produces several artefacts consumed elsewhere
(e.g. `dashboard/data/admin2_predictions.rds`, and the per-outcome rows in
`results/tables/{diagnostics_*,national_estimates_all,single_var_importance,
domain_ablation_all}.csv`). These were intentionally **not** relocated, because
the dashboard and manuscript still read them for the side-by-side
primary-vs-sensitivity comparison. They are labelled as sensitivity outputs in
those consumers rather than physically moved.
