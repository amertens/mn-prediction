# District Factors (model-agnostic SHAP) — staged, NOT yet run

This is the **Option B** batch that populates the Importance tab's **District
factors** sub-panel. It is fully set up but **has not been run**. It explains
the *actual* SuperLearner predictions (including the BART winner), not a
surrogate model.

## Why it's needed
The pipeline's built-in SHAP (`R/shap_explanations.R`) only extracts tree-SHAP
from an **xgboost** base learner. In the production fits the ensemble put weight
1.0 on **BART** and 0 on xgboost, and only the winning model was retained — so
that path returns empty (`no_xgb_learner`). This batch instead computes
model-agnostic SHAP against the real `predict()`.

## Method
Self-contained Monte-Carlo (sampling) SHAP — Strumbelj & Kononenko (2014),
vectorised over instances. **No `fastshap`/`kernelshap` dependency** (neither is
installed). It:
- predicts via `predict(fit$sl_fit$mlr3_fit, newdata)` — the working path used by
  the ablation modules (the SL wrapper's own `$predict` has a known cluster_id bug);
- explains the saved, model-ready `train_data` (already imputed + recipe-baked),
  recovering each row's district by position-aligning with the non-missing-outcome
  rows of `outcome_data$data` (with a hard length check);
- to bound cost, explains the **top-K features** (from `single_var_importance.csv`)
  and a **per-district sample** of individuals. These caps are logged, not silent.

## Files
- `shap_kernel/shap_kernel_lib.R` — SHAP estimator + per-slice compute.
- `shap_kernel/run_shap_kernel.R` — checkpointed driver (one RDS per slice in
  `results/shap_kernel/`), then combine + rebuild the dashboard bundle.
- `run_shap_kernel.ps1` — detached, single-instance, logged launcher.
- Outputs: `results/tables/shap_district_factors.csv`,
  `shap_global_importance.csv`, and a rebuilt `dashboard/data/importance.rds`.

## Cost (why it's "heavy")
Per slice the cost is ~`top_k x nsim x 2` model-prediction calls, each over the
explained rows; BART prediction dominates. With defaults (top_k=20, nsim=20,
~400 explained rows) expect on the order of **minutes per slice**, i.e. roughly
**0.5–2 hours** for all ~24 slices — but verify with the smoke run first.
The machine should be otherwise idle (no competing heavy R jobs).

## How to run (unattended, resumable)

1. **Smoke test one slice first** (seconds, validates the whole path end-to-end):
   ```
   SHAP_SMOKE=1 SHAP_ONLY=gambia_women_iron \
     Rscript --vanilla shap_kernel/run_shap_kernel.R
   ```
   Confirm it prints `[ok ] gambia_women_iron ...` and writes a checkpoint.

2. **Full run, detached + logged** (PowerShell):
   ```powershell
   Start-Process -WindowStyle Hidden powershell.exe `
     -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File',`
     'C:\Users\andre\OneDrive\Documents\mn-prediction\run_shap_kernel.ps1'
   ```
   Or directly: `Rscript --vanilla shap_kernel/run_shap_kernel.R`.

3. **Resume:** just re-run the same command — completed slices are skipped
   (their checkpoint RDS exists in `results/shap_kernel/`). Delete a slice's RDS
   to force recompute.

4. **Tuning** via env vars: `SHAP_NSIM`, `SHAP_TOPK`, `SHAP_NPD`
   (per-district sample), `SHAP_NMAX` (max explained rows), `SHAP_BG`
   (background size). Raise `SHAP_NSIM` for lower-variance SHAP at higher cost.

## After it finishes
`dashboard/data/importance.rds` is rebuilt with real `shap` rows, so the
**District factors** panel populates automatically. Then redeploy:
`Rscript dashboard/deploy.R`. Consider updating the panel's methods note to say
the attributions are model-agnostic (sampling SHAP over the full ensemble),
replacing the current xgboost-specific wording.

## Status
- Scripts written and **syntax-checked**; **not executed** (per request).
- Not validated end-to-end yet — run the smoke step (1) first; it is cheap and
  confirms the predict path + district alignment before the full batch.
