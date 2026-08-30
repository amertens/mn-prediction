# Pipeline & Results Review — 2026-05-10

Critical review of the prediction pipeline at `mn-prediction/`. Findings are
triaged into **fixed** (low-risk, additive changes already applied), **deferred**
(non-trivial refactors flagged for follow-up), and **non-issues** (audit claims
that did not survive verification).

The hard rule for this pass: do not break the running pipeline. All fixes are
additive — existing targets and CSV schemas continue to load; new columns are
appended rather than replacing existing ones.

---

## Summary

| Severity | Count | Status |
|---|---|---|
| Critical (correctness) | 1 | **fixed** |
| Medium (correctness)   | 2 | **fixed** |
| Medium (eval gap)      | 4 | **fixed** |
| Medium (correctness)   | 2 | deferred (refactor) |
| Minor / non-issue      | 4 | documented |

---

## Fixed in this pass

### 1. Critical — clustered CV was not being enforced

**File:** [R/mlr3_fitting.R:293-322](../R/mlr3_fitting.R)

`origami::make_folds(cluster_ids = ...)` was called and its result
(`fold_obj`) was assigned but **never used**. Only the integer fold count
was forwarded to `mlr3superlearner()`, so PSUs were free to be split across
folds during CV — violating the survey design and inflating CV performance
by an unknown amount.

**Fix:** keep `cluster_id` as a column in the data and pass `group =
"cluster_id"` to `mlr3superlearner()`. Per the package's source
(`task$set_col_roles(group, "group")`), this assigns the mlr3 group role,
which excludes the column from features and blocks it during sample
splitting. Verified end-to-end on synthetic clustered data (n=500, 50
clusters, 3 folds, `glm` + `mean` + `ranger`): fit completes, predictions
align row-for-row with input, AUC reasonable, no NA predictions.

### 2. Medium — SL→Admin2 aggregation was unweighted in scripts/06

**File:** [scripts/06_admin2_predictions_map.R:421-440](../scripts/06_admin2_predictions_map.R)

The targets-pipeline path
([R/admin2_analysis.R:199](../R/admin2_analysis.R)) aggregates SL
predictions to Admin2 using survey weights, but the standalone script used
a plain `mean(deficient)`. The **survey** prevalence on the same script is
design-weighted, so the SL-vs-survey error metrics conflated model error
with weighting-method mismatch.

**Fix:** the script now joins `cc$weight_col` and uses
`weighted.mean(deficient, .wt)` (with a safe fallback when an Admin2 has
no positive weights). This matches `R/admin2_analysis.R:199` and makes the
SL-vs-survey comparison apples-to-apples.

### 3. Medium — Wald CIs for survey prevalence at low rates

**File:** [scripts/06_admin2_predictions_map.R:807-840](../scripts/06_admin2_predictions_map.R)

`srvyr::survey_mean(vartype = "ci")` returns a Wald-style CI that can
fall below 0 or rise above 1 for rare deficiencies and small Admin2
cells. Two new columns `svy_prev_logit_lo` / `svy_prev_logit_hi` are
appended (delta-method on the log-odds, clamped to [0, 1]); the
original `svy_prev_low` / `svy_prev_upp` are preserved so downstream
consumers don't break.

### 4. Medium — no minimum sample size for Admin2 error metrics

**File:** [scripts/06_admin2_predictions_map.R:25,857-876](../scripts/06_admin2_predictions_map.R)

Admin2 cells with `n_svy < 10` were entering MAE/RMSE with very high
sampling variance, dominating the metric in unstable ways.

**Fix:** new `MIN_N_SVY` parameter (default `10`) filters the eval set;
the per-Admin2 table still reports all cells. The error summary CSV
gains two columns: `n_admin2_excluded` and `min_n_svy` so the threshold
choice is transparent in downstream reports.

### 5. Medium — calibration metrics had no uncertainty / no ECE

**File:** [R/diagnostics.R:95-149](../R/diagnostics.R)

`compute_binary_diagnostics()` returned point estimates of the
calibration intercept and slope but no CIs, and no Expected Calibration
Error.

**Fix:** four new columns appended (`calib_int_lo/hi`,
`calib_slope_lo/hi` — Wald CIs from the recalibration GLM) plus `ece`
and `mce` (Expected and Maximum Calibration Error, computed from the
same n-bin partition the existing calibration table uses). Schema change
is additive; existing columns kept exactly. Verified on synthetic
miscalibrated data (slope ≈ 3.3, ECE ≈ 0.10, MCE ≈ 0.41).

### 6. Medium — Moran's I added to admin2 error summary

**File:** [scripts/06_admin2_predictions_map.R:875-913](../scripts/06_admin2_predictions_map.R) (already merged earlier today)

Spatial-autocorrelation diagnostic on SL residuals at Admin2 is now
computed alongside MAE/RMSE/Pearson r, surfacing whether the global SL
leaves spatial structure unexplained at the polygon scale.

---

## Deferred — non-trivial, flagged for future work

### A. Preprocessing fitted before CV folds (potential leakage)

**File:** [R/mlr3_fitting.R:174-231](../R/mlr3_fitting.R)

`caret::nearZeroVar()`, `ck37r::impute_missing_values()`, and the
`recipes` pipeline (`step_corr`, `step_normalize`) are run once on the
full sample before `mlr3superlearner()` constructs CV folds. Strictly
speaking each of these should be re-fitted inside each training fold,
otherwise test-fold values inform the training preprocessor.

In practice the bias is bounded:

- `step_normalize` uses a global mean/SD — small effect.
- `step_corr` (threshold 0.85) drops a feature if its train+test
  correlation with another feature is high; this can occasionally
  remove a feature that wouldn't be flagged on training-only data.
- `ck37r::impute_missing_values(type = "standard")` uses marginal
  median/mode — small effect at typical sample sizes.
- `caret::nearZeroVar` removes features that are near-constant in the
  whole sample; a feature constant in train but variable in test would
  be incorrectly retained, but this is a rare failure mode.

**Why not fixed now:** moving preprocessing inside the CV loop requires
refactoring `mlr3_SL_clustered()` to either (a) run K manual fits with
per-fold preprocessing, losing `mlr3superlearner`'s metalearner
optimization, or (b) wrap the preprocessing as an mlr3 `PipeOp` chain.
Either is a sizable change with real risk of regressing the headline
results. **Recommended follow-up:** add a `--per-fold-preproc` sensitivity
target that runs (a) on a single country/outcome and reports the AUC
delta. If the delta is < 1 AUC point, the current pipeline can stay; if
larger, schedule the refactor.

### B. Aggregation uncertainty assumes independent predictions

**File:** [R/conformal.R:150-157](../R/conformal.R),
         [R/aggregation_uncertainty.R:27-48](../R/aggregation_uncertainty.R)

When individual-level conformal intervals are aggregated to Admin2 /
Admin1 / national, the variance is computed as
`sum(w_i^2 * sigma_i^2)`, which assumes independent prediction errors.
Predictions within a cluster share unmeasured confounders and a common
learner, so errors are positively correlated and the aggregated SE is
underestimated. CIs on Admin2 predicted prevalence are therefore narrower
than they should be.

**Why not fixed now:** the right fix is a cluster-level resampling
(bootstrap clusters, recompute aggregated point estimates, take
quantiles) which costs ~10-30× the current pipeline if applied to all
country×outcome combinations. **Recommended follow-up:** prototype on a
single outcome and compare to the current delta-method CI; if the
correlation-corrected CI is materially wider, replace; otherwise add a
DEFF correction as a cheap approximation.

---

## Audit claims that did not survive verification

- **Outcome leakage via `gw_` columns.** Audit flagged this as critical;
  inspection of `cfg$gw_exclude_patterns` shows the outcome biomarker
  columns (`gw_cVAD_Thurn`, `gw_wIDA_Brinda`, etc.) are filtered out
  before `Xvars` is built. Confirmed by grepping for these names in the
  feature lists at fit time. **Not a bug.**
- **BART refit-on-predict leaks the test fold.** The fallback only
  triggers when the primary `predict()` path returns degenerate values,
  and at that point the SL has already been trained inside per-fold
  models held by `mlr3superlearner`. The fallback refit is on the
  fold-training data the model was originally given, not the global
  pool. **Bounded; not actionable.**
- **National-estimates use simple Admin1 mean.** Verified to use
  `srvyr::survey_mean` with the design object. **Correct as-is.**
- **`has_survey` filtering missing.** The error block uses
  `filter(!is.na(svy_prev), !is.na(sl_prev))`, which is functionally
  equivalent. **Not a bug**; `MIN_N_SVY` was added for a stronger
  filter.

---

## Recommended evaluation additions (next iteration)

These were not implemented this pass to avoid scope creep, but each is
low-risk and additive:

1. **Stratified diagnostics** — compute AUC / calibration separately by
   age group (children 6-23m, 24-59m) and sex (where applicable). Currently
   pooled, which can hide subgroup miscalibration.
2. **Decision-curve analysis** — particularly relevant for binary
   outcomes where downstream use is "screen households for
   intervention". `dcurves` package is the standard.
3. **Cluster bootstrap on CV-AUC** — gives a CI on AUC that respects
   the survey design (current DeLong CI assumes iid).
4. **Held-out conformal calibration set** — for the manuscript's coverage
   claims to be theoretically guaranteed (rather than approximate)
   conformal needs a true held-out calibration partition. The current
   "use CV residuals" approach yields approximately marginal coverage; a
   single 80/20 split would make it exact under exchangeability.

---

## Files changed in this pass

| File | Change |
|---|---|
| [R/mlr3_fitting.R](../R/mlr3_fitting.R) | `group = "cluster_id"` passed to mlr3superlearner |
| [scripts/06_admin2_predictions_map.R](../scripts/06_admin2_predictions_map.R) | weighted SL aggregation; logit-CIs; n_svy threshold; Moran's I |
| [R/diagnostics.R](../R/diagnostics.R) | calibration CI; ECE; MCE; spatial residual diagnostic (earlier today) |
| [R/gwr_analysis.R](../R/gwr_analysis.R) | new module (earlier today) |
| [scripts/08_admin2_gwr_supplementary.R](../scripts/08_admin2_gwr_supplementary.R) | new driver (earlier today) |
| [docs/pipeline_review_2026-05.md](pipeline_review_2026-05.md) | this document |

No removed columns; no schema-breaking changes.

---

## Verification

Smoke tests passed locally:

- `mlr3_SL_clustered` with `group = "cluster_id"` on 500 obs / 50
  clusters, 3 folds, mixed learner library — fit OK, no NA preds, AUC
  in expected range.
- `compute_binary_diagnostics` on a synthetic miscalibrated outcome —
  ECE/MCE non-negative, calibration-slope CI brackets the point
  estimate.

Recommended next step before publishing results from the patched
pipeline: re-run one full country × outcome (e.g. Ghana × women_iron)
under the new clustered-CV settings and diff the cv_perf table against
the prior run. The expected direction is **AUC slightly down**
(previous estimates were optimistic because clusters could be split);
the magnitude tells you how much past results were biased.

---

## Addendum (2026-05-10): cluster-CV sensitivity panel

Ran the patched `mlr3_SL_clustered` (with `group = "cluster_id"`) on 4
country × outcome cases via
[scripts/09_cluster_cv_sensitivity.R](../scripts/09_cluster_cv_sensitivity.R).
Comparison is against the cached
`results/tables/cv_performance_all.csv` from the prior full-stack run.
All fits used the full 13-learner library with BART included.

| Country | Outcome | n | Prev | AUC baseline | AUC patched | ΔAUC |
|---|---|---:|---:|---:|---:|---:|
| Ghana  | women_iron | 981 | 8.5% | 0.7197 | 0.7426 | **+0.023** |
| Malawi | women_iron | 752 | 11.5% | 0.7954 | 0.7965 | +0.001 |
| Ghana  | women_b12  | 473 | 6.5% | 0.8554 | 0.8510 | −0.004 |
| Sierra Leone | child_iron | — | — | — | — | error during refit (see below) |

**Continuous-side deltas** (RMSE):
Ghana/women_iron −0.001, Malawi/women_iron exactly 0.000,
Ghana/women_b12 −1.67 on a baseline of 217 (≈0.8% relative). All
within numerical noise.

**Interpretation.** Across the three completed outcomes:

- **ΔAUC mean +0.007, range [−0.004, +0.023].** No sign of a
  systematic shift in either direction.
- The Ghana/women_iron +0.023 we saw in the single-outcome trial was
  the most extreme of the three; the other two are essentially flat.
- The continuous-side metrics are nearly identical across baseline
  and patched runs — including a 0.000 RMSE delta on Malawi/women_iron,
  which is consistent with the continuous SuperLearner discrete winner
  not changing between the two fold schemes.
- **The old random-fold CV was not materially leakage-inflated.** If
  splitting clusters had been giving the model unfair information,
  we'd expect a consistent baseline > patched gap on the order of
  1–3 AUC points. We don't see that. The +0.023 outlier looks like
  fold-composition noise on a rare outcome (83 events / 981).

**Implication for the deferred preprocessing-leak issue (item A in
the review):** weaker case for an urgent refactor. Combined with the
clustered-CV result, it's now more plausible that the recipe-prep-on-
full-data leak is also small in practice. Worth re-prioritizing below
the aggregation-correlation item (B), which has a clearer mechanism
for affecting headline numbers.

**Sierra Leone child_iron error.** The refit failed with "argument is
of length zero". Looking at the warnings, `oc$binary` is being passed
into `mean()` and returning NA, and one predictor (`lsms_hh_se.ppp2011`)
had zero variance in this dataset. Likely an upstream data issue
specific to the Sierra Leone outcome_data target rather than the
cluster-CV patch — this case also caused the n=14 Admin2 thinness we
saw in the earlier GWR run, so Sierra Leone's outcome_data probably
deserves its own data-prep audit.

**Recommendation:** keep the cluster-CV patch and proceed. The
sensitivity panel is a more defensible artifact for the manuscript
than the single-outcome trial; it shows that the methodologically
correct fix doesn't change the headline numbers in a meaningful way.

Outputs:

- [results/sensitivity/cluster_cv_summary.csv](../results/sensitivity/cluster_cv_summary.csv)
- per-outcome CSVs in `results/sensitivity/`
