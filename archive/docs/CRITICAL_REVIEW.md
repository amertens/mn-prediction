# Critical review — micronutrient deficiency prediction pipeline & dashboard

**Date:** 2026-06-14
**Scope:** `R/` modeling code, `_targets.R` factory, `dashboard/`.
**Goals:** (1) improve the models; (2) communicate *if / when / how* the models
are actually useful to **non-technical decision-makers** — beyond MSE/AUC.
**Method:** code reading + methodological critique + lightweight prototypes from
*existing* results (no retraining; a heavy job is running on this machine).

> Convention: findings are tagged **[LEAKAGE] / [CORRECTNESS] / [METHODS] /
> [CODE] / [COMMS]** and reference `file:line` or pipeline target names.
> Severity = how much it can mislead a decision; Effort = S/M/L.

---

## TL;DR — headline findings

1. **Within-country cross-validated performance is optimistic by construction.**
   In `R/mlr3_fitting.R` all preprocessing — imputation, near-zero-variance
   pruning, **supervised feature prescreening against the full outcome
   (`washb_prescreen`, line 307)**, correlation filtering and normalization
   (line 341) — is fit on the **entire dataset before CV folds are created
   (line 409)**. The cluster-blocked SuperLearner (`group="cluster_id"`) then
   cross-validates on already-leaked features. So `cv_perf`, `diagnostics_*`
   (ROC/PR/Brier-skill/calibration) overstate real performance. The supervised
   prescreen is the worst offender (it uses `Y` directly). **This does not
   affect the LOCO/transport path** (which calls the same function on
   training-country data only) — which is precisely why the honest
   leave-one-country-out numbers are so much worse than the within-country
   ones.

2. **Calibration is assessed in-sample, so the "recalibration fixes it" story
   is overstated.** The Platt recalibrator is fit and evaluated on the *same*
   out-of-fold predictions (`R/benchmark_models.R:~1957`, `fit_binary_recalibrator`;
   `R/diagnostics.R:100-112`). Reported `calib_slope≈1`, ECE/MCE and the Brier-skill
   recovery in `diagnostics_binary_calibrated.csv` are optimistic. A
   fold-aware routine (`calibrate_binary_oof`) exists in the codebase but is
   **not the one called**.

3. **Uncertainty intervals likely under-cover.** Conformal intervals
   (`R/conformal.R`) use CV residuals without a held-out calibration split, are
   computed at Admin-1 and **broadcast to Admin-2** (no district-level coverage
   guarantee), and the delta-method aggregation treats observations as i.i.d.,
   **ignoring the survey design** (strata/PSU) → intervals too narrow where it
   matters most (sparse districts).

4. **The model's real decision value is modest and wildly heterogeneous — and
   the dashboard never showed it.** Using the honest out-of-sample (LOCO)
   predictions, model-guided targeting of the worst 20% of the population
   ranges from **~1.9× better than no targeting (child_iron, Ghana)** to
   **worse than chance (women_vitA: Malawi 0.15×, Sierra Leone 0.32×)**. A good
   AUC hid this. The new **Decision value** dashboard tab makes it explicit.

5. **"Truth" carries sampling error that is never accounted for.** Every error
   metric (`admin2_error_all`, `benchmarks_all`, `area_comparison_all`) compares
   predictions to the *direct survey estimate*, whose own variance
   (`svy_prev_se`, ~6–10 pp for n≈50 districts) inflates apparent error and
   biases correlations downward. Sparse districts look "wrong" when they are
   merely noisy.

---

## Priority matrix (do these first)

| # | Recommendation | Part | Impact | Effort |
|---|----------------|------|--------|--------|
| P1 | Move ALL preprocessing (impute/NZV/**prescreen**/corr/normalize) *inside* the CV folds; re-report within-country `cv_perf`/`diagnostics_*` | A/B | ★★★★★ | M |
| P2 | Make calibration & Platt assessment **out-of-fold** (use `calibrate_binary_oof`); re-issue `diagnostics_binary_calibrated.csv` | A/B | ★★★★ | S–M |
| P3 | Lead the dashboard & manuscript with **decision-value** metrics (targeting/enrichment, lift, fitness-for-purpose) not AUC/MSE | D | ★★★★★ | S (prototyped) |
| P4 | Fix interval coverage: split-conformal + **survey-design SEs** (`srvyr`) at Admin-1/2; stop broadcasting Admin-1 CI to Admin-2 silently | A/B | ★★★★ | M |
| P5 | Report all accuracy metrics against **sampling-error-aware** truth (denoised RMSE, or weight/stratify by `n_svy`) | A/B | ★★★ | S–M |
| P6 | Add **out-of-support / trust flags** (covariate-support distance, wide-interval, poor-transport) and surface on the map | C/D | ★★★★ | M |
| P7 | Partial pooling / hierarchical area model (proper SAE) as the primary estimator, with design-based direct estimates as anchors | C | ★★★★ | L |
| P8 | Fix latent leakage in two-stage centering (`benchmark_models.R:577-579`) and OOS imputation fallback (`oos_prediction.R:189`); join on `(country, Admin2)` everywhere | A | ★★★ | S |

---

## A. Codebase review

### A1. Structure & reproducibility (good foundation)
- `_targets.R` is a clean dynamic factory (country × outcome) with `format="file"`
  tracking on the merged `.rds` (lines ~398-417) so on-disk changes invalidate
  downstream — good practice.
- Outcome harmonization to uniform WHO cutoffs before pooling
  (`R/admin2_analysis.R:317-328`, `UNIFORM_TRANSPORT_TAGS`) is correct and
  important for cross-country comparability.
- CV is genuinely **cluster-blocked** (`origami::make_folds(cluster_ids=...)`,
  `R/mlr3_fitting.R:409`; `group="cluster_id"` passed to mlr3superlearner) —
  the right call for survey data.
- **Risk:** two parallel target stores (`_targets/` fast, `_targets_full/` full);
  the dashboard reads `_targets_full`. Document which store is "publication"
  truth and add a guard so stale-store mixing can't happen silently.
- **Risk:** many `_tmp_*.R` / `scripts/*` exist outside the targets graph
  (e.g. `_tmp_loco_cmp.R`, `scripts/extend_loco_micronutrients.R`,
  `results/tables/_loco_new_partial.rds`). Results that feed the dashboard but
  are produced *outside* `tar_make()` are a reproducibility hole. Fold them into
  targets or clearly mark them as ad-hoc.

### A2. Data leakage (highest priority)
- **[LEAKAGE] Preprocessing & supervised prescreen on full data before CV** —
  `R/mlr3_fitting.R:282-343` (prep at 341) then folds at 409. **`washb_prescreen(Y=Y, Ws=cov, …)` at line 307 selects features using the full outcome.** Within-country CV/diagnostics are optimistic. *Fix:* wrap impute+NZV+prescreen+corr+normalize in a `recipes`/`mlr3pipelines` step that is fit per training fold (mlr3 `PipeOp` so it is cross-validated honestly), or compute OOF predictions with a pipeline graph. **Impact ★★★★★.**
- **[LEAKAGE] Class weights from full-data prevalence** — `R/mlr3_fitting.R:375`
  computes `scale_pos_weight`/`class.weights` from the whole sample before CV.
  Minor vs the above; fold-internal computation closes it.
- **[LEAKAGE] (latent) Two-stage centers the held-out country by its OWN mean** —
  `R/benchmark_models.R:577-579` (`country_means_te` from `test`). If
  `two_stage` is ever in the active benchmark set this leaks the held-out
  level. *Fix:* center test by the training-pooled mean. (Default
  `area_transport` uses `center=FALSE`, `R/transportability_area.R:38`, so the
  shipped numbers are safe — but the code path is a trap.)
- **[LEAKAGE] OOS imputation fallback uses OOS medians** —
  `R/oos_prediction.R:189`: if `area_model$X_train` is absent, CIV features are
  imputed with CIV's own medians. Always persist `X_train`; warn on fallback.
- **Cleared (good):** the LOCO/transport individual path fits the recipe on
  `train_data` only (`R/transportability.R:194-214`); area-LOCO imputes test
  with **training** medians (`R/benchmark_models.R:95-106`) and screens on
  training only (`:2117`).

### A3. Correctness risks that can silently mislead
- **[CORRECTNESS] Cross-country `Admin2` name collisions.** Joins keyed on
  `Admin2` alone (`R/oos_prediction.R:462-465`; and I reproduced a hard
  **segfault / cartesian blow-up** when merging the pooled LOCO frame to
  population by `Admin2` only). District names are not unique across countries.
  *Fix:* key every join on `(country, Admin2)`.
- **[CORRECTNESS] Missing covariates padded with 0** — `R/oos_prediction.R:175-183`.
  Zero is a real value after normalization-elsewhere; padding shifts the linear
  predictor. Warn when >5% of model vars are missing; impute to training median.
- **[CORRECTNESS] Silent clip to [0,1]** — `R/benchmark_models.R:109` (`.bound01`).
  No diagnostic when many predictions fall outside the simplex (a model-spec
  smell). Log the clip fraction.
- **[CORRECTNESS] Silent row drops / fallbacks** — NA-outcome rows dropped with
  no count (`R/mlr3_fitting.R:399`); `washb_prescreen` failure silently keeps
  all vars (`:309`); `srvyr` SE failure returns NA SE while the point estimate
  still prints (`R/national_estimates.R:47-78`). Each should emit a visible flag.
- **[CODE] gw_ exclusion is an unanchored `grepl` blacklist** —
  `R/config.R:103-107` & `R/data_prep.R:382`. Substring patterns
  (`"Fol"`, `"RBP"`, `"Zinc"`) risk both over- and under-matching. Anchor the
  patterns or switch to an explicit allow-list of proxy predictors. (Project
  rule: primary models must be proxy-only / exclude `gw_` survey variables —
  verify this holds after anchoring.)

### A4. Code quality (lower priority)
- `R/benchmark_models.R` is 2,391 lines — split SAE methods, calibration, and
  the LOCO driver into separate files.
- Hardcoded `K=5` folds for every country/outcome (`R/config.R:507,519`); some
  outcomes have <500 usable rows. Make `K` data-adaptive.
- Imputation differs between paths (ck37r vs `apply(median)`); unify.
- Add a per-run "failure summary" (which method × hold-out × outcome returned
  NULL) instead of scattered `cat()`.

---

## B. Methods review

### B1. SuperLearner & library
- Library is a reasonable ensemble (mean + lasso + ranger ×2 + lasso→xgb + xgb,
  per `_targets.R` header). **But:** (a) the meta-learner weights and the
  reported within-country performance share the same folds *and* leaked
  features (A2) — not nested/honest; (b) **survey weights are never used in
  fitting** (no `weights=` to `mlr3superlearner`), so the ensemble optimizes
  unweighted loss while the targets are design-based prevalences — an internal
  inconsistency that should at least be documented and ideally fixed via
  observation weights.
- *Better:* an outer LOCO/spatial loop for *honest* performance and an inner CV
  loop for weight selection (nested CV). Add monotone/known-direction priors
  only if defensible.

### B2. Cross-validation honesty
- Cluster-blocking is correct but **insufficient for small-area estimation**:
  spatially adjacent clusters in *different* CV folds still leak spatial signal.
  *Better:* **spatial block CV** (contiguous Admin-1 or k-means-on-centroids
  blocks) as a second, more conservative estimate; report both.
- The within-country numbers should be regenerated after P1; expect AUC/Brier-
  skill to fall and to move toward the LOCO numbers.

### B3. Small-area estimation & benchmarks
- **Fay-Herriot** *is* given proper per-area sampling variances
  (`R/benchmark_models.R:251-254`, `p(1-p)/n_eff` with `deff=1.5`) — good — but
  `deff` is a constant guess and FH is hard-capped to 5 covariates
  (`:241-246`); other methods use all screened vars, so the comparison is not
  fully like-for-like. Use survey-estimated design effects; document the cap.
- **BYM2/INLA** PC-priors are fixed and not sensitivity-tested
  (`:405-407`); isolated/island areas are zeroed with a comment that "INLA
  handles it" but no warning (`:383-391`); adjacency build
  (`build_adjacency_list`, `:2029-2076`) has fragile `match()` indexing.
- **Baseline choice.** The "baseline" is the training-pooled mean — fine for
  transport. But the **right status-quo comparator for a ministry** is *"use the
  national survey average for every district"* (current practice) and *"use the
  direct survey estimate where it exists."* Add those explicitly so the benefit
  claim is against real practice, not a strawman.
- **Area-level vs unit-level** are separate pipelines that never meet in one
  LOCO frame (`R/area_level_comparison.R` vs the individual SL aggregated to
  Admin-2). They cannot be cleanly compared today. Put both through one
  evaluation harness on identical hold-outs.
- **Metrics vs noisy truth** (B-wide): `compute_area_metrics`
  (`R/benchmark_models.R:42-82`) and `compute_admin2_error`
  (`R/admin2_analysis.R:728-749`) treat `svy_prev` as fixed. *Better:* report a
  noise-adjusted error `E[(obs−pred)²] − E[sampling var]`, or inverse-variance
  weight by `n_svy`, or stratify metrics by `n_svy` quintile.

### B4. Conformal / prediction intervals
- **No held-out calibration split** (`R/conformal.R:47-79`) — quantile taken on
  the same CV residuals → coverage not guaranteed for small folds.
- **Admin-1 → Admin-2 broadcast** (`fct_helpers.R`, prep §1) has **no
  district-level coverage guarantee**; conformal validity is per-prediction,
  not preserved under aggregation/broadcast.
- **Design ignored** in the delta-method SE (`R/conformal.R:150-157`): i.i.d.
  weights, no PSU/strata → too narrow. *Better:* compute Admin-1/2 point + SE
  via `srvyr::survey_mean()` (design-based) and layer the conformal width on the
  design SE; validate empirical coverage on the LOCO hold-outs.

### B5. Calibration (Platt)
- **In-sample fit+evaluate** (`R/benchmark_models.R:~1957`; `R/diagnostics.R:100-112`).
  Switch to the existing fold-aware `calibrate_binary_oof`, or fit Platt on a
  calibration fold and evaluate on a disjoint fold. Re-issue
  `diagnostics_binary_calibrated.csv`. Until then, label the "recalibrated"
  metrics as optimistic in the dashboard (the new diagnostics tab should carry
  this caveat).

### B6. Transportability / LOCO design
- LOCO is **methodologically honest** (preprocessing fit on training countries
  only) — the design is sound. But **n=4 held-out countries** cannot support
  statistical inference: report it as **hypothesis-generating**, give per-country
  results (already done), and avoid pooled summaries that hide the
  Malawi/Sierra-Leone failures.
- **CIV external validation** has a documented **feature-scale mismatch**
  (commit `bb221bb`; `scripts/civ_external_validation_pipeline.R`): pre-downloaded
  CIV rasters use different units/aggregation than pipeline-extracted training
  features (e.g. soil-calcium training mean 0.30 vs CIV 87). Predictions there
  are **low-confidence**; keep the honest framing and add an out-of-support flag
  rather than a number on a map.

---

## C. Model-improvement recommendations (prioritized)

| Rec | What & why | Impact | Effort |
|-----|-----------|--------|--------|
| C1 | **Honest pipeline + spatial CV** (P1, B2). Fold-internal preprocessing as an mlr3 graph; add spatial-block CV. Truthful numbers are prerequisite to every claim. | ★★★★★ | M |
| C2 | **Partial-pooling area model as primary** (P7). A hierarchical/Bayesian area model (e.g. BYM2 or a multilevel logistic with random Admin-1 effects + covariates) borrows strength to sparse districts and yields coherent intervals — better than the current "pick a winner among MSE/logit losses." | ★★★★ | L |
| C3 | **Design-aware fitting & uncertainty** (P4/P5). Observation weights in SL; `srvyr` design-based direct estimates as the SAE anchor; sampling-error-aware metrics. | ★★★★ | M |
| C4 | **Out-of-support / transport-risk features** (P6). Per-district covariate-support distance (Mahalanobis/`isolation` to training cloud), interval width, and LOCO calibration → a single "reliability" score. Improves *decision usefulness* even if RMSE is flat. | ★★★★ | M |
| C5 | **Spatial structure in features** (cheap win). Admin-2 centroids already added for prescreened SL (`add_admin2_centroids`); add spatial smooths / neighbor-mean covariates and a spatial random effect to the area model. | ★★★ | M |
| C6 | **Better feature curation** for transport: keep only covariates with consistent cross-country definitions and stable scales (the CIV mismatch shows scale drift breaks transport). Pre-register the transport feature set. | ★★★ | S–M |
| C7 | **Right-size the ensemble per outcome/N.** Rare outcomes (<5–10%) are unstable; consider firth/penalized logistic + isotonic calibration and explicit "insufficient signal" reporting instead of forcing the full stack. | ★★ | S |

> Note on payoff: for *decision usefulness*, C3/C4 likely beat chasing a higher
> AUC. The targeting analysis (Part D) shows accuracy is already good enough to
> help in some outcome×country cells and useless in others — knowing *which* is
> worth more than a small average AUC gain.

---

## D. Stakeholder translation (the core ask)

A ministry planner does not act on AUC. They act on questions like *"Which
districts do I fund first? How much better is that than what I do now? Where
should I not trust this?"* These are answerable from existing results.

### D1. Targeting accuracy / decision value — **prototyped**
*If you rank districts by the model and act worst-first, what share of the truly
deficient population do you reach vs having no sub-national data (where you can
only reach people in proportion to where they live)?* This is a burden-
concentration ("enrichment") curve; the gap above the 45° line is the value the
model adds. Headline numbers (real, from the **out-of-sample LOCO** predictions
in `transport_calibration.rds$adm2_area`, population-weighted):

| Outcome | Scope | n | transport r | reach @ worst-20% pop | lift vs no-targeting |
|---------|-------|---|------|------|------|
| child_iron | Ghana | 75 | 0.53 | 38% | **1.90×** |
| women_iron | Ghana | 75 | 0.35 | 31% | 1.56× |
| women_iron | POOLED | 206 | 0.34 | 28% | 1.42× |
| child_iron | POOLED | 206 | 0.37 | 26% | 1.32× |
| child_vitA | Ghana | 75 | 0.38 | 36% | 1.82× |
| women_vitA | Sierra Leone | 14 | −0.14 | 6% | **0.32×** |
| women_vitA | Malawi | 87 | −0.12 | 3% | **0.15×** |

The message a decision-maker can read directly: **for iron in Ghana the model
roughly doubles targeting efficiency; for vitamin-A in Malawi it is worse than
acting blindly.** This is exactly what AUC hides.

### D2. Status-quo vs model — **prototyped**
Per outcome × country: targeting **lift** (prevalence in the model's worst
quintile ÷ overall) and **top-quintile hit-rate** (share of the truly-worst
quintile the model flags). Green/amber/red so a non-technical reader sees at a
glance where model-guided prioritisation beats "use the national average."

### D3. Fitness-for-purpose ("when NOT to trust it") — **prototyped**
A traffic-light scorecard by outcome × country, from the **out-of-sample**
transported agreement (rank correlation + MAE): *Good for targeting /
Use with caution / Not reliable for targeting.* Tied directly to the LOCO
diagnostics, with the optimistic within-country correlation shown alongside for
contrast. This operationalizes "where the model should not be trusted."

### D4. Calibration in plain language — **prototyped**
Instead of a calibration slope, sentences: *"When the model predicts around
30%, the survey-observed rate averages 27% (about as predicted; n = 1,200)."*
Generated per country × outcome from the binned diagnostics.

### D5. Further ideas (not yet built — recommended next)
- **Decision-curve / net-benefit** at an explicit action threshold (treat a
  district if predicted prevalence ≥ t), showing net benefit vs treat-all /
  treat-none — the clinical-decision-analysis analogue of D1, once intervals are
  honest (P4).
- **Trust map:** shade the choropleth by reliability (out-of-support distance +
  interval width + LOCO verdict), or hatch "low-confidence" districts so the eye
  is not drawn to precise-looking numbers in untrustworthy areas (needs C4).
- **People-level framing everywhere:** "≈X,000 deficient children you would
  reach" rather than percentages (population join already exists).
- **Survey-vs-model toggle** on the map: where a survey exists, show it; where it
  doesn't, show the model *with* its reliability flag — making the "if/when"
  explicit spatially.

### Prototypes delivered (this review)
- **New dashboard tab "Decision value"** — `dashboard/R/mod_decision_value.R`,
  wired in `dashboard/app.R` (after *National burden*). Three panels:
  *Targeting accuracy* (enrichment curve + headline value boxes, realistic/
  optimistic toggle), *Status quo vs model* (lift/hit-rate scorecard), and
  *Fitness for purpose* (traffic-light verdicts + plain-language calibration).
  All computed from existing `transport_calibration.rds`, `admin2_predictions`,
  `admin2_population`, `benchmarks.rds`, `model_diagnostics.rds` — no retraining.
- (Context) The earlier session also added **Model diagnostics** and
  **Benchmarks** tabs surfacing ROC/PR/calibration and SAE comparisons that
  were previously computed but unshown.

---

## Verification
- `model_diagnostics.rds` / `benchmarks.rds` bundles build from existing CSVs
  (`dashboard/data-raw/_build_new_bundles.R`).
- Headless construction: `global.R` sources; every `mod_*_ui()` builds; the full
  `shiny.appobj` constructs (all tabs incl. *Decision value*).
- `shiny::testServer` smoke test on `mod_decision_value_server`: all panels
  render under both the realistic (LOCO) and within-country settings, pooled and
  single-country, with no runtime errors (one NA-bin bug found and fixed in the
  calibration readout).
- Targeting headline numbers reproduced independently
  (`dashboard/data-raw/_calc_headline.R`).

## How to regenerate
```
# refresh dashboard bundles (light; no GADM download)
Rscript --vanilla dashboard/data-raw/_build_new_bundles.R
# full dashboard data refresh after a pipeline rerun
Rscript --vanilla dashboard/data-raw/01_prepare_dashboard_data.R
# launch dashboard
R -e "shiny::runApp('dashboard')"
```

## Top 5 actions, in order
1. **P1** — fold-internal preprocessing/prescreen; regenerate honest
   within-country `cv_perf`/`diagnostics_*`.
2. **P3/D** — lead with decision-value metrics (now in the dashboard); reframe
   the manuscript's results around targeting value and fitness-for-purpose.
3. **P2** — out-of-fold calibration; re-issue calibrated diagnostics.
4. **P4/P5** — design-based, design-aware uncertainty and sampling-error-aware
   error metrics.
5. **P6/C4** — out-of-support & trust flags on the map (the visual "where not to
   trust it").
