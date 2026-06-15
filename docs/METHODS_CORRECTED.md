# Corrected methods: eight methodological additions (P1–P8)

**Status:** parallel pipeline built, validated, and **run at full scale**
(24 slices / 293 targets; `tar_make` completed 2026-06-14 23:12, 0 errored
targets). Full-run numbers are filled in below. **Note:** corrected fits
produced metrics for **18 of 24** slices; the **6 Sierra Leone** slices returned
no corrected fit (only the production reference is available for them) — a
known gap, debuggable by resuming the cached run (see limitations). One-slice
smoke values are retained and labelled *(smoke: Gambia × women_iron)*.

**Cross-references.** Problems and recommendations: [`CRITICAL_REVIEW.md`](CRITICAL_REVIEW.md)
(P1–P8). External justification: [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md).
Implementation/run: [`CORRECTED_PIPELINE_RUN.md`](CORRECTED_PIPELINE_RUN.md);
preliminary results: [`CORRECTED_PIPELINE_RESULTS.md`](CORRECTED_PIPELINE_RESULTS.md).

## Overview

The production pipeline (`_targets.R`) is left untouched. A parallel pipeline
(`_targets_corrected.R`, store `_targets_corrected/`, helpers in `R_corrected/`)
re-derives the analysis with eight methodological corrections, reusing
production data-loading outputs (`outcome_data_*`, `svy_admin2_*`,
`gee_admin2_*`, read from `_targets_full/objects/` via `read_prod()`) so that
only the *methods under study* differ. Each correction emits a target whose
output is placed side-by-side with the production analogue
(`methods_comparison` → `dashboard/data/methods_comparison.rds`,
`results/tables/corrected/*.csv`). Comparisons are *identity-controlled*: the
honest and optimistic estimators use the **same data and the same learner
library**, so any difference isolates the method, not the configuration.

---

## P1 — Leakage-free fitting and deployment-honest cross-validation

**Problem (production).** In `R/mlr3_fitting.R` the entire preprocessing chain is
fit on the full sample *before* resampling: imputation and near-zero-variance
pruning (`:282–293`), **supervised univariate prescreening against the full
outcome** (`washb_prescreen(Y = Y, Ws = cov, …)`, `:307`), correlation filtering
and normalization (`recipes::…|> prep()`, `:341`), with cross-validation folds
created only afterwards (`origami::make_folds`, `:409`) and handed to
`mlr3superlearner`. Because feature selection and scaling have already seen the
held-out rows (the prescreen uses `Y` directly), the cross-validated metrics
(`cv_perf`, `diagnostics_*`) are optimistically biased. Folds are also
cluster-blocked but not spatially blocked, so spatially adjacent clusters in
different folds still leak spatial signal.

**Correction (`R_corrected/p1_fitting.R`, `00_corrected_utils.R`).** A
fold-internal pipeline: for every outer fold, `fit_preproc()` is fit on the
*training rows only* — drop constant/all-NA columns → median imputation
(training medians) → variance filter → **supervised prescreen** (top-K by
|point-biserial correlation with the training Y|) → greedy correlation filter →
standardize (training mean/sd) — and then `apply_preproc()` transforms the
held-out fold with the training-derived parameters. A compact **discrete
SuperLearner** selects, by inner CV on the training fold, the lowest-loss base
learner (library mean + glmnet + ranger [+ rpart at full scope]), refits it on
the full training fold, and predicts the held-out fold, yielding genuinely
out-of-fold predictions. The fit is run under three schemes: **cluster-block CV**
(folds = survey PSU groups), **region/spatial-block CV** (folds = Admin-1 region
groups, `make_group_folds`), and an **optimistic** mirror of the production
leakage (preprocessing fit on the full data + naive random CV,
`run_oof_optimistic`).

**Why it is more honest.** Validation error must be estimated with the data's
structure respected and with no information flowing from test to train;
otherwise predictive error is "seriously underestimated" (Roberts et al. 2017,
Ecography) and apparently-skilful models can have near-null real skill (Ploton
et al. 2020, Nat. Commun.). Fold-internal preprocessing is the standard remedy
for the exact leakage class we have (John, Saurette & Heung 2025, Geoderma);
spatial blocking aligns evaluation with deployment to new areas (Brenning &
Suesse 2026; Linnenbrink, Nowosad & Meyer 2026). See
[`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 1.

**Evaluation.** Per slice: ΔAUC and ΔBrier-skill of *optimistic − honest*, and
*cluster − spatial*; the production `cv_perf` is shown for reference. Honest
schemes should sit below the optimistic/production estimates; spatial-block is
the most conservative.

**Results.** *(smoke: Gambia × women_iron, same library/data):* optimistic AUC
0.592 vs honest cluster 0.579 vs honest spatial-block 0.579 (≈ production 0.579);
Brier-skill 0.026 optimistic vs 0.004 spatial-block. **Full run (18 corrected
slices):** mean OOF ROC-AUC optimistic **0.571**, cluster-block **0.547**,
spatial-block **0.518**; mean leakage inflation **optimistic − spatial-block =
+0.053 AUC** (optimistic exceeds spatial-block in 15/18 slices), spatial-block
the most conservative throughout — confirming that the production within-country
numbers are optimistic and that honest spatial CV is materially lower.

---

## P2 — Out-of-fold (honest) probability calibration

**Problem (production).** The Platt recalibrator is fit and evaluated on the
*same* out-of-fold predictions (`R/benchmark_models.R::fit_binary_recalibrator`,
called in `run_diagnostics_calibrated` ~`:1957`; logit-scale slope/intercept in
`R/diagnostics.R:100–112`). In-sample calibration is optimistic by construction:
the reported post-calibration slope ≈ 1 and the Brier-skill "recovery" in
`diagnostics_binary_calibrated.csv` overstate real calibration. A fold-aware
routine (`calibrate_binary_oof`) exists but is not the one used for the headline.

**Correction (`R_corrected/p2_p6_methods.R::diagnostics_calibrated_oof`).** A
fold-aware Platt recalibrator: for each fold, fit the logistic recalibrator on
the *out-of-fold* predictions and apply it *in-fold* (mirroring
`calibrate_binary_oof`'s logic). We report the calibration slope/intercept and
Brier for three states — raw, in-sample Platt (the production anti-pattern), and
out-of-fold Platt (honest) — so the optimism is visible. (Isotonic can be
swapped in identically; Platt is the default for the low-prevalence binaries.)

**Why it is more honest.** Calibration is a decision-relevant performance domain
and must be assessed out-of-sample, not on the data used to fit the recalibrator
(Lancet Digital Health 2025 guidance; cf. the under-reporting documented by
Ferreira et al. 2020). See [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 3–4.

**Evaluation.** In-sample vs out-of-fold calibration slope/intercept and Brier
per slice; the divergence is the optimism removed.

**Results.** *(smoke):* in-sample Platt slope **1.00**, intercept 0.00; honest
out-of-fold Platt slope **−0.05**, intercept −0.83 — i.e. essentially no
calibratable signal once assessed honestly. **Full run:** the in-sample Platt
slope is **exactly 1.00 in all 18 slices** (the artefact); the honest
out-of-fold slope averages **−0.17** (range −4.39 to 0.83) — the in-sample
calibration illusion is universal across slices.

---

## P3 — Decision-value metrics as first-class targets

**Problem (production).** Model usefulness was reported only as AUC/MSE, which do
not tell a decision-maker whether acting on the model beats the status quo; the
dashboard's decision view was computed ad hoc rather than from pipeline targets.

**Correction (`R_corrected/p2_p6_methods.R::decision_value_corrected`).** From the
corrected Admin-2 predictions (honest OOF aggregated to district means, joined to
survey truth on `(country, Admin2)`) we compute, as a pipeline target: the
**burden-concentration / enrichment** reach at the worst-20% of population
(`reach_at_20pct`), the **targeting lift** versus no sub-national information
(prevalence in the model's worst quintile ÷ overall), and the **top-quintile
hit-rate**. Districts are weighted by survey size (`n_svy`) as an in-pipeline
population proxy; the dashboard uses true subgroup population.

**Why it matters.** Clinical/operational utility at an action threshold — net
benefit relative to treat-all/treat-none — is the decision-relevant evaluation
domain (Vickers & Elkin 2006; Singh, Shah & Vickers 2023 for the
resource-constrained "fund the top-N" case; Lancet Digital Health 2025). See
[`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 4. (A full
decision-curve/net-benefit target is the planned extension once P4 intervals are
in place.)

**Evaluation.** reach@k, lift and hit-rate per slice, contrasted with the
no-targeting baseline (lift = 1).

**Results.** *(smoke):* reach@20% 0.23, lift 1.15×, hit-rate 0.17. Full
**Full run:** mean reach@20% **0.26**, mean lift **1.29×**, mean top-quintile
hit-rate **0.31**; lift > 1.2× in **13/18** slices — model-guided targeting beats
the national-average status quo in most outcome × country cells, but not all.

---

## P4 — Honest prediction intervals: split-conformal + survey-design CIs

**Problem (production).** Conformal intervals use cross-validation residuals with
no held-out calibration split (`R/conformal.R:47–79`), are computed at Admin-1
and **broadcast to Admin-2** (no district-level coverage guarantee;
`dashboard/R/fct_helpers.R`, prep §1), and the delta-method aggregation treats
observations as i.i.d., ignoring strata/PSU (`R/conformal.R:150–157`) — so
intervals are likely too narrow, especially for sparse districts.

**Correction (`R_corrected/p2_p6_methods.R::intervals_corrected`).** Two
district-level interval sources, no broadcast: (i) **split-conformal** — reserve
a random half of districts as a *calibration split* whose absolute residuals set
the conformal half-width `q̂ = Quantile_{(1−α)(n_cal+1)/n_cal}(|pred−direct|)`,
applied as `pred ± q̂` to the remaining districts, with **empirical coverage
reported on the held-out districts**; (ii) **design-based CIs** taken directly
from the survey (`srvyr`-derived `svy_prev_low/upp` per district), with an
explicit `design_available` flag where a district has no survey estimate.

**Why it is more honest.** Split-conformal restores the exchangeability/held-out
property that gives conformal its coverage guarantee; district uncertainty for
survey quantities is properly design-based (Mercer et al. 2015; Wakefield 2020).
See [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 2–3.

**Evaluation.** Empirical held-out coverage vs the nominal level; interval width;
count of districts with design-based CIs (vs broadcast).

**Results.** *(smoke):* split-conformal half-width 37.3 pp, empirical held-out
coverage 0.87 (target 0.90), 30/30 districts with design CIs. Full cross-slice
**Full run:** mean empirical held-out coverage **0.92** (target 0.90; slightly
conservative), mean half-width 29.5 pp, **1,266** district-level design-based CIs
(no Admin-1 → Admin-2 broadcast).

---

## P5 — Sampling-error-aware accuracy metrics

**Problem (production).** Error metrics compare predictions to the direct survey
estimate as if it were exact truth (`R/admin2_analysis.R::compute_admin2_error`
`:728–749`; `R/benchmark_models.R::compute_area_metrics` `:42–82`). The survey
estimate carries its own sampling variance (~`p(1−p)/n_eff`, ≈4–10 pp for n≈50),
which inflates apparent error and biases correlations downward — sparse districts
look "wrong" when they are merely noisy.

**Correction (`R_corrected/p2_p6_methods.R::admin2_error_corrected`).** We report
the naive RMSE/MAE *and* a **sampling-variance-adjusted** RMSE that subtracts the
mean sampling variance of the truth from the mean squared error,
`adj_MSE = max(0, E[(pred−direct)²] − E[Var_sampling(direct)])`, using the
reported `svy_prev_se²` (or `p(1−p)/n_svy` when absent). Pearson r is reported
alongside.

**Why it is more honest.** When the reference is itself an estimate, regressing
out its sampling variance recovers the model's true error (errors-in-variables
logic; the design-based SAE tradition makes the design variance explicit —
Mercer et al. 2015; Wakefield 2020; the micronutrient SAE work of arXiv
2604.14971). See [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 2.

**Evaluation.** naive vs adjusted RMSE (and the sampling RMSE removed) per slice.

**Results.** *(smoke):* naive RMSE 22.4 pp, sampling RMSE 4.6 pp, adjusted
21.9 pp; r 0.63. **Full run:** mean naive RMSE **15.6 pp** vs sampling-adjusted
**15.4 pp** (mean survey sampling RMSE 2.5 pp removed); mean r 0.23.

---

## P6 — Out-of-support / transport-based trust flags

**Problem (production).** No per-area signal of where predictions should *not* be
trusted: an unsurveyed district far outside the surveyed covariate support, or
one with a very wide interval, is shown with the same visual authority as a
well-supported one.

**Correction (`R_corrected/p2_p6_methods.R::trust_flags_corrected`).** For every
Admin-2 polygon (surveyed and not), using the GEE covariate table: standardize a
capped, well-behaved covariate subset, compute each district's **out-of-support
distance** as the mean standardized distance to the centroid of the surveyed
(training) districts, and flag `out_of_support` beyond the 95th percentile of the
training distances. Combine with **small sample** (`n_svy < 15`) and **wide
interval** (P4 width > 1.5× median) into a traffic-light verdict —
*OK to use / Use with caution / Do not rely*.

**Why it matters.** Whether a prediction is an interpolation or an extrapolation
relative to the training support determines which evaluation regime applies and
whether to trust it (Brenning & Suesse 2026; Linnenbrink et al. 2026); flagging
unreliable areas operationalizes the value-of-information caution against acting
on uncertain estimates (Griffin, Claxton & Sculpher 2008). See
[`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Themes 1 & 4.

**Evaluation.** Distribution of districts across trust tiers; concordance between
"Do not rely"/out-of-support flags and large prediction errors where survey
truth exists.

**Results.** *(smoke):* 32 OK / 3 caution / 2 do-not-rely (of 37 polygons). Full
**Full run (all polygons × outcomes):** **2,080 OK to use / 1,294 use with
caution / 278 do not rely** (out-of-support).

---

## P7 — Design-aware partial-pooling area model as a primary estimator

**Problem (production).** Small-area-estimation methods are used only as
*comparators* in a benchmark, not as a primary estimator, and the area machinery
has issues: Fay-Herriot is hard-capped to 5 covariates (`R/benchmark_models.R`
`:241–246`) while other methods use more; BYM2/INLA priors are fixed and
untested, and island/disconnected areas are zeroed without warning (`:383–391`).

**Correction (`R_corrected/p7_area.R::area_partial_pooling_corrected`).** A
design-aware partial-pooling estimator promoted to first-class: a **Fay-Herriot
area-level model** (reusing `R/benchmark_models.R::fit_predict_fh`, which supplies
proper per-area sampling variances `p(1−p)/n_eff`) produces a smoothed
(shrunken) district prevalence, placed beside the raw direct survey estimate and
the corrected ML estimate. When FH is numerically singular (collinear area
covariates on few areas) the target degrades gracefully to **empirical-Bayes
shrinkage** toward the design-weighted global mean with weight `B/(B+v_i)`
(between-area variance `B`, area sampling variance `v_i`) — still design-aware.

**Why it is more correct.** A design-aware partial-pooling cluster/area model
should be the primary estimator (it borrows strength to sparse areas and yields
coherent uncertainty), not merely a benchmark; the directly analogous
micronutrient SAE study finds a cluster-level Beta-binomial best and an
area-level joint-smoothing model best among design-aware options (arXiv
2604.14971; Wakefield 2020; Mercer et al. 2015; the principled workflow of
Dong & Wakefield 2025). See [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 2.

**Evaluation.** Partial-pooled vs direct vs ML estimates; shrinkage magnitude;
out-of-sample MSE to held-out survey truth. *(Planned extension: a cluster-level
Beta-binomial and a properly-conditioned BYM2.)*

**Results.** *(smoke):* FH was computationally singular on 30 collinear-covariate
areas → empirical-Bayes shrinkage fallback engaged and produced smoothed
estimates. **Full run:** Fay-Herriot (REML) converged for **486** area-rows; the
empirical-Bayes shrinkage fallback covered **780** (where FH was singular) — both
design-aware. BYM2 / cluster-level Beta-binomial remain planned upgrades.

---

## P8 — Latent-bug fixes (centering leak, OOS imputation, key joins)

**Problem (production).** Three latent defects flagged in the review: (i) the
two-stage transport model centers the held-out country by its *own* mean
(`R/benchmark_models.R:577–579`) — a forward look at the test level; (ii) the OOS
predictor matrix pads missing covariates with 0 and can fall back to imputing
with the *held-out* country's own medians (`R/oos_prediction.R:175–189`); and
(iii) survey/area joins keyed on `Admin2` alone collide across countries — which
produced a cartesian blow-up (an actual segfault during the review).

**Correction (enforced by construction in `R_corrected/`).** (i) All
standardization uses **training statistics only** (`fit_preproc`/`apply_preproc`
in `00_corrected_utils.R`; no held-out centering anywhere). (ii) Out-of-sample
imputation uses **training medians** only — the corrected path never imputes a
target row from its own values, and never zero-pads silently. (iii) **Every**
survey/area join is keyed on `(country, Admin2)` (`corrected_admin2`,
`build_methods_comparison`, P7 merges), eliminating the cross-country collision.

**Why it is correct.** These remove forward-looking leakage and a data-integrity
fault; the leakage definition and its consequences are documented in John,
Saurette & Heung (2025) and the spatial-CV literature
([`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) Theme 1).

**Evaluation.** Regression-style: the corrected pipeline builds all 24 slices /
293 targets without the cartesian blow-up; spot-check that transported/OOS
estimates do not shift when held-out rows are permuted (invariance to test-set
content). Smoke: full graph built, no crash. **Full run integrity:** all 24
slices built, **0 errored targets** (`tar_errored()` empty), no cartesian
blow-up.

---

## Head-to-head evaluation plan and results skeleton

The `methods_comparison` target assembles all corrected outputs with the
production analogues into `dashboard/data/methods_comparison.rds` and
`results/tables/corrected/*.csv`; the dashboard **"Methods (corrected)"** tab
renders them. Full-run results (`CORRECTED_SCOPE=full`, 24 slices; 18 with
corrected metrics):

| Comparison | Metric | Production / optimistic | Corrected (honest) | Result |
|------------|--------|------------|--------------------|----|
| P1 CV honesty | mean OOF ROC-AUC | optimistic **0.571** | cluster **0.547** / spatial **0.518** | leakage inflation **+0.053 AUC** (opt−spatial), 15/18 slices |
| P2 calibration | Platt slope | in-sample **1.00** (all 18) | out-of-fold mean **−0.17** | in-sample calibration is an illusion |
| P3 decision value | lift vs no-targeting | (ad hoc) | mean **1.29×**; **13/18** > 1.2× | targeting beats status quo in most cells |
| P4 intervals | held-out coverage | broadcast Admin-1 (no guarantee) | split-conformal **0.92** @ target 0.90; **1,266** design CIs | coverage restored, no broadcast |
| P5 error | RMSE (pp) | naive **15.6** | sampling-adjusted **15.4** | survey sampling noise (2.5 pp) removed |
| P6 trust | tier counts | (none) | **2,080** OK / **1,294** caution / **278** do-not-rely | unreliable areas flagged |
| P7 area model | estimator role | comparator-only | primary partial pooling: FH **486** + EB **780** rows | design-aware estimator promoted |
| P8 integrity | builds / errors | n/a (cartesian blow-up risk) | 24 slices, **0 errored**, no blow-up | by construction |

Filled from the completed full run (2026-06-14 23:12). The 6 Sierra Leone slices
lack corrected metrics (production reference only) and can be backfilled by
resuming the cached run; numbers above are over the 18 slices with corrected
fits. The dashboard "Methods (corrected)" tab renders the same bundle.
