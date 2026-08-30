# Analysis / modeling code audit — leakage & correctness

READ-ONLY audit (2026-06-23). Scope: in-sample-vs-OOS leakage, preprocessing
leakage, metric/calibration errors, transportability honesty. Effective n = #
Admin-2 areas (14-87), so optimism inflates metrics badly.

Context established up front (NOT bugs):
- `res$yhat_full` IS genuine out-of-fold in both fitting paths
  (`src/analysis/sl_helpers.R:130` uses `sl_fit$predict_fold(SL_task,"validation")`;
  `R/sensitivity/mlr3_fitting.R:704` uses `yhat_oof`, falling back to in-sample
  only when OOF unavailable). So `diagnostics.R`, `admin1_analysis.R::extract_cv_performance`,
  and the binary metrics in `national_estimates.R` read honest OOF predictions.
- `fit_area_sl` (R/area_level_comparison.R) was already fixed: `loo_preds`/metrics
  now come from an outer-CV loop (lines 135-189), not `predict(fit, train)`.
- The two ablation files correctly compare in-sample baseline vs in-sample
  permuted predictions (`yhat_insample`), so the *drop* is a fair contrast.

---

## HIGH severity

### H1. Conformal "coverage" and interval width are computed in-sample on the same residuals used to set the width — not a held-out check
File: `R/conformal.R:72-119` (and the reported `constant_coverage`/`adaptive_coverage` in `diagnostics`, lines 256-265).
- The conformity scores `residuals <- abs(Y - yhat)` (line 72) are OOF residuals
  (good), but the (1-alpha) quantile (line 79) is taken over the SAME residual
  vector that coverage is then evaluated on (lines 89, 117). By construction the
  empirical coverage of the in-sample quantile is ~1-alpha; it is NOT an
  out-of-sample validity check and should not be reported as evidence the
  intervals "cover."
- The locally-adaptive path is worse: `q_adaptive` (line 106) is the quantile of
  normalized residuals, then `adaptive_width` rescales by the per-point
  `sqrt(p(1-p))` / `median(residuals)`. Coverage (line 117) is again measured on
  the fitting data. Reported adaptive coverage is circular.
Why a bug: any coverage number reported from this function is self-fulfilling;
it cannot detect miscalibration. Severity high because CIs feed the dashboard
and `aggregation_uncertainty.R`.
Fix (one-line intent): split the OOF residuals into a calibration half (set the
quantile) and an evaluation half (measure coverage), or at minimum stop
reporting `constant_coverage`/`adaptive_coverage` as a validity check.

### H2. Conformal Admin-1/Admin-2/national interval reverses the conformal guarantee (delta-method shrinks the interval instead of preserving coverage)
File: `R/conformal.R:153-164, 178-182, 225-241`.
- Individual conformal half-width is a distribution-free *prediction* interval.
  The aggregation converts it to a per-point Gaussian SD
  `sigma_i = half_w / qnorm(0.975)` and then forms `SE(mean) = sqrt(sum(w_i^2 sigma_i^2))`
  (lines 153-157). For ~equal weights over `m` districts this divides the width
  by ~sqrt(m), so a 95% individual PI collapses to a tiny mean-CI. That is a
  *confidence interval for the mean prediction*, not a prediction interval, and
  it silently drops the conformal coverage guarantee the function advertises.
- It also ignores correlation between district errors (spatial autocorrelation),
  so even as a mean-CI it is anti-conservative.
Why a bug: the Admin-1/national CIs are mislabeled and far too narrow; the
`boot_mean`/`ci_*` columns are consumed downstream as if bootstrap PI bounds.
Severity high (numbers reported in manuscript/dashboard).
Fix: aggregate by propagating *prediction* intervals (e.g. average the
individual `ci_lo`/`ci_hi`, or bootstrap district draws) rather than
sqrt(n)-shrinking a Gaussianized half-width.

### H3. `compute_continuous_diagnostics` R^2 (and `extract_cv_performance` r2) use a non-CV baseline for the OOF predictions, but more importantly r2 can be reported as a "CV" fit while the threshold-derived prevalence metrics are computed point-wise on OOF predictions of a continuous model
File: `R/diagnostics.R:209` and `R/admin1_analysis.R:42`.
- `r2 <- 1 - sum((Y-yhat)^2)/sum((Y-mean(Y))^2)` with OOF `yhat` is fine as an
  OOS R^2 *definition*, but note `extract_cv_performance` (admin1_analysis.R:42)
  computes it WITHOUT `na.rm`, so a single NA yhat silently yields NA — low risk,
  flag for verification.
- Real issue: `compute_continuous_diagnostics` derives `prev_pred` by
  thresholding individual OOF continuous predictions (lines 224-227). Thresholding
  a regression prediction at the individual level systematically biases prevalence
  (Jensen-type bias: `mean(1[yhat<c])` != `mean(1[Y<c])` even if yhat is unbiased
  for Y). `prev_bias` will look like model bias when it is largely an artifact of
  per-individual thresholding. Severity high-ish because `national_estimates.R:91`
  does the SAME thing for the continuous model's national `pred_prev`.
Fix: predict prevalence as `mean(P(Y<c))` via a calibrated probability, or fit
the binary model for prevalence rather than thresholding point predictions.

---

## MEDIUM severity

### M1. `calibrate_binary_oof` default path fits the Platt calibrator on ALL data then applies it to the same data (in-sample recalibration)
File: `R/benchmark_models.R:2006-2013`.
- When `fold_col` is NULL (the default and the path `run_diagnostics_calibrated`
  effectively relies on — see below) the calibrator is fit on every row and
  applied to those same rows. The comment admits it is "slightly optimistic."
- `run_diagnostics_calibrated` (lines 1957-1992) fits the Platt recalibrator on
  `res$yhat_full` / `res$Y` (all OOF rows) at line 1963 and then evaluates the
  calibrated predictions' Brier/calibration on the SAME rows (line 1968). The
  2-parameter Platt fit is cheap to overfit at n = number of areas; the resulting
  `diagnostics_binary_calibrated.csv` Brier-skill improvement is optimistically
  biased.
Why a bug: the calibrated diagnostics overstate the benefit of recalibration.
Severity medium (it is a recalibration step, small parameter count, but n is tiny).
Fix: route `run_diagnostics_calibrated` through the fold-aware branch of
`calibrate_binary_oof` (pass the CV fold id) so calibration is out-of-fold.

### M2. `fit_predict_glm` / `fh` / `mixed` / `gam` / `quasibinomial` / `betareg` / `dag` do correlation-based variable PRE-SELECTION using the training labels — honest for LOCO, but `train_pred` (used as the "in-sample" reference everywhere) screens on the full training outcome
File: `R/benchmark_models.R` (e.g. 168-175, 242-247, 471-476, 716-721, 927-932, 1007-1012, 1729-1735).
- For the LOCO test predictions this is fine: screening uses TRAIN-only
  `train$svy_prev`, never the held-out country. Confirmed `run_area_benchmarks_loco`
  (line 2141) splits train/test by country before `.screen_vars` (line 2144) and
  before each method's internal screen. No held-out leakage. GOOD.
- Caveat (not leakage, but optimism): `within`-country k-fold runner
  (`run_area_benchmarks_within`, lines 2238-2273) calls the SAME methods, which
  re-screen by univariate correlation *inside each fold's training set* — correct
  — BUT `.screen_vars` at line 2242 and the top-K correlation caps reduce p using
  only train rows, so within-country CV is honest. Verify only that no method uses
  `test$svy_prev`. None do. Flagging as VERIFIED-CLEAN, listed for completeness.

### M3. `compute_area_metrics` calibration slope/intercept are fit on the held-out points but reported without guarding tiny-n instability; `calib_intercept` is on the raw prevalence scale while several methods report it as if comparable to the logit-scale calibration elsewhere
File: `R/benchmark_models.R:64-67`.
- `lm(obs ~ pred)` on as few as 3 held-out areas gives a slope with enormous
  variance; reported to 3 decimals it looks precise. Not wrong, but misleading at
  n_test = 3-14. Severity medium for interpretation.
- Note the transportability LOCO (`transportability.R:316-320`) computes
  calibration on the logit scale (`glm(Y ~ qlogis(pred))`) while
  `compute_area_metrics` uses the identity scale — the two `calib_slope` columns
  are NOT comparable but share a name across output tables. Flag for the
  manuscript so they are not pooled.
Fix: add `n_eval` gating / CI to calibration slope, and disambiguate the column
name (e.g. `calib_slope_linear` vs `calib_slope_logit`).

### M4. `transportability.R` LOCO calibration intercept/slope swap risk + `null_brier` reference is the prevalence variance, not the prevalence Brier
File: `R/transportability.R:316-320, 349`.
- `glm(Y_test ~ qlogis(preds_clamp))`: `calib[1]` is intercept, `calib[2]` slope —
  used correctly (intercept at line 352, slope at 353). OK.
- `null_brier = mean(Y_test) * (1 - mean(Y_test))` (line 349) is the variance of a
  Bernoulli at the test prevalence, i.e. the Brier of the *test-prevalence*
  constant predictor. But the model is trained on OTHER countries; the honest
  "no-skill" reference for transportability is predicting the *training* mean (or
  a fixed 0.5), not the held-out country's own prevalence (which the model never
  sees). Reporting Brier skill against the held-out prevalence flatters nothing
  but is conceptually the wrong null for a transport claim — it gives the model
  credit for cross-country level shift it cannot know. Severity medium.
Fix: also report Brier vs the training-prevalence constant so the skill score
reflects transportable signal, not in-country base-rate luck.

### M5. `national_estimates.R` `obs_n_events` formula is a no-op / wrong
File: `R/national_estimates.R:44`.
- `obs_n_events <- sum(y_bin*wts)/sum(wts) * sum(valid)` = `weighted_prev * n` =
  expected events under the weighted prevalence, NOT observed event count
  (`sum(y_bin)`). If this column is read as a raw count it is wrong (it is a
  weighted-prevalence-implied count). Severity medium / low depending on use.
Fix: `obs_n_events <- sum(y_bin[valid])` for raw count, or rename to
`obs_expected_events`.

---

## LOW severity / verify

### L1. `bootstrap_loco_ci` resamples held-out areas AFTER fitting (no refit) — CI reflects only sampling of the eval set, not model uncertainty
File: `R/benchmark_models.R:1864-1876`.
- The model is fit once on train; only the held-out `(obs,pred)` pairs are
  resampled to bootstrap the Pearson r CI. This is a legitimate CI for the r
  estimate given a fixed model, but it does NOT include training/model
  variability, so it understates total uncertainty. Documented behavior, but the
  manuscript should state the CI is conditional on the fitted model. Low.

### L2. `bootstrap.R` (individual-level) refits SL on each cluster-bootstrap and predicts on ORIGINAL data — this is the correct optimism-aware bootstrap, but predictions for the CI come from a model that has seen (resampled copies of) the prediction rows
File: `R/bootstrap.R:39-117`.
- Bootstrap resamples clusters with replacement, refits, then predicts on the
  full original `d` (line 67-111). Because resampled clusters overlap the
  prediction set, the per-replicate predictions are partly in-sample. For a
  *prevalence aggregate* CI this is the standard cluster bootstrap and is
  acceptable, but the resulting CI is a stability/aggregation CI, not an
  out-of-sample predictive interval. Low; document.

### L3. `aggregation_uncertainty.R` treats district CIs as independent Gaussians
File: `R/aggregation_uncertainty.R:38-48`.
- `sigma_i = half_w/qnorm(0.975)` then independent `rnorm` per district. Ignores
  (a) the conformal width's non-Gaussian tails and (b) spatial correlation of
  prediction errors. Combined with H2 (the input `ci_lo/ci_hi` are already too
  narrow), the propagated national interval is doubly anti-conservative. Severity
  low on its own, but compounds H2 — fix H2 first.

### L4. `single_var_ablation.R` / `domain_ablation.R` importance is in-sample
File: `R/single_var_ablation.R:43`, `R/domain_ablation.R:202`.
- Both deliberately use `yhat_insample` for baseline AND permuted predictions, so
  the *difference* (AUC drop) is a fair in-sample contrast. This is the standard
  permutation-importance recipe and is internally consistent. Not a leakage bug.
  Only caveat: importances are in-sample magnitudes, so absolute drops are
  optimistic; rankings are the intended output. VERIFIED-CLEAN, noted.

### L5. `transportability_area.R` centered-test path uses held-out country's own covariate means
File: `R/transportability_area.R:211-214`.
- When `recipe$center = TRUE`, test covariates are centered by the held-out
  country's own column means and predictions get `+ level` (training mean logit).
  This uses only the held-out COVARIATES (label-honest) but injects the held-out
  country's covariate centroid — a transductive assumption. The default recipe has
  `center = FALSE`, so the shipped path is unaffected. Verify no production target
  flips `center = TRUE`. Low.

---

## Transportability honesty — summary
- `run_loco_cv` (transportability.R), `run_area_loco` (area_level_comparison.R),
  `run_area_benchmarks_loco` (benchmark_models.R), and `run_area_transport_loco`
  (transportability_area.R) ALL split train/test by country BEFORE screening,
  imputation, standardization, and fitting. Imputation medians (`.impute_train_test`,
  benchmark_models.R:95-106; `.tr_prep_X`, transportability_area.R:143-155) are
  TRAIN-ONLY and applied to test. CORAL (method 22) uses only test covariates for
  covariance/mean alignment, never the held-out label. No held-out-label leakage
  found in any LOCO path. This is the good news.
- The remaining issues are (1) circular conformal coverage (H1), (2) mislabeled /
  too-narrow aggregated conformal CIs (H2), (3) individual-threshold prevalence
  bias (H3), (4) in-sample recalibration optimism (M1), and (5) wrong/ambiguous
  null references and counts (M4, M5).
