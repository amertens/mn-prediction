# Corrected-pipeline patch plan (honest-evaluation audit)

Concrete, developer-executable fixes for seven confirmed methodological /
correctness issues in the `R/corrected/` honest-evaluation pipeline and its
upstream comparators (`R/transportability_area.R`, `R/benchmark_models.R`,
`R/conformal.R`). All findings verified against current `main` (line numbers
below reflect the present files; re-check before editing — they drift).

This is a **patch plan only** — no code is changed here. Each item gives the
file/function/line range, the one-sentence defect, the concrete change, the
downstream targets/CSVs/manuscript tables to regenerate, a verification check,
and an effort/risk rating. A recommended execution order and a list of
manuscript numbers likely to move follow at the end.

Key shared facts used throughout:

- The corrected pipeline is wired in [`_targets.R`](../_targets.R) lines
  ~1400–1492 (per-slice loop building `corrected_sl_*`, `corrected_err_*`,
  `corrected_dec_*`, `corrected_int_*`, `corrected_area_*`, `corrected_recon_*`,
  then `corrected_methods_comparison` and `protocol_reconciliation`).
- CSVs land in `results/tables/corrected/` via
  [`R/corrected/comparison.R`](../R/corrected/comparison.R) `build_methods_comparison()`
  (writes `cv_honesty_compare.csv`, `calibration_oof_compare.csv`,
  `admin2_error_sampling_adjusted.csv`, `decision_value.csv`, `trust_flags.csv`,
  `area_partial_pooling.csv`, `interval_summary.csv`, `methods_comparison.rds`)
  and via [`p9_protocol_reconciliation.R`](../R/corrected/p9_protocol_reconciliation.R)
  `build_protocol_reconciliation()` (`protocol_reconciliation.csv`,
  `protocol_reconciliation_medians.csv`).
- The per-individual survey weight column is `cc$weight_col` (e.g.
  `gw_svy_weight`, `gw_sWeight`, `svy_weight` — see [`R/config.R`](../R/config.R)
  lines 35/132/259/368). It is already correctly threaded in
  `reconcile_protocols()` (p9) as `w_ind`.

---

## 🔴 Issue 1 — FH "partial-pooled" estimate reported in-sample (EB ceiling sold as predictive)

- **File / function / lines:**
  [`R/corrected/p7_area.R:12–75`](../R/corrected/p7_area.R)
  `area_partial_pooling_corrected()`. Covariate selection by `|cor(gee_v, svy_prev)|`
  on the full surveyed set (lines 27–37); FH fit with `train == test == area`
  (lines 43–52, `fit_predict_fh(tr, tr, vars, ...)` then `est <- fh$train_pred`).
- **What's wrong:** `partial_pooled_prev` is an in-sample empirical-Bayes
  smoother whose covariates were chosen using the very `svy_prev` it is
  correlated against, so any reported `partial_pooled_prev` vs `direct_prev`
  agreement is an EB/selection ceiling, not out-of-area predictive skill.
- **Concrete change** — pick (a) for a real number or (b) for honesty-by-labelling:
  - **(a) Leave-one-area-out wrapper (preferred).** Replace the single full-data
    selection+fit (lines 27–61) with a loop over rows of `area`: for each held-out
    Admin2 `i`, (1) rank `gee_cols` by `|cor|` with `svy_prev` **using `area[-i,]`
    only**, (2) `fit_predict_fh(area[-i, ], area[i, , drop=FALSE], vars_i,
    model_type="continuous")`, (3) take `fh$pred` (the synthetic-only out-of-sample
    prediction, which `fit_predict_fh` already returns as `pred`, not `train_pred`)
    into `est[i]`. Keep the EB-shrinkage fallback per-fold when `fh` is `NULL`.
    Add a column `estimate_type = "leave-area-out FH (out-of-sample)"`.
  - **(b) Relabel (minimum).** Keep the in-sample fit but rename
    `partial_pooled_prev` → `partial_pooled_prev_insample`, set
    `attr(out,"method")` to include `"(in-sample EB smoother / ceiling — NOT
    predictive)"`, and add an explicit column `is_in_sample = TRUE`. Then forbid
    any correlation of this column vs `direct_prev` from being reported as skill
    (see downstream).
- **Downstream to regenerate:** target `corrected_area_*` (all slices) →
  `area_pp` field of `corrected_result_*` → `build_methods_comparison()` →
  `results/tables/corrected/area_partial_pooling.csv` and `methods_comparison.rds`
  (and the dashboard copy `dashboard/data/methods_comparison.rds`). The dashboard
  "Methods comparison" tab ([`dashboard/R/mod_methods_comparison.R`](../dashboard/R/mod_methods_comparison.R))
  reads `area_pp`. Note: P9's `sae_glmm_insample` column is the *correct,
  already-labelled* analogue — cross-reference it so the two ceilings are named
  consistently.
- **Verification:** under (a), confirm `est[i]` never uses row `i` — assert that
  re-running with row `i`'s `svy_prev` shuffled does not change `est[i]`. Under
  (b), grep the manuscript/dashboard for any `cor(... partial_pooled ...)` and
  confirm none remains. Sanity: out-of-area FH `pearson_r` should drop well below
  the current in-sample value.
- **Effort / risk:** (a) **M**, medium risk (FH can fail on 1-area-removed folds —
  the existing fallback covers this, but watch degenerate variance, see Issue 4).
  (b) **S**, low risk.

---

## 🔴 Issue 2 — Admin-2 headline tables use UNWEIGHTED prediction means vs design-WEIGHTED survey truth

- **File / function / lines:**
  [`R/corrected/p2_p6_methods.R:60–73`](../R/corrected/p2_p6_methods.R)
  `corrected_admin2()` aggregates OOF predictions with
  `stats::aggregate(yhat_full ~ Admin2, FUN = mean)` (line 63) — an **unweighted**
  district mean — then joins the **design-weighted** `svy_prev`. Consumed by P5
  `admin2_error_corrected()` (76–97), P3 `decision_value_corrected()` (100–123),
  P4 `intervals_corrected()` (131–167).
- **What's wrong:** the predicted district value is a plain mean of individual OOF
  probabilities while the "truth" is the survey-weighted prevalence, so every
  RMSE / Pearson-r / lift / coverage number conflates a real prediction error with
  a weighting mismatch — exactly the confound P9 already isolates
  (`indiv_*_unwt` vs `indiv_*_wt`).
- **Concrete change:**
  1. Give `corrected_admin2()` access to per-individual weights. The OOF frame
     `slfit$honest_*$oof` does **not** currently carry weights, so either:
     - thread weights in at fit time — in
       [`p1_fitting.R:112–117`](../R/corrected/p1_fitting.R) `attach_meta()`, add
       `x$oof$w <- w_ind` where `w_ind` is built exactly as in p9
       (`cc$weight_col`, non-finite/≤0 → 1); **and** pass `cc` through to
       `corrected_admin2()` so it can fall back to weight 1 when absent; **or**
     - re-derive weights in `corrected_admin2()` from `outcome_data` by row index
       (the OOF frame has `row` = original row id; align on it).
     The first is cleaner and reuses the row alignment P1 already owns.
  2. Replace line 63 with a weighted aggregation:
     `agg <- do.call(rbind, lapply(split(d, d$Admin2), function(g)
       data.frame(Admin2 = g$Admin2[1],
                  pred_prev = stats::weighted.mean(g$yhat_full, g$w, na.rm=TRUE))))`.
  3. Keep an **unweighted** variant available for the explicit weighting-confound
     comparison only (mirroring P9), but make `pred_prev` (the value all headline
     tables use) the **weighted** one.
- **Downstream to regenerate:** `corrected_err_*`, `corrected_dec_*`,
  `corrected_int_*`, `corrected_area_*` (it calls `corrected_admin2` at p7_area.R:64),
  `corrected_recon_*` is already weighted-aware (no change), then
  `build_methods_comparison()` → `admin2_error_sampling_adjusted.csv`,
  `decision_value.csv`, `interval_summary.csv`, `area_partial_pooling.csv`,
  `methods_comparison.rds` + dashboard copy. Any manuscript Admin-2 RMSE/r,
  reach@20%, lift, hit-rate@top20, and conformal-coverage numbers.
- **Verification:** for one slice, confirm the new `corrected_admin2()$pred_prev`
  equals `reconcile_protocols()`'s weighted aggregation
  (`indiv_<scheme>_wt` path) on the same OOF object — they should now agree to
  rounding. Confirm P3's `note` no longer claims `weight = n_svy` if the change
  uses individual weights (update the note string at line 121).
- **Effort / risk:** **M**, medium risk — touches the shared helper feeding four
  P-targets; the row-alignment/weight-threading is the only delicate part.

---

## 🟡 Issue 3 — Dormant test-centering leak in area-LOCO

- **File / function / lines:**
  [`R/transportability_area.R:209–214`](../R/transportability_area.R)
  inside `run_area_transport_loco()`. When `recipe$center` is TRUE,
  `Xte <- sweep(Xte, 2, colMeans(Xte), "-")` (line 211) centers the **held-out
  country on its own column means**.
- **What's wrong:** test-set centering on test means is information leakage from
  the held-out country into its own features; it is currently **inert only**
  because both shipped recipes set `center = FALSE`
  ([`AREA_TRANSPORT_RECIPE`](../R/transportability_area.R) line 39,
  `AREA_TRANSPORT_RECIPE_RIDGE` line 47), so a future re-enable would silently
  re-introduce the leak.
- **Concrete change:** the within-country-centering intent for a held-out country
  is undefined out-of-sample. Two safe options:
  - **Preferred:** center the test rows with the **training** column means:
    replace line 211 with `Xte <- sweep(Xte, 2, colMeans(Xtr), "-")` (use the
    *pre-screen, pre-center* `Xtr` means — capture `tr_means <- colMeans(Xtr)`
    immediately after `.tr_prep_X`, before `.tr_center_by` mutates `Xtr` at line
    210). Note this still leaves the train side group-centered, which is the
    documented "predict deviance" recipe; the principled transportable analogue is
    `fit_predict_two_stage()` in benchmark_models.R, so a comment should point there.
  - **Minimum:** delete line 211 (no test centering) and add a guard:
    `if (isTRUE(recipe$center)) stop("center=TRUE for held-out country is not
    well-defined out-of-sample; use fit_predict_two_stage instead")` — making the
    dormant path fail loudly instead of leaking.
- **Downstream to regenerate:** none under the shipped recipes (both `center=FALSE`),
  so **no CSV changes today**. If `AREA_TRANSPORT_RECIPE` is ever switched to
  `center=TRUE`, this gates `area_transport_*` targets →
  `transportability_area_loco_metrics.csv` / `_predictions.csv`.
- **Verification:** with `center=FALSE`, byte-identical CSV outputs before/after
  the patch (regression guard). With a temporary `center=TRUE`, confirm held-out
  predictions change only via train means (or that the guard fires).
- **Effort / risk:** **S**, low risk (inert path; pure hardening).

---

## 🟡 Issue 4 — Inconsistent / degenerate FH sampling variance

- **File / function / lines:**
  [`R/benchmark_models.R:233–295`](../R/benchmark_models.R) `fit_predict_fh()` —
  `deff_default = 1.5` (line 235), `n_e <- pmax(train$n_svy / deff_default, 1)`,
  `sv <- p*(1-p)/n_e` (lines 252–254). Contrast
  [`R/corrected/p7_area.R:54–60`](../R/corrected/p7_area.R) EB fallback, which uses
  `vt <- direct*(1-direct)/pmax(n_svy,1)` — i.e. **deff = 1**. The FH random-effect
  variance (`refvar`/σ²_u from `sae::eblupFH`) can collapse to ~0 for near-constant
  low-prevalence outcomes (e.g. `women_vitA`), shrinking every district to the
  global synthetic mean (~0.001).
- **What's wrong:** the design-effect assumption (and therefore the shrinkage
  weight `γ_i = σ²_u/(σ²_u+v_i)`) is **1.5 in the benchmark FH but 1.0 in the
  corrected EB fallback**, `deff` is meaningless for the many ~1-cluster Admin-2
  areas (sae_area_level / survey_weighting_issues memos), and a σ²_u→0 fit
  produces degenerate near-constant district estimates that are silently reported.
- **Concrete change:**
  1. **Unify deff.** Introduce one constant (e.g. `FH_DEFF_DEFAULT <- 1.5` in
     benchmark_models.R, or 1.0 — pick one and document the choice) and use it in
     **both** `fit_predict_fh()` and the p7_area.R EB fallback (line 56). If the
     project's working assumption is genuinely deff≈1 at Admin-2 (≈1 cluster/area),
     prefer 1.0 everywhere and note it; otherwise pass `deff` explicitly from the
     caller.
  2. **Handle 1-cluster / degenerate `n_svy`.** Floor effective n and the sampling
     variance: after `sv` is computed, `sv <- pmax(sv, <floor>)` and skip/flag
     areas with `n_svy < 2` (cannot estimate a within-area variance).
  3. **Guard σ²_u → 0 collapse.** After `sae::eblupFH` succeeds (line 274),
     read the estimated random-effect variance
     (`fit$fit$refvar` for `sae::eblupFH`) and if it is `< eps` (e.g. `1e-6`) OR
     the convergence flag is bad, **return `NULL` with a logged reason** (so the
     method is skipped rather than emitting ~0.001 everywhere). Alternatively floor
     `refvar` and recompute `γ_i`, but flag+drop is safer for a headline table.
- **Downstream to regenerate:** benchmark FH rows in the benchmark targets
  (`benchmark_targets` in `_targets.R`) and any FH-derived metrics CSV; `women_vitA`
  FH rows specifically (currently degenerate). In the corrected pipeline,
  `corrected_area_*` if its EB fallback deff changes → `area_partial_pooling.csv`.
- **Verification:** for `women_vitA`, confirm FH no longer returns a
  near-constant `train_pred` (sd of `train_pred` should be > the floor, or the
  method should be explicitly skipped with a logged reason). Confirm the EB
  fallback and FH now use the same deff (grep both call sites). Add a unit check
  that `sd(fh$train_pred) == 0` ⇒ flagged.
- **Effort / risk:** **M**, medium risk — `sae::eblupFH`'s slot for σ²_u
  (`fit$fit$refvar`) must be confirmed against the installed `sae` version before
  relying on it; degenerate-area dropping changes which districts appear.

---

## 🟡 Issue 5 — Two conflicting conformal implementations + binary-Y interval mis-scaling

- **Files / functions / lines:**
  [`R/conformal.R:43–282`](../R/conformal.R) `compute_conformal_ci()` —
  `q_level <- ceiling((1-alpha)*(n+1))/n` (line 86), in-sample self-fulfilling
  coverage (lines 97–103, already footnoted), and **binary-Y individual intervals**
  built around `yhat ∈ [0,1]` with an absolute-residual half-width (lines 81,
  94–95, 126–127) that does not respect the binary support of `Y`.
  Versus the correct split-conformal in
  [`R/corrected/p2_p6_methods.R:131–167`](../R/corrected/p2_p6_methods.R)
  `intervals_corrected()` (held-out calibration split, finite-sample level
  `qprob <- min(1,(1-alpha)*(|cal|+1)/|cal|)` line 142, honest held-out coverage
  lines 151–154).
- **What's wrong:** the corrected pipeline ships two different conformal methods;
  `conformal.R` reports an in-sample (mechanically ≥1−α) coverage and constructs
  individual prediction intervals for a **binary** outcome (where a "covered"
  individual interval is not a meaningful predictive statement), while
  `intervals_corrected()` is the methodologically sound district-level
  split-conformal. The two `q_level`/`qprob` formulas also differ
  (`ceiling(...)/n` vs `(...)*(m+1)/m`).
- **Concrete change:**
  1. **Retire `compute_conformal_ci()` from the corrected/honest path.** It is wired
     as `conformal_ci_*` ([`_targets.R:217–228`](../_targets.R)) on the *production*
     `sl_fit`, not the corrected SL — confirm no corrected manuscript table reads it,
     then either (a) leave it for the production pipeline but add a header banner
     "in-sample diagnostic only — not a held-out predictive interval", or (b) drop
     the `conformal_ci_*` targets if unused downstream.
  2. **Single canonical interval estimator** for any honest reporting:
     `intervals_corrected()`. If `conformal.R` must stay, make its `q_level` match
     the finite-sample form `min(1,(1-alpha)*(n+1)/n)` (drop the `ceiling/n`
     variant) so the two never disagree.
  3. **Binary-Y note:** document that for binary `Y`, individual conformal
     intervals are not reported; only **district-aggregated** intervals
     (`intervals_corrected`) are valid, because the predictand is a district
     prevalence (continuous), not an individual label.
  4. **`adjusted_rmse` upward-bias footnote.** In
     [`p2_p6_methods.R:88`](../R/corrected/p2_p6_methods.R),
     `adj_mse <- max(0, naive_mse - samp_mse)` is a truncated (non-negative)
     estimator and is therefore **upward-biased** when the true adjusted MSE is
     near zero; add a one-line code comment and a manuscript/table footnote stating
     that `adjusted_rmse_pp` is a conservative (bias-toward-positive) point estimate.
- **Downstream to regenerate:** `interval_summary.csv` (only if the `qprob`
  formula or the weighting from Issue 2 changes it), `admin2_error_sampling_adjusted.csv`
  (footnote only, no number change), and the `conformal_ci_*` production targets if
  dropped/relabelled. Manuscript conformal-coverage and half-width text.
- **Verification:** grep the codebase + manuscript for `compute_conformal_ci`
  consumers; confirm corrected tables draw only from `intervals_corrected`. Confirm
  `qprob` is identical wherever a conformal level is computed.
- **Effort / risk:** **S–M**, low risk (mostly deprecation + documentation; the
  honest estimator already exists).

---

## 🟡 Issue 6 — No CIs / null reference on any shipped performance metric

- **Files / functions / lines:**
  - [`R/transportability_area.R:198–249`](../R/transportability_area.R)
    `run_area_transport_loco()` — per-country `metrics` data.frame (lines 235–242)
    has `pearson_r`, `spearman_r`, `rmse_pp`, … but **no `ci_lo`/`ci_hi`/`boot`/
    `perm_p`** columns.
  - [`R/corrected/p2_p6_methods.R:76–97`](../R/corrected/p2_p6_methods.R)
    `admin2_error_corrected()` — `pearson_r` with no CI/null.
  - [`R/corrected/p9_protocol_reconciliation.R:95–104`](../R/corrected/p9_protocol_reconciliation.R)
    `reconcile_protocols()` — the matched-protocol r's have no CI/null.
- **What's wrong:** every shipped area/transportability correlation and AUC is a
  point estimate over 14–87 areas (fe_effective_n memo) with no uncertainty and no
  null — so a reader cannot tell signal from sampling noise.
- **Concrete change** — add two columns *pairs* to each correlation/AUC metric:
  1. **Cluster/area bootstrap CI.** Resample the analysis units **with the right
     unit**: for area-level LOCO and Admin-2 metrics, resample **Admin-2 areas with
     replacement** (B ≈ 1000) within the held-out set, recompute the metric, take
     the 2.5/97.5 percentiles → new columns `ci_lo`, `ci_hi`, `n_boot`. For
     individual-aggregated metrics (`admin2_error_corrected`, p9 `indiv_*`), the
     resampling unit is the **cluster** (`cc$cluster_id`) to respect the design.
  2. **Within-country area-label permutation null.** Permute the Admin2 ↔ `svy_prev`
     labels **within each country** (preserving country composition), recompute the
     metric P ≈ 1000 times, and report a one-sided `perm_p` = fraction of permuted
     metrics ≥ observed (and optionally `null_r_95` for the 95th percentile of the
     null). For LOCO this is permutation within the held-out country only.
  3. Wrap both in a single helper, e.g. `area_metric_ci_null(obs, pred, area, country,
     metric = c("pearson","auc"), B = 1000, seed = ...)`, returning a one-row frame
     of `ci_lo, ci_hi, perm_p`, and `cbind` it into each metric frame.
- **Result-producing functions/targets to modify and new columns:**
  - `run_area_transport_loco()` metrics → add `pearson_ci_lo, pearson_ci_hi,
    spearman_ci_lo, spearman_ci_hi, perm_p, n_boot`; written to
    `transportability_area_loco_metrics.csv` by the `area_transport_summary` target
    ([`_targets.R:980–995`](../_targets.R)).
  - `admin2_error_corrected()` → add `pearson_ci_lo, pearson_ci_hi, perm_p`;
    flows to `admin2_error_sampling_adjusted.csv`.
  - `reconcile_protocols()` → add `_ci_lo/_ci_hi/_perm_p` suffixed columns for at
    least the matched `indiv_region_wt` and `area_enet_region`/`sae_glmm_region`
    columns; flows to `protocol_reconciliation.csv`.
  - If any AUC is headline (e.g. `binary_metrics$auc` in
    [`00_corrected_utils.R:189–208`](../R/corrected/00_corrected_utils.R)
    surfaced via `extract_cv_perf_corrected`), add cluster-bootstrap `auc_ci_lo,
    auc_ci_hi` there too → `cv_honesty_compare.csv`.
- **Downstream to regenerate:** all of the above CSVs + `methods_comparison.rds`
  + dashboard copy; every manuscript table reporting an area r / transport r / AUC
  gains CI and a null p-value.
- **Verification:** confirm CIs bracket the point estimate and widen for
  small-area outcomes (women_iron); confirm `perm_p` is ~uniform when `pred` is
  replaced by noise (calibration of the null). Check column presence in each CSV.
- **Effort / risk:** **L**, medium risk — compute cost (B×P resamples per slice ×
  many slices) and the subtlety of choosing the correct resampling unit
  (area vs cluster) per metric; bugs here mis-state uncertainty.

---

## 🟡 Issue 7 — Silent failures + load-order-dependent `%||%`

- **Files / functions / lines:**
  - **`%||%` divergence:** NA-aware definition in
    [`00_corrected_utils.R:16`](../R/corrected/00_corrected_utils.R)
    (`is.null(a) || (length(a)==1 && is.na(a))`) vs is.null-only definitions in
    [`benchmark_models.R:35`](../R/benchmark_models.R) and
    [`transportability_area.R:187`](../R/transportability_area.R). Whichever file is
    sourced **last** wins globally (all define the bare operator in the global env).
  - **Silent failures:** `tryCatch(..., error = function(e) NULL)` followed by
    `next`/`return(NULL)` swallows fits across the pipeline (e.g.
    [`transportability_area.R:174,176`](../R/transportability_area.R);
    [`p9_protocol_reconciliation.R:132,172,187`](../R/corrected/p9_protocol_reconciliation.R);
    benchmark_models.R FH/BYM2/mixed/etc.), and the rollup
    ([`comparison.R:30–35`](../R/corrected/comparison.R),
    [`p9 build_protocol_reconciliation`](../R/corrected/p9_protocol_reconciliation.R)
    `Filter(... nrow>0 ...)`) still writes CSVs from the surviving slices — so a
    degraded run looks complete.
- **What's wrong:** the meaning of `%||%` depends on source order (NA-aware vs
  null-only changes results wherever an operand can be a scalar `NA`, e.g.
  `cc$cluster_id %||% "gw_cnum"`), and failed model fits vanish without a trace, so
  CSVs can silently represent a partial pipeline.
- **Concrete change:**
  1. **One canonical `%||%`.** Define the NA-aware version **once** (keep it in
     `00_corrected_utils.R:16`) and **delete** the duplicates at
     `benchmark_models.R:35` and `transportability_area.R:187`. Confirm sourcing
     order in `_targets.R` so the canonical one is always present; or, safer,
     namespace it (e.g. `.or_default()`) and replace bare-operator call sites.
     Audit call sites where NA-awareness changes behaviour (search `%||%` operands
     that can be scalar `NA`).
  2. **Failure-logging convention.** Replace `error = function(e) NULL` with a
     small logger, e.g. `.log_skip(slice, stage, e)` that appends a row to a
     run-scoped `results/tables/corrected/skipped_fits.csv`
     (`country, outcome, stage, reason, timestamp`) and *then* returns `NULL`.
     Every rollup (`build_methods_comparison`, `build_protocol_reconciliation`,
     `area_transport_summary`) should additionally write a coverage line: number of
     slices attempted vs surviving, so a degraded run is visible at a glance.
- **Downstream to regenerate:** no metric values change **if** the canonical
  `%||%` matches what was effectively in force at the last run (verify!) — but any
  call site where the operative definition flips (null-only ⇄ NA-aware) can change
  a default and therefore a result; treat as a full corrected + transport re-run
  with a before/after diff. New artifact: `skipped_fits.csv` + coverage counts.
- **Verification:** add a startup assertion that `is.na(NA %||% 1)` is `FALSE`
  (NA-aware) and that only one definition exists (`environmentName` / `getAnywhere`
  check, or a test that greps for duplicate definitions). Run the full corrected
  pipeline and confirm `skipped_fits.csv` row count matches the number of `NULL`
  returns observed in logs.
- **Effort / risk:** **S–M**. The `%||%` unification is **higher-risk than it
  looks** — flipping null-only → NA-aware at a call site that intentionally passed
  a scalar `NA` will change a default. Audit every `%||%` operand before deleting
  the duplicates.

---

## Recommended execution order

Ordered to fix correctness-of-the-number issues first, then uncertainty, then
hardening; and to avoid re-running the expensive pipeline more than necessary.

1. **Issue 7 (`%||%` unification + failure logging)** — do first so subsequent
   re-runs are observable and deterministic, and so later changes aren't masked by
   silent skips. (S–M)
2. **Issue 4 (FH sampling variance unify + σ²_u guard)** — fixes a degenerate
   estimator the EB fallback in Issue 1 depends on. (M)
3. **Issue 2 (weighted Admin-2 aggregation)** — the single highest-impact
   correctness fix; feeds P3/P4/P5/P7. Do before Issue 1(a) so the LOAO FH is
   compared against the corrected weighted truth path. (M)
4. **Issue 1 (LOAO FH or relabel)** — depends on 2 and 4. (M, or S if relabel)
5. **Issue 5 (single conformal + footnotes)** — cheap; do after 2 because the
   interval table depends on the corrected aggregation. (S–M)
6. **Issue 6 (bootstrap CIs + permutation nulls)** — do last among substantive
   changes: it adds columns to the *final* corrected/transport metrics, so run it
   after the point estimates are correct to avoid recomputing CIs on numbers that
   will change. (L)
7. **Issue 3 (dormant centering)** — pure hardening, no current output change; can
   land any time, but bundle it with the transport re-run triggered by Issue 6. (S)

After items 2–6, re-run: all `corrected_*` per-slice targets →
`corrected_methods_comparison` and `protocol_reconciliation`; and
`area_transport_*` → `area_transport_summary`. Regenerate every CSV in
`results/tables/corrected/` and `results/tables/transportability_area_loco_*.csv`,
then refresh `dashboard/data/methods_comparison.rds`.

---

## Manuscript numbers that will likely CHANGE

Authors should expect movement in the following (so they are not surprised):

- **Admin-2 error metrics** (`naive_rmse_pp`, `adjusted_rmse_pp`, `pearson_r`) for
  every country×outcome — Issue 2 (weighting) and Issue 5 footnote. Direction not
  guaranteed; r may rise or fall once truth and prediction share a weighting basis.
- **Decision-value table** (`reach_at_20pct`, `lift_vs_no_targeting`,
  `hit_rate_top20`) — Issue 2 re-ranks districts by weighted predictions.
- **Conformal interval summary** (`split_conformal_halfwidth_pp`,
  `empirical_coverage_heldout`) — Issue 2 changes the residuals; possibly Issue 5.
- **Partial-pooling / area-model headline** — Issue 1: the previously reported
  in-sample `partial_pooled` agreement will **drop** to an honest out-of-area
  number (or be relabelled as a ceiling and removed from "predictive" claims).
- **Fay-Herriot benchmark rows**, especially **`women_vitA`** — Issue 4: degenerate
  σ²_u→0 districts will no longer collapse to ~0.001 (or the method will be
  explicitly skipped); the unified deff also shifts shrinkage weights for all FH rows.
- **Transportability area-LOCO metrics** — Issue 6 adds `ci_lo/ci_hi/perm_p`;
  the **point** r's are unchanged by Issue 6 but the *interpretation* (which are
  distinguishable from the null) will change, and several modest r≈0.30 results may
  prove not significant at the area level.
- **Protocol-reconciliation table (P9)** — Issue 6 adds CIs/nulls; the matched
  `indiv_region_wt` vs `area_enet_region` comparison may move from "tied" to
  "indistinguishable from null for both".
- **Any AUC headline** — Issue 6 adds cluster-bootstrap CIs; point AUCs unchanged.

No change expected to: cross-country LOCO **point** estimates from
`run_area_transport_loco()` under the shipped `center=FALSE` recipes (Issue 3 is
inert), nor to upstream data-cleaning / merge outputs.
