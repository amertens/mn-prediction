# Transportability (LOCO) methods

**Project:** Subnational micronutrient-deficiency prediction from proxy indicators
**Countries:** Gambia (2018), Ghana (2017), Sierra Leone (2013), Malawi (2015–16), Tanzania (TDHS 2009–10, vitamin A only — see §4)
**Outcomes:** child & women iron (ferritin), vitamin A (RBP), folate, B12 (women)
**Last updated:** 2026-08-18 (Tanzania added as a fifth country; see §4)

This document catalogs every method tried for the **leave-one-country-out (LOCO)** transportability evaluation — the test of whether a model trained on three countries can predict Admin-2 deficiency prevalence in the held-out fourth.

---

## 1. The LOCO setup

- **Unit of analysis:** Admin-2 area (district). Each country contributes 14–87 areas; pooled n ≈ 200–220.
- **Target:** survey-weighted deficiency prevalence per Admin-2 (`svy_prev`), with sample size `n_svy`.
- **Predictors:** ~150 common GEE/remote-sensing covariates screened to those present across all four countries (`.screen_vars`), optionally augmented with outcome-aware interaction terms (`fe_*`).
- **Protocol:** for each held-out country, train on the other three's Admin-2 rows, predict the held-out country's Admin-2 prevalences. The held-out country's `svy_prev` is **never** used in training (domain-adaptation methods may use the held-out *covariates* but never its labels).
- **Why it's hard:** the effective sample size is the *number of areas* (p ≫ n), and the four surveys differ in biomarker assay, inflammation adjustment, and year — so both signal and comparability are limited.

### Evaluation metrics
| Metric | Meaning |
|---|---|
| Pearson r / Spearman ρ | rank/linear agreement of predicted vs observed Admin-2 prevalence (scale- & shift-invariant → measures **pattern** transfer) |
| MAE / RMSE (pp) | absolute error in percentage points (measures **level** transfer) |
| mean bias (pp) | systematic over/under-prediction |
| calibration slope/intercept | calibration of predicted vs observed |

A method can have decent **r** (ranks districts correctly) yet large **bias/MAE** (wrong absolute level) — and in this project the failures are predominantly of the *level* kind.

---

## 2. Methods tried

### Baseline
- **`baseline`** — country-mean: predict every held-out district at the pooled training mean prevalence. Reference for "no spatial information."

### Frequentist regressions
- **`glm`** — GLM on the covariates (continuous Gaussian and logit variants).
- **`penalized`** — elastic-net (glmnet, α=0.5), `n_svy`-weighted; the standard regularized regression for p ≫ n.
- **`mixed`** — random-effects/mixed model (country random intercept) on the logit of prevalence.
- **`gam`** — generalized additive model (smooth covariate terms).
- **`quantile`** — quantile (median) regression; robust to outliers.
- **`betareg`** — beta regression for proportion outcomes.
- **`quasibinomial`** — quasi-binomial logit GLM with `n_svy` weights (allows over/under-dispersion).

### Variable selection / invariance filters
- **`forward`** — forward-selection regression within a LOCO loop (greedy add-by-CV).
- **`dag`** — DAG-informed covariate selection (uses a biology-driven variable set per outcome).
- **`invariance_filter`** — keep covariates whose covariate→outcome relationship is *stable across training countries* (invariant-prediction idea), to favor transportable signal.
- **`combined_filter`** — combines invariance with cross-outcome information.

### Small-area estimation (SAE)
- **`fh` — Fay–Herriot.** Area-level random-effects model: `svy_prev ~ X·β + u_area`, with the area's known sampling variance. Borrows strength via shrinkage. The canonical survey-statistics method for small areas.
- **`bym2` — BYM2 via INLA.** Bayesian areal model with structured (ICAR, spatial) + unstructured area random effects (Riebler PC priors). For LOCO, the graph is block-diagonal across training countries; held-out areas get the covariate (fixed-effect) prediction. *Fit with `safe = FALSE`* to avoid an INLA rerun deadlock on near-constant rare outcomes.

### Spatial models
- **`spatial_coords`** — covariates + a smooth function of longitude/latitude (GAM spline, k≈30); models broad spatial gradients.
- **`spatial_plus_h6`** — spatial coordinates plus the "H6" engineered cross-outcome feature stack.
- **`spatial_plus_soil`** — spatial coordinates plus top SoilGrids features (soil adds +~0.045 LOCO r over spatial-only).

### Ensembles & two-stage
- **`stacked`** — convex (NNLS) stacker over base methods (baseline, penalized, mixed, gam), with **country-disjoint** inner-LOCO weighting — a transportability-aware SuperLearner stacker.
- **`two_stage`** — non-spatial ensemble in stage 1, then a spatial GAM on the stage-1 residuals (smoothing/mapping).
- **`sl_prescreened`** — the project's optimized SuperLearner: five-stage prescreening + spatial coordinates + a fast learner library. This is the "individual-spirit" ensemble adapted to the area level.

### Domain adaptation
- **`coral` — CORrelation ALignment** (Sun, Feng & Saenko 2016). Aligns the *covariance* of the training covariates to the held-out country's covariate covariance (whiten-source / recolor-to-target), then fits elastic-net on the aligned features. Label-honest (uses only held-out *covariates*). Directly attacks cross-country feature-scale mismatch; size-conditional mean-alignment for small held-out countries. Fixes scale/units shift but **not** an outcome-level offset or covariate→outcome sign flip.

### Dropped for tractability
- **`grouplasso`** (grpreg) — group-lasso. **Removed from the default** (2026-06-17): pathological convergence (multi-hour, single compiled call) on the near-constant women-vitamin-A outcome. Per-method timing confirmed it was the *only* pathological method; pass it explicitly to re-enable.

---

## 3. What the LOCO results show (current)

**Overall: median Pearson r ≈ 0.10, median MAE ≈ 13.7 pp.** Transport is weak across the board.

Median LOCO r by method (best → worst, abbreviated):

| Rank | Method | median r | best r (any split) |
|---|---|---|---|
| 1 | spatial_plus_soil | 0.33 | 0.66 |
| 2 | spatial_coords | 0.27 | 0.70 |
| 3 | spatial_plus_h6 | 0.21 | 0.63 |
| 4 | coral | 0.14 | 0.64 |
| 5 | combined_filter | 0.12 | 0.63 |
| … | gam / two_stage / glm / fh | 0.08–0.11 | 0.6–0.7 |
| low | **sl_prescreened** | **0.06** | 0.46 |
| low | bym2 / penalized / forward / dag / betareg | ≤0.04 | — |

**Key takeaways**
1. **Geography transports best.** Spatial coordinates + soil are the top methods — physically comparable signals (climate/soil) port across countries better than survey-derived structure.
2. **The flexible SuperLearner is among the weakest** for transport, and is beaten by a non-SL method on *every* outcome × held-out split. Flexibility doesn't pay when n = areas.
3. **No method is reliably good** — even the best top out around r ≈ 0.3 median.

### Four distinct failure modes (from the harmonization investigation)
| Mode | Where | Cause | Fixable? |
|---|---|---|---|
| Level offset, ranking OK | Gambia iron | genuinely ~6× lower ferritin (66.6% deficiency corroborated externally) | No — real biology |
| Level artifact | folate / B12 | units/cutoff/assay differences (Malawi folate units bug — **fixed**) | Yes (data) |
| Pattern non-transfer | Malawi / SL iron | genuine West- vs East-Africa difference; r ≈ 0 | No |
| Weak signal | vitamin A (both) | RBP harmonized but proxies don't predict VAD pattern | No |

---

## 4. Adding a fifth country (Tanzania, TDHS 2009-10)

Tanzania contributes **vitamin A only** (RBP-based VAD, women + children). Iron
is excluded because TDHS 2010 assayed sTfR, not ferritin, so it is not
comparable to the definition the other four use; iodine does not fit the
individual-to-Admin-2 prevalence framework.

Two things had to be true before it could enter the LOCO at all, and neither was
true when the merged dataset was first built.

### 4.1 The covariate vocabulary had to be reconciled

Every pooled/LOCO builder takes a strict **intersection** of covariate names
across countries. The four original countries' Admin-2 covariates are zonal
means over local `.tif` exports, named after the **raster filename**
(`gee_ndvi_2013`, `gee_soilzinc_mean_0_20`). Tanzania's came from the Earth
Engine API, named after **EE bands** (`gee_a2_NDVI`,
`gee_a2_SoilZinc_mean_0_20`). The two vocabularies overlap **0 of 149**.

Switching Tanzania on by config alone would therefore have silently cut the
area-level LOCO covariate set from **149 to 0**, and the individual-level pooled
proxy set from **585 to 129** (its `gee_` block from 419 to **0**) — for every
country, with no error raised. This is worth stating plainly because it is
invisible in the outputs: the models still fit, and still report metrics.

**Resolution.** `scripts/build_gee_legacy_parity.R` (built on
`src/GEE/extract_gee_legacy_parity.R`) pulls the same EE assets and emits the
**legacy** column names, so a country with no local rasters joins the existing
vocabulary directly. `extract_gee_admin2()` / `extract_area_covariates()` fall
back to it automatically, `_targets.R` tracks the CSV's checksum so a rewrite
invalidates downstream, and `build_pooled_dataset()` now logs per-domain
per-country covariate counts and warns when a domain contributes zero.

**Validation is the load-bearing part.** The extractor was checked against Sierra
Leone, which has both representations: each column must reproduce the
raster-derived values at **r >= 0.9 AND a mean within 10%**. Rank agreement alone
is not sufficient — a constant scale offset on one country is read by the pooled
model as a genuine country difference, which is precisely the failure mode
catalogued in section 3. The first version of the extractor passed only **22 of
148** columns; the failures were three silent scale errors (a year-over-year
delta band averaged into the base value, halving TRMM; a missing iSDAsoil
`exp(x/10) - 1` back-transform; and mask-vs-zero semantics on WSF). After the
fixes: **135/135 pass**
(`results/sensitivity/gee_legacy_parity_validation_sierraleone.csv`).

**Cost of admitting Tanzania**, to be reported rather than buried:

| Set | Without Tanzania | With Tanzania |
|---|---|---|
| Area-level (Admin-2) LOCO GEE covariates | 149 | 135 |
| Individual-level pooled `gee_` covariates | 419 (Admin-2 + cluster-buffer) | 135 (Admin-2 only) |

A further consideration for the individual-level pool: rows are **not**
country-weighted, and Tanzania's vitamin A sample (6,238 children) is five to
thirteen times any other country's. In every LOCO fold holding out a small
survey, Tanzania would supply roughly two-thirds of the training rows. A
country-weighted sensitivity run is warranted before reading anything into the
individual-level numbers.

The 14 dropped Admin-2 variables are `gee_accessibility_2019` (EE asset
corrupted server-side), the five `gee_esa_worldcereal_2021_*` summaries, seven of
nine `gee_soiltotalcarbon_*`, and `gee_ndvi_2022` (beyond AVHRR CDR coverage).
The 186 **cluster-buffer** variables come from an analyst EE Code Editor export
that is not in the repo, and the `fpp_`/`tpp_` families cannot be identified from
the column names alone, so Tanzania cannot reproduce them. That affects only the
individual-level SuperLearner, which is the documented **sensitivity** analysis
(`sensitivity/README.md`) — the **primary area-level SAE and its LOCO are
unaffected**. Recovering the original export script would close the gap.

### 4.2 The vitamin A outcome is NOT method-identical

The other four countries use two-marker BRINDA (CRP + AGP) on raw RBP. TDHS 2010
assayed **no AGP**, and CRP on only ~27% of the RBP sample. Tanzania had been
falling through a warning path in `apply_brinda_vita_binary()` and silently
keeping the DHS-supplied `rbpadcrp` instead.

It now takes an explicit, declared path: a **CRP-only** BRINDA correction where
CRP exists, falling back to the survey agency's own adjusted RBP where it does
not, rather than leaving those rows raw and mixing corrected with uncorrected
values in one outcome. `brinda_country_method()` reports what each country
actually received, and the method is printed on every run.

Harmonized prevalence: **child VAD 23.3%, women 7.2%** (previous configured
binaries: 23.9% / 7.2%).

**Carry into the manuscript:** Tanzania's vitamin A outcome is a *weaker*
inflammation correction than the other four. Given that section 3's dominant
failure mode is a level offset rather than a ranking failure, a residual level
difference for Tanzania should be interpreted as partly methodological, not
purely biological. Two further caveats compound this: TDHS 2010 has no
micronutrient-subsample weight (the household weight is used, as for Gambia), and
2010 sits four to eight years before the rest of the panel.

---

## 5. Honest conclusion

Cross-country transport of **absolute** prevalence is **not deployable** — the dominant failure is a level/calibration offset driven by real biological differences and cross-survey biomarker comparability, neither of which covariate-based models can recover. Transport of the **spatial pattern** (rank) is modestly better, and best captured by geography/soil features. The defensible deliverable is therefore **within-country** area-level prediction (anchored on the held-out country's own survey), with cross-country results reported as rank-scale, interval-bounded, and explicitly caveated.

*Implementation:* `R/benchmark_models.R` (`run_area_benchmarks_loco` + `fit_predict_*`); results in `results/tables/benchmarks_all.csv` and `transportability_area_loco_metrics.csv`.
