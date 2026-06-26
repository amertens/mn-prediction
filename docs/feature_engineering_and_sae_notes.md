# Feature Engineering, Preprocessing & Area-Level Modelling — Notes & Extensions

**Date:** 2026-06-15
**Companion to:** `docs/pipeline_improvements_todo.qmd` (items here are written to fold into that roadmap)
**Evidence base:** sandbox investigation in `sandbox_fe/` — see `sandbox_fe/FINDINGS.md` (preprocessing & feature
engineering) and `sandbox_fe/SAE_FINDINGS.md` (area-level modelling), with figures `FE_summary.png`,
`FE_confirmation.png`, `FE_sae.png` and raw outputs `sandbox_fe/results_*.csv`.

This note records every idea generated while exploring how to preprocess the proxy predictors and frame the
prediction task, with enough evidence and context to act on each one. It is deliberately exhaustive; the
**TL;DR** and **priority table** at the top are the fast path.

---

## TL;DR

1. **The pipeline is implicitly an n = 14–87 problem, not n ≈ 500–1400.** Proxy predictors are almost all
   Admin2-level, so the effective sample size is the number of areas. This single fact explains the weak,
   unstable, occasionally anti-predictive results and dictates the rest of these recommendations.
2. **Ship now (validated, drafted behind default-off flags):** rank/quantile normalization instead of z-score;
   outcome-specific feature bundles; treat aggressive supervised dimension reduction as mandatory, not optional.
3. **Do *not* ship unsupervised per-domain PCA** — the metadata won't support it and its loadings don't transport.
4. **Reframe area maps as small-area estimation (SAE).** Admin2 is below the surveys' design resolution
   (1–3 clusters/area), so model-based estimation is mandatory; report leave-area-out CV, not in-sample fit
   (the optimism gap is large and outcome-specific).
5. **Variable/domain importance** should be computed out-of-fold, in *both* within-country and transfer regimes,
   with unique-vs-shared decomposition and area-level permutation.

## Priority table

| # | Item | Priority | Status | Effort |
|---|---|---|---|---|
| 1 | Report effective n (# areas) + area-blocked uncertainty | **P1** | Planned | M |
| 2 | Area-level SAE (Fay-Herriot → BYM/INLA) as co-primary | **P1** | Prototyped | L |
| 3 | Leave-area-out CV for reported skill & importance | **P1** | Prototyped | M |
| 4 | Rank/quantile normalization (`step_percentile`) | **P1** | Drafted (flag) | S |
| 5 | Out-of-fold + transfer + unique/shared domain importance | **P1** | Planned | M |
| 6 | Outcome-specific feature bundles | P2 | Drafted (flag) | S |
| 7 | Supervised reduction as first-class step (not PCA) | P2 | Validated | S |
| 8 | Fall back to direct/national estimate where CV skill ≤ 0 | P2 | Planned | S |
| 9 | Resolution-matched features (cluster vs Admin2) | P2 | Validated | M |
| 10 | Design effect for SAE sampling variance | P2 | Planned | S |
| 11 | Seasonal construct summaries (transportable features) | P3 | Drafted (unwired) | M |
| 12 | Continuous-then-threshold for rare binary outcomes | P3 | Hypothesis | M |
| 13 | Expand country coverage / multilevel partial pooling | P3 | Strategic | L |

Status: **Validated** = measured in sandbox; **Drafted** = code written behind a default-off flag;
**Prototyped** = working sandbox prototype, not in pipeline; **Planned** = design clear; **Hypothesis** = needs test.
Effort: S/M/L.

---

## 1. The core structural insight — effective sample size is the number of areas

**Finding.** The number of *distinct predictor rows equals the number of Admin2 units*: Sierra Leone 14, Gambia
30, Ghana 72, Malawi 87 — despite 486–1402 individuals. ~100% of GEE variance is between-cluster; granularity
profiling shows the vast majority of the ~1000–1700 columns are Admin2-level (all of DHS/MICS/IHME/LSMS/FSEC/SOIL),
with only ~70–200 GEE/MAP columns at finer GPS-cluster resolution. (`sandbox_fe/07_effective_n.R`,
`08_granularity.R`.)

**Why it matters.** Fitting hundreds of features against 14–87 effective points is ~30× over-parameterized;
apparent n ≈ 1000 precision is pseudo-replication. This is the root cause behind every downstream symptom.

**Actions (P1):**
- Report **# areas alongside n_individuals** in every results table and the manuscript.
- Use **area-blocked CV and area-level bootstrap** for uncertainty (current cluster-blocked resampling still
  over-credits independent observations). Files: `R/bootstrap.R`, `R/aggregation_uncertainty.R`.

---

## 2. Validated preprocessing changes (ready to adopt)

All measured under leakage-free cluster-blocked CV and confirmed through the **real `sl3` SuperLearner** (fast and
full 6-learner stacks). See `sandbox_fe/FINDINGS.md` §validation and `results_04_*`, `results_14*`.

### 2.1 Rank/quantile normalization instead of z-score — **P1, Drafted**
Heavy-tailed environmental predictors hurt linear learners under z-scoring. Rank/percentile mapping gives a
consistent small gain (+0.003–0.012 AUC where it helps; never a material loss; transfers better). Smaller than
the lasso-only effect because the SL stack already contains transform-invariant ranger — but real on the full stack
(Ghana child_vitA 0.539→0.543; Gambia women_iron 0.604→0.607).
- **Drafted:** `FE_NORMALIZE=rank` swaps `step_normalize` → `step_percentile(outside="both")` in
  `DHS_SL_clustered()`. `outside="both"` is essential — the default returns NA for new-country values outside the
  training range and would break prediction. Bake-aware, so the bootstrap/prediction path reproduces it.
  Files: `R/config.R`, `src/analysis/sl_helpers.R`. Default = legacy.

### 2.2 Outcome-specific feature bundles — **P2, Drafted**
Biology-driven predictor subsets match/beat the full set with ⅓–⅛ the features and improve iron transfer.
- `vitA_env = gee_ + soil_ + fsec_`; `iron_health = MAP_ + ihme_ + fsec_`.
- Evidence: Ghana child_iron 0.723 with **40 features** vs 0.702 with 350; women_iron area-transfer Spearman
  0.274→0.323. Specificity (matched > mismatched) holds in 3/4 within-country cells; the exception (Gambia
  women_iron) is because environment features carry broad within-country signal. So bundles are best framed as
  **principled dimension reduction + interpretability**, not a universal AUC booster.
- **Drafted:** `FE_BUNDLES=true`; `bundle_prefixes_for_outcome()` + `Xvars_bundle` from `build_outcome_dataset()`.
  Files: `R/config.R`, `R/data_prep.R`, `R/sl_fitting.R`. Default = legacy.

### 2.3 Aggressive supervised reduction as a first-class step — **P2, Validated**
~27 components (or the cluster-resolution subset) match the full ~1000 features within-country; the extra
features are redundant given 14–87 effective points. Prefer **supervised selection of named features**
(prescreen + L1) or curated indices — the existing prescreen (p<0.2) + corr filter is a reasonable start; consider
tightening and documenting it as required, not optional.

### 2.4 Imputation — leave as is — **Resolved**
Median ≈ mean ≈ zero ≈ KNN within 0.003 AUC; missingness indicators neutral. Do **not** build missForest/MICE/KNN.

---

## 3. Why NOT unsupervised per-domain PCA (metadata critique) — **decided against**

`sandbox_fe/17_pca_metadata_eval.R`. Per-domain PCA is only as good as its grouping, and the metadata doesn't
support it:
- **Prefix "domains" are incoherent** — "GEE" (270 cols) top-3 PCs hold only 55% of variance (mixes precip,
  vegetation, temperature, elevation, cropland, nightlights).
- **Loadings don't transport** — cross-country |cor| of PC1 loadings is 0.22 (temperature), 0.33 (precip), 0.46
  (malaria); this is the mechanism behind observed negative transfer (one component flipped to Spearman −0.72).
- **Insufficient metadata** — only name strings; the conceptual mapping (`metadata/variable_conceptual_domains.csv`)
  helps but construct parsing still leaves gaps; PCA rotations are estimated from 14–87 points; only 4 countries to
  validate transfer.

**Decision:** ship rank-normalization + bundles + interpretable construct summaries instead. An experimental,
construct-grouped, pooled-fit reducer is provided **unwired** in `R/feature_engineering_constructs.R` for review
only — promote only after a proper variable dictionary and validation on more countries.

---

## 4. Area-level modelling & honest uncertainty (SAE) — **P1, Prototyped**

`sandbox_fe/19_area_sae_prototype.R`, `SAE_FINDINGS.md`, `FE_sae.png`. Directly addresses the effective-n problem.

- **Admin2 is below the surveys' design resolution** (1–3 clusters/area) — design-based per-area variances are
  non-estimable, so **model-based SAE is mandatory** for any Admin2 map.
- **Optimism gap:** in-sample cor(area prev, prediction) overstates leave-area-out CV: Ghana child_iron 0.70→0.62
  (+21% skill vs grand mean), Gambia women_iron 0.72→0.41 (+8%), Ghana child_vitA 0.49→0.35 (+1%), Malawi
  child_vitA 0.45→**0.05** (−13%, worse than the mean). **Report leave-area-out CV; the gap is outcome-specific so
  no flat correction works.**
- **Fay-Herriot borrows strength correctly** — few-cluster areas shrink toward the covariate model
  (γ≈0.3–0.6 sparse vs 0.7–0.85 well-sampled); FH tightens the uninformative direct CIs (±0.17–0.21) by 19–33%
  to ±0.12–0.14.

**Actions:**
- **(P1)** Adopt area-level SAE as **co-primary** for Admin2 maps (Fay-Herriot minimum; **BYM/INLA spatial
  multilevel** is the recommended production form — adds neighbour structure and can use cluster-resolution proxies).
  Keep individual-level SL as sensitivity. Build on `R/area_level_comparison.R`, `R/admin2_analysis.R`.
- **(P1)** Switch reported skill & importance to **leave-area-out CV**.
- **(P2)** **Fall back to direct/national estimate where CV skill ≤ 0** (e.g. Malawi child_vitA) — never publish a
  proxy map worse than the mean.
- **(P2)** **Design effect for `psi`:** the prototype uses a binomial sampling variance (DEFF≈1, non-estimable at
  1 PSU/area), making FH CIs slightly optimistic. Estimate a DEFF at Admin1/national level and apply it.

---

## 5. Variable / domain importance — improvements — **P1, Planned**

Current `R/domain_ablation.R` / `R/conceptual_ablation.R`: permutation (no-refit) importance on **in-sample**
predictions, within-country, grouped via `metadata/variable_conceptual_domains.csv`.

**Improvements:**
1. **Out-of-fold (or transfer) importance** — in-sample permutation rewards overfitting; measure on held-out
   predictions. (`domain_ablation.R` line ~200 explicitly uses in-sample.)
2. **Two regimes: within-country AND LOCO transfer** — they disagree and the disagreement is the point.
   Environment is broadly useful *within* a country; for *transfer* the biology-matched domain dominates
   (iron ← malaria/health, vitamin A ← environment). (`sandbox_fe/12_domain_signal.R`, `15_feature_bundles.R`.)
3. **Unique vs shared contribution (commonality analysis)** — domains are spatially confounded (poor areas have
   high malaria *and* low vegetation *and* food insecurity), so marginal/permutation importance double-counts.
   Run leave-one-domain-out (unique) *and* only-one-domain (total), report the gap.
4. **Permute at the area level & bootstrap over areas** — shuffling an Admin2-constant column across individuals
   creates impossible rows; at 14–87 areas, importance ranks are high-variance, so report rank stability not point
   ranks.

---

## 6. Transportability / external validity — **strategic**

- Only **433 columns are shared across all 4 countries** (GEE/IHME/MAP/FSEC/SOIL) — no survey aggregate transfers;
  this is the deployable feature set for a new country. (`sandbox_fe/02_overlap.R`.)
- Transfer works where countries resemble each other — Gambia/Ghana area-burden is rankable from the others
  (Spearman 0.34–0.69), Sierra Leone (14 areas)/Malawi are not. **This is a hard ceiling set by having only 4
  countries, not a preprocessing failure.**
- **(P3, strategic)** The biggest external lever is **more country surveys**; the biggest internal lever is
  multilevel **partial pooling** (borrow strength across countries with country random effects) rather than
  hard "predict a brand-new country" framing.

---

## 7. Backlog / open questions

- **Continuous-then-threshold (P3, hypothesis):** rare binary outcomes (prevalence 1–3% for some women_vitA/b12)
  are very hard; modelling the continuous biomarker and thresholding may use more signal. Untested here.
- **Seasonal construct summaries (P3, drafted unwired):** the monthly climatology (`gee_<var>_<Mon>_10km`) can be
  compressed to transportable mean/amplitude/phase features (`R/feature_engineering_constructs.R`). Roughly neutral
  within-country in tests; revisit for transfer/spatial models. Note: the older prototype
  `R/feature_engineering.R` matches a `monthlyMM` token that does not occur in these names — confirm matchers
  against the actual schema before relying on it.
- **Resolution-matched features (P2):** if predicting below Admin2, prioritize the ~70–200 cluster-resolution
  GEE/MAP columns; Admin2 aggregates add ~no independent information at finer scales.
- **Count vs rate columns:** IHME/MAP `_count`/`_Count` columns are population-scaled and hurt transfer; prefer
  rate/prevalence columns (modest effect; fold into bundle curation).
- **Drop near-duplicate temporal snapshots** (e.g. MAP 202206 vs 202406) before modelling.

---

## 8. Evidence index

| Topic | Sandbox script | Output |
|---|---|---|
| Effective n / granularity | `07_effective_n.R`, `08_granularity.R` | console |
| Shared-feature overlap | `02_overlap.R` | `_shared_*.rds` |
| Imputation × transform × reduction (within) | `04_within.R` | `results_04_within.csv` |
| Reduction strategies (within / LOCO) | `09_reduction.R`, `10_reduction_loco.R` | `results_09/10_*.csv` |
| Area-level transfer | `11_area_level.R` | `results_11_area.csv` |
| Domain-biology signal | `12_domain_signal.R` | `results_12_domain.csv` |
| Real-SL recipe confirmation | `14_sl_recipe_variants.R`, `14b_fullstack_confirm.R` | `results_14*_*.csv` |
| Feature bundles | `15_feature_bundles.R` | `results_15_*.csv` |
| PCA metadata critique | `17_pca_metadata_eval.R` | console |
| Smoke test of drafted changes | `18_smoke_test.R` | console |
| Fay-Herriot SAE | `19_area_sae_prototype.R`, `20_sae_figure.R` | `results_19_sae*.csv/.rds` |
| Figures | — | `FE_summary.png`, `FE_confirmation.png`, `FE_sae.png` |

**Drafted production code (default-off flags):** `R/config.R`, `R/data_prep.R`, `R/sl_fitting.R`,
`src/analysis/sl_helpers.R`; experimental `R/feature_engineering_constructs.R` (unwired).
A/B usage: `FE_NORMALIZE=rank FE_BUNDLES=true Rscript <pipeline entry>`.
