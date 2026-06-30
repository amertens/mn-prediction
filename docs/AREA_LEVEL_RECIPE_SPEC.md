# Area-level prevalence recipe — production spec

Status: implemented in [`R/corrected/area_recipe.R`](../R/corrected/area_recipe.R);
wired as `area_recipe_*` targets in [`_targets.R`](../_targets.R). Arrived at by a
methods-development sandbox (iterations 1–3 + 2b, June 2026). This doc is the
rationale, the recipe, and the integration/migration plan.

## 1. Why — the modelling unit was the binding choice

The shipped corrected pipeline models **individual binary deficiency** and
aggregates to Admin-2. Because almost every predictor is resolved at the district
level, those individual rows are *pseudo-replicated* (the covariates are constant
within an area), and the honest district-ranking signal looked near-chance
(matched-protocol median r ≈ 0.06).

Modelling the **area prevalence directly** — one row per district, proper count
loss, aggressive p≫n selection — is the dominant lever. Honest out-of-sample
district-rank correlation for the signal cells rises from the individual-SL
~0.06–0.30 to:

| Cell | within-country district | within-country region |
|---|---|---|
| Ghana child iron | 0.64 (p<.005) | 0.55 (p<.005) |
| Gambia women iron | 0.60 (p<.005) | 0.51 (p<.005) |
| Malawi child iron | 0.20 (p=.03) | ~0 (ns) |

The earlier "near-chance" headline was a **median across signal + null cells**
*and* an artifact of the individual-level unit — not a property of the data.

## 2. The recipe

| Element | Choice | Evidence |
|---|---|---|
| **Unit** | one row per Admin-2; predict area prevalence | iter 1: dominant lever |
| **Loss** | **binomial (count)**, weight = survey n | iter 2a: same ranking as logit-gaussian, much better calibration (slope ≈0.85 vs 0.75; logit-gaussian over-shrinks) |
| **Features (within-country)** | **enriched**: `gee_` (soil/climate/veg/terrain/human) + **`MAP_` malaria** + **`ihme_` anaemia/stunting** | iter 2a: malaria/IHME rescue high-burden iron (Malawi null → 0.2–0.38) |
| **Features (transport)** | **universal**: `gee_` only | iter 2b: MAP/IHME are country-calibrated and *hurt* transport (mean 0.30 → 0.20) |
| **Selection** | mild **global near-null pre-filter** (drop bottom 70% by marginal \|cor\|) → in-fold supervised top-K screen → **penalised regression (enet)** | iter 3: pre-filter is **leakage-free** (honest re-screen permutation null stays at 0; observed p ≤ .02) and stabilises screening; penalised regression ≈ RF/GBM and is simpler/safer at n = #areas |
| **SAE layer** | intervals + uncertainty **only**; national-anchor the level | iter 2b: covariates carry *ranking*, not *level*; Fay-Herriot shrinking the level toward the covariate mean *increases* absolute error |
| **Distributional variant** | reserve for **vitamin A** (tail outcome) | iter 2a/earlier: continuous-biomarker integration helps the lower-prevalence cells |
| **Late-fusion stacking** | opt-in per-outcome booster | iter 3: wins on distributed signal (Gambia women-iron 0.76) but **collapses sparse signal** (Malawi → 0.03) — never a default, never PCA |

## 3. Two deployment regimes (different covariate sets)

1. **Within-country gap-filling** — interpolate unsurveyed districts in a country
   with *some* survey data. Use the **enriched** features. Honest r ≈ 0.5–0.8 for
   iron/vitamin-A signal cells. This is the primary, defensible product.
2. **Cross-country transport** — predict a country with *no* biomarker data. Use
   **universal** features and score **ranking only** (Spearman; the level does not
   transport — anchor it to a known national prevalence). Honest Spearman ≈ 0.3–0.5
   for epidemiologically-similar countries (West Africa), ~0 for dissimilar
   (Malawi). Report trust/out-of-support flags front and centre.

## 4. Integration

New, **additive** targets (do not modify the existing corrected pipeline):

- `area_recipe_<country>_<outcome>` → `ar_within_country()` (enriched), one row per
  CV scheme (district + region block) with `pearson_r`, `spearman_r`,
  `calib_slope`, `mae_pp`, and bootstrap CI + permutation null
  (`pearson_ci_lo/_hi/_perm_p`).
- `ar_frame_univ_<country>_<outcome>` → `ar_build_frame(..., "universal")`, reused
  by the transport rollups.
- `area_recipe_transport_<outcome>` → `ar_transport_loco()` for outcomes measured
  in ≥3 countries (iron, vitamin A), per-held-out-country Spearman + CI/null.
- Rollups `area_recipe_within` / `area_recipe_transport` →
  `build_area_recipe_tables()` writes
  `results/tables/corrected/area_recipe_within.csv` and
  `results/tables/corrected/area_recipe_transport.csv`.

All metrics ship CI + permutation-null columns (these estimates are fold-seed
sensitive at the margin — Malawi region bounces 0.06–0.27 — so the CI/`perm_p`
are mandatory for interpretation, not optional).

## 5. Learner-library note

The headline change is the **estimator (area-level), not the library**: at
n = #areas a penalised regression matches RF/GBM/MARS and is more stable
(iter 1/3). The richer convex SL is retained as an *opt-in* (`late-fusion`) for
distributed-signal outcomes only. For the **individual-level** SL sensitivity
analysis, the production stack’s apparent strength is largely the leakage from
full-data prescreening; honest within-fold prescreening removes most of it, so a
richer individual library is not a priority.

## 6. Migration / caveats

- The area recipe **complements** the existing corrected/area targets; the
  manuscript should present it as the primary district-ranking estimator and
  reframe the individual-SL result as the pseudo-replication baseline.
- Soil/climate is the transportable backbone; the malaria/IHME enrichment is
  **within-country only** — do not use it for transport or the no-data product.
- PCA loadings flip sign across countries (project memory); never use unsupervised
  PCA for the headline or for transport.
- The SAE/Fay-Herriot layer is for **intervals**, not point-level denoising; the
  level must come from national-anchor calibration.
