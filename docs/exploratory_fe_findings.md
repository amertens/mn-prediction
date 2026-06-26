# Feature-Engineering & Preprocessing Investigation — Findings

**Scope.** A self-contained sandbox (`sandbox_fe/`) that re-uses the production
config/data-prep but iterates *fast* on preprocessing choices, evaluated with
simple learners (lasso, random forest) under honest, leakage-free, cluster-blocked
CV (within country) and leave-one-country-out / LOCO (transfer). Focus outcomes:
**child_vitA, women_iron, child_iron** (the three present in all four countries).
All preprocessing statistics are computed on training folds only.

> The point of this exercise was not to squeeze out the last 0.01 of AUC, but to
> learn *what kind of preprocessing the data structure actually rewards* — and,
> just as usefully, what it does **not**.

---

## TL;DR — what to change in the pipeline

1. **Stop treating this as an n≈500–1400 problem. It is an n=14–87 problem.**
   The predictors are almost all **Admin2-level**, so the effective sample size is
   the number of areas, not individuals. This single fact explains the weak,
   unstable, and occasionally *anti-predictive* transfer results — and dictates
   everything below.
2. **Aggressively reduce dimensionality before modelling.** ~27 features (3
   PCs × domain) match the full ~1000-feature set within country. The other ~970
   columns are redundant given so few effective points.
3. **Replace z-score normalization with rank/quantile normalization** for the
   linear learners. Consistent +0.01–0.02 AUC across all outcomes, for free.
4. **Don't invest in fancy imputation.** Median ≈ mean ≈ zero ≈ KNN (within 0.003
   AUC). Missingness indicators don't help. Median is fine.
5. **For transfer, select predictor *domains* by outcome biology — don't pool
   everything.** Vitamin A transfers on **environment (GEE)**; iron on
   **malaria (MAP) + modelled health burden (IHME)**. Pooling all domains
   *dilutes* the transferable signal, and the best single domain often beats ALL.
6. **Prefer named features + L1 selection over unsupervised PCA for transfer.**
   PCA loadings are not transportable — the same component flipped sign across
   countries (Spearman −0.72 in one case).

---

## I1 — The headline: effective sample size is the number of areas (Panel A)

`07_effective_n.R`, `08_granularity.R`. The number of **distinct predictor rows**
equals the number of Admin2 units:

| Country | Individuals | Clusters | Admin2 | Distinct predictor rows |
|---|---|---|---|---|
| Sierra Leone | 486–774 | 60 | **14** | **14** |
| Gambia | 1020–1402 | 70 | **30** | **30** |
| Ghana | 981–1165 | 90 | **75** | **72** |
| Malawi | 752–1102 | 103 | **87** | **87** |

~100% of GEE feature variance is *between-cluster* (constant within area).
Granularity profiling shows the split clearly: of ~1000–1700 columns, **the vast
majority are Admin2-level** (DHS, MICS, IHME, LSMS, FSEC, SOIL — all of them) and
only **~70–200 GEE/MAP columns** carry finer GPS-cluster resolution.

**Implication.** Fitting hundreds of features against 14–87 independent points is
~30× over-parameterized. The apparent precision from n≈1000 individuals is
pseudo-replication. This is *the* reason transfer is fragile.

## I2 — ~27 components ≈ ~1000 features (Panel B)

`09_reduction.R` (within-country CV, mean over countries):

| Strategy | child_vitA | women_iron | child_iron |
|---|---|---|---|
| full (~1000) | 0.543 | 0.567 | 0.562 |
| per-domain PCA (3 PCs, ~27) | 0.538 | 0.558 | 0.557 |
| cluster-resolution only (~70–200) | 0.529 | **0.580** | 0.563 |

A handful of components reproduce the full model; for **women_iron the
cluster-resolution remote-sensing subset alone *beats* the full set** — the
Admin2 aggregates add noise, not signal.

## I3 — Imputation method is irrelevant

`04_within.R` (ranger; trees are invariant to monotone transforms so this isolates
imputation/reduction). median, mean, zero, and KNN all land within **0.003 AUC** of
each other on every outcome. Missingness indicators are neutral-to-slightly-negative.
**Don't build a missForest/MICE pipeline — median imputation is not the bottleneck.**

## I4 — Rank/quantile normalization beats z-score for linear learners (Panel C)

`04_within.R` (lasso). The heavy-tailed environmental predictors hurt a linear model
under z-scoring; mapping each feature to its rank fixes it:

| Outcome | z-score (current) | rank/quantile |
|---|---|---|
| child_vitA | 0.538 | **0.549** |
| women_iron | 0.528 | **0.548** |
| child_iron | 0.547 | **0.558** |

Consistent direction across all three → credible. The production recipe currently
ends in `step_normalize()` (z-score); a rank/percentile step would help every linear
learner in the SuperLearner stack at no cost. (Trees are unaffected — expected.)

## I5 — Only the globally-available surfaces transfer

`02_overlap.R`. Of ~1000–1700 columns per country, only **433 share names across all
four countries**, and they are exclusively: **GEE 270, IHME 116, MAP 30, FSEC 10,
SOIL 7**. **No DHS/MICS/LSMS survey aggregate is universal** — they require a recent
survey in the target country. The deployable feature set for a *new* country is the
remote-sensing + modelled-surface set, full stop.

## I6 — Trees overfit and anti-predict in transfer; reduction rescues them

`05_loco.R`, `10_reduction_loco.R`. On the full shared set, random forest transfer is
unstable and sometimes catastrophic — **women_iron worst-case AUC 0.38–0.42** (worse
than chance) when held out on Gambia/Malawi. Per-domain PCA before the forest lifts the
**worst-case women_iron transfer from 0.379 → 0.526**. A linear model with L1 + z-score
is already the most *stable* transfer baseline (no rescue needed). Lesson: if you keep
tree learners in the transfer stack, they must be fed a reduced feature space.

## I7 — Unsupervised PCA loadings are not transportable

`11_area_level.R`. Per-domain PCA → lasso gave **Spearman −0.72** transferring to
Gambia (women_iron): the principal component that meant one thing in the training
countries meant the opposite in the target. Named features under L1 selection keep a
fixed meaning and avoid this failure mode. **For transfer, prefer supervised selection
of interpretable features over unsupervised components.**

## I8 — Match domains to outcome biology; pooling dilutes (Panel D)

`12_domain_signal.R` (area-level LOCO, single-domain lasso, mean Spearman):

| Outcome | Best domain | Best | ALL pooled |
|---|---|---|---|
| child_vitA | **GEE (environment)** | 0.280 | 0.272 |
| women_iron | **IHME (health burden)** | 0.300 | 0.274 |
| child_iron | **IHME (health burden)** | 0.210 | 0.165 |

Malaria (MAP) helps iron (+0.24 women_iron) but *hurts* vitamin A (−0.06).
Environment (GEE) drives vitamin A but is weak for iron. The biology is sensible
(VAD ↔ diet/vegetation/seasonality; iron-deficiency anemia ↔ malaria + health
system), and **the best single domain matches or beats throwing everything in.**

## I9 — The area level is the more honest unit (and has a hard ceiling)

`11_area_level.R`. Aggregating to Admin2 prevalence and pooling **~200 areas** across
countries makes the real n explicit. Transfer then works *where countries resemble each
other* — held-out **Gambia (Spearman 0.61–0.69) and Ghana (0.34–0.37)** are rankable;
**Sierra Leone (14 areas) and Malawi are not** (≈0 / negative). No preprocessing fixes
this — it is a function of having only four countries, two of which are dissimilar.

---

## Concrete recommendations for the production pipeline

- **Recipe:** swap `step_normalize` → a rank/percentile transform for linear learners
  (keep raw scale for trees; SuperLearner can hold both representations).
- **Dimension reduction:** add a per-domain PCA (≈3 PCs/domain) or a supervised
  pre-screen down to ~25–50 features *as a first-class step*, not an afterthought.
  Document that this is required by the area-level effective-n, not optional tuning.
- **Imputation:** keep median; stop here (don't add KNN/missForest/MICE).
- **Transfer / transportability models:** (a) restrict to the 433-column
  globally-available set; (b) prefer L1-selected named features over PCA; (c) if trees
  are used, reduce dimensionality first; (d) consider building **outcome-specific
  domain bundles** (VAD→environment, iron→malaria+health) rather than one pooled model.
- **Honesty in reporting:** report the **number of areas** alongside n_individuals, and
  prefer area-blocked uncertainty. Within-country proxy-only discrimination is modest
  (AUC ≈ 0.54–0.58) and cross-country transfer is weak and country-dependent — this is a
  property of the data structure (few areas, four countries), not a tuning failure.

---

## Validation against the production SuperLearner (follow-up a)

The findings above used fast standalone learners. We re-ran the recipe changes through the
**actual `sl3` SuperLearner** (`sandbox_fe/sl_variant.R` mirrors `DHS_SL_clustered`'s exact
preprocessing, then swaps the final step), cross-validated AUC, both the fast 3-learner and the
production 6-learner stacks.

**Fast stack, 5 cells** (`results_14_sl_recipe_fast.csv`), mean AUC: rank **0.586** ≥ current
0.584 > rankpca 0.578 ≈ pca 0.577. Per-cell Δ vs current:

| cell | current | rank | per-domain PCA |
|---|---|---|---|
| Gambia women_iron | 0.609 | 0.600 | **0.616** (p 163→15) |
| Ghana child_vitA | 0.544 | 0.548 | **0.574** (+0.030, p→21) |
| Ghana child_iron | 0.702 | 0.705 | 0.706 (p 350→24) |
| Malawi child_vitA | 0.540 | **0.552** | 0.491 (−0.049) |
| Gambia child_iron | 0.524 | 0.524 | 0.498 |

**Full 6-learner stack, 2 cells** (`results_14b_fullstack.csv`) — confirms the direction:

| cell | current | rank | per-domain PCA |
|---|---|---|---|
| Ghana child_vitA | 0.539 (213s) | 0.543 (274s) | **0.572** (43s, p=21) |
| Gambia women_iron | 0.604 (114s) | 0.607 (144s) | **0.610** (35s, p=15) |

**Confirmed conclusions:**
- **Rank/quantile normalization is a safe marginal win** (+0.003–0.012 in the cells where it
  helps; never a material loss). The benefit is smaller than in the lasso-only test because the
  SL stack already contains transform-invariant ranger — but it is real and holds on the full stack.
- **Per-domain PCA is a high-variance, high-value reduction.** It cut features ~15× and runtime
  ~3–5× (213s→43s) and *won* on every cell tested on the full stack, with a +0.030 win on Ghana
  child_vitA. But on the fast stack it lost on two weak-signal cells (Malawi child_vitA −0.049),
  partly because the univariate prescreen over-compressed to p≈4. **Recommendation:** adopt PCA
  reduction with a floor on #components (skip or relax the post-PCA prescreen) — it is a
  parsimony/stability/speed win, and within-country it is roughly break-even-to-positive.

## Outcome-specific feature bundles (follow-up b)

We turned the domain-biology finding into two bundles and tested them through the real SL
(within-country) and area-level LOCO (transfer). Crucially we also ran each outcome with the
**mismatched** bundle, so any benefit reflects biology, not just smaller feature count.

```
vitA_env    = gee_ + soil_ + fsec_      (environment, soil, food security)
iron_health = MAP_ + ihme_ + fsec_      (malaria, modelled health burden, food security)
```

**Within-country, real SL** (`results_15_bundles_within.csv`):

| cell | all | MATCHED | mismatched |
|---|---|---|---|
| Ghana child_vitA | 0.544 (p291) | **0.564** (p90) | 0.557 |
| Malawi child_vitA | **0.540** (p135) | 0.519 (p78) | 0.473 |
| Gambia women_iron | 0.609 (p163) | 0.593 (p20) | **0.616** (p70) |
| Ghana child_iron | 0.702 (p350) | **0.723** (p40) | 0.704 |

**Area-level LOCO transfer** (`results_15_bundles_loco.csv`, mean Spearman):

| outcome | all (433 feat) | MATCHED bundle |
|---|---|---|
| child_vitA | 0.272 | 0.263 (287 feat) |
| women_iron | 0.274 | **0.323** (156 feat) |
| child_iron | 0.165 | 0.155 (156 feat) |

**Conclusions:**
- **Biological specificity holds in 3 of 4 within-country cells** (MATCHED > mismatched), strongest
  for iron: Ghana child_iron 0.723 with only **40 features** beats the 350-feature full model
  (0.702). The exception is Gambia women_iron, where the environment bundle wins — environmental
  features (the largest, finest-resolution domain) carry broad within-country signal for everything.
- **For transfer, the iron-health bundle clearly improves women_iron** (Spearman 0.274→0.323) using
  one-third the features, and matches for child_iron. The vitamin-A environment bundle matches the
  full set for transfer.
- **Net:** bundles deliver equal-or-better performance with **⅓–⅛ the features**, plus
  interpretability and a biological prior — exactly what the effective-n constraint
  (see I1) calls for. They are best understood as principled dimension reduction, not a universal
  AUC booster (environment features are broadly useful within country).

**Recommended bundles for the pipeline:**
- **Vitamin A** → environment-forward (`gee_`, `soil_`, `fsec_`); malaria/`MAP_` features can be dropped (they hurt VAD transfer).
- **Iron (women & children)** → health-forward (`MAP_` malaria, `ihme_` burden, `fsec_`); this both shrinks the model and improves women_iron transfer.

## Do we have the right metadata for per-domain PCA? (critical evaluation)

`17_pca_metadata_eval.R`. Per-domain PCA is only as sound as the grouping it rests on. We
tested whether the **source-prefix "domain" labels** are an appropriate grouping, whether a
finer **construct** grouping parsed from column names does better, and whether the resulting
loadings are **transportable across countries**.

**Q1 — prefix domains are incoherent grab-bags.** Variance captured by the top-3 PCs:

| domain | #cols | PC1 | PC1–3 |
|---|---|---|---|
| GEE | 270 | 0.32 | **0.55** |
| IHME | 113 | 0.50 | 0.74 |
| MAP | 24 | 0.33 | 0.64 |
| SOIL | 7 | 0.48 | 0.87 |

"GEE" mixes precipitation, vegetation, temperature, elevation, cropland, nightlights, aerosol —
3 PCs keep only 55% of its variance. Reducing "GEE → 3 PCs" both discards a lot and blends
unrelated constructs into uninterpretable axes.

**Q2 — a construct grouping (parsed from names) is far more coherent for some, and reveals real
multidimensionality in others.** cropland PC1 = 0.91, elevation PC1 = 0.90 (essentially rank-1 —
a single mean would do); but **precip 0.33 / temperature 0.46 / vegetation 0.49** have genuine
multi-dimensional *seasonal* structure that one PC cannot summarize.

**Q3 — the loadings are NOT transportable for the climate constructs.** Cross-country |cor| of
PC1 loading vectors (1 = same meaning everywhere, 0 = country-specific):

| construct | mean &#124;cor&#124; (range) | construct | mean &#124;cor&#124; (range) |
|---|---|---|---|
| **temperature** | **0.22** (0.04–0.43) | cropland | 0.88 (0.82–1.00) |
| **precip** | **0.33** (0.13–0.58) | ihme_health | 0.76 (0.54–0.90) |
| **malaria (MAP)** | **0.46** (0.02–0.89) | vegetation | 0.71 (0.45–0.93) |
| GEE (whole prefix) | 0.59 (0.36–0.72) | elevation | 0.61 (0.17–1.00) |

This is the mechanism behind the earlier −0.72 transfer failure (I7): an unsupervised component
of temperature/precip/malaria means something *different* in each country, so projecting a target
country onto training-country axes is unreliable.

**Verdict — we do NOT have the right metadata, and unsupervised prefix-PCA should not be the
production reduction:**
- The only per-column metadata is the **name string**; there is **no curated variable dictionary**
  mapping the ~600 GEE / ~120 IHME columns to constructs. Construct labels must be heuristically
  parsed (here, 31 GEE cols fell to an "other" bucket) — fragile and incomplete.
- Even *with* construct grouping, the climate-construct loadings (precip/temp/malaria) are
  non-transportable, so per-construct PCs would still transfer poorly.
- PCA rotations for area-level blocks are estimated from **14–87 area-points** (high sampling
  variance), and we have **only 4 countries** to judge which components transfer — too few to
  validate empirically.

**Better-justified reductions (what to ship instead):**
1. **Rank/percentile normalization** (`step_percentile(outside="both")`) — fixed, transportable,
   bake-aware. (Use `outside="both"`: the default `"none"` returns NA for new-country values
   outside the training range and would break downstream learners.)
2. **Outcome-specific bundles** — supervised, biological, fixed-meaning features (no data-driven
   rotations to mis-transport).
3. **Interpretable construct summaries** (seasonal mean / amplitude / phase) where seasonal
   structure exists — transportable *by definition* because the formula doesn't depend on the data.
4. If PCA is still wanted: restrict to the **coherent + transportable constructs** (cropland,
   elevation, vegetation, IHME), **fit loadings on pooled multi-country data**, keep a component
   floor — and treat as experimental, not primary.

## Production changes drafted (behind config flags, default = legacy)

Smoke-tested end-to-end (`18_smoke_test.R`): all files parse, both recipe paths run, and the
rank path bakes new data with **0 NAs** (prediction-path safe). With default flags the behaviour
is unchanged from the current pipeline.

| File | Change | Default |
|---|---|---|
| `R/config.R` | `get_pipeline_params()` gains `normalize_method` (`FE_NORMALIZE`) and `use_outcome_bundles` (`FE_BUNDLES`); new `bundle_prefixes_for_outcome()` helper | `zscore`, `FALSE` → legacy |
| `src/analysis/sl_helpers.R` | `DHS_SL_clustered()` + `one_bootstrap()` gain `normalize_method`; recipe swaps `step_normalize` ↔ `step_percentile(outside="both")` | `"zscore"` → legacy |
| `R/data_prep.R` | `build_outcome_dataset()` returns `Xvars_bundle` (outcome-specific proxy subset) | additive; ignored unless used |
| `R/sl_fitting.R` | `fit_sl_models()` reads the two flags, selects bundle vs full `Xvars`, threads `normalize_method` | flags off → legacy |
| `R/feature_engineering_constructs.R` | **NEW, EXPERIMENTAL, unwired:** construct mapping, transportable seasonal summaries, pooled construct-PCA | not called by `_targets` |

**How to A/B against the current pipeline:**
```bash
FE_NORMALIZE=rank  Rscript run_targets.R     # rank/quantile scaling
FE_BUNDLES=true    Rscript run_targets.R     # outcome-specific bundles
FE_NORMALIZE=rank FE_BUNDLES=true Rscript run_targets.R   # both
```

**Deliberately NOT shipped: unsupervised per-domain PCA.** The metadata evaluation showed the
prefix grouping is incoherent and the climate-construct loadings are non-transportable across the
4 countries. The construct-level reducer in `feature_engineering_constructs.R` is provided for
review only, restricted to transportable constructs and requiring a pooled fit; promote it only
after building a proper variable dictionary and validating loading transfer on more countries.

## Files
- `harness.R` — fast leakage-free CV/LOCO harness (ranger ~2s, fixed-λ lasso ~10s).
- `sl_variant.R` — drives the real `sl3` SuperLearner with swappable recipes (current/rank/pca/rankpca/bundle).
- `fe_features.R` — seasonality + count/rate helpers.
- `0*_*.R`, `1*_*.R` — one experiment each; `results_*.csv` — raw outputs.
- `FE_summary.png` — within-country/transfer findings figure; `FE_confirmation.png` — production-SL confirmation + bundles.

### Side-note (verify before relying on)
`R/feature_engineering.R::engineer_admin2_features` is a **prototype** only invoked by
`scripts/prototype_feature_engineering.R`, not by the main `_targets` pipeline. Its
year/anomaly matcher (`_YYYY`) does match the admin2 `gee_a2_*_YYYY` columns, but its
**seasonality matcher looks for `monthly[0-9]{2}`**, whereas the cluster-level monthly
climatology in the merged data is named `gee_<var>_Jan_10km … _Dec_10km` (no `monthlyMM`
token). If seasonality features are ever promoted into the pipeline, confirm the matcher
against the actual `gee_admin2` schema first (see `fe_features.R::.month_groups` for a
matcher that handles the `_Mon_` naming). In our tests seasonality summaries were roughly
neutral within-country, so this is low priority.
