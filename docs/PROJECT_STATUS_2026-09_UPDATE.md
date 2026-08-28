# Project status update — September 2026

Addendum to `PROJECT_STATUS_2026-08_UPDATE.md`, covering the
`analysis-updates-2026-08` branch. Read the June and August documents first.

This note separates results that **changed** from capabilities that are **new**.
Every number below is read from a produced CSV and cited with its file and the
rows it comes from. Two figures quoted in the August note could not be
reproduced from code, and that is stated rather than papered over.

Per-workstream detail with full source maps is in `docs/findings/`. The
regression gate (`scripts/regression_gate.R`) reports all 52 frozen baselines in
`results/tables/frozen_2026-08/` unchanged after every workstream, so nothing
previously published was regenerated.

---

## Part 1: results that changed

### 1. The cross-country transport signal is largely a selection artifact

`build_best_transportable_predictors()` selects its predictors by scoring
candidates on the same leave-one-country-out folds the selected set is then
reported on. `screen_bivariate_loco()`, `select_stepwise_loco()` and
`add_construct_pca_features()` all see every held-out country.

`R/corrected/p10_nested_loco.R` reruns the entire unmodified selection procedure
per outer fold using only the training countries. Scored by the same estimator
in both arms, over 16 folds
(source: `results/tables/corrected/loco_nested_selection_summary.csv`,
`metric = auc`):

| Scorer | Original selection | Nested selection | Delta |
|---|---|---|---|
| glm | 0.6199 | 0.5380 | -0.0818 |
| SuperLearner | 0.5977 | 0.5214 | -0.0763 |

Nested selection is better in 1 of 16 folds under either scorer, and 0 of 16
under equal-country weighting. Four of 16 nested SuperLearner folds fall below
0.5. Women's vitamin A, the strongest cell under the original scheme at 0.7257
glm AUC, falls to 0.5770.

### 2. The eight-SoilGrids spatial GAM: the answer was already in the repository

`sandbox_parsimony/R/34_rescore_spatial_plus_soil.R` re-scored the manuscript's
model under in-fold selection and its output was never propagated. Mean over the
same 16 cells (source:
`results/tables/frozen_2026-08/sandbox_parsimony_out/spatial_plus_soil_rescored.csv`,
grouped by `variant`):

| Variant | Mean Pearson r |
|---|---|
| `locked_published` | 0.3324 |
| `selected_in_fold` | 0.1426 |
| `spatial_only` (no soil features) | 0.2861 |

Honest in-fold selection performs worse than using no soil features at all.
`docs/manuscript_mcn.qmd:666` reports the 0.330 as "substantially above any
other parsimonious comparator we tested"; that clause does not survive. Five
affected passages are itemised in `docs/findings/WS1_nested_loco.md`. **The
manuscript has not been edited.**

### 3. Covariate hygiene is now on by default

`GEE_COVARIATE_HYGIENE` defaults to true as of 2026-08-27
(`R/gee_band_semantics.R`). Across 16 paired folds the pruned set lowered RMSE in
12, worst case +0.380 pp and best case -13.390 pp, and lowered absolute national
bias in 11; correlation metrics were ambiguous in both directions
(source: `results/sensitivity/covariate_hygiene_paired_summary_ws2b_2026-08.csv`).
The decision rests on the error metrics and on the pruned columns not being
physical quantities, not on the correlations. The flip is an isolated commit.

### 4. The national level offset is not an adjustment artifact

Vitamin A was already harmonized. Iron was not: `UNIFORM_TRANSPORT_TAGS` applies
a uniform WHO cut-point to four differently adjusted ferritin columns.
`brinda_adjust_ferritin()` adds a uniform BRINDA CRP+AGP alternative
(source: `results/tables/corrected/uniform_brinda_loco_level_bias.csv`):

| Outcome | Scheme | Mean absolute national bias (pp) | Mean Pearson r |
|---|---|---|---|
| child_iron | configured | 15.087 | 0.1988 |
| child_iron | uniform_brinda | 14.290 | 0.3150 |
| women_iron | configured | 16.080 | 0.0257 |
| women_iron | uniform_brinda | 14.572 | 0.1748 |

Harmonizing the adjustment moves the level bias by under 1.5 pp on a 15 to 16 pp
bias, and raw child ferritin medians still span a factor of 4.85 across
countries afterwards. It does improve the transported spatial pattern
substantially and is recorded as the `uniform_brinda` scheme, not made the
default.

The decomposition explains why. On the log biomarker scale the pure country
intercept is the larger term for three of four outcomes, and for child iron the
training pool sits between 0.526 and 2.679 times the held-out country's level
(source: `results/tables/corrected/level_bias_decomposition.csv`).

### 5. Anchoring is worth a third of the error, and the anchor does not exist

(source: `results/tables/corrected/anchored_transport_loco.csv`, grouped by
`anchor_source`):

| Anchor | Folds | Absolute national bias (pp) | MAE (pp) | RMSE (pp) |
|---|---|---|---|---|
| unanchored | 16 | 8.679 | 13.530 | 16.336 |
| own survey | 16 | 0.000 | 8.876 | 11.411 |
| external | 3 | 44.141 | 44.791 | 46.323 |

A perfect anchor removes 34 percent of the mean absolute error and improves 10
of 16 folds, with Spearman identical across arms to four decimal places; the
code asserts that monotonicity rather than assuming it.

The external row is the finding. Only **3 of 16** primary country-outcome cells
have any external anchor, there is **no iron entry for any of the four
countries**, and the three vitamin A anchors predate their surveys by 15 to 19
years (source: `metadata/anchors/national_anchors.csv`). Anchoring to them is
far worse than not anchoring. The method is sound and the input is missing.

### 6. Covariates hurt at cluster level

A point-referenced spatial model on the existing cluster coordinates, evaluated
leave-one-Admin-2-out over the 21 cells with a real noise ceiling
(source: `results/tables/corrected/cluster_mbg_within_country.csv`):

| Arm | MAE (pp) | RMSE (pp) |
|---|---|---|
| `spatial_only` | 10.087 | 12.663 |
| `matern_spamm` | 9.947 | 12.591 |
| `national_mean` (null) | 10.599 | 13.019 |
| `covariates_only` | 13.873 | 18.773 |
| `spatial_plus_covariates` | 15.106 | 20.701 |

Both covariate arms are worse than doing nothing, and adding covariates to the
spatial model moves MAE from 10.087 to 15.106. Spatial borrowing beats the null
in 13 of 24 cells, mean -0.48 pp MAE. The Matérn arm produced 391
nearly-singular correlation matrix warnings and is reported but not
recommendable.

### 7. The subsample result reproduces; two caveats are added

Rebuilt from the August note's prose specification because no code for it exists
(source: `results/tables/corrected/subsample_contrasts.csv`):

| Stratum | Comparison | Mean delta (pp) | Median | Percent favouring |
|---|---|---|---|---|
| districts retaining clusters | model minus null | +4.8922 | +3.9533 | 11.3 |
| zero-cluster districts | covariate model minus null | -0.5955 | -0.1230 | 57.3 |
| zero-cluster districts | spatial model minus null | -0.7113 | -0.1642 | 56.1 |

The zero-coverage figure reproduces the note's "about 0.7 percentage points".
The "+7.7pp" figure does not: the rebuild gives +4.8922 pp, same direction and
order but not the same number.

Two caveats the note does not carry. The mean is not the typical case: the
zero-cluster median gain is -0.1642 pp against a mean of -0.7113, and only 56
percent of replicates favour the model. And the gain grows as more of the survey
is retained (-0.38 pp at 50 percent retention, -1.38 pp at 80 percent), which is
the opposite of what a cost-saving argument wants.

---

## Part 2: new capabilities

### 8. Nested selection is available as a reportable scheme

`R/corrected/p10_nested_loco.R`, pipeline targets `nested_loco_result` and
`nested_loco_tables`, driver `scripts/run_nested_loco.R` with a glm scorer for
the full grid and an optional SuperLearner scorer for comparability with
published numbers.

### 9. Uniform ferritin adjustment, assay lineage, and a lineage gate

`brinda_adjust_ferritin()` mirrors the RBP function with the clamp sign reversed,
because inflammation raises ferritin rather than depressing it.
`metadata/assay_lineage.csv` is generated by `scripts/build_assay_lineage.R` so
it cannot drift from the code, and the new `assay_lineage_check` target fails the
pipeline when a modelled outcome has no provenance row. Specimen, instrument and
laboratory are UNKNOWN in all 24 rows because nothing in the repository records
them; each carries a source cell saying where to look. They are not guessed.

### 10. Distributional estimator across every outcome

`R/corrected/p12_distributional.R` runs both arms over all 24 cells. It improves
area correlation (+0.0784) and calibration slope and worsens AUC (-0.0285), Brier
skill, MAE and `r_share`
(source: `results/tables/corrected/distributional_paired_signal.csv`).
Recommendation is to keep it as a declared alternative, not the default; all five
largest gains are cells under 10 percent prevalence, as the prototype predicted.

### 11. Cluster spatial model, anchoring, exceedance, calibration curve

`R/cluster_mbg.R`, `R/corrected/p13_anchored_transport.R`,
`R/corrected/p14_subsample.R`, with drivers in `scripts/`. Exceedance
probabilities are produced against `metadata/who_severity_thresholds.csv`, a new
config carrying a source document per row. **Every row is marked
`verified_against_source = FALSE`** and the flag propagates to every output row;
the file records that the women's vitamin A categories are applied by analogy and
that the 20 percent iron figure is a literature convention, not a WHO category.

### 12. Reproducibility infrastructure

`ANALYSIS_PROFILE` (smoke or full) in `get_pipeline_params()`, so the profile
enters the target hash. `scripts/freeze_baselines.R` and the 52-file frozen
baseline with md5 manifest. `scripts/regression_gate.R`, which recomputes nothing
and classifies every table as unchanged, new scheme rows, changed baseline or
missing. `docs/API_DATA_SOURCES.md` recording every endpoint tested and which
credentials exist.

---

## Part 3: defects found

1. **`ranger_low_mtry` hard-codes `mtry = 8`** (`R/sensitivity/mlr3_fitting.R:53`).
   With fewer predictors ranger exits and aborts the whole SuperLearner fit. This
   cost 28 of 32 cells on the first WS1 attempt. It also affects the production
   `loco_best_*` targets, which feed a 4-to-7 variable curated set into the same
   stack and have never been built in `_targets_full`. WS1 guards locally;
   the production fix is filed separately and production behaviour is unchanged.

2. **`.tr_fit_predict()` silently converts a fit failure into the training mean**
   (`R/transportability_area.R`). On a within-survey-centred response that
   fallback is exactly zero, so a failed fit is indistinguishable from a genuine
   null result. It produced a plausible-looking false negative in WS3 before
   being caught. The root cause there was a named response vector from `tapply()`
   making `cv.glmnet` fail; the wider concern about the silent fallback across
   the production area-transport path is recorded, not changed.

3. **Tanzania.** The brief that generated this work treated Tanzania as the
   active fifth country. It was dropped on 2026-08-26 (`R/config.R:507`) and all
   work here is four-country. This is noted because external documents may still
   describe a five-country pipeline.

---

## Part 4: what was not done, and why

- **Displacement-integrated Earth Engine extraction at cluster points.**
  Dropped before any Earth Engine time was spent. `FINDINGS.md` section 12
  measured the noise ceiling falling as the spatial level gets finer
  (Admin-1 0.59, Admin-2 0.31, cluster 0.22), and the covariate arms it would
  have improved lose to a covariate-free smoother by about 5 pp of MAE
  (section 6 above).
- **WorldPop population weighting.** Not retrieved. No population-weighted
  aggregation exists on this branch; WS5 and WS7 use cluster count or survey n
  and record that on every row.
- **GBD extract.** No GHDx credential is present. The placeholder
  `dashboard/data/gbd_estimates.rds` is left in place rather than deleted, since
  deleting a dashboard input was flagged for confirmation.
- **Admin1+Admin2 join-key migration, Malaria Atlas refresh, SoilGrids refresh,
  GFDx covariates, candidate-country pre-extraction, hemoglobin track.** Not
  attempted. No result on this branch depends on any of them. The hemoglobin
  track is gated on `ENABLE_HB_TRACK`, which is unset.

---

## Part 5: what follows

The evidence now points consistently in one direction. Covariate-based transport
is close to chance once selection is honest (section 1), covariates actively hurt
at cluster level (section 6), and the residual error is dominated by a country
intercept no covariate can supply (section 4). The one intervention that
addresses that intercept is anchoring, and it works: a third of the mean absolute
error (section 5).

The binding constraint is therefore not modelling. It is that no current national
iron estimate exists for these countries in any source available here. Acquiring
a usable anchor is worth more than any further covariate work, and it is the
cheapest of the remaining options.

The August note's reframing toward a survey-design decision aid survives, with
the sharper statement that the model's value is confined to districts you
deliberately do not sample, is worth a few tenths of a percentage point on
average, and is absent in nearly half of replicates.
