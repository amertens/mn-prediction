# WS1. Selection-honest transport numbers (nested LOCO)

## The defect

`build_best_transportable_predictors()` (`R/transportability_best_model.R:247`)
chooses its predictor set by scoring candidates on leave-one-country-out folds,
and the chosen set is then reported on those same folds. Three steps inside one
call see every held-out country:

| Step | Location | What it sees |
|---|---|---|
| `screen_bivariate_loco()` | `R/transportability_best_model.R:84` | ranks candidates by mean absolute AUC departure from 0.5 across all held-out countries |
| `select_stepwise_loco()` | `R/transportability_best_model.R:199` | greedily adds whichever candidate maximises mean LOCO AUC across all held-out countries |
| `add_construct_pca_features()` | `R/transportability_best_model.R:157` | fits PCA loadings on the pooled matrix of all countries |

The third is acknowledged in its own docstring as "a mild leakage in the
SELECTION/SCREENING step". The first two are not. The reported transport number
is therefore not out-of-sample with respect to the selection.

This repository has already measured the size of that bias once, for a
different method. `sandbox_parsimony/R/34_rescore_spatial_plus_soil.R` moved the
selection of the eight locked SoilGrids variables inside the fold and reported
the two side by side. WS1 does the same for the stepwise predictor set and
adopts that file's variant structure.

## What was implemented

`R/corrected/p10_nested_loco.R`. For each held-out country, both the merged data
and the country configs are restricted to the training countries, and the
**entire unmodified selection procedure** is rerun on that subset, so its
internal folds become leave-one-of-the-training-countries-out. The selected set
is then scored once on the held-out country.

The nested path calls `build_best_transportable_predictors()` itself rather than
reimplementing the screen, dedup, PCA and stepwise search. Reimplementing them
would risk the two schemes differing for reasons other than the nesting, which
is the one thing this comparison must not confound.

Pipeline targets `nested_loco_result` and `nested_loco_tables` are wired in
`_targets.R`, depending on `merged_gambia`, `merged_ghana`, `merged_malawi`,
`merged_sierraleone` and `pipeline_params`. Driver:
`scripts/run_nested_loco.R`.

Outputs, all new files:

- `results/tables/corrected/loco_nested_selection.csv`
- `results/tables/corrected/loco_nested_selection_paired.csv`
- `results/tables/corrected/loco_nested_selection_summary.csv`
- `results/tables/corrected/loco_nested_selected_vars.csv`

Scope: 4 outcomes by 4 held-out countries, 16 folds. Two scorers, and for the
glm scorer two weightings.

## Result 1: the stepwise predictor set

Mean LOCO AUC over 16 folds, glm scorer, equal-row weighting
(source: `results/tables/corrected/loco_nested_selection_summary.csv`, row
`scorer = glm`, `weighting = equal_row`, `metric = auc`):

| Scheme | Mean AUC | Folds better | Mean paired delta | Range of delta |
|---|---|---|---|---|
| `original_selection` | 0.6199 | 15 of 16 | | |
| `nested_selection` | 0.5380 | 1 of 16 | -0.0818 | -0.2726 to +0.0027 |

Under equal-country weighting the picture is the same: 0.6201 against 0.5391,
mean delta -0.0810, nested better in 0 of 16 folds (same source, row
`weighting = equal_country`).

Per outcome (source:
`results/tables/corrected/loco_nested_selection_paired.csv`, rows with
`metric = auc`, `scorer = glm`, `weighting = equal_row`):

| Outcome | Original | Nested | Delta |
|---|---|---|---|
| child_vitA | 0.5906 | 0.5229 | -0.0677 |
| women_vitA | 0.7257 | 0.5770 | -0.1487 |
| child_iron | 0.5896 | 0.5298 | -0.0597 |
| women_iron | 0.5737 | 0.5224 | -0.0513 |

The women's vitamin A cell carries the largest correction. Its apparent
discrimination of 0.7257 is the strongest transport signal in the set under the
original scheme, and it is also the one that falls furthest, to 0.5770. Its
worst single fold is Ghana, 0.7970 to 0.5244, a delta of -0.2726 (same source,
row `outcome = women_vitA`, `held_out = Ghana`).

Three of the 16 nested folds fall below 0.5, against one of 16 under the
original scheme (same source, column `nested_selection` and
`original_selection`).

The selected sets themselves differ per fold
(source: `results/tables/corrected/loco_nested_selected_vars.csv`):

| Scheme | Held out | Predictors selected | Candidates screened |
|---|---|---|---|
| `original_selection` | ALL | 5 | 585 |
| `nested_selection` | Gambia | 7 | 803 |
| `nested_selection` | Ghana | 6 | 585 |
| `nested_selection` | Sierra Leone | 4 | 634 |
| `nested_selection` | Malawi | 6 | 589 |

The candidate count varies because the common predictor set is an intersection
across the countries in scope, and dropping a country can enlarge it. A
predictor selected on three countries can be absent from the fourth; those are
dropped at scoring time, which is what deployment in a genuinely new country
would force. Mean `n_scored` is 3.75 against `n_selected` of 4 to 7 (source:
`results/tables/corrected/loco_nested_selection_summary.csv`, row
`metric = n_scored`).

## Result 2: the spatial-plus-soil model already had its answer in the repository

`sandbox_parsimony/R/34_rescore_spatial_plus_soil.R` had already re-scored the
manuscript's eight-SoilGrids spatial GAM under in-fold selection. Its output was
never propagated. Mean over the same 16 outcome-by-holdout cells
(source: `results/tables/frozen_2026-08/sandbox_parsimony_out/spatial_plus_soil_rescored.csv`,
grouped by `variant`):

| Variant | Mean Pearson r | Mean Spearman | Mean RMSE (pp) |
|---|---|---|---|
| `locked_published` | 0.3324 | 0.3164 | 18.262 |
| `selected_in_fold` | 0.1426 | 0.1170 | 18.238 |
| `spatial_only` | 0.2861 | 0.2739 | 16.497 |

Honest in-fold selection of the soil features (0.1426) performs worse than
using no soil features at all (0.2861). The published figure of 0.3324 is
attributable to selecting the eight variables on the metric they are reported
against. The covariate-free spatial smoother also has the lowest RMSE of the
three.

## Which manuscript sentences this bears on

The manuscript is not edited by this workstream. The following passages are
affected.

1. `docs/manuscript_mcn.qmd:666` (paragraph opening at line 654) states that a weighted generalised additive
   model with a thin-plate spline plus the eight soil features "achieved mean
   LOCO Pearson r = 0.330 across 16 outcome-holdout combinations on the four
   main outcomes -- substantially above any other parsimonious comparator we
   tested". The 0.330 corresponds to `locked_published` (0.3324). Under in-fold
   selection the figure is 0.1426, and the comparator with no soil features at
   all reaches 0.2861, so the "substantially above any other parsimonious
   comparator" clause does not survive.

2. `docs/manuscript_mcn.qmd:658` states that univariate leave-one-country-out
   scoring "identified eight SoilGrids variables ... as the most consistently
   transportable single predictors". That ranking is itself computed on the
   evaluation folds. The claim is about which variables the procedure picks, so
   it remains true as a description of the procedure, but it cannot be read as
   evidence that those variables transport.

3. `docs/manuscript_mcn.qmd:671` reports that bootstrap 95 percent confidence
   intervals on the Pearson r "were strictly above zero for seven of the sixteen
   holdouts". Those intervals are computed on the locked, fold-contaminated
   feature set. The equivalent count under in-fold selection is **not yet
   computed**; the rescore file carries point estimates, not bootstrap
   intervals.

4. `docs/PROJECT_STATUS_2026-08_UPDATE.md` section 2 reports "Overall mean LOCO
   AUC with this set: 0.5718 (16 folds = 4 outcomes x 4 held-out countries)" for
   the five-variable parsimonious set, with an aggregate 95 percent interval of
   [0.544, 0.600] said to exclude chance. Under nested selection scored by the
   same estimator family the mean is 0.5214, with 4 of 16 folds below 0.5
   (Result 3). The interval quoted in that note is computed on a predictor set
   chosen using the folds it is evaluated on, so it does not support the
   exclusion-of-chance claim. The nested equivalent of that interval is **not
   yet computed**; producing it needs the same aggregation code applied to the
   per-fold nested rows.

5. The same note's claim that "for this outcome/country set, more data domains
   is not the bottleneck" is unaffected. It compares candidate pools under one
   selection scheme, and nesting moves both arms.

## Result 3: the SuperLearner scorer

The glm scorer above is fast and matches the metric the search optimises, but
the published transport numbers were produced by the production SuperLearner, so
a claim about those numbers has to be scored the same way. All 32 SuperLearner
cells were fitted under `PIPELINE_MODE=full`, with 8 learners per fit and no
failures (source: `results/tables/corrected/loco_nested_selection.csv`, rows
with `scorer = sl`, columns `note`, `sl_learners_used`, `pipeline_mode`).

Mean LOCO AUC over 16 folds (source:
`results/tables/corrected/loco_nested_selection_summary.csv`, row `scorer = sl`,
`metric = auc`):

| Scheme | Mean AUC | Nested better in | Mean paired delta | Range of delta |
|---|---|---|---|---|
| `original_selection` | 0.5977 | | | |
| `nested_selection` | 0.5214 | 1 of 16 folds | -0.0763 | -0.2267 to +0.0041 |

Per outcome (source:
`results/tables/corrected/loco_nested_selection_paired.csv`, rows with
`metric = auc`, `scorer = sl`):

| Outcome | Original | Nested | Delta |
|---|---|---|---|
| child_vitA | 0.5713 | 0.5187 | -0.0527 |
| women_vitA | 0.6461 | 0.5196 | -0.1265 |
| child_iron | 0.5856 | 0.5219 | -0.0637 |
| women_iron | 0.5879 | 0.5256 | -0.0623 |

Four of the 16 nested folds fall below 0.5, against none of 16 under the
original scheme (same source, columns `nested_selection` and
`original_selection`).

The two scorers agree on both the direction and the approximate size of the
correction: -0.0763 under the SuperLearner against -0.0818 under the glm. Under
nested selection the SuperLearner mean of 0.5214 is close enough to 0.5 that the
per-outcome ordering carries little information.

**Comparability with the published figure.**
`docs/PROJECT_STATUS_2026-08_UPDATE.md` section 2 reports mean LOCO AUC 0.5718
over the same 16 folds for the five-variable set. The `original_selection`
figure here is 0.5977, not 0.5718. The two are close but not identical, and
three differences are known:

1. The predictor data has changed since that note was written. The 2026-08-27
   cross-country linkage fixes altered the candidate pool, so the selection
   procedure does not necessarily return the same five variables it did then.
   This run's `original_selection` chose 5 predictors from 585 candidates
   (source: `results/tables/corrected/loco_nested_selected_vars.csv`).
2. `ranger_low_mtry` is dropped from the library in every cell here, for the
   reason described below. The published run either included it or failed; which
   of the two is **not yet computed**, because the `loco_best_*` targets have
   never been built in the store.
3. The published number came from `run_loco_cv()`, which computes the metric
   across all folds of a single fitted object per outcome. This runs one fit per
   held-out country per scheme, which is what a per-fold variable set requires.

None of those differences affect the comparison being made. Both schemes here
are scored by the same estimator on the same test rows, and the delta between
them is the quantity of interest. The 0.5977 is not offered as a replacement for
0.5718.

## An unrelated defect found while doing this

The first attempt at the SuperLearner scorer returned metrics for 4 of 32 cells.
The other 28 failed with "mtry can not be larger than number of variables in
data. Ranger will EXIT now", which aborts the whole SuperLearner fit rather than
the one learner. The cause is `ranger_low_mtry`, which hard-codes `mtry = 8`
(`R/sensitivity/mlr3_fitting.R:53`). Only the cells with 6 or more scored
predictors survived, and only barely.

This is not confined to WS1. The production `loco_best_<outcome>` targets
(`_targets.R:1437`) feed a 4-to-7 variable curated set from
`build_best_transportable_predictors()` straight into the same full stack, and
would fail identically. Those targets have never been built in the
`_targets_full` store, which is consistent with the defect.
`R/transportability_best_model.R`'s header claims `run_loco_cv()` "handles small
predictor counts correctly (see mlr3_fitting.R's SMALL_P_SKIP_CORR_THRESHOLD)";
that threshold only skips `washb_prescreen` and `step_corr` and does nothing
about `mtry`.

WS1 guards against it locally in `.nl_fit_score_sl()` by dropping any learner
whose explicit `mtry` exceeds the predictor count, mirroring how `run_loco_cv()`
drops BART (`R/transportability.R:278`). The learner is dropped rather than
clamped, because silently rewriting `mtry` would change what the learner is
while keeping its id, and this scorer exists to be the published estimator. The
dropped learner ids are recorded per row in the `sl_learners_dropped` column.

With the guard in place all 32 cells fit. Exactly one learner, `ranger_low_mtry`,
was dropped in every cell, leaving 8 (source:
`results/tables/corrected/loco_nested_selection.csv`, rows with `scorer = sl`,
columns `sl_learners_dropped`, `sl_learners_used`, `note`; the `note` column is
empty for all 32).

The production fix is deliberately out of scope here, so that production
behaviour was not changed underneath a running comparison. It is filed
separately.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* The point of the workstream. Under `nested_selection`, selection
  receives a world restricted by `.nl_subset()` to the training countries, in
  both the data and the configs, so there is no code path by which the held-out
  country reaches the screen, the dedup, the PCA fit or the stepwise search.
  Verified by reading `nested_select_features()` and by the per-fold selected
  sets differing across folds, which they could not do if selection were
  fold-independent.
- *Selection on evaluation folds.* `original_selection` is retained precisely
  because it selects on the evaluation folds; it is the comparator, labelled as
  such in the `scheme` column, not a result.
- *Identical scoring.* Both schemes pass through the same `.nl_fit_score()` or
  `.nl_fit_score_sl()` on the same pooled dataset and the same test rows. Only
  the variable list differs.
- *Standardization.* Means and standard deviations come from the training fold
  only (`R/corrected/p10_nested_loco.R`, the loop over `scored`), matching
  `loco_glm_auc_single()`.
- *Fold construction.* Leave-one-country-out over four countries. Countries are
  the unit of training, evaluation and inference.
- *Survey weights.* The glm scorer is run twice, equal-row and equal-country.
  Equal-country is the sensitivity the brief asked for, implemented in
  `.nl_country_weights()` as within-training-set inverse country size,
  renormalized to mean 1 so effective sample size is unchanged. It does not
  change the conclusion. The SuperLearner scorer is equal-row only, because
  `mlr3_SL_clustered()` has no weight role at all; that is the finding recorded
  in `docs/PROJECT_STATUS_2026-08_UPDATE.md` section 6, and it is not fixed
  here.
- *Seeds.* `set.seed(params$seed)` is called before each selection and before
  each SuperLearner fit, and the seed is recorded in the `seed` column. The glm
  scorer and the selection procedure contain no stochastic step, so the seed is
  a guard rather than a requirement.
- *Inference.* Four outer folds. Counts, means and ranges of the paired delta
  are reported. No p-value is computed anywhere in this workstream.
- *Not confounded with WS2b.* `prune_gee_covariates()` is called only from
  `R/area_level_comparison.R:379`, `R/transportability_area.R:139` and
  `R/data_prep.R:405`, none of which are on the `build_pooled_dataset()` path
  this workstream uses. No `[gee_hygiene]` or `[hygiene]` line appears in the
  run log. The covariate-hygiene flip therefore does not affect any number
  above.

**Reproducibility reviewer.**

- `targets::tar_manifest()` parses 843 targets, up from 841, and
  `tar_network()` reports 2217 edges, up from 2211. `nested_loco_result` and
  `nested_loco_tables` appear, with `nested_loco_result` depending on the four
  `merged_*` targets and `pipeline_params`, and `nested_loco_tables` depending
  on `nested_loco_result`.
- `ANALYSIS_PROFILE` was added to `get_pipeline_params()` (`R/config.R`) and
  validated against `smoke` and `full`, so it enters the target hash rather than
  changing behaviour invisibly.
- The smoke profile runs and writes to its own `*_SMOKE.csv` filenames, so a
  development run cannot be mistaken for the reportable grid.
- No new file input was introduced, so no stamp target is required. The
  workstream reads only `targets` store objects and writes only new files.
- Paths resolve through `here::here()`. There is no absolute path and no
  `setwd()`.
- `scripts/regression_gate.R` reports all 52 frozen baselines unchanged, so
  nothing under `results/tables/` was regenerated.
