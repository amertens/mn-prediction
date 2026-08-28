# WS3. Measurement and label harmonization

Two questions, both about whether the dominant leave-one-country-out failure
mode is epidemiology or measurement artifact:

1. Does putting every country on one inflammation-adjustment protocol shrink the
   national level bias?
2. On the continuous biomarker scale, what share of each fold's national offset
   is a pure country intercept that no covariate could have supplied?

## A Stage 0 correction

The Stage 0 audit recorded `brinda_country_method()` as referenced in a comment
but never defined, and the plan therefore listed writing it as a WS3 task. That
was wrong: the function exists at `R/brinda_adjustment.R:137`. The audit matched
only the comment at line 38 and stopped. Nothing else in the plan depended on
it, and the function is used here as it stands.

## What each country-outcome actually receives

`adjustment_inventory()` (`R/brinda_adjustment.R`) reads the configs rather than
asserting from memory, so the table cannot drift from the code.
(source: `metadata/adjustment_inventory.csv`.)

**Vitamin A is already harmonized.** `apply_brinda_vita_binary()` overwrites
every country's VAD binary at runtime with one BRINDA CRP+AGP definition, and
`brinda_country_method()` returns `BRINDA CRP+AGP` for all four active
countries. The configured binaries in `R/config.R` (`gw_cVAD_Thurn` and so on)
are never the ones used. A uniform-adjustment sensitivity for vitamin A would be
a no-op, which is why the sensitivity below is about iron.

**Iron is not.** `UNIFORM_TRANSPORT_TAGS` (`R/admin2_analysis.R:20`) applies a
uniform WHO cut-point to iron, but it applies it to each country's own
configured continuous column, and those are adjusted four different ways
(source: `metadata/adjustment_inventory.csv`, column
`configured_continuous`):

| Country | Configured continuous | Adjustment |
|---|---|---|
| Gambia | `gw_LogFerAdj` | survey-agency adjustment, log scale |
| Ghana | `gw_cFerrAdjThurn` / `gw_wFerrAdjThurn` | Thurnham |
| Sierra Leone | `gw_cFerrAdj` (children) / `gw_wFerAdjBR1` (women) | two different survey-provided adjustments within one country |
| Malawi | `sf_reg` | regression-adjusted serum ferritin |

A uniform cut-point applied to non-uniform adjustments is not a uniform outcome.

## Result 1: uniform adjustment moves the pattern, not the level

`brinda_adjust_ferritin()` (`R/brinda_adjustment.R`) implements the same
regression correction as the existing RBP function with the sign reversed:
inflammation raises ferritin rather than depressing it, so the coefficients are
clamped to be non-negative and the correction lowers the adjusted value.
Clamping in the wrong direction would let collinearity manufacture or erase
deficiency, so the two clamps are not interchangeable.

Both arms are rebuilt through `build_outcome_dataset()` and scored through
`compute_svy_admin2()`, so the survey design, weighting and domain handling are
identical and only the biomarker adjustment differs.

Area-level LOCO, 4 held-out countries per outcome
(source: `results/tables/corrected/uniform_brinda_loco_level_bias.csv`, grouped
by `outcome` and `scheme`):

| Outcome | Scheme | Mean absolute national bias (pp) | Mean RMSE (pp) | Mean Pearson r |
|---|---|---|---|---|
| child_iron | `configured` | 15.087 | 23.727 | 0.1988 |
| child_iron | `uniform_brinda` | 14.290 | 21.975 | 0.3150 |
| women_iron | `configured` | 16.080 | 25.180 | 0.0257 |
| women_iron | `uniform_brinda` | 14.572 | 25.610 | 0.1748 |

The level bias barely moves: 0.797 pp for children and 1.508 pp for women, on
biases of 15 to 16 pp. The correlation moves a great deal: 0.1988 to 0.3150 for
children and 0.0257 to 0.1748 for women.

The same split shows in the biomarker medians
(source: `results/tables/corrected/uniform_brinda_prevalence.csv`, columns
`raw_ferritin_median` and `adj_ferritin_median`):

| Outcome | Raw ferritin median, lowest to highest country | Ratio | After uniform adjustment | Ratio |
|---|---|---|---|---|
| child_iron | 11.33 to 71.35 | 6.30 | 7.97 to 38.65 | 4.85 |
| women_iron | 23.46 to 47.86 | 2.04 | 16.57 to 35.93 | 2.17 |

Uniform adjustment removes part of the children's spread, from a factor of 6.30
to 4.85, and none of the women's. Sierra Leone's children carry a median CRP of
3.98 mg/L against 0.49 to 1.47 elsewhere (same source, column `crp_median`),
which is why the children's correction is the larger one.

**Answer to question 1: no.** Harmonizing the adjustment does not shrink the
national level bias to any useful degree. It does improve the transported
spatial pattern enough to be worth adopting on its own merits, and it is
recorded as the `uniform_brinda` scheme rather than made the default.

## Result 2: the offset is mostly a pure intercept

`R/corrected/p11_level_decomposition.R` works at Admin-2 level on the log
harmonized biomarker. Each survey is standardized within itself, a model is
fitted on the training countries, and the transported prediction is
reconstituted with the training-pooled location and scale, because a deployment
in an unsurveyed country has no other option. The national bias then splits
exactly, by algebra rather than by fitting, into

- `location_offset` = training mean minus held-out mean, the pure intercept, and
- `pattern_term` = whatever net level the covariates themselves assert.

(source: `results/tables/corrected/level_bias_decomposition.csv`, summarised by
`outcome`):

| Outcome | Median absolute location term | Median absolute pattern term | Median bounded location share | Folds where location is larger | Level ratio, min to max |
|---|---|---|---|---|---|
| child_iron | 0.5350 | 0.4271 | 0.449 | 1 of 4 | 0.526 to 2.679 |
| women_iron | 0.2193 | 0.0826 | 0.689 | 2 of 4 | 0.790 to 2.130 |
| child_vitA | 0.1539 | 0.0508 | 0.745 | 3 of 4 | 0.808 to 1.203 |
| women_vitA | 0.0740 | 0.0763 | 0.507 | 2 of 4 | 0.901 to 1.265 |

The terms are on the log scale, so `level_ratio_train_over_holdout` is their
multiplicative reading. For child iron the training pool sits anywhere from
0.526 to 2.679 times the held-out country's level, a swing of a factor of five
across four folds. For child vitamin A the same range is 0.808 to 1.203.

`location_share_bounded` is reported rather than the raw share because the raw
share is unbounded: when the covariate term opposes the intercept the two partly
cancel and the share exceeds 1, reaching 47.4 in one fold where they nearly
cancel exactly. That value is informative about that fold and useless as a
summary, so both are written to the CSV and only the bounded one is averaged.

**Answer to question 2: the intercept is the larger term for three of the four
outcomes**, and for child iron it is comparable to the pattern term while
opposing it in three of four folds. The transported model's remaining error is
dominated by a country-level height that the covariates cannot supply, which is
an argument for anchoring (WS7) rather than for more predictors.

## A masked failure found while doing this

The first decomposition run returned `n_selected = 0` and `pattern_term = 0` in
all 16 folds. Read at face value that says the covariates contribute nothing,
which would have been a publishable-looking negative result. It was wrong.

`.ld_standardize()` built the standardized response through `mu[country]`, which
`tapply()` returns as a NAMED vector. A named response reaches
`glmnet::cv.glmnet()` and fails with "non-conformable arrays".
`.tr_fit_predict()` (`R/transportability_area.R:170`) wraps that call in
`tryCatch(..., error = function(e) NULL)` and on failure returns the training
mean with an empty variable list. Because every country's standardized response
is centred within survey, the pooled training mean is exactly 0, so the fallback
produced exactly the signature a genuine null result would produce.

Stripping the names fixes it, and the fits then select 6 or more predictors with
out-of-sample correlations up to 0.44. The fix and the reason are recorded in
`R/corrected/p11_level_decomposition.R`.

The wider point is about `.tr_fit_predict()`, not about this file: its silent
fallback converts any fit failure into a plausible-looking null result anywhere
in the production area-transport path. That is worth a separate look and is not
changed here.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* The uniform adjustment is fitted within each survey on that
  survey's own rows and involves no outcome and no held-out country. The
  decomposition standardizes within survey before any fold is formed, and the
  transported reconstruction deliberately uses the TRAINING location and scale,
  never the held-out country's, which is the whole point of the exercise.
- *Selection on evaluation folds.* `.tr_screen()` runs inside each fold on the
  training rows only.
- *Identical treatment of the two arms.* Both are rebuilt through
  `build_outcome_dataset()` and aggregated through `compute_svy_admin2()`. The
  stored `outcome_data_*` targets were deliberately not reused: they were built
  by an older code version and no longer agree with the current population
  filter on row count, so using one stored arm against one fresh arm would have
  compared two code versions rather than two adjustments.
- *Survey weights.* The Admin-2 aggregation is design-based through
  `compute_svy_admin2()` and identical in both arms. The decomposition's area
  means use a plain weighted mean because only point estimates enter it and no
  standard error is reported from it; that is stated in the script header.
- *Denominators.* The two arms can differ in usable rows, because a row can have
  a configured adjusted value and no raw marker or the reverse. Both row counts
  are recorded per cell (`n_rows_configured`, `n_rows_uniform`).
- *Seeds.* The only stochastic step is `cv.glmnet` fold assignment, seeded from
  `AREA_TRANSPORT_RECIPE$seed`, identically in both arms.
- *Inference.* Four folds per outcome. Means and ranges only; no p-value.
- *Sign convention.* The ferritin clamp is non-negative and the RBP clamp is
  non-positive, checked against the direction each marker moves under
  inflammation.

**Reproducibility reviewer.**

- `tar_manifest()` parses 845 targets, up from 843, and `tar_network()` reports
  2218 edges, up from 2217. `assay_lineage_stamp` and `assay_lineage_check`
  appear, with the check depending on the stamp.
- `metadata/assay_lineage.csv` is an untracked file input and is guarded by
  the md5-plus-always-cue stamp pattern copied from `gee_parity_stamp_*`.
- The lineage table is generated by `scripts/build_assay_lineage.R`, not
  hand-written, so it cannot drift from the code it describes. 24 rows for 24
  configured country-outcomes; the validator passes.
- Specimen, instrument and laboratory are UNKNOWN in all 24 rows. Nothing in the
  repository records them for any of the four surveys, and the per-cell source
  column says where a reader would have to go. They are not guessed.
- The validator lives in `R/assay_lineage.R` rather than in the script, because
  a pipeline target calls it and `tar_source()` loads `R/` but not `scripts/`.
- All outputs are new files. `scripts/regression_gate.R` reports the frozen
  baselines unchanged.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
