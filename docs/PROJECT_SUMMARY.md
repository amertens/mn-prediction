# Micronutrient prediction project: status summary

Standing overview of the project. Updated 2026-08-28 to fold in the
`analysis-updates-2026-08` branch.

Detail and full source maps live in `docs/findings/` (one file per workstream)
and `docs/PROJECT_STATUS_2026-09_UPDATE.md`. Earlier context is in
`PROJECT_STATUS_2026-06.md` and `PROJECT_STATUS_2026-08_UPDATE.md`.

Every number below is read from a produced CSV. Where a figure quoted in an
earlier version of this document could not be reproduced from code, that is
stated rather than carried forward.

## What we have done

- Built a harmonized 4-country pipeline (Gambia, Ghana, Sierra Leone, Malawi;
  Tanzania archived pending valid biomarker data) spanning individual-level
  survey outcomes, roughly 600 pooled remote-sensing and socioeconomic proxy
  predictors, and both individual-level and Admin-2 area-level modelling tracks.
- Fixed data-infrastructure bugs that were silently zeroing predictor coverage
  for specific countries: WFP food-price admin-name matching (case sensitivity
  for Ghana, wrong admin level for Malawi), a SoilGrids cache-filename bug, a
  Sierra Leone external-predictor filename mismatch, and a `haven_labelled`
  package-load-order crash in outcome recoding.
- Unlocked and harmonized new predictor domains across all four countries: a
  scale-free WFP food-price deviation index, 105 harmonized DHS Admin-2
  covariates, nighttime lights, travel-time and market access, ESPEN helminths,
  MapSPAM crop mix, and IHME modelled disease burden.
- Fixed two SuperLearner reliability bugs: unseeded internal ensemble-weight
  cross-validation, and correlation-based prescreening being misapplied to small
  curated predictor sets where it destroyed signal.
- Compared model classes head to head on identical covariates: regression, full
  SuperLearner ensembles, XGBoost, and random forest, both within country and
  under leave-one-country-out (LOCO) cross-validation.
- Swept predictor-pool size from 150 down to 10 to find where cross-country
  transportability peaks.
- Built the Admin-2 small-area-estimation track with a reliability ceiling and a
  national-benchmark correction.
- Tested rank-based objectives and domain-generalization methods (IRM, Group DRO)
  for few-environment transfer.
- Investigated survey timing and seasonality as a possible confound.
- Ran an external-style critical-review audit of survey weighting, proxy linkage,
  external-data temporality and the ensemble setup.

Added on the `analysis-updates-2026-08` branch:

- Froze 52 headline result tables with an md5 manifest, and built a regression
  gate that classifies every later table as unchanged, new scheme rows, or a
  changed baseline. It has reported the baselines unchanged after every
  workstream.
- Implemented nested LOCO, so predictor selection is a function of the training
  countries only, and reported it against the original scheme with the same
  estimator.
- Implemented a uniform BRINDA ferritin adjustment and a decomposition of each
  fold's national level bias into a country intercept and a covariate-asserted
  component.
- Extended the threshold-integrating distributional estimator to every outcome
  with a continuous biomarker and scored it against the binary path.
- Built a point-referenced cluster-level spatial model on the existing GPS
  coordinates, with no Earth Engine extraction.
- Rebuilt the survey-subsample cost-of-accuracy simulation from its prose
  specification, because no code or output for it existed, then added a spatial
  variant and a sentinel-district calibration learning curve.
- Implemented national-anchor recalibration, new-country prediction intervals,
  and district exceedance probabilities.
- Migrated every Admin-2 join from the name alone to the (Admin1, Admin2) pair,
  and added a regression test asserting unit counts per country.
- Generated an assay-lineage table with a pipeline gate that fails when a
  modelled outcome has no provenance row.

## What we have found

### Cross-country transport is weaker than previously reported

- **The transport signal is largely a selection artifact.** The previously
  reported mean LOCO AUC of about 0.57 came from a predictor set chosen by
  scoring candidates on the same folds it was then reported on. Rerunning the
  entire selection procedure inside each outer fold, and scoring both schemes
  with the same estimator over 16 folds:

  | Scorer | Original selection | Nested selection | Delta |
  |---|---|---|---|
  | glm | 0.6199 | 0.5380 | -0.0818 |
  | SuperLearner | 0.5977 | 0.5214 | -0.0763 |

  Nested selection is better in 1 of 16 folds under either scorer, and 0 of 16
  under equal-country weighting. Four of 16 nested SuperLearner folds fall below
  0.5. The aggregate confidence interval quoted previously, [0.544, 0.600]
  excluding chance, was computed on the contaminated set and does not support the
  exclusion-of-chance claim. The nested equivalent of that interval is not yet
  computed.
  (source: `results/tables/corrected/loco_nested_selection_summary.csv`)

- **The same correction applies to the eight-SoilGrids spatial model**, and the
  answer was already in the repository, unpropagated. Mean Pearson r over the
  same 16 cells: `locked_published` 0.3324, `selected_in_fold` 0.1426,
  `spatial_only` with no soil features at all 0.2861. Honest in-fold selection
  does worse than using no soil features.
  (source: `results/tables/frozen_2026-08/sandbox_parsimony_out/spatial_plus_soil_rescored.csv`)

### The residual error is a country intercept, and anchoring is the remedy

- **Harmonizing the inflammation adjustment does not fix the level offset.**
  Vitamin A was already harmonized; iron was not, with four different
  survey-provided adjustments across four countries. Putting iron on one uniform
  BRINDA protocol moves the mean absolute national bias from 14.637 to 15.637 pp
  for children and from 17.325 to 15.477 pp for women, while raw child ferritin
  medians span a factor of 6.30 across countries and still span 4.85 after
  adjustment. It does improve the transported spatial pattern substantially
  (child iron Pearson 0.1808 to 0.3097; women iron -0.1865 to 0.0207), so it is
  worth adopting on those grounds, and it is recorded as a scheme rather than
  made the default.
  (source: `results/tables/corrected/uniform_brinda_loco_level_bias.csv`)

- **The decomposition says why.** On the log biomarker scale the pure country
  intercept is the larger term for three of four outcomes, and for child iron
  the training pool sits between 0.523 and 2.669 times the held-out country's
  level.
  (source: `results/tables/corrected/level_bias_decomposition.csv`)

- **Anchoring works, and it is the highest-value remaining lever.** A monotone
  logit-scale shift forcing the map's aggregate to a national value, over 16
  folds:

  | Anchor | Absolute national bias | MAE | RMSE |
  |---|---|---|---|
  | unanchored | 9.135 pp | 13.671 pp | 16.476 pp |
  | own survey | 0.000 pp | 8.934 pp | 11.499 pp |
  | external | 44.141 pp | 44.690 pp | 46.386 pp |

  A perfect anchor removes 35 percent of the mean absolute error and improves 9
  of 16 folds, with district ranking untouched.
  (source: `results/tables/corrected/anchored_transport_loco.csv`)

- **But no anchor of the required quality exists.** Only 3 of 16 primary
  country-outcome cells have any external anchor, there is no iron entry for any
  of the four countries, and the three vitamin A anchors predate their surveys
  by 15 to 19 years. Anchoring to them is far worse than not anchoring, because
  a stale number is applied with full confidence. The method is sound and the
  input is missing, which is a data-acquisition problem rather than a methods
  problem.
  (source: `metadata/anchors/national_anchors.csv`)

- **A new country's level is not predictable for iron.** A country
  random-intercept model gives a new country's geometric-mean RBP within roughly
  a factor of 2.3, and its mean child ferritin only within a factor of 79.
  (source: `results/tables/corrected/new_country_intervals.csv`)

### Within-country prediction

- **Admin-2 is finer than these surveys can resolve.** The median district
  contributes few biomarker measurements, the reliability ceiling is low, and
  the ceiling falls further at cluster level (Admin-1 0.59, Admin-2 0.31,
  cluster 0.22). This is the structural constraint behind every weak result
  here: the effective sample size for cross-country learning is the number of
  areas, not the number of individuals.

- **Covariates hurt at cluster level.** A point-referenced spatial model on the
  existing GPS coordinates finds both covariate arms worse than doing nothing,
  and adding covariates to the spatial smoother makes it worse still. Spatial
  borrowing beats the national-mean null in about half of cells by a fraction of
  a percentage point. *Recomputation after the join-key migration is in progress;
  the pre-migration run gave spatial-only MAE 10.087 pp against a null of 10.599
  pp, better in 13 of 24 cells.*
  (source: `results/tables/corrected/cluster_mbg_within_country.csv`)

- **The distributional estimator should not be the default.** Over the 23 cells
  with a reliability ceiling above zero it improves area correlation (-0.0592 to
  -0.0326) and worsens discrimination (AUC 0.5944 to 0.5887), Brier skill, MAE
  and the correlation read against the ceiling. Where it does help is where
  theory predicts, in the tail: the largest gains are all cells under 10 percent
  prevalence. It is kept as a declared alternative.
  (source: `results/tables/corrected/distributional_paired_signal.csv`)

### The survey-subsample result

- **The flagship result had no code.** It existed only as prose, so its figures
  could not be cited and the simulation was rebuilt from its specification. The
  zero-coverage figure reproduces: the model beats the national mean by about
  0.7 pp, against the note's "about 0.7 percentage points". The other figure
  does not: the rebuild gives roughly +4.9 pp for direct-beats-model where
  clusters remain, against the note's "+7.7pp". Same direction and order, not
  the same number, and the difference is unexplained because the prose does not
  record the original covariate set. *Recomputation after the join-key migration
  is in progress; figures here are from the pre-migration rebuild.*
- **The mean is not the typical case.** The zero-cluster gain has a median about
  a quarter of its mean, and only about 56 percent of replicates favour the
  model at all. The gain also grows as more of the survey is retained, which is
  the opposite of what a cost-saving argument wants, and it is negative for the
  lowest-prevalence outcome.
- **No p-value is reported for these contrasts.** Replicates within a cell share
  one underlying survey and are not independent, so the p-values quoted in the
  August note for these comparisons are not supportable.
- **Sentinel calibration is insurance, not improvement.** Calibrating a
  transported map on k of the target country's districts drops the mean held-out
  error and barely moves the median, while cutting the upper tail sharply. No k
  in the grid reaches within-country model error.

### Unchanged from earlier work

These were established before this branch and nothing here contradicts them.

- Complex ensembles beat simple models within a country; simple, heavily
  regularized models transport across countries at least as well as complex ones.
- More predictors hurt cross-country transport. The newly unlocked domains never
  displaced the small climate, soil and malaria set; the best DHS candidate
  ranked 51st of about 590.
- Reframing the evaluation from classification to district ranking does not
  surface hidden signal.
- Seasonality and interview timing do not confound the results, and this is a
  closed question.
- The individual-level SuperLearner was fit unweighted; a weighted refit did not
  change the headline (0.558 against 0.572). That comparison exists only as prose
  and has no committed table, so its figures cannot be cited from code.

## Defects found and fixed

- **Admin-2 join key.** Malawi's GADM layer has four genuine same-named district
  pairs in different regions. Joining on the name alone caused three separate
  failures: covariates averaged across regions, the pooled area table inflated
  from 87 to 90 Malawi districts, and cluster tables inflated from 103 to 107.
  Affected districts carried double weight in every area-level fit and one copy
  of each carried another region's coordinates. Now keyed on (Admin1, Admin2),
  with a regression test and a multiplication warning at every join.
- **`ranger_low_mtry` hard-codes `mtry = 8`**, which aborts the entire
  SuperLearner fit for any curated set with fewer predictors. This cost 28 of 32
  cells on the first nested-LOCO attempt and also affects the production
  `loco_best_*` targets, which have never been built.
- **`.tr_fit_predict()` silently converts a fit failure into the training mean.**
  On a within-survey-centred response that fallback is exactly zero, so a failed
  fit is indistinguishable from a genuine null result. It produced a
  plausible-looking false negative before being caught.
- **`compute_svy_admin2()` filter-before-design ordering** was fixed in August
  2026; the audit item is closed.

Still open from the earlier audit: the 2 km GPS buffer for the cluster-level
comparator undercorrects for DHS displacement of up to 5 to 10 km, and the
food-security domain accepts Admin-2 name matches at a 30 percent rate.

## Plan

- **Reframe the flagship deliverable** from predicting deficiency in a country
  with no data, which the evidence does not support, to helping a country decide
  which districts to skip sampling. This survives the rebuild, with the sharper
  statement that the model's value is confined to districts deliberately not
  sampled, is worth a few tenths of a percentage point on average, and is absent
  in nearly half of replicates.
- **Acquiring a usable national anchor is now the single highest-value action.**
  Anchoring is worth roughly a third of the transport error, and no current
  national iron estimate exists for these countries in any source available. This
  is cheaper than any further covariate work and worth more.
- **Run one or two national case studies** built around the corrected subsample
  result, with the narrow claim rather than the broad one.
- **Bring in outside costing expertise** to convert the zero-coverage edge into a
  dollar framing.
- **Scope any GBD contribution narrowly**, as within-country-year disaggregation
  constrained to reproduce a national total, and state that remote-sensing
  predictors cannot reliably extend back to 1990.
- **Keep pursuing survey data acquisition** without gating the strategy pivot on
  it.
- **Correct the manuscript** before submission. Five passages are affected by the
  nested-selection result, itemised in `docs/findings/WS1_nested_loco.md`. The
  manuscript has not been edited by this work.
- **Finish the two remaining audit items** above before anything goes external.
- **Verify the WHO and IZiNCG severity thresholds** in
  `metadata/who_severity_thresholds.csv` against their source documents. Every
  row is currently marked `verified_against_source = FALSE`, and two carry
  substantive caveats.

## Methods compared that did not improve estimation

- **XGBoost**, shallow and deep, and **random forest**: no LOCO advantage over a
  well-selected simple regression.
- **Full SuperLearner ensembles against a single linear model** on identical
  unscreened covariates: no advantage across countries.
- **Larger candidate pools** of 100 to 150 variables: degraded transport relative
  to the parsimonious set.
- **Newly harmonized predictor domains**, forced in and left to compete: none
  displaced the original small set.
- **Construct-based PCA reduction**: mostly no value beyond existing variables,
  though one vegetation component earned a place.
- **Rank-based objectives**: no method significantly beat another, and the
  comparison was underpowered.
- **Domain generalization (IRM, Group DRO)**: neither improved on naive pooled
  regression with three training countries per fold; Group DRO significantly
  underperformed.
- **The distributional estimator as a default** (new): improves area correlation
  and calibration, worsens discrimination, Brier skill and MAE. Retained as an
  alternative for tail-prevalence outcomes.
- **Cluster-level covariates** (new): both covariate arms are worse than the
  national-mean null, and adding covariates to the spatial smoother makes it
  worse.
- **Cluster-level spatial resolution as such** (new): the reliability ceiling
  falls going finer than Admin-2, so the Earth Engine cluster-covariate
  extraction that would have supported it was dropped before it was run.

## Provenance of the figures above

Current as of the join-key migration rebuild: nested LOCO, covariate hygiene,
uniform BRINDA, level decomposition, distributional estimator, anchored
transport, new-country intervals.

Recomputation in progress at the time of writing: the cluster spatial model and
the subsample simulation. Their figures above are from the pre-migration run and
are labelled as such; only Malawi rows are expected to move.

Not reproducible from code, and therefore not citable: the August note's "+7.7pp"
subsample figure and its weighted-versus-unweighted SuperLearner comparison.
