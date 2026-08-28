# WS6. Survey-subsample cost of accuracy

## Why this is a rebuild

`docs/PROJECT_STATUS_2026-08_UPDATE.md` section 5 reports this project's
flagship result. No code or output for it exists anywhere in the repository: the
Stage 0 audit searched for `*subsample*`, `*cost_of_accuracy*`, `*retention*`
and `*sentinel*` and found nothing, and `results/` contains no retention table.
The figures live only as prose, so under this branch's anti-fabrication rule
they could not be cited. The simulation was rebuilt from the specification in
that section before anything was added to it.

Specification as stated there, and as implemented: retention fractions 50 to 90
percent of clusters, 25 replicates per fraction per country-outcome, a
covariate-based area model against "local direct estimate where you have one,
else the national average", both benchmarked against each country's own full
survey, with genuine k-fold cross-validation so a district's model prediction is
never trained on its own subsampled estimate.

Scope: 16 country-outcome cells, 5 retention fractions, 25 replicates
(source: `results/tables/corrected/subsample_cost_of_accuracy.csv`).

## Does it reproduce?

Yes on the result that matters, approximately on the other.

RMSE against the full survey, split by whether a district retained any clusters
(source: `results/tables/corrected/subsample_summary.csv`):

| Stratum | Strategy | Mean RMSE (pp) | Mean MAE (pp) |
|---|---|---|---|
| has clusters | `direct_or_national` | 5.197 | 3.211 |
| has clusters | `spatial_cv` | 10.044 | 7.711 |
| has clusters | `model_cv` | 10.151 | 7.815 |
| zero clusters | `spatial_cv` | 10.513 | 8.496 |
| zero clusters | `model_cv` | 10.943 | 8.928 |
| zero clusters | `direct_or_national` | 11.225 | 9.343 |

Paired within replicate
(source: `results/tables/corrected/subsample_contrasts.csv`):

| Stratum | Comparison | Replicates | Mean delta (pp) | Median delta (pp) | Percent of replicates favouring the first |
|---|---|---|---|---|---|
| has clusters | `model_cv` minus null | 1945 | +4.8922 | +3.9533 | 11.3 |
| has clusters | `spatial_cv` minus null | 2000 | +4.8475 | +4.2597 | 5.2 |
| zero clusters | `model_cv` minus null | 1582 | -0.5955 | -0.1230 | 57.3 |
| zero clusters | `spatial_cv` minus null | 1655 | -0.7113 | -0.1642 | 56.1 |
| zero clusters | `spatial_cv` minus `model_cv` | 1582 | -0.1886 | -0.2382 | 55.8 |

**The zero-coverage figure reproduces.** The note reports the model beating the
national mean by "about 0.7 percentage points" in districts with no surveyed
clusters. The rebuilt covariate model gives **-0.5955 pp** and the new spatial
variant gives **-0.7113 pp**. Same direction, same magnitude.

**The other figure reproduces in direction but not in size.** The note reports
direct estimates beating the model by "+7.7pp RMSE" wherever any clusters
remain. The rebuild gives **+4.8922 pp**. The conclusion is unchanged and if
anything the rebuild is more conservative, but the numbers are not
interchangeable and the difference is not explained. The likely cause is the
covariate set: this rebuild uses the fixed a-priori eight covariates that WS4
and WS5 use, and the original's covariate set is not recorded in the prose. The
original figure remains **not reproducible**; this one is.

## The mean is not the typical case

The mean and median deltas differ by a factor of four in the zero-cluster
stratum: -0.7113 against -0.1642 for the spatial arm, with only 56.1 percent of
replicates favouring it. The average gain is carried by a minority of replicates
with large gains, and in nearly half of them the national mean is as good or
better. Any statement of the form "the model is worth 0.7 pp in unsurveyed
districts" should carry that, because a decision-maker choosing districts to
skip faces the distribution, not its mean.

The gain also depends on outcome and on how much survey remains
(source: `results/tables/corrected/subsample_cost_of_accuracy.csv`, spatial arm
against the null in the zero-cluster stratum):

| Outcome | Mean delta (pp) | Percent favouring |
|---|---|---|
| women_iron | -1.8265 | 63.9 |
| child_iron | -0.9144 | 55.1 |
| child_vitA | -0.3035 | 58.5 |
| women_vitA | +0.1967 | 46.9 |

| Retention | Mean delta (pp) | Percent favouring |
|---|---|---|
| 0.5 | -0.3811 | 53.9 |
| 0.6 | -0.3844 | 51.1 |
| 0.7 | -0.5106 | 55.2 |
| 0.8 | -1.3808 | 62.3 |
| 0.9 | -1.0618 | 59.2 |

Women's vitamin A, the lowest-prevalence outcome, is actively worse under the
model. And the gain grows as more of the survey is retained, which is the
opposite of what a cost-saving argument would want: the model helps most when
you have cut the least.

## Does spatial borrowing beat covariates for skipped districts?

Marginally. `spatial_cv` beats `model_cv` in the zero-cluster stratum by
**-0.1886 pp** mean RMSE, in 55.8 percent of replicates. So the answer to the
question WS6b posed is yes, the zero-coverage gain does exceed the
covariates-only figure, but by about a fifth of a percentage point. Both remain
small.

## Calibration learning curve

How many sentinel districts does a transported map need before it is usable?
The area model is trained on the other countries, calibrated by a logit-scale
location shift fitted on k of the target country's districts, and scored on the
districts not used for calibration
(source: `results/tables/corrected/calibration_learning_curve.csv`; figure at
`results/figures/calibration_learning_curve.png`):

| k sentinel districts | Replicates | Mean RMSE (pp) | Median RMSE (pp) | 2.5th | 97.5th |
|---|---|---|---|---|---|
| 0 | 400 | 16.336 | 14.658 | 1.829 | 47.694 |
| 3 | 400 | 14.196 | 14.353 | 1.983 | 30.104 |
| 5 | 400 | 13.844 | 14.467 | 1.816 | 26.758 |
| 10 | 300 | 14.167 | 15.467 | 4.231 | 24.739 |
| 15 | 300 | 14.258 | 15.277 | 3.790 | 24.270 |
| 20 | 300 | 14.229 | 15.030 | 3.799 | 24.353 |

The mean falls from 16.336 to 13.844 by k=5 and then flattens. The median
barely moves at all, from 14.658 to 14.467, and is higher at k=20 than at k=0.
What actually changes is the upper tail: the 97.5th percentile falls from 47.694
to 26.758.

So sentinel calibration is insurance, not improvement. It truncates the
catastrophic case, which is the fold where the transported map lands at the
wrong level, and does close to nothing for the typical case. That is exactly
what WS3 and WS7 predict: a location shift fixes a location error, and the
median fold's error is not mostly location.

**The question the plan asked was at which k the transported-plus-calibrated
error matches within-country model error. No k in the grid achieves it.** The
best calibrated transported RMSE is 13.844 pp at k=5, while the within-country
cluster spatial model reaches 12.663 pp
(source: `results/tables/corrected/cluster_mbg_within_country.csv`, arm
`spatial_only`, signal cells). Adding sentinel districts up to 20 does not close
that gap, and beyond k=5 it does not help at all.

## What this workstream supports

The August note's framing survives, and the rebuild sharpens it. Direct
estimation dominates wherever any survey remains, by roughly 5 pp of RMSE and in
89 percent of replicates. The model's value is confined to districts you
deliberately do not sample, it is worth a few tenths of a percentage point on
average, it is absent in nearly half of replicates, and it is negative for the
lowest-prevalence outcome. Spatial borrowing is the better of the two model
options for that narrow use by about 0.19 pp.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* The correction the original note describes is implemented and
  checked. `.ss_model_cv()` assigns district-level folds and a district's
  prediction comes from a fit that excluded it. Zero-cluster districts are
  absent from training entirely, so their prediction from the full fit is
  out-of-sample by construction. `.ss_spatial_cv()` does the same at cluster
  level, holding out every cluster in a fold's districts together.
- *Evaluation target.* All strategies are scored against the FULL survey
  district prevalence, never against the subsampled estimate, which is what
  makes the comparison meaningful.
- *Stratification.* Results are split by whether a district retained clusters
  and never averaged across the two strata. The strata answer different
  questions and pooling them hides the result.
- *Thin districts.* Districts whose own full-survey estimate is below the same
  `DIST_MIN_NSVY` threshold used elsewhere are excluded from the evaluation, so
  the target is not itself noise.
- *Mean against median.* Both are reported for every contrast, along with the
  percent of replicates favouring each side, because the means are skewed and
  the mean alone would overstate the model's value.
- *Replicate structure.* 5 retention fractions by 25 replicates by 16 cells,
  recorded per row. Counts differ slightly between strategies because a fit can
  fail for a replicate; the row counts are reported rather than padded.
- *Seeds.* Seeded from `params$seed` and recorded per row.
- *Inference.* Means, medians, and percent-favouring. No p-value: replicates
  within a cell share the same underlying survey and are not independent, so a
  p-value computed across them would be badly overstated. This is a deliberate
  departure from the original note, which quotes p < 1e-300 and p = 1.2e-20 for
  these contrasts.

**Reproducibility reviewer.**

- No `targets` node added; the simulation is script-driven and expensive, so
  wiring it into the default pipeline would make every build pay for it.
- Smoke profile runs on 2 countries, 1 outcome and 5 replicates and writes to
  `*_SMOKE.csv`.
- `tar_manifest()` parses 845 targets and `tar_network()` 2218 edges, unchanged,
  confirming the new `R/corrected/p14_subsample.R` does not break sourcing.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
- All outputs are new files. The regression gate reports the frozen baselines
  unchanged.
