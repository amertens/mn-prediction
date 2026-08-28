# WS7. Anchored transport, new-country intervals, exceedance

WS3 established that the transported map's remaining error is dominated by a
country-level height the covariates cannot supply. A monotone shift that forces
the map's aggregate to a known national value is the direct remedy for exactly
that error, and it leaves district ranking untouched, so it cannot manufacture
spatial signal. This workstream measures how much that is worth, and whether an
anchor of the required quality actually exists.

The two questions are separable and the answers point in opposite directions.

## Result 1: anchoring works, and it works where WS3 said it would

Area-level LOCO, 4 outcomes by 4 held-out countries, 16 folds
(source: `results/tables/corrected/anchored_transport_loco.csv`, grouped by
`anchor_source`):

| Anchor source | Folds | Mean absolute national bias (pp) | Mean MAE (pp) | Mean RMSE (pp) | Mean Spearman |
|---|---|---|---|---|---|
| `unanchored` | 16 | 8.679 | 13.530 | 16.336 | 0.0490 |
| `own_survey` | 16 | 0.000 | 8.876 | 11.411 | 0.0490 |
| `external` | 3 | 44.141 | 44.791 | 46.323 | -0.1032 |

`own_survey` uses the held-out country's own design-based national value. It is
not available in deployment; it measures the value of a perfect anchor. On that
basis anchoring removes 34 percent of the mean absolute error, 13.530 to 8.876
pp, and 30 percent of the RMSE.

Spearman is identical to four decimal places between the unanchored and anchored
arms, which is the intended behaviour: the shift is monotone. The code asserts
this rather than assuming it, and the assertion passed for every fold
(`ranking_preserved` is TRUE throughout).

The gain concentrates exactly where WS3's decomposition said the level offset
was largest (source: same file, per-fold `mae_pp` by `anchor_source`):

| Outcome | Held out | MAE unanchored | MAE anchored | Change |
|---|---|---|---|---|
| child_iron | Gambia | 45.73 | 10.64 | -35.08 |
| women_iron | Gambia | 41.77 | 14.38 | -27.39 |
| women_iron | Ghana | 15.33 | 9.21 | -6.13 |
| child_iron | Sierra Leone | 7.82 | 4.75 | -3.06 |
| child_vitA | Gambia | 12.93 | 11.29 | -1.64 |
| women_vitA | Malawi | 2.12 | 2.56 | +0.44 |

Anchoring lowers MAE in 10 of 16 folds. The six where it does not are all
women's vitamin A or already well-calibrated cells, where the correction is
between +0.08 and +0.44 pp. Gambia's iron cells, which WS3 identified as
carrying a training-to-holdout level ratio of 2.68, account for most of the
total gain. This is the same finding arriving by a second route.

## Result 2: no anchor of the required quality exists

`scripts/build_anchors.R` consolidates every national anchor available for these
countries and outcomes. Coverage of the four primary outcomes
(source: `metadata/anchors/national_anchors.csv`, and the coverage report the
script prints):

**3 of 16 country-outcome cells have an external anchor at all.** There is no
iron entry for any of the four countries. The three that exist are vitamin A,
and each predates its survey by 15 to 19 years.

Anchoring to them is worse than not anchoring
(source: `results/tables/corrected/anchored_transport_loco.csv`, rows with
`anchor_source = external`):

| Outcome | Held out | Anchor value | Year gap | MAE (pp) | National bias (pp) |
|---|---|---|---|---|---|
| child_vitA | Gambia | 0.473 | 19 | 29.064 | +27.521 |
| child_vitA | Malawi | 0.597 | 15 | 50.594 | +50.519 |
| women_vitA | Malawi | 0.557 | 15 | 54.715 | +54.384 |

Mean MAE 44.791 pp against 13.530 unanchored. A stale anchor does not merely
fail to help; it drags a roughly correct map to a badly wrong level, because the
shift is applied with full confidence in a number that is two decades old.

So the method is sound and the input is missing. That is a data-acquisition
problem, not a methods problem. The sources checked were WHO VMNIS (the only
per-nutrient national source in the repository, with no reachable machine
endpoint) and Stevens et al. 2022 (a composite "any micronutrient deficiency"
figure, flagged `anchor_type = "composite"` so it cannot be joined onto a
single-nutrient outcome by accident). GBD was not retrieved: no GHDx credential
is present in `.Renviron`, and the placeholder
`dashboard/data/gbd_estimates.rds` has been left in place rather than deleted,
because deleting a dashboard input is not reversible from here and was flagged
for confirmation.

## Result 3: how wide is a new country?

A country random-intercept model on the Admin-2 mean log biomarker, giving the
interval a genuinely new country's national level would fall in
(source: `results/tables/corrected/new_country_intervals.csv`; the
`*_ratio_scale` columns are `exp()` of the log-scale interval, so they are in
the marker's own units):

| Outcome | Between-country SD (log) | Low | Centre | High | Span |
|---|---|---|---|---|---|
| child_vitA (RBP, umol/L) | 0.1192 | 0.6600 | 1.0088 | 1.5419 | 2.3x |
| women_vitA (RBP, umol/L) | 0.1149 | 0.9603 | 1.4454 | 2.1755 | 2.3x |
| women_iron (ferritin, ug/L) | 0.3997 | 6.1199 | 25.3723 | 105.1905 | 17x |
| child_iron (ferritin, ug/L) | 0.6149 | 2.2277 | 19.8641 | 177.1262 | 80x |

The `t_interval` and `lmer_random_intercept` estimators agree closely and
neither hit a variance boundary, so the width is a real feature of four
countries rather than an artifact of one estimator.

Vitamin A level transports: a new country's mean RBP is predicted within roughly
a factor of two. Iron does not: a new country's mean ferritin is predicted only
within a factor of 80 for children. Any claim that this pipeline can place an
unsurveyed country's iron prevalence on the right level is not supported.

## Result 4: exceedance probabilities

`results/tables/corrected/exceedance_probabilities.csv`, 1,648 rows across 32
country-outcome-threshold combinations, plus
`dashboard/data/exceedance_probabilities.rds` for the dashboard to consume. No
dashboard module was modified. That RDS is not committed: `dashboard/data/` is
gitignored in this repository, as generated dashboard inputs are throughout, and
the file is regenerated by rerunning the script.

Thresholds come from `metadata/who_severity_thresholds.csv`, which carries a
source document and URL per row. **Every row has
`verified_against_source = FALSE`**, and that flag is carried onto every output
row. The values are transcribed from the cited WHO and IZiNCG documents but were
not checked against those documents in this session, and the file records two
substantive caveats: the women's vitamin A categories are applied by analogy to
the preschool-child scale, which WHO does not publish separately, and the 20
percent iron figure is a convention in the literature rather than a WHO
category. Neither should be presented as a WHO threshold without verification.

## What this means together

Anchoring is the highest-value remaining lever on transport accuracy, worth
about a third of the mean absolute error, and it is the intervention WS3's
decomposition predicts. It cannot be used today for the outcome that needs it
most, because no current national iron estimate exists for these countries in
the sources available. The practical consequence is that acquiring a usable
anchor is worth more than any further covariate work.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* The `own_survey` arm uses the held-out country's own national
  value, which is by construction not available in deployment. It is labelled as
  a value-of-anchoring measurement throughout and never presented as achievable
  performance. The `external` arm uses only published values.
- *Monotonicity.* Asserted, not assumed: Spearman is compared across anchor
  sources within each fold and a warning is raised if it moves. It did not.
- *Weights.* Aggregation uses survey n, not population. WorldPop ingestion was
  not run, and a survey-n-weighted aggregate is not the same quantity as a
  design-based national prevalence. Every row records this in
  `aggregation_weight`, and the anchor target is computed with the same weights
  as the aggregate so the two are the same functional.
- *Anchor matching.* Where multiple external anchors exist for a cell the one
  with the smallest year gap is used, and the gap is recorded per row.
- *Interval estimators.* Two are reported because with four countries they can
  disagree; `boundary_fit` flags a collapsed variance component. Neither
  collapsed here.
- *Inference.* Sixteen folds for the own-survey arm, three for the external arm.
  Counts and means only; no p-value.
- *Reused estimator.* The shift delegates to `benchmark_area_predictions()`
  (`R/benchmark_area.R:62`), the project's existing implementation, rather than
  a second one that could drift from it.

**Reproducibility reviewer.**

- No `targets` node added. The workstream is script-driven; wiring it in is
  deferred until an anchor source exists to make it meaningful.
- `metadata/anchors/national_anchors.csv`, `metadata/external_provenance.csv`
  and `metadata/who_severity_thresholds.csv` are all in the tracked `metadata/`
  tree, not under the gitignored `data/`.
- The anchors table is generated by `scripts/build_anchors.R` and records, per
  row, source, URL, licence, year gap and whether it is usable as a
  single-nutrient anchor.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
- All outputs are new files. The regression gate reports the frozen baselines
  unchanged.
