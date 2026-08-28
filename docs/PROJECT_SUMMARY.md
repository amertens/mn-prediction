# Micronutrient prediction project: status summary

Standing overview of the analysis, its findings, and where it can go next.
Per-workstream detail with full source maps is in `docs/findings/`.

Every number below is read from a produced CSV and cited with its file.

## What the analysis does

- **Scope.** Four countries with biomarker microdata: Gambia, Ghana, Sierra
  Leone, Malawi. Outcomes are biomarker-defined deficiency for vitamin A, iron,
  folate, B12 and zinc, in preschool children and women of reproductive age.
- **Predictors.** Roughly 600 pooled proxy variables that require no survey:
  remote-sensing (SoilGrids, iSDAsoil, climate, vegetation, nighttime lights),
  Malaria Atlas, IHME modelled disease burden, harmonized DHS Admin-2
  indicators, food prices, travel time and market access, helminths and crop mix.
- **Two modelling tracks.** An individual-level SuperLearner ensemble aggregated
  to Admin-2, and an Admin-2 area-level model. Both are evaluated within country
  and under leave-one-country-out (LOCO) transport.
- **A leakage-free evaluation layer.** Honest spatial cross-validation,
  out-of-fold calibration, design-aware Admin-2 error metrics, split-conformal
  and design-based intervals, out-of-support trust flags, and a reliability
  ceiling that says how much correlation sampling noise permits.
- **Selection nested inside the evaluation folds.** Predictor selection is a
  function of the training countries only, so a reported transport number is
  out-of-sample with respect to the selection as well as the data.
- **Harmonized outcomes.** One uniform WHO cut-point per outcome, and one
  inflammation-adjustment protocol per biomarker, so a cross-country comparison
  is not comparing definitions.
- **Reproducibility.** Frozen baseline tables with an md5 manifest, a regression
  gate that classifies every table as unchanged, new scheme rows or a changed
  baseline, an assay-lineage gate that fails when a modelled outcome has no
  recorded measurement provenance, and md5 stamps on untracked file inputs.

## Findings

### Cross-country transport is close to chance

Predictor selection, rerun inside each outer fold so it sees only the training
countries, over 16 folds:

| Scorer | Mean LOCO AUC |
|---|---|
| SuperLearner | 0.5214 |
| glm | 0.5380 |

Four of 16 SuperLearner folds fall below 0.5. The eight-SoilGrids spatial model
behaves the same way: selecting its features inside the fold gives mean Pearson
r of 0.1426, while a covariate-free spatial smoother on the same cells reaches
0.2861. Honest selection of soil covariates performs worse than using none.

The practical reading is that this pipeline cannot place an unsurveyed country's
subnational map on the right footing from proxies alone.
(source: `results/tables/corrected/loco_nested_selection_summary.csv`,
`results/tables/frozen_2026-08/sandbox_parsimony_out/spatial_plus_soil_rescored.csv`)

### The residual error is a country intercept, not a covariate failure

On the log biomarker scale, each fold's national bias splits by algebra into a
pure country location offset and a covariate-asserted component. The intercept
is the larger term for three of four outcomes, and for child iron the training
pool sits between 0.523 and 2.669 times the held-out country's level.

Raw child ferritin medians span a factor of 6.30 across the four countries. One
uniform BRINDA CRP+AGP protocol still leaves a factor of 4.85, and moves the
mean absolute national bias by roughly a percentage point on a 15 pp bias. It
does improve the transported spatial pattern substantially (child iron Pearson
0.1808 to 0.3097, women iron -0.1865 to 0.0207), which is why it is available as
a scheme.

Measurement harmonization is therefore worth doing for pattern and does not
solve the level.
(source: `results/tables/corrected/level_bias_decomposition.csv`,
`results/tables/corrected/uniform_brinda_loco_level_bias.csv`,
`results/tables/corrected/uniform_brinda_prevalence.csv`)

### Anchoring to a national estimate is the effective remedy

A national anchor is a known country-level prevalence figure used to fix the
height of a transported map. A single constant shift on the logit scale makes
the map's weighted aggregate equal the anchor. Because the shift is monotone,
district ranking is untouched, so it corrects level without manufacturing
spatial signal. It targets exactly the error term above.

| Anchor | Absolute national bias | MAE | RMSE |
|---|---|---|---|
| none | 9.135 pp | 13.671 pp | 16.476 pp |
| the country's own survey | 0.000 pp | 8.934 pp | 11.499 pp |
| a published external estimate | 44.141 pp | 44.690 pp | 46.386 pp |

A perfect anchor removes 35 percent of the mean absolute error and improves 9 of
16 folds. The country's own national value is not available when deploying to an
unsurveyed country; it measures what a perfect anchor is worth.

The external row is the binding constraint. Only 3 of 16 primary
country-outcome cells have any external anchor, there is no iron entry for any
of the four countries, and the three vitamin A anchors predate their surveys by
15 to 19 years. Anchoring to a stale figure is far worse than not anchoring,
because the shift is applied with full confidence in it.
(source: `results/tables/corrected/anchored_transport_loco.csv`,
`metadata/anchors/national_anchors.csv`)

### A new country's level is predictable for vitamin A and not for iron

A country random-intercept model on the Admin-2 mean log biomarker gives the
range a genuinely new country's national level would fall in:

| Outcome | Interval, marker units | Span |
|---|---|---|
| child vitamin A (RBP, umol/L) | 0.662 to 1.536 | 2.3x |
| women vitamin A (RBP, umol/L) | 0.961 to 2.174 | 2.3x |
| women iron (ferritin, ug/L) | 6.12 to 105.02 | 17x |
| child iron (ferritin, ug/L) | 2.23 to 176.40 | 79x |

With four countries these intervals are wide, and that width is the honest
product. A factor of 79 for child ferritin means no claim about an unsurveyed
country's iron level is supportable.
(source: `results/tables/corrected/new_country_intervals.csv`)

### Admin-2 is finer than these surveys can resolve

The reliability ceiling, the correlation attainable given sampling noise alone,
falls as the spatial unit gets finer: 0.59 at Admin-1, 0.31 at Admin-2, 0.22 at
cluster level. The median district contributes few biomarker measurements, and
in three of four countries the median district holds a single survey cluster.

This is the structural constraint behind the weak results. The effective sample
size for cross-country learning is the number of areas, between 14 and 89 per
country, not the number of individuals.

### Covariates do not help at cluster level

A point-referenced spatial model on the survey GPS coordinates, evaluated
leave-one-district-out so the smoother must extrapolate to an unsurveyed
location, finds both covariate arms worse than giving a district the national
average, and finds that adding covariates to the spatial smoother makes it
worse than the smoother alone. Spatial borrowing beats the national-mean null in
about half of cells, by a fraction of a percentage point of MAE.
(source: `results/tables/corrected/cluster_mbg_within_country.csv`)

### The distributional estimator helps only in the tail

Modelling the continuous biomarker and integrating past the cut-point, rather
than dichotomising first, improves area-level correlation (-0.0592 to -0.0326)
and worsens discrimination (AUC 0.5944 to 0.5887), Brier skill, MAE and the
correlation read against the ceiling, over the 23 cells with a ceiling above
zero. Its gains concentrate exactly where the efficiency argument predicts, in
outcomes whose cut-point sits in a tail: every one of the largest gains is a
cell under 10 percent prevalence.

It is available as a declared alternative rather than the default.
(source: `results/tables/corrected/distributional_paired_signal.csv`)

### Modelling helps only for districts a survey skips entirely

Simulating reduced survey budgets, 50 to 90 percent of clusters retained with 25
replicates per fraction, and scoring against each country's full survey with
district-level cross-validation:

- For any district that keeps even a few clusters, the direct survey estimate
  beats the model by roughly 5 percentage points of RMSE, and the model wins in
  about a tenth of replicates.
- For districts with zero retained clusters, the model beats guessing the
  national average by roughly 0.6 to 0.7 pp, and a spatial smoother is slightly
  better than a covariate model.
- The mean overstates it. The median gain is about a quarter of the mean, and
  only about 56 percent of replicates favour the model at all.
- The gain grows as more of the survey is retained, and is negative for the
  lowest-prevalence outcome.

No p-value is reported for these contrasts: replicates within a cell share one
underlying survey and are not independent.
(source: `results/tables/corrected/subsample_contrasts.csv`)

### Sentinel calibration is insurance, not improvement

Calibrating a transported map on k of the target country's own districts lowers
the mean held-out error, barely moves the median, and sharply cuts the upper
tail. It removes the catastrophic fold, where the map lands at the wrong level,
and does close to nothing for the typical one, which is what a location-only
correction should do. No k in the tested grid reaches within-country model error.
(source: `results/tables/corrected/calibration_learning_curve.csv`)

## Conclusions

1. **Predicting subnational deficiency in a country with no survey is not
   supported.** Transport discrimination is near chance once selection is
   honest, and a new country's biomarker level is unpredictable for iron.
2. **The obstacle is level, not pattern.** The transported map has roughly the
   right shape and the wrong height, and the height is a country intercept that
   no covariate in this pool supplies.
3. **Anchoring is the one intervention that addresses it**, worth about a third
   of the transport error for a single parameter, with ranking untouched.
4. **The anchor does not exist.** No current national iron estimate is available
   for these countries in any source checked. This is a data-acquisition
   constraint, not a modelling one.
5. **More or finer covariates are not the bottleneck.** Wider pools degrade
   transport, cluster-level covariates lose to doing nothing, and a
   covariate-free spatial smoother beats honestly-selected soil covariates.
6. **Within a country, the survey beats the model wherever the survey exists.**
   The model's defensible value is confined to districts deliberately not
   sampled, where it is worth a few tenths of a percentage point on average and
   nothing at all in nearly half of replicates.
7. **Admin-2 is below the design resolution of these surveys.** Much of the
   between-district spread being modelled is sampling noise, which caps what any
   method can achieve.

## Next steps and extensions

**Highest value.**

- **Acquire a usable national anchor.** This is the single highest-value action
  available: it is what converts a correctly-shaped map into a usable one, it
  costs a data request rather than a modelling programme, and it is currently
  missing for the outcome that needs it most. Candidate routes are a GBD
  extract, a current national survey figure, or a commissioned national estimate.
- **Acquire additional biomarker surveys.** More countries is the only thing
  that widens the four-country base the transport estimates rest on, and it is
  also what would narrow the new-country prediction intervals.

**Worth doing.**

- **Publish at Admin-1 alongside Admin-2.** The reliability ceiling nearly
  doubles at Admin-1. A more reliable estimate of a coarser unit is a real
  product, and the Admin-2 map should carry its reliability explicitly.
- **Run one or two national case studies** on the survey-design question, with
  the narrow claim: trust the direct estimate wherever you sample, and use the
  model only for districts you deliberately skip.
- **Bring in external costing expertise** to convert the zero-coverage edge into
  a dollar figure, which is a smaller and more honest ask than a general
  accuracy claim.
- **Adopt the uniform inflammation adjustment for iron** on pattern grounds, and
  the distributional estimator for tail-prevalence outcomes.
- **Verify the WHO and IZiNCG severity thresholds** in
  `metadata/who_severity_thresholds.csv` against their source documents before
  any exceedance map is presented. Every row is currently marked
  `verified_against_source = FALSE`, and two carry substantive caveats.

**Scoped narrowly if pursued.**

- **A GBD contribution** should be within-country-year disaggregation of an
  existing national estimate, constrained to reproduce the national total,
  rather than free-standing subnational prediction. Remote-sensing predictors
  cannot reliably extend back to 1990.

**Open technical items.**

- The 2 km GPS buffer for the cluster-level comparator undercorrects for survey
  confidentiality displacement of up to 5 to 10 km in rural clusters.
- The food-security domain accepts Admin-2 name matches at a 30 percent rate.
- Reconcile survey weighting with the pipeline's rare-outcome class weighting so
  the individual-level ensemble can be weighted without destabilising rare
  folds.

## Methods tested that do not improve estimation

- **XGBoost**, shallow and deep, and **random forest**: no transport advantage
  over a well-selected simple regression.
- **Full SuperLearner ensembles against a single linear model** on identical
  unscreened covariates: no advantage across countries. The ensemble's advantage
  is within country.
- **Larger candidate pools** of 100 to 150 variables: degrade transport relative
  to a parsimonious set.
- **Additional predictor domains** (DHS wealth and health, food prices,
  nighttime lights, market access, IHME): tested both forced into the search and
  left to compete freely, none displaced the small climate, soil and malaria
  set; the best DHS candidate ranked 51st of about 590.
- **Construct-based PCA reduction**: little value beyond the existing variables,
  though one vegetation component earns a place.
- **Rank-based objectives** in place of classification loss: no method
  significantly beat another, and the comparison is underpowered.
- **Domain generalization (IRM, Group DRO)** for few-environment transfer:
  neither improves on naive pooled regression with three training countries per
  fold, and Group DRO significantly underperforms. These methods need more
  environments than a four-country study supplies.
- **Cluster-level spatial resolution**: the reliability ceiling falls going
  finer than Admin-2, so a finer map is less reliable rather than more.
- **Reframing the evaluation** from classifying individuals to ranking
  districts: the same order of magnitude either way, because the constraint is
  the data rather than the metric.
- **Seasonality and interview timing** as a confounder: each survey window is
  short, interview month shows no association with outcomes, and timing is
  nearly collinear with district. Closed question.

## Sources

Figures above are read from `results/tables/corrected/`,
`results/sensitivity/` and `metadata/`, with the specific file cited under each
section. The cluster spatial model and subsample figures are being refreshed at
the time of writing; their conclusions are stable but the decimals may move
slightly.
