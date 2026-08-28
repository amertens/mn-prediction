# Micronutrient prediction project: status summary

What the analysis does, what it found, and where it can go next. Workstream
detail is in `docs/findings/`.

## What the analysis does

- **Scope.** Four countries with biomarker microdata: Gambia, Ghana, Sierra
  Leone, Malawi. Outcomes are biomarker-defined deficiency for vitamin A, iron,
  folate, B12 and zinc, in preschool children and women of reproductive age.
- **Predictors.** Roughly 600 pooled proxy variables that need no survey:
  remote sensing (SoilGrids, iSDAsoil, climate, vegetation, nighttime lights),
  Malaria Atlas, IHME modelled disease burden, harmonized DHS Admin-2
  indicators, food prices, travel time and market access, helminths, crop mix.
- **Two modelling tracks.** An individual-level SuperLearner ensemble aggregated
  to Admin-2, and an Admin-2 area-level model. Both are scored within country
  and under leave-one-country-out (LOCO) transport.
- **A leakage-free evaluation layer.** Spatially blocked cross-validation,
  out-of-fold calibration, design-aware Admin-2 error metrics, split-conformal
  and design-based intervals, out-of-support trust flags, and a reliability
  ceiling that says how much correlation sampling noise permits.
- **Selection nested inside the folds.** Predictor selection sees only the
  training countries, so a reported transport number is out of sample with
  respect to the selection as well as the data.
- **Harmonized outcomes.** One WHO cut-point per outcome and one
  inflammation-adjustment protocol per biomarker, so a cross-country comparison
  compares prevalences and not definitions.

## Findings

### Cross-country transport is close to chance

**Individual level.** The production predictor set comes from a three-step search
over roughly 590 pooled candidates: rank each candidate by its own cross-country
transport strength, collapse near-duplicate variants of the same construct, then
add predictors greedily for as long as leave-one-country-out AUC improves. It
settles on five variables: soil calcium variability, November rainfall, malaria
incidence, a vegetation principal component, and MAP parasite rate. A small,
interpretable set of that kind is the natural thing to report.

Scoring it on the same folds the search consulted gives the first column below.
Rerunning the entire search inside each fold, so it sees only the three training
countries, gives the second. Both columns use the same estimator on the same 16
folds, so the gap between them is the search, not the model.

| Scorer | Search sees every fold | Search nested in the fold |
|---|---|---|
| SuperLearner | 0.5977 | 0.5214 |
| glm | 0.6199 | 0.5380 |

Nesting the search costs about 0.08 AUC and it costs it almost everywhere:
the nested arm is worse in 15 of the 16 folds under either scorer, and four of
its SuperLearner folds fall below 0.5.

The five variables are not the problem. The search that found them ranked
candidates on the same folds it was later scored on, so most of the apparent
transport signal is a property of the search rather than of the predictors. Any
selection procedure run this way will produce a set that looks transportable.

**Area level.** Eleven model families predict district prevalence from the proxy
covariates, trained on three countries and scored on the fourth, over the same
16 outcome-by-country cells. Each map is placed at the training countries' mean
level, which is what a deployment without local data would have to do.

Three metrics, each answering a different question:

- **Pearson r** measures the *relative ranking*: whether the map puts districts
  in the right order. Correlation is unchanged by a constant shift, so it says
  nothing about the absolute level.
- **r_share** rescales that correlation by the reliability ceiling, the
  correlation attainable if the model were perfect and only sampling noise in
  the survey estimates stood in the way. The ceiling for these cells is 0.305,
  so a model at r = 0.15 has captured about half of what the data permits.
  Without this, a low correlation looks like a modelling failure when it is
  often a data limit.
- **RMSE** and **level bias** together separate the two ways a map goes wrong.
  RMSE is total error and absorbs both a wrong ranking and a wrong absolute
  level. Level bias isolates the absolute level. A model with modest RMSE and
  large level bias ranks districts about right and places all of them too high or
  too low.

| Model | What it is | Pearson r | r_share | RMSE (pp) | Level bias (pp) |
|---|---|---|---|---|---|
| ridge, all predictors, z-scored | every predictor, ridge penalty; outcome z-scored within country | 0.2242 | 0.735 | 17.81 | 8.73 |
| decorrelated 20, country-centred | 20 predictors, one per correlated cluster; covariates centred per country | 0.1715 | 0.562 | 15.48 | 9.53 |
| ridge, all predictors, country-centred | every predictor, ridge; covariates centred per country | 0.1694 | 0.555 | 19.70 | 13.13 |
| elastic net, top 30, country-centred | top 30 by correlation screen; covariates centred per country | 0.1439 | 0.472 | 21.93 | 15.43 |
| curated 16-variable set | 16 variables chosen in advance, one per subject-matter construct | 0.1349 | 0.442 | 18.17 | 9.08 |
| correlation-screened top 30 | top 30 predictors by correlation with the outcome | 0.1345 | 0.441 | 18.57 | 8.67 |
| **production parsimonious recipe** | the pipeline default: elastic net on the top 30 screened predictors | **0.1228** | **0.403** | **14.51** | **8.90** |
| **covariate-free spatial spline** | thin-plate spline on district coordinates; no covariates at all | **0.1151** | **0.377** | **16.53** | **11.52** |
| decorrelated 20, z-scored | the 20 de-duplicated predictors, outcome z-scored within country | 0.1031 | 0.338 | 18.56 | 8.98 |
| production ridge | the production covariate set with a ridge penalty | 0.0964 | 0.316 | 16.13 | 9.20 |
| spatial spline plus decorrelated 20 | the spline with the 20 de-duplicated predictors added | 0.0936 | 0.307 | 20.05 | 12.12 |

Every family reaches between a third and three quarters of the ceiling, so none
of them is close to exhausting what the data holds, and none is far off the
others either.

The ordering rewards a design choice as much as a model. The z-scored variants
standardise the outcome inside each country, which removes between-country level
and spread from the target and leaves the model fitting relative rankings alone.
That is why they lead on Pearson r while their level bias stays ordinary. They
bound how much of the relative ranking is learnable, and they are not a
deployable recipe.

Among the variants that keep the full problem, the parsimonious production recipe
and the covariate-free spline land within 0.008 of each other, and the production
recipe posts the lowest RMSE in the table at 14.51 pp. A spline that knows only
where a district is does as well as the covariates do. Adding those covariates to
the spline pushes it to the bottom.

**Selection honesty accounts for much of the apparent signal.** The
eight-SoilGrids spatial model, re-scored with its features chosen inside each
fold:

| Variant | Pearson r | Spearman | RMSE (pp) |
|---|---|---|---|
| features locked from a full-data ranking | 0.3324 | 0.3164 | 18.26 |
| covariate-free spatial spline, same cells | 0.2861 | 0.2739 | 16.50 |
| features selected inside each fold | 0.1426 | 0.1170 | 18.24 |

Choosing the soil features honestly more than halves the correlation. It also
falls below using no soil features at all.

The two tables come from separate experiments whose ceilings differ, 0.205 here
against 0.305 above. Read each within itself. No model family in the pool
changes the conclusion: proxies alone cannot place an unsurveyed country's map
at the right absolute level.

### The residual error is a country intercept, not a covariate failure

On the log biomarker scale, each fold's national bias splits exactly into two
terms, a country location offset and whatever level the covariates assert. The
intercept is the larger term in three of the four outcomes. For child iron the
training pool sits between 0.523 and 2.669 times the held-out country's level.

Raw child ferritin medians span a factor of 6.30 across the four countries. A
single uniform BRINDA CRP+AGP protocol leaves a factor of 4.85 and moves the
mean absolute national bias by about one percentage point on a bias of fifteen.
It does sharpen the transported relative rankings: child iron Pearson rises from
0.1808 to 0.3097, women iron from -0.1865 to 0.0207. The pipeline offers it as a
scheme on those grounds.

Harmonizing measurement improves the relative rankings. It does not correct the
absolute level.

### Anchoring to a national estimate is the effective remedy

A national anchor is a known country-level prevalence figure that fixes the
absolute level of a transported map. One constant shift on the logit scale makes
the map's weighted aggregate equal the anchor. The shift is monotone, so the
relative ranking of districts survives untouched: it corrects the absolute level
without inventing spatial signal. It attacks precisely the error term above.

| Anchor | Absolute national bias | MAE | RMSE |
|---|---|---|---|
| none | 9.135 pp | 13.671 pp | 16.476 pp |
| the country's own survey | 0.000 pp | 8.934 pp | 11.499 pp |
| a published external estimate | 44.141 pp | 44.690 pp | 46.386 pp |

A perfect anchor removes 35 percent of the mean absolute error and improves 9 of
16 folds. A country's own national value is of course unavailable when deploying
to an unsurveyed country; that row measures what a perfect anchor would be
worth.

The external row is the binding constraint. Three of 16 primary
country-outcome cells have an external anchor at all. None of the four countries
has an iron entry. The three vitamin A anchors predate their surveys by 15 to 19
years. A stale anchor is worse than no anchor, because the shift trusts it
completely.

### A new country's level is predictable for vitamin A and not for iron

A country random-intercept model on the Admin-2 mean log biomarker gives the
range a genuinely new country's national level would fall in:

| Outcome | Interval, marker units | Span |
|---|---|---|
| child vitamin A (RBP, umol/L) | 0.662 to 1.536 | 2.3x |
| women vitamin A (RBP, umol/L) | 0.961 to 2.174 | 2.3x |
| women iron (ferritin, ug/L) | 6.12 to 105.02 | 17x |
| child iron (ferritin, ug/L) | 2.23 to 176.40 | 79x |

Four countries make these intervals wide, and the width is itself the finding. A
factor of 79 for child ferritin will not support any claim about an unsurveyed
country's iron level.

### Admin-2 is finer than these surveys can resolve

The reliability ceiling falls as the unit gets finer: 0.59 at Admin-1, 0.31 at
Admin-2, 0.22 at cluster level. The median district contributes few biomarker
measurements, and in three of the four countries it holds a single survey
cluster.

This constraint sits behind every weak result here. What limits cross-country
learning is the number of areas, 14 to 89 per country, not the number of
individuals.

### Covariates do not help at cluster level

A point-referenced spatial model on the survey GPS coordinates, held out one
district at a time so the smoother must extrapolate to an unsurveyed location,
ranks both covariate arms below simply giving a district the national average.
Adding covariates to the spatial smoother makes it worse than the smoother
alone. Spatial borrowing beats the national-mean null in about half the cells,
by a fraction of a percentage point of MAE.

### The distributional estimator helps only in the tail

Modelling the continuous biomarker and integrating past the cut-point, instead
of dichotomising first, lifts area-level correlation from -0.0592 to -0.0326 and
costs discrimination: AUC drops from 0.5944 to 0.5887, and Brier skill, MAE and
the ceiling-adjusted correlation all fall. That is across the 23 cells whose
ceiling exceeds zero.

The gains concentrate where the efficiency argument predicts, in outcomes whose
cut-point sits in a tail. Every one of the largest gains is a cell under 10
percent prevalence. The pipeline offers it as a declared alternative, not the
default.

### Modelling helps only for districts a survey skips entirely

Reduced survey budgets, 50 to 90 percent of clusters retained, 25 replicates per
fraction, scored against each country's full survey under district-level
cross-validation:

- Where a district keeps even a few clusters, the direct survey estimate beats
  the model by about 5 percentage points of RMSE. The model wins a tenth of
  replicates.
- Where a district keeps none, the model beats the national average by 0.6 to
  0.7 pp, and a spatial smoother edges out a covariate model.
- The mean overstates that gain. The median is a quarter of it, and 56 percent
  of replicates favour the model at all.
- The gain grows as more of the survey is kept, and turns negative for the
  lowest-prevalence outcome.

These contrasts carry no p-value. Replicates within a cell share one survey and
are not independent.

### Sentinel calibration is insurance, not improvement

Calibrating a transported map on k of the target country's own districts lowers
the mean held-out error, barely moves the median, and cuts the upper tail
sharply. It rescues the catastrophic fold, where the map lands at the wrong
absolute level, and leaves the typical fold alone. A location-only correction should do
exactly that. No k in the tested grid reaches within-country model error.

## Conclusions

1. **Predicting subnational deficiency in a country with no survey is not
   supported.** Transport discrimination sits near chance once selection is
   honest, and a new country's biomarker level is unpredictable for iron.
2. **The obstacle is absolute level, not relative ranking.** The transported map
   ranks districts roughly correctly and places all of them at the wrong
   absolute level, which is a country intercept no covariate in this pool
   supplies.
3. **Anchoring is the one intervention that addresses it.** One parameter buys
   about a third of the transport error and leaves ranking untouched.
4. **The anchor does not exist.** No current national iron estimate covers these
   countries in any source checked. That is a data problem, not a modelling one.
5. **More or finer covariates are not the bottleneck.** Wider pools degrade
   transport, cluster-level covariates lose to doing nothing, and a
   covariate-free spatial smoother beats honestly selected soil covariates.
6. **Within a country, the survey beats the model wherever the survey reaches.**
   The model earns its place only in districts nobody sampled, where it is worth
   a few tenths of a percentage point on average and nothing at all in half the
   replicates.
7. **Admin-2 sits below the design resolution of these surveys.** Much of the
   between-district spread being modelled is sampling noise, which caps what any
   method can achieve.

## Next steps and extensions

**Highest value.**

- **Acquire a usable national anchor.** This is the single highest-value action
  available. It converts a correctly ranked map into a usable one, it costs a
  data request instead of a modelling programme, and it is missing for the
  outcome that needs it most. Candidates are a GBD extract, a current national
  survey figure, or a commissioned national estimate.
- **Acquire additional biomarker surveys.** More countries is the only thing
  that widens the four-country base the transport estimates rest on, and the
  only thing that narrows the new-country intervals.

**Worth doing.**

- **Publish at Admin-1 alongside Admin-2.** The reliability ceiling nearly
  doubles at Admin-1. A more reliable estimate of a coarser unit is a real
  product, and the Admin-2 map should carry its reliability on its face.
- **Run one or two national case studies** on the survey-design question, with
  the narrow claim: trust the direct estimate wherever you sample, and use the
  model only for districts you deliberately skip.
- **Bring in external costing expertise** to price the zero-coverage edge, a
  smaller and more honest ask than a general accuracy claim.
- **Adopt the uniform inflammation adjustment for iron** on pattern grounds, and
  the distributional estimator for tail-prevalence outcomes.
- **Verify the WHO and IZiNCG severity thresholds** before any exceedance map is
  presented. Every threshold is currently unverified against its source
  document, and two carry substantive caveats.

**Scoped narrowly if pursued.**

- **A GBD contribution** should disaggregate an existing national estimate
  within a country-year, constrained to reproduce the national total, and not
  attempt free-standing subnational prediction. Remote-sensing predictors cannot
  reliably reach back to 1990.

**Open technical items.**

- The 2 km GPS buffer for the cluster-level comparator undercorrects for survey
  confidentiality displacement, which reaches 5 to 10 km in rural clusters.
- The food-security domain accepts Admin-2 name matches at a 30 percent rate.
- Survey weighting and the pipeline's rare-outcome class weighting need
  reconciling before the individual-level ensemble can be weighted without
  destabilising rare folds.

## Methods tested that do not improve estimation

- **XGBoost**, shallow and deep, and **random forest**: no transport advantage
  over a well-selected simple regression.
- **Full SuperLearner ensembles against a single linear model** on identical
  unscreened covariates: no advantage across countries. The ensemble earns its
  keep within a country.
- **Larger candidate pools** of 100 to 150 variables: worse transport than a
  parsimonious set.
- **Additional predictor domains** (DHS wealth and health, food prices,
  nighttime lights, market access, IHME): forced into the search and also left
  to compete freely, none displaced the small climate, soil and malaria set. The
  best DHS candidate ranked 51st of about 590.
- **Construct-based PCA reduction**: little value beyond the existing variables,
  though one vegetation component earns a place.
- **Rank-based objectives** in place of classification loss: no method beat
  another, and the comparison was underpowered.
- **Domain generalization (IRM, Group DRO)** for few-environment transfer:
  neither improves on naive pooled regression with three training countries per
  fold, and Group DRO does markedly worse. These methods need more environments
  than a four-country study supplies.
- **Cluster-level spatial resolution**: the ceiling falls below Admin-2, so a
  finer map is less reliable, not more.
- **Reframing the evaluation** from classifying individuals to ranking
  districts: the same order of magnitude either way, because the data constrains
  the answer and the metric does not.
- **Seasonality and interview timing** as a confounder: each survey window is
  short, interview month shows no association with outcomes, and timing is
  nearly collinear with district. Closed question.

---

The cluster spatial model and subsample figures are being refreshed. Their
conclusions are stable; the decimals may move slightly.
