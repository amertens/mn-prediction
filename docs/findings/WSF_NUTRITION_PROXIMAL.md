# WS-F. The nutrition-proximal predictors

`results/tables/nutrition_proximal_confirmatory.csv`,
`results/figures/nutrition_proximal_confirmatory.png`.

## The hypothesis this was built to give its best chance

WS-A tested 294 harmonised predictors for consistent association across cells
and found **two** surviving family-wise correction. The obvious objection is
that the predictor set is the wrong one: climate, soil texture and vegetation
are distal to what a person eats, so a null result there says little about
whether *nutrition-proximal* covariates would work.

WS-F answers that objection on its own terms. Thirty-one predictors were added,
chosen for mechanistic closeness to intake rather than for availability:

| Block | Examples |
|:---|:---|
| Soil micronutrients | `soil_zinc_mean_0_20`, `soil_zinc_mean_20_50`, `soil_zinc_stdev_20_50` |
| Soil chemistry | `soilgrids_ph`, `soilgrids_organic_carbon` |
| Crop composition | `spam_share_cereals`, `spam_share_pulses`, `spam_share_roots`, `spam_share_oilcrops` |
| Programme coverage | `np_vita_capsule`, `np_deworming`, `np_iron_pregnancy` |
| Antenatal contact | `dhs_w_anc1` |

Nine came from the rdhs cache, which supplied them offline: `h33` for vitamin A
capsule receipt, `h43` for deworming, `m45_1` for iron in pregnancy, `hv234a`
for salt.

## The result

**Zero survivors, in every family and every cell set tested.**

| Family | Cell set | Predictors | Testable | q analytic < .05 | q permutation < .05 | FWER < .05 |
|:---|:---|---:|---:|---:|---:|---:|
| marginal | all | 31 | 29 | 3 | **0** | 0 |
| marginal | shared | 31 | 29 | 4 | **0** | 0 |
| region-partialed | all | 31 | 29 | 2 | **0** | 0 |
| region-partialed | shared | 31 | 29 | 6 | **0** | 0 |

The analytic false-discovery rate flags two to six predictors depending on the
family. **Every one of them is extinguished once multiplicity is calibrated
against the data's own permutation null.** The strongest analytic q in the whole
scan is 3.9e-07, for a predictor whose permutation q is 0.405.

The nearest miss is `spam_share_cereals` under region-partialing, at z 3.50 and
family-wise p 0.764.

## Why this matters more than the WS-A null

WS-A can be dismissed as a screen over predictors nobody expected to work. WS-F
cannot. These were nominated in advance on mechanism, they include the direct
programme-coverage measures a nutritionist would ask for first, and they were
tested on the same cells under the same protocol. The result is a cleaner null
than WS-A's, and it is the one that belongs in an abstract.

## The two survivors elsewhere are worth naming against this

WS-A's two family-wise survivors were `soilgrids_ph` (z 7.34, p 0.020) and
`dhs_w_skilled_attendant` (z 7.02, p 0.028). Soil pH governs micronutrient
bioavailability and skilled birth attendance indexes health-system access;
those are the soil-chemistry and access channels one would nominate in advance.
So the project's position is not that nothing is associated with deficiency. It
is that **two predictors out of 294 survive, none of the 31 chosen for
mechanistic proximity do, and no combination of them beats the national mean
under a deployable fold protocol.**

## Reviewer pass

**The completeness rule was fixed mid-analysis and the fix is recorded.** An
earlier version dropped Ghana to four usable predictors through a strict
completeness requirement combined with six unmatched districts. The scan now
uses an 80 percent completeness threshold with median imputation, and records
`imp_rate` per cell: 0.0 percent for Sierra Leone, 0.9 for Gambia, 1.9 for
Malawi, 6.4 for Ghana.

**Cell coverage.** 24 cells, 29 of 31 predictors testable in each; the two
dropped are constant within at least one country.

**Empirical ceilings are reported per cell** and vary from 0.000 to 0.848, so a
null in a cell with `r_max` near zero carries no information and is not counted
as evidence against the predictor.

**Additive.** `nutrition_proximal_confirmatory.csv` is new.
