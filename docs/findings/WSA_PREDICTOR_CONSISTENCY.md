# WS-A. Does any predictor carry a consistent association?

The central question of the session, asked with pooling rather than cell by cell.
`results/tables/predictor_consistency_meta.csv`, 294 predictors, two families,
two cell sets, 500 permutations, seed 20260920.

## Why this is not the FDR screen already in the project

Section 3 reports that no predictor survives FDR control in any of the 24 cells.
That screen tests each predictor separately **inside** a single country-outcome
with 14 to 87 districts, so a predictor with a real but modest association
everywhere fails it in every cell while being exactly what the project wants.
This pools across cells and countries and tests the pooled quantity.

## The answer in one paragraph

**Two predictors survive a family-wise permutation calibration, both in the
region-partialed family, and neither is what the mechanism hypothesis
predicted.** They are `soilgrids_ph` (meta z 7.34, family-wise permutation
p 0.020) and `dhs_w_skilled_attendant` (z 7.02, p 0.028), each estimated from
three countries and twelve cells with nine of twelve cells agreeing in sign.
**Nothing survives in the marginal family.** The analytic normal-theory
q-values are badly anti-conservative at four clusters, claiming up to 22
survivors where the permutation supports at most two, which is why the
calibration was specified.

## The two calibrations, and why they disagree

| Family | Cell set | Testable | q analytic < 0.05 | q permutation < 0.05 | Family-wise p < 0.05 |
|:---|:---|---:|---:|---:|---:|
| marginal | shared | 281 | **22** | 0 | 0 |
| marginal | all | 281 | 4 | 0 | 0 |
| region-partialed | shared | 281 | **16** | 0 | **2** |
| region-partialed | all | 281 | 10 | 0 | **1** |

**The analytic column is not trustworthy.** It treats the random-effects z as
standard normal with four clusters, and the permutation null shows that
assumption fails badly: the smallest analytic q is 5.9e-11 for a predictor whose
permutation family-wise p is 0.020.

**The permutation BH column is resolution-limited and should not be read as
evidence of absence.** With 500 permutations the smallest attainable
per-predictor p is 1/501 = 0.002, and only one to three predictors per family
reach that floor, so the smallest attainable BH q is 0.187 to 0.421. **BH on
permutation p-values cannot reach 0.05 in this design regardless of the truth.**
Reporting "zero survive under permutation" without that caveat would overstate
the negative.

**The family-wise max-|z| calibration is the sound statement.** It compares the
observed |z| against the distribution of the maximum |z| over all 294 predictors
under the null, needs no independence assumption, and is adequately estimated at
the 0.05 level by 500 draws.

## The survivors

| Predictor | Family | Cell set | z | Countries | Cells | Sign agreement | Family-wise p |
|:---|:---|:---|---:|---:|---:|---:|---:|
| `soilgrids_ph` | region-partialed | shared | 7.34 | 3 | 12 | 9 of 12 | **0.020** |
| `dhs_w_skilled_attendant` | region-partialed | shared | 7.02 | 3 | 12 | 9 of 12 | **0.028** |
| `soilgrids_ph` | region-partialed | all | 6.62 | 3 | 18 | 12 of 18 | **0.044** |

**Soil pH is mechanistically interpretable and the sign is the expected one.**
Zinc and iron bioavailability to crops falls as soil pH rises, so alkaline soils
should be associated with more deficiency. The positive z is that direction. It
is the only predictor in the whole scan for which a mechanism and a surviving
association point the same way.

**Skilled birth attendance is not a nutrition variable.** It is a health-service
access measure and sits in the same family as the antenatal-care, improved-water
and improved-sanitation predictors immediately below it in the ranking
(z 5.76, 4.19, 4.09). That is the socioeconomic axis WS4b identified, and it
survives even after region-partialing.

**Soil zinc has the right sign and does not survive.** `soil_zinc_mean_0_20` and
`soil_zinc_mean_20_50` reach z of -3.37 in the region-partialed shared family,
meaning more soil zinc is associated with less deficiency, which is the
mechanistically correct direction. Their family-wise p-values are 0.948 and
0.946. The direction is encouraging and the evidence is not there.

## Why the region-partialed family has survivors and the marginal one does not

This looks backwards and is not. Partialing on Admin-1 removes the dominant
between-region variance from both sides. Where a predictor's association lives
mostly within region, removing the between-region component raises the
signal-to-noise of what remains rather than lowering it. The marginal family
carries every predictor that tracks the between-region socioeconomic gradient,
which makes the marginal null harder to beat because the maximum |z| under
permutation is itself larger there.

## Heterogeneity

Median tau across predictors is 0.128 to 0.399 depending on family and cell set.
Between-country heterogeneity in the association is therefore of the same order
as the associations themselves, which is the quantitative form of the transport
problem this project keeps meeting.

## Reviewer pass, statistical

**Shared-noise check.** Zero for every arm. The predictors are geospatial
rasters, agro-ecological zones, soil chemistry, and DHS aggregates from a
**prior survey round**. No respondent contributes to both a predictor and the
outcome. The one channel that is not zero is shared district-level confounding
between prior-round DHS aggregates and the current survey, which is a different
concern from shared sampling noise and is exactly what the region-partialed
family is designed to reduce.

**Clustering.** Outcomes within a survey share respondents, fieldwork and assay
batch, so the 16 cells are not 16 independent tests. Country is the clustering
unit: z values are averaged within country first, then combined across the four
countries with DerSimonian-Laird. Four clusters is few, which is the entire
reason for the permutation calibration.

**The permutation preserves what matters.** Outcomes are shuffled **within
region within country**, so the between-region structure that WS4a showed
carries most of the signal is preserved and only the within-region link is
destroyed. A null permuting across regions would be far easier to beat and would
manufacture significance out of the socioeconomic gradient.

**Disattenuation.** Correlations are divided by the cell's empirical `r_max`
before the Fisher transform, capped at 0.99, with an undisattenuated column
carried alongside. Disattenuation inflates variance and a reader should be able
to see both.

**Predictor availability differs by cell.** 281 of 294 predictors are testable,
and a predictor complete in Sierra Leone (293 usable) may be incomplete in Ghana
(176 usable). `k_countries` and `n_cells` are carried on every row so a z from
three countries is never read as a z from four.

**Multiplicity is handled three ways and all three are reported.** The
family-wise statement is the one to quote.

## Reviewer pass, reproducibility

**Targets graph.** No new target.

**Seeds.** 20260920, covering the 500 permutations in all four family-by-cellset
combinations.

**Joins.** Survey to covariates through `admin2_join_by()`.

**A bug found and fixed before any result was read.** The first run failed
because `sd()` returns NA on a column with missing values, so the
constant-column mask contained NAs and the subscripted assignment errored. The
mask now applies the per-cell completeness filter first and coerces NA to FALSE.

**Additive.** `predictor_consistency_meta.csv`,
`predictor_consistency_top30.csv` and `results/figures/predictor_consistency.png`
are new files.

**Runtime.** About 55 minutes for four family-by-cellset combinations at 500
permutations. `PROFILE=smoke` runs one cell at 50 permutations in about three
minutes.
