# Multilevel regression and poststratification (MRP)

Added as a comparator in the admin-2 benchmark suite. Implementation:
[`R/mrp.R`](../R/mrp.R). Standalone comparison run:
`Rscript scripts/covariates/06_mrp_comparison.R` →
[`results/tables/mrp_comparison.csv`](../results/tables/mrp_comparison.csv).

## Why it is here

MRP is the reference small-area method a survey-statistics reviewer expects, and
its absence from a 21-method benchmark table was conspicuous. It is also the
individual-level analogue of the area-level models already compared: instead of
modelling a district's direct estimate and its sampling variance (Fay–Herriot,
BYM2), it models **individuals** and rebuilds the district from a
poststratification table.

It addresses a bias the manuscript already concedes — aggregating shrunken
person-level predicted probabilities up to a district mean is not an unbiased
district estimator. Poststratifying to a frame is the principled fix.

## Specification

1. **Model.** `glmer(deficient ~ poststratification variables + k screened area
   covariates + (1 | Admin1) + (1 | Admin1:Admin2), family = binomial)`, fitted
   **unweighted**. That is standard MRP practice: the design information enters
   through the poststratification table, not the likelihood.
2. **Area covariates.** At most three, screened by absolute correlation with the
   district survey prevalence. With 14–87 areas, more than a few area terms in
   the fixed part overfit immediately.
3. **Poststratification.** Survey-weighted cell counts per district over the
   poststratification variables; predict every cell, take the cell-count-weighted
   mean.

## The frame, and its limitation

Textbook MRP poststratifies to a **census** frame. No census microdata is
available for these countries in this project, so the frame is built from the
survey's own design weights: within each district, the weighted cell
distribution estimates the population distribution.

This is the standard fallback when no census frame exists, and it is not
vacuous — the biomarker subsample is not compositionally representative of the
district the full survey describes, and this corrects for that. But it cannot
correct for coverage error in the survey frame itself, which makes MRP here a
**model-assisted smoother of the design-weighted estimate** rather than a fully
independent estimator. Quote the numbers with that attached.

Obtaining a real frame (IPUMS International census microdata: Ghana 2010, Malawi
2008, Sierra Leone 2015) would upgrade this from a smoother to a genuine MRP and
is the obvious next step if the comparator is to carry weight in print.

**Out of sample**, a district with no survey data takes its cell distribution
from its Admin-1 parent (then the national distribution), and its unobserved
district random effect is set to the prior mean of zero — so unsurveyed
districts are separated only by their area covariates. That is the honest
behaviour of MRP out of sample, not a limitation of this implementation.

## Poststratification variables are discovered by pattern

The individual-level survey vocabulary drifts between countries exactly as the
covariate vocabulary does:

| | residence | wealth | sex |
|---|---|---|---|
| Gambia | `gw_urban` | `gw_hWealthquintile` | `gw_cSex` |
| Ghana | `gw_urban/ rural` | `gw_hWealthquintile` | `gw_cSex` |
| Sierra Leone | — | `gw_hWealth**Q**uintile` | `gw_cSex` |
| Tanzania | — | — | `gw_sex` |

A hardcoded name list reduced Ghana and Sierra Leone to *no* usable variables.
`MRP_PS_PATTERNS` resolves one variable per conceptual dimension by regex,
rejects anything with more than 8 levels or more than 40% missing, and logs the
per-country resolution so it is never silent.

## Coverage

MRP runs in **16 of 24** country × outcome cells.

**Malawi (8 cells) cannot be fitted.** Its analytic table carries 1,403 columns
of which exactly one is individual-level (`gw_cnum`); the individual
demographics were not carried through the Malawi merge. This is a data
limitation, not a model failure — recovering it means going back to the raw
Malawi survey. Reported as skipped with the reason logged.

**Tanzania** is not yet in `get_country_configs()`; harmonised covariates and
AlphaEarth already cover it, so MRP will run there as soon as it is activated.

## Results (in-sample, against the district survey estimate)

Metrics match how the other rows of `area_comparison_all.csv` are computed —
they compare methods like-for-like and are **not** out-of-sample accuracy. The
honest out-of-sample numbers come from the block-CV corrected-methods analysis.

| | median Pearson r | median MAE (pp) |
|---|---|---|
| MRP (legacy covariates) | 0.752 | 5.84 |
| MRP (harmonised covariates) | 0.738 | 5.01 |
| National mean (null) | — | 9.55 |

Re-run against the verified 208-predictor build (removing 11 duplicate
population-density columns did not move any figure, because the duplicates were
never among the screened area covariates).

MRP beats the national-mean null on MAE in **16 of 16** fitted cells with
harmonised covariates (15 of 16 with legacy). Between the two covariate sets it
is close to a wash: harmonised has the lower MAE in 10 of 16 cells and the lower
median MAE, legacy the marginally higher median correlation. That is the
expected result given the noise ceiling — a better covariate set cannot buy much
when the target's own sampling error dominates.

The comparator is wired into `run_area_comparison()` and appears as the `MRP`
row in `area_comparison_all.csv` once the pipeline is rebuilt.
