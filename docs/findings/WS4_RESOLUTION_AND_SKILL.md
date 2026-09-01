# WS4. The reliability-skill curve

## Scope actually delivered

**WS4b is delivered in full.** WS4a, the aggregation-level sweep across Admin-1,
Admin-1 split into two and three contiguous parts, and Admin-2, is **not yet
computed**. The reason is compute: the machine was running the repository
owner's own jobs throughout this session, and WS3's individual-level fits
progressed at roughly six of 192 rows per half hour under that load. WS4a needs
an empirical reliability and a model fit at each of four aggregation levels for
24 cells. It is recorded here as not yet computed rather than approximated.

## The finding, and a defect removed before it

**Removing the leak first.** The first run of `ws4b_skill_curve.R` returned
r 0.946 for Gambia anaemia with `r_share` 1.05. Ninety-seven of the 294
harmonised Admin-2 covariates are DHS-derived aggregates, including
`dhs_AN_ANEM_W_ANY`, so predicting the DHS indicator `AN_ANEM_W_ANY` was partly
the outcome predicting itself. All 97 are excluded, from **both** target families,
because excluding them only where they leak would have left the two families
fitted on different predictor sets and the curve would compare models rather than
targets.

### The micronutrient cells are not the problem

| Family | Targets | Median `r_max` | Median achieved `r` | Median `r_share` |
|:---|---:|---:|---:|---:|
| DHS indicator | 81 | 0.784 | **0.071** | 0.109 |
| Micronutrient | 24 | 0.613 | **0.200** | 0.359 |

(source: `results/tables/reliability_skill_curve.csv`.)

Fitting a line through the DHS indicators and placing the micronutrient cells on
it, the micronutrient cells sit **+0.219 above** that line on average and below
it in only **5 of 24** cells.

**This reverses the direction of Section 11 claim 2.** The claim is that low
micronutrient skill is a property of the target rather than the predictors, and
that a well-measured district quantity would be predicted well. Measured on 81
DHS indicators with design-based reliabilities at a median of 0.784, the same
pipeline achieves a median out-of-fold correlation of **0.071**. Well-measured
district quantities are **not** predicted well by these covariates either. The
micronutrient cells are among the better-performing targets, not the worse.

### Reliability does not predict skill

The correlation between `r_max` and achieved `r` across all 105 targets is
**0.043**. Whatever limits achieved skill at Admin-2, it is not the reliability
of the target. That is the opposite of the framing in Sections 2.4 and 11.1,
which treat the ceiling as the binding constraint.

Taken with WS1, the picture is coherent and different from the published one:
the ceiling is much higher than reported (median 0.613, not 0.098), models reach
far less of it than reported, and the shortfall is not explained by how noisily
the target is measured. The constraint sits in the predictors, or in the
Admin-2 aggregation, and not in the target's reliability.

### What this does not establish

The DHS indicators come from a different survey round than the micronutrient
outcomes (Ghana 2014 against 2017, and so on), so the covariates are contemporary
with neither in the same way. The comparison is between targets measured on
different samples; it is not a within-sample contrast. And the reliability used
for the DHS family is design-based from `direct.var`, while the micronutrient
family uses the WS1a split-half estimate. The two are different estimators of the
same quantity, and WS1 showed estimator choice can move `r_max` by 0.4. The
median `r_max` comparison across families should therefore be read loosely; the
achieved-`r` comparison, which needs no reliability at all, is the solid part.

The design-based reliability also inherits the single-PSU problem: a district
with one cluster returns a variance near zero, which would overstate reliability.
Districts whose reported variance is below 1e-8 are excluded and the count is
carried in `n_degenerate_var`.

## Reviewer pass, statistical

**Identical models across families.** Both families use the same 197 non-DHS
covariates, the same `.ds_fit` ridge with a top-20 in-fold prescreen, and the
same leave-one-region-out blocking. The only difference between a DHS row and a
micronutrient row is the target.

**Prescreen placement.** Inside the fold, via `.ds_fit(k_screen = 20)`. Section
2.3 measures +0.182 optimism from screening on all data; that path is not taken.

**Leakage.** Addressed above and it was the first thing found. The residual risk
is a non-DHS covariate that encodes a DHS indicator, which the WS7a leakage
report would surface for the micronutrient targets but does not cover for the
DHS pseudo-targets.

**Seeds.** 20260907. The ridge is deterministic; no stochastic call affects the
result.

**Joins.** Through `admin2_join_by()`, so Malawi is keyed on the pair.

## Reviewer pass, reproducibility

**Targets graph.** No new target; the script writes to `results/tables/` outside
the graph, as `scripts/covariates/` scripts do.

**Stamp targets.** Reads `data/DHS/clean/*_dhs_admin2_direct.rds` and
`data/covariates/harmonized/predictors_admin2_harmonized.csv`, both tracked.

**Paths.** `here()` throughout, no absolute path, no `setwd()`.

**Figure.** `results/figures/reliability_skill_curve.png`, base graphics, no
ggplot dependency.

**Runtime.** About six minutes for 105 targets.
