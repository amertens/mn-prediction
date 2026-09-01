# WS-E. Skill-curve symmetry and decomposition

## E1. The D8 exclusion was symmetric, so no rerun was needed

Stage 0 read `scripts/accuracy_impact/ws4b_skill_curve.R` and found `COVS`, the
reduced set with all 97 DHS-derived columns removed, used at **both** the
DHS-target block and the micronutrient block. The run confirms it at execution:
197 columns for both families.

**No rerun was performed and the figure is not republished.** The workstream's
conditional did not fire.

One consistency note that follows from the symmetry rather than contradicting
it: the skill curve's micronutrient correlations are computed on the reduced 197
set, while `anchor_controls.csv` uses the full 294. Two different unanchored
leave-one-region-out numbers for the same cells therefore exist in committed
tables, and they are not in conflict; they are different predictor sets.

## E2. The interpretation becomes a measured column

`results/tables/skill_curve_decomposition.csv`. Each target's out-of-fold
correlation is decomposed by residualising prediction and target on Admin-1.

| Family | Targets | Median r overall | Median r **between** region | Median r **within** region | Within positive | Variance share between |
|:---|---:|---:|---:|---:|---:|---:|
| DHS indicator | 81 | 0.076 | **0.118** | **0.073** | 64.2% | 0.401 |
| Micronutrient | 24 | 0.200 | **0.367** | **0.111** | 66.7% | 0.382 |

**At the median the claim holds, and for the micronutrient cells it holds
strongly.** The between-region correlation is **0.367** against a within-region
**0.111**, a factor of 3.3. About 40 percent of district-level variance is
between regions in both families, so the between-region component is doing more
with less.

**As a per-target regularity it does not hold.** `r_between` exceeds `r_within`
in **56 of 105** targets, which is 53 percent and not distinguishable from a
coin flip. The median tells one story and the per-target ordering does not
support it as a rule.

**Within-region signal is weak but not zero.** It is positive in 64 to 67
percent of targets, with a median of 0.073 to 0.111. The stated interpretation
that these covariates carry "almost nothing within region" is closer to right
than wrong for the typical target and overstates the case as an absolute.

## How this sits with WS-A and WS4a

The three lines of evidence agree in direction and differ in strength:

- **WS4a** found skill declining monotonically as resolution refines, with
  Admin-1 best in 9 of 14 cells. Consistent with a between-region signal.
- **WS-E2** measures that signal directly: 0.367 between against 0.111 within
  for micronutrients.
- **WS-A** found that the two predictors surviving a family-wise permutation
  correction do so in the **region-partialed** family, which is the within-region
  one. Soil pH and skilled birth attendance carry something after the
  between-region component is removed.

The reconciliation is that the within-region signal is weak on average and not
absent, and that pooling across cells is what makes it detectable. A single cell
cannot see it; sixteen cells across four countries can.

## Reviewer pass, statistical

**Shared-noise check.** Zero for every target. Predictions come from geospatial
and prior-round covariates; targets are the current survey's district estimates
or the DHS direct estimates. No respondent contributes to both.

**Identical treatment across families.** The same 197 covariates, the same ridge
with a top-20 in-fold prescreen, the same Admin-1 blocking, the same
decomposition. Only the target changes.

**The decomposition is descriptive, not inferential.** `r_between` is computed on
region means, so its sample size is the number of regions, which is 4 to 27. A
between-region correlation from four points is not a stable quantity and the
median across 105 targets is doing the work. `n_regions` is carried on every row.

**Degenerate-variance districts are excluded** from the DHS family on the same
rule the skill curve uses: a reported design variance below 1e-8 indicates a
single-PSU district and is dropped.

**Seeds.** 20260923; the ridge is deterministic.

## Reviewer pass, reproducibility

**Targets graph.** No new target.

**Additive.** `skill_curve_decomposition.csv` is new. `reliability_skill_curve.csv`
is not modified: the workstream proposed adding columns to it, and writing a
separate table instead preserves guardrail 4 and keeps the frozen snapshot
byte-comparable.

**Joins.** Through `admin2_join_by()`.

**Runtime.** About nine minutes for 105 targets.
