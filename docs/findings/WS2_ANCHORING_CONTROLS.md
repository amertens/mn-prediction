# WS2. Anchoring: controls and circularity tests

Section 4 is the project's only strongly positive finding. Section 4.4 states the
circularity concern and argues against it. This workstream tests the argument
rather than restating it.

All arms are scored on the **same 24 cells** under the **same leave-one-region-out
folds**, from `results/tables/anchor_controls.csv`. The published arms reproduce
to within 0.004 in mean `r` (hard anchor 0.409 here against 0.413 published,
MAE 8.91 against 8.85), the difference being rounding at four decimal places
rather than three before averaging.

## The two-paragraph answer

**The anchoring gain does not survive a jackknife.** When each district's
regional anchor is computed from the region's *other* districts only, so that no
district contributes to its own correction, mean `r` falls from **0.409** to
**0.147**, against **0.156** for no anchor at all. The gain over the unanchored
model goes from **+0.253** (better in 22 of 24 cells) to **-0.001** (better in 8
of 24), and MAE moves from 1.86 pp better to 1.56 pp **worse**. The published
gain is an artifact of each district contributing roughly a quarter to a third
of the number used to correct it.

**A predictor with no covariates at all beats every covariate model.** Assigning
each district its region's design-based survey estimate, and nothing else,
reaches mean `r` **0.516**, MAE **7.38 pp** and mean absolute bias **0.77 pp**,
against **0.409 / 8.91 / 1.61** for the anchored covariate model and
**0.156 / 10.77 / 3.21** for the unanchored one. It is better than the anchored
covariate model on every metric reported. The information in Section 4 is the
survey's regional estimates; the covariate pattern subtracts from it.

## The arms

| Arm | Cells | Mean r | Median r | MAE pp | Mean abs bias pp | Median r_share (empirical ceiling) |
|:---|---:|---:|---:|---:|---:|---:|
| No anchor (LORO) | 24 | 0.156 | 0.179 | 10.77 | 3.21 | 0.33 |
| **Flat REGIONAL mean, no covariates** | 24 | **0.516** | **0.547** | **7.38** | **0.77** | 0.86 |
| Flat NATIONAL mean, no covariates | 24 | undefined | undefined | 11.12 | 4.25 | not applicable |
| Admin-1 anchor (hard) | 24 | 0.409 | 0.378 | 8.91 | 1.61 | 0.67 |
| Admin-1 anchor (shrunk) | 24 | 0.310 | 0.303 | 9.42 | 2.70 | 0.53 |
| National anchor | 24 | 0.161 | 0.173 | 12.07 | 5.85 | 0.32 |
| **Admin-1 anchor (hard, JACKKNIFE)** | 24 | **0.147** | 0.169 | 12.33 | 4.35 | 0.38 |
| Admin-1 anchor (hard, split-sample) | 24 | 0.321 | 0.249 | 10.47 | 4.09 | 0.55 |

Paired against the unanchored arm on identical cells:

| Arm | delta r | better on r | delta MAE pp | better on MAE | delta abs bias pp |
|:---|---:|:---|---:|:---|---:|
| Flat regional mean | **+0.360** | 19/24 | **-3.39** | 20/24 | **-2.44** |
| Admin-1 anchor (hard) | +0.253 | 22/24 | -1.86 | 18/24 | -1.60 |
| Admin-1 anchor (shrunk) | +0.154 | 16/24 | -1.35 | 15/24 | -0.51 |
| Admin-1 anchor (hard, split-sample) | +0.165 | 18/24 | -0.30 | 14/24 | +0.88 |
| National anchor | +0.005 | 15/24 | +1.30 | 4/24 | +2.64 |
| **Admin-1 anchor (hard, JACKKNIFE)** | **-0.001** | **8/24** | **+1.56** | 8/24 | +1.14 |

The flat national mean has no defined correlation: it predicts one value for
every district, so its `r` is undefined by construction rather than poor. Only
its MAE and bias are interpretable, and both are worse than the flat regional
mean.

## WS2a. The controls Section 4.4 lacked

Section 4.4's counterfactual is the national anchor: if the gain were leakage,
the national anchor would gain too, and it does not. That argument fails because
the two anchors differ in resolution as well as in how much of a district's own
data they contain. A national anchor built from 1,000 respondents contains about
0.1 percent of any one district's data; a regional anchor built from three to
four districts contains a quarter to a third. The absence of a gain at national
resolution is what a leakage account predicts, not evidence against it.

The correct control is an arm at the **same resolution** that uses **no
covariates**. The flat regional mean is that arm, and it wins.

## WS2b. The jackknife, and what the split-sample adds

The jackknife recomputes, for each district, its region's design-based
prevalence from the region's other districts, solves the logit shift on those
other districts, and applies that shift to the held-out district. No quantity
entering a district's prediction has seen that district.

The result is the finding above: the gain vanishes. Both level and pattern
metrics agree, which matters because Section 4.4 concedes `r` may be the wrong
summary. Under the jackknife, MAE and absolute bias are both **worse** than no
anchor, so no reading of the metrics rescues the arm.

The split-sample arm retains about two thirds of the `r` gain (+0.165 against
+0.253). It is a **partial** control and should be read as one: splitting
clusters within a region leaves half of the scored district's own clusters in
the anchoring half, so a district still contributes to its own anchor, at half
weight. That it sits between the hard anchor and the jackknife is consistent
with the leakage account and is not independent evidence for it.

## WS2c. The implied-shift audit, and a defect it uncovered

For each cell the un-anchored population-weighted national aggregate was
compared against the design-based survey national estimate. Section 3 reports
the pipeline recovers national prevalence to a mean absolute error of 0.96 pp.

**The mean absolute gap is 9.60 pp over 24 cells**
(source: `results/tables/anchor_implied_shifts.csv`, rows with
`arm == "national"`). The largest are:

| Country | Outcome | Un-anchored aggregate | Survey estimate | Gap pp |
|:---|:---|---:|---:|---:|
| Sierra Leone | child_vitA | 0.896 | 0.120 | **+77.57** |
| Sierra Leone | women_folate | 0.917 | 0.792 | +12.52 |
| Ghana | women_folate | 0.292 | 0.546 | -25.42 |
| Malawi | child_iron | 0.256 | 0.103 | +15.29 |
| Gambia | child_iron | 0.506 | 0.375 | +13.04 |

Sierra Leone child vitamin A is predicted at **89.6 percent** against a survey
estimate of **12.0 percent**. This is not a level offset; it is a broken fit. The
mechanism is visible in the design: Sierra Leone has 14 districts in 4 regions,
so a leave-one-region-out fold trains a 20-covariate ridge on about 10 districts
and extrapolates to the held-out region.

**This does not contradict Section 3's 0.96 pp claim**, because that claim is
about the individual-level SuperLearner's national estimate
(`national_estimates_all.csv`), not about the area-level leave-one-region-out
ridge that the anchoring arms use as their base model. The two are different
estimators and the documents do not distinguish them.

It does change what Section 4 measures. The anchoring arms are correcting a base
model whose level is wrong by a mean of 9.6 pp and by 77 pp in the worst cell.
A large apparent benefit from anchoring is what any repair of a badly mis-levelled
model would produce, and it says nothing about whether regional information is
valuable relative to a competent base model. The flat regional mean makes the
same point from the other direction: it has no base model at all and beats them.

**Section 4.2's claim that national anchoring is worse than useless was withheld
pending this audit and can now be released, with its explanation changed.**
National anchoring does make absolute bias worse (+2.64 pp against no anchor,
better on MAE in only 4 of 24 cells). The reason is not that a single national
number displaces districts that were already correct, as Section 4.2 states. It
is that the base model's level error varies by region, and a single national
shift cannot correct a region-varying error; it moves every district by the same
amount and leaves the regional structure of the error intact.

## WS2d. Hard against shrunk, under the jackknife

Under the published anchors, hard beats shrunk (0.409 against 0.310), which
Section 4.3 reports. Under the jackknife the hard anchor is not better than no
anchor at all, so the hard-against-shrunk comparison no longer identifies a
preferable estimator. The correct statement is that the choice between hard and
shrunk anchoring is not a live question until an anchoring scheme is found that
beats the unanchored baseline out of sample.

## Replacement interpretation for Section 4, for the claims register

The following is drafted for the later manuscript edit and is entered in
`docs/findings/CLAIMS_REGISTER.md`.

> A district model's level can be corrected by the survey's own regional
> estimates, and doing so improves apparent accuracy substantially. That
> improvement is not evidence that the covariate model contributes: assigning
> every district its region's survey estimate, with no covariates at all,
> outperforms the anchored covariate model on correlation, absolute error and
> bias (0.516 / 7.38 pp / 0.77 pp against 0.409 / 8.91 pp / 1.61 pp over the same
> 24 cells). Nor is the improvement fully out of sample: when each district's
> regional anchor is recomputed from the region's other districts, the gain over
> an unanchored model disappears (mean r 0.147 against 0.156, better in 8 of 24
> cells). What Section 4 establishes is that **the survey's regional estimates
> are the useful quantity**, and that a district map is best read as a regional
> estimate applied to its districts until a covariate model is shown to add to
> one out of sample.

## Reviewer pass, statistical

**Identical cells and folds.** Every arm is produced inside one loop from one
base out-of-fold prediction vector per cell, so the arms differ only in the
anchoring applied afterwards. The paired comparison joins on country and
outcome and reports the pair count for each row.

**Leakage.** The jackknife is the leakage control and is the point of the
workstream. Its own construction was checked for the same defect it tests: the
regional prevalence excludes the scored district's respondents
(`a1_prev_excluding`), the shift is solved on the region's other districts, and
the scored district enters neither. The split-sample arm is explicitly labelled
a partial control because it does not have this property.

**Survey-weight handling.** Regional and national anchors use
`admin1_design_based()` and `national_design_based()`, both survey-weighted. The
jackknife's regional prevalence uses the same weighted mean with a row filter,
so the anchor and its jackknife are the same estimator on different rows.

**Seeds.** Seed 20260905 for the split-sample cluster assignment, 25 splits per
cell. The base ridge fit is deterministic.

**Metric choice.** MAE and signed bias are primary, per the workstream
specification and per Section 4.4's own concession that `r` is the wrong summary
for an arm that moves only the level. The conclusion does not depend on the
choice: the jackknife arm is worse than no anchor on `r`, MAE and absolute bias
simultaneously.

**What this does not test.** Whether a *different* base model, one that is not
mis-levelled by 9.6 pp on average, would show a genuine anchoring gain. WS2c
shows the base model used in Section 4 is badly mis-levelled, and every anchoring
result here is conditional on it. That is a limitation of the published design
which this workstream inherits rather than removes.

## Reviewer pass, reproducibility

**Targets graph.** No new target. The script follows the pattern of
`scripts/covariates/08_admin1_arms.R`, which it re-implements with additional
arms and writes to a new file rather than overwriting `admin1_arms.csv`.

**Additive results.** `anchor_controls.csv` and `anchor_implied_shifts.csv` are
new files. `admin1_arms.csv` is untouched and its frozen copy is in
`results/tables/frozen_2026-09/`.

**Joins.** Scoring joins survey to prediction through `admin2_join_by()`, so
Malawi is keyed on the pair. The unit-count check in WS7a covers the output.

**Stamp targets.** The script reads
`data/covariates/harmonized/predictors_admin2_harmonized.csv`, the same tracked
input `08_admin1_arms.R` reads, and the `_targets_full` store. No untracked
input.

**Paths and determinism.** `here()` throughout, no absolute path, no `setwd()`.
The one stochastic element is seeded.

**Smoke profile.** `PROFILE=smoke` runs Ghana only, six cells, about two
minutes. The full run over 24 cells took about twenty minutes.
