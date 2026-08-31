# WS6. National level model (VMNIS) repairs

## WS6a. The sampling term was wrong by a factor of 4.7

Section 6.3 reports `sd_sampling` 0.816 on the logit scale for the Vitamin A
preschool panel. At a prevalence near 25 percent the delta-method variance of
`logit(p)` is `1 / (n p (1-p))`, so 0.816 implies an effective sample size of
about **13**. National nutrition surveys do not have thirteen respondents.

Three mechanisms, all measured
(source: `results/tables/vmnis_sampling_audit.csv`):

**1. The arithmetic mean of a reciprocal.** `national_noise_ceiling()` computes a
per-survey variance and then takes `v_bar <- mean(v_s)`. The mean of `1/n` is set
by the smallest surveys, not the typical one. The Vitamin A preschool panel holds
**34 surveys under n = 50** and a minimum of **n = 8**, against a median of
**373.5**. Replacing the mean with the median of the same per-survey variances
takes `sd_sampling` from **0.816 to 0.172**.

**2. The prevalence clamp.** `p` is clamped at 0.005, so `p(1-p)` can fall to
0.005 and one survey can contribute up to `301/(n-1)`. The panel holds **37
surveys with prevalence under 2 percent**.

**3. Missing sample sizes.** 80 of the 528 rows have no `Samplesize`. They drop
out of `v_bar` but stay in the `lmer` fit, so the variance components and the
sampling term are computed on different sets of surveys.

A fourth mismatch is worth recording though it is not a defect in the formula:
the ceiling is fitted on **528 un-aggregated VMNIS records** while the LOCO model
is scored on **108 country-year means**. The two describe different populations.

### The corrected ceiling

`national_noise_ceiling_v2()` uses the median per-survey variance, moves the
clamp to 0.02, and requires a recorded sample size for inclusion.

| Panel | Version | Surveys | sd country | sd method | sd resid | sd sampling | r_max report | r_max standardised |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| Vitamin A / preschool | published | 528 | 1.411 | 0.564 | **0.000** | **0.816** | 0.818 | 0.866 |
| Vitamin A / preschool | revised | 448 | 1.436 | 0.383 | **0.703** | **0.171** | **0.869** | **0.893** |
| Zinc / preschool | published | 219 | 1.130 | 0.468 | 0.254 | 0.207 | 0.892 | 0.960 |
| Zinc / preschool | revised | 187 | 1.126 | 0.496 | 0.322 | 0.116 | 0.882 | 0.957 |
| Vitamin A / NPW | published | 93 | 1.232 | 1.996 | **0.000** | 0.409 | 0.517 | 0.949 |
| Vitamin A / NPW | revised | 90 | 1.173 | 2.101 | **0.312** | 0.195 | 0.482 | 0.954 |
| Folate / NPW | published | 184 | 1.764 | 0.950 | 0.442 | 0.331 | 0.849 | 0.954 |
| Folate / NPW | revised | 161 | 1.746 | 0.949 | 0.563 | 0.166 | 0.843 | 0.948 |

(source: `results/tables/national_vmnis_ceiling_revised.csv`.)

**Saturation restated.** Vitamin A preschool: best model r 0.655 against a
revised `r_max_report` of **0.869** is about **75 percent**, not the 80 percent
Section 6.3 states against 0.818. The Section 6.3 correction to the earlier "98
percent saturated" claim stands, and the headroom is larger than it reported.

## WS6b. The boundary fits were not degenerate fits

Section 6.3 flags `sd_resid` at exactly 0.000 in two panels and attributes it to
a degenerate `lmer` fit, concluding `r_max` is untrustworthy there.

**That diagnosis is wrong.** The raw residual variance from `lmer` is non-zero in
all four panels. The reported zero is produced by subtracting the over-large
sampling term from it and flooring at zero. With the corrected sampling term,
`sd_resid` is 0.703 and 0.312 in the two panels concerned, `resid_at_boundary` is
FALSE everywhere, `sampling_exceeds_resid` is FALSE everywhere, and **all four
panels are usable**.

No refit with weakly informative priors was needed. The workstream anticipated
using `blme`, which is not installed; the reframe would have been INLA, which is.
Neither was required, because the problem was in the subtraction and not in the
fit.

## WS6c. Method covariates: not yet computed

Adding assay, adjustment and cut-point metadata as covariates in the LOCO VMNIS
model, and comparing its level error against outcome standardisation, is **not
yet computed**. The variance decomposition already speaks to the question:
method variance is `sd_method` 0.383 for Vitamin A preschool after correction
and **2.101 against a country variance of 1.173** for Vitamin A NPW, so for that
panel method dominates country. Section 6.3's flag that "that panel is measuring
survey methodology more than country" is **confirmed and strengthened** by the
correction.

## WS6d. Bias and error separated, and the degenerate cells marked

| Arm | Cells | MAE pp | Signed bias pp | Absolute bias pp | Cells excl. degenerate | MAE excl. | Signed bias excl. |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Transported, no level | 8 | 5.81 | **-3.14** | 3.35 | 4 | 9.10 | **-5.31** |
| + national null level | 8 | 8.22 | +4.40 | 4.69 | 4 | 10.51 | +4.16 |
| + VMNIS predicted level | 8 | 12.70 | +10.03 | 10.03 | 4 | 19.12 | +15.05 |
| + true national (oracle) | 8 | 5.58 | **+0.35** | 0.73 | 4 | 8.62 | **+0.98** |

(source: `results/tables/national_composition_revised.csv`.)

**Four of the eight cells are near-degenerate.** The women's vitamin A cells sit
at 1.3 to 2.5 pp true national prevalence, so there is almost no level to get
wrong and an arm that fixes the level can barely improve on one that does
nothing.

**The restated oracle result.** On the four non-degenerate cells, an oracle
national level takes MAE from **9.10 to 8.62 pp**, a gain of 0.48 pp, while
taking signed bias from **-5.31 to +0.98 pp**. The correct statement is that a
perfect national level **removes almost all of the bias and almost none of the
error**. The residual error is pattern, and no national quantity can reach it.
Section 6.4's conclusion survives; its framing conflated two things that move
very differently.

## Reviewer pass, statistical

**Like-for-like.** The published and revised ceilings are computed on the same
panels from the same source file in the same run, differing only in the three
changes named. The revised row uses fewer surveys because rows without a
`Samplesize` are excluded; that is one of the three changes and is visible in the
`surveys` column rather than hidden.

**A limit on the median.** Replacing the mean of per-survey variances with the
median is defensible because the quantity wanted is a typical survey's sampling
variance, not an average over a set dominated by tiny ones. It is still a summary
of a heterogeneous set, and a fully correct treatment would weight each survey's
contribution to the variance decomposition by its own precision rather than
subtracting a single pooled term. That is not done here and is recorded as a
limitation.

**Seeds.** No stochastic call; `lmer` with REML is deterministic given the data.

**The composition table is unchanged in its inputs.** WS6d re-presents
`national_composition.csv` with bias separated and degenerate cells marked. The
arm MAEs reproduce the published values exactly (5.81, 8.22, 12.70, 5.58), which
is the check that the re-presentation did not move anything.

## Reviewer pass, reproducibility

**Targets graph.** No new target. `national_noise_ceiling_v2()` is defined inside
`scripts/accuracy_impact/ws6_vmnis_repairs.R` rather than replacing
`national_noise_ceiling()` in `R/national_vmnis.R`, so no committed result
changes and the two can be compared side by side. Promoting it into `R/` is a
follow-up decision for the project owner, not something this branch did
silently.

**Additive results.** `vmnis_sampling_audit.csv`,
`national_vmnis_ceiling_revised.csv` and `national_composition_revised.csv` are
new files. `national_vmnis_ceiling.csv` and `national_composition.csv` are
untouched, and the WS9 regression gate confirms it.

**Stamp targets.** Reads `data/national/vmnis_national_rep.rds`, a tracked file.

**Paths.** `here()` throughout; no absolute path; no `setwd()`.
