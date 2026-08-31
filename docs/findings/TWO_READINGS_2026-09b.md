# Two readings of the same results — late September 2026 revision

This supersedes `TWO_READINGS.md`, which is kept unchanged. That version was
written before the estimator tournament, the predictor-consistency
meta-analysis, the sixteen-cell arm comparison and the nutrition-proximal scan.
Three of its claims have moved, and they are listed at the end.

Both abstracts describe the **same** measurements. Neither contains a claim the
other contradicts. They differ in which measured facts they place first and in
what they take the project to be for. Every number is source-mapped in
`docs/findings/CLAIMS_REGISTER.md`.

---

## Reading A: the optimistic abstract

**A shippable district product exists, and it is the one that admits how little
the covariates know**

*Background.* Sub-national targeting of micronutrient programmes is limited by
survey cost. The standing hope has been that geospatial covariates — climate,
soil, accessibility, remotely sensed agriculture — could interpolate deficiency
prevalence between surveyed districts. Progress against that hope has been hard
to judge because every estimator was scored against the survey's own district
estimates, a yardstick whose noise was never quantified.

*Methods.* Across 24 country-outcome cells in four national surveys (Gambia,
Ghana, Malawi, Sierra Leone) we measured the reliability of the district
yardstick by split-half resampling, then built a simulation in which the true
district prevalence is known by construction and covariates carry a controlled
share ρ of its variance. Six estimators were scored against that truth at four
levels of ρ, paired within replicate. We tested 294 harmonised predictors for
sign- and magnitude-consistency across cells by random-effects meta-analysis
with permutation calibration, added 31 nutrition-proximal predictors covering
soil zinc, crop mix, helminth burden and programme coverage, and compared three
nested individual-level arms under region-blocked folds.

*Results.* The district yardstick has a median reliability of **0.613**, so
roughly forty percent of its variance is sampling noise, and every prior
estimator comparison was scored against a target that partly rewards
reproducing that noise. Against simulated truth an **empirical Bayes blend** —
each district's own estimate shrunk toward a jackknifed regional mean, with the
shrinkage weight set by the measured reliability rather than by a moment
estimator — attains the **lowest error of any estimator at every level of
covariate signal**, 6.05 pp mean absolute error at ρ 0.35 against 6.84 pp for
the unshrunk district estimate and 8.75 to 9.28 pp for every covariate arm. It
uses no predictor at all. We ship it: **19 of 24 Admin-2 layers released**, five
suppressed where reliability is below 0.30, six flagged by a calibration gate,
each district carrying a bootstrap rank interval. Field haemoglobin, where a
survey has already collected it, adds **+0.174** to district ranking and the
gain is specific to iron outcomes (+0.758, +0.697, +0.394) and absent for
vitamin A, which is the pattern physiology predicts. A previously blocked
individual-level linkage for Malawi was solved this session: the survey's
cluster, household and line identifiers match the DHS household roster on
**3097 of 3099 records**, admitting a fifth country's worth of questionnaire
data to the individual-level arms.

*Conclusion.* The project's deliverable is real, defensible and available now.
It is an estimator, not a model: the honest product of this evidence is a
shrunken district prevalence with explicit rank intervals and explicit
suppression, and it beats every covariate alternative on the metric that matters
for targeting. The covariate programme's failure is itself the finding that made
the shippable product identifiable.

---

## Reading B: the pessimistic abstract

**Geospatial covariates do not predict micronutrient deficiency at district
resolution, and four independent lines of evidence now say so**

*Background.* A substantial literature reports that environmental and remotely
sensed covariates predict micronutrient deficiency sub-nationally. Those reports
are almost always validated against the same survey's district estimates, under
folds that permit spatial and administrative leakage, and without correction for
the number of predictors screened.

*Methods.* As above, with the emphasis on falsification: permutation-calibrated
family-wise correction across predictors, leave-one-region-out and
leave-one-country-out folds, decontaminated and strictly nested predictor arms,
and a pre-specified mechanistic hypothesis (soil and agronomic determinants of
dietary iron) tested as a named nine-covariate bridge rather than as a screen.

*Results.* Of **294 harmonised predictors**, **two** survive family-wise
correction for consistent association across cells. Adding **31
nutrition-proximal predictors** chosen for mechanistic proximity to intake —
soil zinc, crop composition, vitamin A capsule and deworming coverage, iron
supplementation in pregnancy — yields **zero survivors** under permutation
calibration in every family and cell set tested; the two to six predictors
flagged by the analytic false-discovery rate are all extinguished once the
correction is calibrated against the data's own null. The pre-specified iron
bridge attains **0.027** under leave-one-region-out and **−0.022** under
leave-one-country-out. In a simulation where covariates are *guaranteed* to
carry signal, covariate estimators lose to a covariate-free shrinkage rule at
every level of that signal, and carry a systematic **−3.1 to −3.3 pp** bias; a
sixth arm that shrinks toward the covariate prediction rather than the regional
mean buys correlation at low signal and imports that bias intact. The published
questionnaire gain of **+0.099** falls to **+0.036** on the same twelve cells
once the arms are decontaminated and nested — roughly two thirds of it was
leakage and non-nesting. The apparent gap between the project's optimistic and
pessimistic internal validation numbers is **+0.262** and is attributable to
fold construction alone, with no difference in learner or predictor set. The
released district rankings have a median rank-interval width of **42 to 66
ranks** in the 87- and 75-district countries.

*Conclusion.* The evidence does not support a covariate-driven micronutrient
map at Admin-2 resolution in these settings. The effective sample size is the
number of districts, not the number of individuals, and at 14 to 87 districts
per country the design cannot support the predictor counts routinely screened.
The one product that survives falsification uses no covariates. Published
sub-national micronutrient prediction accuracies should be treated as
unreliable unless they report the reliability of their own validation target,
use administratively blocked folds, and calibrate multiplicity against a
permutation null — three conditions that, applied here, remove the entire
effect.

---

## What moved since the first version

Three claims in `TWO_READINGS.md` no longer hold as written.

**1. The questionnaire gain does not reverse sign.** The first version reported
it reversing to **−0.013**. That was measured on a four-cell parity subset. On
all twelve computable cells it is **+0.036** — much reduced from the published
+0.099, but positive. The four-cell subset was not representative.

**2. The blend does not overtake the direct estimate on ranking.** A one-cell
smoke run suggested it did at high covariate signal. The eight-cell run does not
reproduce that: the direct estimate wins on correlation at every ρ. The blend's
advantage is on **error**, not on ranking, and the abstracts above say so.

**3. Malawi's individual-level linkage is no longer blocked.** The first version
recorded it as requiring a codebook request. The crosswalk was found and
verified this session at 3097 of 3099 records. Constructing the questionnaire
block from the DHS recodes remains outstanding data engineering, but the
blocking uncertainty is resolved.

## A limitation that biases the shipped comparison

The simulator generates truth as `sqrt(rho) * z + sqrt(1-rho) * e` with **no
explicit regional block structure**, while about **40 percent** of real
district-level variance sits between regions. An estimator that shrinks toward a
regional mean is therefore handicapped in the simulation relative to reality.
This biases *against* the shipped blend, so Reading A's central result is
conservative and Reading B's is not affected. Adding a region random effect to
the generator is the fix and is not done.
