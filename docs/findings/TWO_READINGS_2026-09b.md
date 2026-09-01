# Two readings of the same results — 1 September 2026

This supersedes the late-September revision, which is preserved in git history.
Both abstracts describe the **same** measurements. Neither contains a claim the
other contradicts. Every number is now produced under **one protocol**:
region-blocked folds, a null that is the mean of the training areas only, and
the shared 373-predictor set across 19 domains.

That protocol matters more than any modelling choice in the project. Fold
construction alone moves results by 0.262, and the choice of null moves them by
2.80 pp — enough to flip the central conclusion on its own.

---

## Reading A: the optimistic abstract

**Micronutrient surveys can be extended for targeting, and the extension
transports to countries that have never run one**

*Background.* Sub-national targeting of micronutrient programmes is limited by
survey cost. The standing hope has been that geospatial covariates could
interpolate deficiency between surveyed districts, and progress against that
hope has been hard to judge because estimators were scored against yardsticks
whose noise was never quantified.

*Methods.* Across 24 country–outcome cells in four national surveys we measured
the reliability of the district yardstick by split-half resampling, fixed a
single evaluation protocol, and scored two outcome types — deficiency
prevalence and the underlying biomarker concentration — under region-blocked
folds against an out-of-sample baseline. Predictors were assembled into a
shared set of 373 covariates spanning household surveys, soil chemistry, crop
maps, malaria layers and satellite climate and vegetation, each carrying
coverage and domain metadata. Targeting performance was scored as the
deficiency burden reached by visiting the worst-ranked fifth of districts.

*Results.* **Targeting works, and it transports.** Ranking districts by
predicted burden reaches **1.23 times** the burden of untargeted allocation
when a model is applied to a country held out of its training entirely, above
chance in **12 of 16** cells — statistically indistinguishable from the 1.21
measured within a country. Correlation over the same three protocols runs
0.416, 0.154 and 0.184, so the ranking survives exactly the fold change that
destroys the correlation. **Modelling the biomarker concentration rather than
the deficiency indicator improves prediction in 18 of 23 cells**, median gain
+0.114 and median r 0.206 against 0.058, with fewer sign reversals. A
decomposition that takes the level from a national anchor and the ranking from
covariates cuts cross-country error from 17.5 to **8.9 pp** and bias from 14.2
to **1.0 pp**. And roughly a fifth of a survey, spread across regions, recovers
district accuracy within 1 pp of the full survey in three of four countries.

*Conclusion.* The deliverable is a **ranked priority list with an anchored
level**, and it is available now. It does not require the target country to
have been surveyed, provided a national estimate exists. The covariate
programme's limits are what made this product identifiable.

---

## Reading B: the pessimistic abstract

**Geospatial covariates do not support sub-national micronutrient estimates,
and the failure survives every attempt to rescue it**

*Background.* A substantial literature reports that environmental and remotely
sensed covariates predict micronutrient deficiency sub-nationally. Those
reports are almost always validated against the same survey's estimates, under
folds permitting spatial leakage, against baselines that contain the held-out
unit, and without correction for the number of predictors screened.

*Methods.* As above, with the emphasis on falsification, and with five distinct
rescue attempts pre-specified: more predictors, mechanistically chosen
predictors, a continuous rather than dichotomised outcome, a stacked estimator,
and a model-predicted national anchor.

*Results.* Under region-blocked folds against an out-of-sample baseline, median
correlation is **0.058** for prevalence and **0.206** for the biomarker level,
and neither beats the baseline in more than half the cells. Every rescue
failed. Of **294** harmonised predictors, **two** survive family-wise
correction; of **31** chosen in advance for mechanistic proximity to intake,
**none** do. **No domain is load-bearing** — the median cost of deleting an
entire domain is **exactly zero**, deleting one *helps* in 72 percent of cases,
and **every single domain performs better alone than all 373 together**, which
is overfitting stated without euphemism. The continuous outcome's correlation
advantage **does not reach the decision**: on classifying districts into
public-health-significance bands it wins in only 3 of 10 cells, and **both
outcome types lose to assigning every district the commonest band** in 6 of 10,
by margins reaching 0.30 for child vitamin A. With targeting scored, the
**unshrunken survey estimate outranks every model** at all covariate-signal
levels this project has measured. A model-predicted national anchor is **worse
than assuming the regional average** (17.98 against 13.70 pp).

*Conclusion.* The evidence does not support a covariate-driven micronutrient
map at Admin-2 resolution. The effective sample size is the number of areas —
14 to 87 — not the number of people, and no predictor set repairs that.
Published sub-national accuracies should be treated as unreliable unless they
report the reliability of their validation target, use administratively blocked
folds, calibrate multiplicity against a permutation null, **and state whether
their baseline saw the held-out unit**. Applied here, those four conditions
remove the effect.

---

## What moved since the previous version

**1. The optimistic reading's newest card does not reach the decision.**
Biomarker levels predict better than prevalence in 18 of 23 cells. That
advantage **disappears** when the outcome is which band a district falls in: 3
of 10 cells, median −0.011. A statistical gain that does not survive
translation into the decision it is meant to inform is worth less than its
effect size suggests.

**2. Both outcome types lose to the modal band.** In 6 of 10 cells with a
verified band, assigning every district the commonest band beats the model —
Gambia child vitamin A 0.167 against 0.467, Malawi 0.264 against 0.448. This is
the null comparison applied to the decision the bands actually drive.

**3. `eb_stack` is withdrawn.** It won error at every measured signal level and
**wins targeting lift at none**. When the endpoint moved from estimating
prevalence to ranking districts, the estimator selected under the old objective
stopped being the right answer. The direct district estimate wins lift at 1.72
to 1.74.

**4. No domain is load-bearing, and the full set is worse than its parts.**
Median drop cost exactly zero; every domain beats all 373 together.

**5. The WHO bands are mostly not WHO's.** Verified against source documents:
only the vitamin A bands for **preschool children** are WHO's for the
population they are applied to. The iron bands are WHO's **anaemia**
classification measured by **haemoglobin** (Table 3, p.17) applied to
ferritin-based iron deficiency — a different biomarker and a different
condition. Folate and B12 have **no WHO banding at all**; the anaemia bands were
reused. Zinc has one real threshold (IZiNCG, 20 percent); the 10 and 30 percent
bands have no source.

**6. Every outcome uses an individual-specific cut-off.** Vitamin A's is
age-specific, zinc's is age, sex and time-of-day specific. Measured, 6.5 to 13.4
percent of non-deficient children sit below the highest deficient
concentration. **No district-level concentration threshold recovers any
binary**, so "predict the level, then threshold" cannot work at district
resolution. An out-of-fold isotonic level-to-prevalence map is used instead.

**7. Two cells are withdrawn rather than reported.** Sierra Leone women B12 has
**four deficient women in 768**, district prevalence SD 0.0066 — there is
nothing to predict, and its −0.730 correlation is a model fitting noise. Malawi
zinc's binary cannot be recovered from its continuous variable for the reason in
(6).

**8. A predictor-breadth claim of mine was wrong.** Gambia's GEE extraction was
doubled from 222 to 435 columns. Every downstream result was **bit-identical**.
The binding constraint on the vocabulary is the harmonisation rules file, not
any country's data.

## Limitations that bias the comparison

**Against Reading A.** The simulator used to rank estimators contains **no
regional block structure**, while about 40 percent of real district variance
sits between regions, so shrinkage estimators are handicapped there relative to
reality.

**Against Reading B.** The band-agreement result rests on **10 cells**, and the
child vitamin A cells that fail hardest have four bands over as few as 14
areas, where a modal-band baseline is unusually strong. And Gambia beats the
null in 4 of 4 cells on both outcome types, so the negative is a statement about
the pooled median, not about every country.

**Against both.** Four countries is a small n for any cross-country claim, and
transportability is not estimable at that n. WHO threshold citations are now
verified, and the verification is what demoted most of them.
