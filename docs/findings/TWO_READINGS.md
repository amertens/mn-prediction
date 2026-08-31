# Two readings of the same results

Both abstracts below describe the **same** measurements from the September 2026
accuracy-and-impact work. Neither contains a claim the other contradicts. They
differ in which measured facts they place first and what they take the project to
be for.

They are written to be set side by side, because the honest position is that the
evidence supports both and the choice between them is a judgement about purpose,
not about data.

Every number is source-mapped in `PROJECT_STATUS_2026-09_UPDATE.md` and
`docs/findings/CLAIMS_REGISTER.md`.

---

## Reading A: the optimistic abstract

**Micronutrient deficiency can be mapped where surveys reach, and the ceiling on
district prediction is far higher than previously believed**

*Background.* Geospatial prediction of micronutrient deficiency has been
constrained by an apparent ceiling: district survey estimates were thought so
noisy that a median correlation above 0.098 was unattainable, implying models
were already near-saturated and further effort was futile.

*Methods.* We revalidated that ceiling by split-half resampling within districts
across 24 country-outcome cells in four national surveys, corrected it against a
simulation in which the true district prevalence was known by construction, and
re-scored every estimator against it. We added covariate-free controls, a
jackknife anchor, a third predictor arm admitting field haemoglobin, and a
survey-design experiment varying the anchoring budget.

*Results.* The ceiling was substantially underestimated. Measured empirically it
is **0.613**, not 0.098; the analytic estimator understates the attainable
correlation by a mean of **0.161**, and reports no recoverable signal in 10 of 24
cells where measurement finds it in only **3**. Most cells therefore contain real
district-level signal that had been written off. Regional estimates are excellent
and cheap: full-survey Admin-1 error is **0.297 pp**, and a covariate-free
regional predictor reaches **r 0.516, MAE 7.38 pp, bias 0.77 pp** at district
level. Field haemoglobin, already collected in every DHS, adds **+0.394** to
women's iron prediction in The Gambia, reaching **r 0.848**. At national level
the VMNIS model reaches **75 percent** of a corrected ceiling of 0.869, leaving
real headroom. Micronutrient targets sit **above** the skill-reliability line
defined by 81 standard DHS indicators, so they are not unusually hard.

*Conclusions.* A defensible, low-cost deliverable exists now: survey-anchored
regional prevalence, district rankings with suppression where reliability is
below 0.30, and a survey-design rule to buy clusters within sampled regions
rather than wider region coverage. The constraint is not the target.

---

## Reading B: the pessimistic abstract

**Every positive result in geospatial micronutrient prediction failed its own
control, and the covariates predict wealth rather than nutrition**

*Background.* Geospatial covariates are widely proposed for subnational
micronutrient estimation where surveys are absent. We subjected one such
pipeline, spanning 294 harmonised Admin-2 predictors and four national surveys,
to the controls its positive findings had not previously faced.

*Methods.* We added covariate-free baselines at matched resolution, a jackknife
anchor excluding each district from its own regional target, a nested three-arm
predictor design, an automated leakage report, and a label-derived assay guard
built from original Stata files. We ran the identical pipeline on 81 standard
DHS indicators as pseudo-targets.

*Results.* The reported anchoring gain does not survive its control: mean r falls
from **0.409 to 0.147** against **0.156** for no anchor, and both error and bias
worsen. A predictor using **no covariates at all**, assigning each district its
region's survey estimate, outperforms every covariate model tested
(**0.516 / 7.38 pp / 0.77 pp**). The reported questionnaire advantage reverses
from **+0.075 to -0.013** once fifteen blood-derived columns are removed and the
arms are properly nested. Across 81 DHS indicators the pipeline achieves a median
correlation of **0.071**; the correlation between target reliability and achieved
skill is **0.043**, so the ceiling framing does not explain performance. The
covariates predict education (**0.679 to 0.808**) and wealth (**0.62**) but not
stunting (**-0.122**), water (**0.014**) or sanitation (**0.007**). Fifteen
leaked columns were found in a codebase that had already fixed three, and one
fitted model predicted **89.6 percent** prevalence where the survey measured
**12.0 percent**.

*Conclusions.* District micronutrient prevalence from geospatial covariates is
not currently supportable. What the covariates recover is a socioeconomic
gradient, and what the district maps recover is the regional survey mean. Absent
a predictor set validated against health outcomes rather than wealth proxies,
the defensible output is the survey estimate itself.

---

## What separates them

| Question | Reading A | Reading B |
|:---|:---|:---|
| Is there district signal? | Yes: ceiling 0.613, only 3 of 24 cells empty | Irrelevant: models reach 0.200, and reliability does not predict skill (0.043) |
| Do the covariates work? | They recover a third of an attainable maximum | They are beaten by having no covariates |
| Is the target the problem? | No, and that is good news | No, and that is bad news: the predictors are |
| What should ship? | Regional prevalence plus suppressed district rankings | The survey estimate |

**Both readings agree on the practical output.** Reading A arrives at
survey-anchored regional estimates with cautious district rankings; Reading B
arrives at the survey estimate. The gap between those two deliverables is
narrow, and `docs/FITNESS_FOR_USE.md` sits inside it.

**They disagree about the research programme.** Reading A says the headroom is
real and worth investing in. Reading B says the headroom is real and there is no
evidence any available covariate can close it. Deciding between them needs the
work this branch did not compute: the aggregation-level sweep (WS4a), the full
sixteen-cell protocol comparison, and a covariate set built for nutrition rather
than inherited from poverty mapping.

---

## Correction, 2026-09-21

Both abstracts above were written before two controls were run. Neither abstract
is rewritten; this note records what supersedes them and where.

**1. The covariate-free regional arm was not jackknifed.** Reading B quotes it at
`0.516 / 7.38 pp / 0.77 pp` and calls it the arm that "outperforms every
covariate model tested". Reading A quotes the same figures. Those numbers were
produced with each district contributing to the regional estimate used to
predict it, which is the same mechanism that withdrew the hard-anchor claim
(register row 4.6). Under a symmetric jackknife, where a district's regional
estimate is computed from the region's other districts, the arm scores

| Arm | mean r | MAE pp | mean absolute bias pp |
|:---|---:|---:|---:|
| Flat regional mean, as quoted above | 0.516 | 7.38 | 0.77 |
| **Flat regional mean, jackknifed** | **0.076** | **10.87** | **3.41** |
| No anchor (LORO), for comparison | 0.156 | 10.77 | 3.21 |

(source: `results/tables/anchor_controls_B1.csv`, arms
`2a flat REGIONAL mean (no covariates)` and `2b flat REGIONAL mean (JACKKNIFE)`.)

**The claim that a covariate-free predictor beats every covariate model does not
survive its own control and is withdrawn.** It appears in Reading B's abstract
and conclusions, in Reading A's results paragraph, and in the comparison table's
"Do the covariates work?" row. The corrected statement is that under a matched
out-of-sample control the covariate model and the covariate-free regional mean
are close, with the covariate model ahead on correlation and the two
indistinguishable on error.

This correction was produced by applying to the covariate-free arm the same test
that had already been applied to the anchored arms. The asymmetry was an
oversight in the earlier session, not a difference in the estimators.

**2. The field-haemoglobin result rested on four cells.** Both abstracts cite
`+0.394` for Gambia women's iron and `0.848` for that cell. Those come from a
four-cell subset. The sixteen-cell run supersedes them; see
`docs/findings/WS3_INDIVIDUAL_ARMS.md` and
`results/tables/individual_arms_2026-09_16CELL.csv` when present.

**What does not change.** The reliability findings, the withdrawal of the
anchoring gain, the skill curve across 81 DHS indicators, the resolution sweep,
and the VMNIS repairs are unaffected by either correction.
