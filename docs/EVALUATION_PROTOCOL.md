# Evaluation protocol

Rules for submissions to the shared district-prediction benchmark. Written for
the Proxy Modeling Alliance. The harness that enforces them is in `harness/`.

The rules exist because this project has measured what happens without them, on
its own data, more than once.

---

## 1. What is being predicted

The **district (Admin-2) population prevalence** of a biomarker-defined
deficiency, as defined in `docs/TARGET_ESTIMAND.md`. Submissions are scored
against the survey's design-based district estimate, which is a noisy estimate
of that quantity. Section 3 of `docs/TARGET_ESTIMAND.md` sets out which ceiling
bounds which estimator, and why a ratio to a ceiling is not comparable across
estimators that have and have not seen the survey.

## 2. The files

| File | Contents |
|:---|:---|
| `harness/folds_loro.csv` | region-blocked leave-one-region-out assignment per district |
| `harness/folds_loco.csv` | leave-one-country-out assignment |
| `harness/targets_public.csv` | the district yardstick for **training** cells only |
| `harness/heldout_cells.csv` | the locked test cells, without their values |
| `harness/heldout_regions.txt` | which regions were drawn, and the seed |
| `harness/ceilings.csv` | the validated empirical reliability per cell |

The held-out values are not published. Submit predictions for the cells in
`heldout_cells.csv` and run the scorer.

**The test set is 13 regions across 4 countries, 30 percent of each country's
regions with a floor of two, drawn under seed 20260910 among regions holding at
least three districts.** Held-out cells carry 9 to 17 districts, median 17.

An earlier design held out **one** region per country. It was discarded after
the acceptance test showed it unscorable: one region gives three or four
districts, and any predictor that is constant within a region has zero variance
there, so every correlation returned NA. That included the best arm this project
has found. A test set that cannot score the leading method is not a test set.

## 3. The rules

**3.1 Prescreening happens inside the fold.** Any covariate selection, ranking,
filtering, imputation or scaling that uses the outcome must be recomputed within
each training fold.

*Measured cost of breaking this rule:* two otherwise-identical pipelines, both
holding out whole regions, differing only in whether the covariate prescreen saw
the held-out regions, differ by a mean of **+0.182 in rank correlation, positive
in all 20 measurable cells**, reaching +0.51 (source: Section 2.3 of
`docs/SESSION_FINDINGS_FOR_REVIEW.md`; the honest-against-optimistic scheme
comparison is in `results/tables/frozen_2026-09/corrected/cv_honesty_compare.csv`).
That figure is **carried forward and has not been re-verified on this branch**;
it is cited here as the project's own recorded measurement, not as a
re-established result.

**3.2 No access to held-out cells, for any purpose.** Not for scaling, not for
choosing a hyperparameter, not for deciding which cells to submit. The scorer
**refuses** a file containing predictions for training cells rather than
filtering them out, because a submitter who has predictions for training cells
may have used the test cells too, and silently dropping the extra rows would
conceal that.

**3.3 Report reliability-adjusted metrics.** `r_share` is `r` divided by the
**empirical** ceiling in `harness/ceilings.csv`, suppressed where that ceiling
is below 0.05.

*Do not use an analytic ceiling.* The analytic estimator used earlier in this
project understates the attainable correlation by a mean of **0.161** and
reports exactly zero in **10 of 24** cells against **3 of 24** for the empirical
one (source: `results/tables/reliability_simulation.csv` and
`results/tables/reliability_analytic_vs_empirical.csv`).

**3.4 Report signed bias separately from absolute error.** They answer different
questions and this project has conflated them: an oracle national level in
Section 6 takes signed bias from -5.31 to +0.98 pp while moving MAE only from
9.10 to 8.62 pp (source: `results/tables/national_composition_revised.csv`,
excluding the four near-degenerate cells).

**3.5 Declare what the model saw.** In particular, whether it used the survey's
own regional or national totals. An arm that anchors on the survey's regional
estimate is not comparable with one that does not, and the ceiling does not
bound it (`docs/TARGET_ESTIMAND.md` section 3).

**3.6 Beat the covariate-free baseline or say you did not.** Assigning each
district its region's design-based survey estimate, with no covariates, scores
mean r **0.516**, MAE **7.38 pp** and absolute bias **0.77 pp** across the 24
cells (source: `results/tables/anchor_controls.csv`, arm
`2a flat REGIONAL mean (no covariates)`). It beats every covariate arm this
project has tested. It is the number to beat.

## 4. Metrics returned

`r`, `spearman`, `mae_pp`, `bias_pp` (signed), `r_share` against the empirical
ceiling, and `topk_capture`, the share of the truly worst quarter of districts
that the submission also ranks in its worst quarter.

## 5. Reference result

The project's own covariate-free regional arm, scored through this harness on
the locked cells: **24 cells, median r 0.336, mean MAE 7.19 pp, mean signed bias
-0.13 pp, median `r_share` 0.669 on the 21 cells with a usable ceiling.**

## 6. How to run it

```bash
Rscript harness/score_predictions.R my_predictions.csv results.csv
```

The prediction file needs `country`, `outcome`, `Admin1`, `Admin2`, `pred`, with
`pred` a proportion in [0, 1]. The scorer exits with status 2 and prints the
offending rows if a rule in section 3.2 is broken.
