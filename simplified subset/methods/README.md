# Methods toolkit (simplified subset)

Runnable, self-contained implementations of three of the priority methods to-dos,
working on the simplified-subset data. Each is a teaching/prototyping reference
Andy can run and adapt; all build on `_helpers.R` (data loading, the 16-predictor
elastic-net area learner, and metrics). See
[`ANDY_KIM_PROJECT_PLAN.docx`](../ANDY_KIM_PROJECT_PLAN.docx) and
[`../../docs/METHODS_TODO_IMPLEMENTATION_PLAN.md`](../../docs/METHODS_TODO_IMPLEMENTATION_PLAN.md)
for the full plan and the honesty flags.

## Run (from the repo root)

```
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla "simplified subset/methods/national_anchor_calibration.R"
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla "simplified subset/methods/binomial_loss.R"
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla "simplified subset/methods/aggregate_uncertainty.R"
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla "simplified subset/methods/aggregate_inference.R"
```
Needs `glmnet`. Each script has `OUTCOME` (and other knobs) at the top.

## What each script shows (reference results, women_iron)

- **`national_anchor_calibration.R`** (to-do 2). Leave-one-country-out, then a
  one-parameter logit shift so the predicted national mean matches a known
  anchor. Reference run: mean MAE 20.9 to 9.6 pp, mean |bias| 19.5 to 2.1 pp,
  Spearman ranking unchanged. Standard SAE benchmarking; low novelty.
- **`binomial_loss.R`** (to-do 3b). Within-country 5-fold CV comparing a
  weighted-Gaussian area model with a binomial-counts model. Reference run: RMSE
  10.9 (Gaussian) vs 10.5 pp (binomial), and the binomial keeps every prediction
  in [0,1] by construction. Standard; a default-setting change. (Note the glmnet
  two-column response is (failures, successes), so deficiency counts go second.)
- **`aggregate_uncertainty.R`** (to-do 3c, uncertainty part). Split-conformal
  per-district intervals (reference coverage ~0.96 at the 90% target) plus an
  area bootstrap CI for the national aggregate, and a heuristic
  risk-times-uncertainty value-of-information ranking. Honest result worth noting:
  the national bootstrap CI is narrow and usually misses the truth under
  transport, because it captures sampling noise but not the transport level bias
  that the national-anchor calibration corrects.
- **`aggregate_inference.R`** (to-do 3c, deep dive). The cheap, defensible pieces
  of loss-based aggregate inference: a country-block bootstrap CI on the
  cross-validated aggregate error (Option A, honest and wide at ~4 effective
  units) and a partial-identification band for an unsurveyed country's national
  prevalence built from the across-country range of level-shifts (Option D).
  Reference run: between-country national error 18.5 pp (90% CI [12.4, 24.5]);
  the partial-ID band covers the truth in 2/4 held-out countries versus 1/4 for
  the naive bootstrap. See the To-do 3c deep dive in
  `../../docs/METHODS_TODO_IMPLEMENTATION_PLAN.md` for the full option set
  (including the TMLE / transport-EIF routes that need statistician sign-off).

- **`qa_reporting.R`** (standard output QA, prompted by the WFP/MIMI materials).
  Two reusable, post-prediction checks that attach to any surface:
  `skill_over_baseline()` reports a model as improvement over a naive baseline
  (WFP's "normalized difference vs a dummy") — transport skill vs the training
  mean, and ranking skill vs the country's own mean after removing the level
  bias; and `gradient_sanity_check()` tests whether a predicted surface
  reproduces an expected structural gradient (deficiency higher in the
  less-developed tertile), returning a per-country pass/fail QA gate for the
  data-constrained case with no ground truth. Reference run: under LOCO the
  transport skill vs the training mean is mostly <= 0 for the women's outcomes
  (the effective-n story as one number) and partially positive for child_iron;
  gradient pass rates are 25-75% and correctly flag near-flat transported
  surfaces. Standard reporting conventions; no methodological novelty, but they
  make the honesty story legible. Natural to lift into `R/corrected/` as a
  post-prediction QA step.

## Honesty flags (short version)

- Calibration, binomial loss, split-conformal, and the bootstrap are standard and
  packaged. The VOI ranking is an explicit heuristic, not formal optimal design
  (EVSI). The full novelty discussion, including what needs statistician sign-off
  for the other to-dos (decision-focused weighted ensembling, the transport
  version of sequential multi-outcome, loss-based aggregate inference), is in
  `docs/METHODS_TODO_IMPLEMENTATION_PLAN.md`.

## Promotion path

These prototypes operate on the simplified subset for speed and clarity. Once
validated, the calibration and conformal/bootstrap pieces are natural additions
to `R/corrected/` in the main pipeline (the split-conformal already lives there);
the binomial family is a config change for the area-level learners.
