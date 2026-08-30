# Implementation plan for the four priority methods to-dos

A detailed, step-by-step plan Andy can follow to implement the four
highest-value methodological to-dos. For each one this gives the goal, a few
brainstormed approaches with a recommendation, a concrete step sequence, the R
packages or existing methods to use, and an explicit honesty flag on what is
standard versus what is methodologically novel and should be validated by a
statistician (Alan) before it goes in the manuscript.

Read the honesty flags carefully. Several of these are well-trodden and packaged;
a few are genuinely novel in this setting and must not be presented as settled.

## Shared setup (do this once)

- **Data.** Work from the simplified subset (`simplified subset/data/mn_*.csv`):
  one row per area, with `prev_<outcome>`, `n_<outcome>`, `ndef_<outcome>`,
  `var_<outcome>` per outcome and 16 shared proxy predictors. Levels are
  `mn_cluster.csv`, `mn_admin2.csv`, `mn_admin1.csv`. The full pipeline data is in
  `_targets_full/objects/` if you need more covariates than the 16.
- **Code to build on.** The corrected methods in `R/corrected/` already implement
  honest cluster/spatial-block CV (`p1_fitting.R`), out-of-fold calibration and
  split-conformal intervals (`p2_p6_methods.R`), sampling-error-aware admin-2
  error, decision-value scoring, and partial-pooling SAE (`p7_area.R`). Reuse
  these rather than re-deriving them.
- **Evaluation discipline.** Always report leave-one-country-out (LOCO) for
  transport and cluster/spatial-block CV for within-country. Never use in-sample
  or random-fold numbers as headline. This is settled in the project.

---

## To-do 1: Decision-focused ensembling / top-N targeting

**Goal.** Stop optimizing average squared error across all districts and instead
optimize the operational task: correctly identifying the highest-burden districts
(above a WHO threshold, or the top quintile). This is your sensitivity-constrained
/ quintile work applied to areas.

### Brainstormed approaches

- **A. Binary high-burden reframing (recommended starting point).** Define a
  binary label `high = prev_<outcome> > threshold` (WHO public-health cutoff) or
  `high = top 20% within country`. Fit the existing SuperLearner as a classifier
  of `high`, and choose the ensemble by a ranking/decision metric rather than NLL.
- **B. Continuous prediction, decision-focused model selection.** Keep predicting
  `prev_<outcome>`, but select the discrete SuperLearner (the single best library
  member) by a cross-validated ranking metric (Spearman, recall@top-k, or
  area-under-the-precision-recall-curve for the high-burden label). This is the
  cleanest mapping of Alan's "do it in the ensembling stage" comment.
- **C. Learning-to-rank.** Fit a model whose objective is directly a ranking loss
  (xgboost `objective = "rank:pairwise"` or `"rank:ndcg"`, grouped by country).
  More direct, but harder to calibrate and less familiar to reviewers.

Recommended path: A and B together. Report recall@top-k and a
sensitivity-at-fixed-specificity curve, within-country and under LOCO, against the
MSE-optimal model as the baseline.

### Step by step

1. On `mn_admin2.csv`, for a chosen outcome (start with `women_iron`), build the
   high-burden label two ways: above the WHO threshold, and top 20% of
   `prev_<outcome>` within each country. Keep both; they answer slightly different
   policy questions.
2. Reuse the corrected SL fit (`fit_corrected_sl` in `R/corrected/p1_fitting.R`)
   to get cross-validated base-learner predictions. The base learners stay as they
   are; only the combination/selection step changes.
3. Define the decision metric as a function of the cross-validated predictions:
   recall@top-k (share of true top-k districts that the model ranks in its top-k),
   subject to a cap on false positives. Also compute Spearman and PR-AUC.
4. **Discrete super learner (low risk):** pick the single library member that
   maximizes the decision metric under LOCO and cluster-block CV. This is a
   standard, defensible use of the super learner (cross-validated risk
   minimization with a custom risk).
5. **Weighted ensemble (higher risk):** if you want a convex combination rather
   than a single winner, optimize the weights to maximize the decision metric.
   Because recall@top-k is non-smooth and non-decomposable, use a derivative-free
   optimizer (`optim` with Nelder-Mead, or a small grid / `GenSA`) on the
   simplex.
6. Evaluate within-country (cluster/spatial-block CV) and LOCO, reporting the
   decision metrics against the MSE-optimal ensemble.

### R packages / existing methods

- `SuperLearner`: has `method.AUC` for AUC-optimized weights and
  `method.template` to define a custom metalearner risk. `CV.SuperLearner` for
  honest evaluation. The project also uses `mlr3superlearner`/`mlr3`, which let
  you tune to any `mlr3` measure.
- `xgboost` / `lightgbm`: built-in ranking objectives for approach C.
- `precrec` or `yardstick` for PR-AUC and recall@k.

### Honesty flags

- **Standard / packaged:** classifying a high-burden label, AUC-optimized SL
  weights (`method.AUC`), and selecting the discrete SL by a custom
  cross-validated risk. Learning-to-rank objectives in xgboost/lightgbm.
- **Needs statistician validation:** optimizing the convex ensemble weights for a
  **non-decomposable, area-level, asymmetric** criterion (recall@top-k subject to
  a false-positive bound). The discrete-SL version is on solid theoretical
  footing (van der Laan's oracle results hold for any bounded loss); the
  weighted-combination version for a non-smooth area-level objective is novel
  here. Confirm with Alan that the chosen risk is a valid loss and that the CV
  selection is honest (no leakage of the top-k definition across folds, especially
  when "top 20% within country" is defined using the held-out fold). The
  sensitivity-constraint formulation is your existing area, which helps.

---

## To-do 2: National-anchor calibration

**Goal.** When transporting to a new country, shift the predicted prevalences so
their population-weighted mean equals a known national prevalence (from a past
survey, VMNIS, BRINDA, or DHS anaemia). This fixes the documented level bias
without changing the district ranking.

### Brainstormed approaches

- **A. Logit shift (recommended).** Find one scalar delta so that the
  weighted mean of `expit(logit(p_i) + delta)` equals the national target. Solve
  with `uniroot`. Strictly monotone, so ranking is unchanged.
- **B. Ratio or difference benchmarking.** Multiply (ratio) or add (difference) a
  constant so the weighted mean hits the target. This is classic small-area
  estimation benchmarking. Ratio keeps values in [0,1] only if the target is
  below 1/max; difference can push values outside [0,1] and needs clipping.
- **C. Distribution matching.** If more than the mean is known (for example a
  national prevalence and a rough spread), match with isotonic regression or
  quantile mapping. Usually overkill given we only reliably know the mean.

Recommended: A (logit shift). It is simple, bounded, and rank-preserving.

### Step by step

1. For each held-out country in the LOCO loop, compute the "known" national
   prevalence. For the demo, use that country's own sample-size-weighted mean of
   `prev_<outcome>` (weighting by `n_<outcome>`) as a stand-in for the external
   anchor you would have in practice.
2. Solve for `delta` with `uniroot` so that
   `sum(w_i * expit(logit(phat_i) + delta)) / sum(w_i) == target`, where `phat_i`
   are the transported predictions and `w_i` is area population (or `n_<outcome>`
   as a proxy).
3. Apply the shift to get calibrated predictions. Confirm the ranking metrics
   (Spearman, recall@top-k) are unchanged and the level metrics (MAE, bias) before
   versus after.
4. Report how much of the LOCO level error a single national anchor removes, per
   outcome and per held-out country.
5. Extension for honesty about the anchor: if the national anchor itself has a
   confidence interval, propagate it by repeating the shift across draws of the
   target and widening the prediction intervals accordingly.

### R packages / existing methods

- Base R `uniroot` and `stats::plogis`/`qlogis` for the logit shift.
- `emdi` provides small-area benchmarking functions; `survey::calibrate` /
  `postStratify` implement calibration to known margins. The hand-coded logit
  shift is a few lines and easy to audit.

### Honesty flags

- **Standard / packaged and low novelty.** Benchmarking small-area estimates to a
  reliable aggregate is well established in the SAE literature (ratio and
  difference benchmarking; calibration estimators in survey statistics). This is
  the most ready-to-ship of the four.
- **Minor choices to confirm with a statistician:** logit shift versus ratio
  versus difference (they differ slightly in how they redistribute the
  adjustment), whether to benchmark per outcome, and how to propagate the anchor's
  own uncertainty. None of these are research-novel.

---

## To-do 3a: Sequential multi-outcome modeling

**Goal.** Borrow strength across deficiencies that co-occur. Predict the
best-measured outcome first (iron) and feed the predicted upstream deficiencies in
as covariates for the next outcome, exploiting their covariance.

### Brainstormed approaches

- **A. Regressor chain (the transcript's idea).** Order outcomes (iron, then
  vitamin A, then B12/folate, then zinc). Model k adds the predicted outcomes
  1..k-1 as covariates. Within-country: train on observed upstream, predict with
  predicted upstream. For transport: must use the predicted upstream (no truth in
  the target country). Add a missing-indicator plus imputation when a country did
  not measure an upstream outcome.
- **B. Multivariate Fay-Herriot (recommended for rigor).** Model all outcomes
  jointly with a multivariate area model that estimates the cross-outcome residual
  covariance and borrows strength through it. This is the statistically principled
  version of the same idea and it is packaged.
- **C. Multi-task regularized regression.** `glmnet(family = "mgaussian")` fits a
  joint multi-response lasso that shares variable selection across outcomes. Light
  and easy, captures shared predictors but not residual covariance.

Recommended: do A first (it directly tests the transcript hypothesis and is easy
on the wide data), then B as the rigorous comparison.

### Step by step (approach A)

1. Fix an order by data richness: `women_iron`, `women_vitA`, `women_b12`,
   `women_folate`, then child outcomes, then zinc.
2. Fit outcome 1 from the 16 proxies (reuse the corrected SL).
3. For outcome k: build the design as the 16 proxies plus the cross-validated
   predictions of outcomes 1..k-1, each with a 0/1 indicator for "was this
   upstream outcome measured in this country" and a median imputation where not.
4. Critical for honesty: within-country, the upstream covariate in the test fold
   must be the cross-validated prediction, not the observed value, to avoid
   leakage. For transport, it is always the predicted upstream.
5. Compare, per outcome, proxies-only versus proxies-plus-upstream, separately for
   within-country and LOCO. Expect the clearest gains for the data-poorest
   outcomes (zinc, B12).

### R packages / existing methods

- Approach A: implement with the existing SL plus base R for the chaining. In
  Python the pattern is `sklearn.multioutput.RegressorChain` / `ClassifierChain`
  if you prefer a reference implementation.
- Approach B: `msae` (multivariate Fay-Herriot small area estimation on CRAN) or
  `emdi`. Approach C: `glmnet(family = "mgaussian")`.

### Honesty flags

- **Established components:** regressor/classifier chains are a known multi-output
  ML pattern; multivariate Fay-Herriot is standard SAE with an R package
  (`msae`); multi-task lasso is standard (`glmnet`).
- **Needs statistician validation:** the **transport** version of the chain, where
  you condition on a *predicted* upstream outcome that itself transports poorly.
  The estimand is not obviously well defined and errors compound. Alan framed the
  sequential idea as a forward factorization of the joint density (the same logic
  as G-computation); confirm with him whether the chained transport prediction is
  coherent or whether the multivariate-FH joint model is the defensible route. Be
  explicit in the manuscript that within-country gains and transport gains are
  different claims.

## To-do 3b: Proper binomial loss for prevalence

**Goal.** Model the 0-to-1 prevalence with a binomial/quasi-binomial likelihood so
predictions stay in range and the mean-variance relationship is respected.

### Step by step

1. Refit the area learners on `prev_<outcome>` with `family = quasibinomial` (or a
   binomial likelihood on `ndef_<outcome>` out of `n_<outcome>`), weighting by
   `n_<outcome>` or by precision `1 / var_<outcome>`.
2. Compare calibration and stability against the Gaussian-link version, especially
   for low-prevalence outcomes (vitamin A, B12 in women).

### R packages / honesty flags

- `stats::glm(family = quasibinomial)`, `betareg`, `mgcv` with `family = betar`,
  or `brms` for a Bayesian binomial. For the SuperLearner, set `family = binomial`
  (most learners accept a non-integer 0-1 response).
- **Low novelty and widely used.** The project already has a quasi-binomial
  benchmark; this is a default-setting change, not new methodology.

## To-do 3c: Aggregate uncertainty and value-of-information map

**Goal.** Put honest uncertainty on the aggregate product (national or admin
prevalence), then rank where a new survey would be most informative.

### Brainstormed approaches (uncertainty)

- **A. Split-conformal (already implemented).** Reuse `intervals_corrected` in
  `R/corrected/p2_p6_methods.R`. Distribution-free, finite-sample valid coverage.
- **B. Cluster/area bootstrap.** Resample clusters within areas (or areas within
  countries), refit, and take percentile intervals. Standard.
- **C. Loss-based cross-validated inference for the aggregate target.** This is
  the van der Laan style of getting a confidence interval for a final aggregate
  measure (the approach Alan referred to). Most rigorous, most involved.

### Brainstormed approaches (value of information)

- **D. Heuristic VOI ranking (recommended starting point).** Rank candidate
  survey sites by `predicted risk x prediction uncertainty`. Simple and
  defensible as a triage heuristic.
- **E. Formal expected value of sample information (EVSI).** Decision-theoretic
  optimal design: the expected reduction in decision loss from sampling a site.
  Principled but heavy.

### Step by step

1. Attach per-area uncertainty using A (reuse the corrected split-conformal) and
   cross-check with B (a cluster bootstrap) on a couple of outcomes.
2. Build the heuristic VOI ranking D: order areas by predicted prevalence times
   interval width (or times posterior variance if you use a Bayesian area model).
   Map it as "where to sample next."
3. Sanity-check that the VOI ranking concentrates on districts that are both
   plausibly high-burden and poorly determined (often out-of-support districts;
   cross-reference `trust_flags` from the corrected methods).

### R packages / existing methods

- Conformal: the project's own split-conformal, or `conformalInference`.
  Bootstrap: `boot` or base resampling. For the formal route, `voi` implements
  health-economic EVPI/EVSI but is not a drop-in for spatial sampling design.

### Honesty flags

- **Standard / packaged:** split-conformal and cluster bootstrap intervals. The
  heuristic risk-times-uncertainty VOI ranking is a reasonable, common heuristic,
  but be explicit that it is a heuristic, not optimal design.
- **Needs statistician validation (genuinely novel here):** (i) loss-based
  cross-validated inference for the aggregate target (approach C) is
  efficient-influence-function / TMLE territory and should be designed with Alan
  or Mark; (ii) formal EVSI / optimal sampling design (approach E) is a real but
  involved decision-theory method that has not been applied in this project.
  Ben Arnold and colleagues' work on sampling for neglected-tropical-disease
  targeting is a useful reference point. Do not present the heuristic VOI as
  optimal design.

### Deep dive: loss-based cross-validated inference for the aggregate target

This is the hardest sub-problem, so it gets its own treatment. The first move is
to be precise about what "the aggregate target" is, because the answer splits
sharply by case and each case needs different machinery:

1. A performance summary (the cross-validated misclassification rate of the
   high-burden label, or the weighted aggregate prediction error). This is a risk.
2. A population mean in a labelled population (national prevalence in a surveyed
   country). Identifiable from data.
3. A population mean in an unlabelled population (national prevalence in an
   unsurveyed country). This is the transport case and is where valid inference is
   genuinely hard, because it is partially or non-identifiable under the very
   assumption (stable outcome model across countries) that the results show fails.

Options, from most established to most novel:

- **A. Cross-validated risk as the parameter (the literal version of Alan's
  comment).** Define the target as the expected loss of the prediction procedure
  for a chosen loss (squared error on prevalence, misclassification of the
  high-burden label, or a top-N recall loss). The out-of-fold loss values are a
  sample; their mean is the CV-risk estimate and an influence-function or CLT
  standard error gives the interval. Established (van der Laan and Dudoit CV-risk
  asymptotics). The catch is the cluster structure: areas are nested in only four
  countries, so a country-clustered SE or country-block bootstrap has about four
  effective units and the between-country interval is wide. Packages: `origami`
  for the CV structure; `boot` for the block bootstrap.
- **B. Efficient-influence-function or (CV-)TMLE for a mean in a labelled
  population (case 2).** Target the population-mean deficiency in a surveyed
  country using the proxy model as a nuisance, with a CI from the efficient
  influence function; CV-TMLE permits SuperLearner nuisances without Donsker
  conditions. Established and packaged (tlverse: `tmle3`, `sl3`, `origami`; also
  `drtmle`, `AIPW`). Caveat: for a surveyed country the design-based survey mean
  already gives a valid CI, so the main value is efficiency and a clean framework
  to extend to case 3.
- **C. Transport EIF / transport-TMLE for the unsurveyed-country aggregate
  (case 3).** Frame the target via the transport formula: the target-country mean
  equals the source-fitted outcome model integrated over the target covariate
  distribution, with EIF-based inference that adds a source-to-target density-ratio
  weight (Dahabreh, Robins, Hernan on transportability; Rudolph and van der Laan
  on transported TMLE). The method exists in causal inference, but applying it
  here is genuinely novel and risky: the identifying assumption (the same
  E[Y given X] across countries) is exactly what the manuscript shows fails, and
  the cross-country density ratio is hard to estimate from four source countries.
  The honest likely conclusion is that the transported point is non-identifiable
  without extra structure, and demonstrating that (an EIF whose covariate-shift
  term explodes) is itself a result. Needs statistician design.
- **D. Partial identification / sensitivity bounds for the transported
  aggregate.** Instead of a point plus a falsely tight CI, report bounds under a
  bounded violation of transportability: let the country-level logit shift lie
  within the range observed across the training countries, and propagate that to a
  band on the target national prevalence. Pairs naturally with the national-anchor
  calibration (if an anchor exists the mean is identified and inference shifts to
  within-country shape; if not, bounds are the honest output). Novel in this
  project; sensitivity analysis for transportability is an active research area.
- **E. Hierarchical Bayesian country random effect.** A multilevel model with a
  country-level random intercept (extending the existing BYM2 / Fay-Herriot work)
  gives a posterior-predictive distribution for a new country's aggregate that
  automatically includes between-country variance estimated from the four
  countries. Packaged and standard (`brms`, `INLA`). Caveat: with four countries
  the random-effect variance is barely estimable, so the prior does much of the
  work and the interval is honest about ignorance rather than frequentist-calibrated.

**Recommendation.** Use A for a risk CI and B for a surveyed-country mean (both
established and packaged). For the transported aggregate that the funder actually
wants, the honest path is D (partial-identification bounds) and/or E (Bayesian
posterior predictive), with C as the rigorous but assumption-heavy comparison.
The binding constraint throughout is that four countries carry almost no
information about between-country variation, so any honest aggregate interval for
an unsurveyed country will be wide, and that is a feature. The specific item to
take to Alan or Mark is whether the transported aggregate is identifiable at all,
and if not, to make the partial-identification framing the deliverable.

**Prototype available.** The cheap, defensible pieces (the country-block bootstrap
CI on the cross-validated aggregate error from option A, and the
partial-identification band from the across-country shift range in option D) are
implemented in `simplified subset/methods/aggregate_inference.R`.

---

## To-do 4: Anaemia validation track

**Goal.** Rerun the whole methodology on DHS haemoglobin (anaemia), which is
measured at large sample size in nearly every sub-Saharan DHS, to validate the
pipeline and transport behaviour on a trusted, abundant outcome and to expand the
cross-country panel far beyond four countries.

### Step by step

1. Acquire DHS haemoglobin data for a set of countries (the project already has
   DHS access; use the `rdhs` R package to pull the recode files programmatically).
2. Define anaemia with the standard WHO haemoglobin thresholds by age, sex, and
   pregnancy status, applying the established altitude and smoking adjustments to
   haemoglobin before thresholding. These adjustment formulas are standard (CDC /
   WHO); no special package is needed.
3. Aggregate to admin-2 with survey weights and design (the same machinery as the
   micronutrient outcomes), and link the same proxy predictors. The country
   ingestion stubs in `R/ingest_new_country.R` are the entry point.
4. Run the existing pipeline (within-country cluster/spatial-block CV and LOCO)
   and the corrected methods (calibration, intervals, decision value, trust).
5. Compare: does the method's within-country and transport behaviour on anaemia
   match what we see for iron? This separates "the method works, the micronutrient
   biomarkers are noisy" from "the method does not work."

### R packages / existing methods

- `rdhs` for DHS data access; `survey`/`srvyr` for design-based aggregation (the
  pipeline already uses these); the existing `_targets.R` pipeline for everything
  downstream.

### Honesty flags

- **Low methodological novelty.** This is applying the existing, validated pipeline
  to a new outcome. DHS anaemia, WHO thresholds, and altitude/smoking adjustment
  are all standard. The work is mostly data ingestion and engineering.
- **Interpretation caveat to state plainly:** anaemia is multifactorial (iron is
  only one cause, alongside malaria, haemoglobinopathies, other deficiencies), so
  it validates the **pipeline and transport methodology**, not micronutrient
  biology. Frame it as a methods-validation and panel-expansion exercise.

---

## Suggested sequence for Andy

1. **National-anchor calibration (to-do 2)** and **proper binomial loss (3b)**
   first: low novelty, packaged, quick wins, and calibration directly addresses
   the documented transport failure. Good for the September Ghana talk.
2. **Decision-focused targeting (to-do 1)**: high value, builds on your existing
   work; do the discrete-SL version first, then discuss the weighted version with
   Alan.
3. **Aggregate uncertainty + heuristic VOI (3c, parts A/B/D)**: reuse the
   corrected conformal intervals; ship the heuristic VOI map.
4. **Sequential multi-outcome (3a)** within-country, then the multivariate-FH
   comparison; treat the transport version as exploratory pending Alan's review.
5. **Anaemia validation (to-do 4)** as the larger engineering track once the above
   are demonstrated on the micronutrient outcomes.

Items to put in front of Alan before they become manuscript claims: weighted
ensemble optimization for a non-decomposable area-level criterion (1), the
transport version of the sequential chain (3a), and loss-based aggregate inference
plus formal EVSI (3c, C and E).
