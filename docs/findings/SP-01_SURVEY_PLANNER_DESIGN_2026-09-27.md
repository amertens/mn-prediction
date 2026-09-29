# SP-01. Does model-guided district selection preserve the map when a survey can only reach k districts?

27 September 2026. Design written before the script ran (the repo's rule after
the September audit: objectives and metrics first, numbers second). Script:
`scripts/protocol_v2/64_survey_planner_validation.R`. Results:
`results/tables/protocol_v2/survey_planner_validation.csv` (+ `_summary`).

## The question

The dashboard's Plan-a-survey tab shows what a *small* survey buys (AR-01:
anchor + ranking against district/regional surveys of the same size). It does
not answer the question a planner actually asks next: **if the next survey can
afford biomarkers in only k districts, which k should they be** — and does
choosing them with the transported model beat choosing them the way surveys
are usually spread?

This is the design question the adaptive-geostatistical-design literature
answers with sequential surveys (Chipeta et al. 2016, Spatial Statistics;
Kabaghe et al. 2017, PLOS ONE; Andrade-Pacheco et al. 2020, Sci Rep — batch
selection against a threshold). We cannot run a sequential field pilot, but the
four surveys support the retrospective version: hide most of a country's
surveyed districts, "survey" only k of them, fit the in-country index on those
k, and score the resulting map on the districts held back.

## Design

Per cell (country × outcome, level target, cells with ≥ 14 surveyed
districts):

- Universe = the cell's surveyed districts (only they carry truth).
- The **selection score** available before any blood is drawn: the transported
  climate + soil index (the pre-registered candidate), fitted on the other
  three countries at district level exactly as in
  `scripts/policy_deck/10_viz_tables.R` block C, applied to the target
  country. Selection never sees the target country's outcomes.
- Fractions f ∈ {0.2, 0.35, 0.5, 0.75}; k = ceiling(f · n_surveyed), k ≥ 5.
- Selection arms, 40 replicates each:
  1. `random` — simple random sample of k (the reference).
  2. `pps` — k districts drawn with probability proportional to population
     (how sample tends to be spread when EAs are drawn PPS).
  3. `spread_model` — districts sorted by the transported score, cut into k
     equal strata, one drawn per stratum (model-stratified spread; the
     covariate-balance idea of GRTS / model-assisted design).
  4. `extremes_model` — k/2 from the top and k/2 from the bottom of the
     transported score (anchor the gradient at its ends; the optimal-design
     analogue for a linear index).
- Fit: the protocol's zero-tuning index (`arm_domain_index_v2`) on the k
  selected districts; domain-PC basis built once per cell from the country's
  full predictor matrix (the basis is unsupervised — X only — and every
  district's X is public, so a real planner has it).
- Score on the held-out surveyed districts (n − k):
  - `spearman` — ranking accuracy (primary);
  - `capture` — share of the held-out worst fifth (by survey) that the fitted
    map places in its predicted worst fifth (top-k concordance, k = fifth);
  - `mae_pp` — anchored prevalence error: predictions calibrated (IS-01 map)
    and anchored so the k sampled districts' weighted mean matches their
    survey mean, in percentage points (secondary).
- Comparator row `transport_only`: the transported score itself, no in-country
  fit, scored on the same held-out districts (what k = 0 gives).

## Pre-stated readings

- Primary: at f = 0.35 and 0.5, `spread_model` beats `random` on held-out
  Spearman in the majority of cells, and the paired mean delta is positive.
- Secondary: `extremes_model` behaves like `spread_model` or better at small
  f; `pps` sits at or below `random` (population is not an information
  design).
- Whatever the sign, the measured deltas go on the dashboard's Plan-a-survey
  tab as "what choosing districts well buys", with random as the stated
  reference and this note linked. If `spread_model` does not beat random, the
  tab says so and the planner is presented as a prioritisation aid only
  (which districts to *confirm*), not a precision gain.

## Amendment (same day, after the first run, before any dashboard use)

The first run scored only the district map. These surveys are commissioned for
the **national prevalence**, so the design must also say what k-district
selection does to it. Added metrics, per cell × arm × fraction:

- `nat_bias_pp` — mean over draws of (national estimate from the k visited
  districts − the same estimator over all surveyed districts), in percentage
  points. The estimator is the population-weighted mean of district
  prevalences; for `spread_model` a second, stratified estimator weights each
  model-score stratum by its population (spread selection is one-per-stratum
  stratified sampling — a probability design — so this is its design
  estimator).
- `nat_rmse_pp` and `nat_ci_pp` (1.96 × the across-draw SD) — the
  between-district-selection component of the national estimate's error and
  precision. Within-district sampling noise is the observed surveys' own and
  is held fixed; the reported precision is therefore the component that the
  choice and number of districts moves.

Pre-stated reading for the amendment: `random` and the stratified
`spread_model` estimator are unbiased for the all-district figure;
`extremes_model` (and any burden-first targeting) is an informative sample
and biases the national number, which is why the dashboard presents
model-guided targeting as a top-up to a probability design, never as the
sampling frame.

## Results (recorded 27 September, after the amendment run)

Mean over cells, level target (`survey_planner_validation_summary.csv`):

| Share visited | Random | Population (PPS) | Spread across model | Model extremes |
|---:|---:|---:|---:|---:|
| 0.20 | 0.283 | 0.276 | 0.287 | 0.306 |
| 0.35 | 0.324 | 0.310 | 0.329 | 0.337 |
| 0.50 | 0.357 | 0.312 | 0.362 | 0.343 |
| 0.75 | 0.456 | 0.394 | 0.497 | 0.289 |

- Primary reading met on its letter, not in substance: spread beats random
  in 16/22 cells at 0.35 and 14/22 at 0.5, but the mean gain is +0.005,
  a tie under the project's 0.03 rule. At 0.75 the gain is +0.04 (11/16).
- PPS (population-first) is worse than random at every share.
- Extremes help only at small shares and collapse at 0.75 (range restriction).
- **Transport crossover (not pre-stated):** the transported climate + soil
  score alone ranks the held-out districts at 0.343. A model fitted on the
  visited districts only beats it from about half the districts upward
  (0.357 at 0.5); at 0.20 and 0.35 it is worse (0.283, 0.324). With few
  surveyed districts, the other countries' model is the better ranking;
  combining the two is the obvious next arm.
- National estimate (amendment): random 0.33 pp mean |bias|, ±3.96 pp
  (95% half-width) at half the districts; spread with Hajek weights 0.30 pp,
  ±3.44; extremes 1.37 pp; PPS 0.75 pp, ±1.83.

## What this is not

Not a cluster-level design (selection is at the district grain the pipeline
predicts at); not a sequential/adaptive pilot (one batch, no refit between
rounds); not a sample-size calculation (AR-01 remains the size story). It is
the missing middle: *which* districts, at a fixed k.
