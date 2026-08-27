# Project status update — August 2026

This is an addendum to `PROJECT_STATUS_2026-06.md`, covering what changed in
the two months since: new predictor domains, a validated parsimonious LOCO
predictor set, area-level survey-design fixes, several negative results
worth recording, and a critical-review audit with its follow-up
verification. Read the June document first for full context (data sources,
harmonization conventions, evaluation design); this note assumes it.

## 1. New predictor domains unlocked

Three cross-country linkage bugs were fixed, each of which had been silently
zeroing out a whole predictor domain for one or more countries rather than
producing a visible error:

- **WFP food prices**: case-sensitive admin-name matching failed for Ghana
  (WFP uses ALL-CAPS names; GADM uses title case) and an admin-level
  mismatch failed for Malawi (WFP's "admin2" field is actually GADM
  Admin1-equivalent). Fixed via fuzzy, case-insensitive matching with
  admin-level auto-detection (`R/external_data.R`). Replaced raw price
  levels (not comparable across countries/markets) with a scale-free
  price-deviation/spike index.
- **DHS Admin2 covariates** (wealth, health system, sanitation, vaccination):
  each country already had these, but they never pooled across countries
  because the pipeline's own year-embedded column-naming convention (not the
  underlying indicators, which are DHS-standard and comparable by design)
  broke poolability. Fixed by stripping the year prefix to a canonical name
  (`R/dhs_harmonization.R`); unlocks 105 common columns across all 4
  countries.
- **Nighttime lights**: previously extracted for Ghana only, due to three
  independent bugs (a token gate that misfired for token-free extraction
  methods, an undefined variable, and country-matching that only
  coincidentally worked for Ghana). Fixed for all 4 countries.

## 2. Validated parsimonious LOCO predictor set

A stepwise, LOCO-native search (screen candidates by cross-country bivariate
strength, dedupe near-duplicate constructs, greedy forward selection scored
by LOCO AUC) landed on a 5-variable set: soil calcium variability, November
rainfall, malaria incidence rate, a vegetation-domain PCA component, and MAP
parasite rate. This set was tested against every predictor domain above,
both forced into the candidate pool and left to compete freely — **it won
both times**. The best-ranked newly-unlocked candidate (a DHS indicator)
placed only 51st of ~590 candidates by raw cross-country predictive
strength. This is a real negative result: for this outcome/country set, more
data domains is not the bottleneck.

Overall mean LOCO AUC with this set: **0.5718** (16 folds = 4 outcomes × 4
held-out countries), aggregate 95% CI [0.544, 0.600], excludes chance
(p ≈ 6×10⁻⁵). Only 6/16 individual folds are independently significant;
rare-outcome folds (e.g. women's vitamin A) are underpowered on their own.

## 3. Methods tried that did not beat this ceiling

Recorded because negative results here are as informative as positive ones:

- **XGBoost (shallow and deep), random forest (ranger)** — no LOCO advantage
  over a well-selected linear model on identical covariates.
- **Full SuperLearner ensembles vs. a single glm on identical, unscreened
  covariates** — no advantage across countries (SL's advantage shows up
  within-country, not across).
- **Wider candidate pools (100–150 variables)** — degraded LOCO
  transportability relative to the 5-variable set; more predictors added
  in-pool overfitting, not transportable signal.
- **Rank-based / learning-to-rank objectives** (XGBoost `rank:pairwise` by
  country) for predicting district-level prevalence — no method
  significantly beat regression or SL loss; underpowered (16 folds, 12–83
  areas/fold) to confirm even the numerically-best result.
- **Domain-generalization methods for few-environment transfer** (IRM,
  Group DRO) — neither improved on naive pooled regression with only 3
  training countries per fold; Group DRO significantly *underperformed*
  (−0.030 AUC, p=0.02) by overfitting to whichever training country
  currently had the highest loss. These methods need more domains than a
  4-country study supplies.
- **Seasonality/interview-timing reconstruction** — investigated as a
  possible confound and closed out: each country's survey window is short
  (61–92 days) and shows no significant within-country association with any
  outcome; timing is nearly collinear with Admin2 (fieldwork moves
  district-by-district), so admin fixed effects already absorb it. The
  parsimonious set's fixed-calendar-month variables (e.g. "November
  rainfall") turned out to be a data-driven search artifact rather than a
  deliberate design choice, but rebuilding them as interview-date-relative
  would not plausibly move the results.

## 4. Area-level SAE fixes

- **`compute_svy_admin2()` survey-design ordering**: the function used to
  filter to the outcome's domain (non-missing biomarker, non-missing
  Admin2) *before* declaring the survey design, rather than declaring the
  design on the full weighted sample and restricting to the domain
  afterward (the textbook-correct order for domain estimation). Filtering
  first can make a stratum that genuinely has 2 PSUs look like it has only
  1 whenever one PSU has zero valid reads for a given outcome, misstating
  its variance contribution. Fixed; verified against the pre-fix version
  across all 4 countries × 4 outcomes — identical point estimates and SEs
  everywhere in this project's actual data (a real fix, but it did not
  change any previously reported number).
- **Fay-Herriot degenerate-SE detection**: most Admin-2 areas hold a single
  PSU, so `survey::svymean` cannot form a between-cluster variance and
  returns SE ≈ 1e-17 — finite and positive, so it slipped past a `D <= 0`
  guard. The EBLUP shrinkage weight went to ~1 for exactly the districts
  that most needed shrinking (mean gamma 0.71 as coded vs. 0.26 with a real
  sampling variance; 29% of districts had gamma > 0.95). Fixed by also
  flagging any variance implausibly small relative to the simple-random
  -sample variance.
- **Reliability ceiling** (`R/area_reliability.R`, new): a median Admin-2
  district contributes 6–36 biomarker measurements, so much of the
  between-district spread in survey prevalence is sampling noise, not
  geography. Added a variance-decomposition ceiling on achievable Pearson r
  so a model's correlation can be read against what's actually attainable.
- **National-benchmark correction** (`R/benchmark_area.R`, new): these
  surveys are weighted to deliver a national estimate; Admin-2 is an
  unplanned domain. A single logit-scale shift constrains a model's
  population-weighted Admin-2 aggregate to reproduce the survey's own
  design-based national figure — monotone, so district ranking is preserved
  exactly.
- See `sandbox_parsimony/FINDINGS.md` for the full critical-review sandbox
  these fixes came out of, including the geography-left-out-of-the-model and
  logit-scale findings.

## 5. Survey-subsample "cost of accuracy" result

Tested the idea (raised in the August 2026 Gates strategy discussion) of
whether a country could cut survey costs and let covariate-based modeling
recover the lost subnational accuracy. Simulated reduced survey budgets
(50–90% of clusters retained, 25 replicates per fraction per country/outcome)
and compared a covariate-based area model against "use the local direct
estimate where you have one, else the national average," both benchmarked
against each country's own full survey.

**An initial pass had in-sample leakage** — the area model was trained on a
district's own subsampled estimate and then scored on that same district,
making it look artificially good. Corrected with genuine k-fold
cross-validation (every district's model prediction is held out from its own
training). Result, honestly scored:

- For any district that retains even a few clusters, the direct estimate
  beats the model by a wide, highly significant margin (+7.7pp RMSE,
  p < 10⁻³⁰⁰) — **don't use the model there**.
- The model's only real, statistically reliable value is for districts with
  **zero** surveyed clusters, where it beats guessing the national average
  by about 0.7 percentage points (≈6% relative reduction on an ~12pp error,
  p = 1.2×10⁻²⁰), consistently across retention levels.

This is a real but narrow result. The defensible framing is "a smarter
guess for the specific districts you deliberately decide not to sample," not
"cut your survey budget and the model recovers the loss" — the broader claim
does not survive honest cross-validation.

## 6. Critical-review audit and weighted-SL follow-up

An external-style audit of survey weighting, proxy-data linkage, external-
data temporality, and the SuperLearner ensemble setup surfaced:

- **The individual-level SuperLearner was fit completely unweighted** — no
  survey design at all (`mlr3_SL_clustered()` never set a weight role, and
  `mlr3superlearner()`'s CV-risk computation has no weight parameter to
  extend). Re-verified with a genuinely weighted refit (monkeypatched
  `make_mlr3_task`/`compute_loss`, unweighted reproduction checked to match
  production exactly first): **weighted mean LOCO AUC = 0.558 vs.
  unweighted 0.572 (p = 0.50, not significant)**. The common outcomes
  (child_vitA, child_iron) were essentially unchanged; the rare outcomes
  became substantially noisier under weighting (one fold dropped to 0.40),
  most likely from survey weights stacking with the pipeline's existing
  rare-outcome class-weighting. **The 0.57 ceiling is not an artifact of
  missing weights**, but production weighting needs the two schemes
  reconciled before it's turned on by default — `.sl_weight` is now carried
  through `build_pooled_dataset()` (within-country normalized to mean 1) for
  whoever picks this up.
  - Caught in the process: Malawi's weight column was on the raw DHS ×10⁶
    scale while the other 3 countries' were already normalized — invisible
    before because nothing pooled weights across countries. Fixed at the
    source.
- **`compute_svy_admin2()` filter-before-design ordering** — see §4, fixed.
- **Not yet acted on**: a 2km GPS buffer for the cluster-level comparator
  likely undercorrects for DHS's confidentiality displacement (up to
  5-10km rural); the food-security domain accepts Admin2 name-matches at
  just a 30% match rate. Neither confirmed to change reported numbers.

## 7. Proposed direction for the next year

Per the August 2026 Gates strategy discussion, reframe the flagship
deliverable from "predict deficiency in a country with zero data" (not
supported by current evidence) to "use a country's own survey plus modeling
to make a smarter choice about which districts to skip sampling" — a
survey-design decision aid, not a general accuracy multiplier. Concretely:

- Run 1-2 national case studies/workshops built around the corrected
  subsample result (§5), with the narrow, defensible framing above.
- Bring in outside costing expertise to convert the ~0.7pp zero-coverage
  edge into an honest $-value framing.
- Scope any GBD contribution as within-country-year disaggregation of an
  existing national estimate (constrained to reproduce the national total),
  not free-standing subnational/zero-data prediction — and note explicitly
  that remote-sensing predictors can't reliably extend back to 1990.
- Fix the two audit-confirmed issues most likely to be challenged by a
  skeptical reviewer (§6) before anything goes external.
