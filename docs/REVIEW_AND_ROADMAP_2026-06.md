# Critical review + methodological roadmap (June 2026)

Companion to [`CRITICAL_REVIEW.md`](CRITICAL_REVIEW.md) (the P1–P8 code/methods
audit) and the manuscript ([`manuscript_mcn.qmd`](manuscript_mcn.qmd)). This
document (a) re-reviews the codebase, analysis approach, and results with fresh
eyes, and (b) turns the Mertens–Hubbard–Kim meeting (transcript, June 2026) into
a prioritized methodological roadmap, with an honest read on what each idea can
and cannot fix given the data.

---

## 0. The one thing to keep in front of every decision

Andrew's own framing in the meeting is the correct prior: **the binding
constraint is the data, not the model.** Four surveys, ~206 admin-2 polygons,
biomarkers whose cross-country comparability is itself disputed, and proxy
predictors that are individually weak. The honest out-of-fold results already in
the manuscript bear this out: within-country admin-2 ranking is modest (median
Pearson r ≈ 0.53 for the individual-level aggregate, ≈ 0.12 for the directly-fit
area ensemble), LOCO transport is near chance on most outcome–country splits
(median r ≈ 0.21), and a plain OLS / Fay–Herriot beats the ensemble on ~10 of 13
head-to-head LOCO splits.

**Implication for this roadmap:** almost none of the meeting's ideas will move
within-country AUC much. Their value is in (1) *reframing the deliverable* to
something the data can actually support (ranking / top-N targeting, value-of-
information sampling), (2) *reducing variance* (multi-outcome borrowing,
PCA learner, decision-focused ensembling), and (3) *honest inference* so the
product is "credible, appropriately-uncertain estimates" rather than false
precision. That reframing is itself the most fundable result. Prioritize it.

---

## 1. Meeting summary (for the record)

Three headline directions, in Alan's own end-of-call summary:

1. **Borrow strength across outcomes.** Predict the micronutrients *sequentially*
   in an arbitrary order (most-informative first: iron → vitamin A → B12 →
   folate → zinc), feeding the predicted upstream outcome in as a covariate for
   the next — the classic forward G-computation factorization of a joint density.
   Handle the (frequent) missingness of an upstream outcome with a
   missing-indicator basis function plus imputation, exactly as inside the SL.
2. **Optimize the ensemble for an aggregate, decision-relevant criterion.**
   Base learners minimize MSE/NLL, but the *ensembling / discrete-SL selection
   step* can target any loss — e.g. correctly classifying which admin units are
   above a public-health threshold, or capturing the top-N highest-risk regions
   with controlled false-negatives. This is structurally identical to Andy's
   sensitivity-constrained / quintile work.
3. **Inference on the predictions.** Get valid uncertainty on the *aggregate*
   product (national or admin prevalence), via van der Laan's loss-function-
   based CV inference, conformal, or bootstrap — and use it to drive
   value-of-information: "where are we both high-risk and most uncertain → sample
   there next."

Supporting ideas raised:
- **Calibrate to known national prevalences** (constrain the predicted mean to a
  known/historical national value). Pure level-calibration; does not change
  ranking but directly attacks the documented LOCO level-offset failure.
- **Binomial family** for the 0–1 area-prevalence outcome (proper loss) instead
  of a Gaussian link.
- **PCA / principal-components learner** given many predictors across clear
  domains (rainfall, poverty, soil, …).
- **Within- vs between-country error decomposition** from a single pooled CV.
- **Anemia as a validation outcome** — DHS-measured everywhere, large n, trusted
  — to validate the *methodology* even if the micronutrient biomarkers are noisy.
- Andy's parsimonious-DHS manuscript (boil DHS down to ~4 highly-explanatory
  proxies) is the same "cheap proxy indicator" goal as this project's Aim 1.

---

## 2. Codebase review (delta since CRITICAL_REVIEW.md)

The structure is sound: a `targets` production pipeline (`_targets.R` →
`_targets_full/`), a parallel corrected pipeline (`_targets_corrected.R`,
P1–P8), a modular Shiny dashboard, and now a clean simplified subset. The
recent cleanup archived 15 throwaway scripts and both pipelines still validate
byte-identically. What still deserves attention:

- **Two parallel pipelines is a maintenance liability.** The corrected pipeline
  (leakage-free in-fold preprocessing, spatial/cluster-block CV, out-of-fold
  calibration, partial-pooling SAE) implements what the manuscript already
  treats as the analysis of record. Plan a path to **fold the corrected
  methods back into the main pipeline** and retire the duplicate, so there is a
  single source of truth before Andy and external ML teams build on it.
- **`src/` legacy tree** (~60 country/DHS/GEE scripts with clear `_V2`
  duplicates) is upstream data-prep provenance that no longer matches the
  `targets` data flow. It needs a dedicated, reviewed archival pass (flagged in
  the cleanup report) — important hygiene before sharing the repo with external
  modelers.
- **The feature-engineering / linkage layer is the highest-leverage code to
  invest in**, per Andrew: time-ordered, space-aware linkage of proxies to
  clusters. This is where real predictive signal would come from, more than any
  learner change. Treat "better linkage + feature engineering" as a first-class
  workstream, not a chore.
- **Design weights / effective n.** The pipeline (and the simplified subset)
  assume a fixed DHS design effect of 1.5. With ~1 cluster per admin-2 in places
  this is a real approximation; worth replacing with survey-package design
  variances where the cluster structure supports it (see
  [`survey_weighting_critique.md`](survey_weighting_critique.md)).

## 3. Analysis-approach review

What is right and should be protected:
- Honest out-of-fold evaluation, cluster-blocked CV, LOCO for transport, and
  benchmarking the ensemble against SAE/OLS. This rigor is the manuscript's main
  credibility asset — the finding "ensembles don't beat SAE here" is only
  believable *because* of it.
- Proxy-only predictor discipline (excluding `gw_` biomarker items) so the model
  is deployable where no biomarker survey exists.

Where the approach is fighting the data (and the meeting points to better
framing):
- **Predicting absolute prevalence level is the wrong target for transport.**
  The manuscript already shows LOCO failure is a *level-calibration* error
  driven by between-country prevalence heterogeneity, not (mostly) a ranking
  error. Continuing to score transport by absolute MAE will keep producing
  pessimistic headlines. **Re-center the transport deliverable on ranking /
  top-N classification + an explicit national-level anchor.** (Meeting ideas 2
  and the calibration point.)
- **Per-outcome independent fitting leaves covariance on the table** (meeting
  idea 1), though see §5 for why I'd temper expectations.
- **Within-country individual models have no within-area predictor variation**
  (Andrew's observation) — so the individual vs area distinction is partly
  illusory for the proxy-only set, and the area-level framing is the honest one.

## 4. Results review (sanity check of the headline numbers)

The numbers in the manuscript are internally consistent and I would not change
them. Three I would foreground differently:
- **National estimates are good** (within ±5 pp for all 24 country–outcomes,
  median |Δ| 0.4 pp). The product that already works is *national* prevalence
  for a country with a recent-ish survey. The Davis collaborators' view in the
  meeting — "accurate national prevalences would itself be a huge accomplishment"
  — is supported by the results and is the safest headline.
- **Transport ranking is salvageable where level is not.** Folate preserves
  ordering (slope 0.37) once cutoffs are harmonized; the Sierra Leone admin-3
  re-analysis moved Spearman from negative toward positive (B12 +0.30). This is
  evidence the *ranking* deliverable is more robust than the *level* one — which
  is exactly what the meeting's top-N reframing exploits.
- **Vitamin-A non-transport and B12 instability under rare-outcome class-
  weighting** are real and should stay as cautions; the decision-focused
  ensembling (§5.2) is partly a fix for the class-weighting instability.

---

## 5. Methodological roadmap (prioritized, with honest expectations)

Priority key: **P1** = do first / highest value-for-effort; **P3** = worthwhile
but lower or more speculative.

### P1 — Decision-focused evaluation & ensembling (meeting idea 2)
**What:** Reframe the primary deliverable as *correctly classifying high-risk
admin units* (above a WHO public-health threshold) and *capturing the top-N
highest-prevalence regions*, with an asymmetric loss that penalizes missing a
truly high-risk region (false negative) more than a false alarm. Implement it in
the **SL ensembling / discrete-SL selection step**, leaving base learners as-is
(they still do MSE/NLL). This is Andy's sensitivity-constrained / quintile
machinery applied to areas.
**Why P1:** matches the actual policy need (Andy/Alan), is statistically *easier*
than level estimation, directly uses the part of the signal that survives
(ranking), and is cheap (reuses cached base-learner fits with new ensemble loss).
**Honest expectation:** this is the single most likely thing to turn a
"disappointing" result into a fundable one. It reframes, it does not manufacture,
signal.

### P1 — National-anchor calibration (meeting calibration point)
**What:** When transporting to a country, constrain the predicted prevalence
distribution so its (population-weighted) mean equals a known/historical national
prevalence from VMNIS/BRINDA/DHS-anemia. A one-parameter recalibration on the
logit scale (or isotonic) that shifts level without touching rank.
**Why P1:** the manuscript localizes transport failure to *level*, and national
anchors are obtainable even where individual data is not. Low effort, directly
targets the documented failure mode.
**Honest expectation:** fixes the level bias by construction wherever an anchor
exists; does nothing for ranking (which is the point — pair it with P1 ranking).

### P2 — Proper loss for area prevalence (meeting binomial point)
**What:** Fit area-prevalence learners with `family = binomial` / quasi-binomial
(any 0–1 response), weighted by effective n, instead of a Gaussian link. Already
partly present (quasi-binomial GLM benchmark); make it the default link for the
SL on prevalence.
**Why P2:** correctness and mild gains; bounds predictions to [0,1], respects
mean-variance. Low risk.

### P2 — Sequential multi-outcome borrowing (meeting idea 1)
**What:** Order outcomes by information (iron → vit A → B12/folate → zinc); at
each step add the *predicted* upstream deficiencies as covariates, with
missing-indicator + imputation for outcomes absent in a country.
**Why P2 not P1:** clever and worth testing **within-country first**, where the
covariance is real and measured. But for *transport* you would condition on a
*predicted* upstream outcome that itself transports poorly, compounding error —
so I expect within-country gains and limited transport gains. The simplified
subset now carries all outcomes (wide format), so this is directly testable there.
**Honest expectation:** modest within-country variance reduction; do not expect
it to rescue transport.

### P2 — Inference on aggregate predictions + value-of-information (meeting idea 3)
**What:** Put valid uncertainty on the *aggregate* product (national / admin
prevalence): van der Laan loss-based CV inference for the aggregate target,
conformal for area predictions, bootstrap for LOCO r (already in the manuscript).
Then compute a VOI map: rank candidate survey sites by `(risk × uncertainty)` to
answer "where should Gates fund the next biomarker survey."
**Why P2 (but strategically P1):** the VOI deliverable ("here's where to sample
next for maximum information") is politically the strongest pitch to the funder
and is robust to the data being weak. Implementation is moderate effort.

### P3 — PCA / principal-components learner (meeting PCA point)
**What:** Add a PCA-regression learner (or domain-block PCA) to the SL library.
**Why P3:** cheap to add, theoretically apt for many correlated predictors across
domains; unlikely to be a game-changer but a reasonable library member. Low risk.

### P3 — Within/between error decomposition from pooled CV (meeting decomposition)
**What:** A single pooled (non-spatial) CV whose validation folds mix
same-country and other-country rows, reporting within- vs across-country error
separately. **Caveat (Alan flagged it mid-call):** the training fold is always a
country mixture, so this is not a clean estimand — it is a descriptive
diagnostic, not the transport estimator. Keep LOCO as the transport estimate.
**Why P3:** useful diagnostic framing, but don't oversell it.

### P3 — Anemia validation track (meeting anemia point)
**What:** Run the whole pipeline on **anemia** (haemoglobin, DHS-measured in
nearly every SSA country, large n) to validate the methodology on a trusted,
abundant outcome and to expand the LOCO panel far beyond four countries.
**Why P3 by effort, P1 by strategic value:** large lift (new ingestion), but it
is the cleanest way to separate "the method works, the biomarkers are noisy" from
"the method doesn't work." Strong candidate for the scale-up phase.

### Cross-cutting — single source of truth
Fold the corrected pipeline into the main pipeline before the external ML teams
and Andy build on it (see §2).

---

## 6. Suggested sequence against real deadlines
- **Before the September Ghana talk:** P1 decision-focused evaluation + P1
  national-anchor calibration + the VOI map. These are the high-value,
  low-compute reframes that present well to a non-technical nutrition audience
  ("we can reliably flag the highest-risk districts and tell you where to
  sample next") and don't require new data.
- **Fall (with Andy):** P2 sequential multi-outcome (within-country), proper
  binomial loss, aggregate inference; start the anemia validation scoping.
- **Scale-up phase:** anemia track + East/Southern-Africa expansion (the
  external-validity step the manuscript already flags), then revisit transport.

See [`ANDY_KIM_PROJECT_PLAN.md`](ANDY_KIM_PROJECT_PLAN.md) for the concrete
starter project on the simplified subset.
