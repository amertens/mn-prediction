# Transport calibration ladder and the MIMI intake-vs-biomarker divergence

Two additions to the methods, both prompted by the WFP / MIMI partner materials
(Voukelatou et al. 2025, *Sci Rep*; Martini et al. 2022, *Nature Food*) reviewed
in July 2026. Neither replaces existing work; each fills a specific gap.

- **Part A — the calibration ladder.** A three-rung framing of how we fix the
  documented transport *level* bias, graded by how much the target country tells
  us. Rungs 1 and 3 already exist; this adds the missing middle rung (a
  *gradient-anchored* calibration for a country with no national prevalence) and
  ties all three into one section.
- **Part B — the MIMI divergence track.** Use WFP's modelled dietary-inadequacy
  risk both as a candidate *predictor* and as the subject of a *scientific
  comparison* (where measured biomarker deficiency and modelled intake
  inadequacy disagree), the direct analogue of our GBD zinc diet-supply-vs-
  biomarker story.

Honesty flags are carried through in the style of
`METHODS_TODO_IMPLEMENTATION_PLAN.md`: what is standard versus what is
novel-in-this-setting and must go to Alan before it becomes a manuscript claim.

---

## Part A. The transport calibration ladder

### The problem it addresses

Under leave-one-country-out (LOCO) transport, predictions keep the district
*ranking* but are pulled toward the training-country mean — the documented
*level* bias (see `exploratory_loco_findings.md` and the `fe_transport_level_offset`
finding). Everything below is a way to restore the level without disturbing the
ranking, differing only in **what the target country lets us anchor to**.

| Rung | What the target country gives you | Method | Status |
|---|---|---|---|
| **1** | A known national prevalence (past survey, VMNIS/BRINDA, DHS anaemia) | One-parameter logit shift to that mean | **Done** — `simplified subset/methods/national_anchor_calibration.R` (to-do 2) |
| **2** | *No* national number, but a usable **structural gradient** (wealth / urbanicity / development) | One-parameter logit shift to a **borrowed reference-stratum** level | **New** — `sensitivity/13_gradient_anchored_calibration.R` (this note) |
| **3** | Nothing | **Partial-identification bounds** from the across-country shift range | **Done** — roadmap 3c option D, `simplified subset/methods/aggregate_inference.R` |

All three are single-parameter and strictly rank-preserving, so they never touch
the part of the prediction that *does* transport (the ranking); they only move
the level, which is the part that does not.

### Rung 2: gradient-anchored calibration (the new piece)

**Origin.** WFP's MIMI-ML calibrates a transported model *with no ground truth*
by tuning one parameter (a classification threshold) so the output reproduces an
expected structural gradient — deficiency falls as wealth or urbanicity rises —
verified with "sanity checks" on wealth-quintile and urban/rural contrasts.

**Adaptation to our data.** The simplified subset has no wealth or urban/rural
field; its proxies are environmental/malaria. We therefore use the available
**development / urbanicity axis** (Global Human Modification index, built-up
surface, cropland population density) as the structural gradient, and we anchor
the *level* to a **reference stratum** rather than to the whole distribution:

- Reference stratum = the **most-developed tertile** of areas. Its deficiency
  level is the most cross-country-comparable part of the distribution
  (development converges at the top; the national mean is dominated by each
  country's idiosyncratic rural burden). We **borrow** the top-tertile level from
  the *surveyed (training)* countries and shift the transported target
  predictions to match it.
- The shift is one logit delta → rank-preserving, exactly like rung 1. Only the
  anchor changes: a borrowed subgroup level instead of a known national mean.

**Why it is genuinely deployable.** Rung 1's anchor, in the demo, uses the target
country's *own* true mean (an oracle). Rung 2 borrows only from training
countries, so it needs **zero target-country biomarker data** — the actual use
case. The script reports all three (raw, gradient-anchored, national-anchored
oracle) so the question is quantitative: *how much of the oracle's level
correction does the deployable anchor recover?*

### Prototype results (admin-2, LOCO, elastic-net base learner)

Mean MAE (percentage points) across the four held-out countries:

| Outcome | Raw transport | Gradient-anchored (deployable) | National-anchored (oracle) | Gradient sign stable? |
|---|--:|--:|--:|--:|
| **women_vitA** | 3.73 | **2.26** | 2.53 | **Yes (4/4)** |
| **women_iron** | 20.79 | 15.61 | 9.62 | No (0/4) |
| **child_iron** | 13.56 | 17.52 (worse) | 10.95 | No (mostly) |

The result is nutrient-specific, and that is the point:

- **Vitamin A** — the structural gradient is **stable in sign** across all four
  countries, and the deployable gradient anchor **recovers essentially all** of
  the oracle's level correction (it even edges it, 2.26 vs 2.53 pp), with no
  target-country data. This is the success case.
- **Iron** — the development→deficiency gradient is **weak and sign-unstable**
  (iron deficiency is confounded by malaria, inflammation, and haemoglobin-
  opathies, not cleanly by development), so the anchor is unreliable: it helps
  women_iron partly by luck but *hurts* child_iron. The built-in **sign check
  flags exactly these cases** (`sign_ok = FALSE`).

This maps onto the existing `fe_*` finding that predictor domains should be
matched to outcome biology (vitamin A transfers on environment; iron on malaria
and modelled health burden). The gradient anchor inherits that: it works for the
outcome whose burden really does track the structural axis.

### The two failure modes the script surfaces (use both diagnostics)

1. **Sign instability** — `grad_train` and `grad_target` have opposite signs
   (the gradient does not transport). Then the borrowed anchor pulls the level
   the wrong way. Flagged as `sign_ok = FALSE`.
2. **Level non-transport** — the sign is stable but the borrowed top-tertile
   level (`borrow_top`) badly misses the target's true one (`true_top`), e.g.
   child_iron/Gambia (0.15 borrowed vs 0.35 actual). The sign check passes but
   the anchor still overshoots.

Report both columns; a country needs *both* to pass before the gradient anchor
is trustworthy.

### Honesty flags (rung 2)

- **Standard / low-novelty:** the logit shift and benchmarking of small-area
  estimates to an aggregate.
- **Novel-in-this-setting — NOT yet a manuscript claim (take to Alan):**
  anchoring the transported *level* to a borrowed structural-gradient stratum is
  a non-standard identification assumption. It is valid only when (a) the
  top-tertile level is transportable and (b) the gradient sign is stable — both
  testable, both printed. The honest headline is conditional: *"gradient-anchored
  calibration recovers most of the level error for outcomes with a stable
  structural gradient (vitamin A), and a built-in sign check identifies the
  outcomes/countries where it must not be used (iron, Malawi)."* Do not present it
  as a general fix.

### Does the outcome resolution (cluster / individual) change any of this?

Short answer: **no, and that invariance is itself the effective-*n* argument.**
Tested in `sensitivity/14_outcome_resolution_calibration.R`, which runs the same
LOCO transport at three resolutions — admin-2→admin-2, cluster→cluster, and
cluster→admin-2 (aggregate the cluster predictions back up) — and reports raw
metrics plus the oracle national-anchor correction. Means over the four held-out
countries:

| Outcome | Pipeline | Spearman (rank) | \|level bias\| pp | MAE after anchor (noise floor) |
|---|---|--:|--:|--:|
| women_iron | admin2→admin2 | −0.04 | 15.9 | 9.78 |
| | cluster→cluster | −0.09 | 11.7 | 11.00 |
| | cluster→admin2 | −0.16 | 13.9 | **9.77** |
| child_iron | admin2→admin2 | 0.16 | 10.0 | 11.00 |
| | cluster→cluster | 0.17 | 10.8 | 11.80 |
| | cluster→admin2 | 0.17 | 11.8 | **10.62** |
| women_vitA | admin2→admin2 | 0.04 | 4.0 | 2.52 |
| | cluster→cluster | −0.01 | 1.4 | 2.79 |
| | cluster→admin2 | 0.01 | 1.6 | **2.54** |

Four things hold across all outcomes:

1. **Ranking is a property of the area-level signal, not the outcome
   resolution.** Spearman is essentially the same at every resolution
   (child_iron ~0.16–0.17 throughout; iron weak/negative throughout). Scoring at
   cluster level adds noise but neither creates nor destroys rank signal.
2. **The transport level bias survives at every resolution** — it is an aggregate
   property, so you cannot escape it by going finer, and the calibration ladder
   applies unchanged. The oracle anchor removes it at every resolution.
3. **Cluster-level *evaluation* just adds a sampling-noise floor.** Cluster `n` is
   ~9–30 people, so a cluster "prevalence" is mostly binomial noise; the residual
   after removing the level is higher (women_iron 11.0 vs 9.8 pp). Scoring the
   *estimand* — an area/national aggregate — against noisy cluster prevalences
   penalises you for noise, not model error.
4. **Modelling at cluster level and aggregating back recovers the admin-2 result**
   (bold column: 9.77 ≈ 9.78; 2.54 ≈ 2.52; 10.62 vs 11.00). No free lunch from a
   finer outcome — consistent with effective *n* = number of areas.

**The one documented exception is the base learner, not the evaluation.** The
project's finding that an *individual-level* SuperLearner aggregated to admin-2
ranks districts best (Pearson ~0.53 vs the area ensemble ~0.12) is about *where
you fit*, not *where you score*: a proper class-weighted individual classifier
captures the within-area link and rare-outcome structure, then aggregates. That
is not contradicted here — the naive cluster-prevalence elastic-net in script 14
is not an individual-level classifier, so it does not show that gain. The clean
separation for the manuscript:

- **Estimand / evaluation resolution — always the aggregate (area or national).**
  Scoring at cluster/individual only injects sampling noise. Invariant to
  resolution; the calibration ladder is untouched.
- **Modelling / fitting resolution — can matter for the base learner.**
  Individual-level fitting may improve the aggregated prediction (documented);
  cluster-prevalence fitting does not. Either way it feeds the same
  aggregate-level evaluation and the same calibration.

*Caveats:* run on the simplified-subset elastic-net, not the corrected SL; true
individual-level rows live in `_targets_full/` and were not re-run here, but the
cluster level (small `n`) is a faithful proxy for the noise-injection effect and
the project's own individual-level SL result covers the fitting-resolution
question.

### How rung 2 slots into the manuscript

Present the ladder as one figure: raw → rung 1 (if an anchor exists) → rung 2 (if
a gradient exists) → rung 3 (bounds otherwise), with the per-nutrient recovery
table above. It converts the bare "transport fails on level" finding into a
graded, operational answer: *here is exactly how much we can repair, given
exactly what a data-poor country can supply.*

---

## Part B. The MIMI intake-vs-biomarker divergence track

### The idea in one line

WFP's MIMI produces a modelled **risk of inadequate dietary micronutrient
intake** (from HCES food-consumption data vs Harmonized Average Requirements)
for the same five nutrients we model as **biomarker deficiency**. Intake
inadequacy is *biologically upstream* of biomarker deficiency, so it is both a
candidate predictor and a scientifically interesting comparator.

### B1. MIMI intake risk as a proximal predictor

**Motivation.** Our most transportable feature is soil micronutrient content —
because it sits on the causal pathway (soil → crop → intake → biomarker). MIMI
intake inadequacy sits one step *further down that same pathway* (intake →
biomarker) and should, in principle, carry more direct signal than any of our
environmental proxies.

**Plan.**
1. Obtain, under the WFP collaboration/DUA, MIMI's modelled intake-inadequacy
   risk per nutrient at the finest shared administrative unit (admin-1, possibly
   admin-2) for our four countries.
2. Add it as one candidate covariate domain (`mimi_intake_<nutrient>`), linked by
   the standard `(country, admin)` key used everywhere in the pipeline.
3. Test incremental value under honest evaluation: within-country (cluster/
   spatial-block CV) and LOCO, as (a) an add-on to the parsimonious spatial-
   spline+soil model and (b) a standalone predictor, per nutrient. Expect the
   clearest gains for zinc and vitamin A, where intake is the dominant pathway.

**Honesty flag.** Straightforward feature engineering and evaluation (low
novelty). One caveat to state: MIMI risk is itself a *modelled* quantity, so
using it as a predictor imports its modelling assumptions; report results with
and without it, and treat it as a proxy, not ground truth.

### B2. The intake-vs-biomarker divergence analysis (the scientific payoff)

This is the analogue, with a *measured biomarker on both sides*, of the GBD
divergences catalogued in `PROJECT_STATUS_2026-06.md` §7 — most directly the zinc
case, where GBD uses a diet-supply proxy and we have serum zinc.

**Estimand.** Per nutrient and area, the signed discrepancy
`D_a = biomarker_deficiency_a − intake_inadequacy_a` (both on the prevalence
scale, harmonised to a common population group and reference period). Positive
`D_a` = more biological deficiency than intake alone predicts (points to
absorption/inflammation/infection losses — e.g. malaria and haemoglobinopathies
for iron); negative `D_a` = intake looks inadequate but biomarkers are
adequate (points to adaptation, fortification, or supplementation coverage).

**Plan.**
1. Harmonise the two prevalence surfaces (population group, thresholds, admin
   level, survey-year alignment) — reuse the existing linkage machinery.
2. Map and model `D_a`: where do the two agree, and where do they diverge
   systematically? Regress `D_a` on the confound/pathway proxies we already have
   (malaria burden, inflammation markers where available, fortification/
   supplementation coverage, soil) to *explain* the divergence rather than just
   describe it.
3. Tie back to the biomarker literature: iron's positive divergence should track
   malaria/inflammation (the same drivers that make ferritin hard to interpret);
   vitamin A's should track supplementation coverage (the driver GBD's modelled
   decline leans on).

**Why it matters.** It reframes "intake vs biomarker" from a *limitation* of
WFP's method (their pipeline is intake-based) into the **most valuable joint
product** of the collaboration: a measured-biomarker check on modelled intake
inadequacy, and a map of *where diet-based targeting would misallocate* because
the biological deficiency is driven by absorption/loss, not intake. It also
strengthens the GBD-comparison chapter with an independent, non-GBD intake
benchmark.

**Honesty flags.**
- Descriptive mapping and regression of `D_a` on observed proxies is standard.
- **Needs care (flag to Alan/nutrition leads):** any *causal* reading of the
  divergence (e.g. "malaria explains the iron divergence") is an ecological,
  correlational claim among confounded area-level variables — state it as
  hypothesis-generating. Comparability of the two prevalence definitions
  (biomarker threshold vs H-AR intake cut-point, different population groups) is
  the main threat to validity and must be audited before the discrepancy is
  interpreted.

### Dependency

B1 and B2 both require MIMI outputs for our four countries via the WFP
collaboration. Scope this as a **joint deliverable** with the WFP team rather
than something the core team builds alone — it is the concrete reason to have
them in the alliance, and it is the analysis that most directly needs both
sides' data.

---

## Suggested sequence

1. **Rung 2 gradient anchor** — prototype done (`sensitivity/13_*`); take the
   identification assumption to Alan, then fold the ladder figure into the
   transport section. Quick.
2. **B1 (MIMI as predictor)** — as soon as MIMI outputs are available under the
   DUA; low-novelty feature test.
3. **B2 (divergence analysis)** — the larger, higher-value scientific track;
   co-owned with WFP; slots alongside the GBD-comparison chapter.

## Pointers

- Rung 1: `simplified subset/methods/national_anchor_calibration.R`
- Rung 2: `sensitivity/13_gradient_anchored_calibration.R`
  (writes `results/sensitivity/gradient_anchored_calibration.csv`)
- Rung 3 / bounds: `simplified subset/methods/aggregate_inference.R`
- Roadmap context: `docs/METHODS_TODO_IMPLEMENTATION_PLAN.md` (to-do 2, 3c)
- GBD divergence framing: `docs/PROJECT_STATUS_2026-06.md` §7
