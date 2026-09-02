# Slide deck outline — rebuild, September 2026

For `docs/mn_prediction_slides.qmd`. Written after the 2026-09-01 audit,
protocol v2, and the data-source expansion. Companion to
`docs/findings/TWO_READINGS_2026-09d.md` (abstracts + paper outline) and
`docs/findings/PROTOCOL_V2.md` (numbers and the plain-language arm guide).

## What "protocol" means here

Throughout this outline and in `PROTOCOL_V2.md`, **protocol** means the set of
decisions about HOW A MODEL IS EVALUATED, as distinct from the model itself:

- how the data is split into training and test, and whether that split is
  replicated or reported from a single draw;
- what the model is compared against, and **what information that comparator is
  allowed to see** - a baseline that reads the held-out district's own survey is
  not comparable with a model that cannot;
- which quantity is predicted (deficiency prevalence, or the underlying
  biomarker concentration);
- how districts are weighted when scoring - equally, or by how much information
  each actually carries;
- which null distribution a significance claim is judged against;
- and which deployment question is being answered: filling gaps inside a
  surveyed country, reaching an unvisited region, or transporting to a country
  with no survey at all.

The analogy is a clinical trial, where the protocol - endpoint, comparator,
analysis plan - is registered before enrolment, and the drug is the
intervention. In prediction work the model is the intervention and the
evaluation design is the protocol, but nobody pre-registers it, so it usually
gets settled after the results are in, often without anyone noticing a choice
was made.

The project's central empirical finding is that **the protocol moved results
more than any model choice tested**. Fold draw alone moved a headline from 0.058
to 0.217. The choice of baseline moved a comparison from "covariates lose to a
survey mean" to "covariates win". No estimator in the tournament moved anything
by that much.

## The problem with the current deck

It is 74 slides in four accreted colour-coded eras — the original build, then
"Methods Review 2026-08" (red), "Full Rebuild 2026-08" (green), and "Covariate
Rebuild & Honest Baselines 2026-08-30" (teal). Each phase was **appended rather
than integrated**.

The consequence is not that some numbers are stale. It is that the deck is a
**changelog rather than an argument**: someone reading front-to-back meets
claims in their original confident form and only later, in a differently
coloured section, learns they were withdrawn. Several are withdrawn twice. A
presenter who stops early — which is what happens — presents retracted results.

Three specific slides now argue the opposite of the current evidence:

| Slide | Says | Now |
|:---|:---|:---|
| "Anchoring to regional totals is the design that pays" | anchoring 0.164 → 0.413 | **withdrawn** by the project's own jackknife (0.147 vs 0.156) |
| "Nothing survives FDR once space is accounted for" | 0 of 294 | the zero was **unattainable by construction** at 999 permutations |
| VMNIS "98% saturated, r_max 0.66" | near-ceiling | **two corrections behind** (r_max 0.869, ~75%) |

## The restructure

**From changelog to argument.** One narrative, ~30 slides, in which the
corrections are not a confession section but the **methodological contribution**
— which is also the honest read, since the protocol work is the most
transferable thing the project has produced.

**Three organising principles.**

1. **Lead with the current position.** The first ten slides should be defensible
   on their own if the presenter never reaches slide eleven.
2. **One estimand per claim.** Every performance number gets labelled with which
   of the three deployment questions it answers. Most historical confusion in
   this project came from mixing them.
3. **Corrections as method, not apology.** The framing is "here is what the
   field's DEFAULT evaluation choices cost, measured on our own data" - not
   "here is what we did carelessly." That distinction is both kinder and more
   accurate. A single unreplicated fold draw, a design effect assumed at 1.5, a
   baseline computed from all the data: these are the ordinary choices almost
   every paper in this literature makes, which is exactly why measuring their
   cost is worth a talk. This project caught several of them itself - the
   jackknife control that withdrew the 0.516 baseline was run here, before any
   external audit.

---

## Proposed structure

### Section 1 — The problem (3 slides)
1. **The biomarker data gap.** Keep the existing slide; it still works.
2. **What a programme actually decides.** Which districts to reach, with what.
   Frames every later metric as a targeting decision rather than an R².
3. **Why this is hard to evaluate honestly.** Four surveys, 14–87 districts per
   country, a noisy yardstick. Sets up the protocol section.

### Section 2 — Data (3 slides)
4. **What we hold.** Four national biomarker surveys, 24 country–outcome cells,
   451 harmonised predictors across 19 domains.
5. **The vocabulary, by domain and native resolution.** Include the admin-level
   column — it is where the IHME defect lived.
6. **What was added this year.** 373 → 451: IHME sub-national surfaces, WFP food
   prices, FAOSTAT supply, WorldPop demographic composition, Relative Wealth
   Index, agro-ecological and climate zones. A delivered asset, not a promise.

### Section 3 — The protocol (6 slides) — *the contribution*
7. **Three estimands, three baselines.** In-fill, region extrapolation, country
   transport, each against a comparator that saw the same information.
8. **What one fold draw costs.** Median range 0.174 across ten draws, maximum
   0.705. The published median r 0.058 came from a bottom-decile draw;
   replicated, it is 0.217.
9. **What an unfair baseline costs.** The old number-to-beat (0.516) read the
   held-out district's own survey. Jackknifed it scores 0.076 and loses.
10. **What the effective sample size really is.** Measured design effects
    1.0–5.6, median 2.4, not the assumed 1.5. 77% of Malawi districts have
    n_eff < 5.
11. **What dichotomising costs.** +0.10 to +0.11 for the continuous biomarker,
    consistently across all three estimands.
12. **The four conditions.** What any group should report. This is the slide
    people photograph.

### Section 4 — Results (7 slides)
13. **The leaderboard, in-fill.** Covariates 0.398 vs the jackknifed survey mean
    0.320 on the biomarker level; 0.286 vs 0.193 on prevalence.
14. **The leaderboard, region extrapolation.** Same ordering, wider margin over
    the raw-column fit.
15. **Transport.** Rank correlation 0.288, **positive in 21 of 22 cells**, up
    from 0.200 and 11 of 15.
16. **Targeting, in programme terms.** The worst-ranked fifth reaches 26% more
    burden than the survey's own regional averages — which are themselves no
    better than random.
17. **A named example.** Malawi child iron: targeted districts carry 43%
    deficiency against 26% nationally, capturing 49% of national burden.
    *Always name the country and biomarker; the pooled average is not a
    coherent quantity.*
18. **Simplest wins.** The zero-tuning domain index beats a penalised fit over
    all 451 columns by up to 0.38, and the raw fit goes negative under region
    extrapolation. Explain why: marginal weights need p estimates, a joint fit
    needs p²/2, and n is 14–87.
19. **What geography already does.** Spatial alone 0.391 vs 0.398 with
    covariates. Stated plainly — it is the honest frame for where covariates
    earn their place.

### Section 5 — Where covariates earn their place (3 slides)
20. **Two situations, precisely.** Against a survey-derived baseline, and where
    no smoother can be fitted because there is no survey.
21. **Which predictors replicate.** Night land-surface temperature, soil
    chemistry, land-use embeddings, prior-round wasting and women's education —
    with sign agreement in all four countries.
22. **What is genuinely null.** Malaria, in every probe. Kept because a
    credible deck retires its own hypotheses.

### Section 6 — Limits (3 slides)
23. **Ranking, not level.** Cross-survey biomarker offsets mean a transported
    prevalence is not validatable. Say it before anyone asks.
24. **Where it fails.** Sierra Leone at 14 districts; women's vitamin A, exactly
    zero in 47–87% of districts; transport burden-targeting at lift 1.007, which
    does **not** reproduce the earlier 1.227.
25. **Four countries.** Transportability is asserted on an n that cannot support
    it. This is the argument for the extension.

### Section 7 — What next (2 slides)
26. **Headroom.** Empirical ceiling 0.47–0.61; current models sit well below it.
    The gap is what further work buys.
27. **The queue.** More countries; HCES/LSMS consumption data; livestock density;
    the protocol as a shared benchmark for the Proxy Modeling Alliance.

### Appendix — kept, not presented
Everything retired below, plus per-cell tables, so questions can be answered
without the main line carrying them.

---

## Disposition of the current 74 slides

**Retire from the main line** (move to appendix, clearly dated):
- The entire "Methods Review 2026-08", "Full Rebuild 2026-08" and "Covariate
  Rebuild 2026-08-30" section structure. Their *content* survives in Section 3,
  reframed as method rather than chronology.
- "Anchoring to regional totals is the design that pays" — withdrawn.
- "Nothing survives FDR once space is accounted for" — the zero was unattainable
  at 999 permutations; the honest statement is about power.
- "The area-level models do not beat doing nothing" — superseded; they beat the
  honest baseline.
- "The ceiling, not the learner, is binding" — rests on the analytic ceiling,
  ~4.7× biased low.
- VMNIS "98% saturated" — two corrections behind.
- Binary AUC slides and the forest plot — the deck's headline metric should be
  rank correlation and burden capture, not AUC, because neither maps to the
  targeting decision.

**Keep and rewrite:**
- "A covariate-free spatial smoother matches 294 covariates" → becomes slide 19,
  with the fold-matched comparison rather than the mismatched one.
- "But at Admin-1 the covariates clearly win" → folds into slide 13.
- "Rank transports across borders. Level does not." → slide 23. This one has
  held up throughout and should be given more room.
- "Model the biomarker, or the binary?" → slide 11.
- "Data we hold and do not use" → slide 6, now largely resolved.

**Keep as is:** the biomarker data gap, the objective, the outcome/population
structure, and the negative-results slide.

## Note for whoever builds it

Every number above is in `results/tables/protocol_v2/`. Before rebuilding the
deck, re-run `02b_merge_and_loco.R` against the freshly rebuilt store so the
slides and the DAG agree — the current protocol-v2 tables were built partly from
the pre-rebuild cache.
