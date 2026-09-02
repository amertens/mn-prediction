# Slide deck outline — nutrition policy audience, September 2026

For `docs/mn_prediction_slides.qmd`. Supersedes the methods-first draft of this
outline. Companion to `docs/findings/TWO_READINGS_2026-09d.md` and
`docs/findings/PROTOCOL_V2.md`.

## Who this is for, and what changed

The earlier version of this outline made **evaluation protocol** the centrepiece.
That is the most transferable thing the project has produced and it belongs in a
methods paper — but it is not what a nutrition policy audience came to hear.
They want to know **what travels with deficiency, whether it differs by
nutrient, and what that means for programmes.**

So the protocol work drops to a short credibility section — three slides that
answer "why should I believe these numbers" — and the centre of the deck becomes
**variable and domain importance, outcome by outcome.**

## The finding that should lead

Across four countries, the factors associated with deficiency are
**nutrient-specific and mechanistically coherent**, and they **agree in sign in
all four countries**:

| Outcome | What tracks it most strongly | Direction | Countries agreeing |
|:---|:---|:---|:---|
| **Child vitamin A** | legume consumption; fruit & vegetable consumption | more → less deficiency | 4 / 4 |
| **Child iron** | cereal-dominant cropping; cattle ownership | cereal-dominant → more | 4 / 4 |
| **Women's iron** | soil chemistry (aluminium heterogeneity, CEC); low maternal BMI | → more | 4 / 4 |
| **Women's vitamin A** | child wasting in the same district; rainfall | wasting → more; rain → less | 4 / 4 |
| **Women's B12** | household assets; education | lower → more | 3 / 4 |

**Why this matters more than a p-value.** Under a correction for having examined
451 predictors, few of these clear conventional significance thresholds — the
multiplicity penalty at 14–87 districts per country is severe. The evidence that
should carry the talk is **replication**: the same variable, with the same sign,
in four independent national surveys, in countries with different diets,
different agriculture and different survey teams. That is a stronger and more
honest claim than a starred p-value, and a nutrition audience will recognise it
as such.

---

## Proposed structure (~30 slides)

### Section 1 — The question (3 slides)
1. **The biomarker data gap.** Keep the existing slide.
2. **What a programme actually decides.** Which districts to reach, with which
   intervention. Frames everything after as a targeting decision.
3. **Two questions we can answer.** *What is associated with deficiency?* and
   *can we predict where it is?* — kept separate all the way through, because
   they have different answers.

### Section 2 — What travels with deficiency (9 slides) — *the centre of the talk*
4. **How we looked.** 451 indicators across 19 domains, four countries, six
   nutrients. Plain language: we asked which indicators line up with deficiency,
   and whether the same ones line up in every country.
5. **The headline table.** The five-row table above. One slide, no statistics.
6. **Vitamin A in children is a diet story.** Legume consumption is the single
   strongest signal in the whole scan (consistent in 4/4 countries); fruit and
   vegetable consumption follows. Mechanistically exactly what provitamin A
   carotenoid intake predicts, and the clearest case in the project of the data
   recovering known biology.
7. **Iron in children is a cropping-system story.** Cereal-dominant cropping
   tracks *more* deficiency; cattle ownership tracks *less*. Consistent with
   phytate inhibition in staple-heavy diets and animal-source iron availability.
   Connect explicitly to fortification and dietary-diversification programming.
8. **Iron in women is a soil and maternal-nutrition story.** Soil chemistry
   (aluminium heterogeneity, cation exchange capacity) and low maternal BMI.
   Flag honestly: soil is a *marker* of agro-ecology here, not a demonstrated
   causal pathway.
9. **Vitamin A in women tracks general undernutrition.** District-level child
   wasting is the strongest single predictor across the entire binary scan
   (4/4 countries); rainfall is protective. Programme reading: where children
   are wasted, women are vitamin A deficient — the two co-locate and may be
   targetable together.
10. **B12 and folate are socioeconomic.** Household assets and education
    dominate; this is the outcome set where the wealth gradient is clearest.
11. **What is consistent across nutrients.** Agro-ecology (land use, climate,
    soil) appears for nearly every outcome; night-time land surface temperature
    tracks worse status for vitamin A and women's iron in all four countries.
12. **What is consistently null — and worth saying.** Malaria indicators show
    essentially nothing, in every test we ran, for every outcome. A credible
    talk retires its own hypotheses.

### Section 3 — Domains, and how much each carries (4 slides)
13. **Domain importance by outcome — the matrix.** Rows = 19 domains, columns =
    outcomes, shaded by strength and signed by direction. The single most
    information-dense slide in the deck, and the one people will photograph.
14. **No single domain is decisive.** The signal is *distributed*: many
    indicators each carrying a little, rather than a handful of strong ones.
    Programme implication: there is no shortcut indicator to collect instead of
    a survey.
15. **Diet and agriculture indicators punch above their weight** relative to how
    few of them we have — 11 dietary-diversity indicators against 93
    agricultural ones, and the dietary ones carry more per indicator.
16. **What we still cannot measure well.** Food prices and market access reach
    the model only thinly; household consumption data (HCES/LSMS) is the largest
    unexploited source and is the strongest candidate for the next phase.

### Section 4 — Can we predict where deficiency is? (6 slides)
17. **Three different questions.** Filling gaps inside a surveyed country;
    reaching a region never visited; a country with no survey at all. They have
    different answers and conflating them has caused confusion in this field.
18. **Inside a surveyed country.** Model ranks districts better than the survey's
    own regional averages — the only sub-national alternative currently
    available.
19. **In programme terms.** Visiting the worst-ranked fifth of districts reaches
    26% more of the national deficiency burden than the same effort guided by
    regional survey averages, which are themselves no better than choosing at
    random.
20. **A worked example — Malawi child iron.** Targeted districts carry 43%
    deficiency against 26% nationally; that fifth of districts holds 49% of the
    country's burden in children. *Always name the country and nutrient.*
21. **Transporting to a country with no survey.** Rankings hold — the right
    districts come out in the right order in 21 of 22 tests — but the *level*
    does not transport, so the product is a priority list, not a prevalence map.
22. **What geography alone already gives you.** Neighbouring districts resemble
    each other, and that alone gets much of the way. Covariates earn their place
    in two specific situations: beating a survey-derived baseline, and where
    there is no survey to smooth from.

### Section 5 — Why you can trust these numbers (3 slides) — *was the centrepiece, now support*
23. **We tested our own results hard.** An internal audit re-examined every
    headline claim; several were withdrawn, in both directions.
24. **The two things that mattered most.** Results reported from a single split
    of the data are unstable — one headline moved fourfold on the luck of the
    split. And comparisons are only meaningful when the benchmark has seen the
    same information; one benchmark had already been shown the answers.
25. **What we now report as standard.** Repeat every split and show the spread;
    state what the comparison benchmark could see; show the test could have
    detected an effect if one existed.

### Section 6 — What this means for programmes (3 slides)
26. **Where the model is ready.** Ranked district priority lists, in countries
    with or without their own survey.
27. **Where it is not.** Not a prevalence map; not a replacement for a survey;
    weak where districts are few or deficiency is very rare.
28. **What would move it.** More countries; household consumption data; and a
    modest survey — roughly a fifth of a full one — turns a ranking into
    calibrated estimates.

### Appendix
Per-nutrient predictor tables, per-country breakdowns, the full domain matrix,
methods detail, and the withdrawn results with dates.

---

## Disposition of the current 74 slides

**Retire from the main line** (appendix, clearly dated): the three accreted
review/rebuild sections as *structure*; "Anchoring to regional totals is the
design that pays" (withdrawn); "Nothing survives FDR once space is accounted
for" (the zero was unattainable by construction); "The area-level models do not
beat doing nothing" (superseded); "The ceiling, not the learner, is binding"
(rests on a ceiling ~4.7x biased low); VMNIS "98% saturated" (two corrections
behind); the binary AUC slides and forest plot (AUC does not map to a targeting
decision).

**Keep and rewrite:** "Cross-Outcome Domain Importance" and "Individual Variable
Importance" become Section 2 — they were always the most audience-relevant
slides in the deck and were buried at position 385 and 438. "Rank transports
across borders. Level does not." → slide 21, given more room. "A covariate-free
spatial smoother matches 294 covariates" → slide 22, with the fold-matched
comparison. "Data we hold and do not use" → slide 16.

**Keep as is:** the biomarker data gap, the objective, the outcome/population
structure, and the negative-results slide.

## Note for whoever builds it

Numbers for Section 2 are in `results/tables/signal_probes/`
(`p1_admin1_scan_domains.csv`, `p4_admin1_continuous_domains.csv`, and the
matching `_predictors` files, which carry the per-outcome sign agreement).
Section 4 comes from `results/tables/protocol_v2/`. Present sign agreement
across countries, not p-values: under correction for 451 predictors at 14-87
districts, few individual indicators clear conventional thresholds, and
replication is both the stronger evidence and the more honest claim.
