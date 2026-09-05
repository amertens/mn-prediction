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
**what the indicators actually track.** That turned out to be a single
cross-cutting axis rather than a set of nutrient-by-nutrient stories, so Section
2 leads with the axis and treats the per-nutrient detail as the layer on top of
it. The audience leaves with one idea, not six.

## The finding that should lead

**One axis explains most of what we found, and it is not diet.**

Across four countries, six nutrients and 451 indicators, the district-level
markers of micronutrient deficiency converge on a single recognisable thing:
**rural subsistence agro-ecology.** The districts that look agriculturally
self-provisioning — households keeping cattle and goats, cereal-dominant
cropping, legumes grown and eaten, warm night-time temperatures, heterogeneous
soils — are the districts where deficiency is highest, for every nutrient we
measured, in every country we measured it in.

| Outcome | What tracks it most strongly | Direction | Countries agreeing |
|:---|:---|:---|:---|
| **Child vitamin A** | legume consumption; warm night-time temperatures | more legumes → **more** deficiency | 4 / 4 |
| **Child iron** | cattle ownership; cereal-dominant cropping | more cattle → **more** deficiency | 4 / 4 |
| **Women's iron** | warm night-time temperatures; soil aluminium heterogeneity | → more deficiency | 4 / 4 |
| **Women's vitamin A** | child wasting in the same district; rainfall | wasting → more; rain → less | 4 / 4 |
| **Women's B12** | vegetable production; goat ownership | → more deficiency | 3 / 3 |

**Two of these run backwards from individual nutrition, and that is the
finding.** At the household level, eating legumes and owning cattle are
*protective* — they are exactly what dietary diversification programmes promote.
At the district level both point at **more** deficiency, consistently, in every
country. The sign does not flip because the biology is wrong. It flips because
the districts where subsistence agriculture dominates are the poor, rural,
remote, deficient ones, and that overwhelms the dietary contribution when you
aggregate.

This is a textbook **ecological fallacy** signature. It matters practically:
these indicators are excellent for saying **where** to look and useless for
saying **what** to change. The predictive and targeting results in Section 4 are
unaffected — a reliable marker of a deficient district works whether or not it
is causal.

**Why we can see it at all.** We required the same indicator to show the same
sign in four independent national surveys, in countries with different diets,
different agriculture and different survey teams. Under multiplicity correction
for 451 predictors at 14–87 districts, almost nothing clears a conventional
p-threshold — the per-cell screen has median power 0.005 — so a starred p-value
was never going to be the evidence here. Replication is both the stronger claim
and the one that made an inverted sign visible: a wrong-direction result that
holds 4/4 is conspicuous in a way that one at p = 0.03 in a single country is
not.

> **Note on a correction, for the record — not a slide.** An earlier draft of
> this outline had the legume and cattle directions the *right* way round for
> individual biology and therefore the *wrong* way round for the data, and
> called that "the clearest case of recovering known biology." The source scans
> orient the outcome so positive means more deficiency (p4 negates the
> biomarker), and legumes (+4.87) and cattle (+4.40) are both positive in both
> scans. Do not put this correction on a slide; lead with the finding above.

---


## Proposed structure (~30 slides)

### Section 1 — The question (3 slides)
1. **The biomarker data gap.** Keep the existing slide.
2. **What a programme actually decides.** Which districts to reach, with which
   intervention. Frames everything after as a targeting decision.
3. **Two questions we can answer.** *What is associated with deficiency?* and
   *can we predict where it is?* — kept separate all the way through, because
   they have different answers.

### Section 2 — What travels with deficiency (10 slides) — *the centre of the talk*

*Arc: state the axis, show it replicates, show the two indicators that point the
"wrong" way, explain why, then the nutrient-specific layer on top.*

4. **How we looked.** 451 indicators across 19 domains, four countries, six
   nutrients. Plain language: we asked which indicators line up with deficiency
   district by district, and whether the same ones line up in every country.
   State the evidence standard here, once: **same sign, four countries** — not a
   p-value, and say why in one sentence (451 indicators against 14–87 districts).
5. **The finding, in one line.** *The districts that look agriculturally
   self-provisioning are the deficient ones — for every nutrient, in every
   country.* No table yet. Let the room hold one sentence.
6. **The evidence.** The five-row table. Read the last column first: four out of
   four means the pattern held in countries with different diets, different
   farming and different survey teams. Note that the *same domains* — livestock,
   cropping, climate, soil — recur down every row.
7. **Two of these point the wrong way.** Legume consumption is the single
   strongest signal in the whole scan (+4.87, 4/4) and it points at **more**
   child vitamin A deficiency. Cattle ownership does the same for child iron
   (+4.40, 3/3). Put the household-level expectation on the slide next to the
   district-level result so the contradiction is explicit. **Pause here.**
8. **Why the sign flips.** The ecological-fallacy slide, and the most important
   one in the deck. A household that eats legumes is better off; a *district*
   where legume-growing dominates is poorer, more remote and more deficient, and
   that swamps the dietary contribution. Same data, opposite sign, different unit
   of analysis. Say "ecological fallacy" out loud. Then state the rule the rest
   of the talk obeys: **these indicators say where to look, not what to change.**
9. **How we know it is not a coding error.** Two independent scans — the
   deficiency-prevalence scan and the biomarker-concentration scan — give the
   same signs. Four countries agree. We checked the sign convention in the source
   code after the result surprised us. One slide, because the audience will
   privately wonder and it is better to answer than to be asked.
10. **The nutrient-specific layer on top.** Having established the shared axis,
    what genuinely differs: cereal-dominant cropping tracks *more* child iron
    deficiency while root-crop share tracks *less* (4/4 both ways) — the one
    result that reads the way a nutritionist expects, and worth flagging as such
    precisely because the others do not. Women's iron and vitamin A load on
    night-time temperature and soil heterogeneity. B12 on vegetable production
    and goats, three countries only, so mark it the weakest row.
11. **The exception that behaves.** District child wasting tracks women's vitamin
    A deficiency (+4.62, 4/4) and rainfall tracks less of it (−3.92). This one is
    interpretable exactly as stated, because *both* sides are genuine district
    properties — no aggregation trap. Programme reading: where children are
    wasted, women are vitamin A deficient; the two co-locate and may be
    targetable together.
12. **What is consistently null — and worth saying.** No malaria *burden*
    indicator replicates against any outcome. The strongest malaria-domain signal
    is indoor residual spraying *coverage* (family-wise p 0.14), which marks
    where control programmes operate, not where transmission is. A credible talk
    retires its own hypotheses.
13. **What this licenses, and what it forbids.** Licensed: ranking districts,
    choosing where to send a survey team or a programme, and doing it without a
    biomarker survey in hand. Forbidden: reading any row of that table as an
    intervention target. The bridge into Section 3 and 4.

### Section 3 — Domains, and how much each carries (4 slides)
14. **Domain importance by outcome — the matrix.** Rows = 19 domains, columns =
    outcomes, shaded by strength and signed by direction. The single most
    information-dense slide in the deck, and the one people will photograph.
15. **No single domain is decisive.** The signal is *distributed*: many
    indicators each carrying a little, rather than a handful of strong ones.
    Programme implication: there is no shortcut indicator to collect instead of
    a survey.
16. **Diet and agriculture indicators punch above their weight** relative to how
    few of them we have — 11 dietary-diversity indicators against 93
    agricultural ones, and the dietary ones carry more per indicator.
17. **What we still cannot measure well.** Food prices and market access reach
    the model only thinly; household consumption data (HCES/LSMS) is the largest
    unexploited source and is the strongest candidate for the next phase.

### Section 4 — Can we predict where deficiency is? (6 slides)
18. **Three different questions.** Filling gaps inside a surveyed country;
    reaching a region never visited; a country with no survey at all. They have
    different answers and conflating them has caused confusion in this field.
19. **Inside a surveyed country.** Model ranks districts better than the survey's
    own regional averages — the only sub-national alternative currently
    available.
20. **In programme terms.** Visiting the worst-ranked fifth of districts reaches
    26% more of the national deficiency burden than the same effort guided by
    regional survey averages, which are themselves no better than choosing at
    random.
21. **A worked example — Malawi child iron.** Targeted districts carry 43%
    deficiency against 26% nationally; that fifth of districts holds 49% of the
    country's burden in children. *Always name the country and nutrient.*
22. **Transporting to a country with no survey.** Rankings hold — the right
    districts come out in the right order in 21 of 22 tests — but the *level*
    does not transport, so the product is a priority list, not a prevalence map.
23. **What geography alone already gives you.** Neighbouring districts resemble
    each other, and that alone gets much of the way. Covariates earn their place
    in two specific situations: beating a survey-derived baseline, and where
    there is no survey to smooth from.

### Section 5 — Why you can trust these numbers (3 slides) — *was the centrepiece, now support*
24. **We tested our own results hard.** An internal audit re-examined every
    headline claim; several were withdrawn, in both directions.
25. **The two things that mattered most.** Results reported from a single split
    of the data are unstable — one headline moved fourfold on the luck of the
    split. And comparisons are only meaningful when the benchmark has seen the
    same information; one benchmark had already been shown the answers.
26. **What we now report as standard.** Repeat every split and show the spread;
    state what the comparison benchmark could see; show the test could have
    detected an effect if one existed.

### Section 6 — What this means for programmes (3 slides)
27. **Where the model is ready.** Ranked district priority lists, in countries
    with or without their own survey.
28. **Where it is not.** Not a prevalence map; not a replacement for a survey;
    weak where districts are few or deficiency is very rare. And — the hardest
    one to say to this audience — **not a list of things to fix.** The
    indicators that rank districts are markers of rural subsistence, not levers.
    Repeat the Section 2 rule here so nobody leaves the room planning to
    discourage legumes.
29. **What would move it.** More countries; household consumption data; and a
    modest survey — roughly a fifth of a full one — turns a ranking into
    calibrated estimates.

### Appendix
Per-nutrient predictor tables, per-country breakdowns, the full domain matrix,
methods detail, and the withdrawn results with dates. Include the sign
convention explicitly (positive = more deficiency; the concentration scan negates
the biomarker) — it is the kind of thing a careful audience member will ask
about after slide 7, and the kind of thing a re-analyst will get wrong.

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
slides in the deck and were buried at position 385 and 438. Rebuild both around
the shared agro-ecology axis rather than nutrient by nutrient, and add the two
new slides the axis needs: the wrong-way-round result (7) and the
ecological-fallacy explanation (8). "Rank transports
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

**Sign convention — read this before building any slide.** In both
`p1_admin1_scan_predictors.csv` and `p4_admin1_continuous_predictors.csv` the
outcome is oriented so that **positive `meta_z` = associated with MORE
deficiency**. In p4 this is achieved by negating the biomarker concentration
(`y = -wmean(t, w)` in the script), because a lower concentration is worse
status. Getting this backwards is exactly the error the first draft of this
outline made.
