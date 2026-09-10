# 15-minute deck for a nutrition policy audience

**Draft outline for review. Not the deck.**
Date: 9 September 2026 · Branch: `signal-audit-and-protocol-v2`
Revision 4. Slide 9 confirmed visible; the NCE request read and folded in;
slides 8, 11 and 12 substantially rewritten; the Côte d'Ivoire section updated
to completed status.

Built from the motivation and introduction slides of
`MN-proxy-BMGF-short presentation.pptx` and `MN-proxy-Ghana-presentation.pptx`
(Dropbox/MN prediction), with every result number refreshed from
`results/tables/protocol_v2/` and `results/tables/signal_probes/` as of the
8 September re-runs (RR-09 / VA-01).

Framing, wording and next-steps language follow the NCE request
(`NCE_Proxy_2026-09-02 ANM.docx`, INV-070689) together with
`docs/findings/NCE_IMPLICATIONS_2026-09-03.md`, which quotes the draft sentence
by sentence and gives agreed replacement wording. Where the NCE has a settled
phrasing for a claim, this outline uses it rather than inventing one.

**The deck should read as the public face of the NCE.** Its Section 6 is the
NCE's stated central challenge, and its Section 8 is the NCE's activity list.

Companion documents: `docs/findings/TWO_READINGS_2026-09f.md` (what survives,
stated two ways) and `docs/findings/SLIDE_OUTLINE_2026-09.md` (the longer
30-slide version).

---

## Two things in the NCE draft you should look at, unrelated to the slides

**1. A sentence is corrupted.** In the second paragraph:

> *"…produce national and regional-level vitamin and mineral deficiency (VMD)
> prevalence estimates substantially more accurate than the current default of
> 2) it gives programs information for deciding which areas to prioritize…"*

Text has been lost between "the current default of" and "2)", and the "1)" that
"2)" refers to is gone with it. As written the sentence does not parse. This is
worth fixing before submission, because it sits in the achievements paragraph,
which is the most-read part of the request.

**2. The cost-analysis bullet may overreach, and it is the one place where the
NCE and the evidence disagree.** The bullet proposes to *"compare the
anticipated cost of conducting a planned national nutrition survey using the
standard sampling approach with the cost of using the reduced sample size
identified through the proxy modeling approach."*

However, G4-02 tested this symmetrically and found that the model **does not
reduce the sample a survey needs**. District error degrades in lockstep with the
survey's as the sample shrinks (+3.36 vs +3.33 pp from full to 15%), and
NCE_IMPLICATIONS item 2 says in terms: *do not frame savings as coming from a
smaller sample.*

There **is** a genuine reduced-sample story, but it is narrower than the bullet
implies. Under the **anchor-and-rank design**, a country measures one national
prevalence at roughly 5% of a full sample and then lets the transported ranking
order the districts; that gives district estimates a conventional survey only
matches at about eight times the sample. It is a *different survey design*, not a
smaller version of the same one. **I would rewrite the bullet around that design
specifically.** Otherwise the cost analysis is scoped to answer a question we
have already answered in the negative. Slide 16 uses the narrower framing.

---

## What changed in revisions 2, 3 and 4

| Slide | Change | Why |
|:---|:---|:---|
| **8** | **Rewritten.** No longer "an ensemble competition." | The NCE says explicitly not to call the ensemble the model: four meta-learner losses were tested and none beat the zero-tuning index at these sample sizes. The old decks' ensemble cartoon is now actively misleading. |
| **9** | **Kept visible, hedge removed.** | Your call. |
| **10** | Rewritten with the NCE's scoping language and the per-cell honesty it requires. | The "26%" flag from revision 1 is resolved. See the note on that slide. |
| **11** | Reframed: the axis is agro-ecological, and the survey-derived indicators are the *weak* ones. | Follows from the NCE's "what to collect" item. |
| **12** | **Substantially rewritten.** No longer asserts ecological fallacy. | Your challenge was right. See "On slide 12" below. |
| **16** | Uses the NCE's settled sentence verbatim. | It is better than mine. |
| **18** | Adds the NCE's survey-size correction. | The draft was softened for a reason, and the slide should not undo it. |
| **19** | **Rewritten again** around the NCE's actual activity list. | Revision 2 guessed; the NCE names the activities and the money behind them. |

**Revision 3 (after reading the NCE itself)** additionally changed: slide 4 (five
country databases, four with biomarkers, Côte d'Ivoire being the fifth and the
no-survey test case), slide 6 (reconciles ">1,000 variables" in the NCE with the
454 that reach the models), slide 7 (the software pipeline is a named NCE
deliverable and deserves a bullet), slide 14 (elevated, because it *is* the
NCE's central challenge in the NCE's own words), and slide 17 (the dashboard is
"ready for field testing" with invited users, which answers the deployment
question).

**Revision 4** replaces the Côte d'Ivoire "in progress" section with what
actually happened, and records the two guards the result passed.

---

## On slide 12: you were right to push, and here is where I landed

**No, it is not for sure ecological fallacy, and we should not say it is.** I
over-committed in revision 1. There are four candidate explanations, they carry
different programmatic implications, and, most importantly, **we have the
individual-level data to distinguish two of them.**

### What we can actually test

I checked the individual participant data. Availability is uneven but not empty:

| Country | Individual legume/bean intake | Individual cattle / milk | Testable? |
|:---|:---|:---|:---|
| **Sierra Leone** | `gw_cBeans`, `gw_wBeans`, `gw_FoodGrp2Beans` | `gw_hCattleNumb`, `gw_cMilk`, `gw_cMilkTimes` | **Both, fully** |
| **The Gambia** | none found | `gw_cattle`, `gw_cow`, `gw_goat`, `gw_cMilk`, `gw_cDairy` | **Cattle only** |
| Ghana | none | formula milk only | No |
| Malawi | district aggregates only | district aggregates only | No |

The household-level association can therefore be estimated directly in Sierra
Leone for both indicators, and in The Gambia for cattle. That is a real test
rather than a thought experiment.

*(A note on scope: these are `gw_` variables, drawn from the biomarker survey's
own questionnaire. Using them as **predictors** is barred by the proxy-only
rule, but using them for a **descriptive cross-level check** is a different
thing and does not touch the prediction models. Worth stating explicitly so it
does not look like a violation.)*

### The four explanations

**1. Ecological confounding, the "fallacy" reading.** Legume-heavy and
cattle-keeping districts are poorer and more remote, with lower service access,
so the district-level association may be that rather than the food. This is
*untested for these two indicators specifically.* We tested urbanicity for the
climate+soil **index**, which moves it by less than 0.02, but never partialled
it out of legumes or cattle on their own.

**2. Cross-level substitution, where the aggregate variable means something
different from the individual one.** This is the explanation I now think most
likely for legumes, and it is **not** a fallacy. At the individual level, "the
child ate legumes yesterday" is one tick on a diversity score, and protective.
At the district level, "legumes are what the food system produces and people
eat" means **a plant-dominant, animal-source-food-poor food system.** Vitamin A
is exactly the nutrient that would expose this, because preformed retinol comes
almost entirely from liver, dairy and eggs, and provitamin A carotenoids need
dietary fat to absorb. A high district legume share is therefore a marker for
*low animal-source food availability*. Reading it as "legumes cause deficiency"
is the error; reading it as "this is an ASF-poor food system" is correct.

**3. Genuine household-level biology running the same direction, specifically
for cattle and iron.** This one is worth taking seriously. Fresh cow's milk is
a well-established iron-deficiency risk in 6–59 month olds: it is iron-poor, its
calcium and casein inhibit non-heme iron absorption, and in infants it can cause
occult gastrointestinal blood loss. In agro-pastoral settings cattle are wealth
and traction stores that are rarely slaughtered, so **cattle means milk, not
meat.** If that is what is happening, the association is causal at the
individual level too, is not an aggregation artefact at all, and is
*programmatically actionable* in a way nothing else in this talk is.

**4. Agro-ecological zone confounding.** Legume cropping and cattle keeping both
concentrate in the drier savanna zones: northern Ghana, northern and eastern
Sierra Leone, the Gambian interior. Those zones also carry higher malaria,
higher helminth burden, lower service density and lower dietary diversity. On
this reading, legumes and cattle are noisy, survey-dependent proxies for the
same climate-and-soil axis that carries all the transport signal, which would
neatly explain why the two-domain index transports *better* than the full
vocabulary that contains them.

Explanations 2 and 4 are cousins and point the same programmatic direction.
Explanation 3 points somewhere completely different.

### My hypothesis for these four countries, stated so it can be wrong

- **Legumes and child vitamin A: explanation 2 plus 4.** The cowpea and
  groundnut belts are the ASF-poor drylands. **Prediction:** the association
  attenuates sharply when conditioned on district animal-source-food consumption
  or on the climate+soil index, and the Sierra Leone individual-level
  association is **null or protective**.
- **Cattle and child iron: explanation 3 is live.** **Prediction:** unlike
  legumes, cattle shows a **positive (harmful) association at individual level
  too**, mediated by fresh cow's milk intake in the 6–23 month olds. That is
  testable directly in Sierra Leone through `gw_cMilk` and `gw_cMilkTimes`, and
  partially in The Gambia.

**These two predictions differ, which is what makes the check worth running.**

### What I recommend for the deck

Do not say "ecological fallacy" as a conclusion. Slide 12 below now presents the
contradiction, names the explanations, and says which one we are testing, which
is both more honest and, for this audience, more interesting. It also protects
you from the one bad outcome, namely someone acting on "discourage legumes."

**Offer:** the Sierra Leone and Gambia individual-level check is a scoped
analysis, survey-weighted, cluster-robust and child-age-stratified with the
biomarker as the outcome, and it would settle explanation 3 in particular. It is
maybe half a day. It is not in the current pipeline and involves choices (age
bands, whether to condition on wealth, how to handle the milk-frequency coding)
that I would rather you make than guess. **Say the word and I will run it before
the talk.** If it confirms the cattle-milk mechanism, that becomes the most
quotable slide in the deck.

---

## Three framing choices carried forward from revision 1

1. **The talk is about *where*, not *what*.** Every result supports ranking
   districts, and none supports "fix X to reduce deficiency." This is stated on
   slides 12 and 18. *(Slide 12's version is now softer and better: the rule
   holds whether or not the cattle mechanism turns out to be real, because a
   ranking indicator does not become an intervention target either way.)*
2. **Slide 9 stays visible.** Geography does most of the work inside a surveyed
   country, and saying so unprompted buys the credibility that the
   transportability section spends.
3. **Numbers from the old decks that no longer reproduce are listed at the end**
   rather than silently dropped.

**Timing.** Nineteen slides for 15 minutes, workable if slides 5, 8, 13 and 16
run at 20 seconds. Two slides are marked **PAUSE**. If you must cut, cut 6 and
13 (*cuttable*), which takes you to 17.

---

## Section 1. Motivation (3 slides, ~2.5 min)

### 1. Title

Proxy modelling for vitamin and mineral deficiencies.
Andrew Mertens, UC Berkeley. Gates Foundation funding acknowledgement.
Collaborators line from the existing decks (GroundWork, University of Ghana,
UW–Madison, KEMRI–Wellcome Trust, UNICEF, Global Affairs Canada).

*Reuse the existing title slide as-is.*

### 2. The problem: programmes decide sub-nationally, data exists nationally

**Message:** biomarker surveys are rare, expensive, and, where they exist,
usually only powered to the region rather than the district.

- Blood-based micronutrient surveys cost millions and happen roughly once a
  decade per country. Most of sub-Saharan Africa has none in the last ten years.
- Meanwhile the decisions are district-level: where fortification, or
  supplementation, or a survey team, actually goes.
- Concretely, in our four surveys: **206 districts, and the survey design puts
  1–4 sample clusters in each one.** In Ghana, 83% of districts hold a single
  cluster. The survey was never designed to answer the district question.
- **And the obvious fallback does not work.** For a district the survey did not
  reach, using its region's average is **no better than choosing at random**.
  That is the gap this project is trying to fill.

*Source: `targets_v2.csv` and `variance_components_ceiling.csv`
(`mean_clusters_per_unit`, `share_single_cluster`). The regional-average point
is the NCE's own framing (NCE_IMPLICATIONS item 1), backed by
`nce_targeting_summary.csv`.*

*Keep the photo and framing from Ghana deck slide 2 and BMGF deck slide 20.*

### 3. The idea, in one sentence

**Message:** other things *are* measured everywhere, every year. Can they tell
us where deficiency is?

- Satellites photograph every district monthly, and soil has been mapped
  globally. Household surveys, disease models, livestock censuses, crop maps and
  market prices all exist at district resolution, for every country, updated
  regularly.
- **The question:** do these "proxy" measures line up with measured deficiency
  well enough to rank districts in a country we have never sampled?
- Two questions we keep separate throughout, because they have different
  answers:
  - *What travels with deficiency?* (Section 5)
  - *Can we predict where it is?* (Sections 4 and 6)

*Carry over the Goals block from BMGF deck slide 2, compressed to these two
questions. The three-bullet goals list is too abstract for 15 minutes.*

---

## Section 2. The four surveys (2 slides, ~1.5 min)

### 4. Five country databases; four with blood to check against

**Message:** these are the only places on earth where we can check whether the
proxies are right, and the fifth country is the test of whether we need them.

| Survey | Year | Districts | Outcomes used here |
|:---|:---|---:|:---|
| Sierra Leone | 2013 | 14 | iron, vitamin A, folate, B12 |
| Malawi (DHS-linked) | 2015–16 | 87 | iron, vitamin A, zinc, folate, B12 |
| Ghana | 2017 | 75 | iron, vitamin A, folate, B12 |
| The Gambia (MICS-linked) | 2018 | 30 | iron, vitamin A |
| **Côte d'Ivoire** | none | 33 | **no biomarker survey; the prediction target** |

- Children 6–59 months and women 15–49 y, modelled separately.
- **24 country × outcome combinations** with biomarkers, which is our unit of
  evidence. Every number later in this talk is an average over those 24, or over
  the 22 that support a cross-country test.
- Gold-standard blood biomarkers, inflammation-adjusted, population-based
  sampling. This is the ground truth everything else is scored against.
- **Côte d'Ivoire has the full proxy database and no biomarker survey.** It is
  what the whole method is for, and it comes back on slides 14 and 17.

*Source: `targets_v2.csv`. The five-country database framing is the NCE's own,
from the achievements paragraph. Simplify from Ghana deck slide 4: drop the
per-survey prose and keep the district counts, which are the number that matters
for everything after.*

### 5. What "deficient" means here (20 seconds)

**Message:** standard definitions, nothing invented.

Compressed from BMGF deck slide 25 and Ghana slide 5: ferritin <12 / <15 µg/L,
BRINDA-adjusted; retinol <0.70 µmol/L after each survey's own RBP calibration;
folate <10 nmol/L; B12 <150 pmol/L; zinc by the time-of-day cutoffs.

> **Flag for you, and it is bigger than revision 2 said.** Vitamin A changed in
> the 8 Sep re-run (VA-01). Each survey's RBP-to-retinol calibration differs
> enough that a fixed 0.70-of-RBP cutoff was worth 0.61–0.91 of retinol across
> the four countries. Child vitamin A prevalence is now 20 / 28 / 7 / 0.4%
> (Gambia / Ghana / Sierra Leone / Malawi).
>
> **Checking women's vitamin A turned up that six of the 24 cells are too sparse
> to carry a district signal, not two.** I would put a one-line version of that
> on this slide, along the lines of "some outcomes are too rare to map at
> district level, and we say which," because it is the kind of limitation this
> audience respects and the dashboard makes it visible anyway.

#### Diagnosis: women's vitamin A, and how many cells are actually usable

**The Côte d'Ivoire women's vitamin A panel is not broken.** It predicts
0.8–2.1% because that is what the training data looks like: observed district
prevalence averages 3.5% (Gambia), 3.1% (Ghana), 0.9% (Sierra Leone) and
**0.2%** (Malawi). The model is faithfully reproducing a near-absent outcome.

The real problem is upstream, and it is a **zero-inflation** problem rather than
a level problem. Counting districts with *exactly zero* measured deficiency:

| Country | Outcome | Districts | Zeros | % zero | Mean prevalence |
|:---|:---|---:|---:|---:|---:|
| Malawi | women's vitamin A | 87 | 85 | **98%** | 0.22% |
| Malawi | child vitamin A | 87 | 82 | **94%** | 0.30% |
| Ghana | women's vitamin A | 75 | 55 | 73% | 3.07% |
| Sierra Leone | women's B12 | 14 | 10 | 71% | 0.42% |
| Ghana | women's B12 | 75 | 47 | 63% | 9.37% |
| Sierra Leone | women's vitamin A | 14 | 7 | 50% | 0.90% |

**Two cells are degenerate (≥94% zeros) and six are sparse (≥50% zeros).** At
the other end, Gambia's two iron cells, Malawi women's zinc and Sierra Leone
women's folate have no zero districts at all and are well behaved.

**What VA-01 did to it.** Comparing targets before and after commit `041ec7e`:

| Cell | Mean prevalence before → after | Zero districts before → after |
|:---|:---|:---|
| Malawi women's vitamin A | 1.60% → **0.22%** | 76 → **85** of 87 |
| Malawi child vitamin A | 8.68% → **0.30%** | 38 → **82** of 87 |
| Ghana women's vitamin A | 1.70% → 3.07% | 63 → 55 of 75 |
| Ghana child vitamin A | 13.6% → 26.1% | 21 → 8 of 75 |

The retinol rebuild therefore moved Ghana's vitamin A cells *up* and made them
more usable, while moving Malawi's down into degeneracy. That is a real
consequence of a defensible correction rather than a bug, but it means
**Malawi's vitamin A rows should carry a "too rare to rank" marker instead of a
number**, both in the deck and in the dashboard.

**Three things I would do about it, none of which is on a slide:**

1. Mark the six sparse cells explicitly in `targets_v2.csv` and have the
   benchmark tables report them separately, so that a mean over "24 cells" is
   not quietly a mean over 18 informative ones plus six near-constants.
2. Check Ghana child vitamin A, where one district now reads **100%**
   prevalence. It is almost certainly a small-denominator district, and it will
   drag any mean.
3. Log this as a sandbox entry. It belongs with VA-01 rather than only in a
   slide outline.

*Source: `results/tables/protocol_v2/targets_v2.csv` against
`git show 041ec7e~1:results/tables/protocol_v2/targets_v2.csv`.*

---

## Section 3. Data sources and building the proxy dataset (2 slides, ~2 min)

### 6. What we assembled *(cuttable)*

**Message:** over a thousand variables assembled per country, of which 454
survive into the models. All are routinely collected and none requires a blood
draw.

> **Reconcile with the NCE before this slide is built.** The NCE says each
> country database holds **"more than 1,000 variables."** That is right for the
> assembled database. The number that reaches the models is **454**, after
> dropping near-constant and heavily missing variables, anything that could leak
> the answer, and duplicates across sources. **Say both numbers on the slide**,
> as ">1,000 assembled, 454 into the models", or the deck and the proposal will
> appear to disagree.

| Source | Indicators | What it is |
|:---|---:|:---|
| DHS / MICS district aggregates | 130 | health, WASH, assets, diet, livestock, education |
| Google Earth Engine | 76 | climate, vegetation, night lights, built surface, water |
| AlphaEarth satellite embedding | 64 | a learned 64-number summary of what each district *looks like* |
| SoilGrids / iSDA | 45 | soil nutrients, texture, carbon, and their heterogeneity |
| IHME disease surfaces | 17 | modelled disease and nutrition burden |
| WorldPop / GHS | 14 | population density, urbanisation |
| MapSPAM | 12 | crop mix and production |
| Malaria Atlas | 11 | transmission, interventions, blood disorders |
| Livestock (GLW4) | 8 | cattle, goat, poultry density |
| Köppen/AEZ, JRC water, ESPEN helminths, WFP prices | 18 | agro-ecology, coast and water proximity, worms, market prices |

*Source: `source_ablation_loco_summary.csv` (`n_cols`) and
`domain_representation_summary.csv`.*

### 7. Building the dataset is most of the work

**Message:** the hard part is not the model. It is making four countries' worth
of incompatible data line up on the same map, in the same year.

- Every layer is matched to **each survey's own fieldwork window**, so Gambia
  Jan–Apr 2018 and Malawi 2016 rather than the publication year. One table
  (`metadata/survey_years.csv`) governs this; getting it wrong silently degrades
  everything downstream, and did, twice.
- District boundaries change between census rounds and between sources.
  Everything is reconciled onto one common administrative rung per country.
- **Anything that could leak the answer is removed by rule rather than by
  judgement.** DHS anaemia prevalence and nutrient-named indicators are
  excluded, because a model that predicts iron deficiency from measured anaemia
  has learned nothing transferable.
- Keep the "common challenges" bullet from the existing decks (fragmentation
  across platforms, access hurdles, mismatched geographic and temporal
  resolution, no standard metadata). It resonates most with the data-producer
  half of the room.
- **End on the deliverable rather than the difficulty.** All of this now runs as
  **a software pipeline that automates the aggregation and can be pointed at a
  new country.** That is a named achievement in the NCE, and it is the reason
  adding Ethiopia or Pakistan is a matter of weeks rather than the year it took
  the first time. *A policy audience cares about this more than about the merge
  logic; it is the difference between a study and a capability.*

*Reuse Ghana deck slide 7, adding the fieldwork-matching, leakage-rule and
pipeline bullets, all new since those decks.*

---

## Section 4. Model approach and performance within countries (3 slides, ~3 min)

### 8. How the model works, and how the complicated version lost (20 seconds)

**REWRITTEN. Do not reuse the ensemble cartoon from the old decks as the
description of the model.**

**Message:** we tried the sophisticated approach and the simple one won. That is
a finding, not a shortcut.

- We built what the earlier decks promised: a twelve-learner ensemble, machine
  learning competing against regression, weighted by out-of-sample performance.
- **Under corrected evaluation, it does not win.** Four different ensemble
  meta-learners were tested, and **none beat a simple index built from summary
  scores of each data domain** at these sample sizes.
- The reason is sample size, and it is worth one sentence to a policy audience:
  **we have 14 to 87 districts per country, not thousands of people.** Flexible
  methods need data to be flexible with, and with 87 rows they overfit.
- So the model we report is **a parsimonious index over domain summaries**:
  simple, inspectable, and better than the complicated alternatives on data none
  of them had seen.

*Source: NCE_IMPLICATIONS "One thing not to say"; `sl_domain_index_scores.csv`,
`sl_rank_loss_scores.csv`, `estimator_tournament_truth.csv`.*

> **This replaces the ensemble slide in both old decks.** If you want to keep
> the ensemble cartoon visually, use it as the "what we tried" half of a
> two-panel slide, with the index as the "what won" half. Do not present the
> ensemble as the current model; the NCE is explicit about this and the deck
> should match the proposal.

### 9. Inside a country with a survey: how good is it? **KEEP VISIBLE**

**Message:** we can fill gaps in a surveyed country, and here is exactly what is
doing the work.

Ranking districts on the biomarker level, averaged over the 24 country × outcome
combinations (correlation with truth; higher is better):

| Approach | Score |
|:---|---:|
| Model (covariates + geography) | **0.40** |
| Model (covariates only) | 0.40 |
| **Geography alone, no covariates at all** | 0.39 |
| The survey's own regional averages | 0.31 |
| National average (no information) | −0.23 |

- **The good news:** covariates beat the best survey-based alternative available
  (0.40 vs 0.31), in 15 of 18 testable cells.
- **The honest news:** *inside a surveyed country*, a simple map smoother that
  knows only which districts are neighbours reaches 0.39. Neighbouring districts
  resemble each other, and that alone gets most of the way.
- The covariates' real value is therefore where there is nothing to smooth
  *from*, which is Section 6, and which is the whole point.

*Source: `benchmarks_v2_summary.csv`, `estimand == "infill"`, `target ==
"level"`. On prevalence the ordering is more favourable to covariates: 0.29
model against 0.26 geography-only and 0.22 regional averages.*

*Speaker note: say the third row out loud rather than letting people find it.
"You should ask whether the satellites are doing anything a map couldn't. Inside
a surveyed country, barely. Here is where they earn it." Then go to Section 6.*

### 10. What that buys a programme

**Message:** translate a correlation into a targeting decision, and scope the
claim the way the NCE scopes it.

- Suppose you can only reach the worst-ranked fifth of districts.
- Using the model's ranking, those districts carry **32% average deficiency
  against 24% nationally**, and hold **22% of the national burden**. Using the
  survey's own regional averages: 20%. At random: 14%. A perfect oracle: 47%.
- **Worked example, Malawi child iron:** targeted districts carry **43%
  deficiency against 26% nationally**, and that fifth of districts holds **49%
  of the country's burden in children.** *Always name the country and nutrient
  when quoting a number this good.*
- **The scoped claim, in the NCE's words:** *for districts a survey has not
  reached, regional averages are no better than chance, while the model is.*

*Source: `nce_targeting_summary.csv`, `estimand == "infill"`.*

> **The revision-1 "26%" flag is resolved.** The NCE draft's figure is 24%
> against 19% burden capture, the current tables give 22% against 20%, and the
> "26% more" in the long outline was the ratio of an earlier pair. **More
> important than the drift is what NCE_IMPLICATIONS item 1 says about it:** per
> cell the model beats the jackknifed regional mean in only **8–11 of 18**, the
> pooled margin comes from countries whose regions hold just 3–5 districts, and
> a regional mean that *includes* the district reaches 30%. **The robust claim
> is the ranking claim (0.29 against 0.19, 15 of 18), not the burden ratio.** I
> have written the slide with raw capture percentages and the scoped sentence,
> and I would not put an "N% more than regional averages" headline on a slide at
> all.

---

## Section 5. The variables driving the best model (3 slides, ~3 min)

### 11. What matters is agro-ecological, and the survey indicators are the weak ones. **PAUSE**

**Message:** across four countries, six nutrients and 454 indicators, the signal
that travels is the remotely sensed one.

Domains that measurably improve cross-country prediction (drop-one test, 22
cells; the number is how much prediction falls when the domain is removed):

| Domain | On biomarker level | On prevalence |
|:---|---:|---:|
| Soil characteristics | +0.024 | +0.036 |
| Climate and weather | +0.015 | +0.021 |
| Agriculture, crops, land use | +0.011 | +0.026 |
| Malaria burden and control | +0.009 | +0.023 |
| Livestock density | +0.007 | none |
| *Everything else* | *≈ 0 or negative* | *≈ 0 or negative* |

- **The headline for a data-strategy audience:** dropping the entire DHS/MICS
  block **improves** cross-country prediction (0.28 → 0.33). Ten survey-derived
  domains are dead weight across borders, and **the indicators that travel are
  the ones nobody has to collect**, namely soil and climate, free, global and
  updated.
- Individual indicators holding the **same direction in all four countries**, a
  stronger standard than a p-value with 454 indicators and 14–87 districts:
  night-time temperature, soil aluminium heterogeneity, legume consumption,
  cattle ownership, cereal-dominant cropping, district child wasting.
- Two conspicuous absences: **malaria burden**, where only spraying *coverage*
  shows up, marking where programmes operate rather than where transmission is;
  and **direct nutrition indicators** such as stunting, wasting and dietary
  diversity as a domain. A credible talk retires its own hypotheses.

*Source: `domain_ablation_loco_summary.csv`, `source_ablation_loco_summary.csv`,
`signal_probes/p1_*_predictors.csv` and `p4_*_predictors.csv` (`sign_agree`);
framing per NCE_IMPLICATIONS item 6.*

### 12. Two of these point the wrong way, and we can find out why. **PAUSE**

**REWRITTEN. No longer asserts ecological fallacy.**

**Message:** present the contradiction, name the competing explanations, and say
which one we are testing. Do not resolve it on the slide.

- **Legume consumption is the strongest replicated signal in the entire scan
  (4 of 4 countries), and it points at MORE child vitamin A deficiency.**
- **Cattle ownership does the same for child iron deficiency.**
- At the household level both are protective, and both are what dietary
  diversification programmes promote. At the *district* level both mark
  deficiency, in every country.

**Three things this could be.** Put them on the slide; they are the interesting
part.

1. **The variable means something different once you aggregate.** "This child ate
   legumes" is dietary diversity. "This district eats legumes" is a
   **plant-dominant, animal-source-food-poor food system**, and vitamin A comes
   almost entirely from liver, dairy and eggs. Same word, different quantity.
2. **Both are markers of dryland agro-ecology.** Legume cropping and cattle
   keeping concentrate in the same savanna zones as high malaria, low service
   density and low dietary diversity, which are the zones the soil and climate
   layers already identify.
3. **For cattle, it may be real biology at the household level too.** In
   agro-pastoral settings cattle are wealth rather than meat, so cattle means
   **milk**. Fresh cow's milk is an established iron-deficiency risk in young
   children: it is iron-poor, it inhibits absorption of the iron in other foods,
   and in infants it can cause occult gut blood loss.

**And we can test it.** Sierra Leone's survey recorded individual bean intake,
household cattle numbers *and* child milk frequency alongside the biomarkers,
while The Gambia recorded cattle and milk. **If the household-level association
is protective, explanations 1 and 2 are right. If it is harmful, explanation 3
is, and the cattle finding is programmatically actionable.**

**Either way, the rule the rest of the talk obeys is the same:**

> **These indicators tell you where to look. Do not read them as a list of
> things to change until we have shown which of the three this is.**

*Source: `p4_admin1_continuous_predictors.csv`, where `dhs_c_fg_legumes` has
meta_z +4.87 (4/4) and `dhs_hh_cattle` +4.40 for child iron (3/3). Positive
means more deficiency, because the concentration scan negates the biomarker.
**Do not rebuild this slide without checking the sign convention**, since the
first draft of the long outline got it backwards. Individual-level variables
verified present: `gw_cBeans`, `gw_hCattleNumb`, `gw_cMilk`, `gw_cMilkTimes`
(Sierra Leone); `gw_cattle`, `gw_cow`, `gw_cMilk` (Gambia).*

> **If the individual-level check gets run before the talk, this slide gets
> better and shorter:** state the contradiction, then state the answer. See the
> offer at the top of this document.

### 13. The one that behaves (20 seconds) *(cuttable)*

**Message:** where children are wasted, women are vitamin A deficient.

- District child wasting tracks women's vitamin A deficiency (4 of 4 countries),
  and rainfall tracks less of it (4 of 4).
- This is interpretable exactly as stated, because both sides are genuine
  district properties, with no aggregation question to answer first.
- Programme reading: the two co-locate and may be targetable together.

---

## Section 6. Transportability: countries with no survey (3 slides, ~2.5 min)

### 14. The central challenge, and where we have got to

**Message:** train on three countries, predict the fourth, having never seen it.

Open with the NCE's own framing, because this is the project's stated core
problem and the audience should hear it as such:

> *"Transferring models to countries without survey data remains the central
> challenge."*

- **How we test it:** remove one country entirely from training, fit on the
  rest, then rank the excluded country's districts. That simulates the target
  use case, a country with no VMD survey of its own.
- **Name the constraint the NCE names.** Only a handful of countries have
  national nutrition survey data at all, and their survey-estimated national
  prevalences differ enormously. That second half is not an aside; it is exactly
  why levels do not transport, below.

Then the result, in the NCE's settled sentence:

> *District and regional rankings transport to a country never used in training,
> with a mean rank correlation of about 0.3 at both tiers, far outside a
> permutation null (95th percentile 0.16).*

- **District level: 0.28, positive in 17 of 22 tests. Region level: 0.30,
  positive in 17 of 22.** Nulls top out at 0.08 and 0.16, so *p* ≈ 0.02.
- **But the level does not transport.** Sierra Leone's absolute prevalence cannot
  be predicted from Ghana's, so the honest product is a **priority list rather
  than a prevalence map.**

*Source: `transport_null_calibration.csv` and `benchmarks_v2_summary.csv`
(`estimand == "country"`); wording per NCE_IMPLICATIONS item 4.*

> **Do not quote "12 of 12" or "0.50–0.56" anywhere**; the NCE is explicit. A
> population-file join silently dropped Sierra Leone from every
> population-weighted run (BUG-01).

#### The Côte d'Ivoire figure: built and validated

**This is the slide's visual, and it should be Côte d'Ivoire**, because CIV is a
real country with a full proxy database and no biomarker survey. It makes the
abstract claim concrete. The figure now exists:
`results/figures/policy_deck/fig10_civ_ranking.png`.

**What it shows.**

- **Left, a rank map.** All 33 CIV districts, shaded by **predicted rank**
  rather than by predicted prevalence, captioned *"which districts to reach
  first."* The worst-ranked districts are Tchologo, Bounkani, Bagoué, Poro and
  Hambol, which are the northern savanna belt: the drier, poorer, more pastoral
  half of the country, and the gradient slide 11 describes.
- **Right, how firmly each district is placed.** The training set was resampled
  400 times and the index refitted each time, so every district has a
  distribution of ranks rather than a single one. The map shows the width of the
  90% range, in places out of 33. The median width is 5 places, the five
  worst-ranked districts stay in the worst third in every refit, and 10 of the
  33 districts have at least an 80% chance of being in the worst third while 21
  have at most 20%. Caption: *"how firmly each district is placed."*

**Be precise about what that uncertainty is.** It is estimation uncertainty: how
much the ranking depends on which training districts we happened to learn from.
It is not a test of whether the model transports to Côte d'Ivoire at all, and no
CIV-internal quantity can be, because there is no CIV ground truth. The external
bound on that is the held-out transport accuracy of 0.37. Say both sentences if
anyone asks; the second one is the honest half.

**Methods available for ranking uncertainty, since you asked.** The conformal
machinery in `R/conformal.R` produces intervals on a *level*, not a rank, and
for the stored CIV output those intervals are degenerate. Two resampling
approaches do work on a ranking, and both are implemented in
`scripts/policy_deck/06_civ_rank_uncertainty.R`: a bootstrap over training
districts stratified within country (400 refits, used for the map), and
leave-one-training-country-out (4 refits), which answers the different question
of whether the ranking depends on *which* countries we learned from rather than
on how many districts they contain. The per-district table carries both.

Two maps, and together they make the section's argument: the ordering is
informative at the extremes, and the level is not available at all. Slide 16 is
the answer to the missing level.

**How it was built, and the two guards it passed.** The index was fitted on the
four surveyed countries using the estimand-C harness from
`scripts/protocol_v2/02b_merge_and_loco.R`, then applied to CIV. Two checks ran
first, because a new country silently joining a vocabulary is exactly how this
pipeline has broken before.

1. **Reproduction.** The same code with each surveyed country held out gives
   mean Spearman **0.3690** over 22 cells, 22 of 22 positive, against the
   published **0.369, 22 positive** in `nested_domain_selection.csv` (arm
   `fixed_cs`). Exact to three decimals, so this is the code path that produced
   the deck's numbers rather than a lookalike.
2. **Intersection cost.** Adding CIV to the column pool costs 7 of the 95
   climate-and-soil columns, because its SoilGrids block covers only 8 of 33
   districts and the 0.70 coverage screen drops it. Transport does not degrade:
   **0.386, still 22 of 22**, on 84 common columns. The iSDA soil block (38
   columns) and all 50 climate columns came through complete, which is why the
   loss is harmless.

**Do not build a "predicted prevalence map plus uncertainty map" pair.** Three
reasons, in order of severity:

1. **The uncertainty map would carry no information.** In the stored production
   predictions, the intervals are split-conformal, constant-width, then clamped
   at zero, so `ci_lo` is **exactly 0 for all 132 district-outcomes** and
   `ci_hi` is `pred_prev` plus a per-outcome constant (0.422 child vitamin A,
   0.230 women's vitamin A, 0.508 child iron, 0.652 women's iron, exact to three
   decimals across all 33 districts). An interval-width choropleth is therefore
   **the prediction map re-coloured**: two maps, one piece of information, and an
   audience member who works it out will not trust the rest of the deck.
2. **A predicted-prevalence map contradicts slide 14.** We say the level does
   not transport. Putting an unanchored prevalence map on the next slide undoes
   that in one image, and it is the image people photograph.
3. **Those stored predictions are not from the model the deck describes.** They
   come from `R/oos_prediction.R`, the production area-level elastic net on GEE
   rasters, whereas slides 14–16 report the protocol_v2 domain index and the
   climate+soil index. **They are different models.** The figure now on the slide
   avoids this by being built from the protocol_v2 index directly.

**Two related resolutions.**

- **Women's vitamin A is resolved and is not a bug.** It predicts 0.8–2.1%
  because women's vitamin A deficiency really is that rare in the training data,
  and because the cell is degenerate in Malawi (98% zero districts) and sparse in
  Ghana and Sierra Leone. See the diagnosis under slide 5. **Do not use women's
  vitamin A for the CIV figure.**
- **Child iron is the right outcome for the figure.** It has the widest
  predicted spread and a well-behaved target in every training country, with no
  zero districts in Gambia or Sierra Leone and 10–19% elsewhere, so the rank map
  is legible and defensible.

**On the vintage question, now settled.** No new Earth Engine acquisition was
needed. Every climate and soil layer CIV requires was already on disk as
exported rasters in `data/Cote_dIvoire_GEE_rasters/` (83 files), in the same
form the four training countries use.

| Block | Shared columns | CIV raster on disk | Training-country vintage |
|:---|---:|:---|:---|
| iSDA soil, 11 properties × mean/stdev × 2 depths | 38 | `Soil{Aluminium,CEC,Calcium,Iron,Magnesium,Nitrogen,Phosphorus,Potassium,Sulfur,TotalCarbon,Zinc}_Cote_dIvoire.tif` | same, static |
| Night-time LST, monthly and annual | 18 | `LST_Night_{2014,2015}_Monthly` | Gambia and Malawi also 2014/2015 |
| TerraClimate | 12–13 | `TerraClimate_Cote_dIvoire_{2017,2018}.tif` | Gambia 2017/2018 |
| Precipitation, annual | 10 | `TRMM_Cote_dIvoire_2010…2019.tif` | same series |
| Aerosol optical depth | 1 | `AerosolOptical_Cote_dIvoire_{2017,2018}.tif` | Gambia 2017/2018 |
| Köppen / AEZ16 | 9 | offline, `data/Koppen_geiger_tif/1991_2020` via script 10 | same |
| SoilGrids | 7 | `data/external_cache/cote_divoire_external_predictors.rds` | same, 8 of 33 districts |

**Why 2018 is the right t0 on the evidence rather than by preference.** The `_t0`
suffix does not mean "a raster from the survey year"; it means the nearest
available vintage. The training countries prove it. Gambia's survey is 2018 and
Malawi's is 2016, yet **both** use LST_Night 2014/2015, because that is what was
exported. TerraClimate and aerosol *are* survey-matched, and CIV's are
2017/2018, identical to Gambia's. At t0 = 2018 every CIV layer therefore sits on
exactly the vintage a training country already uses; at 2016 or 2022 it would
not, and the rasters would have to be re-exported. The decision and the data
agree.

*Scripts: `scripts/policy_deck/04_civ_climate_soil_prediction.R` (fit and
guards) and `05_civ_map.R` (figure). Tables:
`results/tables/policy_deck/civ_climate_soil_ranking.csv`, which carries all six
outcomes for all 33 districts, and `civ_transport_guards.csv`.*

### 15. Two things worth knowing about that ranking

**Message:** the transportable signal is remotely sensed, and it is not simply an
urban–rural map.

- **Climate and soil alone transport better than all 454 indicators**, at 0.37
  district level and 0.46 region level against 0.28 and 0.30 for the full set.
  The two domains available for any country on earth, free, are the ones that
  travel.
- **It is not just urbanicity.** Partialling out night lights, population
  density, built surface and travel time moves the result by less than 0.02, and
  urbanicity on its own does not transport at all (−0.07).
- **Caveat, on the slide:** the climate+soil result was *identified* on these
  four countries, so quoting it as a finding would be selection on the test.
  **We have pre-registered seven quantitative predictions for the next
  countries**, including that this two-domain index will transport district
  rankings at least as well as the full model.

*Source: `climate_soil_admin1.csv`, `urbanicity_conditioning.csv`,
`nested_domain_selection.csv`, and
`docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md`; framing per
NCE_IMPLICATIONS item 6.*

### 16. One national number turns a ranking into a map (20 seconds)

**Message:** the cheapest useful survey we can design.

Use the NCE's sentence verbatim, because it is better than anything I would
write:

> *For a country without a micronutrient survey, a small national biomarker
> sample combined with the transported district ranking gives district-level
> estimates that a conventional survey only matches at roughly eight times the
> sample.*

- Concretely: **one national prevalence, roughly 5% of a full survey's sample,
  50–70 blood draws.** District error is **10 percentage points**, against
  **20** for a district survey of the same size and **12** for a regional one,
  with crossover at about 40% of a full sample.
- **Honest counterpoint, and say it:** the anchored design loses to a district
  survey on *burden captured* at every size, because burden concentrates in
  populous districts, which a survey measures precisely and a ranking does not.

*Source: `anchor_and_rank_summary.csv` and `anchor_and_rank_crossover.csv`;
NCE_IMPLICATIONS item 4b. This is the most directly actionable slide in the deck
for a funder audience, so give it more than 20 seconds if you cut elsewhere.*

---

## Section 7. The interactive dashboards (1 slide plus demo, ~1.5 min)

### 17. Everything above is a click away

**Message:** we built the delivery layer, not just the model.

Screenshot the map explorer and list the rest.

- **Shiny dashboard, 15 tabs.** The four a policy audience will use:
  - **Map explorer**, giving predicted prevalence, uncertainty width, population
    at risk and WHO classification by district, with a detail card on click.
  - **District profiles**, showing every nutrient side by side for one district,
    with intervals.
  - **Decision value**, answering "if I can reach N districts, what do I
    capture."
  - **Plan a survey**, showing where a new survey would be most informative.
    *This is the tab that operationalises slide 16.*
- The rest are for the technical audience: benchmarks, diagnostics, what drives
  the estimate, transportability, resolution and anchoring, methods comparison.
- **Côte d'Ivoire tab**, a full worked prediction for a country with **no
  biomarker survey**. This is the deliverable, demonstrated.
- **Country briefs**, a standing 4-page PDF and HTML per country plus
  district-level CSVs, regenerated whenever the pipeline re-runs
  (`dashboard/report/out/`).
- **Status, in the NCE's words: the dashboard is *ready for field testing*, open
  to invited users.** Say that rather than implying it is a public product. It
  also sets up the ask, which is that **we want countries to test it**, with two
  use-case countries and two regional workshops budgeted for exactly that.

> **The deployment question is settled by the NCE.** It is an invited-user tool
> ready for field testing rather than a public link, so say "available to
> invited users; we are looking for field-test partners" and do not put a URL on
> the slide. A 45-second live demo still beats a screenshot if the app runs
> reliably in the room, but a screenshot is the safe default.

*This slide is also the natural place to name the **January 2026 Accra
meeting**, which is where the partnerships behind the field testing came from,
and this audience may include people who were there.*

---

## Section 8. Conclusions, limitations, next steps (2 slides, ~2 min)

### 18. What is ready, and what is not

**Message:** be precise about the boundary, because the failure mode here is
overclaiming, and the NCE has already been through one round of softening that
this deck should not undo.

**Ready now**

- Ranked district priority lists, in countries with a survey *and* countries
  without one.
- A calibrated map for any country willing to measure one national prevalence.
- Delivered through a dashboard and standing country briefs.

**Not ready, and the reasons matter**

- **Not a prevalence map without an anchor.** Rankings transport; levels do not.
- **It does not reduce the sample a survey needs.** Tested symmetrically: as the
  survey shrinks, the model's district error degrades in lockstep with the
  survey's (+3.36 vs +3.33 pp from full sample to 15%). In the NCE's words, *at
  any given survey size the model improves the ordering of districts and
  modestly reduces district error; it does not reduce the sample a survey
  needs.* **Do not frame savings as coming from a smaller sample.**
- **Not a list of things to fix.** *(Repeat the slide 12 rule.)*
- **Weak where districts are few or deficiency is rare.** Sierra Leone has 14
  districts, Malawi's child vitamin A prevalence is 0.4%, and child zinc shows
  no measurable district-level variation at all.
- **Less headroom than we thought.** Most districts are represented by a single
  survey cluster, so earlier ceiling estimates counted cluster noise as
  geography. Correcting for it, the best achievable correlation averages
  **0.44** rather than 0.54. We are already at that ceiling in 3 of 24
  combinations, and 6 of 24 have a ceiling below 0.30. *(NCE_IMPLICATIONS item 3
  suggests this sentence may cost more to defend than it buys; consider dropping
  it and keeping the single-cluster point on its own.)*
- **We tested ourselves hard and withdrew results in both directions.** An
  internal audit re-examined every headline claim, and several did not survive.
  *One line only. It buys more credibility than it costs.*

### 19. Next steps: this is the extension, in plain language

**REWRITTEN against the NCE's actual activity list.** These are funded, scoped
activities with money behind them rather than aspirations. Say so.

- **More countries, and we know what each one buys.** Every additional training
  country adds roughly **0.05** of transported accuracy (**+0.03–0.04** for the
  two-domain index, from a much higher starting point). That is the concrete
  argument for expansion rather than an appeal to generality.
  - **Ethiopia and Pakistan**, with communication under way with both. They are
    also the countries the **seven pre-registered predictions** are written
    against, including that a climate-and-soil-only index will transport at
    least as well as the full model. These were registered before any of their
    data was seen, precisely because that result was found on the current four.
  - **Tanzania** is in progress on our side as an additional biomarker survey.
- **The Proxy Modeling Alliance.** Other modelling groups take **the same
  curated dataset** and apply their own approaches, comparing strategies,
  identifying the most predictive indicators, and testing transportability
  independently of us. Three groups are identified and two more are in
  discussion, with a six-month engagement from November. *For a policy audience
  the point is not the collaboration; it is that our results will be checked by
  people who did not produce them.*
- **A working group with MIMI and GBD** to compare which predictors each of the
  three approaches relies on, and what each can learn from the others, alongside
  continued work with IHME, WFP and the Nutrition Modeling Consortium. *These
  are the other groups producing sub-national micronutrient numbers, and the
  field does not need three unreconciled maps.*
- **Two country use cases and a cost analysis.** Working with countries that
  have run a national survey, we will show what the model gives for the areas
  the survey did not reach, and cost it against the current approach. *Note for
  the speaker: the honest version of the cost story is the anchor-and-rank
  design from slide 16, not "you can run a smaller survey." See the note at the
  top of this document.*
- **Two regional workshops**, building on the January 2026 Accra meeting, on
  what it would take to trust modelled estimates enough to act on them. *This is
  the closing ask, and it is the right one for this room.*
- **On the data side, household consumption and expenditure surveys (HCES and
  LSMS)** are still the largest unexploited source. Everything else tried
  recently moved nothing: food prices matched to the fieldwork month, night-time
  temperature, a 64-dimension satellite embedding, and engineered
  climate/terrain/soil blocks. **The current vocabulary is saturated, and more
  of the same will not help.**

**Closing question to the room:** *What predictive performance is good enough to
be worth acting on?* We can tell you what the model does. We cannot tell you
where your threshold is, and that is what the workshops are for.

---

## Appendix (have ready, do not present)

- Full outcome and cutoff table (BMGF slide 25).
- The 24-domain × 6-outcome importance matrix, the slide people photograph.
- Per-country transport results, all 22 cells.
- WHO category accuracy: with an anchored transported ranking, **54% exact and
  92% within-one-category** at district level, and 77% exact against the WHO 20%
  threshold (`risk_category_accuracy_summary.csv`).
- The reliability-ceiling analysis (`variance_components_ceiling.csv`).
- The ensemble-versus-index comparison, for anyone who asks why the ensemble is
  not the headline (`estimator_tournament_truth.csv`).
- The Côte d'Ivoire guard table (`civ_transport_guards.csv`), for anyone who
  asks how we know the CIV map is on the training scale.
- Withdrawn results with dates (`docs/findings/CLAIMS_REGISTER.md`).
- **The sign convention**, explicitly: positive means more deficiency, and the
  concentration scan negates the biomarker. Somebody will ask after slide 12.

---

## Numbers retired from the old decks: check each before it disappears

| Old claim | Where | Status now |
|:---|:---|:---|
| "Predicts national deficiency prevalences accurately" | both decks, conclusions | Not re-verified under the current protocol. Not carried forward. |
| "Transportable across countries" (unqualified) | Ghana deck 25 | **Qualified.** Rankings transport (0.28–0.30); levels do not. |
| The 12-learner ensemble as *the model* | both decks, methods | **Retired by the NCE.** Four meta-learners tested; none beats the zero-tuning domain index at these sample sizes. See slide 8. |
| Admin-2 risk-category accuracy 44% / 31% exact | both decks, open questions | **Superseded.** Now 54% exact and 92% within-one, with an anchored transported ranking. |
| B12 AUC 0.689, PR gain 3.09×, Brier skill 0.095 | both decks | Individual-level model, now the *sensitivity* analysis rather than the headline. Recommend dropping, since AUC does not map onto a targeting decision. |
| Spatial ensemble comparison A/B/C/D, Ghana child iron | BMGF 7 and 28 | Superseded by slide 9. |
| "26% more burden than regional averages" | long outline, slide 20 | **Do not use as a headline** (NCE_IMPLICATIONS item 1). Current tables: 22% against 20%; NCE draft: 24% against 19%; per cell the model wins in only 8–11 of 18. Use the ranking claim instead. |
| Regional transport "12 of 12 at 0.50–0.56" | earlier drafts | **Withdrawn (BUG-01).** All 22 cells: 0.30. The NCE says do not quote it anywhere. |
| "Models capture two-thirds of attainable accuracy" | NCE draft, Edit 6 | **Weakened.** Upper bound only, since most districts are a single cluster. Consider omitting. |
| Tanzania listed among the four surveys | BMGF 3 | In progress, not yet in results. Slide 4 lists the four that are. |
| "~100–200 DHS variables" etc. (source table) | BMGF 4 | **Replaced** with actual current counts (slide 6). |

## Open questions for you

1. **Run the individual-level legume and cattle check before the talk?** See the
   top of this document. It would upgrade slide 12 from a puzzle to an answer,
   and for cattle it might produce the one actionable finding in the deck.
2. **The two NCE issues at the top**, namely the corrupted sentence in the
   achievements paragraph and the cost-analysis bullet that proposes to price a
   "reduced sample size" the evidence says does not exist. Both concern the
   request itself rather than the slides.
3. **Audience and venue?** Written for a general nutrition-policy audience. If it
   is a specific funder or one country's nutrition secretariat, slides 10 and 16
   should be rebuilt around that country's numbers.
4. **Slide 5: do we say on the slide that Malawi's vitamin A cells are
   near-empty?** I recommend yes.
5. **Slide 18: keep or drop the headroom bullet?** The NCE note suggests it now
   costs more to defend than it buys.
6. **Dashboard: live demo or screenshots?** Depends on whether there is a
   reachable deployment.
7. **Any national-level (VMNIS/BRINDA) content?** Both old decks had a slide. At
   15 minutes I dropped it, and it is the easiest thing to add back.
