---
title: "A 15-minute talk for the Micronutrient Forum: critical evaluation and a brainstormed outline"
subtitle: "Working document for editing — Andrew Mertens and Sonja Hess, MNF 2026, Accra"
date: "18 September 2026"
---

# How to use this document

This is a **working draft for you to cut into**, not a deck and not a recommendation you have to accept whole. It has two halves.

**Part 1 is a critical evaluation** of the analysis pipeline, the proxy dataset, the results and the 92-slide full talk. It exists because several things the full talk currently *says* about its own numbers do not survive contact with the tables, and those sentences would be the ones a statistician in Accra picks up. Read sections 1.1 and 1.4 if you read nothing else.

**Parts 2 to 8 are the brainstormed 15-minute talk.** At every branch point you get two or three labelled options with the trade-off, and a recommendation. Nothing is presented as settled. Part 8 lists the decisions only you can make; several of them gate figure-building work, so they are worth settling first.

## Contents

- **Part 0. Report the measurable combinations, not the average of everything**
    - The problem with both current framings
    - A screen that uses only survey properties
    - The eight exclusions, each with a reason a programme officer accepts
    - Why this is safer than the 0.64 cut
    - How to use it in the talk
    - One consequence for the rest of this document
- **Part 0B. Reading across all the models: what kind of thing is predictable?**
    - 0B.1 The headline: three findings that hold together
    - 0B.2 Which deficiencies are predictable?
    - 0B.3 Which country did we learn the most from, and which was easiest to learn?
    - 0B.4 Are women or children more predictable?
    - 0B.5 What does each nutrient's model actually reach for?
    - 0B.6 How to put this in the talk
- **Part 1. Critical evaluation of the pipeline, dataset, results and deck**
    - 1.1 Overall verdict
    - 1.2 The analysis pipeline and protocol
    - 1.3 The proxy dataset
    - 1.4 The results
    - 1.5 The national-level (VMNIS) track
    - 1.6 The full talk deck as a communication object
    - 1.7 Findings raised and refuted in verification
- **Part 2. The recommended 15-minute talk**
    - B0. Title options and the one sentence the room should repeat
    - B1. Recommended spine at a glance
    - B2. Slide-by-slide, with options and talking points
    - B3. Alternative spines
- **Part 3. Communicating accuracy to a policy audience**
    - C1. Plain-language accuracy statements
    - C2. Six devices for saying "how good is it" without a correlation
- **Part 4. Visualisation brainstorm**
    - D1. Visuals to keep from the 90-slide deck
    - D2. New visuals to build
    - D3. Visuals to retire, and why
- **Part 5. The national-level (VMNIS) slide: include or appendix**
- **Part 6. Next steps and augmenting nutrition policy**
    - 6.0 The conclusion point to lead with: more surveys make every map better, and we are asking for them
    - F1. Three alternative framings of the closing
    - F2. The where / when / how-to-survey evidence, in plain language, with its honest limits
    - F3. Which proxy data to collect next
    - F4. The ask to the room
    - F5. The fifth survey and pre-registration
- **Part 7. Appendix plan**
- **Part 8. Open questions for you**

---

## The three decisions that gate everything else

1. **Which Ghana cell is the worked example.** The current deck uses child vitamin A, which is the one Ghana cell where the survey's own regional average beats the model on both targets. Child iron is the cell that makes the point. Switching requires rebuilding `fig7_ghana_map.png` and producing a `deploy_ghana_child_iron.csv`.
2. **Whether the headline accuracy is reported over all 24 country-nutrient combinations or over the combinations that are measurable at district level at all.** This is the subject of the next section, and it is the single biggest change to the talk's headline numbers.
3. **Eight slides or twelve.** The Forum's one-slide-per-two-minutes rule gives eight including the title and Sonja's opening. The 12-slide variant in section 2 (B1) is a 19–20 minute talk, not a 15.

## Conventions used throughout

- **"Held out"** means the model never saw that district's, region's or country's blood results while it was being fitted. Every accuracy number in this document is held out.
- **"The regional average"** means the survey's own Admin-1 mean *computed without the district itself* (jackknifed, respondent-weighted). It is the number a ministry would use today for an unsurveyed district, and it is the comparator that matters.
- **Ranking accuracy** is Spearman correlation: 0 is chance, 1 is perfect. Where possible it is translated into a share of district pairs ordered correctly, where 50% is a coin toss.
- Numbers carry their source table. Nothing on the withdrawn list appears — no "12 of 12 regions at 0.50–0.56", no "26% more burden than regional averages", no "two-thirds of attainable accuracy", no "the SuperLearner chose the index".

---

# Part 0. Report the measurable combinations, not the average of everything

*This section responds to a specific instruction: highlight the best-working models rather than an average dragged down by outcomes nobody could predict. It is new analysis, run on 18 September 2026, and it is the recommended headline framing for the talk.*

## The problem with both current framings

The full talk currently offers two numbers and neither is quite right for a policy room.

- **The 24-cell average (0.40 in-country, 0.30 transported).** Honest, but it averages over four combinations where the deficiency is so rare that no method could rank districts (Malawi child and women's vitamin A at 0.4% and 0.2% national prevalence) and four where the survey's own district estimates carry no district-level signal at all (both Malawi zinc cells, Sierra Leone child vitamin A and women's iron). Reporting those in the mean understates the method wherever it can actually be used.
- **The six "strong cells" (0.64).** This is the number you want to be true, but the six were selected *on the same held-out scores they are reported with*. Three of the six are The Gambia. A reviewer will say so, and they will be right.

## A screen that uses only survey properties

Keep a country-nutrient combination if **both** of these hold, each computable from the survey alone before any model is fitted:

1. **National prevalence ≥ 2%** — WHO's own threshold for a deficiency being a public-health problem.
2. **The survey's reliability ceiling ≥ 0.30** — i.e. the survey's own district estimates carry district-level signal. This comes from the variance-component model (respondents in clusters in districts) and is a property of the survey design, not of our model.

Neither condition looks at how well the model did. Sixteen of the 24 combinations pass; 14 of those have in-fill scores.

| Set | Cells | In-country | Regional avg | Better in | Pairs right | Transported |
|---|---|---|---|---|---|---|
| **Measurable** | 16 (14 scored) | **0.49** | 0.39 | 11 of 14 | **64%** vs 59% | **0.35** |
| Excluded | 8 | 0.26 | 0.17 | — | 45% | 0.27 |
| All combinations (current headline) | 24 (18 scored) | 0.40 | 0.31 | 13 of 18 | 60% vs 56% | 0.30 |
| Post-hoc cut at ρ ≥ 0.5 (do not use) | 6 | 0.64 | 0.48 | 6 of 6 | 70% | 0.60 |

*Sources: `results/tables/protocol_v2/targets_v2.csv` (prevalence), `variance_components_ceiling.csv` (ceiling), `benchmarks_v2_cells.csv` (scores), `nce_targeting_metrics.csv` (pairs). Recomputed 18 September 2026.*

## The eight exclusions, each with a reason a programme officer accepts

**Too rare to rank at district level:** Malawi child vitamin A (0.4% national prevalence), Malawi women's vitamin A (0.2%), Sierra Leone women's B12 (0.6%), Sierra Leone women's vitamin A (1.0%). At these levels most districts record zero deficient respondents; there is no order in the survey for any model to recover.

**The survey itself finds no district-level signal:** both Malawi zinc combinations, Sierra Leone child vitamin A, Sierra Leone women's iron. The variance model estimates zero district-level geography — what looked like between-district variation is cluster and respondent noise. For zinc there is a known mechanism: serum zinc is unadjusted, fasting rules vary between clusters, and an afternoon blood draw lowers it by 3.5–7.3%.

## Why this is safer than the 0.64 cut

1. **It is not selected on the outcome.** The screen uses prevalence and the survey's own resolution. The 0.64 cut uses the model's score.
2. **It is insensitive to where you put the thresholds.** Every combination tried — prevalence at 2%, 5% or 10%, ceiling at 0, 0.30 or 0.40 — lands between 0.46 and 0.50. That insensitivity is itself the argument: *any* reasonable definition of "measurable" gives the same answer.
3. **It keeps a loss in.** Ghana child vitamin A passes the screen and the model loses there (0.28 against 0.41 for the regional average). Retaining a visible loss is what makes 0.49 credible.
4. **It answers a question the room finds interesting on its own.** "Which deficiencies can be mapped at district level at all, and which cannot?" is a real finding about micronutrient surveys, not an apology for a model. Answering it *before* showing accuracy converts a caveat into a contribution.

## How to use it in the talk

**Recommended:** make the screen a visible step on slide 5 rather than a footnote. Two sentences:

> "Before we show you accuracy, a question worth asking on its own: which deficiencies can a survey even see at district level? Eight of our twenty-four country-nutrient combinations fail that test before any model is involved — either the deficiency is too rare to rank, or the survey's own district numbers are noise. We report those separately, and we can tell you in advance which combinations will fail."

Then the headline: **0.49 against 0.39 for the regional average, 64% of district pairs against 59%.**

**Caveat to state aloud, and to act on:** the 2% and 0.30 thresholds were chosen after seeing these data. For the fifth survey the screen belongs in `docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md` alongside the two transport candidates. Written down in advance, 0.49 becomes a prediction rather than a description — which is exactly the move the project already makes for climate-and-soil.

**Alternative A — report both.** Give the 24-combination mean and the measurable-combination mean side by side, with the screen shown as a filter. *Trade-off:* maximally transparent, costs 20 seconds and splits the room's attention across two numbers when it will only remember one.

**Alternative B — keep the 24-combination mean as the headline and put the screen in the appendix.** *Trade-off:* safest possible, and it is what the current deck does; but it systematically understates the method in every combination where it can be used, and the room never learns that some deficiencies are unmappable, which is a finding.

## One consequence for the rest of this document

Where numbers below are quoted over 18 or 24 combinations (0.40, 60% of pairs, 0.30 transported), the measurable-combination equivalents are 0.49, 64% and 0.35. Both are correct; they answer different questions. The slide-by-slide scripts in Part 2 are written with the 24-combination numbers because that is what the current tables and figures carry — swap them if you adopt the screen, and note that the figures listed in Part 4 will need regenerating on the screened subset.

---


# Part 0B. Reading across all the models: what kind of thing is predictable?

*New analysis, 18 September 2026, run across every country, outcome and population, and across the back-projected variable importances. It answers four questions the talk currently has no slide for: which deficiencies are predictable, which country taught us most, whether women or children are easier, and what each nutrient's model actually reaches for. This is the "story beyond the statistics" layer, and it is the strongest candidate for the slide the room remembers.*

## 0B.1 The headline: three findings that hold together

1. **Predictability tracks how directly a nutrient's status is written on the landscape and on other measured surfaces — not how common the deficiency is.** B12 and vitamin A are the best-predicted nutrients; iron is middling; folate ranks at home but does not travel; zinc is unpredictable everywhere.
2. **Country differences are survey-design differences, not data differences.** The ordering of how well each country is predicted follows the ordering of what its survey could resolve, and in every country the model reaches roughly three-quarters to nine-tenths of what was attainable.
3. **Each nutrient reaches for a different, biologically legible second signal on top of a shared environmental base** — and the nutrient with the clearest dietary marker (B12 → flesh-food consumption) is the best-predicted one. That is the sentence that turns a table of correlations into a finding.

## 0B.2 Which deficiencies are predictable?

Measurable combinations only (the Part 0 screen), in-country ranking on the biomarker level, with the survey's own reliability ceiling beside it:

| Nutrient | Cells | Prevalence | Ceiling | **In-country** | Regional avg | **Transported** | Share attainable |
|---|---|---|---|---|---|---|---|
| **B12** (women) | 2 | 10% | 0.65 | **0.63** | 0.53 | **0.52** | ~0.97 |
| **Vitamin A** | 4 | 13% | 0.53 | **0.53** | 0.44 | **0.51** | ~1.0 |
| **Iron** | 7 | 27% | 0.61 | 0.44 | 0.31 | 0.33 | ~0.72 |
| **Folate** (women) | 3 | 52% | 0.64 | 0.40 | 0.38 | **0.07** | ~0.61 |
| **Zinc** (Malawi) | 2 | 60% | not estimable | **−0.07** | −0.07 | — | — |

*Sources: `benchmarks_v2_cells.csv`, `variance_components_ceiling.csv`, `targets_v2.csv`. "Share of attainable" is the achieved score over the survey's ceiling and is indicative only — the ceiling estimates are fragile (26 of 44 variance fits are singular) and this ratio should not be quoted as a headline percentage.*

**Three things to say about this table.**

**Prevalence does not predict predictability.** Folate and zinc are the two most common deficiencies in the set (52% and 60%) and are the two the model handles worst. Vitamin A at 13% and B12 at 10% are the best. A programme officer's intuition — "surely the common ones are easier" — is wrong, and saying so is a good 15 seconds.

**Folate is the interesting failure, and it is a measurement failure, not a model failure.** Folate ranks districts inside a country at 0.40, essentially as well as iron, and then collapses to 0.07 when the model crosses a border. The reason is documented: Malawi's folate was measured with a microbiologic assay and the others with a Roche immunoassay, which reads systematically lower, so "folate deficiency" is not the same quantity in the training and test countries. This is the clearest illustration in the whole project of the point the transport slide makes abstractly — *what does not cross a border is the measurement, not the geography.* It is also an argument the Forum audience can act on, because assay harmonisation is a live agenda item there.

**Zinc is the honest floor.** Serum zinc is unadjusted, fasting rules vary between clusters, an afternoon draw lowers it by 3.5–7.3%, and the survey's own variance model finds no district-level geography at all. Both zinc combinations score below zero. The right framing is not "our model failed on zinc" but "the survey's own district zinc numbers are noise, and we can show you that before anyone fits a model."

## 0B.3 Which country did we learn the most from, and which was easiest to learn?

These are two different questions and they have different answers.

| Country | Surveyed units | Clusters/district | Ceiling | **In-country** | **Held out** | **Teaches others** |
|---|---|---|---|---|---|---|
| **The Gambia** | 30 districts | 2.33 | 0.70 | **0.62** | **0.60** | 0.16 |
| **Ghana** | 75 of 260 districts | 1.20 | 0.54 | 0.44 | 0.35 | **0.29** |
| **Malawi** | 87 Traditional Authorities | 1.18 | 0.57 | 0.43 | 0.20 | 0.18 |
| **Sierra Leone** | 14 districts | 4.29 | 0.62 | not scoreable | 0.16 | — |

*"Teaches others" is the mean transported ranking achieved in the other countries when this country is the sole training country (`training_country_curve.csv`, level target). "Predicted (held out entirely)" is leave-one-country-out (`benchmarks_v2_cells.csv`).*

**The Gambia is the most learnable country; Ghana is the best teacher.** Trained on Ghana alone, the model ranks The Gambia's districts at 0.61 — as well as when it has The Gambia's own survey. Trained on The Gambia alone it reaches Ghana at 0.30 and Malawi at 0.09.

**Why The Gambia is easiest: its survey could see districts.** 2.33 clusters per district against 1.2 in Ghana and Malawi, and the highest reliability ceiling (0.70). Note what this is *not*: The Gambia has the fewest districts and the narrowest environmental range in the set. Being easy to predict is a property of the survey, not of the country.

**Why Ghana teaches best: it spans the widest gradient.** Ghana's 260 districts run from Sahelian north to forested south — the widest spread of greenness (NDVI SD 0.16 against 0.04–0.09 elsewhere) and plant productivity (GPP SD 1.24 against 0.36–0.91) in the panel. A training country that contains the whole gradient teaches the model the whole gradient.

**A hypothesis I tested and had to drop.** The obvious explanation for Sierra Leone's failure is extrapolation — that it sits outside the environmental range the model was trained on. Measured across all 105 climate and soil columns, that is false: **91% of Sierra Leone's district values fall inside the other three countries' range, the highest share of any country**, against 80% for Malawi. Extrapolation does not explain the country differences and should not be offered as an explanation. What remains for Sierra Leone is more mundane and better evidenced: 14 districts is too few to run the district test at all (a Spearman on 14 units has a standard error near 0.28), four of its six combinations fail the measurability screen, and its child file as supplied contains only the survey's anaemic children (532 of 654 assayed).

**What this means for a fifth country, stated so it can be wrong:** a country will be *easy to predict* if its survey puts two or more clusters in a district; a country will *teach* if it has many districts spanning a wide environmental gradient. With three donor countries these are hypotheses, not findings — but both are pre-registerable for Ethiopia and Pakistan, and both are checkable before a single biomarker is measured.

## 0B.4 Are women or children more predictable?

**No consistent difference, and the honest answer is more interesting than a difference would be.** Paired within country and nutrient, the model reaches 0.39 for women's outcomes and 0.34 for children's, with women ahead in 5 of 7 pairs — well inside noise at this sample size. Across borders the two are identical (0.30 each, women ahead in 3 of 8).

The pattern is not population but *combination*: Gambia women's iron 0.68 against child iron 0.41, but Ghana child iron 0.53 against women's iron 0.40 — the ordering reverses between countries. Whatever drives predictability, it is not whose blood was drawn.

One asymmetry is worth a clause: women's vitamin A is far rarer than children's (2.6% against 18.6% in The Gambia; 3.1% against 27.6% in Ghana) and is nonetheless predicted as well or better. That cuts against the "rare outcomes cannot be ranked" rule at the margin — the rule bites below about 2%, not below 5%.

## 0B.5 What does each nutrient's model actually reach for?

Every model rests on the same environmental base — climate, soil and the satellite embedding together take 46% of the index weight, in every outcome. But those are also the three biggest blocks, and a domain's share is partly just its size (share correlates 0.79 with the number of components a domain contributes). **Dividing by block size shows what each nutrient reaches for that the others do not**, and the result is unexpectedly legible:

| Nutrient | Distinctive signal (weight per component) | Leading individual layers (+ = more deficiency) |
|---|---|---|
| **Women's B12** | **Infant and young child feeding, 2.5% per component** — the highest of any domain for any outcome | **− flesh-food consumption**, + exclusive breastfeeding, + diarrhoea, + livestock units per person |
| **Women's iron** | **Anaemia surfaces, 8.0% on a single component** | + modelled anaemia, + moderate anaemia, + ruminant share, − overweight, − pigs, − malaria mortality |
| **Women's folate** | **Anaemia surfaces, 8.6% on a single component**; fortification 3.6% | + modelled anaemia, − oil crops, − iodised salt, − tetanus-vaccine coverage |
| **Child vitamin A** | Child anthropometry and SES, ~2.1% per component | **+ child wasting**, + grassland cover, + under-5 population share, + dependency ratio, − plant productivity |
| **Women's vitamin A** | Child anthropometry, SES, household assets | + malaria prophylaxis in pregnancy, **+ child wasting**, + underweight, + livestock units, + own-production food share, − wealth |
| **Child iron** | Livestock and agriculture, ~1.6–2.2% per component | + cereal-dominant cropping, **+ ruminant share**, + cattle density, − root crops, − overweight, − relative staple price |

*Source: `index_importance_domains.csv` and `index_importance_top.csv`, pooled four-country fits, biomarker level. Every layer named holds its sign in all four leave-one-country-out fits.*

**The story this tells, in three sentences the talk can use.**

**Each nutrient's model reaches for something a nutritionist would recognise.** B12 — the nutrient available almost only from animal-source food — leans hardest on flesh-food consumption, and districts where more infants eat meat or fish have better B12 status. Iron and folate lean on the modelled anaemia surfaces, which is their clinical consequence. Vitamin A leans on child wasting and underweight, the other undernutrition that co-locates with it. Child iron leans on the agro-pastoral gradient: cereal-dominant cropping, ruminant share and cattle density mark worse status, root crops and market-integrated food prices mark better.

**And the nutrient with the clearest marker is the best predicted.** B12 has the most specific dietary signal in the dataset and is the best-predicted outcome (0.63 in-country, 0.52 transported, essentially at the survey's ceiling). That is the cleanest available answer to "why does this work at all": not because satellites see deficiency, but because the food system, the disease environment and the landscape that produce a deficiency are all visible from space and from other people's surveys.

**Three cautions that must travel with this slide.**

- **These are district-level associations, never causes and never things to change.** The sharpest illustration is Malawi's B12 model, where the modelled anaemia and malaria surfaces carry *negative* weight — districts with more malaria have *better* B12 — because they mark the fish-eating lakeshore. The model finds a place-marker and cannot tell it from a mechanism.
- **Two of the strongest replicated signals point opposite to household-level nutrition advice.** Cereal-dominant cropping and ruminant/cattle density mark *more* child iron deficiency at district level, while at the household level owning cattle and eating legumes are protective. The district-level reading is "this is a subsistence agro-pastoral food system", not "discourage cattle". Say the rule aloud: *these indicators tell you where to look, not what to change.*
- **The anaemia surfaces are the one dependency worth disclosing.** They carry 8% of the iron and folate models on a single component, and they are themselves modelled from earlier DHS haemoglobin. They are admissible as predictors under the project's leakage rule, but "none of this needs a blood sample" is not literally true of them — say "no *new* blood sample", and note that dropping the anaemia domain entirely costs 0.02 of transported accuracy.

## 0B.6 How to put this in the talk

**Recommended — one slide, four rows, as the answer to "what kind of thing is predictable?"** This is a better fifth slide than a second accuracy slide, because it is the only place the talk says something about *micronutrient epidemiology* rather than about a model. Suggested structure:

> **Headline:** *What we can map depends on the nutrient, and on what the survey could see*
>
> Row 1 — **B12 and vitamin A map best; folate travels worst; zinc not at all.** And it is not about how common the deficiency is: folate and zinc are the commonest and the least predictable.
> Row 2 — **Each nutrient's model reaches for something recognisable:** B12 for flesh-food consumption, iron and folate for the anaemia surfaces, vitamin A for child wasting, child iron for the agro-pastoral gradient.
> Row 3 — **Folate's failure across borders is an assay difference, not a model failure** — a harmonisation problem this room can fix.
> Row 4 — **Country differences are survey-design differences:** where a survey put two clusters in a district we reach 0.62; where it put one, 0.43.

**Alternative A — fold rows 1 and 3 into the accuracy slide and give the whole slide to the variable-importance story** (row 2 expanded into five plain-language maps). *Trade-off:* the drivers are what people ask about and the maps are beautiful; but it spends two minutes on associations the talk then has to disown as non-causal, which is a strange use of the budget.

**Alternative B — lead the talk with this slide instead of with the map trio.** *Trade-off:* it is the most scientifically interesting content in the deck and it earns the room's attention honestly; but it delays the deliverable past minute five, and a policy audience came to see a map.

**What to build:** a single figure with nutrients on the y-axis, two dots each (in-country and transported) on a 0–0.7 ranking axis, the survey's ceiling drawn as an open bar behind each, and the distinctive driver printed as text at the right-hand margin. Everything needed is in `benchmarks_v2_cells.csv`, `variance_components_ceiling.csv` and `index_importance_domains.csv`; no new computation is required.

---


# Part 1. Critical evaluation of the pipeline, dataset, results and deck

Sources: the six phase-1 reports (`scratchpad/phase1_{pipeline,dataset,results,deck,national,survey_design}.md`), the adversarial verification round (survived / refuted findings), and `scratchpad/accuracy_statements.md`. Every number traces to a table under `C:/Users/andre/OneDrive/Documents/mn-prediction/results/tables/`. Sign convention: positive = more deficiency. No withdrawn claim is used.

---

## 1.1 Overall verdict

What is solid is the evaluation machinery and the dataset. Folds are cut at or above the district, in-fill is replicated over 10 draws, region and country hold-outs are exhaustive, every covariate-free comparator is information-matched and jackknifed (`arm_region_mean_jk_v2()` never sees the held-out district), the cross-country claim carries a country-block permutation null (95th percentile of the mean rho = 0.079 at Admin-2), and there is a written claims register with withdrawn statuses. The dataset (575 columns, 554 districts, 28 sources, 29 domains) is the best-documented part of the project and is under-sold: one metadata record per column with source, tier, alignment rule, year used and worst-case gap; 23 exclusion rules each with an evidence file; leakage defined by survey instance rather than analyte name. Exactly one product survives all of this: **a district ordering**. In-country it scores Spearman 0.399 on the biomarker level and 0.282 on prevalence against 0.314 / 0.220 for the survey's own jackknifed regional averages (better in 13 of 18 scored cells), and orders 60% of district pairs as the survey does against 56% for the regional average and 50% for a coin. Across borders it scores 0.302 (18 of 22 cells positive, null 0.079), which translates into a worst-third call that is right 52% of the time against 33% by chance.

What is fragile is nearly everything the deck says *about* those numbers. Six verified defects, in order of Q&A damage: (1) the prevalence-error slide plots the **un-calibrated** ranker, whose held-out district prevalence error (14.25 pp, population-weighted) is *worse* than handing every district one national number (11.49 pp) in 13 of 18 cells, while the speaker note says the opposite; (2) the transport figures the closing slide leads with (0.378 climate+soil, 0.414 five-domain) come from domain sets chosen on the same 22 cells they are scored on — the only honest nested selection gives 0.312, indistinguishable from the full index's 0.302; (3) the anchor-and-rank slide is framed as substitution ("a district survey needs 25% of full size to match") when the same table shows a regional survey at 10–15% does the same or better, and the anchored design is at chance on burden; (4) "zero-tuning" and "the SuperLearner chose the index" overstate independence — the discrete SuperLearner picked the index in 3.1% of 1,304 honest fits under squared error and 18.1% under rank loss; (5) the reliability-ceiling defence rests on 26 of 44 singular variance-component fits, six exact zeros and a "59% of attainable" ratio whose numerator and denominator come from different cell sets; (6) burden targeting is at chance and the deck still lists "that fifth holds a fifth of the deficient people" among the gains.

What the talk **must** claim: the model produces an order, tested against the best a statistician could do with the survey alone, and it beats that comparator on the order in most cells and in a country never surveyed. What it **must not** claim: a district prevalence from the un-shrunk index; that the model finds burden (it is at chance, and so are the regional averages); that a survey can be made smaller; that 0.38 is a validated transport number; that the SuperLearner selected the winner.

---

## 1.2 The analysis pipeline and protocol

### Strengths worth saying aloud

Blocked, replicated folds with a fold-draw audit (in-fill draw SD 0.028, median over cells, so draw luck is no longer the story). Information-matched comparators, including a covariate-free spatial smoother (`arm_spatial_v2`, 0.390 in-fill) correctly barred from transport. A permutation null that preserves within-country outcome correlation (500 reps, same permutation across outcomes; p < 0.001). PC orientation learned from training countries only under LOCO (`sign_rows = tr`), after per-country orientation was found to give −0.06 transport. A pre-registration document with seven quantitative predictions dated 3 Sept 2026, and an honest nested-selection arm that labels the post-hoc domain sets as such. Rank-interval calibration reported as stability rather than truth — rarer than it should be in this literature.

### Verified weaknesses, ranked

| # | Weakness | Evidence | Consequence for the talk |
|---|---|---|---|
| 1 | The prevalence slide plots the un-calibrated ranker and the notes have the sign backwards | `benchmarks_v2_summary.csv` infill/prev: index 14.25 pp vs `null_train_mean` 11.49, `region_mean_jk` 11.77, `spatial` 10.93 (18 scored cells; Sierra Leone has no in-fill MAE). Worse than the national number in 13 of 18, worse than the regional average in 15 of 18. Cause: `arm_domain_index_v2()` rescales with rho = 1 (`R/protocol_v2.R` 471–535), over-dispersing by 1/rho | Every level or planning figure must come from the calibrated arm or be dropped. **Editorial correction (verified 18 Sept against the committed tables; the deck chunk was fixed the same day):** `domain_index_cal` is already present in both `benchmarks_v2_summary.csv` and `benchmarks_v2_cells.csv`, at **10.74 pp** mean wMAE over the 18 scored cells — better than the national number (11.49), the regional average (11.77) and the un-calibrated index (14.25), and better than the national number in 9 of 18 cells. **No rerun is needed**; the deck simply does not plot that arm — the chunk at `MN-proxy-full-talk-2026-09.qmd:588` filters `arm %in% c("domain_index","region_mean_jk")` and, although its `case_when` names "One national number", never draws it. So this is a one-line plotting fix, not an analysis task. (For reference, the same arm also appears in the 18 Sept shards `benchmarks_v2_raw_{gambia,ghana,malawi,sierraleone}.csv` at **10.69 pp** wMAE (unweighted 11.43), beating the national number (11.49) and the regional average (11.77) and matching the smoother (10.93) — better than national in 10 of 18, than regional in 15 of 18; child iron 10.8 / 12.6 / 11.8 vs national 12.8 / 14.7 / 11.8.) Cost: pooled Spearman falls (level 0.40 → 0.35, prev 0.29 → 0.21), so **rankings stay with `domain_index`; levels, MAE and bias come from `domain_index_cal`** |
| 2 | Transport domain sets selected on the test cells | `nested_domain_selection.csv`: honest inner-LOCO selection 0.312 level / 0.174 prev vs full index 0.302 / 0.223; fixed climate+soil 0.378; five-domain 0.414. The nested arm recovers climate+soil together in only 4 of 44 selections. The post-hoc-minus-honest gap (+0.07 level, +0.10 prev) is the entire margin over the full index | Quote **0.30** (full vocabulary, 18 of 22 positive, null 0.08) as the transport number; call 0.38 / 0.41 two simpler candidates found on these four countries and pre-registered for the fifth. Note the Côte d'Ivoire map, its rank-interval maps and the anchor design are computed **only** under the climate-soil set, and that `civ_candidates_agreement.csv` (0.88–0.97) compares the two post-hoc candidates with each other, not the full index with climate+soil — a CIV robustness statement against the selection problem needs a full-index CIV ranking that does not exist in `results/tables/policy_deck/` |
| 3 | Anchor-and-rank framed as substitution | `anchor_and_rank_summary.csv` (prev, climate_soil, median of 22 cell means): A1 national anchor + ranking 10.05 pp at f = 0.05 and 8.47 at f = 1; district survey B 9.47 at f = 0.25; **regional survey C 9.88 at f = 0.10 and 9.07 at f = 0.15**, and below A1 at equal budget from f = 0.25 (8.27 vs 8.63) to f = 1 (6.49 vs 8.47). A1's Spearman (0.252) and burden capture (0.217; chance 0.20, ceiling 0.411, B 0.280 even at f = 0.05) are f-invariant by construction. One noiseless national number for every district gives a median **8.62 pp** on the same 22 cells, better than A1-at-5% in 18 of 22 cells. Design B shrinks the existing survey, so at f = 0.05 it has ~1 effective respondent per district | Do not say "replace" or "matches a survey four to five times its size". Defensible: with only a national estimate the model supplies an order (Spearman ~0.25, at chance on burden) at the national number's level; with a regional survey at 15–25% or more of full size the survey wins on level and the model's contribution is the within-region order (A2 0.350 vs C 0.274 at f = 0.25; 0.516 vs 0.488 at full size) |
| 4 | "Zero-tuning" and "the SuperLearner chose the index" | Over 1,304 honest SL fits (SL-06, `sl_selection_sl.csv`) the squared-error discrete SL selected the index in **3.1%** (4.9 in-fill / 1.2 region / 9.1 transport), the rank-loss discrete SL in 18.1%; NNLS weight on the index 0.08–0.13. The `pcvar` representation was adopted because LOCO on these cells rose 0.151 → 0.255; the soil exclusions were reverted because climate+soil fell 0.369 → 0.341 | "No parameter is tuned inside a fit, so a held-out district cannot be overfit; the design choices were made looking at four countries, which is why the fifth survey is pre-registered. We compared the index against 13 alternatives and four ensemble rules on identical folds; nothing beat it by more than 0.03 on any estimand, nothing beat it on prevalence within a country, and the ensembles scored below it because their inner cross-validation runs on 30–70 districts." Drop "the SuperLearner chose it". Honest nuance: in The Gambia (30 districts) the index beats every tuned learner in 8 of 8 in-fill cells by 0.06–0.16; in Ghana (75) ridge, ranger and the spatial GAM tie or edge it (Ghana child iron level: ridge 0.578–0.587 vs index 0.497–0.500); in Malawi (87) the index is ahead again. "Simplest wins" is a sample-size finding, not a monotone law |
| 5 | Reliability ceiling: sound construction, fragile estimate, defensive wording | `variance_components_ceiling.csv`: **26 of 44** Admin-2 fits singular; `ceiling_vc` exactly 0 in six rows; per-cell headroom ratio runs −0.05 to 1.07 with two cells at 0.99–1.00; no ceiling carries an interval. The deck's 59% divides an 18-cell prevalence numerator by a 21-cell denominator (16 cells overlap). On the level target the SAY actually quotes (0.40) the matched share is **0.66**, not 0.59. TWO_READINGS g's "85–90% of a 0.62–0.78 ceiling" is wrong on both ends (ceilings 0.62–0.87, ratios 0.80–1.07) and those six cells were selected on the same achieved score that forms the numerator | Show the ceiling as a **by-country range** (prevalence: Gambia 0.60–0.76, Ghana 0.34–0.75, Malawi 0.45–0.69 non-zero, Sierra Leone 0.62 in two cells and 0 in three) with the single-cluster share beside it (Ghana 83%, Malawi 85%, Gambia 57%, Sierra Leone 0%). Do **not** quote a share of attainable (59% or 66%) — it revives the withdrawn "two-thirds" framing. Close with the design result: two to three clusters per district raises the bar every model is judged against |
| 6 | Burden targeting is at chance, and prevalence is the wrong ranking for a burden decision | `nce_targeting_summary.csv`, 18 in-fill cells: index captures 0.219 of deficient people, regional average 0.203, smoother 0.198, one national number 0.136, oracle 0.482; lift > 1 in 11 of 18, permutation p ≈ 0.5. Under transport 0.209 (10 of 22). Dropping the two Malawi vitamin A cells below 0.3% prevalence gives 0.246 vs 0.206 vs oracle 0.417 (p ≈ 0.10). Structural cause: burden sits in populous districts (Gambia: Kombo Saint Mary 26% of children, Kanifing 16%) while high prevalence is small and rural, so even the true-prevalence fifth captures 0.13–0.49 | Say plainly that on *cases* the model and the regional averages are both at chance and only population data moves that number. `deck_prose.md:364` still lists "that fifth holds a fifth of the deficient people" among gains labelled "real, modest, and better than the number a ministry has today"; fix it. A predicted-cases arm exists in no script; if built it must be scored against a **population-only** ranking, not the prevalence oracle |
| 7 | Cell-count arithmetic | The in-fill mean is over **18** cells, not 24 — Sierra Leone's 14 districts cannot be 5-folded (a ≥ 12-training-district floor in `12_nce_targeting_metrics.R`), and `nce_targeting_summary.csv`'s `cells=24` is an na.rm artefact. `deck_prose.md:191` already says eighteen; line 366 says 24 | Fix every "24 combinations" label attached to an in-fill mean |
| 8 | Raised in phase 1, not put through the verification round | (a) In-country PCs are learned on all surveyed districts including held-out ones (`02_run_benchmarks_v2.R` line 92, `sign_rows = NULL`); no outcome information leaks and transductive PCA is deployment-realistic, but the deck's sentence "the rotations use training districts only" is false for the in-fill and region estimands. (b) Index component weights are implicitly proportional to component variance: domain importance share correlates 0.79 (Spearman) with the number of PC axes a domain contributes, and the three largest blocks take 46% of the index | (a) reword to "the predictor summaries use no outcome information; the weights are learned on the training districts only". (b) describe importance as "where the index's weight sits", not "what drives deficiency", and lead with the drop-one ablation (any single domain costs ≤ 0.025 of 0.30) |

### Reviewer questions, one-sentence honest answers

1. *"Your district prevalence is worse than the national number."* — "Yes for the un-shrunk ranking (14.3 vs 11.5 pp); shrunk to its out-of-sample correlation it is 10.7 and beats both the national number and the survey's regional breakdown, which at district level is itself no better than the national number."
2. *"0.38 was picked on the 22 cells you score it on."* — "Correct; the honest nested number is 0.31, the same as the full index at 0.30, and climate-and-soil is a pre-registered candidate for the fifth survey, not a validated result."
3. *"Does the model let a country run a smaller survey?"* — "No; fitted symmetrically on a reduced survey it degrades in step with the survey (+3.36 vs +3.33 pp from full sample to 15%), and a regional survey at 15% beats the national-anchor design on level; what the model adds at any sample size is the ordering."
4. *"A fifth of districts holding a fifth of the burden is chance."* — "It is, for us and for the regional averages, because burden sits in populous districts; the honest targeting statements are 32% deficiency in the flagged fifth against 24% nationally, and 60% of district pairs ordered as the survey does."
5. *"Isn't 55% of pairs across borders a coin toss?"* — "Nearly, at the pair level; the usable statement is that a worst-third call in a country never surveyed is right 52% of the time against 33% by chance, with a permutation null of 0.08 against our 0.30."
6. *"You chose the estimator, the representation and the domain set on these data."* — "No parameter is tuned inside a fit; the design choices were made on four countries, which is exactly why the fifth survey is pre-registered."
7. *"The ceiling is an excuse."* — "It is a measurement of the survey: 83% of Ghana's and 85% of Malawi's surveyed districts hold one cluster, and the fix is two or three clusters per district in the next survey."
8. *"How much does geography alone get you?"* — "Inside a surveyed country almost everything: a neighbour smoother scores 0.39 to our 0.40; the covariates earn their keep in unsurveyed regions and countries, where no smoother can be fitted."
9. *"Aren't the IHME anaemia surfaces the outcome by another name?"* — "They are modelled from earlier DHS haemoglobin, not from our surveys; dropping all 29 IHME columns changes transport by +0.010 on the level and −0.007 on prevalence."
10. *"Do the map's intervals contain the truth?"* — "No; they cover the held-out survey's rank 38% of the time, so they are stability under retraining; an honest 90% interval is 55% of the list."

### Plain-language methods sentences

1. **Honest testing.** "Every number comes from predicting districts, regions or whole countries the model was not allowed to see, and we compare it against the best a statistician could do with the survey alone — the region's average computed without that district."
2. **Order, not level.** "The model tells you which districts are worse than which. It cannot tell you the prevalence unless the country measures one national number."
3. **Same-size surveys, better placed.** "The model does not shrink the survey you need; it tells the survey where to go, and it fills in the districts a survey cannot afford to visit."
4. **The bar is the survey, and the survey is thin.** "Most districts here hold one cluster of about a dozen people, so even a perfect model would agree only about 0.5 to 0.7 with the survey's own district numbers."
5. **Simple by necessity, pre-registered by design.** "With 14 to 87 districts per country nothing fancier than a weighted sum of public layers beats a weighted sum — and because we chose that sum on four countries, we have written down in advance what it should score on the fifth."

---

## 1.3 The proxy dataset

### Strengths

One metadata record per column (13 fields): domain, source, tier, subnational flag, modelled-surface flag, alignment rule, year used, worst gap. 53 alignment rules; 23 exclusion rules each with an evidence file; coded p-code joins with a build audit that fails on domain-key collision or fan-out. The leakage rule (LK-02: same *survey instance*, not analyte name) is the right rule and is implemented where it can be — 70 GMNS-hosting MICS clusters dropped in The Gambia, 105 MNS clusters dropped from every Malawi DHS aggregate, with tests. 518 of 575 columns exist in all four countries; 78 are static; 136 sit at the exact survey year everywhere.

Two under-sold facts a policy audience needs more than any of the above. First, **the load-bearing layers are the free, static, permission-less ones**: the climate-only arm scores 0.335 under leave-one-country-out and soil-only 0.327, each above the full 575-column index at 0.302. Second, **the household-survey microdata blocks that cost weeks of RA time make cross-border prediction slightly worse** (MICS drop delta −0.013, HCES −0.009) and add under 0.01 in-country; the fully open tier transports best of the three (0.316 open / 0.302 headline / 0.281 with DHS; in-fill 0.400 with DHS vs 0.399 without). A ministry with no survey partner needs none of the survey data.

### Verified weaknesses

**The headline set is nearly blind to programmes.** Of the 35 columns in "Food fortification and supplementation", 25 are national constants and 7 are DHS (excluded from the headline arm). What actually enters the design matrix is **2 columns in every LOCO and pooled fit** (`mics_salt_iodised_15ppm`, `mics_salt_any_iodine`) and 3 in Malawi in-country (adding `mics_c_vas_6mo`), out of ~329–362. Lead with the standalone arm, not the drop-delta (the delta is 0.000, but 8 of 21 domains are ≤ 0.000 so it proves little): the domain's "only" arm scores **0.032** on the level and **−0.108** on prevalence against climate alone 0.335 and soil alone 0.327, and its median importance share is 0.007–0.008, 19th of 21 domains. The precise statement: no district-level reach for wheat flour, maize flour, oil, sugar or rice; none for IFA or deworming; VAS only in Malawi; salt only at Admin-1 granularity in Ghana (10 distinct values over 260 districts) and Malawi (31 over 243). Consequence: the ranking is need-as-geography, so it cannot separate "high risk because unreached" from "high risk despite coverage", and cannot evaluate a programme. One honest qualifier for the ask: where a district reach variable does exist it tracks status hard — district salt iodisation against low-UIC prevalence is Spearman −0.48 in The Gambia (27 districts) and −0.77 in Sierra Leone (14) in `iodine_in_country.csv` — though from the survey's own salt measure, and iodine is outside the 24 headline cells. The build is feasible: LSMS item modules are on disk for Gambia, Malawi and Sierra Leone with item maps already tagging flour, oils, sugar and salt (Ghana is aggregate-only).

**"District" is not one unit across countries.** From `metadata/admin2_spine.csv`: The Gambia's 30 surveyed units are Districts (70 clusters); Ghana's 75 are 37 Districts + 36 Municipalities + 2 Metropolises of 260 (90 clusters); Sierra Leone's 14 are its whole Admin-2 (60 clusters); Malawi's 87 are Traditional Authorities and sub-chiefdoms nested inside its 27 Admin-1 units — **which are what a Malawian programme officer calls a district** (103 clusters). "Traditional Authority" appears nowhere in the deck. Two sub-claims do *not* hold and must not be repeated: Sierra Leone does not contaminate the in-country headline (all 12 of its in-fill rows and all 60 of its in-fill targeting rows are NaN, so the 0.598 pair statistic is identical with or without it), and Malawi's predictors are not coarse (61.6% of shared predictors vary at full Admin-2 resolution in Malawi against 54.4% in Ghana; `consistent_district_tier.csv` shows the tier choice moves Malawi by ≤ 0.03). Where Sierra Leone **does** distort is cross-border: its 6 cells average 0.402 pair concordance — below a coin toss — against Gambia 0.621, Ghana 0.580, Malawi 0.608, so the 22-cell mean of 0.546 becomes 0.601 without it, and transport Spearman 0.302 / 0.223 becomes 0.355 / 0.293.

**The CONSORT numbers describe the wrong tier and are stale.** The exclusions slide (qmd 1405–1416) reads 614 / 575 / 575 / 525 / 464, but the 525/464 pair is the `open+survey_public+survey_dhs` row of a 15 Sept design CSV while the deck's results come from `open,survey_public`; the same table's headline row is 375/326, and live metadata adds the 8 ACLED columns. The 614 literal is unsupported — the cited log reads 588 → 575. One chain, one tier: **588 built → 575 after 23 written exclusion rules → 383 vary within a country and need no DHS → ~333 present in all four countries at ≥ 70% coverage with variation → 94 / 125 in the two pre-registered deployment candidates**.

**Temporal misalignment does not explain the weak cells**, and the deck should stop implying it might: Ghana has the best-aligned headline set (mean gap 0.58 y, 242 exact) and mid-table results; The Gambia has the worst microdata gaps and the best results; Sierra Leone's 145 columns ≥ 3 years off are AlphaEarth, MICS and HEAT, none of them load-bearing. The honest causes of weak cells are 1.2 clusters per district in Ghana and Malawi, 14 units in Sierra Leone, and the Sierra Leone child file being the anaemic subset only (532 of 654).

Lower-severity items worth one appendix line each: 72 of Ghana's 383 headline columns are region-broadcast (MICS at 10 old regions, HEAT, MIMI), so the district detail in the Ghana worked example comes from rasters, not household surveys; `source_registry.csv` is stale against the live set (ESPEN, `ghsl_smod`, GLW4 year); `variable_sheet.csv` still reads "needs_RA" for 288 of 608 rows, so "every column carries a definition" is true for source and year, not for a written definition.

### The one-slide description for policy

> **575 things we already know about every district, none of them from a blood sample taken for this project.**
> - **What it is:** 575 measurements for each of 554 districts in four countries, from 28 public sources — rainfall, temperature, greenness, night lights and a 64-number satellite fingerprint; soil chemistry; crops and livestock; malaria and anaemia maps built from earlier surveys; helminth treatment records; market prices; and the public household surveys (MICS, budget surveys) summarised to districts.
> - **What we deliberately left out:** anything measured on the same people whose blood we predict. The DHS is out of the headline models entirely; where a MICS or DHS was fielded inside the biomarker survey's own clusters, those clusters are dropped.
> - **What carries the ranking to a new country:** the free, static layers. Climate normals alone rank a never-seen country's districts at 0.34 and soil chemistry alone at 0.33, against 0.30 for all 575 columns. No household survey is required.
> - **What it cannot see:** programme reach below the national level — vitamin A rounds, iron-folate distribution, fortified-vehicle coverage. Two district-level programme columns reach the model; the rest is one number per country.
>
> Speaker note: "Every layer is downloadable today at no cost; three need a free account. What the ministry adds is the thing no satellite sees: where the programmes actually are."

### What a new country needs

| Tier | Contents | Common columns in Côte d'Ivoire | Effort | Measured return |
|---|---|---|---|---|
| 1 — the deployed product | TerraClimate normals + survey-year anomaly, MODIS LST, SoilGrids v2 / iSDAsoil, SRTM, Köppen + AEZ. Needs an Admin-2 polygon file and an Earth Engine account | **94** | 1–2 days (the CIV climate normals were "a morning's work") | the pre-registered climate+soil candidate |
| 2 — second candidate | + IHME anaemia surfaces, Malaria Atlas, MapSPAM, GLW4, ESPEN, Copernicus land cover | **125** | 2–3 days (CIV went 132 → 212 columns in one session) | the five-domain candidate |
| 3 — in-country only | + MICS microdata (registration; GPS needed for district resolution), HCES diet modules (~a week per survey for item coding), WFP/RTFP prices, ACLED | — | 2–4 weeks, RA-level | < 0.01 in-fill, **negative** across borders; keep for interpretability, not accuracy |
| Not needed at any tier | DHS microdata (in-fill 0.400 vs 0.399; transport −0.02), national series, LSMS region means, FluNet | — | — | — |
| No public source has it | DHIS2/HMIS district VAS and IFA counts; the biomarker survey itself if a scored in-country ranking is wanted | — | — | — |

---

## 1.4 The results

### Headline-number verification

| Quantity | Deck value | Verified | Verdict |
|---|---|---|---|
| In-fill index, level / prev | 0.399 / 0.282 | 0.399 / 0.282 over **18** cells (summary prints `cells=24`) | OK; fix the label |
| In-fill jackknifed regional mean | 0.314 / 0.220 | same | OK |
| In-fill spatial smoother | 0.390 / 0.257 | same | OK — it is the stronger survey-only comparator and the deck compares against the weaker one |
| Index beats regional average | 13 of 18 | 13 of 18 | OK |
| Pairs: index / regional / spatial / transport | 59.8 / 55.5 / 59.1 / 54.6 % | same; transport 60.1% excluding Sierra Leone | OK; add the SL split |
| Strong-cell pairs | 72.2% | **69.9%** exact (0.722 is an arcsin approximation) | Fix |
| In-fill prev wMAE: index / regional / national / smoother | index plotted alone | **14.25 / 11.77 / 11.49 / 10.93**; calibrated index **10.69** | **Speaker text wrong in direction** |
| Transport: full / climate+soil / five-domain / honest nested | 0.302 / 0.378 / 0.414 / 0.312 | all reproduce; null 0.079 | OK; the middle two are post hoc |
| Burden capture: index / regional / null / oracle | 21.9 / 20.3 / 13.6 / 47.5 % | same; lift > 1 in 11 of 18, p ≈ 0.5 | Reproduces, but at chance |
| Prevalence in flagged fifth | 32.2 vs 24.3, oracle 41.3 | 32 vs 24, regional 30, **oracle 47 on the matched 18 cells** | Use the matched oracle |
| WHO vitamin A bands | 64 / 91% vs 66 / 94% | same | OK — the regional average wins |
| Worst-third hit rate | 52.4% vs 33.0% | same, 1,176 held-out districts | OK — this is a **leave-one-country-out** statistic and must not be mixed with in-fill figures |
| Rank-interval coverage / honest width | 38% (16–70); 55% of list | same | OK |
| Anchor A1 at 5% | 9.3 pp | 9.32 row-median / **10.05** cell-median; a flat national number **8.62** | Inconsistent aggregation, and the flat number is better |
| Ceiling / share of attainable | 0.480 / 59% | 0.480 over 21 rows including four zeros; matched prevalence share 0.58, **level share 0.66** | Do not quote a share |
| `tr_enet` in quantities.json | 0.266 | **0.260** after the 18 Sept run | Re-render the deck |

Six tables are stale (pre-RR-12, 1–9 Sept) and structurally excluded from `rerun_downstream.sh`: `headroom_by_cell.csv` — which has **no producer anywhere in the repo**, so it must be rewritten rather than re-run — plus `consistent_district_tier.csv`, `model_augmented_survey_summary.csv`, `survey_size_symmetric_summary.csv`, `comparator_fairness.csv`, `fold_draw_risk.csv`. Drift is real and like-for-like (Gambia child vitA 0.745 stale vs 0.694 current), caused by the 15 Sept outcome reconciliation (iron binaries off the IDA columns, VITA_RULE → rbp070). None currently feeds a slide; cite them as "earlier run, direction only" if at all.

### The ten policy numbers, with framing and caveat

1. **60% of district pairs vs 56% for the regional average and 50% by chance.** "Ask the model which of two districts is worse and it agrees with the survey six times in ten; the number a ministry uses today manages five and a half." Caveat: a covariate-free neighbour smoother gets 59.1% — much of the in-country signal is geography, and the covariates earn their keep where no smoother can be fitted.
2. **0.40 against 0.31, better in 13 of 18 cells.** "Better than the survey's own regional averages at ordering districts." Caveat: on level, burden and severity bands it is not better; the sentence needs the word *order*.
3. **Six of 18 cells at 0.53–0.70** (Gambia child vitA 0.69, Gambia women's iron 0.68, Gambia women's vitA 0.70, Ghana child iron 0.53, Ghana women's B12 0.57, Malawi women's B12 0.70; 65–76% of pairs). Caveat: the 0.5 cut was chosen after seeing the scores and three of six are one country — describe, do not promise.
4. **Transport 0.30, 18 of 22 cells positive, chance 0.08.** Caveat: this is the honest number; 0.38 and 0.41 are post-hoc candidate sets. Spread by held-out country is large (0.60 / 0.35 / 0.20 / 0.16) but each is a single realisation with no interval and reflects a training-set-and-target pair, so state it as spread, not as a per-country deployability rating.
5. **Worst-third call right 52% vs 33% by chance; a top-five call right 62%.** The plainest sentence the project has, and it is the *new-country* number, so it is the right caption for Côte d'Ivoire. Caveat: "right" means the survey's own noisy rank; by held-out country 56 / 52 / 52 / 50%.
6. **The flagged fifth runs at 32% deficiency against 24% nationally; a perfect map 47%.** Caveat: a mean over outcomes spanning 0.2% to 81% prevalence.
7. **Burden: the flagged fifth holds 22% of deficient people; the regional averages 20%; chance 20%; a perfect map 48%.** Frame as the honest negative — on where the cases are, the model and the regional averages are both at chance, because burden sits in populous districts. Do not headline the 2-point margin.
8. **District prevalence: the calibrated index is off by about 10.7 points, against 11.5 for one national number and 11.8 for the survey's regional averages.** "About eleven points off, a little better than the national figure or the survey's regional averages — which at district level are themselves no better than the national figure." Caveat: only from the calibrated arm; the un-shrunk ranker is 14.3 and must never appear on a prevalence axis.
9. **Anchor-and-rank: with only a national estimate, the model supplies an order (Spearman ~0.25) at the national number's level (~8.5–10 pp, depending on the anchor's own precision).** Caveat: a regional survey at 15–25% of full size wins on level; the anchored ranking is at chance on burden; never say "replace".
10. **Every added training survey raises every other country's transported ranking: 0.20 → 0.26 → 0.30 with one, two and three training surveys (13 of 16 paired cells improve); climate+soil 0.33 → 0.38 → 0.40, positive in 95% of single-training-country fits.** The pooling argument, and the best next-steps figure. Caveat: three points from four countries on one continent; the gain is in Gambia and Ghana as held-out countries and flat for Malawi and Sierra Leone.

### Oversold vs undersold

Oversold: the level SAY ("only a little better than one national number" — it is worse); "better than the survey's own regional averages" without the word *order*; "close to the best any model could be" (uncertain ceiling, cells selected on the same score); "a whole unsurveyed region: the ranking holds at 0.39" (the covariate-free smoother reaches 0.377, so the region estimand does not demonstrate covariate value — only the country estimand does); "59% of attainable"; the SuperLearner race the index wins by 0.01; transport pairs 54.6%, contaminated by Malawi women's vitamin A at 0.1% national prevalence (concordance 0.965) — say 53%, or say 60% for the three countries where it transports and at chance in the fourth.

Undersold: The Gambia as a whole (0.41–0.70 in-fill on all four outcomes; 0.59–0.77 when held out entirely) — the demonstration that the method works when the survey lets it, and the argument for the survey-design ask; the class-probability statements (52 vs 33; 62%); the five-domain candidate's 22 of 22 positive cells; climate+soil already positive in 95% of single-training-country fits; fold-draw stability of the strong cells (SD ≤ 0.03).

### Underused results worth surfacing

Per-country in-fill **with the smoother row**: Gambia 0.62 / jk 0.41 / smoother 0.54; Ghana 0.44 / 0.39 / **0.50**; Malawi 0.26 / 0.21 / 0.23 — the covariates' in-country value is Gambia and Malawi, not Ghana. Malawi at its true district rung (`consistent_district_tier.csv`, stale): women's iron 0.24 → 0.60 on 27 real districts, which is the Malawi number a programme officer would recognise. Worst-fifth probability and its calibration (16 / 27 / 34 / 36% across four bands; 1,266 district-outcome rows) — a ready-made "how much to trust a flag" figure with ground truth, where the deck shows only the CIV version without it. The model-augmented-survey test: putting the model's within-region pattern on top of a survey's regional level made district estimates *worse* (12.2 vs 9.3 pp at full sample) — the only tested way of blending survey and model failed, and it is not in the deck, though it is the honest answer to "models together with survey data". Urbanicity strata (Malawi rural-tercile child iron −0.35): the "where not to trust it" slide.

### Which models won, plainly

The zero-tuning principal-component **domain index** — rank-normalise within country, PCA per domain, weight components by training-fold Spearman, sum. Against the survey's own regional averages it wins the order in 13 of 18 cells; against a covariate-free neighbour smoother it ties in-country (0.399 vs 0.390; ahead in 13 of 18 paired, behind in Ghana) and wins where no smoother can be fitted; against a full SuperLearner library of 14 candidates on identical folds it has the best mean Spearman on 4 of 6 estimand × target combinations and trails by 0.00–0.03 on the other two, with every ensemble rule at or below it; against the DHS-style geostatistical model it ranks better (0.384 vs 0.271 over 24 cells) and measures worse (9.23 vs 10.66 pp **in that comparison's own scoring run** — do not set 9.23 against the 14.25 from a different run).

Where it fails: levels without calibration (14.3 pp); burden at every estimand; zinc (−0.07) and any outcome below 1% national prevalence; held-out Malawi (0.20) and Sierra Leone (0.16) with the full vocabulary; Malawi's rural tercile (−0.35 on child iron); uncertainty (90% rank intervals cover 38%); and combination with survey data beyond the national anchor.

---

## 1.5 The national-level (VMNIS) track

**What it found.** On the WHO VMNIS panel for preschool vitamin A (108 country-years, 69 countries), a random forest on 17 World Bank indicators with the country held out reaches Spearman 0.625 / Pearson 0.655 and MAE 11.75 pp against a no-covariate null of 16.23 (ridge 12.88 / 0.562). Folate in non-pregnant women reaches ~0.51 (32 country-years), zinc 0.36 by ridge; women's vitamin A and B12 show nothing. But composing that predicted national level onto a transported district map makes the district estimates **worse**: over the four non-degenerate cells, MAE is 9.10 pp with no level, 10.51 with a pooled-null level, **19.12** with the VMNIS-predicted level and 8.62 with the country's own true national number; signed bias −5.31 / +4.16 / **+15.05** / +0.98. Per country the predicted national levels are Gambia 20.27 (true 19.78), Ghana 20.77 (14.93), Sierra Leone **40.90** (12.19), Malawi **31.21** (9.18).

**Why it fails, in the checkable form.** The vitamin A / preschool panel contains **none of the four project surveys** (Ghana and Sierra Leone have zero rows; The Gambia's six are 1999, Malawi's six 2001), 515 of 528 rows are serum retinol — never RBP, the project's analyte — inflammation adjustment is recorded for only 33.5% of rows, and the panel holds **30 African country-years, 22 of them pre-2005 averaging 50.3 pp and only 4 after 2010 (25.0 pp)**, against 79 non-African country-years at 16.9 pp. The ridge therefore maps Sierra Leone 2013 and Malawi 2015 onto that pre-2005 African cloud. Scope correction, for accuracy in Q&A: the four surveys **are** present in the VMNIS folate and B12 panels (Ghana 2017, Sierra Leone 2013, Malawi 2016) and Malawi 2016 zinc; the absence is specific to vitamin A.

**Does it earn a slide?** No main-talk slide; one sentence in the "models alone vs models with survey data" section, and an appendix slide as the backup for "why not model the national number too?". The sentence is a design rule, and it is the rule the anchor-and-rank slide depends on: *"A model cannot replace the national blood sample: predicted national levels miss by about 12 points on average and by 22–29 points in two of our four countries, and a wrong anchor moves every district by the same amount — district error doubles, from 9.1 to 19.1 points."* A second, more interesting framing for this audience is the measuring-stick result: for women's vitamin A the between-method variance (sd 2.10 logit) exceeds the between-country variance (1.17), i.e. how a survey measured explains more of the gap between two countries' reported numbers than which country it is — a harmonisation ask, on 20 countries, worth at most one spoken sentence.

**Caveats that must travel with any use.** The composition tables were written 31 Aug 2026 from the **retired pre-protocol-v2 elastic-net transport pattern** (`transportability_area_loco_predictions.csv`), not the v2 index; they carry stale survey years (Gambia 2021, Malawi 2015 against `metadata/survey_years.csv`'s 2018 and 2016), which is load-bearing because the target country's covariate row is picked at the nearest year; the script header's "null beats the covariate model in 6 of 8 cells" is **5 of 8** in the table; nothing in the deck reads these tables; and the deck's appendix slide titled "National-level track" actually renders `national_estimates_all.csv` — the within-survey person-level SuperLearner aggregate, i.e. self-consistency, not national prediction — so that slide is mislabelled and should be relabelled or replaced. Re-pointing script 19 at v2 is blocked: no v2 script writes per-district LOCO predictions, so that export must be built first.

---

## 1.6 The full talk deck as a communication object

92 slides: 66 before the appendix (57 content + 9 dividers) + 26 appendix. At the Forum's one-slide-per-two-minutes rule the 15-minute budget is 7–9 slides, so the full talk is a **source library**, not a draft to trim — about one content slide in eight survives. On its own terms it is an excellent technical record: every number is an expression over one `Q` list, the SAY/DETAIL split is disciplined, and the withdrawn claims are genuinely absent.

**What works.** Slide 3 (national number → 16 regional averages → all 260 districts) is the entire value proposition in one image and its note already says so. Slide 5's "big-data opportunity, small-data problem" is the best framing sentence in the deck. Slides 22–23 (three programme questions: a missed district, a missed region, a country with no survey) are the right organising device, and the held-out cartoon is the only method figure a policy room needs. Slide 35 (Côte d'Ivoire, named districts, WFP/MIMI corroboration) is the most quotable result. Slide 41 is a genuine graded decision inventory with the status-quo comparator on each row. Slide 62 costs four designs against each other in survey-size units. Slide 65's yes/partly/no checklist is the best credibility device in the deck. The SAY notes are 40–70 words, present tense, one idea each — a usable script spine after a register change.

**What loses a policy room.**
1. **Structure.** The Results divider is slide 24 and the first per-cell accuracy slide is 25, after 15 slides of inputs and protocol; the Côte d'Ivoire deliverable is 35, the decision grid 41, the costed design 62, the close 65–66. Everything the room came for is in the back half. (Fair caveat: the qmd declares this "the 30–45 minute talk", so this is a hazard for the distillation, not a defect of the source.)
2. **Repetition.** Eight main-body slides translate the same in-fill result — 30, 37, 38, 39, 40, 41, 65, 66 — with 36 and 61 restating it in other units, and appendix 84 repeating 40. The notes read "DRAFT A / B / C / D of four candidate translations", and slide 37's note says keep *this* one, which conflicts with the natural instinct to keep 41: a live branch point for the PI. Slides 31 (Ghana) and 42 (Malawi) are the same four panels on the same outcome.
3. **Jargon**, counted on-slide: "in-fill" 10×, "fold(s)" 15×, "held out" 12×, "ceiling" 7×, "correlation" 6×, "Spearman" 4×, "geostatistical" 4×, "leave-one-country-out" 3×, "SuperLearner" 2×, "permutation null" 2×; "outcome-country combination" ~20× on-slide and 60×+ in notes.
4. **Text-processing damage in the notes**, which matters because SAY is the script: "inherited red-outcome-country combination disorders", "country-outcome outcome-country combinations", "optimiztic", "Status of the January analyzes".
5. **The worked example is the wrong cell.** Ghana child vitamin A carries the opener trio, the exceedance map and the policy deck's three-map fig7 — and it is a cell where the index **loses** to the survey's own regional average on both targets (0.277 / 0.265 vs 0.410 / 0.314; prevalence error 21.2 pp vs 15.2), with the regional average winning in 9 of 10 fold draws. fig7's caption benchmarks 0.31 against chance (0.08) and never against the 0.314 regional average that beats it. The Goals slide's note even calls Ghana child *iron* "the worked example that recurs through the talk"; the next slide switches. Ghana women's B12 is the only Ghana cell that beats the regional average on both targets (0.566 / 0.486 level; 0.317 / 0.235 prev); Ghana child iron is defensible if the claim is stated on the level/ranking target (0.53 vs 0.36) and the prevalence tie (0.502 vs 0.498, with worse error) is said aloud. Scope: across the 36 in-fill country-outcome-target rows the regional average beats the index in 12, so this is a presentational, localized defect — the deck reports the aggregate honestly.
6. **`fig3_top_predictors.png` is not presentable, and it carries the 3-minute short oral.** Four of 30 bars show raw codes ("mics salt iodised 15ppm", "mics salt any iodine", "mics heat vtetprot sy", "u5 diarrhoea prev"); 9 of 30 are in the wrong colour group (MICS and HCES columns coloured "Remotely sensed environment"); the "Household survey (DHS)" legend entry is dropped entirely; the speaker note describes bars that no longer exist ("fewer women owning their home or in paid work, less schooling"); and the same teal/purple that mean "our model / DHS-style model" on figs 1, 2, 5 and 9 mean data *type* here. Minimum fix is code-level: add `mics_`/`hces_` branches to `group_of()` in `scripts/policy_deck/01_figures_main.R`, rename the group "Public household surveys (DHS, MICS, budget surveys)", add the four missing `PLAIN` entries. If a single-outcome redraw is chosen, use women's B12 or child vitamin A — child iron's own top five leads with `ihme_overweightprevalence` at beta −4.33, which is harder to defend than the salt codes.

**The 12 strongest slides for policy, with the fix each needs:**

| # | Slide | Why it carries weight | Fix before use |
|---|---|---|---|
| 1 | 3 — what a country has today vs what the model adds | The whole value proposition in one image | Move off Ghana child vitamin A; caption "185 of 260 districts never surveyed" |
| 2 | 35 — Côte d'Ivoire ranking | A real country with no survey, ranked from free data; corroborated on vitamin A (rho ≈ 0.45 over 32 districts vs the WFP/MIMI map) | Say "stable, not verified"; drop the near-black P(worst third) row |
| 3 | 41 — which decisions is it good enough for | The only graded decision inventory, with the status quo on each row | Re-grade "sequencing a roll-out": on burden and bands the model is not better than the regional averages |
| 4 | 62 — anchor-and-rank | The one result a funder can cost, in survey-size units | Reframe off substitution; add the regional-survey row and the burden loss |
| 5 | 33 — does it work in a country it has never seen | The central challenge, with a permutation null | Quote 0.30; label 0.38 / 0.41 as candidates; replace "Spearman" with "ranking accuracy" |
| 6 | 5 — big-data opportunity, small-data problem | Explains why the simple model won without saying "shrinkage" | Keep "14 to 87 districts per country is the sample size of every model" |
| 7 | 36 / fig6 — burden captured by the worst fifth | The only slide that speaks in people | Say the honest version: at chance, and so are the regional averages |
| 8 | 65 — the can / cannot checklist | The best credibility device in the deck | One number per row; rewrite the coverage row as plain odds |
| 9 | 21 — how the model works, four steps | Enough method to be trusted, no more | "A few summary scores per family of data" |
| 10 | 54 / fig3 — what it looks at | Answers "what is it looking at" | Rebuild the figure (above); keep the footer "these weights say where deficiency is, not what causes it" |
| 11 | 52 / fig5 — learning curve | The pooling ask: every survey improves every other country's map | Keep as is |
| 12 | 25 — per-cell spread | The truth that it works where the survey resolves districts | Too dense at 18 rows; compress to a two-bar statement or say the sentence |

**Missing slides:** a named user and a named decision (the deck names no institution as user; "ministry" appears only in speaker notes); any money figure or timeline; a "with survey / without survey" grid (rows: surveyed country with unsurveyed districts / a region never reached / a country with no survey / planning the next survey); a one-district before/after card; a single plain-accuracy panel; a map of where the next survey should sample; and the ask — the 9 Sept outline's funded activity list (Ethiopia, Pakistan, Tanzania, the Proxy Modeling Alliance, the MIMI/GBD working group, two regional workshops) and its closing question, "what predictive performance is good enough to be worth acting on?", appear nowhere in the deck, and `git log -S` shows they were never in it.

**Jargon → plain (the dozen that matter):** Spearman → "ranking accuracy: 0 is chance, 1 is perfect" and better still the pair count; estimand → "the three questions we tested"; in-fill → "filling gaps inside a surveyed country"; leave-one-country-out → "we hid a whole country and predicted it"; principal components → "each family of data summarised into a few scores"; fold / ten draws → "every district predicted with itself hidden, ten times over"; jackknifed regional average → "the region's average computed without that district"; permutation null → "chance (0.08)"; reliability ceiling → "the most any model could score against a survey that samples one community per district"; population-weighted MAE → "how many percentage points off, on average"; SuperLearner → "a competition of algorithms scored on data they had not seen"; transport → "works in a country it has never seen"; outcome-country combination → "each nutrient in each country".

**Visuals — keep / cut / missing.** Keep: `fig10_civ_ranking.png`, `fig1_model_comparison.png`, `fig6_targeting.png`, `fig5_learning_curve.png`, `fig11_worst_fifth_probability.png`, the Ghana source-map panels (one set, not both), the Hess 2023 framework with the coverage overlay, and `fig13` top row only. Cut or send to the appendix: `fig3` as currently drawn, `fig2_simpler_wins`, `fig4_domain_scatter`, `fig8_twenty_layers`, `fig9_geostatistical`, `fig12_vim_forest_*`, every `index_importance_top10_*.png` (variable codes), and the January ensemble cartoons (`ensemble_orchestra`, `ensemble_vfold`) — showing them as the method contradicts the deck's own finding. Build: a today-vs-model trio on a cell the model wins; a one-district card; the decision grid; a plain-accuracy panel; the survey-design cost curve with plain labels; and a "sure vs right" scatter (model firmness against rank error) that shows honestly what the pale/deep legend does and does not mean — because it does not mean what the current legend says. Verified: the worst-fifth probability is only weakly calibrated (35.7% of the 157 districts at p > 0.8 are truly in the worst fifth, against 15.8% at p ≤ 0.2 and a 20.4% base rate), and within country-outcome the correlation between uncertainty and rank error has a median of **−0.018** across 18 cells. Reword the legend to "deep = stable under retraining; pale = would move with different training data", and make model-guided cluster placement an explicit next-step test rather than an instruction.

---

## 1.7 Findings raised and refuted in verification

Each was investigated and did **not** survive. Residues are noted because they are cheap to fix.

1. **"Modelled surfaces are recycled surveys; Malawi has an unchecked same-household path through IHME; the drop-modelled sensitivity was never run."** Refuted. The Gambia has the same nested design (70 MICS6 clusters, not only Malawi's 105 DHS clusters); the second-order path is written down in three places in the repo and both nested designs are named on-slide; and the substantive question has been answered three times for LOCO — dropping all 29 IHME columns costs +0.010 on the level and **improves** prevalence by 0.007, dropping Malaria Atlas costs +0.009 / +0.013, and adding the whole anaemia block to climate+soil gives −0.004 / −0.003. Live exposure is 40 of 383 subnational headline columns, not 48 (8 are national constants dropped at fit time, 6 are Ghana-only MIMI). Residue: the 11 Malaria Atlas columns are **not** flagged `modelled_surface` although the deck's note says the switch covers them; there is no ablation for the in-fill or region estimands; and the "open tier" label overclaims, since IHME and MAP surfaces are themselves fitted to earlier DHS rounds.
2. **"Per-country transport numbers should replace the 0.30 / 0.38 mean on a slide."** Refuted as a recommendation. The ceiling column that would sit beside it is an artefact (3 of Sierra Leone's 5 prevalence rows are singular zeros; dropping them makes Sierra Leone the *second highest* ceiling, not the lowest); the proposed mechanism is contradicted (Sierra Leone has the most clusters per district, 4.26, and no single-cluster district, yet the worst transport; Ghana has 1.20 and 83% single-cluster districts and the second best); country explains 43% of the sum of squares against outcome's 26%, and within-country SD (0.194) ≈ between-country SD of means (0.199); Malawi's low in-fill is largely its two Malawi-only zinc cells (0.260 → 0.369 without them); Sierra Leone has no in-fill or region values at all, so the table has a structural hole; and the LOCO rows have reps = 1 and no intervals. The deck already carries a per-cell three-facet forest with the regional average on every row. Residue: the spread is real — say it in one sentence, do not build a causal four-row table.
3. **"Nothing in the deck answers 'what would a programme officer do differently' — no user, decision, cost or ask."** Refuted on three of four. Slide 41 is a five-row graded decision inventory carrying the status-quo comparator; slide 62 costs four designs in survey-size units; three uses are quantified against current practice (pair order, burden capture, WHO bands), not one. Residue, and it is real: no named institution, no currency figure, no timeline, no ask — plus a sharper point the original finding missed, that on the two uses actually scored against current practice the model does **not** beat the regional average (burden 22 vs 20%, bands 64 vs 66% exact), so slide 41's "yes" for "sequencing a roll-out" is a mis-grade. Note also that converting slide 62 into a dollar pair would revive the savings framing the project has already withdrawn.
4. **"Withdrawn WS-document numbers (0.516 flat regional mean, +0.253 anchor gain, 0.7 pp subsample benefit, WS4's r = 0.043) may leak into the talk."** Refuted. The deck loads only protocol-v2, cluster-level, policy-deck and national tables and computes all 287 quantities at render time; grep for every withdrawn literal returns zero hits in the qmd, the concept yaml, the quantities file and the prose. Two of the cited items invert authorship — WS2 and WS6 are the documents that *withdrew* those numbers. Residue: `WS2_ANCHORING_CONTROLS.md` and `WS5_ANCHORING_BUDGET.md` still assert the 0.516 flat-regional-mean claim with no withdrawal banner, and they are the files a drafter would naturally open for "how does the model compare with the regional average". Add a banner; keep sourcing numbers from the v2 tables.


# Part 2. The recommended 15-minute talk

**Scope note.** Everything below is designed AT the Forum budget: 16:9, one slide per two minutes of speaking, main body = title + 2 Sonja + 5 Andrew = **8 slides, 14.5 minutes**. Sonja 3.0 min, Andrew 11.25 min, 0.25 min close. Every number is traceable to a named table; none is on the withdrawn list. Comparator convention used throughout: **`region_mean_jk`** = the survey's own Admin-1 mean with the district itself left out, respondent-weighted (`n_eff`). Where an unweighted variant changes a number, it is flagged.

---

## B0. Title options and the one sentence the room should repeat

The submitted abstract title stays in the programme and in 14-pt grey on slide 1; the room has read it. The on-screen headline should be plainer, and should not promise a prevalence.

**Title options**

| | Title | Why / trade-off |
|---|---|---|
| **A (recommended)** | **"Which districts first? Ranking every district for micronutrient deficiency from free public data"** | Names the decision and the product in one line. "Ranking", not "predicting prevalence", so the title itself cannot be quoted against the LEVEL result. 12 words. |
| B | "Between surveys, every district gets its region's average. We can do better than that." | Strongest rhetorical hook; sets the comparator on screen before any number. Risk: it is a criticism of survey practice in a room full of survey teams, and it is also the thing Sonja has just thanked them for. |
| C | "A district map for countries between surveys — an order, not a prevalence" | Most honest, most self-limiting; best if the PI expects a hostile methods room. Loses the "which districts first" decision framing that policy people respond to. |

**The one sentence the room should repeat (three options)**

| | Sentence | Notes |
|---|---|---|
| **A (recommended)** | **"An order, not a number."** | Four words, survives the corridor and the coffee queue, exactly matches what the tables support. Say it three times: slide 4 build 2, slide 5 grade strip, closing line. |
| B | "Three yeses, one partly, two noes." | From the decision-first proposal. Excellent *as a caption for slide 7*, weaker as the talk's thesis because it means nothing without the ledger on screen. Use it on slide 7 only. |
| C | "A shortlist, not a measurement." | Slightly warmer than A, and "shortlist" is programme language. Loses the explicit contrast between order and number that the LEVEL sentence needs. |

**Do not use as the thesis** any sentence longer than about 12 words. All five judged proposals wrote 45-55-word theses and every judge said nobody will repeat them.

---

## B1. Recommended spine at a glance

**Which proposal this is based on, and why.** The spine is **BEFORE / AFTER / PRODUCT** (`scratchpad/proposal_before-after-product.md`, best-judged, mean 39.3/50), because it is the only structure that puts the deliverable on screen in the first three minutes, fits eight slides without apology, and prints its own negative result on the slide rather than leaving it in the speaker's mouth. Four grafts, each taken because every judge named it as the best single object in its proposal:

1. **The decision ledger** from DECISION-FIRST (proposal_decision-first.md slide 7) replaces BEFORE/AFTER/PRODUCT's thinner "how to use it" grid. It is the most-praised artefact across all fifteen judge reports and the only one that answers "what would a programme officer do on Monday" row by row.
2. **"We deleted Ghana's survey"** from ONE COUNTRY'S STORY (proposal_one-country-story.md slide 5) becomes build 3 of slide 5. Caption verbatim: *"This is the Côte d'Ivoire experiment, run where we have the answer key."* It converts an untestable deployment claim into a test with ground truth.
3. **The ceiling staircase** from DATA-STRATEGY (proposal_data-strategy.md slide 7) becomes the left half of slide 8, converting the reliability-ceiling defence into a costed survey-design ask.
4. **The two printed sentences (ORDER ✓ / LEVEL ✗)** are kept from BEFORE/AFTER/PRODUCT and are the one thing on slide 5 that may never be cut.

From HONEST SCIENCE the withdrawal ledger is taken to the **appendix** (Part G, tier 1) with one spoken clause, not a main-body slide: at eight slides it costs 20 seconds the talk does not have, and its credibility value is realised in Q&A.

### The 8-slide Forum spine

| # | Speaker | Claim-form headline | Min | Visual | Section |
|---|---|---|---|---|---|
| 1 | Sonja | *Predicting district micronutrient deficiency from public data — with four national survey partnerships* | 0.5 | Title, 19 authors, partner logos, "Unpublished data" | — |
| 2 | Sonja | *A programme picks districts. The evidence stops at the region.* | 1.25 | Ghana: one national number → 16 regional means → 260-cell coverage waffle (75 filled, 62 of them single-cluster, 185 empty) | Problem |
| 3 | Sonja | *Four national surveys, four partnerships — and a fifth country with none* | 1.25 | Four-country strip (districts / clusters / nutrients) + greyed Côte d'Ivoire column | Surveys and outcomes |
| 4 | Andrew | *The model adds a third map: every district, including the 185 no survey reached* | 2.75 | Map trio (shrunk before-pair + large Ghana **child iron** all-260 map), three named districts incl. one miss, bottom strip = what it is built from | Data sources / proxy dataset |
| 5 | Andrew | *The order is better than the region's average. The level is not.* | 2.75 | Odds ruler; the two printed sentences (✓/✗); "we deleted Ghana" twin panel; per-country grade strip | Modelling approach and what won |
| 6 | Andrew | *Côte d'Ivoire has no survey. Here is its ranking, and how firmly we hold it.* | 2.25 | CIV **women's vitamin A** ranking map + firmness strip + external-check corner + dashboard thumbnail | What drives it / product |
| 7 | Andrew | *Good enough to order districts and place survey clusters; not to count cases or judge a programme* | 2.25 | The decision ledger: 6 rows × Decision / What you use today / Verdict, revealed one row at a time | Using the models alone and with surveys |
| 8 | Andrew | *What limits the map is the survey's design — and every survey added improves everyone's map* | 1.5 | Ceiling staircase (clusters per district) + pooling curve with "your survey here" + three asks | Next steps and the ask |
| — | Andrew | Close, spoken over slide 8 | 0.25 | — | — |
| | | **Total** | **14.5** | | |

### The relaxed 12-slide variant (if the PI decides the budget is soft)

Each of Andrew's five compound slides splits at a real seam. Sonja's slides do not split; her 3.0 minutes stay as they are.

| Forum slide | Splits into | Min | What the split buys |
|---|---|---|---|
| 4 | **4a** *The third map* — map trio + three-district callout | 1.75 | The named-district moment gets room to land; the miss (Atiwa East) gets its own clause |
| | **4b** *What it is made of* — 575 layers / 28 sources, the access column (download / free account / agreement), the cost-to-acquire ladder (0.32 / 0.30 / 0.28) | 1.5 | Restores the PI's "proxy dataset construction" section and the single best dataset finding: the free tier travels best |
| 5 | **5a** *How we tested it* — the three hold-outs, the comparators, the odds ruler, the ✓/✗ sentences | 1.75 | "Held-out" gets defined properly; the comparator (regional average, neighbour map, chance) gets named |
| | **5b** *Where it works and where it does not* — per-country grade strip, the six strong cells labelled "chosen after the fact", Sierra Leone named as the failure | 1.5 | Pre-empts "will it work in MY country" with four numbers instead of one mean |
| 6 | **6a** *Côte d'Ivoire ranking + firmness + external check* | 1.5 | The product gets a clean beat |
| | **6b** *What the model is looking at* — climate normals, soil chemistry, livestock, wasting, prices; "where, not what to change" | 1.5 | Restores the PI's "proxy variables driving them" section, currently one spoken clause |
| 7 | **7a** *The three yeses* | 1.25 | |
| | **7b** *The partly and the two noes* — level, cases, programme evaluation, with the rates-vs-cases twin map | 1.25 | The burden "no" becomes a picture and a fix (rank by rate × population) rather than a confession |
| 8 | **8a** *The survey sets the bar* — ceiling staircase + clusters per district | 1.25 | |
| | **8b** *The ask* — pooling curve, pre-registration, three asks, the question back to the room | 1.25 | The ask is the last thing on screen, which is where asks belong |
| | **Total (12 main slides + title)** | **~19.5 min** | i.e. the 20-minute version, not 15 |

**Promote-first order** if the PI gains one slide only: **5b** (per-country heterogeneity) first — it closes the single most likely hostile question. Then **4b** (the cost-to-acquire ladder), then **7b**.
**Cut-first order** if the session runs late: slide 8's pooling curve (speak it), then slide 6's dashboard thumbnail, then slide 4's bottom strip. **Never** cut the ✗ LEVEL sentence on slide 5 or the two "No" rows on slide 7.

---

## B2. Slide-by-slide, with options and talking points

### Section: Problem

---

### Slide 1 — Title (Sonja, 0.5 min)

**Headline:** *Predicting district micronutrient deficiency from public data — with four national survey partnerships*

**Visual.** Forum 16:9 template. Submitted abstract title in 14-pt grey above the headline. 19 authors in two columns. Logos: UC Berkeley, UC Davis, the four national partner institutions, funder. **Red footer, on this and every result slide: "Unpublished data — please do not photograph the result slides."** (Forum guidance requires unpublished data be labelled; put the footer in the template so it cannot be forgotten on the USB upload.)

**Talking points (Sonja).**
- "Good morning. I am Sonja Hess; this is joint work with Andrew Mertens at Berkeley and with the teams who ran four national micronutrient surveys."
- "I will take three minutes on why we did this and who we did it with. Andrew then shows you the maps and what they are good for."
- "Some of what you will see is unpublished, so the slides say so."

**Cut if short:** no.

---

### Slide 2 — What a country has today (Sonja, 1.25 min)

**Headline (claim form):** *A programme has to pick districts. The evidence stops at the region.*

**Options for this slide**

- **Option A (recommended) — the three-part coverage build.** Left: Ghana filled with one flat colour and one number. Middle: the same outline broken into regional means. Right: a **260-cell waffle** — 75 filled, of which 62 half-tone (single survey cluster) and 13 solid (two or more), 185 outline-only. *Trade-off:* three objects in 75 seconds is tight, but the waffle does the entire "what a survey can and cannot see" job with no sentence, and it is reused on slide 8 as the survey-design ask.
- **Option B — the two maps only** (national number, regional means), coverage stated in words. *Trade-off:* cleaner and safer for Sonja's pace; loses the image that makes the ask on slide 8 land.
- **Option C — a fieldwork photograph with the two numbers overlaid.** *Trade-off:* warmer opening for a policy room, but "the region is all you get" is much weaker in words than in a picture of 16 flat polygons.

**Visual (Option A).** New figure; buildable today from `results/tables/policy_deck/viz/oof_child_iron.csv` (Ghana rows: `n_psu`) and `metadata/admin2_spine.csv`. Ghana child iron throughout.

**Talking points (Sonja).**
- "A national micronutrient survey happens roughly once a decade, and it is designed to report regions, not districts."
- "Ghana's 2017 survey did exactly what it was designed to do. It put 90 clusters into 75 of Ghana's 260 districts and it reports its regions well."
- "But the decisions are not regional. Somebody has to say which districts get the next round, which get screened first, where the next survey team goes."
- "Today the honest answer for an unsurveyed district is: it gets its region's average. In Ghana that is 185 districts out of 260 with no measurement of their own."
- "And in 62 of the 75 districts that were reached, the district's whole number comes from one cluster — one community, about a dozen children."
- "That is the gap. Andrew will show you what free public data can and cannot do about it."

**Numbers and sources.** 260 districts / 75 surveyed / 185 unsurveyed — `quantities.json` `exc_n`, `districts_Ghana`, `exc_unsurveyed`. 90 clusters — `quantities.json clusters_Ghana`. 82.7% single-cluster → 62 of 75 — `results/tables/protocol_v2/variance_components_ceiling.csv` `share_single_cluster`. National child iron **20%** — `results/tables/protocol_v2/nce_targeting_metrics.csv` `prev_national` = 0.2009 (population-weighted; **do not** use the 22.6% n_eff-weighted mean of surveyed districts, which three of the five proposals used and which a Ghana Health Service statistician will recompute differently).

**Cut if short:** no — this is the whole BEFORE argument.
**split_if_relaxed:** 2a = the decision and what exists today; 2b = the coverage arithmetic across all four countries (206 surveyed districts; three in four hold a single cluster).

---

### Section: Surveys and outcomes

---

### Slide 3 — The four surveys and the fifth country (Sonja, 1.25 min)

**Headline (claim form):** *Four national surveys, four partnerships — and a fifth country with none*

**Options**

- **Option A (recommended) — the four-column strip + greyed fifth column.** Per country: outline with cluster dots; survey name and year; three figures — *districts surveyed / clusters / nutrients*. A fifth greyed column, **Côte d'Ivoire — no survey**, which Sonja hands to Andrew.
- **Option B — partnerships and stakeholder engagement** (the January 2026 Accra meeting, the dashboard co-design, the country briefs), with the survey table reduced to one line. *Trade-off:* this is the slide funders like and it is missing from every current deck version; but the room then has no sample-size anchor for anything Andrew says, and "14 to 87 districts per country" is the honest reason the simplest model wins.
- **Option C — both, at 1.5 min**, taking 15 seconds from slide 8. *Trade-off:* the ask slide is already the thinnest.

**Visual (Option A).** New; from `quantities.json` `districts_*`, `clusters_*`, `outcomes_*`, `survey_*`, `year_*`.

| | The Gambia 2018 | Ghana 2017 | Sierra Leone 2013 | Malawi 2015-16 |
|---|---|---|---|---|
| Districts surveyed | 30 | 75 of 260 | 14 | 87 |
| Clusters | 70 | 90 | 60 | 103 |
| Nutrients | vitamin A, iron | + folate, B12 | + folate, B12 | + folate, B12, zinc |

**Talking points (Sonja).**
- "These are the four national biomarker surveys and the teams who ran them."
- "Between them: 206 districts, 323 clusters, 24 country-nutrient combinations. This is the only ground truth of this kind for this question on the continent, and it exists because these teams collected it."
- "A word on what 'district' means, because it differs: 30 districts in The Gambia, 75 of Ghana's 260, 14 in Sierra Leone, and in Malawi the 87 units are Traditional Authorities inside 27 districts."
- "Ghana is the worked example today, because the Ghana team is in this room and can check us."
- "And the fifth column is Côte d'Ivoire, which has no biomarker survey at all. That is the country this project exists for, and Andrew will come back to it."

**Numbers and sources.** `quantities.json` `n_districts` 206, `n_clusters` 323, `n_cells` 24; unit definitions from `results/tables/protocol_v2/targets_v2.csv` and `scratchpad/phase1_pipeline.md` F10.

**Cut if short:** trim the unit-definition line to a footnote (it is the most likely floor correction, so keep it somewhere on the slide).
**split_if_relaxed:** 3a = the four surveys; 3b = stakeholder engagement and the workshop.

---

### Section: Data sources and the proxy dataset

---

### Slide 4 — The third map (Andrew, 2.75 min) — **four-step build, one slide**

**Headline (claim form):** *The model adds a third map: every district, including the 185 no survey reached*

**Options for the worked cell — decide this first, it constrains everything**

- **Option A (recommended) — Ghana CHILD IRON.** In-fill level Spearman **0.529** vs regional average **0.356**; prevalence 0.502 vs 0.498 (a tie, concede it if asked). **DONE 18 September 2026.** `scripts/policy_deck/02_figure_ghana_map.R` is now parameterised (`FIG7_OUTCOME`, `FIG7_TARGET`) and both variants are built: `fig7_ghana_map_child_iron_level.png` (**ranking accuracy 0.52 against 0.36 for the survey's own regional average** — the version copied to `fig7_ghana_map.png`) and `fig7_ghana_map_child_iron_prev.png` (0.52 against 0.50 — a tie, so do not use it to claim a margin). The child vitamin A original is preserved as `fig7_ghana_map_child_vitA_prev.png`. The subtitle now names the regional-average comparator instead of the old "0.08 for chance", which was the cross-country permutation null and the wrong baseline for an in-fill map; the comparator is computed with the pipeline's own `arm_region_mean_jk_v2` on the same folds, and reproduces `benchmarks_v2_cells.csv` exactly (0.356 level, 0.498 prevalence).
- **Option B — Ghana women's B12** (level 0.566 in-fill, transport 0.594 — Ghana's best cell on the level target). *Trade-off:* stronger numbers, but B12 is not a programme the room runs; iron is.
- **Option C — keep child vitamin A** (what the current deck and `fig7_ghana_map.png` use). **Do not.** Child vitamin A is the one Ghana cell where the regional average beats the index on *both* targets (level 0.277 vs 0.410; prevalence 0.265 vs 0.314) and on district error (21.2 vs 15.2 pp, `benchmarks_v2_cells.csv`). Presenting the losing cell as the worked example in Accra is self-inflicted.

**Visual — four builds.**
1. The two BEFORE maps from slide 2 shrink to the left third.
2. The third map lands large on the right: Ghana child iron, all 260 districts, shared warm scale. Corner tag in the same type size as the legend: **"Order is what is validated. Not a prevalence."**
3. A red ring drops on three districts with a callout table (ranks among the 75 districts the survey reached — **print that caption at 10 pt or the room hears "1st of 260"**):

   | District | Its region's average says | Model says (district held out) | Survey found |
   |---|---|---|---|
   | Binduri, Upper East | 13th worst | **1st worst** | 3rd worst |
   | Talensi, Upper East | 10th worst | 6th worst | 5th worst |
   | Atiwa East, Eastern | 29th worst | 64th | 11th worst |

   The third row is the deliberate miss and stays on screen.
4. Bottom strip: four thumbnails — *30-year climate normals · soil chemistry · satellite land cover and livestock · public household-survey summaries* — with one line: **"575 public layers, 28 sources. No new blood sample, nothing to license. The map you will see next runs on 94 of them."**

**Talking points (Andrew) — a two-minute script.**
- "Here is the whole talk in one picture. One national number. Sixteen regional numbers. And then a value for every district in Ghana, including the 185 that no survey cluster ever reached."
- "Take three districts, and remember these are held-out predictions — Binduri was hidden from the model that ranked it."
- "Binduri: its region's average puts it thirteenth worst; the model puts it first; the survey found it third. That is a move a programme would act on."
- "Talensi: regional average tenth, model sixth, survey fifth."
- "And Atiwa East, where we are wrong: the survey found it eleventh worst, the model put it sixty-fourth. I will tell you in ninety seconds how often that happens."
- "The right-hand map is built from public layers: thirty-year climate normals, soil chemistry, satellite land cover and livestock, and summaries of other people's household surveys. Five hundred and seventy-five of them, from twenty-eight sources."
- "Nothing here required a new blood sample and nothing had to be licensed or purchased. Three of the sources need a free account. And when we cross a border, the map runs on just ninety-four of those layers — the climate and soil ones."

**Numbers and sources.** Ranks recomputed from `results/tables/policy_deck/viz/oof_child_iron.csv`, Ghana rows, respondent-weighted Admin-1 jackknife (Binduri 3/1/13 under both weightings; **Talensi is 10 under the weighted jackknife used everywhere else in this project and 12 unweighted — use 10**; Atiwa East survey rank 11 under min-ranking, 12 under average-ranking because three districts tie at 0.40 — use 11 and be ready to say why). 575 columns / 28 sources — `data/covariates/harmonized/predictors_admin2_shared_metadata.csv`, `quantities.json n_pred`, `n_src`. 94 deployed columns — `results/tables/policy_deck/civ_candidates_summary.csv` `n_common`. Access wording from `scratchpad/phase1_dataset.md` F9.

**Forbidden here.** "None of it is a blood test" (the modelled anaemia surfaces are fitted to earlier DHS haemoglobin — say **"no new blood sample"**). "575 layers built the map you are looking at" (deployment is 94 or 125). The stale chain "614 assembled / 525 / 464".

**Cut if short:** build 4's drop-one numbers (keep the four thumbnails and the one line).
**split_if_relaxed:** see B1, 4a / 4b. In 4b, add the **cost-to-acquire ladder**: leave-one-country-out ranking accuracy by what a ministry must do to obtain the data — *download today, no account* **0.316** (19 of 22 cells beat chance) · *+ public household-survey microdata, free account* **0.302** (18/22) · *+ DHS, data-use agreement* **0.281** (17/22), with a dashed chance line at **0.079**. Sources: `benchmarks_v2_summary_open.csv`, `benchmarks_v2_summary.csv`, `benchmarks_v2_summary_withdhs.csv`, `transport_null_calibration.csv`. This is the single most original chart available and it does not exist anywhere in the project. One line with it: *"the blocks that cost us weeks and registrations add one-thousandth inside a country — 0.400 with DHS, 0.399 without — and make cross-border prediction slightly worse."*

---

### Section: Modelling approach and what won

---

### Slide 5 — Is the new map any good? (Andrew, 2.75 min) — **three builds**

**Headline (claim form):** *The order is right more often than the region's average. The level is not.*

**Options for the top half**

- **Option A (recommended) — the odds ruler.** One horizontal axis, 50% to 75%, four ticks: *coin toss 50 · the survey's own regional averages 56 · model 60 · model in its six best country-nutrient combinations 70*. Greyed at the chance end. Replaces the current deck's seven separate restatements of 0.40. **Add a fifth tick: a covariate-free neighbour map, 59.** *Trade-off:* the fifth tick shrinks the apparent margin, and it is the first question a methods-literate listener asks; owning it is worth more than the two points it costs. The clause that makes it a win rather than a loss is in the talking points below.
- **Option B — four paired bars** (pairs, worst-third precision, flagged-fifth prevalence, level), each with its own comparator. *Trade-off:* more information, but it puts four metrics in a policy room in 40 seconds and the room retains one.
- **Option C — the two-districts pictogram** (ten pairs of icons, six lit for the model, five and a half for the regional average). *Trade-off:* the most intuitive object in any proposal; but half a lit icon reads as noise, and a 4-point gap drawn as ten icons lands as "these are the same".

**Options for the bottom half**

- **Option A (recommended) — the "we deleted Ghana" twin panel.** Left: Ghana's 75 districts, survey prevalence on x, held-out prediction on y (in-fill, ρ = 0.50; hollow markers = single-cluster districts; marker area ∝ respondents). Right: identical axes, but the model has never seen a single Ghanaian record — trained on The Gambia, Sierra Leone and Malawi only (ρ = 0.57). Caption spanning both: **"This is the Côte d'Ivoire experiment, run where we have the answer key."**
- **Option B — the per-country grade strip** (four tiles: Gambia 0.60, Ghana 0.35, Malawi 0.20, Sierra Leone 0.16 when each is held out entirely). *Trade-off:* it is the single most decision-relevant breakdown and pre-empts "will it work in my country"; but it is four numbers with no picture and it ends the slide on heterogeneity rather than on the demonstration.
- **Option C (recommended compromise) — both**, the twin panel large and the four country tiles as a small strip beneath it. This is what the 2.75 minutes buys.

**Middle of the slide, non-negotiable, 28 pt, printed not spoken:**

> **ORDER — yes.** Pick any two districts: the model says which is worse off **60%** of the time; the survey's own regional averages **56%**; a coin toss **50%**.
> **LEVEL — no.** Asked for a district's percentage, the model is off by about **14 points**; handing every district the national figure is off by **11.5**. It is a ranking, not a measurement.

**Talking points (Andrew) — a two-and-a-half-minute script.**
- "We tested it the hard way. We hid one district in five and predicted them, ten times over. Then whole regions. Then whole countries."
- "'Held out' means the model never saw that district's blood results while it was being fitted — so every number here is a prediction, not a fit."
- "And we did not test against zero. We tested against the number a programme officer actually has on her desk: her district's regional average, computed with her own district left out."
- "Here is the summary in two sentences, and they are printed so I cannot skip them."
- "On the order of districts, we beat that regional average: six district pairs in ten in the survey's order, against five and a half, against five for a coin. In our six strongest country-nutrient combinations, seven in ten."
- "On the level — the actual percentage in a district — we are worse than simply writing the national number on every district. So we do not put a percentage on a district without a survey behind it. An order, not a number."
- "One honest tick on that ruler: a map made only from neighbouring districts, with no satellite data at all, scores 59. Inside a surveyed country, geography gets you most of the way. What the covariates buy you is the border — and there is no neighbour map in a country with no survey to smooth."
- "Which brings the harder test. We deleted Ghana's survey entirely and trained on the other three countries. The model had never seen a Ghanaian blood sample, and it ordered Ghana's districts at 0.57 — about the same as when it had Ghana's own data."
- "That is the Côte d'Ivoire experiment, run where we can check the answer. And the spread matters: held out entirely, The Gambia scores 0.60, Ghana 0.35, Malawi 0.20, Sierra Leone 0.16. Chance is 0.08."
- "Sierra Leone is the honest failure. Fourteen districts of half a million people each — too few to run the district test at all, and below chance when we transport into it. We can tell you in advance which countries look like that one."

**Numbers and sources.** Pairs 0.598 / 0.555 / 0.591 / 0.50 — `results/tables/protocol_v2/nce_targeting_summary.csv`, `nce_targeting_metrics.csv` (**say "eighteen combinations", never 24**: all six Sierra Leone in-fill cells are NaN). Strong-cell pairs **0.699 exact**, regional 0.639 — `accuracy_statements.md` #2 (**never** the deck's 0.722 asin approximation). Level 14.25 vs 11.49 pp, index worse in **13 of 18** cells — `benchmarks_v2_summary.csv` rows 14/21, `benchmarks_v2_cells.csv`. Ghana child iron in-fill prevalence 0.502 / level 0.529 vs regional 0.498 / 0.356; leave-one-country-out 0.568 prev / 0.554 level — `benchmarks_v2_cells.csv`. Transport by held-out country 0.598 / 0.350 / 0.198 / 0.160, null 0.079 — same file, `transport_null_calibration.csv`.

**One clause on method, spoken, no slide element** (this is how the PI's "modelling approach / which models won" section is paid for at eight slides): *"The method is deliberately plain — rank every layer within the country, summarise each family of layers into a few scores, weight them by how well they tracked deficiency in the training districts, add them up. Nothing is tuned inside a fit. We ran it against fourteen alternatives — penalised regressions, random forests, boosting, geostatistical smoothers and four ways of ensembling them — on identical held-out folds, and none of them beat it by more than 0.03. With 14 to 87 districts per country the simplest estimator wins; that is a sample-size finding, not a law."* Sources: `SANDBOX_LOG_2026-09.md` SL-06 (1,304 fits), `results/tables/protocol_v2/sl_selection_sl.csv`. **Do not say "the SuperLearner chose the index"** — the discrete SuperLearner picks it in 5-9% of fits under squared-error loss and 16-30% under rank loss.

**Cut if short:** the strong-cell tick on the ruler (it is the post-hoc number anyway). **Never cut the LEVEL sentence.**
**split_if_relaxed:** see B1, 5a / 5b.

---

### Section: Using the models alone and with surveys to find high-risk areas — part 1, the product

---

### Slide 6 — A country with no survey (Andrew, 2.25 min)

**Headline (claim form):** *Côte d'Ivoire has no survey. Here is its ranking, and how firmly we hold it.*

**Options for the flagship panel**

- **Option A (recommended) — women's vitamin A.** Worst five: **Tchologo, Bounkani, Bagoué, Poro, Folon**. This is the one outcome with external corroboration that survives multiplicity correction and a latitude control: Spearman **0.495** with the WFP/MIMI (Tang et al., *Nature Food* 2026) vitamin A map over all **33** districts, BH q = 0.043 in a pre-specified 10-test family, **0.479** after partialling out centroid latitude. *Trade-off:* iron is the nutrient a programme room cares most about; but for CIV child iron the same external check gives 0.178 (not significant) and women's iron 0.053, and `docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md` says in terms *"do not present iron transport as externally supported"*. Putting a vitamin A corroboration next to an iron map is the trap three of five proposals fell into.
- **Option B — child iron as the flagship, external check moved to its own spoken sentence with "for vitamin A" said aloud.** *Trade-off:* keeps the nutrient the room wants; costs 15 seconds of caveat and invites the mismatch question anyway.
- **Option C — show both nutrients side by side** and say "the two maps are 0.93 correlated with each other — this is largely one map of environmental deprivation, and we say so". *Trade-off:* the most honest, and it pre-empts the sharpest criticism in `TWO_READINGS_2026-09g.md` Reading B; but it spends the slide on a limitation rather than on the product.

**Visual (Option A).** Left, large: CIV ranking map, 33 districts, women's vitamin A, from `results/figures/policy_deck/fig10_civ_ranking.png` (regenerate for this outcome). Right: a **firmness strip** — 33 dots on one horizontal rank axis with 90% rank-interval whiskers; **10 districts solid** (firmly in the worst third across refits), **21 hollow** (firmly out), 2 open circles in the middle. Legend, written as behaviour and reused on every map in the deck and dashboard: **"deep = the model does not change its mind when retrained · pale = it would move with different training data."** Corner box: the external check. Small thumbnail bottom-right: the dashboard.

**Talking points (Andrew) — a two-minute script.**
- "Côte d'Ivoire is the case this project exists for: all the public layers, no biomarker survey."
- "Trained on the four surveyed countries, the model ranks its 33 districts for six nutrients. The northern savanna belt comes out worst — Tchologo, Bounkani, Bagoué, Poro, Folon."
- "Two independently chosen versions of the model agree at 0.88 to 0.97 and share 83% of their worst fifth, so this is not a knife-edge."
- "And there is one independent map for Côte d'Ivoire: the WFP and MIMI vitamin A work published this year, built from a household budget survey, nothing to do with us. It puts the same belt worst — rank agreement 0.50 across all 33 districts, and it survives adjusting for latitude. For iron there is no independent map, and where we can check it, iron does not corroborate."
- "How much to trust it. With a whole country held out, a district we call worst-third really is in the worst third 52% of the time, against 33% by chance. When we name our five worst, 62% of those calls land."
- "The shading on the right is important and easy to misread. Deep means the model does not change its mind when we retrain it on different districts. It does not mean we have verified that district. Ten districts are firmly in the worst third, twenty-one firmly out, and two we will not call."
- "So: a shortlist, not a measurement. And it exists today, for a country with no survey, from data anyone can download."

**Numbers and sources.** 33 districts, 6 outcomes — `quantities.json civ_n`, `civ_outcomes`. Firm counts (women's vitamin A, climate-soil candidate): p80 = 10, p20 = 21 — `results/tables/policy_deck/civ_candidates_summary.csv` (**for child iron it is 10 / 20; match the counts to the outcome you map**). Median rank width 4 of 33 for this outcome, 5 for child iron — same file. Candidate agreement 0.876-0.972, worst-fifth overlap 83.3% — `civ_candidates_agreement.csv`. Worst-third call 52.4%, chance 33.0%, top-five 61.9% — `results/tables/policy_deck/viz/rank_interval_districts.csv`, `quantities.json hit3`, `hit3_base`, `hit5`. External check — `docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md`, `results/tables/tang_lsff_prespecified.csv`. CIV database 212 columns, deployed 94 (climate+soil) or 125 (five-domain) — `civ_candidates_summary.csv`.

**Forbidden here.** "0.45 over 32 districts" (it is 0.495 over 33 for women's vitamin A; 0.395 for children). "Permutation p = 0.01" (there is no permutation test in that table; p = 0.004, q = 0.043, family of 10). "Pale = survey here first" — `scratchpad/phase1_survey_design.md` §1b tested it and the undecided districts are *not* where the model is most wrong.

**Cut if short:** the dashboard thumbnail and the candidate-agreement clause.
**split_if_relaxed:** 6a as above; 6b = **"what it is looking at"** — five or six plain-language Ghana layer maps with their direction (+ = more deficiency), footer *"these weights say where deficiency is; they are not causes and not things to change."* Build from `results/tables/protocol_v2/index_importance_columns.csv` (Ghana / child_iron / level: the largest column shares are IHME overweight and wasting prevalence, ruminant livestock share, malaria parasite rate and mortality, relative staple and vegetable prices, distance to coast) and the drop-one domain ablation (`domain_ablation_loco_summary.csv`: climate +0.025, soil +0.025, anaemia surfaces +0.020). **Before using `results/figures/policy_deck/fig3_top_predictors.png`, rebuild it**: it currently prints raw variable codes ("mics salt iodised 15ppm") and colours MICS and budget-survey columns as remotely sensed. Two honest clauses to pair with it: *"about a third of the Ghana weight sits in satellite embedding dimensions we cannot name in words"*, and *"hold urbanicity constant and the cross-border result barely moves — 0.378 to 0.364 — so this is not just a map of cities"* (`urbanicity_conditioning.csv`).

---

### Section: Using the models alone and with surveys — part 2, the decision

---

### Slide 7 — Which decisions is this good enough for? (Andrew, 2.25 min) — **rows revealed one at a time**

**Headline (claim form):** *Good enough to order districts and to place survey clusters; not good enough to count cases or to judge a programme*

**Options**

- **Option A (recommended) — six rows, three columns, numbers spoken not printed.** *Decision · What you use today · Verdict* with a large tick / half-tick / cross. At 16:9 in a plenary hall, a four-column grid of sentence-length cells is ~180 words and unreadable from row 15. Cap each cell at five words and put the numbers in the mouth.
- **Option B — four columns including "how often right vs today vs chance"** (the decision-first original). *Trade-off:* far more defensible when photographed, which is the point of the slide; but it is the densest object in the talk.
- **Option C — a single Ghana case card**: one named district, what each decision would have been under today's information versus the model, truth revealed. *Trade-off:* the most memorable and most policy-native; but an n-of-1 card invites the charge that the district was chosen after seeing the answer, and the two "No" rows do not fit on it.

**Visual (Option A).**

| Decision | What you use today | Verdict |
|---|---|---|
| Which districts to visit, screen or enrol first | The region's average | **Yes** |
| Where the next survey puts its clusters | Convenience, equal spread | **Yes** |
| Sequencing a roll-out across districts | Region by region | **Yes** |
| A planning percentage for a district | The national number | **Only with a measured national sample** |
| How many deficient people a campaign reaches | Population × the national rate | **No** |
| Whether a programme worked | A repeat survey | **No** |

Caption under the grid: *"Three yeses, one partly, two noes."*

**Talking points (Andrew) — a two-minute script.**
- "This is the slide to photograph. Three yeses, one partly, two noes."
- "Yes to an order: which districts to reach first, where the next survey should sample, how to sequence a roll-out. All three need a ranking, and the ranking is the thing that works — sixty per cent of pairs against fifty-six for the regional averages you use today."
- "Partly on a planning percentage. The ranking has no level of its own. Measure one national number — about five per cent of a survey's sample, roughly forty blood draws for one nutrient — hang the ranking on it, and district figures land within about nine to ten points. Spend the same five per cent on a regional survey and you get about twelve."
- "And that anchor has to be measured, not modelled. We tried predicting a country's national level from world indicators across sixty-nine countries. It ranks countries reasonably. But it missed Sierra Leone by twenty-nine points and Malawi by twenty-two, and every district inherits that miss."
- "Two flat noes. First, cases. If your question is how many deficient people a campaign reaches, our worst fifth holds twenty-two per cent of them; the regional averages twenty; a perfect map forty-eight. In Ghana, simply going to the most populous fifth of districts reaches forty-six per cent. Burden sits where the people are. Rank by rate times population — you already have the population map — and read ours as the rate."
- "Second, evaluation. This is one point in time, and we are nearly blind to programme reach: in the whole dataset there are three district-level columns on supplementation and fortification. We cannot tell you whether a district is high-risk because it is unreached or high-risk despite being reached."
- "And one thing this does not do: it does not let you field a smaller survey. We tested that symmetrically — shrink the survey and our error degrades faster than the survey's. It tells the survey where to go; it does not reduce the sample you need."

**Numbers and sources.** Pairs 0.598 / 0.555 — `nce_targeting_summary.csv`. Anchored design: median district MAE **9.3 pp** at a 5% national anchor (mean over cell means 10.05), regional survey at the same 5% **11.7 pp**, district survey **17.3 pp**, district survey matches the anchored design at **25%** of a full sample (per-cell crossover median 40%) — `results/tables/protocol_v2/anchor_and_rank.csv`, `anchor_and_rank_summary.csv`, `anchor_and_rank_crossover.csv`; 5% = 24-71 respondents, median ≈ 42 (`accuracy_statements.md` #17). VMNIS: Sierra Leone predicted 40.9 vs true 12.2 (28.7 pp), Malawi 31.2 vs 9.2 (22.0 pp), model MAE 11.75 pp vs a no-information 16.23 over 108 surveys in 69 countries — `results/tables/national_composition_levels.csv`, `national_vmnis_loco.csv`; a national sample of ~1,000 children measures the level to ±2.4-4.2 pp at 95% (`national_estimates_all.csv` `obs_se` 1.2-2.1 pp). Burden 21.9% / 20.3% / 47.5% — `nce_targeting_summary.csv`; **a truly random fifth captures about 20% by arithmetic — do not print the table's 13.6% null arm as "chance"**, it is the training mean with ties. Most-populous fifth, Ghana child iron 46% — `viz/oof_child_iron.csv` `cap_pop`; across the six iron cells 44-69% (`accuracy_statements.md` #11). Three district-varying programme columns in the headline no-DHS set, ten counting the DHS block, 35 in the domain of which 25 are national constants — `predictors_admin2_shared_metadata.csv`. Symmetric survey-size test — `survey_size_symmetric_summary.csv`: survey +3.33 pp, spatial smoother +3.36 pp, **the index +6.66 pp** from full sample to 15%. **This table is dated 3 September and predates the 15-September outcome reconciliation; either re-run it or say the claim qualitatively ("it does not let you run a smaller survey") without the numbers.**

**Honest concession to have ready, because Ghana is the worked example.** In the Ghana child-iron cell specifically the regional-average list beats us on burden (34.6% vs 22.4%) and on the prevalence of the flagged fifth (45.8% vs 36.6%), and the pipeline's own pair-concordance for that cell is a tie (0.683 index vs 0.680 regional). The win there is the worst-third list (72% vs 60%) and the 185 districts the regional rule cannot reach at all. Say "at the third, not at the tenth" and never slice finer than a third in public.

**Cut if short:** the symmetric survey-size clause. **Never cut the two "No" rows.**
**split_if_relaxed:** see B1, 7a / 7b; in 7b use the **rates-versus-cases twin Ghana map** (predicted rate worst-fifth beside predicted cases worst-fifth, lines joining the districts in both lists) — it turns the burden "no" into an instruction rather than a confession.

---

### Section: Next steps and the ask

---

### Slide 8 — What limits the map, and what we are asking for (Andrew, 1.5 min + 0.25 close)

**Headline (claim form):** *What limits the map is how the surveys were designed — and every survey added improves everyone's map*

**Options**

- **Option A (recommended) — the ceiling staircase + the pooling curve + three asks.** Left: a step chart, x = clusters per district (1, 2, 3, 5), y = *"the best score any model could get against this survey's own district map"*, steps at **0.54 → 0.66 → 0.73 → 0.81**, with each country pinned at its real design (Sierra Leone 4.3, Gambia 2.3, Ghana 1.2, Malawi 1.2) and a hollow marker for what the model actually achieved (Gambia 0.54, Ghana 0.29, Malawi 0.20). Right: the pooling curve, three measured points and a dashed fourth labelled **"your survey here"** with no value attached, and three boxed asks.
- **Option B — the asks alone, large type, spoken over the decision ledger** (buying a slide back against the budget). *Trade-off:* asks are remembered when they are the last thing on screen; and the staircase is the only recommendation in the talk that this particular room can act on.
- **Option C — close on the checklist** ("what we can and cannot do yet", yes/partly/no). *Trade-off:* the most honest closer, and reviewers love it; but the checklist content already lives in slide 7's ledger two slides earlier, and ending a funder talk on limitations leaves no ask on screen.

**Talking points (Andrew) — a ninety-second script, then the close.**
- "Before you judge any district model, ask what the survey can see. Across our four surveys, three districts in four rest on a single cluster — one community standing for a whole district. In Ghana and Malawi it is more than eight in ten."
- "We measured what that costs, using the surveys' own variance structure. With one cluster per district, a model that knew each district's true prevalence exactly would still only score about 0.54 against the survey's district numbers. That is the bar, not one."
- "And the bar moves with the design. Two clusters per district takes it to 0.66; three to 0.73. The Gambia already sits near two, and The Gambia is where both the survey's own district map and ours are best."
- "So the clearest recommendation we have is a survey-design one, and it raises the bar for your own district estimates as much as for ours: put two or three clusters in a district rather than one district in two more places."
- "The second thing is pooling. With one training survey a country we have never seen ranks at 0.20; with three, 0.30. It has not flattened — though I should say the gain so far is in The Gambia and Ghana, and Malawi and Sierra Leone have not improved. A survey is a public good: yours improves everybody's map."
- "Three asks. If your country has a micronutrient survey — recent, planned or in a drawer — talk to us. We want two implementation partners willing to take one real targeting decision twice, once with the regional averages and once with our map, and let us publish what changed. And to survey designers: two clusters in a district."
- **Close (0.25 min, spoken over the slide):** "We have written two candidate models down in advance, dated the third of September, so the fifth survey is a real test and not another description of these four countries. And the question we would like this room to answer is one we cannot: how accurate does a district ranking have to be before you would use it? We can tell you what we get. Only you can tell us what is good enough. An order, not a number — come and find us and we will show you your country's map on a laptop."

**Numbers and sources.** Single-cluster share: **74% over all 206 surveyed districts** (Gambia 57%, Ghana 83%, Malawi 85%, Sierra Leone 0%) — `variance_components_ceiling.csv`. **`quantities.json single_share` = 0.797 excludes Sierra Leone; if you use 80%, say "of the districts we could run the in-fill test in".** Clusters per district 2.33 / 1.20 / 1.18 / 4.29 — same file. Staircase 0.539 / 0.663 / 0.731 / 0.806 (prevalence target, 17 non-degenerate cells; level target 0.582 / 0.698 / 0.759 / 0.825) — `scratchpad/ceiling_by_clusters_projection.csv`, recomputed and confirmed. Achieved by country (prevalence) 0.541 / 0.292 / 0.201 — same file. Pooling curve 0.198 / 0.259 / 0.302 (44, 54 and 16 fits respectively; on the fixed 16-cell set 0.184 / 0.253 / 0.302; climate+soil 0.315 / 0.368 / 0.403) — `results/tables/protocol_v2/training_curve_climate_soil.csv`, recomputed. Pre-registration, 7 predictions, 3 September 2026 — `docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md`.

**Honest caveats to carry.** The staircase is a projection from measured design effects on the cells whose variance model converged (26 of 44 Admin-2 fits are singular); say "projected from these surveys' own variance structure". Sierra Leone has the most clusters per district (4.3) and the weakest results — the answer is that 14 districts of ~500,000 people is the wrong *unit*, not the wrong design, and the recommendation is clusters per district **at a sensible number of districts**. Someone will raise this; have it ready.

**Cut if short:** the pooling curve becomes one spoken sentence, and ask 3 drops.
**split_if_relaxed:** see B1, 8a / 8b.

---

## B3. Alternative spines

**1. DECISION-FIRST** (`proposal_decision-first.md`; mean 39.0, effectively tied for best). Opens on "Ghana has 260 districts; this year you can reach about a third — which third?" and hangs every result on a row of the decision ledger. *Prefer it when* the room is known to be programme managers and district planners rather than donors, or when the PI wants the ledger on screen for a full 2.5 minutes. *It does better:* the opening frame is the best in the set; the two refusals (level, cases) are structural rather than appended; "three yeses, one partly, two noes" is a better repeatable line than most theses. *It loses:* the product arrives at slide 8 of 9, so the room spends eight minutes hearing about a map it has not seen; it never says what the model is or why it works; and it runs nine slides against a budget of eight.

**2. ONE COUNTRY'S STORY — Ghana end to end** (`proposal_one-country-story.md`; mean 37.3). Ghana's 2017 survey, the free layers over all 260 districts, two hold-out tests in Ghana (district-out 0.50, Ghana deleted entirely 0.57), the full-country map with grades, then one slide generalising to four countries plus Côte d'Ivoire. *Prefer it when* the PI wants maximum narrative memorability in Accra and is confident the child-iron rebuild will be done. *It does better:* the single strongest device in the whole set ("we deleted Ghana's survey and it still put the right districts on top"); one continuous thread instead of six; the grade strip. *It loses:* five of Andrew's minutes rest on one cell of 24, against a programme title that promises a multi-country framework; it stakes everything on the country whose external check is explicitly labelled "uninformative and latitude-dominated" (`EXTERNAL_CHECK_TANG_LSFF_2026.md` §3); and in Ghana specifically a covariate-free neighbour smoother essentially ties the model (0.526 vs 0.529 on the level), which is the first Q&A question and has no good answer inside a Ghana-only frame.

**3. THE DATA-STRATEGY TALK** (`proposal_data-strategy.md`; mean 36.0). The talk is an argument about what data to buy next, settled by a four-country experiment: the free layers travel, the expensive ones do not, our surveys set the ceiling, each survey improves everyone's map. *Prefer it when* the audience is donors and survey funders rather than programme officers — a WHO/UNICEF/Gates-heavy session. *It does better:* the most original argument in the set and the only one that turns the project's modest headline margin into a large, clean, actionable finding; the cost-to-acquire ladder and the ceiling staircase are its inventions and both are grafted into the recommended spine. *It loses:* it addresses a budget line that no national nutrition directorate holds; it never shows what a ministry would be handed; and its title ("the data that limit us are our surveys") is delivered in Accra immediately after Sonja thanks the four survey teams.

**4. HONEST SCIENCE AS THE ASSET** (`proposal_honest-science.md`; mean 34.3). Lead with the test design, then what survived, then a ledger of the project's own withdrawn claims (three that flattered it, one that did not), then the product with its uncertainty described as stability, then the ask. *Prefer it when* the session is a methods session, or when the PI expects the sharpest statisticians in the room to be the ones he must convince. *It does better:* the withdrawal ledger is the single most differentiating image available and nothing at the Forum will look like it; the comparator discipline ("what we had to beat") is the best in the set; it gets more disputed numbers right than any other proposal. *It loses:* five of nine slides lead with or end on a negative, so a policy room leaves with "not ready yet"; the product arrives at slide 8 of 9; and it spends about 2.5 of Andrew's minutes on meta-honesty (ceiling thermometer, withdrawal ledger, calibration strip) that a programme officer cannot act on.

**5. A hybrid the judges implied but nobody proposed.** Open on the decision (DECISION-FIRST slide 4), deliver the product at minute 4 (BEFORE/AFTER/PRODUCT slide 4), evaluate once (slide 5), close on the data strategy (DATA-STRATEGY slides 7-8). That is the recommended spine above. Its distinctive risk is that it has no single authorial voice — mitigate by saying "an order, not a number" at slides 4, 5 and 8, which is the only thread that has to survive.

---

# Part 3. Communicating accuracy to a policy audience

## C1. Plain-language accuracy statements

Drawn from `scratchpad/accuracy_statements.md`; only A- and B-rated statements are listed. C-rated statements (levels, WHO bands, burden in absolute people) are excluded from the main body and appear in Part G. "Perfect" = the oracle arm where one exists.

| # | Statement (say it this way) | Model | Today (regional avg) | Chance | Perfect | Scope | Caveat | Rating |
|---|---|---|---|---|---|---|---|---|
| **1** | "Pick any two districts. The model says which is worse off **60%** of the time; the survey's own regional averages **56%**; a coin toss **50%**." | 59.8% | 55.5% | 50% | — | In-fill, 18 country-nutrient combinations (Gambia, Ghana, Malawi) | A covariate-free neighbour map scores **59.1%** — say so; the covariates earn their keep across borders, where no neighbour map exists | **A** — headline |
| **2** | "Per country: The Gambia **69%** against 60; Ghana **62%** against 60; Malawi **54%** against 50." | — | — | 50% | — | Same, by country | Ghana's margin is 2 points. Sierra Leone has no in-fill rows at all (14 districts) | **B** — the honest breakdown |
| **3** | "In the six combinations where it works best, **70%** of pairs, against **64%** for the regional averages." | 69.9% | 63.9% | 50% | — | 6 of 18 cells with in-fill level ρ ≥ 0.5 | Chosen after seeing the scores; three of the six are The Gambia. Say "in the best cases" out loud. **Use 0.699, never the deck's 0.722 approximation** | **B** |
| **4** | "In a country the model has never seen, it orders **53%** of pairs correctly; chance is 50%." | 52.6% | none exists | 50% | — | Leave-one-country-out, 21 cells | The 22-cell figure is 54.6% but is inflated by Malawi women's vitamin A (0.1% national prevalence, one deficient district, concordance 0.965). **Say 53%.** Per country: Gambia 62%, Ghana 58%, Malawi 54% (60% with the artefact), Sierra Leone **40% — below chance** | **B** |
| **5** | "A district we call worst-third in a country we have never surveyed really is in the worst third **52%** of the time; chance is **33%**. When we name our five worst, **62%** of those calls land." | 52.4% / 61.9% | none exists | 33% | — | 1,176 held-out districts, 22 cells, climate+soil | Remarkably flat by country: Gambia 56%, Ghana 52%, Malawi 52%, Sierra Leone 50% — a better fact than the correlation spread | **A** — the right caption for Côte d'Ivoire |
| **6** | "Of the districts the model names Ghana's worst third for child iron, **18 of 25** really are — 72%. The regional-average list gets **60%**. Chance is 33%." | 72% | 60% | 33% | 100% | Ghana child iron, in-fill, k = 25 of 75 | 60% is the respondent-weighted regional jackknife used everywhere in this project; the unweighted variant gives 68%. **Print the comparator definition on the slide.** At the worst *fifth* the regional list wins in this cell (67% vs 47%) — never slice finer than a third in public | **B** |
| **7** | "The fifth of districts the model flags runs at **32%** deficiency against **24%** nationally; a perfect map would find **41%**." | 32.2% | 30.1% | — | 41.3% | 18 cells, in-fill, population-weighted | The regional-average list's flagged fifth is 30.1%, so the margin is 2 points, not 8. **Quote the regional arm or do not quote the national one.** | **B** |
| **8** | "A district the model calls high-priority has on average **38.5%** iron deficiency, against **30%** nationally and **22%** in the districts it calls low-priority." | 38.5% | — | — | — | 6 iron cells, in-fill, unweighted district means | The regional-average list's top fifth is similar or higher in Ghana child iron (44% vs 36%). Use as a gradient statement, not as a win | **A** — the most intuitive sentence in the set |
| **9** | "When the model is sure — a district lands in its worst fifth in at least 80% of forty retrainings — it is in the survey's worst fifth **35%** of the time; when it is sure a district is not, **16%**. The base rate is 20%." | 34.9% / 15.8% | — | 20% | — | 1,266 district-outcome rows, 18 cells | **"Sure" means stable, not verified.** Nearly twice the base odds, and still wrong more often than right | **B** — the honest way to caption any priority map |
| **10** | "If a programme can reach a fifth of districts, the model's list reaches **22%** of the deficient people; the regional averages **20%**; a random fifth about **20%**; a perfect list **48%**." | 21.9% | 20.3% | ~20% | 47.5% | 18 cells, in-fill | **REGIONAL AVERAGE ESSENTIALLY TIES.** In Ghana child iron the regional list wins outright (35% vs 22%), and going to the most populous fifth reaches 46%. The table's null arm reads 13.6% — that is the training mean with ties, not a random fifth; **do not print 14% as chance** | **B** — say it as a "no" |
| **11** | "Asked for a district's percentage rather than its position, the model is off by about **14 points**; giving every district the national figure is off by **11.5**." | 14.25 pp | 11.77 pp | 11.49 pp (national) | — | 18 cells, in-fill, population-weighted | Worse than the national number in **13 of 18** cells. Shrunk properly it reaches 10.74 pp — better than both, and better than the national number in 9 of 18 cells — but its ranking then falls to 0.206, below the regional average's 0.220. **The calibrated arm is already in the committed tables; plotting it is a one-line change.** **Two settings, two products; say so if you quote the calibrated number** | **A** — the ✗ sentence |
| **12** | "One measured national number — about five per cent of a survey's sample, roughly forty blood draws — plus the ranking gives district figures within about **nine to ten points**. The same five per cent spent on a regional survey gives **twelve**; on a district survey, **seventeen**." | 9.3 pp | 11.7 / 17.3 pp | — | — | 22 cells, leave-one-country-out, climate+soil | A district survey matches it at 25% of a full sample (per-cell crossover median 40%), and beats it on burden at *every* size. The anchored design's error barely improves with more anchor (8.5 pp at 100%) because the ranking, not the anchor, is the error | **B** — the costable result |
| **13** | "With one cluster per district, a model that knew the truth exactly would still only score about **0.54** against the survey's own district numbers. Two clusters: **0.66**. Three: **0.73**." | — | — | — | — | 17 non-degenerate cells, prevalence target | A projection from these surveys' variance structure; 26 of 44 fits are singular. **Never quote a share-of-attainable percentage** (withdrawn, `CLAIMS_REGISTER` 9.3) | **B** — the ask, not a defence |

**Where the regional average wins, and must be conceded:** district prevalence (11.77 vs 14.25 pp), WHO vitamin A severity bands (66% vs 64% exact, 94% vs 91% within one), burden capture in the Ghana child-iron cell (35% vs 22%) and its worst-fifth precision (67% vs 47%), and the whole of Ghana child vitamin A and women's folate. The defensible line is: *"on the order it beats the regional average; on the number and on the cases it does not, and neither does the regional average beat one national figure."*

## C2. Six devices for saying "how good is it" without a correlation

| Device | Form | Numbers | When to use it |
|---|---|---|---|
| **Pairs of districts** | "Pick any two. We say which is worse off six times in ten; you get five and a half today; a coin gets five." | 59.8 / 55.5 / 50 | The default. It is the only device that carries a comparator, a chance line and a unit the room already has |
| **Worst-fifth precision** | "Of the ten districts we name, seven really are in the worst third." | Ghana child iron 18 of 25 at the third | When the audience is a programme with a fixed shortlist size. **State the cut at the third, never at the tenth** — the comparison flips at fine slices |
| **Severity gradient** | "The districts we send you to run at 38% deficiency; the ones we send you away from, 22%." | 38.5 / 30.0 / 22.0, six iron cells | Warmest available statement; converts a rank into a consequence without claiming a level |
| **The ceiling** | "Even a model that knew the truth would only score 0.54 against a survey with one cluster per district." | 0.54 → 0.66 → 0.73 at 1/2/3 clusters | Only as a *survey-design ask*, never as a defence. Delivered cold in a room of survey producers it reads as "the survey is bad, don't blame us" |
| **The decision checklist** | Three yeses, one partly, two noes, with the comparator in the second column | Slide 7 | The closer, and the thing the room photographs. Also the safest way to deliver the two negatives |
| **"Changes its mind" stability** | "Deep = the model does not change its mind when retrained. Pale = it would move. That is stability, not truth." | 166 firm-in districts at 35% true vs a 20% base rate; 90% rank intervals cover the truth 38% of the time | Any map with shading. **Never let a legend imply that deep = verified**; the honest gloss is the frequency, "when we call a district like this, it is in the survey's worst fifth about one time in three, against one in five by chance" |

**Two devices to avoid.** (a) Any share-of-attainable percentage ("59%", "two-thirds") — withdrawn, and the ratio moves from 48% to 83% depending on which ceiling and which target you pick. (b) Any transform of a Spearman into a percentage via `(1 + 2/π·asin ρ)/2` — it is what produced the deck's 0.722 where the exact pair count is 0.699, and it is why the national-track "71% of country pairs" should be dropped in favour of "ranks countries at 0.63".

---

# Part 4. Visualisation brainstorm

## D1. Visuals to keep from the 90-slide deck

| File | Why |
|---|---|
| `results/figures/policy_deck/fig10_civ_ranking.png` (top row only) | The deliverable made concrete: a real country with no survey, ranked, with a firmness panel. Regenerate for **women's vitamin A** if Option A on slide 6 is taken. Do **not** show the bottom row of `fig13_civ_candidates.png` — near-black P > 0.75 fill, unreadable at slide size |
| `results/figures/policy_deck/fig5_learning_curve.png` | Best next-steps figure in the project; 18-pt fonts, chance band drawn. Add the dashed "your survey here" fourth point |
| `results/figures/policy_deck/fig6_targeting.png` | The only figure that speaks in people. Keep for the appendix; the slide-7 "No" row is a sentence, not this chart |
| `results/figures/policy_deck/fig1_model_comparison.png` | Policy-grade already (≤6 categories, "needs a survey here" in place of missing bars). Appendix, answers "what did you actually fit" in one image |
| `docs/slides/img/ghana_sources_land.png`, `ghana_sources_health.png` | The right way to show "what the layers look like" in 20 seconds. Crop four tiles for slide 4's bottom strip; use one panel, not both |
| `docs/slides/img/goals_ghana_map.png` | Already built on child iron — the correct cell. Candidate for Sonja's slide 2 |
| `results/figures/policy_deck/fig11_worst_fifth_probability.png` | The "how sure" map with ground truth behind it; better than the exceedance slide for the appendix |

## D2. New visuals to build

Grouped by slide, each with its source table and the reason a policy audience reads it faster than what the deck has now.

**Slide 2 (Sonja)**
1. **The 260-cell coverage waffle.** 75 filled (62 half-tone single-cluster, 13 solid), 185 outline. Source: `viz/oof_child_iron.csv` (`n_psu`), `metadata/admin2_spine.csv`. *Faster because* the deck currently states coverage in three sentences across two slides; the waffle needs none, and it is reused on slide 8 as the ask.

**Slide 4**
2. **The three-district callout with leader lines**, including the deliberate miss. Source: `viz/oof_child_iron.csv`. *Faster because* it is the only named, single-district, checkable claim anywhere in the project; a room remembers one district and no correlations.
3. **The rebuilt Ghana child-iron trio** (survey 75 / held-out prediction / all 260). Requires a new `deploy_ghana_child_iron.csv`; `dashboard/data-raw/05_build_protocol_v2_bundles.R` already fits the deployment ranking for every cell. *Faster because* the existing `fig7_ghana_map.png` is child vitamin A, the cell where the regional average wins.
4. **(4b, relaxed) The cost-to-acquire ladder.** Three bars labelled by what a ministry must *do* to obtain the data (download / free account / data-use agreement), 0.316 / 0.302 / 0.281, chance line 0.079, with the cells-beating-chance count on each (19 / 18 / 17 of 22). Sources: `benchmarks_v2_summary_open.csv`, `benchmarks_v2_summary.csv`, `benchmarks_v2_summary_withdhs.csv`, `transport_null_calibration.csv`. *Faster because* the deck currently reports predictor tiers by domain name; a ministry reads acquisition burden, not taxonomy. **Do not mix taxonomies on one axis**: climate 0.335 and soil 0.327 are *domain*-only arms (`domain_ablation_loco_summary.csv`), while GEE 0.267 and IHME 0.206 are *source*-only arms (`source_ablation_loco_summary.csv`) and the climate layers live inside the GEE source.

**Slide 5**
5. **The odds ruler.** One horizontal axis, 50-75%, five ticks (coin 50 · regional 55.5 · neighbour map 59.1 · model 59.8 · six best cells 69.9), chance end greyed. Source: `nce_targeting_summary.csv`, `accuracy_statements.md` #2. *Faster because* it replaces the seven separate restatements of 0.40 on slides 30, 36-41 and 66 of the full talk.
6. **The "we deleted Ghana" twin scatter.** Identical axes; hollow markers = single-cluster districts, marker area ∝ respondents; one caption across both. Source: `viz/oof_child_iron.csv` for the left panel. **The right panel needs a table that does not exist** — no protocol-v2 script writes per-district leave-one-country-out predictions (`rank_interval_districts.csv` holds ranks, not predictions). Either commission that export or make the right panel rank-versus-rank and say so on the slide.
7. **The four-country transport strip.** Four tiles, Gambia 0.60 / Ghana 0.35 / Malawi 0.20 / Sierra Leone 0.16, chance line 0.08. Source: `benchmarks_v2_cells.csv` (estimand `country`). *Faster because* it answers "will it work in my country" with four numbers instead of one mean that hides a 0.16-to-0.60 range.

**Slide 6**
8. **The firmness strip.** 33 CIV districts as dots on one rank axis with 90% whiskers; solid = firmly worst-third (10), hollow = firmly out (21), open = undecided (2). Source: `civ_rank_uncertainty_all.csv`. *Faster because* it replaces the bivariate choropleth, whose legend nobody decodes in 20 seconds, and because the undecided band doubles as a survey-siting cue without claiming one.
9. **The travelling grade legend.** One legend card, reused on every map in the deck, the dashboard and the country briefs: *deep = does not change its mind when retrained · pale = would move · hatched = no survey cluster.* Source: the A/B/C convention already in the full talk at `.qmd` line 948. *Faster because* the room learns it once and it travels outside the talk.
10. **(6b, relaxed) Five plain-language driver maps** with + / − direction. Sources: `index_importance_columns.csv`, `domain_ablation_loco_summary.csv`. *Faster because* `fig3_top_predictors.png` currently prints raw variable codes and mis-colours MICS/HCES items as remotely sensed.

**Slide 7**
11. **The decision ledger** (six rows, three columns, ticks). *Faster because* it is the only artefact in the talk a programme officer can photograph and use in a budget meeting without a statistician.
12. **(7b, relaxed) Rates-versus-cases twin Ghana map** with the two ten-district lists and lines joining the few districts in both. Source: `viz/oof_child_iron.csv` × `dashboard/data/admin2_population.rds` (CSV already at `scratchpad/admin2_population.csv`; **normalise the "Sierra Leone" spelling**). *Faster because* it converts the burden negative into an instruction — rank by rate × population.

**Slide 8**
13. **The ceiling staircase.** Steps at 0.54 / 0.66 / 0.73 / 0.81 for 1/2/3/5 clusters per district, each country pinned at its real design with a hollow achieved marker. Source: `scratchpad/ceiling_by_clusters_projection.csv` (**not a committed project table — commit and re-run the projection script before the deck is built**). *Faster because* the deck's current ceiling slide is a static bar that reads as an excuse; this one is a costed lever.
14. **The pooling curve with "your survey here."** Three measured points and a dashed fourth with no value. Source: `training_curve_climate_soil.csv`. *Faster because* a visual ask beats a spoken one and it does not over-promise a number.
15. **The unit-definition footer strip** (country / what "district" means / surveyed units / clusters per unit / median respondents: Gambia 30 / 2.3 / 22 · Ghana 75 / 1.2 / 14 · Malawi 87 TAs inside 27 districts / 1.2 / 12 · Sierra Leone 14 / 4.3 / 36). Usable on any map slide. Sources: `variance_components_ceiling.csv`, `targets_v2.csv`. *Faster because* it pre-empts the single most likely correction from the floor.

## D3. Visuals to retire, and why

| Retire | Why |
|---|---|
| `fig3_top_predictors.png` **as it stands** | Six panels × five bars; raw codes leak ("mics salt iodised 15ppm", "u5 diarrhoea prev"); the two-colour legend labels MICS survey items "Remotely sensed environment". Rebuild as one outcome, five bars, plain names, three colour groups — or do not show it |
| `fig7_ghana_map.png` **as it stands** | Built on Ghana child vitamin A, the one Ghana cell where the regional average beats the index on both targets and on district error. Rebuild on child iron |
| `ensemble_orchestra.png`, `ensemble_vfold.png`, `ensemble_no_overfit.png`, `old_s10_ensemble_*` | The January ensemble cartoons, and slide 20's "the best team wins" caption, contradict the project's own finding that the untuned index won |
| `fig13_civ_candidates.png` bottom row | Near-black P > 0.75 fill on a 0.00-1.00 legend; unreadable at 24 pt |
| Slides 26-28, 46-47, 50, 57-59, 61, 78, 89 of the full talk | Correct answers to questions this room will not ask: population-weighted MAE by cell, RMSE/SD ratios, Spearman-vs-prevalence scatters, predicted-vs-observed clouds, rank-interval coverage, ablation deltas, meta-analytic z, weighting arms. All appendix |
| The three near-identical four-panel map slides (3, 31, 42) plus the three Ghana map slides (45-47) | Six slides competing for one visual slot |
| The Malawi B12 fish/anaemia/malaria forest (55-56) | The most interesting science in the deck and the most dangerous here: a non-statistician's take-home is "more malaria means better B12". Appendix, with the maps rather than the forest, and only if the speaker has 40 seconds |
| The "two readings" framing as a slide title | It invites the room to pick the pessimistic reading. Keep the discipline in the notes; on screen use slide 7's ledger form |

---

# Part 5. The national-level (VMNIS) slide: include or appendix

**Decision: appendix only, plus one spoken sentence inside slide 7's "planning percentage" row.**

**The one sentence (recommended wording, already in the slide-7 script):** *"And that anchor has to be measured, not modelled. We tried predicting a country's national level from world indicators across sixty-nine countries. It ranks countries reasonably, but it missed Sierra Leone by twenty-nine points and Malawi by twenty-two — and every district inherits that miss."*

**Why not a main-body slide, in order of force.**
1. It answers a between-country question this room is not asking; countries already have national estimates (VMNIS itself, GBD).
2. Its one usable result is negative, and the negative is already carried by the slide-7 sentence.
3. It is not re-audited: the tables are dated 31 August / 1 September from an April 2026 sub-pipeline, the composition uses the *pre-protocol-v2* elastic-net transport pattern, and it used the wrong survey years for The Gambia (2021 for a 2018 survey) and Malawi (2015 for 2016).
4. The panel does not contain any of the four project surveys. The child vitamin A panel is 515 of 528 rows serum retinol, median survey year **2001**, adjustment unrecorded in 290 of 528 rows — a different era and a different assay from the BRINDA-adjusted RBP the composition is scored against. That, not the covariates, is why Sierra Leone comes back at 41% against a true 12%.
5. The deck's *current* appendix slide titled "National-level track" does not show VMNIS at all — it tabulates `national_estimates_all.csv`, the within-survey person-level SuperLearner aggregate, whose near-zero errors are self-consistency. **Relabel or replace it regardless of this decision.**

**Exact appendix-slide content** (title in claim form): **"A model cannot replace the national blood sample: predicted national levels miss by 12 points, and every district inherits the miss."**
- Left: dumbbell chart, four countries, child vitamin A — survey national prevalence with its ±2 SE bar (Gambia 19.8 ± 4.3, Ghana 14.9 ± 2.5, Sierra Leone 12.2 ± 3.4, Malawi 9.2 ± 3.2), the no-covariate "other countries' average" (17.1-17.5), and the World-Bank-indicator model (20.3, 20.8, **40.9**, **31.2**).
- Right: four bars, district-map mean absolute error when the transported map is anchored to nothing **9.1**, the training countries' average **10.5**, the model's national level **19.1**, the country's own measured national level **8.6** pp (4 non-degenerate cells; the 8-cell versions are 5.8 / 8.2 / 12.7 / 5.6 with a +10.0 pp bias on the modelled arm).
- Footer: "WHO VMNIS, vitamin A in children: 108 national surveys, 69 countries, leave-one-country-out; 17 World Bank indicators; random forest. MAE 11.75 pp against 16.23 for a no-information guess. **Earlier sub-pipeline, not re-audited.**"
- Sources: `results/tables/national_vmnis_loco.csv`, `national_composition_levels.csv`, `national_composition_revised.csv` (`MAE_pp_excl`, `mean_signed_bias_excl`), `national_estimates_all.csv` (`obs_se`).

**Numbers not to use.** "71% of country pairs" (an asin approximation of ρ = 0.625 — say "ranks countries at 0.63" or compute the exact pair count). Any African-subset panel description ("30 country-years, 22 before 2005 at about 50%") — it is not reproducible from the audited material; the verified panel facts are 528 rows, 70 countries, 109 country-years, 515 serum retinol, median year 2001, era means 29.8% pre-2000 / 21.1% 2001-10 / 20.5% 2011+.

**Two alternatives, if the PI wants a national slide anyway.**
- **Alt 1 — "the measuring stick problem."** For women's vitamin A, between-method variance exceeds between-country variance (sd 2.10 vs 1.17 in logit; r_max 0.48 report-level vs 0.95 standardised). Plain language: *how* a survey measured explains more of the gap between two countries' reported numbers than *which* country it was. Policy hook: harmonised assays and BRINDA adjustment before any cross-country comparison — a live Forum theme, and it is the same mechanism that defeats level transport at district level. Caveat: 90 surveys, ~20 countries, "method" is a post-hoc string class. One spoken sentence at most.
- **Alt 2 — "rank countries to decide where to survey next."** ρ = 0.63 over 69 countries, MAE 28% below the null for preschool vitamin A; folate ≈ 0.5; nothing for women's vitamin A or B12. Weakest option: countries already have national numbers, the panel's median year is 2001, and the project's deliverable is subnational.

---

# Part 6. Next steps and augmenting nutrition policy

## 6.0 The conclusion point to lead with: more surveys make every map better, and we are asking for them

*Added 18 September 2026 at the PI's direction. This is the project's own result, it is the only next-step claim backed by a measured curve rather than a plan, and it should be the last thing the room hears.*

**The claim.** What limits this method today is not the predictor data — it is the number of biomarker surveys it can learn from. With one training survey a country the model has never seen is ranked at **0.20**; with three, **0.30**; about **+0.05 per country added**, and the curve has not flattened at three. On the pre-registered climate-and-soil set the same curve runs 0.31 → 0.37 → 0.40, and a single training survey already gives a positive ranking in 95% of fits.

**Why it is not a generic "more data would help" plea.** Two independent results say the predictors are saturated and the surveys are not:

- Eleven candidate predictor blocks were scored under one harness this year and **none moved in-fill or transport by more than 0.01 at the median**.
- The **open, permission-free tier transports best** (0.316) — better than the full headline set (0.302) and better than the set with DHS microdata added (0.281). More layers is demonstrably not the lever; more countries is.

**The refinement that makes the ask specific — and this is new.** It is not only how many countries, but which. Trained on **Ghana alone**, the model ranks The Gambia's districts at **0.61**; trained on **The Gambia alone**, it reaches Ghana at 0.30 and Malawi at 0.09. Ghana's 260 districts span the widest environmental gradient in the panel (NDVI SD 0.16 against 0.04–0.09 elsewhere; GPP SD 1.24 against 0.36–0.91). **A training country that contains the whole gradient teaches the whole gradient** — and whether a candidate country extends the range can be checked from free rasters before a single blood sample is drawn. That turns "please send surveys" into a specific, testable request.

**Three qualifiers to carry, because someone will find them.**

1. The gain is uneven. As held-out countries, The Gambia improves 0.40 → 0.60 and Ghana 0.30 → 0.39 as training surveys are added; **Malawi (0.16) and Sierra Leone (0.06) do not improve at all.** "Every survey improves every map" is a statement about the average.
2. Three points on a curve from four countries on one continent. The slope is a description, not a law, and it is on the biomarker-level target; **levels do not improve with more training countries**, because the obstacle there is cross-survey measurement offsets, not sample size.
3. Say the falsifier. If the pre-registered predictions fail in Ethiopia or Pakistan, the transport result was a property of West and Southern Africa rather than of the proxies.

**Suggested closing wording (about 35 seconds).**

> "The last thing I want to leave you with is what actually limits this. It is not the satellite data — we added eleven new blocks of it this year and not one moved the answer by more than a hundredth. It is the number of biomarker surveys we can learn from. One survey gets a new country to 0.20; three get it to 0.30; the line has not flattened. And it matters which country: trained on Ghana alone we rank The Gambia's districts at 0.61, because Ghana spans the whole gradient from the Sahel to the forest. So the ask is simple. If your country has a micronutrient survey — recent, planned, or sitting in a drawer — put it in the pool. It improves the map of every other country in it, and yours gets a ranking for the districts your survey never reached."

**Where it goes.** Slide 8 in the recommended spine, as the second half, with the pooling curve and the "your survey here" dashed point. If the talk is cut to seven slides, this is the one closing point that survives; the pre-registration and the cluster-per-district ask can both be spoken over it.

## F1. Three alternative framings of the closing

**Framing A (recommended) — "The survey is the bottleneck, and every survey is a public good."**
Ceiling staircase, pooling curve, three asks. *Why:* it is the only closing in which the audience — survey teams, funders, ministries — holds the lever. It converts the project's most defensive result (the reliability ceiling) into its most actionable one, and it visibly ends the talk asking for *more* measurement, which is the necessary antidote to "so we can stop surveying". *Risk:* it is a criticism of designs made by people in row two. Mitigation: attribute the one-cluster design to standard region-level power calculations, which is what it is, and put the benefit to them first — *"two clusters gives your own district estimate a number worth having."*

**Framing B — "Three uses now, three later."** The decision ledger repeated as the closing image, with the three yeses on the left and the bridge on the right (two pre-registered indices scored on the fifth survey; a retrospective test of model-guided cluster placement on the four surveys we already have; two country pilots; two regional workshops). *Why:* it is the most directly actionable and the ledger is the slide people photograph. *Risk:* the room has seen it two slides earlier; repetition at 1.5 minutes reads as running out of material.

**Framing C — "Come and get your country's map."** Product-forward: the dashboard, the country brief, an offer to run any country in three weeks with a grade attached. *Why:* it guarantees follow-up and it is what a funder room responds to. *Risk:* it writes a cheque the evidence cannot clear — transport ranges 0.60 to 0.16 and Côte d'Ivoire has no ground truth. **If used, the refusal condition must be on the slide:** *"and we will tell you when the grade is too low to act on — that happened in one of our four countries."*

## F2. The where / when / how-to-survey evidence, in plain language, with its honest limits

**WHERE — two answers, one strong and one not yet.**
- *Strong (design):* put two or three clusters in a district rather than one district in two more places. The Gambia's 2.33 clusters per district supports an attainable bar of ~0.70 and achieved 0.41-0.70 on all four of its country-outcome combinations; Ghana and Malawi at 1.2 clusters sit at a bar of ~0.54 and ~0.56 and achieve about 0.34. Projected: one cluster 0.54, two 0.66, three 0.73, five 0.81. *Limit:* this is a correlation across four countries plus a variance-component projection, not a reallocation experiment; nobody has moved a fixed budget from districts to clusters and re-scored. And Sierra Leone has the most clusters per district (4.3) and the worst results, because 14 units of half a million people is the wrong unit — the rule needs "and enough districts with real spread between them" attached.
- *Not yet (model-guided placement):* the maps show where the model is **sure**, not where it is **right**. Tested directly: binned by worst-fifth firmness, the undecided districts' median rank error is 25% of the list against 21% for firm-out and **27% for firm-in** — the stability signal does not find the model's errors. So *"pale = survey here first"* is a design proposal, not a validated rule, and must be said that way. What can be done before anyone fields anything: a retrospective sentinel-selection experiment on the four surveys — choose k = 3, 5, 10 districts to measure directly by four rules (random, most populous, model-undecided, largest model-vs-regional disagreement), combine with the anchored ranking for the rest, and score against the full survey. That is a one-day script and it has not been run.

**WHEN — the honest answer is "we cannot say, and here is why."** All four surveys are single cross-sections; three of four sat in the dry or post-harvest window (Gambia Jan-Apr 2018, Ghana Apr-Jun 2017 with 82 of 90 clusters in May, Malawi Dec 2015-Feb 2016 in the lean season of an El Niño year, Sierra Leone Nov-Dec 2013), so season cannot be separated from country. Every fieldwork-timed covariate was tested and added under 0.01: fieldwork-window food prices (+0.005 in-fill, ±0.003 transport), fieldwork-month climate including night LST and CHIRPS (+0.001 / −0.002), survey-year rainfall and temperature anomalies (0.341 with vs 0.346 without). The reason is structural: *a covariate defined by when the survey happened cannot rank districts that have never been surveyed.* **Recommendation to say aloud:** time the survey to the programme's calendar, and record the fieldwork window per cluster so the next panel can test seasonality. A stated ignorance is worth more to a policy room than a silence they will discover later.

**HOW — three levers, priced in survey-sample units, never in currency.**
1. *One measured national number plus the ranking* — median district error 9.3 pp at a 5% national anchor (≈ 42 respondents per outcome), against 11.7 for a same-cost regional survey and 17.3 for a same-cost district survey; a district survey matches it at 25-40% of a full sample. Honest limits: the anchored design's error floors at ~8.5 pp however large the anchor, because the ranking is the error; it loses to a district survey on burden captured at *every* size (0.22 vs 0.28-0.41); and at 15% or more of full size a regional survey beats it on level.
2. *More clusters per district* — the ceiling argument above.
3. *Do not expect a smaller survey.* Shrink the survey and the model's district error degrades faster than the survey's (index +6.7 pp from full sample to 15%, survey +3.3 pp). The sentence the partnership adopted: *"at any given survey size the model improves the ordering of districts; it does not reduce the sample a survey needs."* **The table behind this (`survey_size_symmetric_summary.csv`) is dated 3 September and pre-dates the 15-September outcome reconciliation — re-run it or make the claim without the numbers.**

## F3. Which proxy data to collect next

**The headline finding is a negative worth saying out loud: the vocabulary is saturated.** Eleven predictor blocks added this year each moved the score by less than 0.01 at the median (two engineered soil blocks moved it by 0.014-0.019, still inside noise). The open, permission-less tier transports *better* than the full set (0.316 vs 0.302 vs 0.281 with DHS); climate alone (0.335) and soil alone (0.327) each beat all 383 headline columns together (0.302); and DHS microdata add 0.001 inside a country (0.400 vs 0.399). The honest advice to a ministry: *you already have the predictor data. Do not spend the next grant on more layers.*

**The one real gap, and it is the one that matters for policy: programme reach.** Of 575 columns, the "food fortification and supplementation" domain holds 35, of which 25 are single national constants and only 10 vary between districts — seven of those are DHS and excluded from the headline set, leaving **three** (two salt-iodine indicators and one Malawi-only VAS indicator). So the model cannot distinguish "high risk because unreached" from "high risk despite good coverage", and that is precisely why the answer to "did the programme work" is No. Highest-value new block, in order:
1. **DHIS2 / HMIS district VAS-round and IFA-distribution counts** — the first non-public block worth requesting from ministries.
2. **HCES fortified-vehicle purchase shares** (wheat flour, oil, sugar, salt) — computable now from microdata already on disk for The Gambia, Malawi and Sierra Leone; about a day per country.
3. **MICS cluster GPS for Ghana 2017-18 and Malawi 2013-14** (request with UNICEF) — turns the MICS block from region-level to district-level in the two countries with the most districts; 72 of Ghana's 383 headline columns currently carry only 10 or 16 distinct values across 260 districts.
4. **HCES diet modules for Ghana (GLSS7) and Sierra Leone (SLIHS item codes)** — in-country value only.
5. **Not worth collecting** on present evidence: more remotely sensed variants of what is there (engineered climatology, terrain, soil-bioavailability, the satellite embedding as its own domain), FluNet, LSMS region means, time-matched market prices.

**And the honest framing of "saturation":** it is saturation *of this estimator at 14-87 districts per country*. The effective sample size is the number of districts, not the number of people. A block can carry real, replicated signal and still not move a five-fold in-fill by 0.01. More districts — which means more countries — is the way to make new blocks scoreable, which is the same argument as pooling.

## F4. The ask to the room

Three asks, in this order, each with a check the room can hold the project to:

1. **"If your country has a micronutrient survey — recent, planned, or sitting in a drawer — talk to us."** Every survey added measurably improves every other country's map: 0.20 with one training survey, 0.30 with three, and the curve has not flattened. *Honest qualifier to say:* the gain so far is in The Gambia and Ghana; Malawi and Sierra Leone as held-out countries have not improved. That is itself a reason to want a fifth country.
2. **"We want two implementation partners willing to make one real targeting decision twice — once with the regional averages, once with our map — and let us publish what changed."** This is the only route by which the worst-third precision numbers ever get validated prospectively, and it converts the project's weakest scientific position (no external validation) into a request rather than a hole.
3. **"To survey designers: two or three clusters in a district rather than one district in two more places."** Framed as a shared ceiling — *"that is the bar every district estimate is judged against, ours and your own survey's."*

**The question back to the room, and the workshop agenda item:** *"How accurate does a district ranking have to be before you would use it? We can tell you what we get. Only you can tell us what is good enough."*

## F5. The fifth survey and pre-registration

Seven quantitative predictions were written down and dated **3 September 2026** (`docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md`), for **Ethiopia and Pakistan**. The central prediction (P1) is a district Spearman above zero in at least 75% of outcomes with a mean above the country-block permutation null's 95th percentile; the recipe test (P2) is that climate+soil lands within 0.03 of, or above, the full index in at least 90% of cells. Tanzania is a fifth country whose data are already in hand (TDHS 2009-10 vitamin A, 164 Admin-2 units, merged dataset built, RBP unit trap fixed, not yet in any protocol-v2 table) — **so name Tanzania as an incorporation, not as a pre-registration target**, or the credibility device is inverted. The falsifier, stated plainly: if P1 fails in Ethiopia or Pakistan, the transport result on these four countries was a property of West and Southern Africa rather than of the proxies.

This is also the correct answer to the sharpest methods question the talk will face — *"you chose climate and soil on the same 22 cells you score them on."* The honest reply is on the record: **yes; the honest nested-selection figure is 0.31, the same as the full index at 0.30; climate+soil at 0.38 and the five-domain candidate at 0.41 are descriptions of these four countries and pre-registered candidates for the fifth, not validated results.** Quote **0.30** as the transport number on every main-body slide.

---

# Part 7. Appendix plan

The main body is 8 slides; the appendix is long and is the Q&A instrument, grouped so one question is answered by jumping to one group. DETAIL notes from the 90-slide deck carry over verbatim.

**A. "Is it better than our regional estimates in *our* country?"**
Per-cell forest (in-fill / region / country) with the jackknifed regional average as a cross on every row, including the cells we lose (Ghana child vitamin A, Ghana women's folate, Malawi child iron, both zincs). Per-country transport dot plot with the permutation null shaded at 0.08. The six strong cells labelled "chosen after the fact". Sierra Leone's missing in-fill rows explained (14 districts, 12-district floor). Fold-draw stability (strong cells move by sd 0.01-0.03 across the ten draws).

**B. "What did you actually fit, and did you try the complicated things?"**
Fourteen candidates on identical folds, 1,304 fits: index 0.294, rank-loss ensemble 0.285, random forest 0.265, squared-error ensemble 0.246, PC-HAL 0.215, elastic net 0.189 (prevalence, 17 cells). The geostatistical (DHS-style) comparator: it estimates levels better (9.2 vs 10.7 pp) and ranks worse (0.27 vs 0.38, 20 of 24 cells) — *"if you want an order use ours; if you want a number use theirs."* The covariate-free neighbour smoother (0.390 in-fill against the index's 0.399; it wins in Ghana; it cannot exist in an unsurveyed country). What "zero tuning" does and does not mean.

**C. "How good could any model be, and how wrong are your error bars?"**
Reliability ceiling by country on non-degenerate fits (Gambia 0.699, Ghana 0.538, Malawi 0.564, Sierra Leone 0.622, prevalence), with clusters per district beside it and the explicit statement that 26 of 44 Admin-2 variance fits are singular and the zero-ceiling cells are excluded rather than plotted as zero. **Never a share-of-attainable percentage.** Rank-interval coverage: 38% against a nominal 90%; the honestly calibrated 90% interval is ±55% of the list against ±68% for a random ordering; median miss 20% of the list against 29% random. The worst-fifth calibration ladder (16 / 27 / 34 / 36% across the four probability bands, base rate 20%).

**D. "Where did the data come from, and did your survey predict itself?"**
One CONSORT chain on the headline tier: **575 built → 383 that vary within a country and need no DHS → 344 present in all four countries → 94 (climate+soil) or 125 (five-domain) deployed in Côte d'Ivoire**, with the 23 written exclusion rules and their evidence files. **Drop the deck's 614 literal and the 525/464 pair, or label the latter "including DHS".** Alignment rules and worst-case year gaps (78 static layers, 136 at the exact survey year in every country; the worst gaps fall on layers the model gives little weight). The leakage policy: same survey *instance*, not analyte name — 70 GMNS-hosting MICS clusters dropped from The Gambia, 105 MNS clusters from every Malawi DHS aggregate. The honest qualifier: IHME and Malaria Atlas inputs cannot be audited from here, Malawi is the one country where a same-household path exists through IHME, and dropping the anaemia domain costs 0.02 of transport. The unit definition per country.

**E. "What is it actually looking at?"**
Domain drop-one ablation (climate +0.025, soil +0.025, anaemia surfaces +0.020, agriculture +0.012, infection +0.010; 8 of 21 domains cost nothing or help when dropped) with the caveat that weight share tracks how many axes a domain contributes (Spearman 0.79 with the number of PC axes). Top-20 layers for child iron with plain names and corrected colour groups. The Malawi B12 fish/anaemia/malaria example as a **surrogate, not a cause** — with the maps, not the forest. The urbanicity check: hold urbanicity constant and transport moves 0.378 → 0.364; urbanicity alone is −0.074.

**F. "Can we survey less, and where should we survey?"**
The four-design cost curve with the 5% anchor, the 25-40% crossover, and the burden loss drawn on the same panel. The symmetric survey-size test (re-run first). The clusters-per-district ceiling analysis. The explicit statement that we have **no** evidence on *when* to survey, and that the "pale = survey here" heuristic is a proposal the data currently contradict (firm-in districts have the *largest* rank error).

**G. "Give me my country's map."**
Côte d'Ivoire, all six outcomes, both candidates, agreement and rank uncertainty, the transport guards, the external WFP/MIMI check with its family size (10 pre-specified tests; in the wider 28-test cut nothing survives BH, min q = 0.138) and the explicit note that iron corroborates nowhere. The dashboard and the country-brief template — **and if the country brief does not exist by 28 September, label the card "in preparation" rather than showing a mock as a deliverable.**

**H. "What about the national number?"**
The VMNIS slide as specified in Part E, labelled "earlier sub-pipeline, not re-audited", replacing the current mislabelled "National-level track" slide.

**I. "What did you try that did not work?"**
The person-level SuperLearner (retained as a sensitivity, not the headline). The eleven added predictor blocks that each moved the score by under 0.01. The blend of model and survey estimate that made district numbers *worse* (12.2 vs 9.3 pp at full sample). Burden capture at chance under transport. The **withdrawal ledger** — four rows with dates and direction arrows: "12 of 12 regions at 0.50-0.56" (withdrawn, fold bug), "26% more burden than regional averages" (withdrawn), "two-thirds of attainable accuracy" (withdrawn, the ceiling was a diagnostic not a bound), and the covariate-free baseline at 0.516 that beat us (**also withdrawn — one in our own favour**). This is the single highest-credibility object in the appendix and the first thing to promote if the budget is ever relaxed beyond 12 slides.

**J. Presenter's own cards.**
The "statements that must not be made" list (Part H below), the jargon-to-plain table, the definitions card ("held-out", "rank correlation", the district unit per country, BRINDA adjustment, the fieldwork windows), and the Q&A ambush list with one-sentence honest answers.

---

# Part 8. Open questions for you

1. **Worked cell for Ghana — child iron, women's B12, or keep child vitamin A?** Recommended: child iron, which requires rebuilding `fig7_ghana_map.png` and producing `deploy_ghana_child_iron.csv`. **This is the binding build decision**; if it slips, the talk falls back to the cell where the regional average beats us. Confirm by when the rebuild can be done.
2. **Flagship Côte d'Ivoire outcome — women's vitamin A (externally corroborated, q = 0.043, survives latitude) or child iron (the nutrient the room cares about, externally uncorroborated at 0.178)?** You cannot show an iron map and a vitamin A corroboration on the same slide without saying so aloud.
3. **Does the odds ruler carry the neighbour-smoother tick (59.1%)?** It shrinks the visible margin from 4 points to under 1 and is the first methods question from the floor. Recommended yes, with the clause "geography gets you most of the way inside a country; the covariates are what let you cross a border."
4. **Eight slides or twelve?** The relaxed variant in B1 is a 19-20 minute talk, not a 15. If the session chair is strict, the 12-slide version cannot be delivered; decide before building.
5. **Is the country brief real by 28 September?** If not, remove it from any product shelf or label it "in preparation". Showing a mock deliverable to a funder room is a commitment.
6. **Is `survey_size_symmetric_summary.csv` re-run before the talk?** If not, the "it does not let you run a smaller survey" claim has to be made qualitatively, without the +3.3 / +6.7 numbers.
7. **Is the ceiling-by-clusters projection committed and re-run?** It currently exists only as `scratchpad/ceiling_by_clusters_projection.csv`; slide 8's centrepiece depends on it.
8. **Does slide 5 get a per-district leave-one-country-out prediction export** (for the right-hand panel of the "we deleted Ghana" twin scatter), or does that panel become rank-versus-rank? No protocol-v2 script currently writes those predictions.
9. **How hard should Sierra Leone's failure be said?** Recommended: one clause on slide 5 ("in one of our four it failed, and we can tell you in advance which countries look like that one"). It will be asked from the floor; volunteering it costs 8 seconds and buys the room.
10. **Cost language.** There is no currency figure anywhere and the recommended spine keeps it that way. But the single most decision-relevant fact for a director deciding whether to invite the team — *"the deployed climate-and-soil tier is 94 columns, about one to two days of a GIS analyst including QC, nothing to license"* — is currently on no slide. Does it go on slide 4's bottom strip, or into Q&A only?
11. **Sonja's slide 3: does the greyed Côte d'Ivoire column stay with her** (a tease delivered by the speaker who does not pay it off, across a handover), or move to Andrew's slide 4 build 1? If it stays, she needs the explicit handoff line.
12. **Does the talk name a user?** No version currently names who acts — Ghana Health Service, GroundWork, Côte d'Ivoire's Programme National de Nutrition, UNICEF country offices. One sentence on slide 6 ("for the PNN in Abidjan, this list says which ten of thirty-three districts to reach first while they wait for a survey") would make the product concrete; it also commits the project publicly. Your call.
13. **Ghana national child iron: 20.1% (population-weighted, the pipeline's own figure) or the GMS 2017 published figure?** They should be reconciled against the published report before a number goes on a slide in Accra. The 22.6% used in three of the five proposals is an n_eff-weighted mean of surveyed districts and should not be used.
14. **Does the withdrawal ledger stay in the appendix, or is it worth a 20-second beat in the main body?** It is the most differentiating object the project owns and the only one that would make the talk memorable to methodologists. At eight slides it does not fit; at twelve it does.
