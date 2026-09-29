---
title: "MNF15 joint talk: critique of the 27 September draft and a revised outline and script"
subtitle: "Micronutrient Forum 2026, Accra. Sonja Hess and Andrew Mertens. Session on national nutrition surveys, Thursday 1 October."
date: "2026-09-27"
---

# 0. The two constraints the current draft does not meet

**Time.** In the 18 September check-in Sonja said the slot is 15 minutes, "plan for like 13", with questions pooled at the end of the session (20 minutes for all speakers). The current script is 2,302 spoken words over 15 slides (Sonja 459 including the title, Andrew 1,843), before Sonja's two new slides. With those it is roughly 2,450 to 2,500 words: 16.5 minutes at 150 words a minute, nearer 18 at the 135 to 140 a minute an international room needs. The revised script below is 1,474 words: about 11 minutes at 135 wpm, 11.5 with clicks and the handoff, which leaves a pause after each number and a minute in hand.

**Slides.** The Forum's presenter guidelines say "limit the number of slides to a maximum of one slide per two minutes of speaking time". For 15 minutes that is 7 or 8. The current draft is 17 (15 plus Sonja's two). The revision below is 11 including the title; section 5 shows the three merges that bring it to 8 if you want to meet the guideline exactly.

# 1. Critique of the current draft

What works and should survive: the thesis line ("Modelling extends a survey's reach. It does not replace the survey."), comparing against the regional average a programme uses today rather than against zero, the honest hold-out design, the Binduri worked example, the decisions slide as a concept, the clean icon style of the concept slides, and a deep appendix for the pooled Q&A.

## 1A. Accuracy problems an expert in the room could catch

These are worth fixing even if the structure stays as it is.

| # | Slide | Problem | Fix |
|---|---|---|---|
| 1 | 5, Ghana map | The caption says "Adding environmental, dietary and health data predicts district deficiency better." For this exact example (Ghana, children's iron, average status) the covariate-free neighbour smoother scores 0.526 and the model 0.529: a tie (`benchmarks_v2_cells.csv`). Across all 18 in-country combinations the means are 0.390 and 0.399, and the smoother's median is higher (0.459 vs 0.407). A geostatistician will ask. | Say it on the slide: regional average 0.36, neighbour map 0.53, model 0.53. Then use the result you are not showing: with Ghana's survey removed entirely, the model trained on the other three countries scores **0.55** for the same outcome. That is the best argument in the project for public data, and it is set in Ghana, in front of the Ghana team. |
| 2 | 11 and 12, Côte d'Ivoire | The map is **children's iron** (script 06 writes `civ_rank_uncertainty.csv` only for child_iron x climate+soil). The "2007 survey agrees" check is **women's B12**. The room will assume the map they just saw was validated. The notes of appendix slide 30 also say "where iron can be checked it does not corroborate", which comes from the Tang dietary-inadequacy comparison and is now contradicted for biomarker data by XV-01 (iron positive in 10 of 11 African tests, and on average in every country). | Label the map with its nutrient, and back it with XV-01's iron results. Keep the 2007 B12 check as one spoken sentence. Update the slide 30 note. |
| 3 | 11, speaker note | "Dark means the model does not change its mind when we retrain it." The legend says the opposite: light is "exact", dark is "could move 9 places". | Invert the palette so dark = firm and relabel "How firmly each district is placed" (the 18 Sep ANM copy used that label), or fix the note. Do not call it "uncertainty": your own appendix shows these 90% rank intervals cover 38% under a held-out country. It is stability. |
| 4 | 14 and 15 | "Every survey added improves every other country's map", in the title, the footnote, the notes ("It improves the map of every other country in it") and the conclusions. Your backup note says the gain is in The Gambia and Ghana as held-out countries; Malawi and Sierra Leone have not improved. The 18 Sep ANM copy had that qualifier on the slide; the current render replaced it with the overclaim. The hollow "your survey here" point is a projection drawn like data. | "On average, each survey added improved the ranking for countries held out (0.20 with one survey, 0.30 with three)." Drop the projected point or label it "projected". |
| 5 | 13, decisions | (a) "District prevalence in a surveyed country: yes." Today's CP-01 90% bands are a median of +/-26 points, +/-30 to 37 for iron, and the dashboard the audience will open shows them. (b) "Today: an equal spread" for survey sampling: national surveys draw clusters with probability proportional to size within strata, and a survey designer in this session will say so. (c) SP-01 (today) shows model-guided district choice is about equal to random choice (0.362 vs 0.357 at half the districts; the gain is over population-weighted choice, 0.312). | (a) Move district prevalence to "use with survey data, planning figure, wide bands". (b) "Today: probability sampling within regions." (c) Frame the model's survey role as "which districts to confirm first", not a precision gain. |
| 6 | 9, by nutrient | The title says accuracy varies "not [by] prevalence", but prevalence is not on the figure. The folate note says the immunoassay "reads systematically lower": a within-country rank is unchanged by a constant offset, so a level shift cannot by itself collapse a ranking, and a biomarker person will say so. | Retitle ("Some deficiencies can be mapped from public data; others cannot yet"), or add prevalence to the row labels. For folate: "different assays and cut-offs, so the three surveys' folate is not the same measurement". |
| 7 | 8 and 9, vitamin A | Vitamin A is presented as best-predicted alongside B12. Externally it is the weakest nutrient: positive in 12 of 14 African tests but at a mean of 0.23 against 0.42 for iron and 0.70 for B12, and both clear wrong-way results (Zambia children, Ethiopia women, each -0.43) are vitamin A. In Ghana children's vitamin A the regional average beats the model (0.41 vs 0.28). This matters to a room that runs vitamin A supplementation. Note also that the dashboard's "Vitamin A fails in the African deposits" is stronger than its own dots show. | Say it plainly on the nutrient slide: rank well inside the four countries, weakest externally, do not lean on it for vitamin A supplementation yet. Soften the dashboard line to "weakest in the African deposits". |
| 8 | 5, 13, notes | "Geographic interpolation" is used for the jackknifed regional average (slide 13 footnote, notes), while slide 8 lists "a map of neighbouring districts" separately. Geospatial experts will read "interpolation" as the smoother. | "Regional average" everywhere; "neighbour map" for the smoother. |
| 9 | Sonja's slide 1 | "Targeted machine learning." In a Berkeley co-authored talk, the ML and statistics people will hear TMLE. The estimator is an untuned weighted index. | "Machine learning guided by conceptual frameworks." |
| 10 | 5 notes | "A programme working from the regional average sends its first team somewhere else." The regional average ranks Binduri 12th of 75, inside the worst sixth; any top-15 list includes it. | Keep the example, drop that sentence. |
| 11 | talk vs dashboard | The dashboard's start tiles read 0.40 in-country (regional 0.31) and 0.38 for a new country (climate and soil alone); the talk says 60% of pairs and 0.30. Some of the audience will open the app during the talk. | One scale, the dashboard's (see 1B), and quote the same numbers. |

## 1B. Structure and clarity

1. **The payoff comes too late.** The deck runs problem, data, framework, map, predictors, methods, three accuracy slides, importance, transport, and only then "what decisions can this support" at slide 13 of 17, eleven minutes in. A policy audience needs "what this does for me" early, with the evidence attached to each use.
2. **The accuracy section loses people.** The odds ruler pools 18 nutrient-country combinations and shows the model at 60 against the neighbour map at 59. In the check-in Sonja said she could not follow it: "our model is 1% better than looking at the neighbouring district? That's not particularly good." The slide is unchanged. The problem is pooling two different situations. Split by situation and the answer is clear: with a survey, geography does most of the work (say so); without one, there is nothing to smooth and the model still ranks districts (Ghana held out: 0.55).
3. **Five number scales.** 0.52 (slide 5), 60% of pairs (8), 0.63 (9), 0.95 (12), 0.20 to 0.30 (14). Pick one. Recommendation: the dashboard's "match score, 0 = guessing, 1 = perfect match", shown on the same small ruler each time, with one plain-language translation on the no-survey slide ("flags a district as worst-third and is right about half the time, against a third by guessing"). The pairs ruler goes to backup.
4. **The strongest trust evidence is missing.** XV-01/02 (22 September): six countries never used in training (Zambia, Ethiopia, Sudan, Nigeria, Pakistan, India), 36 of 41 tests pointing the right way, 0.33 to 0.40, no worse than inside the panel. The talk instead rests its external claim on nine 2007 zones for one nutrient.
5. **The thing a policymaker can use is hidden.** The dashboard appears only in the appendix and there is no QR code.
6. **Two talks, one story.** Your 3-minute short oral carries the predictor-importance story and the Malawi lakeshore warning. The main talk can drop its domains slide and point to the short oral, which saves 70 seconds and avoids duplication.
7. **Session context.** Three talks in the session cover why national micronutrient surveys matter (Sonja's worry about duplication). Open by building on them ("how do we get more out of each one?") and make one use explicitly about survey design; Fabian Rohner, whose GroundWork team partnered on surveys in this panel, speaks in the same session.

## 1C. Sonja's two slides

- **Introduction and objectives.** Text-only, and the concrete version of its first bullet ("information on prevalence is limited") is already slide 2's three Ghana maps. Suggest folding her objective into slide 2 as the question under the maps, and speaking the workshop sentence. "Targeted machine learning" per row 9 above.
- **Conceptual framework of iron deficiency.** The right content for her purpose (showing the room this is "not a wild goose chase"). But about 25 boxes and three footnotes are unreadable in a hall. The figure is a 949 x 693 px screenshot, soft when projected at full width, and it has a text cursor baked in after "Ecology". The footnote box runs past the bottom edge of the slide (in the LibreOffice render the source line is cut off; check it in PowerPoint).
- **The coverage slide does not match it.** Slide 4 re-draws the framework with different tiers and boxes: physiological vulnerability moves from intrinsic risk to "distal drivers"; health policies and culture move up a tier; food insecurity and health services move from intermediate to underlying; infectious disease, inflammation and gynaecological conditions move from direct causes to intermediate; obesity, family planning, health knowledge, built environment and other micronutrient deficiencies disappear, and genetic risk survives only as inherited red-cell disorders. The script calls it "the same framework". Sonja will notice, and so will anyone who knows the paper.
- **Fix, and it saves two slides:** one slide. Her figure as published, then on a click the same boxes get green, amber or grey badges. The green/amber/grey call for each of her boxes needs her sign-off; the missing boxes above need a call (for example, obesity is green via IHME overweight and DHS BMI; genetic risk is green via Malaria Atlas sickle-cell surfaces).

## 1D. Visual design

- **Titles are topics** ("How accurate is it?", "Where the predictors come from"). Make every title the sentence you want remembered.
- **Captions are too small.** Six figure slides end in 9 to 10 pt grey captions that no one past row five can read. Move the sentence into the title or the notes; nothing on screen below about 16 pt.
- **The worked example lives only in speech.** Circle and label Binduri and Atiwa East on the Ghana maps.
- **Palettes change meaning.** Teal is prevalence on slide 2 and rank on slide 5; Côte d'Ivoire's rank map is red-orange; its stability map is grey with dark meaning less stable. Use one "darker = worse" rank palette on every rank map, and one firmness channel where dark always means firm.
- **The decisions slide is a text table.** It is the slide people will photograph; give it three columns with icons and the QR code.

# 2. Proposed spine: problem, approach, three uses, limits, ask

11 slides including the title, 1,474 spoken words (Sonja 310, Andrew 1,164), about 11.5 minutes at 135 wpm including clicks and the handoff.

| # | Speaker | Title (the takeaway) | Visual | Time | Words |
|---|---|---|---|---|---|
| 1 | Sonja | Title (abstract title unchanged) | Existing title slide | 0:15 | 29 |
| 2 | Sonja | A national survey measures regions. Programmes act on districts. | figL, enlarged, with two big-number callouts and Sonja's objective as the question under the maps | 1:00 | 127 |
| 3 | Sonja | We started from what causes deficiency, not from what data exist | Sonja's framework, then coverage badges on click; four-survey strip at the bottom | 1:10 | 154 |
| 4 | Andrew | Learn where blood was drawn, apply everywhere, test only on what the model never saw | **New** flow diagram (replaces slides 6 and 7) | 1:10 | 152 |
| 5 | Andrew | Use 1. With a survey, the model ranks all 260 of Ghana's districts | fig7 with Binduri and Atiwa East circled; three-chip score strip | 1:05 | 139 |
| 6 | Andrew | Without a single Ghanaian blood sample, the model ranks Ghana's districts just as well | **New** two Ghana maps (survey, model trained elsewhere) and per-country strip | 1:10 | 152 |
| 7 | Andrew | Use 2. No survey at all? A first map, checked in six countries the model never saw | fig10 left panel (children's iron) and **new** "where it was checked" map | 1:10 | 156 |
| 8 | Andrew | Some deficiencies can be mapped from public data; others cannot yet | **New** traffic-light table by nutrient | 0:55 | 119 |
| 9 | Andrew | Use 3. Make the next survey count for more | **New** three-panel: clusters-per-district staircase, small national sample, planner screenshot | 1:10 | 159 |
| 10 | Andrew | What the model can and cannot do for you | Decisions slide redrawn as three columns, with QR code | 1:15 | 162 |
| 11 | Andrew | Modelling extends a survey's reach. It does not replace the survey. | Conclusions two-panel, pooling inset, QR code | 1:00 | 125 |

Moved to backup: four-survey table (now a strip on slide 3), predictor groups (slide 4 carries the icons), methods two-panel (slide 4), odds ruler, accuracy by nutrient bar chart (replaced by slide 8), domains (the short oral), the 2007 CIV B12 scatter (one sentence on slide 7), the stand-alone pooling curve (inset on slide 11).

# 3. Slide by slide

Scripts are drafts in the deck's voice: short sentences, British spelling, numbers only where they carry weight. Sonja's are suggestions for her to rewrite. Every number is traced in the "Sources" line.

## Slide 1. Title (Sonja, 0:15)

**Visual.** The existing concept title slide. Optional one-line subtitle in the policy register: *Ranking every district from free public data, so each national survey reaches further.*

> Good morning. I'm Sonja Hess from UC Davis. This is joint work with Andrew Mertens at Berkeley, and with the four national survey teams whose data made it possible.

## Slide 2. A national survey measures regions. Programmes act on districts. (Sonja, 1:00)

**Visual.** figL's three Ghana maps (one national number, regional averages, what the survey measured), larger. Two callouts at 28 pt or more: **185 of 260 districts: no survey cluster** and **62 of the 75 reached: one cluster, about a dozen children**. Across the bottom, Sonja's objective recast as the question: *Can free public data rank every district, even in a country with no survey?* This absorbs Sonja's introduction slide.

> You have just heard why national micronutrient surveys matter, so I won't repeat it. Our question is how to get more out of each one.
>
> Suppose you had to choose ten districts in Ghana for an iron programme tomorrow. Ghana's 2017 survey did exactly what it was designed to do: a good national number and good regional numbers. But Ghana has 260 districts. The survey reached 75. The other 185 have no measurement of their own; today they get their region's average. And 62 of the 75 it reached rest on a single cluster, about a dozen children.
>
> So we asked: can public data that already exist for every district, for free, rank them? And could the same model help a country with no survey at all?

*The opening sentence assumes the survey talks come first; drop it if the order changes.* Sources: `quantities.json` (districts, clusters, unsurveyed); `variance_components_ceiling.csv` (single-cluster share, 82.7% of 75 = 62).

## Slide 3. We started from what causes deficiency, not from what data exist (Sonja, 1:10)

**Visual.** Sonja's framework as published, high resolution, footnotes moved to the notes, source line visible. On click, a coloured badge on each of her boxes: green = a public district-level layer, amber = national or modelled only, grey = nothing public. Bottom strip: four country outlines with name, year and districts (The Gambia 2018, 30; Ghana 2017, 75 of 260; Sierra Leone 2013, 14; Malawi 2015-16, 87 Traditional Authorities). Replaces slides 3 and 4 of the current deck.

> We did not start from satellites. We started from a workshop of nutritionists, epidemiologists and data scientists, and from what causes deficiency. This is our conceptual framework for iron.
>
> [click] Then we went box by box and asked: is there a public, district-level measure of this? Green means yes: climate and ecology from satellites; poverty, schooling, water and sanitation, and infant feeding from household surveys; malaria and worms from disease atlases. Amber means only a national figure or a modelled surface, such as health policy or inflammation. Grey means nothing public reaches it: cultural norms, absorption.
>
> So the predictors are not a fishing expedition. They are the parts of this framework that public data happen to cover, and the model can only see those parts.
>
> To learn from, we had four national biomarker surveys: The Gambia, Ghana, Sierra Leone and Malawi, 206 districts where blood was drawn. They exist because these teams collected them. Andrew.

*Badge colours for each of Sonja's boxes need her sign-off (see 1C).*

## Slide 4. Learn where blood was drawn, apply everywhere, test only on what the model never saw (Andrew, 1:10)

**Visual (new).** A left-to-right flow. (1) Four country outlines: "206 districts with blood samples". (2) A stack of eight layer icons (the existing predictor-group icons): "575 public layers, 28 sources, every one of 554 districts, none bought". (3) "A simple weighted index; nothing tuned". (4) A small ranked map: "every district, worst to best". Band underneath, three icons: *hide one district in five*, *hide a whole region*, *hide a whole country*, then "compared with the regional average programmes use today". One small line: "14 machine-learning alternatives tried; none did meaningfully better at 14 to 87 districts per country". Replaces slides 6 and 7.

> Thank you, Sonja. The method in one picture. On the left, the 206 districts where blood was drawn. In the middle, 575 public layers we assembled for every district in the four countries: climate and soil, crops and livestock, malaria, food prices, household surveys. None of it had to be bought.
>
> The model learns which combination of layers tracks deficiency where we have blood, and applies it everywhere. We tried fourteen machine-learning methods, including random forests, boosting and ensembles. With 14 to 87 districts per country, none did meaningfully better than a simple weighted index with nothing tuned, so that is what we use.
>
> Every number I show you comes from hiding data and predicting it: one district in five, a whole region, a whole country. And we compare against what a programme uses today, the regional average. The score runs from zero, guessing, to one, a perfect match with the survey.

Sources: `predictors_admin2_shared_metadata.csv` (575, 28); SL-06 (14 alternatives, 1,304 fits, none better by more than 0.03).

## Slide 5. Use 1. With a survey, the model ranks all 260 of Ghana's districts (Andrew, 1:05)

**Visual.** fig7 (survey, model with each district hidden, model on every district), with Binduri and Atiwa East circled on the first two panels and named. Under the maps, three chips on the match-score ruler: *Regional average 0.36*, *Neighbour map 0.53*, *Model 0.53*. No caption paragraph.

> Ghana, children's iron. Left: what the survey measured in the 75 districts it reached. Middle: the model's prediction for each of those districts with that district's blood results hidden. Right: the model applied to all 260.
>
> Take Binduri, in Upper East. The survey measured it 6th worst of 75. Its regional average puts it 12th. The model, which never saw Binduri's results, puts it first. And we are wrong in places: Atiwa East, which the survey put 24th, we put 56th.
>
> Overall the model scores 0.53 here, against 0.36 for the regional average. But I want to be straight with you: a map built only from neighbouring districts, with no public data at all, also scores 0.53. Inside a surveyed country, geography does most of the work. So why bother? Because most countries have no survey to borrow from.

Sources: `fig7_ghana_district_ranks_child_iron_level.csv` (Binduri 6/12/1, two clusters; Atiwa East 24/28/56); `benchmarks_v2_cells.csv`, Ghana child_iron level in-fill (index 0.529, spatial 0.526, region_mean_jk 0.356). Backup, if asked about prevalence rather than average status: in-fill prevalence is a tie with the regional average (0.50 vs 0.50).

## Slide 6. Without a single Ghanaian blood sample, the model ranks Ghana's districts just as well (Andrew, 1:10)

**Visual (new).** Two Ghana maps side by side, same palette: *What the 2017 survey measured* and *Model trained only on The Gambia, Sierra Leone and Malawi*. A large **0.55**, with "0.53 for the neighbour map, which needed Ghana's survey" under it. Along the bottom, the four held-out countries as chips: The Gambia 0.60, Ghana 0.35, Malawi 0.20, Sierra Leone 0.16, guessing below 0.08.

> So here is the harder test, run where we can check the answer. We removed Ghana's survey entirely. The model learned only from The Gambia, Sierra Leone and Malawi, and never saw a Ghanaian blood sample. It still scores 0.55 for children's iron: as well as the neighbour map that needed Ghana's own survey.
>
> Ghana's children's iron is one of our better cases, and I don't want to oversell it. Holding out each country in turn, averaged over nutrients, we get 0.60 for The Gambia, 0.35 for Ghana, 0.20 for Malawi and 0.16 for Sierra Leone, where guessing stays below 0.08. In practical terms: when the model puts a district in the worst third of a country it has never seen, it is right about half the time, against a third by guessing. That is a shortlist, not a measurement. But it is a shortlist for a country that today has nothing.

Sources: `benchmarks_v2_cells.csv`, Ghana child_iron level, estimand country, domain_index 0.554; `cell_master.csv` column `tr` (per-country means over the 22 transport cells); `transport_null_calibration.csv` (0.079); `quantities.json` hit3 52%, hit3_base 33%. **Needs a new figure**: the per-district predictions of that fit are not saved with district names (`rank_interval_districts.csv` has no Admin2), so re-run the Ghana-held-out child_iron fit and keep them.

## Slide 7. Use 2. No survey at all? A first map, checked in six countries the model never saw (Andrew, 1:10)

**Visual.** Left: fig10's ranking panel only, titled *Côte d'Ivoire, children's iron, no Ivorian data used*, top three districts named. Right (new): a map of Africa with a South Asia inset. Grey = trained on (4 countries), teal = checked against WHO-deposited survey results (Zambia, Ethiopia, Sudan, Nigeria, Pakistan, India), orange = Côte d'Ivoire. One line: *36 of 41 tests pointed the right way; guessing stays below 0.20.* Small print: regions and provinces, WHO VMNIS deposits.

> Most countries in this region have no recent biomarker survey. Côte d'Ivoire is one. Trained on our four countries, the model ranks its 33 districts for children's iron. The northern savanna comes out worst: Tchologo, Bounkani, Bagoué.
>
> Can you trust a map like this? We looked for surveys the model had never seen. The WHO micronutrients database holds regional results from national surveys in Zambia, Ethiopia, Sudan and Nigeria, and in Pakistan and India. None were used to build the model. In 36 of 41 tests its ranking of their regions pointed the right way, scoring 0.33 to 0.40, no worse than inside our own four countries. Iron pointed the right way in every one of the six. And Côte d'Ivoire's own 2007 survey ranks its nine zones for vitamin B12 almost exactly as the model does.
>
> These are regions, not districts, and the checks are published summaries. But they are data we did not touch.

Sources: `civ_rank_uncertainty.csv` (child_iron, climate+soil); `xv_transport_pooled.csv` (Africa level 0.402, 12/12; Africa prevalence 0.332, 17/21; off-continent 0.394, 7/8; country-block null 0.18 to 0.23); `xv_transport.csv` (iron country means all positive); `civ_b12_vs_2007_rho.txt` (0.95, nine zones). Backup: Tang et al. (WFP/MIMI) vitamin A map for Côte d'Ivoire agrees at 0.50 over 33 districts and survives latitude adjustment; against their modelled dietary inadequacy, iron does not corroborate, which is expected because iron deficiency is not only an intake problem.

## Slide 8. Some deficiencies can be mapped from public data; others cannot yet (Andrew, 0:55)

**Visual (new).** A traffic-light table. Rows: B12 (women); iron (women, children); vitamin A (women, children); folate (women); zinc. Columns: *inside a surveyed country*, *country held out*, *six external countries*, *use it?* Cells coloured green, amber or red with the number inside.

| | Inside a surveyed country | Country held out | Six external countries | Use it? |
|---|---|---|---|---|
| B12, women | 0.63 | 0.52 | 0.70, 4 of 4 | Yes |
| Iron, women / children | 0.44 / 0.43 | 0.39 / 0.28 | 0.42, 10 of 11 | Yes |
| Vitamin A, women / children | 0.58 / 0.49 | 0.54 / 0.49 | 0.23, 12 of 14, both clear misses | Not yet for supplementation decisions |
| Folate, women | 0.40 | 0.07 | 0.30, 3 of 4 | Within a country only |
| Zinc | not mappable | n/a | n/a | No |

> Not every deficiency can be mapped this way. Vitamin B12 maps best, inside a country and across borders, probably because it follows animal-source foods. Iron works and held up externally. Vitamin A ranked well inside our four countries but was the weakest in the external check, so for vitamin A supplementation we would not lean on it yet. Folate ranks districts within a country but not across borders, most likely because the surveys measured it with different assays. And zinc, measured only in Malawi, cannot be mapped by district: the survey's own district zinc numbers carry no signal. None of this follows how common a deficiency is. Folate and zinc are the commonest in our data, and the hardest.

Sources: figB labels and `cell_master.csv` (in-fill and transport, level); `xv_transport.csv`, iSDA arm, domain_index, by nutrient (B12 0.70 4/4; iron 0.42 10/11; folate 0.30 3/4; vitamin A 0.23 12/14); Ghana child vitamin A in-fill 0.28 vs regional average 0.41; `figA` measurability screen (zinc).

## Slide 9. Use 3. Make the next survey count for more (Andrew, 1:10)

**Visual (new).** Three panels. (a) The ceiling staircase, simplified to three bars: one cluster per district 0.54, two 0.66, three 0.73, labelled "the best any model could score against the survey". (b) An icon of a small national sample: "about 40 blood draws per nutrient, district figures typically within about 10 points". (c) A thumbnail of the dashboard's Plan-a-survey tab.

> The third use matters most right now, with survey budgets under pressure: making each survey go further.
>
> First, how districts are sampled. Across our four surveys, three districts in four rest on a single cluster. Against numbers that noisy, even a perfect model would score only about 0.55. Two clusters per district lifts that ceiling to about 0.66, three to 0.73, for our maps and for the survey's own district estimates. If district estimates matter, put two or three clusters in a district.
>
> Second, what to do before a full survey. A small national sample, about forty blood draws per nutrient, is enough to anchor the ranking; district planning figures then typically land within about ten percentage points.
>
> Third, where to look first. The ranking tells a survey or a programme which unsurveyed districts to confirm first, and the dashboard's planning tab lets you try a design.
>
> None of this replaces a survey. It helps each one reach further.

Sources: `variance_components_ceiling.csv` and the staircase projection (0.539 / 0.663 / 0.731, projected from the surveys' own variance structure; 74% of 206 districts single-cluster); `anchor_and_rank_summary.csv` (median district error 9.3 pp at a 5% anchor, about 42 respondents; regional survey of the same size 11.7). Backup, have ready: Sierra Leone has the most clusters per district (4.3) and the weakest results, because 14 districts of half a million people is the wrong unit, not the wrong design. SP-01: visiting half the districts keeps about 90% of the ranking accuracy, but model-guided choice is about equal to random (0.362 vs 0.357); both beat population-weighted choice (0.312), and design weights keep the national estimate unbiased. Do not claim the model makes a smaller survey possible: the symmetric sample-size test says the opposite.

## Slide 10. What the model can and cannot do for you (Andrew, 1:15)

**Visual.** The decisions slide redrawn as three columns with icons. **Use it to:** put districts in order to reach first; sequence a roll-out; choose which unsurveyed districts to confirm first. **Use it with survey data for:** a district prevalence for planning (anchored to a survey; 90% bands about +/-25 points); reaching the most deficient people (multiply by population). **Don't use it to:** replace a survey; judge whether a programme worked; decide what to change. QR code to the dashboard, bottom right, with its URL in 18 pt.

> So what can you use it for? Use it to put districts in order: who to reach first, how to sequence a roll-out, which unsurveyed districts to check first. That is where it beats what programmes use today.
>
> Use it together with survey data for anything that needs a number. The ranking has no level of its own. Anchored to a survey, it gives planning prevalences with honest, wide bands. And if the goal is to reach the most deficient people rather than the highest rates, multiply by population: in Ghana, simply going to the most populous fifth of districts reaches almost half of the deficient children.
>
> Don't use it to replace a survey, to judge whether a programme worked, or to decide what to change. The model says where deficiency is, not why.
>
> All of this is in a public dashboard: every district in the four countries and Côte d'Ivoire, with one-page country briefs. The QR code is on the screen.

Sources: `conformal_prev_cells.csv` (median half-width 26 pp, LOO coverage 0.90 to 0.93); `nce_targeting_summary.csv` (worst fifth by rate reaches 21.9% of deficient people, regional averages 20.3%, random about 20%, perfect 47.5%); `viz/oof_child_iron.csv` cap_pop (most populous fifth 46%). Backup: the Malawi lakeshore surrogate (anaemia marks where people eat fish) if asked why weights are not causes.

## Slide 11. Modelling extends a survey's reach. It does not replace the survey. (Andrew, 1:00)

**Visual.** The existing conclusions two-panel, rewritten: three findings on the left, three asks on the right, the thesis as the title rather than a banner. Small inset: the pooling curve (0.20, 0.26, 0.30, "so far, on average"). QR code repeated.

> Three things to take home. Free public data carry real information about micronutrient deficiency, much of it in climate and soil, measured everywhere, every year. That information ranks districts better than regional averages, and it works in countries the model has never seen. And the limit today is surveys, not data: on average, each survey we added improved the ranking for the countries held out.
>
> So three asks. If your country has a micronutrient survey, recent, planned or sitting in a drawer, share it with the pool. If you are planning one, let's design it together. And tell us how accurate a district ranking needs to be before you would use it.
>
> Modelling extends a survey's reach. It does not replace the survey. Thank you.

Sources: `training_curve_climate_soil.csv` (0.198 / 0.259 / 0.302; gains in The Gambia and Ghana as held-out countries, not yet in Malawi or Sierra Leone); eleven new predictor blocks this year, none moved the median by more than 0.01. Optional fourth line if Sonja wants it: "Next, with the MIMI and GBD teams, we will compare our predictors framework by framework."

# 4. Visual build list, in priority order

**P1, before anything else (accuracy, under an hour).** Fixes 1 to 11 in section 1A, whatever structure you keep: slide 5 caption and chips; Côte d'Ivoire nutrient label, firmness palette and note; the pooling overclaim in four places; decisions-slide wording; folate and vitamin A wording; "regional average" terminology; the dashboard's vitamin A sentence; QR code.

**P2, restructure (reuses existing figures).**

| Slide | Status | Work |
|---|---|---|
| 2 | Modify figL slide | Two callouts, the objective line, drop the caption |
| 3 | New build on Sonja's figure | High-resolution framework image from the paper; badge overlay; four-survey strip. Needs Sonja's box-by-box calls |
| 4 | New concept slide | A new `flow` layout in the concept YAML, or native shapes; reuse the Font Awesome icons already in the predictor grid |
| 5 | Modify fig7 | Callout circles and labels (02_figure_ghana_map.R); three chips |
| 7 | Modify fig10, new map | Left panel only; new "where it was checked" map from `xv_transport.csv` and the GADM outlines already on disk |
| 10, 11 | Modify concept YAML | Checklist to three columns; conclusions rewritten; QR image |

**P3, new figures.**

| Slide | Work |
|---|---|
| 6 | Re-run the Ghana-held-out children's iron fit keeping district names; two-map figure; four-country chips |
| 8 | Traffic-light table from `cell_master.csv` and `xv_transport.csv` (native pptx table or ggplot tiles) |
| 9 | Three-bar staircase from the projection CSV; anchor icon; dashboard thumbnail (figS_dashboard_planner.png exists) |

# 5. If you want 8 slides (the Forum guideline)

Merge 5 and 6 into one four-map slide (survey, regional average, neighbour map, model trained elsewhere, each with its score). Fold slide 8's table into slide 7 as a strip. Put the three asks on slide 10 and end on it. Words stay about the same; the talk gets denser, so I would keep 11 unless the session chair enforces it.

# 6. Backup, reordered for the pooled Q&A

Order the appendix by the question most likely to come, with a one-slide index at its top.

1. "Isn't it just a north-south gradient?" Spatial smoother result; the Côte d'Ivoire vitamin A ranking correlates 0.73 with latitude, yet its agreement with the Tang map survives latitude adjustment (0.495 to 0.479). *Not yet computed: a latitude-only transported ranking under hold-out. Worth an hour if time allows; someone will ask.*
2. "How do I get a prevalence, not a rank?" Anchor design; CP-01 bands.
3. "Why not machine learning?" 14 alternatives on identical folds (fig1).
4. "What is the ceiling?" figJ staircase.
5. "Will it work in my country?" figC by country; the external-check map.
6. "What drives it?" figD domains; Malawi lakeshore (figM). Point to the short oral.
7. "What does the uncertainty on the map mean?" Stability versus calibrated; the 38% check; CP-01.
8. "Which deficiencies can be mapped at all?" figA screen.
9. "Does targeting reach more deficient people?" fig6.
10. "How does this compare with DHS geostatistics?" fig9.
11. "Numbers we withdrew."
12. "What does the data cost?" figI (with MICS moved into the registration tier, per the 18 Sep check-in).

# 7. Decisions needed

1. Confirm the spoken budget with Sonja: 13 minutes for both of you?
2. 11 slides (recommended) or the Forum's 8?
3. Will Sonja fold her introduction into slide 2 and put the coverage badges on her own framework figure?
4. One scale: the dashboard's 0-to-1 match score (recommended) or pairs out of ten?
5. XV-01/02 and the Ghana-held-out result in the main body (recommended)?
6. Does the session order put the survey talks before yours? It sets Sonja's first sentence.
