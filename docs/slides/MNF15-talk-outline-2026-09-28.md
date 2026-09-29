---
title: "MNF15 joint talk: outline and script for deck v6"
subtitle: "Micronutrient Forum 2026, Accra. Sonja Hess and Andrew Mertens. Session on national nutrition surveys, Thursday 1 October."
date: "2026-09-28"
---

# 0. What this version is

**Deck.** `docs/slides/MNF15-talk-2026-09-v6.pptx`: 23 main slides including the title, then the appendix (37 slides). It follows the order Andrew set on 28 September: motivation, data assembly, methods, national prevalence, district prevalence, ranking, out of country, variable importance, policy relevance and next steps, conclusions. It replaces the 27 September outline (an 11-slide spine built around three uses); the deck has since followed Andrew's edits.

**Time.** 2,075 spoken words (Sonja 421, Andrew 1,654): about 15 minutes at 140 words a minute; the timings in the notes add up to 15:25. The slot is 15 minutes, "plan for like 13" (Sonja, 18 September check-in), with questions pooled at the end of the session.

**Slides.** The Forum's guideline is at most one slide per two minutes of speaking, 7 or 8 for this slot. v6 has 23: 22 once Andrew keeps one of the two targeting slides (13 and 14), including the two figures added on 29 September (the population cartogram and the Côte d'Ivoire applicability check). The pace (about 40 seconds a slide) is the main risk. Section 5 gives the cuts to about 10 minutes and 16 slides.

**Honesty rule for this version.** Overall performance is reported for every combination, including the weak ones; single-country or single-nutrient figures show cases where the model works (Ghana children's iron held out, women's B12), and each one says which setting its number comes from.

# 1. What changed since the 27 September outline

| Change | Why |
|---|---|
| Order: motivation, data, methods, national, district prevalence, ranking, out of country, importance, policy, conclusions | Andrew's structure, 28 September |
| Andrew's edits and 17 comments on v4 applied (v5) | Titles, legends, "n of N" districts, the national slide with every outcome, the prevalence-error chart with an X for the national estimate, no "full talk" tags, specific objectives (Sonja) |
| Numbers rebuilt with the 15 September survey fixes (27 to 28 September) | Vitamin A on RBP < 0.70; non-pregnant women and children 6 to 59 months only; iron and zinc on the surveys' own flags. The measurable set is still 16 combinations: Ghana women's vitamin A dropped out (1.6%), Malawi children's vitamin A came in (8.7%) |
| The survey's regional average is scored fairly: 65 of 100 pairs, not 61 | The earlier scoring left each district out of its own region's figure, which reverses two surveyed districts of the same region by construction. Scored as a planner uses it (every district in a region gets the same figure), the regional figures reach 65 against the model's 66, model ahead in 10 of 14. The model's value is where there is no survey |
| Targeting and B12 claims say where they hold | Inside a surveyed country the flagged fifth holds 25% of deficient people (B12: 39 to 46%). With the country's survey removed it holds 17% (B12: 11 to 12%), below a random fifth, because the model ranks by rate and the high-rate districts are small. B12 transports on average level (0.60, 0.44), not on the share deficient (0.24, 0.22) |
| B12 prevalence error: about 9 points, about as close as the survey's regional figures | Post-fix; the earlier 6 was a pre-fix median |
| Three slides to the appendix | "Both a big-data opportunity...", "How do we choose among algorithms?" (it said the ensemble wins; slide 9 says it does not), "How much of what a survey can see..." (its ceiling table is from 16 September) |
| Two merges | Variables: the figure plus the three strongest importance findings (the full findings slide is in the appendix). Closing: "Conclusions" and "What have we learned" became one slide that mirrors the outline |
| Old closing slide from v4-ANM (59) dropped | Pre-fix numbers and the "every survey improves every map" overclaim |
| Outline as five questions; the closing slide answers them | Andrew, 28 September: the earlier outline read as slogans |
| Two targeting slides side by side (13 text, 14 figure) | Andrew will keep one |
| Population cartogram (after the targeting slides) and Côte d'Ivoire applicability check (after the Côte d'Ivoire slide) | Andrew, 29 September: both in the main talk |
| Slide 'What variables drive model predictions' reworded | Leads with the drop-a-domain evidence (dropping climate or soil makes the held-out ranking worse in all four countries); 21 of 24 single-country models, not every; 'satellite imagery summary', not land cover; no 'about half the weight', which mostly reflects how many layers each family has |
| Measurability screen redone with the post-fix survey ceilings (28 Sep, 05:34) | Malawi children's vitamin A left the measurable set (ceiling 0.15, under 0.30) and Malawi women's zinc entered (0.37; in-country only). Every figure and number that averages over measurable combinations was redrawn: in-country pairs 65.6, fair regional figures 65.3 (model ahead in 9 of 14, was 10), held out 62 over 15; targeting 31% per person, most populous 49% |
| Slide 9 on the post-fix SuperLearner runs | In-country and region hold-outs from the 28 September runs (rank-normalised predictors, the production default), country hold-out from a post-fix rerun |
| Speaker notes: why average status, zinc, survey precision | Dichotomising at the clinical cut-off lowers the survey's own district reliability by 0.20 (continuous level more reliable in 25 of 27 combinations); children's zinc has no district-level variance (district share 0.000, cluster 0.172); 86% of direct district prevalence estimates have a CV above 33% |
| Targeting by budget rule (TC-01) | With a budget per person the model's highest-rate districts reach 31% of the deficient with a fifth of the people (random 20, regional figures 27); with a budget in districts the most populous fifth reaches 48% and the model adds nothing. Slide 13 and the appendix targeting figure |
| Only B12 has a map of its own (NX-01) | Swapping weights between nutrients under a country hold-out: B12 specific in 2 of 3 countries, iron 2 of 8, vitamin A 0 of 8. Slides 16, 17, 20 |
| Other checks folded into notes and appendix | Borrowing other countries' weights in-country fails (BO-01); on pairs the survey clearly separates the model gets 83, regional figures 80, neighbour map 83 (SEP-01); district levels in a new country are no better than the national figure (LV-02, redrawn appendix figure); about half of the prevalence error is survey noise (NZ-01) |
| Appendix refreshed on the fixed data | Survey ceiling (0.47, about two-thirds reached), vitamin A severity bands (model 50% exact, regional figures 53, neighbour map 54; was 64 and 66), the small-national-sample figure, the pooled-pairs ruler, the in-country against held-out figure |

# 2. The spine

| # | Speaker | Section | Title | Time | Words |
|----|---------|-------------|----------------------------|-------|-------|
| 1 | Sonja | Motivation | Title | 0:15 | 29 |
| 2 | Sonja |  | Introduction and study objectives | 0:40 | 96 |
| 3 | Sonja |  | Talk outline | 0:20 | 48 |
| 4 | Sonja | Data assembly | Conceptual framework of iron deficiency | 0:20 | 41 |
| 5 | Sonja |  | Which of these causes can public data measure? | 0:45 | 93 |
| 6 | Sonja |  | 575 proxy variables from 29 conceptual domains | 0:25 | 55 |
| 7 | Sonja |  | Which surveys did we use to train the models? | 0:25 | 59 |
| 8 | Andrew | Methods | How do we get honest estimates of model performance? | 0:40 | 82 |
| 9 | Andrew |  | Why not a more complex machine-learning model? | 0:45 | 90 |
| 10 | Andrew | National prevalence | But can proxy models recover national-level prevalence estimates? | 0:45 | 114 |
| 11 | Andrew | District prevalence | How far off is the predicted prevalence, outcome by outcome? | 0:40 | 69 |
| 12 | Andrew | Ranking | How well does it rank districts, nutrient by nutrient? | 0:50 | 124 |
| 13 | Andrew |  | What does that accuracy mean on the ground? | 0:50 | 117 |
| 14 | Andrew |  | Which districts should a programme go to first? | 0:45 | 90 |
| 15 | Andrew |  | Where do the young children live? | 0:35 | 84 |
| 16 | Andrew | Out of country | Can it rank Ghana's districts without Ghana's survey? | 0:50 | 99 |
| 17 | Andrew |  | What about countries with no survey? | 0:55 | 118 |
| 18 | Andrew |  | Is Côte d'Ivoire like the places we learned from? | 0:35 | 86 |
| 19 | Andrew | Variable importance | What variables drive model predictions | 1:00 | 146 |
| 20 | Andrew | Policy relevance and next steps | Which deficiencies can it map? | 0:45 | 93 |
| 21 | Andrew |  | What does this mean for policy? The vitamin B12 model | 0:50 | 128 |
| 22 | Andrew |  | Next steps | 0:30 | 69 |
| 23 | Andrew | Conclusions | Conclusions | 1:00 | 145 |
| | | | **Total** | **15:25** | **2,075** |

# 3. Slide by slide

The script is the speaker notes of v6, as they will appear in PowerPoint (18 pt). Sonja's lines are suggestions for her to rewrite. "Backup, if asked" carries the sources and the numbers for questions.

## Motivation

### Slide 1. Title (Sonja, 0:15, 29 words)

**Visual.** Title, authors, the three logos and the Gates acknowledgement.

> Good morning. I'm Sonja Hess from UC Davis. This is joint work with Andrew Mertens at Berkeley, and with the four national survey teams whose data made it possible.

### Slide 2. Introduction and study objectives (Sonja, 0:40, 96 words)

**Visual.** Sonja's bullets with the two specific objectives on the left; on the right the WHO VMNIS map of the latest national biomarker survey in each country, with counts by nutrient. Source as a footnote.

> Information on how common vitamin and mineral deficiencies are is limited. The map shows how limited for sub-Saharan Africa: 17 of 47 countries have no national biomarker survey on record, and only 10 have one from 2015 or later. A survey gives national and regional numbers, while programmes are planned district by district. So we brought together nutritionists, epidemiologists and data scientists to ask which approaches and data could help. Our objectives: to see whether machine learning on aggregated proxy data can estimate deficiency prevalence, nationally and by district, and find the districts at greatest risk.

*Backup, if asked:* map and bars from the WHO VMNIS Micronutrients Database (export of 25 February 2025), nationally representative surveys measuring ferritin, retinol or RBP, zinc, folate or vitamin B12 (scripts/policy_deck/33_mnf15_v4_vmnis_coverage.R); surveys not yet deposited are not shown, and all four of our surveys are in it. For zinc, only five countries have a survey from 2015 or later. Worldwide, in low- and middle-income countries over 1988 to 2018: vitamin A data in 77 countries, iron 53, folate 24, zinc 21, B12 7 (Brown, Moore, Hess et al., Am J Clin Nutr 2021, CC BY 4.0). The global estimate of 56 per cent of young children and 69 per cent of women with at least one deficiency rests on 24 surveys in 22 countries (Stevens et al., Lancet Glob Health 2022). Fortification is set nationally (Nyumuah et al., Food Nutr Bull 2012); coverage of fortified foods is lower among poor and rural households (Aaron et al., J Nutr 2017). Why model at all: by the conventional survey bands for the coefficient of variation, 86 per cent of the 1,350 direct district prevalence estimates in our four surveys have a CV above 33.3 per cent and 5 per cent are under 16.6 per cent (median CV: Ghana 100 per cent, Malawi 96) (explore/out/19_cv_bands.csv). Two guards: these are the survey’s own direct district figures, not the model’s, so this says modelling is needed, not that the model’s outputs are unreliable; and CV (standard error over prevalence) penalises rare outcomes by construction, so the estimates that pass are simply the high-prevalence ones.

### Slide 3. Talk outline (Sonja, 0:20, 48 words)

**Visual.** The five questions the talk answers, in talk order. The closing slide answers them in the same order.

> We’ll try to answer five questions. Can public data tell us how common a deficiency is, nationally or by district? Can they rank a country’s districts? Does that work where there has never been a survey? Which data matter? And how could programmes and survey planners use it?

## Data assembly

### Slide 4. Conceptual framework of iron deficiency (Sonja, 0:20, 41 words)

**Visual.** Sonja's published framework (Hess et al. 2023) with its footnotes.

> We did not start from satellites. We started from what causes deficiency. This is our conceptual framework for iron, adapted from our 2023 paper: intrinsic risk, the fundamental drivers, underlying and intermediate risk factors, and the direct causes of iron deficiency.

*Backup, if asked:* Hess et al., Ann N Y Acad Sci 2023, adapted.

### Slide 5. Which of these causes can public data measure? (Sonja, 0:45, 93 words)

**Visual.** The same framework with every box coloured by public availability: dark green intrinsically ecological, light green individual data linkable as a district aggregate, amber partly, grey not public, teal the survey outcome.

> So we went back to the framework, box by box, and asked: is there a public, district-level measure? Dark green: yes, and the data describe the place itself: climate, ecology, built-up land. Light green: yes, as a district average standing in for individuals: poverty, schooling, water and sanitation, infection. Amber: only a national figure, a modelled surface, or part of the box. Grey: nothing public. Everything here is a district average, never the person whose blood was drawn. So these layers mark where deficiency is likely; they cannot explain any one person’s status.

*Backup, if asked:* Inherited red blood cell disorders: the only layers are the Malaria Atlas 2012 modelled allele-frequency surfaces (sickle haemoglobin, haemoglobin C, G6PD deficiency), so amber, not green; population gene frequencies, not anyone’s genotype, filed with infection because they largely follow historical malaria, and about 1 per cent of the four-country model’s weight. Genetic risk beyond these (for example variants affecting iron absorption): nothing public, so grey. Fundamental drivers: climate, ecology and geography are measured everywhere; politics is conflict events only (ACLED), economy district wealth averages, and inequity one layer, the spread of relative wealth within the district, so that half is amber. Chronic disease: nothing beyond HIV. Health and nutrition knowledge: the only layers are schooling and literacy, already counted under education, so grey. Family planning, obstetric history, pregnancy and women’s BMI come from DHS: public, so coloured, though the headline model leaves DHS out. Food insecurity: IPC phase, food consumption score and coping index each cover two or three of the countries, none of them Ghana, which has household food spending and market prices. Blood loss: hookworm and schistosomiasis prevalence from WHO ESPEN, nothing on menstrual or obstetric loss. Inflammation: infection burden only, no CRP or AGP. Calls still marked draft for Sonja: built environment, family planning, obesity, health and nutrition knowledge, other micronutrient deficiencies (scripts/policy_deck/30_mnf15_v4_framework_drawn.py).

### Slide 6. 575 proxy variables from 29 conceptual domains (Sonja, 0:25, 55 words)

**Visual.** Andrew's slide: predictors by source (bars) and the eight domain groups with their counts.

> From that exercise we assembled 575 public variables from 28 sources, in 29 conceptual domains, grouped here into eight. Every one is a district average, and all of it is free. Because there are far more variables than districts, each domain enters the model as a few summary scores, not as hundreds of raw columns.

*Backup, if asked:* 575 columns against 14 to 87 districts per country cannot be fitted column by column, and the columns are heavily correlated within a source. Grouping them by what they measure and taking principal components inside each group keeps most of each group's variance, treats a source with 64 columns and a source with 3 as one construct each, and keeps the outcome out of the representation: the rotations are learned from the training districts' predictors alone, so nothing is selected on the biomarker. The grouping is written down: each source block brings a label, and 10 override rules (metadata/covariates/domain_overrides.csv) regroup columns by construct (infection and inflammation, immunisation, anaemia, child growth, infant feeding, supplementation, water and sanitation, education, child mortality, assets), so that Malaria Atlas, IHME, DHS and MICS versions of the same thing are weighed together. The headline model uses 383 of the 575 (no DHS microdata).

### Slide 7. Which surveys did we use to train the models? (Sonja, 0:25, 59 words)

**Visual.** The four countries on a map and a table: survey, districts with blood samples as n of N, nutrients measured.

> The models learn from four national micronutrient surveys: The Gambia 2018, Ghana 2017, Sierra Leone 2013 and Malawi 2015 to 16. Together, 206 districts where blood was drawn: 30 of The Gambia’s 37 districts, 75 of Ghana’s 260, all 14 of Sierra Leone’s, and 87 of Malawi’s 243 areas. These data exist because these teams collected them. Thank you.

*Backup, if asked:* outcome definitions, adjustments and cut-offs are in the appendix. Sierra Leone’s child file holds only anaemic children (532 of 654 assayed). Malawi’s areas are Traditional Authorities. Ghana’s survey reached 90 clusters in its 75 districts. Survey partners for Ghana’s 2017 survey: University of Ghana, GroundWork, University of Wisconsin-Madison, KEMRI-Wellcome Trust, with UNICEF and Global Affairs Canada.

## Methods

### Slide 8. How do we get honest estimates of model performance? (Andrew, 0:40, 82 words)

**Visual.** Cross-validation schematic (five rounds) and three maps: one district in five hidden, a whole region hidden, a whole country hidden. Orange is the held-out test data.

> Thank you, Sonja. Every number I show comes from data the model never saw. On the left, cross-validation: we split the surveyed districts into five parts, fit on four, predict the fifth, and rotate until every district has been predicted without its own data; then we repeat with ten random splits. Next, the same with a whole region held out. And a whole country: Ghana predicted from the other three countries. Every score is on the orange districts, the held-out test data.

*Backup, if asked:* folds are cut at the district or above, never through a district’s respondents; component weights and every other choice are learned inside the training folds. District folds are random, not grouped in space; the region hold-out is the spatially grouped test. The middle map is the benchmark’s first fold draw (make_folds_v2, rep 1).

### Slide 9. Why not a more complex machine-learning model? (Andrew, 0:45, 90 words)

**Visual.** Three lines of text; below, six methods compared under the three hold-outs (the index, a SuperLearner of twelve, random forest, boosting, lasso, elastic net), post-fix.

> The obvious question for the machine-learning people here: why such a simple model? This is a big-data project in predictors, 575 layers, but a small-data project in places: 14 to 87 districts per country, each measured in about one community. We compared 14 methods on the same held-out districts, including random forests, boosting and a SuperLearner ensemble that combines twelve. None did more than slightly better on average across nutrients, and none was better everywhere. So we use the simplest one, a weighted sum of principal components, with nothing tuned.

*Backup, if asked:* average status, the post-fix SuperLearner runs of 28 September (predictors rank-normalised within country, the production default; one district in five and a region from the NS-01 runs, a country from its own run), every method on the same folds. One district in five hidden: index 0.39, SuperLearner 0.37, random forest 0.38, boosting 0.37, elastic net 0.32, lasso 0.32. A region hidden: 0.39, 0.37, 0.41, 0.39, 0.28, 0.25. A country hidden: 0.28, 0.28, 0.24, 0.19, 0.29, 0.29. Per outcome the winner changes; folate is the index’s worst case (0.27 behind the best method with the country hidden). On the share deficient the index leads in all three: 0.30 against at most 0.25 inside a country, 0.28 against at most 0.24 with a region hidden, 0.17 against at most 0.16 with the country hidden. A sample-size finding, not a claim about machine learning in general.

## National prevalence

### Slide 10. But can proxy models recover national-level prevalence estimates? (Andrew, 0:45, 114 words)

**Visual.** Left: 15 country-outcomes, each country predicted from national indicators with its own survey left out, the gap labelled. Right: typical error in WHO's national survey panels, model against the average of the other countries.

> A question from our January meeting: can proxy models give a national prevalence? Inside a surveyed country the survey already gives it. Without the country’s survey, not reliably. On the left, every outcome we can test, each of our countries predicted from national indicators with its own survey left out: close for some, but 28 points too high for Sierra Leone’s children’s vitamin A, 38 too low for its women’s folate, 34 too low for Malawi’s zinc. On the right, WHO’s national survey data from 21 to 69 countries: the model beats the average of the other countries only for children’s vitamin A and folate. Levels do not travel between surveys; district order does.

*Backup, if asked:* national track: WHO VMNIS national survey panel with World Bank indicators; a SuperLearner over mean, ridge, lasso, elastic net and random forest, folds grouped by country (national_vmnis_loco_sl.csv, _pred.csv, national_levels_sl.csv; scripts/covariates/19b; figure script 37). Left panel: vitamin A at each survey’s year against the survey’s own national figure; folate, B12 and zinc against the country’s own VMNIS row. Median miss 7 points, largest 38. Iron has no national panel. The January version compared national predictions made inside each country, which match by construction. Why levels do not travel: assays, cut-offs, inflammation adjustment and RBP-to-retinol calibration differ between surveys (ferritin levels differ several-fold).

## District prevalence

### Slide 11. How far off is the predicted prevalence, outcome by outcome? (Andrew, 0:40, 69 words)

**Visual.** Typical error in each district's predicted prevalence, each district hidden, for 18 combinations: the model calibrated for level, the survey's regional average, and the national estimate (X).

> Same combinations, now the level: how many percentage points each district's predicted prevalence is from the survey's, each district held out. At district level nobody has a good level: the survey's own regional averages are no better than one national number for every district. The model, shrunk for this job, is slightly better than both, and that is all. It is built to order districts, not to measure them.

*Backup, if asked:* population-weighted mean absolute error of the held-out district prevalence (benchmarks_v2_cells.csv, prevalence target, post-fix). Means over 18 combinations: model (calibrated for level) 10.6 points, survey's regional average 12.0, national estimate 11.7. Read the error against the outcome: 10 points on an outcome at 40% differs from 3 points on one at 4%, which is why rare outcomes look accurate here. The ranking uses the uncalibrated index; any level uses the calibrated one.

Added 28 September: about half of these errors is the survey's own sampling noise (NZ-01: 47 per cent of the squared error over all 27 combinations, 60 over the 22 main ones; typical error against the true prevalence about 14 points rather than 19, root mean square). In a country with no survey, district levels anchored to the national figure are no closer than the national figure itself (LV-02: tie in 11 of 22 combinations at every sample size; external 11 of 29).

## Ranking

### Slide 12. How well does it rank districts, nutrient by nutrient? (Andrew, 0:50, 124 words)

**Visual.** One panel, per nutrient and country: inside the country (circle, 95% interval), whole country held out (triangle), the survey's regional figures scored fairly (diamond). National prevalence in brackets. Zinc (Malawi women) has no held-out test.

> Here is every nutrient and country, with a 95 per cent interval, and the national prevalence in brackets. Circles: inside a surveyed country. Triangles: the whole country held out. Diamonds: the survey’s own regional figures, which inside a country do about as well as the model. B12 is the most mappable: about 70 of every 100 pairs right, at home and across borders, probably because it follows animal-source foods. Vitamin A works well in The Gambia, less in Ghana. Iron works in most places, not all: across borders, Malawi’s children are no better than a coin. Folate works inside a country but not across borders, most likely because the surveys measured it differently. Zinc, measured only in Malawi, is no better than a coin.

*Backup, if asked:* measurable combinations only; Sierra Leone is not shown (14 districts leave too few to train on inside the country). In-country intervals clear 50 per cent in 11 of 14; held out in 9 of 15 (all measurable with a held-out test). Averages over the 14 shown: 66 per cent inside a country (65.6), 63 held out (the 13 with a held-out test), the survey’s regional figures 65 (65.3; model ahead in 9 of 14). Diamonds: every district in a region gets the figure of the region’s other surveyed districts, so two districts in the same region count as a coin toss (scripts/policy_deck/39_fair_regional_pairs.py). The earlier scoring, each district left out of its own region’s figure, gave 61 because it reverses two surveyed districts of the same region by construction. The regional figures need the survey; the model does not. On the 17 per cent of pairs the survey clearly separates (difference larger than 1.96 standard errors), the model gets 81.0 in 100, the regional figures 80.5, the neighbour map 81.5: every method does better on easy pairs (SEP-01, sep01_summary.csv). After the fixes Ghana’s women’s vitamin A dropped out (prevalence 1.6 per cent, under the 2 per cent screen) and, with the post-fix variance components, so did Malawi’s children’s vitamin A (ceiling 0.15, under the 0.30 screen); Malawi’s women’s zinc came in (ceiling 0.37) and has no held-out test, zinc being measured in one country only. Folate: immunoassay in Ghana and Sierra Leone, microbiologic assay in Malawi.

### Slide 13. What does that accuracy mean on the ground? (Andrew, 0:50, 117 words)

**Visual.** Three icon rows: pairs in the survey's order (model, regional figures, no survey); with a budget per person, rank by the model; with a budget in districts, start with the most populous. ALTERNATIVE to the next slide: keep one.

> What does that accuracy mean on the ground? Take any two districts: the model orders them as the survey does 66 times in 100. The survey's own regional figures do about as well, 65, but only where there is a survey; in a country with none, the model still gets 62. For a programme, it depends on the budget. If it is per person, say supplements, rank by the model: covering a fifth of the people in its highest-rate districts reaches 31 per cent of the deficient, against 20 at random. If the budget is a number of districts, start with the most populous: that fifth holds almost half the deficient people, and the model adds nothing.

*Backup, if asked:* pairs, post-fix, 14 in-country combinations (v3_percell_pairs.csv, v6_fair_regional_pairs.csv): model 66, the survey's regional figures 65 scored as a planner would use them (model ahead in 9 of 14); held out 62 (15). Budget rules (TC-01, tc01_targeting_by_cases_summary.csv; appendix figure): a fifth of the people in the model's highest-rate districts reaches 31 per cent of the deficient inside a surveyed country (regional figures 28, random 20, perfect 47) and 25 with no survey; a fifth of the districts: most populous 49 and 45 per cent, model expected cases 49 and 44, model rate 25 and 18, perfect 63 and 60. The earlier line (the worst fifth by rate holds 25 per cent) answered a question no programme asks. The next slide shows the budget rows as a figure; keep one of the two.

### Slide 14. Which districts should a programme go to first? (Andrew, 0:45, 90 words)

**Visual.** The targeting figure: share of deficient people reached under each rule, inside a surveyed country and with no survey, for a budget of districts and a budget of people. ALTERNATIVE to the previous slide's budget rows: keep one.

> Which districts should a programme go to first? It depends on the budget. Top half: a budget of districts. Ranking by the model’s rates reaches only a quarter of the deficient people; simply taking the most populous fifth reaches almost half, and the model adds nothing to that. Bottom half: a budget per person, say supplements. Covering a fifth of the people in the model’s highest-rate districts reaches 31 per cent of the deficient inside a surveyed country, and 25 in a country with no survey, against 20 at random.

*Backup, if asked:* the same result as the budget rows of the previous slide, as a figure; keep one of the two slides. TC-01 (scripts/protocol_v2/74, tc01_targeting_by_cases_summary.csv), post-fix measurable combinations (14 in-country, 15 held out). Budget of districts: most populous 49 and 45 per cent, model expected cases 49 and 44, model rate 25 and 18, perfect 63 and 60. Budget of people: model rate 31 and 25, the survey’s regional figures 28 (in-country only), perfect 47 and 47, random 20. Grey dots: each combination.

### Slide 15. Where do the young children live? (Andrew, 0:35, 84 words)

**Visual.** Ghana, children's iron: left, districts coloured by the model's predicted priority (highest-rate third: 33% of districts, 63% of the land); right, the same districts as circles sized by young children (that third holds 27%; the most populous fifth 41%).

> Why does going where the people are win when the budget is a number of districts? Here is Ghana. On the left, the model’s highest-rate third of districts: most of the north, 63 per cent of the land. On the right, the same districts sized by the number of young children. That third holds 27 per cent of them; the 52 most populous districts hold 41 per cent. Choose districts by population when you fund districts; rank by rate when you fund per person.

*Backup, if asked:* Ghana, children’s iron: the in-country model fitted on the 75 surveyed districts and applied to all 260 (prevalence target; Spearman 0.95 against the average-status map in the appendix); children aged 6-59 months from admin2_population.rds; circles laid out so they do not overlap (script 48, ghana_child_iron_priority_all_districts.csv). Across the 14 in-country combinations, measured against survey prevalence in surveyed districts, a fifth of districts chosen by population reaches 49 per cent of deficient people and by the model’s rate 25; per person, the model’s rate reaches 31 per cent with a fifth of the people (TC-01). Do not quote case shares computed from the model’s own predicted rates: they overstate the spread between districts (the highest-rate fifth would seem to hold 43 per cent of expected deficient children), the same level problem as on the national slide.

## Out of country

### Slide 16. Can it rank Ghana's districts without Ghana's survey? (Andrew, 0:50, 99 words)

**Visual.** Ghana's children's iron: the survey's ranking and the model's with no Ghanaian data. Bars: each country held out in turn.

> Now the test that matters for most of this region: we removed Ghana’s survey entirely. The model learned only from The Gambia, Sierra Leone and Malawi. On the left, the survey’s ranking of districts for children’s iron; next to it, the model’s, with no Ghanaian data at all. It puts 69 of every 100 pairs of Ghana’s districts in the survey’s order. This is one of our better cases. On the right, each country held out in turn, all nutrients: 70 for The Gambia, 62 for Ghana, 57 for Malawi, 56 for Sierra Leone. A shortlist, not a measurement.

*Backup, if asked:* correlation 0.54 held out against 0.53 in-country for this outcome, within each other’s intervals. Why about as good as the model that saw Ghana’s survey: it learned from 131 districts in three countries instead of about 60 in Ghana; over all six of Ghana’s outcomes it does a little worse, 62 pairs in 100 against 65. The survey’s own regional figures, which need the survey, put 63 in 100 in order (scored as a planner would use them; 61 the older way, which under-rates them). Ghana held out does better than Ghana in-country for children’s iron, children’s vitamin A and B12, worse for women’s vitamin A, iron and folate.

### Slide 17. What about countries with no survey? (Andrew, 0:55, 118 words)

**Visual.** Left: Côte d'Ivoire's children's iron ranking from the climate-and-soil model. Right: the six countries whose WHO regional survey results checked the model.

> Most countries in this region have no recent biomarker survey. Côte d’Ivoire is one. Trained on our four countries, the climate-and-soil version of the model ranks its 33 districts for children’s iron; the northern savanna comes out worst. Can you trust a map like this? WHO’s micronutrients database holds regional results from national surveys in Zambia, Ethiopia, Sudan, Nigeria, Pakistan and India. None were used to build the model. In Africa its ranking of their regions pointed the right way in all 12 tests of average status and 18 of 21 tests of prevalence; in South Asia, 7 of 8. And Côte d’Ivoire’s own 2007 survey ranks its nine zones for B12 almost exactly as the model does.

*Backup, if asked:* XV-01/02, rerun on the fixed targets: climate-and-soil index at admin-1; Africa average status 0.40 (country-block null 0.20), prevalence 0.34; off the continent 0.38. Iron points the right way in all 14 iron tests across the six countries (weakest Ethiopia, 0.01). Regions, not districts; not the pre-registered test, which needs a new country’s own microdata at district level. Côte d’Ivoire 2007 B12: Spearman 0.95 over nine zones. The Côte d’Ivoire map is post-fix (ranking rebuilt 28 September, 05:07); the survey fixes barely moved it (rank agreement 0.999 with the 8 September map; the same seven worst districts, led by Tchologo, Bounkani and Bagoué).

### Slide 18. Is Côte d'Ivoire like the places we learned from? (Andrew, 0:35, 86 words)

**Visual.** One row per country: each district's distance from the climate and soil conditions of the surveyed districts, with the edge of the training range. Côte d'Ivoire 33 of 33 inside; each surveyed country, held out, below it.

> Can we trust a map of a country with no survey? One check we can make: is it like the places the model learned from? Each dot is a district, placed by how far its climate and soil are from the nearest surveyed district; the line is the edge of what the model has seen. Every Ivorian district sits well inside. This checks the inputs, not the accuracy: Malawi is mostly outside and Sierra Leone mostly inside, yet the model orders their districts about equally well.

*Backup, if asked:* area of applicability (Meyer and Pebesma 2021): unweighted dissimilarity on the 94 climate and soil columns the Côte d’Ivoire model uses, pooled raw values standardised on the training districts, threshold from country-blocked cross-validation (CAST’s rule: Q3 + 1.5 IQR, capped at the maximum), computed in base R (script 47, civ_applicability.csv). Côte d’Ivoire 33 of 33 inside (the furthest, Abidjan, at 0.51 of the edge). Each surveyed country held out: Ghana 74 of 75, The Gambia 24 of 30, Sierra Leone 9 of 14, Malawi 19 of 87 (8 under the older boxplot rule; half of Malawi’s districts are 1.14 of the edge or beyond). Pairs in the survey’s order with the country held out: Ghana 62, The Gambia 70, Sierra Leone 56, Malawi 57. It is not a probability of being right, and four countries cannot calibrate it against accuracy.

## Variable importance

### Slide 19. What variables drive model predictions (Andrew, 1:00, 146 words)

**Visual.** Three findings on the left (climate and soil carry the cross-border ranking: dropping either makes it worse in all four countries; climate, soil and a satellite imagery summary are in the top five of all 22 cross-country and 21 of 24 single-country models; only B12 has a map of its own). Right, children's iron: top-ten domains in Ghana's, Malawi's and the four-country model, and the four-country model's top ten layers with their direction.

> What does the model rely on? On the right, children’s iron: the top ten families of layers in Ghana’s model, Malawi’s, and the four-country model we use for a new country, and that model’s top ten layers. The strongest evidence is what happens when we take a family of layers away: without climate, or without soil, the ranking of a country the model has never seen gets worse, in all four countries, and climate alone does about as well as all 575 layers. Climate, soil and a satellite imagery summary are also in the top five of almost every model we fitted. B12 has markers of its own: meat and fish in children’s diets, and household wealth. The iron and vitamin A maps mostly mark the same poorer districts: swap their weights and the rankings barely change. These are markers of where deficiency is, not causes.

*Backup, if asked:* post-fix importance and ablation tables (index_importance_domains.csv, index_importance_columns.csv, domain_ablation_loco_summary.csv). Pooled weight shares, mean over six outcomes: climate 18 per cent, satellite embedding 16, soil 13, greenness 6; climate + soil + satellite + greenness 45 to 59 per cent by outcome. Country held out: climate alone 0.33, soil alone 0.33, all layers 0.30; removing climate or soil costs 0.03, anaemia 0.02, crops and livestock about 0.01 each. Layers in the top 20 of four outcomes, same sign every time: malaria mortality, child overweight, relative staple price, plant productivity (less deficiency), grassland cover (more). Malaria marks LESS deficiency for iron and vitamin A, most likely because inflammation raises ferritin and lowers RBP differently from deficiency: a feature of the biomarkers. Schooling, water and sanitation, household diets and infant feeding help inside a country but make cross-border ranking slightly worse. Fertility and reproductive health are in no top list. Drop-a-domain test (domain_ablation_loco.csv, average status, 22 cross-country combinations): dropping climate costs 0.030 (The Gambia 0.030, Ghana 0.014, Malawi 0.030, Sierra Leone 0.045), soil 0.029 (0.043, 0.014, 0.023, 0.040); climate alone minus all layers +0.028 (better in 3 of 4 countries), soil alone +0.029 (2 of 4). The satellite embedding carries weight but adds almost nothing across borders: dropping it costs 0.002, and alone it scores 0.10 lower. Top five by weight: all 22 cross-country and all 6 pooled models, and 21 of 24 single-country models (not Malawi children’s iron, Malawi women’s zinc or Sierra Leone women’s B12). Weight shares largely follow how many layers a family has (climate, soil and the embedding are 44 per cent of the model’s 383 columns and get 47 per cent of the weight), so weight alone is not evidence of importance; the drop-a-domain test is. Nutrient specificity (NX-01, nx01_summary.csv): with the country held out, each outcome scored with another nutrient’s weights. B12’s own weights beat every other nutrient’s in Ghana and Malawi (own 0.60 and 0.44, best other 0.54 and 0.31), not in Sierra Leone; iron is specific in 2 of 8 combinations, vitamin A in 0 of 8, folate in 1 of 3; one generic index over all outcomes scores 0.29 against 0.30 for each outcome’s own. Children’s iron figure: weight shares from index_importance_domains.csv (level target), layers from index_importance_columns.csv, pooled fit; every layer keeps its sign in all four leave-one-country-out fits. The full findings slide is in the appendix.

## Policy relevance and next steps

### Slide 20. Which deficiencies can it map? (Andrew, 0:45, 93 words)

**Visual.** Traffic-light table by nutrient: inside a surveyed country, whole country held out, four African countries in WHO data, and whether to use it. One line underneath: only the B12 map is specific to its nutrient.

> Putting it together. B12 and iron can be ranked from public data: inside a country, across borders, and in the external check. Vitamin A ranks well inside our countries but was the weakest externally, so not yet. Folate works only within a country. Zinc cannot be mapped: the survey itself shows no real difference between districts. One caution, in the line at the bottom: only the B12 map is its own. Rank iron districts with the vitamin A weights and you do about as well, because both maps mark the same poorer districts.

*Backup, if asked:* external by nutrient, African deposits, post-fix: B12 0.66 (4 of 4 tests), iron 0.44 (11 of 11), folate 0.28 (3 of 4), vitamin A 0.24 (12 of 14; both clear misses are vitamin A). Vitamin A is the strongest outcome in India and Pakistan. The table’s numbers are correlations (0 = chance); measurable combinations only. Nutrient specificity (NX-01): with the country held out, B12’s own weights beat every other nutrient’s in 2 of 3 countries; iron 2 of 8, vitamin A 0 of 8, folate 1 of 3; the iron indices rank vitamin A districts at least as well as vitamin A’s own (0.39 and 0.34 against 0.33). Zinc: nested variance components put Malawi children’s zinc at district share 0.000 and cluster share 0.172, the largest cluster share in the study (explore/out/22_zinc_variance_components.csv): the survey’s zinc differences sit between communities, not districts.

### Slide 21. What does this mean for policy? The vitamin B12 model (Andrew, 0:50, 128 words)

**Visual.** Four cards: 70 to 75 of 100 pairs (65 to 71 held out, average B12 level); about 2x capture inside a surveyed country; about 9 points typical prevalence error; 4 of 4 WHO checks.

> What does this mean for a programme? Take our best model, women’s B12. It puts 70 to 75 of every 100 district pairs in the survey’s order of average B12 level, and 65 to 71 with the country’s own survey removed. Inside a surveyed country, the fifth of districts it ranks worst hold 39 to 46 per cent of deficient women: about twice a random choice. Across borders the ranking of average level holds, but the share deficient does not, so there it needs a survey to anchor it. Its district prevalences are typically about 9 points from the survey’s, about as close as the survey’s own regional figures. And in WHO survey data from Zambia, Ethiopia and Nigeria it pointed the right way in all four checks.

*Backup, if asked:* pairs on average B12 level (v3_percell_pairs.csv, post-fix): Ghana 70 in-country, 71 held out; Malawi 75 and 65 (Sierra Leone 71 held out, but its B12 prevalence, 0.6 per cent, is under the 2 per cent screen). Worst fifth, inside a surveyed country (nce_targeting_metrics.csv, post-fix): Ghana 46 per cent (regional figures 45, perfect 65), Malawi 39 (27, 65). With the country’s own survey removed the worst fifth holds 12 per cent (Ghana) and 11 (Malawi) of deficient women, less than a random fifth, and its districts are less deficient than average (Ghana 4.9 against 7.8 per cent): on the share deficient B12 transports at 0.24 and 0.22, against 0.60 and 0.44 on average level. Prevalence error, each district hidden, population-weighted, model calibrated for level (benchmarks_v2_cells.csv, post-fix): Ghana 9.3 points (regional figures 9.1), Malawi 9.1 (8.9). The 90 per cent bands were rerun on 28 September (median half-width 23 points over all combinations), but script 66 still anchors on a national-estimates file that predates the survey fixes (Malawi B12 anchor 2.3 against 10.9 per cent), so no band is quoted. B12 is the one nutrient whose map is its own: its weights beat every other nutrient’s in Ghana and Malawi (NX-01). External B12, post-fix: 4 of 4 in Africa (mean 0.66); Côte d’Ivoire’s nine 2007 survey zones in almost the same order (0.95); Pakistan’s B12 is a miss.

### Slide 22. Next steps (Andrew, 0:30, 69 words)

**Visual.** Three asks as cards, the thesis banner and the QR code.

> Three asks. If your country has a micronutrient survey, recent, planned or in a drawer, share it with the pool: on average, each survey we added improved the ranking for the countries held out. If you are planning one, let’s design it together. And tell us how accurate a district ranking must be before you would use it. Modelling extends a survey’s reach. It does not replace the survey.

*Backup, if asked:* Pooling curve: 0.20 with one training survey, 0.26 with two, 0.30 with three; the gain so far is in The Gambia and Ghana, not yet Malawi or Sierra Leone. Two candidate models were written down on 3 September 2026, so a fifth survey is a real test. Other surveys help a country with no survey, not a surveyed one: adding the other countries’ weights to a surveyed country’s own made its ranking slightly worse (BO-01: -0.02, better in 8 of 18). Optional line for Sonja: next, with MIMI and GBD, we compare predictors framework by framework.

## Conclusions

### Slide 23. Conclusions (Andrew, 1:00, 145 words)

**Visual.** The five questions answered, each with its number, and the QR code.

> To answer our five questions. Can public data tell us how common a deficiency is? Not on their own: without a survey, national estimates missed by up to 38 points, and district levels need the survey’s national figure. Can they rank a country’s districts? Yes, modestly: two pairs in three in the survey’s order, about as well as the survey’s own regional figures. In a country with no survey? Mostly: six pairs in ten, and the right direction in 37 of 41 checks in six more countries. Which data? Climate and soil carry the cross-border ranking; and only B12 has a map of its own. And for programmes: with a budget per person, rank by the model; with a budget in districts, start with the most populous. It does not replace a survey. Everything is on the dashboard; the QR code is on the screen.

*Backup, if asked:* national levels without the country’s survey: median miss 7 points, largest 38 (national_levels_sl.csv, national_vmnis_loco_sl_pred.csv). Pairs, post-fix, 14 in-country combinations: model 66, survey’s regional figures 65 (scored fairly; model ahead in 9 of 14); held out 62 over 15. External (WHO VMNIS regions): Africa 12 of 12 positive at level (0.40), prevalence 18 of 21 (0.34), South Asia 7 of 8. Nutrient specificity: NX-01. Targeting (TC-01): a fifth of the people covered in the model’s highest-rate districts reaches 31 per cent of the deficient inside a surveyed country (regional figures 28, random 20, perfect 47) and 25 with no survey; with a budget of a fifth of the districts, the most populous reach 49 and 45 per cent, and the model’s expected cases add nothing. District levels in a new country: the model’s ranking adds nothing to the national figure (LV-02). About half of the district prevalence error is the survey’s own sampling noise (NZ-01).

# 4. Where each number stands

| Number | Slide | Status | Source |
|---|---|---|---|
| Pairs: in-country 66, regional figures 65, held out 62 | 12, 13, 20 | Post-fix | `policy_deck/v3_percell_pairs.csv`, `v6_fair_regional_pairs.csv` |
| National levels: median miss 7 points, largest 38 | 10 | Post-fix | `national_levels_sl.csv`, `national_vmnis_loco_sl_pred.csv` |
| District prevalence error: model 10.6, regional 12.0, national 11.7 points | 11 | Post-fix | `benchmarks_v2_cells.csv` (prevalence, calibrated arm) |
| Ghana children's iron held out: 69 of 100 pairs (0.54); countries 70, 62, 57, 56 | 14 | Post-fix | script 25, `v3_percell_pairs.csv` |
| WHO external check: Africa 12 of 12 (0.40), prevalence 18 of 21 (0.34), South Asia 7 of 8 | 15, 20 | Post-fix | `external_validation/xv_transport_summary.csv` |
| Targeting: 31% against 22%; 25% of deficient people; 17% with no survey | 13 | Post-fix | `nce_targeting_metrics.csv` |
| B12: 70 to 75, 65 to 71 pairs; 39 to 46%; about 9 points; 4 of 4 | 18 | Post-fix | `v3_percell_pairs.csv`, `nce_targeting_metrics.csv`, `benchmarks_v2_cells.csv`, `xv_transport.csv` |
| Budget rules: 31% per person (random 20, regional 28), most populous fifth 49%, model cases 49 | 13, 20, appendix | Post-fix | `tc01_targeting_by_cases_summary.csv` |
| Only B12 nutrient-specific (2 of 3); iron 2 of 8; vitamin A 0 of 8 | 16, 17, 20 | Post-fix | `nx01_summary.csv` |
| Confident pairs: model 81.0, regional 80.5, neighbour map 81.5 | 12 notes, appendix | Post-fix | `sep01_summary.csv` |
| Anchored levels = national figure (tie 11 of 22) | 11 notes, appendix | Post-fix | `lv02_internal_summary.csv`, `lv02_external_summary.csv` |
| Survey noise share of prevalence error: 47% (all), 60% (22 main) | 11 notes | Post-fix | `nz01_overall.csv` |
| Ceiling 0.47, model 0.30; vitamin A severity bands 50 / 53 / 54 | Appendix | Post-fix | `variance_components_ceiling.csv` (28 Sep), `risk_category_accuracy.csv` |
| Importance findings | 16 | Post-fix | `index_importance_domains.csv`, `index_importance_columns.csv`, `domain_ablation_loco_summary.csv` |
| Method comparison figure | 9 | Post-fix | `weight_sources_raw_ns01_sl_rank_[a-d].csv`, `weight_sources_raw_v6_sl_country.csv` |
| Côte d'Ivoire ranking map | 16 | Post-fix (rebuilt 28 Sep 05:07; rank agreement 0.999 with the 8 September map) | `civ_rank_uncertainty.csv` |
| Prevalence bands (90%, CP-01) | 18 notes | **Partly pre-fix** | Rerun 28 Sep (median half-width 23 points), but script 66 still anchors on a national-estimates file from 1 September; no band is quoted |
| Appendix: the geostatistics comparison (8-9 September cluster targets) and the survey-planner line in the planning slide's notes | Appendix | **Pre-fix** | Every other appendix figure was redrawn on post-fix data on 28 September (scripts 43 to 46); the geostatistics rerun needs 4 to 6 hours |

# 5. If the time is tighter

**About 10 minutes, 16 slides.** Drop the outline (3); fold the 575 variables (6) into slide 5's last sentence; keep one of the two targeting slides (13, 14), or fold it into slide 12's script (one sentence: 66 in 100, the regional figures 65, no survey 62); fold "Which deficiencies can it map?" (18) into the first sentence of the B12 slide. About 300 words saved: roughly 10.5 minutes.

**The Forum's 8.** Also merge the framework and its coverage (4, 5) as one slide with a click, the surveys into the honest-estimates slide (7 into 8), national and district prevalence (10, 11), Ghana held out and countries with no survey (15, 16), and end on the conclusions with the asks as its last line (20 into 21). Denser; keep 16 unless the chair enforces the guideline.

# 6. Appendix, in deck order

1. How does the model work, and how was it tested?
2. How does the best-performing model work?
3. Both a big-data opportunity and a small-data problem
4. How do we choose among algorithms?
5. How well does it rank districts inside a surveyed country?
6. Malawi: survey, held-out prediction, every area ranked, and uncertainty
7. What happens when a whole region is held out?
8. Where do the rankings work?
9. How much ranking accuracy is lost without the country's own survey?
10. Accuracy by nutrient, against the best possible
11. How much of what a survey can see does the model recover?
12. How do proxy models compare with the survey's own information?
13. How far off is the predicted mean concentration, outcome by outcome?
14. Is it just geography? Comparison with DHS-style geostatistics
15. Does it match the survey at the regional level?
16. Could a small national survey plus the model replace a district survey?
17. How can it help plan the next survey?
18. How does it compare with other methods, including machine learning?
19. Which data sources contribute most to model prediction performance?
20. What travels across borders, and what does not?
21. What can we learn about the consistently most important predictor domains?
22. What can we learn about unused predictor domains?
23. What drives the transport models?
24. Top 20 individual predictors for child iron
25. What are the weights behind one of the best rankings? Malawi, women's B12
26. Causal driver or surrogate marker? Fish, anemia and B12 in Malawi
27. Which outcomes and cut-offs does the analysis use?
28. Which deficiencies do we predict?
29. How are the layers matched to districts in space and time?
30. Does targeting reach more deficient people?
31. Does the model put districts in the right severity band?
32. Accuracy pooled over all nutrients, in district pairs
33. What the models can and cannot do yet
34. Every combination of training countries
35. Does the model track the survey, district by district?
36. The 2007 Côte d'Ivoire survey and the B12 prediction
37. The live dashboard

# 7. Open items before 1 October

1. **Main talk:** every figure is now post-fix. The Côte d'Ivoire map was rebuilt on 28 September and barely moved (rank agreement 0.999).
2. **Appendix:** every figure was redrawn on post-fix data on 28 September except the geostatistics comparison (a 4 to 6 hour rerun) and the survey-planner line in the planning notes. The prevalence bands need script 66's national anchor file rebuilt on the fixed targets (a dashboard input too, so coordinate with the dashboard session).
3. **Slide 6** (Andrew's): the title overlaps the chart in the render; check in PowerPoint.
4. **Sonja's slides**: her introduction and framework-coverage notes were shortened (to about 40 and 45 seconds); she should read them.

# 8. Decisions

Made 28-29 September: keep 20 slides, plus the two figures added on 29 September; keep the five-question outline; show both targeting versions (13 and 14) for Andrew to choose; the nutrient finding stands (the session covers all micronutrients, not iron); slide 9 redrawn on the post-fix runs.

Still open:

1. Slide 13 (the budget rule in words) or slide 14 (the same result as a figure)?
