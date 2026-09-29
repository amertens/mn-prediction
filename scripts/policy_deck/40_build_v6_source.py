"""Build docs/slides/MNF15-talk-2026-09-v6.qmd from the v5 qmd: the 15-minute order
Andrew asked for on 28 September (motivation, data assembly, methods, national
prevalence, district prevalence, ranking, out of country, variable importance, policy
and next steps, conclusions), about 13 minutes spoken.

    python scripts/policy_deck/40_build_v6_source.py
    bash scripts/render_deck.sh docs/slides/MNF15-talk-2026-09-v6.qmd --forum --no-label --text 18,12 \\
         --concept docs/slides/MNF15-talk-2026-09-v6.concept.yaml
    python scripts/policy_deck/41_merge_v6_slides.py

Changes from v5:
  - Ghana held out moves after the ranking slides (the out-of-country section)
  - "Both a big-data opportunity...", "How do we choose among algorithms?" and "How much of
    what a survey can see..." go to the appendix (script 41)
  - "What variables drive model predictions" is the figure plus the three strongest
    findings; the full findings slide goes to the appendix
  - "Conclusions" and "What have we learned..." become one closing slide that mirrors the
    outline; "Next steps" comes before it; v4-ANM's old closing slide is dropped
  - the survey's regional average is scored fairly everywhere it is compared (script 39):
    66 against 65 pairs inside a country, not 66 against 61
  - B12 district prevalence error is the post-fix figure (about 9 points), not the
    8 September conformal median (6)
  - notes of the introduction and the framework-coverage slide trimmed
"""
import os
import re

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
V5 = ROOT + "docs/slides/MNF15-talk-2026-09-v5.qmd"
V6 = ROOT + "docs/slides/MNF15-talk-2026-09-v6.qmd"
F5 = "../../results/figures/mnf15_v5/"
F6 = "../../results/figures/mnf15_v6/"

src = open(V5, encoding="utf-8").read()
blocks = {}
for m in re.finditer(r"^## ([^\n]+)\n(.*?)(?=^## |^# |\Z)", src, re.S | re.M):
    blocks[m.group(1).strip()] = m.group(2).strip() + "\n"


def notes(text):
    return "::: {.notes}\n" + text.strip() + "\n:::\n"


def slide(title, body):
    return f"## {title}\n\n{body.strip()}\n\n"


def keep(title, subs=(), new_title=None):
    """A v5 slide, with exact text replacements (each must be found)."""
    b = blocks[title]
    for a, c in subs:
        assert a in b, (title, a[:70])
        b = b.replace(a, c)
    return slide(new_title or title, b)


header = src[:src.index("```{r setup}")].replace("VERSION 5", "VERSION 6")
setup = '''```{r setup}
# MNF 2026 joint talk, VERSION 6 (28 September 2026). v1 to v5 are untouched.
# v6 = v5 in the 15-minute order Andrew asked for (motivation, data assembly, methods,
# national prevalence, district prevalence, ranking, out of country, variable importance,
# policy and next steps, conclusions), about 13 minutes spoken; three slides to the
# appendix, two pairs merged; the survey's regional average scored fairly (script 39).
# Render: python scripts/policy_deck/40_build_v6_source.py; bash scripts/render_deck.sh
#   docs/slides/MNF15-talk-2026-09-v6.qmd --forum --no-label --text 18,12
#   --concept docs/slides/MNF15-talk-2026-09-v6.concept.yaml; python scripts/policy_deck/41_merge_v6_slides.py
suppressPackageStartupMessages({library(dplyr); library(knitr)})
```

'''
out = [header, setup]

# ------------------------------------------------------------ motivation
out.append(slide("Introduction and study objectives", notes('''
SONJA, about 40 seconds (92 words).

Information on how common vitamin and mineral deficiencies are is limited. The map shows how limited for sub-Saharan Africa: 17 of 47 countries have no national biomarker survey on record, and only 10 have one from 2015 or later. A survey gives national and regional numbers, while programmes are planned district by district. So we brought together nutritionists, epidemiologists and data scientists to ask which approaches and data could help. Our objectives: to see whether machine learning on aggregated proxy data can estimate deficiency prevalence, nationally and by district, and find the districts at greatest risk.

Backup, if asked: map and bars from the WHO VMNIS Micronutrients Database (export of 25 February 2025), nationally representative surveys measuring ferritin, retinol or RBP, zinc, folate or vitamin B12 (scripts/policy_deck/33_mnf15_v4_vmnis_coverage.R); surveys not yet deposited are not shown, and all four of our surveys are in it. For zinc, only five countries have a survey from 2015 or later. Worldwide, in low- and middle-income countries over 1988 to 2018: vitamin A data in 77 countries, iron 53, folate 24, zinc 21, B12 7 (Brown, Moore, Hess et al., Am J Clin Nutr 2021, CC BY 4.0). The global estimate of 56 per cent of young children and 69 per cent of women with at least one deficiency rests on 24 surveys in 22 countries (Stevens et al., Lancet Glob Health 2022). Fortification is set nationally (Nyumuah et al., Food Nutr Bull 2012); coverage of fortified foods is lower among poor and rural households (Aaron et al., J Nutr 2017). Why model at all: by the conventional survey bands for the coefficient of variation, 86 per cent of the 1,350 direct district prevalence estimates in our four surveys have a CV above 33.3 per cent and 5 per cent are under 16.6 per cent (median CV: Ghana 100 per cent, Malawi 96) (explore/out/19_cv_bands.csv). Two guards: these are the survey's own direct district figures, not the model's, so this says modelling is needed, not that the model's outputs are unreliable; and CV (standard error over prevalence) penalises rare outcomes by construction, so the estimates that pass are simply the high-prevalence ones.
''')))
out.append(slide("Talk outline", notes('''
SONJA, about 20 seconds (50 words).

We'll try to answer five questions. Can public data tell us how common a deficiency is, nationally or by district? Can they rank a country's districts? Does that work where there has never been a survey? Which data matter? And how could programmes and survey planners use it?
''')))

# ------------------------------------------------------------ data assembly
out.append(keep("Conceptual framework of iron deficiency"))
cov = blocks["Which of these causes can public data measure?"]
cov_backup = cov[cov.index("Backup, if asked:"):]
out.append(slide("Which of these causes can public data measure?", notes('''
SONJA, about 45 seconds (98 words).

So we went back to the framework, box by box, and asked: is there a public, district-level measure? Dark green: yes, and the data describe the place itself: climate, ecology, built-up land. Light green: yes, as a district average standing in for individuals: poverty, schooling, water and sanitation, infection. Amber: only a national figure, a modelled surface, or part of the box. Grey: nothing public. Everything here is a district average, never the person whose blood was drawn. So these layers mark where deficiency is likely; they cannot explain any one person's status.

''' + cov_backup.replace(":::", "").strip())))
out.append(keep("Which surveys did we use to train the models?"))

# ------------------------------------------------------------ methods
out.append(keep("How do we get honest estimates of model performance?"))
out.append(keep("Why not a more complex machine-learning model?", subs=[
    ("Backup, if asked: average status, SL-06 run (8 September targets; the survey fixes moved the index by less than 0.005). One district in five hidden: index 0.39, SuperLearner 0.38, random forest 0.38, boosting 0.37, elastic net 0.32, lasso 0.29. A region hidden: 0.38, 0.36, 0.39, 0.40, 0.23, 0.21. A country hidden: 0.28, 0.27, 0.25, 0.20, 0.29, 0.31. Per outcome the winner changes; folate is the index's worst case. On the share deficient the index leads inside a country (0.27 against at most 0.24) and with a region hidden (0.27 against at most 0.23); across borders the ensemble is marginally ahead (0.20 against 0.19).",
     "Backup, if asked: average status, the post-fix SuperLearner runs of 28 September (predictors rank-normalised within country, the production default; one district in five and a region from the NS-01 runs, a country from its own run), every method on the same folds. One district in five hidden: index 0.39, SuperLearner 0.37, random forest 0.38, boosting 0.37, elastic net 0.32, lasso 0.32. A region hidden: 0.39, 0.37, 0.41, 0.39, 0.28, 0.25. A country hidden: 0.28, 0.28, 0.24, 0.19, 0.29, 0.29. Per outcome the winner changes; folate is the index's worst case (0.27 behind the best method with the country hidden). On the share deficient the index leads in all three: 0.30 against at most 0.25 inside a country, 0.28 against at most 0.24 with a region hidden, 0.17 against at most 0.16 with the country hidden.")]))

# ------------------------------------------------------------ national, then district prevalence (script 41 copies the error slide)
out.append(keep("But can proxy models recover national-level prevalence estimates?"))

# ------------------------------------------------------------ ranking (script 41 copies "What does that accuracy mean on the ground?")
out.append(slide("How well does it rank districts, nutrient by nutrient?", "![](" + F6 + "v4_percell_pairs_combined.png){width=\"9.6in\"}\n\n" + notes('''
ANDREW, about 50 seconds (125 words).

Here is every nutrient and country, with a 95 per cent interval, and the national prevalence in brackets. Circles: inside a surveyed country. Triangles: the whole country held out. Diamonds: the survey's own regional figures, which inside a country do about as well as the model. B12 is the most mappable: about 70 of every 100 pairs right, at home and across borders, probably because it follows animal-source foods. Vitamin A works well in The Gambia, less in Ghana. Iron works in most places, not all: across borders, Malawi's children are no better than a coin. Folate works inside a country but not across borders, most likely because the surveys measured it differently. Zinc, measured only in Malawi, is no better than a coin.

Backup, if asked: measurable combinations only; Sierra Leone is not shown (14 districts leave too few to train on inside the country). In-country intervals clear 50 per cent in 11 of 14; held out in 9 of 15 (all measurable with a held-out test). Averages over the 14 shown: 66 per cent inside a country (65.6), 63 held out (the 13 with a held-out test), the survey's regional figures 65 (65.3; model ahead in 9 of 14). Diamonds: every district in a region gets the figure of the region's other surveyed districts, so two districts in the same region count as a coin toss (scripts/policy_deck/39_fair_regional_pairs.py). The earlier scoring, each district left out of its own region's figure, gave 61 because it reverses two surveyed districts of the same region by construction. The regional figures need the survey; the model does not. On the 17 per cent of pairs the survey clearly separates (difference larger than 1.96 standard errors), the model gets 81.0 in 100, the regional figures 80.5, the neighbour map 81.5: every method does better on easy pairs (SEP-01, sep01_summary.csv). After the fixes Ghana's women's vitamin A dropped out (prevalence 1.6 per cent, under the 2 per cent screen) and, with the post-fix variance components, so did Malawi's children's vitamin A (ceiling 0.15, under the 0.30 screen); Malawi's women's zinc came in (ceiling 0.37) and has no held-out test, zinc being measured in one country only. Folate: immunoassay in Ghana and Sierra Leone, microbiologic assay in Malawi.
''')))

out.append(slide("Which districts should a programme go to first?", "![](" + F6 + "v6_targeting_rules.png){width=\"10.4in\"}\n\n" + notes('''
ANDREW, about 45 seconds (101 words).

Which districts should a programme go to first? It depends on the budget. Top half: a budget of districts. Ranking by the model's rates reaches only a quarter of the deficient people; simply taking the most populous fifth reaches almost half, and the model adds nothing to that. Bottom half: a budget per person, say supplements. Covering a fifth of the people in the model's highest-rate districts reaches 31 per cent of the deficient inside a surveyed country, and 25 in a country with no survey, against 20 at random.

Backup, if asked: the same result as the budget rows of the previous slide, as a figure; keep one of the two slides. TC-01 (scripts/protocol_v2/74, tc01_targeting_by_cases_summary.csv), post-fix measurable combinations (14 in-country, 15 held out). Budget of districts: most populous 49 and 45 per cent, model expected cases 49 and 44, model rate 25 and 18, perfect 63 and 60. Budget of people: model rate 31 and 25, the survey's regional figures 28 (in-country only), perfect 47 and 47, random 20. Grey dots: each combination.
''')))

out.append(slide("Where do the young children live?", "![](" + F6 + "a6_targeting_cartogram.png){width=\"10.4in\"}\n\n" + notes('''
ANDREW, about 35 seconds (84 words).

Why does going where the people are win when the budget is a number of districts? Here is Ghana. On the left, the model's highest-rate third of districts: most of the north, 63 per cent of the land. On the right, the same districts sized by the number of young children. That third holds 27 per cent of them; the 52 most populous districts hold 41 per cent. Choose districts by population when you fund districts; rank by rate when you fund per person.

Backup, if asked: Ghana, children's iron: the in-country model fitted on the 75 surveyed districts and applied to all 260 (prevalence target; Spearman 0.95 against the average-status map in the appendix); children aged 6-59 months from admin2_population.rds; circles laid out so they do not overlap (script 48, ghana_child_iron_priority_all_districts.csv). Across the 14 in-country combinations, measured against survey prevalence in surveyed districts, a fifth of districts chosen by population reaches 49 per cent of deficient people and by the model's rate 25; per person, the model's rate reaches 31 per cent with a fifth of the people (TC-01). Do not quote case shares computed from the model's own predicted rates: they overstate the spread between districts (the highest-rate fifth would seem to hold 43 per cent of expected deficient children), the same level problem as on the national slide.
''')))

# ------------------------------------------------------------ out of country
out.append(keep("Can it rank Ghana's districts without Ghana's survey?", subs=[
    ("The survey's own regional averages, which need the survey, put 61 in 100 in order.",
     "The survey's own regional figures, which need the survey, put 63 in 100 in order (scored as a planner would use them; 61 the older way, which under-rates them).")]))
out.append(keep("What about countries with no survey?", subs=[
    ("The Côte d'Ivoire map is from the 8 September build until the downstream run finishes.",
     "The Côte d'Ivoire map is post-fix (ranking rebuilt 28 September, 05:07); the survey fixes barely moved it (rank agreement 0.999 with the 8 September map; the same seven worst districts, led by Tchologo, Bounkani and Bagoué).")]))
out.append(slide("Is Côte d'Ivoire like the places we learned from?", "![](" + F6 + "a6_civ_applicability.png){width=\"10.4in\"}\n\n" + notes('''
ANDREW, about 35 seconds (86 words).

Can we trust a map of a country with no survey? One check we can make: is it like the places the model learned from? Each dot is a district, placed by how far its climate and soil are from the nearest surveyed district; the line is the edge of what the model has seen. Every Ivorian district sits well inside. This checks the inputs, not the accuracy: Malawi is mostly outside and Sierra Leone mostly inside, yet the model orders their districts about equally well.

Backup, if asked: area of applicability (Meyer and Pebesma 2021): unweighted dissimilarity on the 94 climate and soil columns the Côte d'Ivoire model uses, pooled raw values standardised on the training districts, threshold from country-blocked cross-validation (CAST's rule: Q3 + 1.5 IQR, capped at the maximum), computed in base R (script 47, civ_applicability.csv). Côte d'Ivoire 33 of 33 inside (the furthest, Abidjan, at 0.51 of the edge). Each surveyed country held out: Ghana 74 of 75, The Gambia 24 of 30, Sierra Leone 9 of 14, Malawi 19 of 87 (8 under the older boxplot rule; half of Malawi's districts are 1.14 of the edge or beyond). Pairs in the survey's order with the country held out: Ghana 62, The Gambia 70, Sierra Leone 56, Malawi 57. It is not a probability of being right, and four countries cannot calibrate it against accuracy.
''')))


# ------------------------------------------------------------ variable importance (drawn by the concept builder)
imp = blocks["What can we learn about the consistently most important predictor domains?"]
imp_backup = imp[imp.index("Backup, if asked:"):].replace(":::", "").strip()
out.append(slide("What variables drive model predictions", notes('''
ANDREW, about 60 seconds (135 words).

What does the model rely on? On the right, children's iron: the top ten families of layers in Ghana's model, Malawi's, and the four-country model we use for a new country, and that model's top ten layers. The strongest evidence is what happens when we take a family of layers away: without climate, or without soil, the ranking of a country the model has never seen gets worse, in all four countries, and climate alone does about as well as all 575 layers. Climate, soil and a satellite imagery summary are also in the top five of almost every model we fitted. B12 has markers of its own: meat and fish in children's diets, and household wealth. The iron and vitamin A maps mostly mark the same poorer districts: swap their weights and the rankings barely change. These are markers of where deficiency is, not causes.

''' + imp_backup + " Drop-a-domain test (domain_ablation_loco.csv, average status, 22 cross-country combinations): dropping climate costs 0.030 (The Gambia 0.030, Ghana 0.014, Malawi 0.030, Sierra Leone 0.045), soil 0.029 (0.043, 0.014, 0.023, 0.040); climate alone minus all layers +0.028 (better in 3 of 4 countries), soil alone +0.029 (2 of 4). The satellite embedding carries weight but adds almost nothing across borders: dropping it costs 0.002, and alone it scores 0.10 lower. Top five by weight: all 22 cross-country and all 6 pooled models, and 21 of 24 single-country models (not Malawi children's iron, Malawi women's zinc or Sierra Leone women's B12). Weight shares largely follow how many layers a family has (climate, soil and the embedding are 44 per cent of the model's 383 columns and get 47 per cent of the weight), so weight alone is not evidence of importance; the drop-a-domain test is. Nutrient specificity (NX-01, nx01_summary.csv): with the country held out, each outcome scored with another nutrient's weights. B12's own weights beat every other nutrient's in Ghana and Malawi (own 0.60 and 0.44, best other 0.54 and 0.31), not in Sierra Leone; iron is specific in 2 of 8 combinations, vitamin A in 0 of 8, folate in 1 of 3; one generic index over all outcomes scores 0.29 against 0.30 for each outcome's own. Children's iron figure: weight shares from index_importance_domains.csv (level target), layers from index_importance_columns.csv, pooled fit; every layer keeps its sign in all four leave-one-country-out fits. The full findings slide is in the appendix.")))

# ------------------------------------------------------------ policy and next steps
out.append(slide("Which deficiencies can it map?", notes('''
ANDREW, about 45 seconds (104 words).

Putting it together. B12 and iron can be ranked from public data: inside a country, across borders, and in the external check. Vitamin A ranks well inside our countries but was the weakest externally, so not yet. Folate works only within a country. Zinc cannot be mapped: the survey itself shows no real difference between districts. One caution, in the line at the bottom: only the B12 map is its own. Rank iron districts with the vitamin A weights and you do about as well, because both maps mark the same poorer districts.

Backup, if asked: external by nutrient, African deposits, post-fix: B12 0.66 (4 of 4 tests), iron 0.44 (11 of 11), folate 0.28 (3 of 4), vitamin A 0.24 (12 of 14; both clear misses are vitamin A). Vitamin A is the strongest outcome in India and Pakistan. The table's numbers are correlations (0 = chance); measurable combinations only. Nutrient specificity (NX-01): with the country held out, B12's own weights beat every other nutrient's in 2 of 3 countries; iron 2 of 8, vitamin A 0 of 8, folate 1 of 3; the iron indices rank vitamin A districts at least as well as vitamin A's own (0.39 and 0.34 against 0.33). Zinc: nested variance components put Malawi children's zinc at district share 0.000 and cluster share 0.172, the largest cluster share in the study (explore/out/22_zinc_variance_components.csv): the survey's zinc differences sit between communities, not districts.
''')))
out.append(slide("What does this mean for policy? The vitamin B12 model", notes('''
ANDREW, about 50 seconds (118 words).

What does this mean for a programme? Take our best model, women's B12. It puts 70 to 75 of every 100 district pairs in the survey's order of average B12 level, and 65 to 71 with the country's own survey removed. Inside a surveyed country, the fifth of districts it ranks worst hold 39 to 46 per cent of deficient women: about twice a random choice. Across borders the ranking of average level holds, but the share deficient does not, so there it needs a survey to anchor it. Its district prevalences are typically about 9 points from the survey's, about as close as the survey's own regional figures. And in WHO survey data from Zambia, Ethiopia and Nigeria it pointed the right way in all four checks.

Backup, if asked: pairs on average B12 level (v3_percell_pairs.csv, post-fix): Ghana 70 in-country, 71 held out; Malawi 75 and 65 (Sierra Leone 71 held out, but its B12 prevalence, 0.6 per cent, is under the 2 per cent screen). Worst fifth, inside a surveyed country (nce_targeting_metrics.csv, post-fix): Ghana 46 per cent (regional figures 45, perfect 65), Malawi 39 (27, 65). With the country's own survey removed the worst fifth holds 12 per cent (Ghana) and 11 (Malawi) of deficient women, less than a random fifth, and its districts are less deficient than average (Ghana 4.9 against 7.8 per cent): on the share deficient B12 transports at 0.24 and 0.22, against 0.60 and 0.44 on average level. Prevalence error, each district hidden, population-weighted, model calibrated for level (benchmarks_v2_cells.csv, post-fix): Ghana 9.3 points (regional figures 9.1), Malawi 9.1 (8.9). The 90 per cent bands were rerun on 28 September (median half-width 23 points over all combinations), but script 66 still anchors on a national-estimates file that predates the survey fixes (Malawi B12 anchor 2.3 against 10.9 per cent), so no band is quoted. B12 is the one nutrient whose map is its own: its weights beat every other nutrient's in Ghana and Malawi (NX-01). External B12, post-fix: 4 of 4 in Africa (mean 0.66); Côte d'Ivoire's nine 2007 survey zones in almost the same order (0.95); Pakistan's B12 is a miss.
''')))
out.append(keep("Next steps", subs=[("Two candidate models were written down on 3 September 2026, so a fifth survey is a real test.", "Two candidate models were written down on 3 September 2026, so a fifth survey is a real test. Other surveys help a country with no survey, not a surveyed one: adding the other countries' weights to a surveyed country's own made its ranking slightly worse (BO-01: -0.02, better in 8 of 18).")]))

# ------------------------------------------------------------ conclusions (mirrors the outline)
out.append(slide("Conclusions", notes('''
ANDREW, about 60 seconds (138 words).

To answer our five questions. Can public data tell us how common a deficiency is? Not on their own: without a survey, national estimates missed by up to 38 points, and district levels need the survey's national figure. Can they rank a country's districts? Yes, modestly: two pairs in three in the survey's order, about as well as the survey's own regional figures. In a country with no survey? Mostly: six pairs in ten, and the right direction in 37 of 41 checks in six more countries. Which data? Climate and soil carry the cross-border ranking; and only B12 has a map of its own. And for programmes: with a budget per person, rank by the model; with a budget in districts, start with the most populous. It does not replace a survey. Everything is on the dashboard; the QR code is on the screen.

Backup, if asked: national levels without the country's survey: median miss 7 points, largest 38 (national_levels_sl.csv, national_vmnis_loco_sl_pred.csv). Pairs, post-fix, 14 in-country combinations: model 66, survey's regional figures 65 (scored fairly; model ahead in 9 of 14); held out 62 over 15. External (WHO VMNIS regions): Africa 12 of 12 positive at level (0.40), prevalence 18 of 21 (0.34), South Asia 7 of 8. Nutrient specificity: NX-01. Targeting (TC-01): a fifth of the people covered in the model's highest-rate districts reaches 31 per cent of the deficient inside a surveyed country (regional figures 28, random 20, perfect 47) and 25 with no survey; with a budget of a fifth of the districts, the most populous reach 49 and 45 per cent, and the model's expected cases add nothing. District levels in a new country: the model's ranking adds nothing to the national figure (LV-02). About half of the district prevalence error is the survey's own sampling noise (NZ-01).
''')))

# ------------------------------------------------------------ appendix
out.append("# Appendix\n\n")
out.append(keep("How does the model work, and how was it tested?"))
out.append(keep("How well does it rank districts inside a surveyed country?", subs=[
    ("Overall 68 of 100 pairs in the survey's order against 61 for the regional average (correlations 0.53 and 0.36).",
     "Overall 68 of 100 pairs in the survey's order against 61 for the regional average as drawn (correlations 0.53 and 0.36). The ruler scores the regional average with each district left out of its own region's figure, which reverses two surveyed districts of the same region; scored as a planner would use it (every district in a region given the same figure), it is 63 (scripts/policy_deck/39).")]))
out.append(keep("What happens when a whole region is held out?", subs=[
    ("Backup, if asked: SL-06 run, every arm on the same folds, 18 nutrient-country combinations (17 on the share deficient, where one lacks the all-layer elastic net). The benchmark table gives the index 0.39 on average status in this estimand; this run 0.38.",
     "Backup, if asked: the post-fix SuperLearner run of 28 September (NS-01, rank-normalised predictors), every method on the same folds, 18 nutrient-country combinations. Average status: index 0.39, neighbour map plus index 0.38, neighbour map 0.37, SuperLearner 0.37, elastic net on all layers 0.28. Share deficient: 0.28, 0.26, 0.26, 0.21, -0.01. The elastic net on the domain components is not in this run and is no longer shown.")]))
out.append(slide("How much ranking accuracy is lost without the country's own survey?", "![](" + F6 + "v6_incountry_vs_heldout.png){width=\"10.2in\"}\n\n" + notes('''
Each measurable combination: district pairs in the survey's order inside the country (each district hidden in turn, circle) and with the whole country held out (triangle), with the survey's regional figures scored fairly (tick). Over the 14 with an in-country test: 66 inside, the survey's regional figures 65; held out 63 over the 13 with both tests (62 over all 15); held out as good or better in 4 of 13 (Ghana children's vitamin A and iron, The Gambia women's vitamin A, Ghana B12). The large losses are Malawi children's iron and folate in Ghana and Malawi. Sierra Leone has no in-country test (14 districts); zinc, measured only in Malawi, has no held-out test. Borrowing the other countries' weights inside a surveyed country did not help (BO-01). Script: scripts/policy_deck/43.
''')))
out.append(keep("Accuracy by nutrient, against the best possible"))
out.append(keep("How do proxy models compare with the survey's own information?"))
out.append(keep("Is it just geography? Comparison with DHS-style geostatistics", subs=[
    ("(0.495 to 0.479).", "(0.495 to 0.479). Pre-fix: this comparison rests on the 8-9 September cluster targets, before the survey fixes; a post-fix rerun needs 4 to 6 hours of geostatistical fits.")]))
out.append(keep("How can it help plan the next survey?"))
out.append(keep("How does it compare with other methods, including machine learning?", subs=[
    ("mnf15_v3/fig1_model_comparison.png", "mnf15_v6/a6_model_comparison.png"),
    ("Inside a surveyed country: the Domain-PC index 0.40, neighbour smoother with proxies 0.39, neighbour smoother alone 0.39, the full SuperLearner 0.38, a twenty-layer composite 0.36, penalised regression 0.33, the regional average 0.31. A region held out: 0.39, 0.38, 0.38, 0.36, 0.34, 0.26. A country with no survey: index 0.30, SuperLearner 0.27, twenty layers 0.23, penalised regression 0.28; the smoothers and the regional average need a survey. The SuperLearner row comes from the SL-06 run with the same design, in which the index scored 0.39, 0.38 and 0.28.",
     "Post-fix. Inside a surveyed country: the Domain-PC index 0.39, neighbour smoother with proxies 0.38, neighbour smoother alone 0.37, the full SuperLearner 0.37, a twenty-layer composite 0.35, penalised regression 0.27, the regional average 0.31. A region held out: 0.39, 0.38, 0.37, 0.37, 0.34, 0.20. A country with no survey: index 0.28, SuperLearner 0.28, twenty layers 0.24, penalised regression 0.24; the smoothers and the regional average need a survey. The index, SuperLearner, smoothers and penalised regression come from the 28 September SuperLearner runs on the same folds (the penalised regression there is survey-weighted; the benchmark version scores 0.33, 0.27 and 0.26); the regional average and the twenty layers come from separate runs, in which the index scored 0.40, 0.39 and 0.30. The ensemble trails the index by 0.015 inside a country (better in 6 of 18), 0.014 with a region hidden (9 of 18) and 0.015 across borders (5 of 21) (a6_model_comparison_values.csv).")]))
out.append(keep("Which data sources contribute most to model prediction performance?"))
out.append(slide("What can we learn about the consistently most important predictor domains?", notes('''
The findings behind the variable-importance slide, across every model we fitted (post-fix).

''' + imp_backup + " Drop-a-domain test (domain_ablation_loco.csv, average status, 22 cross-country combinations): dropping climate costs 0.030 (The Gambia 0.030, Ghana 0.014, Malawi 0.030, Sierra Leone 0.045), soil 0.029 (0.043, 0.014, 0.023, 0.040); climate alone minus all layers +0.028 (better in 3 of 4 countries), soil alone +0.029 (2 of 4). The satellite embedding carries weight but adds almost nothing across borders: dropping it costs 0.002, and alone it scores 0.10 lower. Top five by weight: all 22 cross-country and all 6 pooled models, and 21 of 24 single-country models (not Malawi children's iron, Malawi women's zinc or Sierra Leone women's B12). Weight shares largely follow how many layers a family has (climate, soil and the embedding are 44 per cent of the model's 383 columns and get 47 per cent of the weight), so weight alone is not evidence of importance; the drop-a-domain test is.")))
out.append(keep("What can we learn about unused predictor domains?"))
out.append(keep("Which outcomes and cut-offs does the analysis use?"))
out.append(slide("Does targeting reach more deficient people?", "![](" + F6 + "v6_targeting_rules.png){width=\"10.4in\"}\n\n" + notes('''
Share of a country's deficient people in the districts chosen, under two budget rules (TC-01, scripts/protocol_v2/74; design written down before the run; script 12's units and folds, reproduced exactly). Grey dots: each measurable combination; large dots: the mean.

Budget of districts (a fifth of them): the most populous districts reach 49 per cent inside a surveyed country and 45 with no survey; the model's expected cases (rate times population) 49 and 44, no gain; ranking by the model's rate 25 and 18. Budget of people (a fifth of the target population): the model's highest-rate districts reach 31 per cent inside a surveyed country (the survey's regional figures 28, random 20, perfect 47) and 25 with no survey (above 20 in 9 of 15). So rank by rate when the budget is per person, by population when it is a number of districts. In a new country the anchored rate is often flat, because the training countries' own cross-border ranking was at or below zero for B12 and several iron combinations. Findings note: docs/findings/TC-01_TARGETING_BY_CASES_2026-09-28.md.
''')))
out.append(slide("Accuracy pooled over all nutrients, in district pairs", "![](" + F6 + "v6_pooled_pairs.png){width=\"10.0in\"}\n\n" + notes('''
Inside a surveyed country, each district hidden, the 14 measurable combinations (post-fix): the model puts 65.6 of 100 district pairs in the survey's order; a map of neighbouring districts built from the survey 65.9; the survey's own regional figures 65.3 (scored fairly: every district in a region gets the figure of the region's other surveyed districts); a coin 50. The model for B12 alone: 73 (Ghana 70, Malawi 75). On the 17 per cent of pairs the survey clearly separates: model 81.0, neighbour map 81.5, regional figures 80.5 (SEP-01). Inside a surveyed country geography does much of the work; the model's advantage is in a country with no survey, where the regional figures and the neighbour map do not exist. Replaces the full talk's odds ruler (pre-fix: 56, 59, 60, 70).
''')))
out.append(keep("Every combination of training countries", subs=[
    ("../../results/figures/full_talk_extracts/ft_training_curve.png){width=\"9.2in\"}", "../../results/figures/mnf15_v6/a6_training_curve.png){width=\"6.8in\"}"),
    ("the gain is in The Gambia and Ghana.", "the gain is in The Gambia and Ghana. Post-fix (training_curve_climate_soil.csv, 28 Sep): average status, climate and soil 0.31, 0.36 and 0.40 with one, two and three training countries, the full set 0.19, 0.25 and 0.30; prevalence 0.24, 0.27, 0.31 and 0.08, 0.14, 0.23. The Gambia 0.38 to 0.59 and Ghana 0.30 to 0.39; Malawi and Sierra Leone flat. Over all 22 combinations, 7 of them not measurable; on the measurable ones the full set gives 0.23, 0.31 and 0.40.")]))
out.append(keep("Does the model track the survey, district by district?", subs=[
    ("../../results/figures/full_talk_extracts/ft_pred_vs_observed.png){width=\"9.2in\"}", "../../results/figures/mnf15_v6/a6_pred_vs_observed.png){width=\"10.0in\"}"),
    ("the model ranks, and compresses the range.", "the model ranks, but its spread is as wide as the survey's or wider, so levels are off."),
    ("much of the vertical scatter is the survey's noise.", "much of the vertical scatter is the survey's noise. Post-fix (viz/oof_child_iron.csv, 28 Sep), children's iron: Ghana 0.50 and 16.0 points average error, Malawi 0.25 and 14.4, The Gambia 0.37 and 12.9; 68, 59 and 63 of 100 district pairs in order.")]))
for t in ["The 2007 Côte d'Ivoire survey and the B12 prediction", "The live dashboard"]:
    out.append(keep(t))

txt = re.sub(r"\n{3,}", "\n\n", "".join(out))
# figures redrawn on the post-fix measurability screen (28 Sep evening) live in mnf15_v6: use them wherever they exist
import glob as _glob
for _f in sorted(os.path.basename(x) for x in _glob.glob(ROOT + "results/figures/mnf15_v6/*.png")):
    txt = txt.replace("mnf15_v5/" + _f, "mnf15_v6/" + _f)
open(V6, "w", encoding="utf-8").write(txt)
print("wrote", V6, "with", len(re.findall(r"^## ", txt, re.M)), "slides")
