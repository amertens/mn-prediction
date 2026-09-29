"""Build docs/slides/MNF15-talk-2026-09-v5.qmd from the v4 qmd, following Andrew's
edits and comments in docs/slides/MNF15-talk-2026-09-v4-ANM.pptx (28 September; read
only) and the post-fix results of the RR-14 rebuild (27-28 September).

    python scripts/policy_deck/36_build_v5_source.py

The v5 deck is made in three steps:
    python scripts/policy_deck/36_build_v5_source.py            # this file: qmd
    bash scripts/render_deck.sh docs/slides/MNF15-talk-2026-09-v5.qmd --forum --no-label --text 18,12 \\
         --concept docs/slides/MNF15-talk-2026-09-v5.concept.yaml
    python scripts/policy_deck/38_merge_v5_slides.py            # slides copied from v4-ANM

Only slides that pandoc or the concept builder draws live in the qmd; the full-talk
slides and Andrew's hand-edited slides are copied from v4-ANM by script 38, next to
the slide named in its plan. Notes of unchanged slides are taken from the v4 qmd.
"""
import re

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
V4 = ROOT + "docs/slides/MNF15-talk-2026-09-v4.qmd"
V5 = ROOT + "docs/slides/MNF15-talk-2026-09-v5.qmd"
F5 = "../../results/figures/mnf15_v5/"
F4 = "../../results/figures/mnf15_v4/"

src = open(V4, encoding="utf-8").read()
blocks = {}
for m in re.finditer(r"^## ([^\n]+)\n(.*?)(?=^## |^# |\Z)", src, re.S | re.M):
    blocks[m.group(1).strip()] = m.group(2)


def notes_of(title):
    b = blocks[title]
    return b[b.index("::: {.notes}"):].rstrip() + "\n"


def image_of(title):
    m = re.search(r"^!\[\]\((.+?)\)\{(.+?)\}", blocks[title], re.M)
    return m.group(1), m.group(2)


def slide(title, notes, image=None, width="9.6in"):
    img = f"![]({image}){{width=\"{width}\"}}\n\n" if image else ""
    return f"## {title}\n\n{img}{notes.strip()}\n\n"


def keep(title, new_title=None, image=None, width=None):
    """Reuse a v4 slide: same notes, optionally a new title and image."""
    img = None
    if "![](" in blocks[title]:
        path, attr = image_of(title)
        img = image or path
        width = width or re.search(r'width="([^"]+)"', attr).group(1)
    elif image:
        img = image
    return slide(new_title or title, notes_of(title), img, width or "9.6in")


def notes(text):
    return "::: {.notes}\n" + text.strip() + "\n:::\n"


header = src[:src.index("```{r setup}")]
setup = '''```{r setup}
# MNF 2026 joint talk, VERSION 5 (28 September 2026). v1 to v4 are untouched.
# v5 = Andrew's edits, deletions and comments in docs/slides/MNF15-talk-2026-09-v4-ANM.pptx
# (his working copy; read, never written), on the results rebuilt with the 15 September
# survey-report fixes (RR-14, scripts/protocol_v2/rerun_survey_fixes_2026-09-27.sh):
#   - the order and deletions of v4-ANM; slides he marked "appendix" moved there, in order
#   - no "full talk" tags; slides copied from v4-ANM keep his hand edits (script 38)
#   - figures redrawn to his comments in results/figures/mnf15_v5/ (scripts 26 FIG_VER=v5,
#     29 / 31 / 34 with FIG_OUT=results/figures/mnf15_v5, 37): surveys "n of N", the test
#     schematic plain, national prevalence for every outcome with a national panel,
#     prevalence error with an X for the national estimate, Ghana held out without text,
#     per-nutrient pairs in one panel with prevalence in brackets, no Cote d'Ivoire in the
#     check map, no text baked under the appendix charts
#   - numbers in the main slides are post-fix; appendix figures from the full talk, fig1,
#     fig9 and the full-talk extracts are from the 8 September build
# Render: python scripts/policy_deck/36_build_v5_source.py; bash scripts/render_deck.sh
#   docs/slides/MNF15-talk-2026-09-v5.qmd --forum --no-label --text 18,12
#   --concept docs/slides/MNF15-talk-2026-09-v5.concept.yaml; python scripts/policy_deck/38_merge_v5_slides.py
suppressPackageStartupMessages({library(dplyr); library(knitr)})
```

'''

out = [header, setup]

# ---------------------------------------------------------------- main body
out.append(keep("Introduction and study objectives"))
out.append(slide("Talk outline", notes('''
SONJA, about 25 seconds (60 words).

Here is what we will tell you. First, aggregate public data carry real information about micronutrient deficiency; the top predictors are climate, soil and land cover. Second, a simple model does best, because there are few places to learn from. Third, it ranks districts no survey reached. Fourth, it works for some nutrients more than others. And fifth, it extends surveys; it does not replace them.
''')))
out.append(keep("Conceptual framework of iron deficiency"))
out.append(keep("Which of these causes can public data measure?"))
out.append(slide("Which surveys did we use to train the models?", notes('''
SONJA, about 25 seconds (58 words).

The models learn from four national micronutrient surveys: The Gambia 2018, Ghana 2017, Sierra Leone 2013 and Malawi 2015 to 16. Together, 206 districts where blood was drawn: 30 of The Gambia's 37 districts, 75 of Ghana's 260, all 14 of Sierra Leone's, and 87 of Malawi's 243 areas. These data exist because these teams collected them. Thank you.

Backup, if asked: outcome definitions, adjustments and cut-offs are in the appendix. Sierra Leone's child file holds only anaemic children (532 of 654 assayed). Malawi's areas are Traditional Authorities. Ghana's survey reached 90 clusters in its 75 districts. Survey partners for Ghana's 2017 survey: University of Ghana, GroundWork, University of Wisconsin-Madison, KEMRI-Wellcome Trust, with UNICEF and Global Affairs Canada.
'''), F5 + "v4_surveys.png"))
out.append(slide("How do we get honest estimates of model performance?", notes('''
ANDREW, about 40 seconds (90 words).

Thank you, Sonja. Every number I show comes from data the model never saw. On the left, cross-validation: we split the surveyed districts into five parts, fit on four, predict the fifth, and rotate until every district has been predicted without its own data; then we repeat with ten random splits. Next, the same with a whole region held out. And a whole country: Ghana predicted from the other three countries. Every score is on the orange districts, the held-out test data.

Backup, if asked: folds are cut at the district or above, never through a district's respondents; component weights and every other choice are learned inside the training folds. District folds are random, not grouped in space; the region hold-out is the spatially grouped test. The middle map is the benchmark's first fold draw (make_folds_v2, rep 1).
'''), F5 + "v4_cv_schematic.png"))
out.append(slide("Why not a more complex machine-learning model?", notes('''
ANDREW, about 45 seconds (100 words).

The obvious question for the machine-learning people here: why such a simple model? This is a big-data project in predictors, 575 layers, but a small-data project in places: 14 to 87 districts per country, each measured in about one community. We compared 14 methods on the same held-out districts, including random forests, boosting and a SuperLearner ensemble that combines twelve. None did more than slightly better on average across nutrients, and none was better everywhere. So we use the simplest one, a weighted sum of principal components, with nothing tuned.

Backup, if asked: average status, SL-06 run (8 September targets; the survey fixes moved the index by less than 0.005). One district in five hidden: index 0.39, SuperLearner 0.38, random forest 0.38, boosting 0.37, elastic net 0.32, lasso 0.29. A region hidden: 0.38, 0.36, 0.39, 0.40, 0.23, 0.21. A country hidden: 0.28, 0.27, 0.25, 0.20, 0.29, 0.31. Per outcome the winner changes; folate is the index's worst case. On the share deficient the index leads inside a country (0.27 against at most 0.24) and with a region hidden (0.27 against at most 0.23); across borders the ensemble is marginally ahead (0.20 against 0.19). A sample-size finding, not a claim about machine learning in general.
''')))
out.append(slide("But can proxy models recover national-level prevalence estimates?", notes('''
ANDREW, about 45 seconds (105 words).

A question from our January meeting: can proxy models give a national prevalence? Inside a surveyed country the survey already gives it. Without the country's survey, not reliably. On the left, every outcome we can test, each of our countries predicted from national indicators with its own survey left out: close for some, but 28 points too high for Sierra Leone's children's vitamin A, 38 too low for its women's folate, 34 too low for Malawi's zinc. On the right, WHO's national survey data from 21 to 69 countries: the model beats the average of the other countries only for children's vitamin A and folate. Levels do not travel between surveys; district order does.

Backup, if asked: national track: WHO VMNIS national survey panel with World Bank indicators; a SuperLearner over mean, ridge, lasso, elastic net and random forest, folds grouped by country (national_vmnis_loco_sl.csv, _pred.csv, national_levels_sl.csv; scripts/covariates/19b; figure script 37). Left panel: vitamin A at each survey's year against the survey's own national figure; folate, B12 and zinc against the country's own VMNIS row. Median miss 7 points, largest 38. Iron has no national panel. The January version compared national predictions made inside each country, which match by construction. Why levels do not travel: assays, cut-offs, inflammation adjustment and RBP-to-retinol calibration differ between surveys (ferritin levels differ several-fold).
'''), F5 + "v5_national_all.png"))
out.append(slide("Can it rank Ghana's districts without Ghana's survey?", notes('''
ANDREW, about 50 seconds (112 words).

Now the test that matters for most of this region: we removed Ghana's survey entirely. The model learned only from The Gambia, Sierra Leone and Malawi. On the left, the survey's ranking of districts for children's iron; next to it, the model's, with no Ghanaian data at all. It puts 69 of every 100 pairs of Ghana's districts in the survey's order. This is one of our better cases. On the right, each country held out in turn, all nutrients: 70 for The Gambia, 62 for Ghana, 57 for Malawi, 56 for Sierra Leone. A shortlist, not a measurement.

Backup, if asked: correlation 0.54 held out against 0.53 in-country for this outcome, within each other's intervals. Why about as good as the model that saw Ghana's survey: it learned from 131 districts in three countries instead of about 60 in Ghana; over all six of Ghana's outcomes it does a little worse, 62 pairs in 100 against 65. The survey's own regional averages, which need the survey, put 61 in 100 in order. Ghana held out does better than Ghana in-country for children's iron, children's vitamin A and B12, worse for women's vitamin A, iron and folate.
'''), F5 + "v2_ghana_heldout.png", "9.4in"))
out.append(slide("How well does it rank districts, nutrient by nutrient?", notes('''
ANDREW, about 50 seconds (110 words).

Here is every nutrient and country, with a 95 per cent interval, and the national prevalence in brackets. Circles: inside a surveyed country. Triangles: the whole country held out. B12 is the most mappable: about 70 of every 100 pairs right, at home and across borders, probably because it follows animal-source foods. Vitamin A works well in The Gambia, less in Ghana and Malawi, and was the weakest in the external check. Iron works in most places, not all: across borders, Malawi's children are no better than a coin. Folate works inside a country but not across borders, most likely because the surveys measured it differently.

Backup, if asked: measurable combinations only; Sierra Leone is not shown (14 districts leave too few to train on inside the country). In-country intervals clear 50 per cent in 12 of 14; held out in 10 of 16 (all measurable). Averages over the 14 shown: 66 per cent inside a country, 63 held out, regional average 61. Orange diamonds: the survey's regional average, computed without the district, which is what a survey that missed the district would know. Ghana's women's vitamin A dropped out after the fixes (prevalence now 1.6 per cent, under the 2 per cent screen); Malawi's children's vitamin A came in (8.7 per cent). Folate: immunoassay in Ghana and Sierra Leone, microbiologic assay in Malawi.
'''), F5 + "v4_percell_pairs_combined.png", "9.6in"))
out.append(keep("What does the model rely on?", "What variables drive model predictions", F5 + "v4_domains_child_iron.png"))
out.append(slide("What can we learn about the consistently most important predictor domains?", notes('''
ANDREW, about 50 seconds (112 words).

Across every model we fitted, the same picture. Climate, satellite land cover and soil are in the top five domains of all 22 cross-country models and of every single-country model; with vegetation greenness they carry about half of the weight, for every nutrient. On its own, climate or soil ranks a country held out about as well as all the layers together. Then a few markers recur, always in the same direction: grassland, ruminants and cereal cropping mark more deficiency; plant productivity, market access and child overweight mark less. And each nutrient adds its own: meat and fish in children's diets for B12, anaemia for iron and folate, child wasting for vitamin A. Markers of where, not causes.

Backup, if asked: post-fix importance and ablation tables (index_importance_domains.csv, index_importance_columns.csv, domain_ablation_loco_summary.csv). Pooled weight shares, mean over six outcomes: climate 18 per cent, satellite embedding 16, soil 13, greenness 6; climate + soil + satellite + greenness 45 to 59 per cent by outcome. Country held out: climate alone 0.33, soil alone 0.33, all layers 0.30; removing climate or soil costs 0.03, anaemia 0.02, crops and livestock about 0.01 each. Layers in the top 20 of four outcomes, same sign every time: malaria mortality, child overweight, relative staple price, plant productivity (less deficiency), grassland cover (more). Malaria marks LESS deficiency for iron and vitamin A, most likely because inflammation raises ferritin and lowers RBP differently from deficiency: a feature of the biomarkers. Schooling, water and sanitation, household diets and infant feeding help inside a country but make cross-border ranking slightly worse. Fertility and reproductive health are in no top list.
''')))
out.append(slide("What about countries with no survey?", notes('''
ANDREW, about 55 seconds (120 words).

Most countries in this region have no recent biomarker survey. Côte d'Ivoire is one. Trained on our four countries, the climate-and-soil version of the model ranks its 33 districts for children's iron; the northern savanna comes out worst. Can you trust a map like this? WHO's micronutrients database holds regional results from national surveys in Zambia, Ethiopia, Sudan, Nigeria, Pakistan and India. None were used to build the model. In Africa its ranking of their regions pointed the right way in all 12 tests of average status and 18 of 21 tests of prevalence; in South Asia, 7 of 8. And Côte d'Ivoire's own 2007 survey ranks its nine zones for B12 almost exactly as the model does.

Backup, if asked: XV-01/02, rerun on the fixed targets: climate-and-soil index at admin-1; Africa average status 0.40 (country-block null 0.20), prevalence 0.34; off the continent 0.38. Iron points the right way in all 14 iron tests across the six countries (weakest Ethiopia, 0.01). Regions, not districts; not the pre-registered test, which needs a new country's own microdata at district level. Côte d'Ivoire 2007 B12: Spearman 0.95 over nine zones. The Côte d'Ivoire map is from the 8 September build until the downstream run finishes.
'''), F5 + "v2_civ_checked.png", "9.4in"))
out.append(slide("Which deficiencies can it map?", notes('''
ANDREW, about 40 seconds (88 words).

Putting it together. B12 and iron can be mapped from public data: inside a country, across borders, and in the external check. Vitamin A ranks well inside our countries but was the weakest externally, so not yet. Folate works only within a country. Zinc, measured only in Malawi, cannot be mapped: the model is no better than chance, because the survey itself shows no real difference between districts. How common a deficiency is does not decide this: folate and zinc are the commonest in our data, and the hardest.

Backup, if asked: external by nutrient, African deposits, post-fix: B12 0.66 (4 of 4 tests), iron 0.44 (11 of 11), folate 0.28 (3 of 4), vitamin A 0.24 (12 of 14; both clear misses are vitamin A). Vitamin A is the strongest outcome in India and Pakistan. The table's numbers are correlations (0 = chance); measurable combinations only.
'''), F5 + "v2_nutrients.png", "9.4in"))
out.append(slide("What does this mean for policy? The vitamin B12 model", notes('''
ANDREW, about 45 seconds (98 words).

What does this mean for a programme? Take our best model, women's B12. It puts 70 to 75 of every 100 district pairs in the survey's order, and 65 to 71 with the country's own survey removed. The fifth of districts it ranks worst hold 39 to 46 per cent of deficient women: about twice a random choice. Its district prevalences are typically within about 6 points of the survey's. And in WHO survey data from Zambia, Ethiopia and Nigeria it pointed the right way in all four checks, as it did for Côte d'Ivoire's 2007 survey zones.

Backup, if asked: pairs (v3_percell_pairs.csv, post-fix): Ghana 70 in-country, 71 held out; Malawi 75 and 65; Sierra Leone 71 held out, 14 districts. Worst fifth (nce_targeting_metrics.csv, post-fix): Ghana 46 per cent (regional averages 45, perfect 65), Malawi 39 (regional 27, perfect 65). Prevalence error, each district hidden (conformal_prev_cells.csv, 8 September build until rerun): Ghana 6.6 points, Malawi 5.5; 90 per cent bands 21 and 27 points either side. External B12, post-fix: 4 of 4 in Africa (mean 0.66); Pakistan's B12 is a miss.
''')))
out.append(keep("Conclusions"))
lrn = notes_of("What have we learned, and what comes next?").replace("inside a surveyed country, about two pairs in three in the survey's order", "inside a surveyed country, about two pairs in three in the survey's order").replace(
    "pairs 67 against 62 per cent inside a country (14 measurable combinations), 62 held out (22)",
    "pairs 66 against 61 per cent inside a country (14 measurable combinations), 62 held out (22), post-fix")
out.append(slide("What have we learned, and what comes next?", lrn))
out.append(keep("Next steps"))

# ---------------------------------------------------------------- appendix
out.append("# Appendix\n\n")
out.append(keep("How does the model work, and how was it tested?"))
out.append(slide("How well does it rank districts inside a surveyed country?", notes('''
Ghana, children's iron, each district's own blood results hidden in turn. Left: what the survey measured in the 75 districts it reached; middle: the model's held-out prediction; right: the model for all 260 districts. Central Gonja: 13th worst in the survey, 51st by its regional average, 21st by the model. Prestea-Huni Valley: 73rd in the survey, 39th by its regional average, 69th by the model. Overall 68 of 100 pairs in the survey's order against 61 for the regional average (correlations 0.53 and 0.36).

Honest caveats: on the share deficient rather than average status, the model ties the regional average here (0.50 each), and a map of neighbouring districts built from the survey does as well as the model (69 of 100 pairs); inside a surveyed country, geography does much of the work. The model is wrong in places: its median district is 13 places from the survey's rank, the regional average's 18; the largest miss is Atiwa East (24th in the survey, 65th in the model). The regional average is the district's own region averaged over the other surveyed districts there (the district itself left out).
'''), F5 + "v2_ghana_infill.png", "9.4in"))
out.append(keep("What happens when a whole region is held out?", image=F5 + "v4_region_heldout.png"))
acc = notes_of("Accuracy by nutrient, against the best possible").replace("Inside a country we reach 63 to 73.", "Inside a country we reach 63 to 75.")
out.append(slide("Accuracy by nutrient, against the best possible", acc, F5 + "v4_accuracy_by_nutrient.png", "9.4in"))
skill = notes_of("How do proxy models compare with the survey's own information?").replace(
    "Public data alone improve on it by about 11 per cent on average: as much as a map of neighbouring districts built from the survey, and more than the survey's regional averages, at 5 per cent.",
    "Public data alone improve on it by about 11 per cent on average: as much as a map of neighbouring districts built from the survey (10 per cent), and more than the survey's regional averages, at 3 per cent.").replace(
    "public data 11.4 per cent, neighbour map 11.4, regional averages 4.6", "public data 10.9 per cent, neighbour map 10.2, regional averages 3.1 (post-fix)")
out.append(slide("How do proxy models compare with the survey's own information?", skill, F5 + "v4_skill_by_outcome.png"))
out.append(keep("Is it just geography? Comparison with DHS-style geostatistics"))
out.append(keep("How can it help plan the next survey?"))
out.append(keep("How does it compare with other methods, including machine learning?"))
out.append(keep("Which data sources contribute most to model prediction performance?", image=F5 + "v4_sources.png"))
out.append(keep("What can we learn about missing predictor domains?", "What can we learn about unused predictor domains?"))
cut = notes_of("Which outcomes and cut-offs does the analysis use?")
cut = re.sub(r"IMPORTANT caveat for Andrew:.*?propagated to the result tables\.",
             "The main-body results in this deck use the targets rebuilt on 27 September with these corrections (vitamin A on the RBP < 0.70 rule, non-pregnant women and children 6-59 months only, Sierra Leone women's iron on Thurnham, Malawi repeated Traditional Authority names, Malawi zinc on the survey's own flag); appendix figures copied from the full talk are from the 8 September build.", cut, flags=re.S)
tbl = blocks["Which outcomes and cut-offs does the analysis use?"]
table = tbl[:tbl.index("::: {.notes}")].strip()
for a, b in [   # the table as the rebuilt (post-fix) targets define the outcomes
        ("Children under 5; women 15-49", "Children 6-59 months; non-pregnant women 15-49"),
        ("| Women 15-49 |", "| Non-pregnant women 15-49 |"),
        ("Retinol-binding protein, BRINDA-adjusted, converted to retinol with each survey's calibration", "Retinol-binding protein, BRINDA-adjusted"),
        ("| Under 0.70 µmol/L |", "| RBP under 0.70 µmol/L |")]:
    assert a in table, a
    table = table.replace(a, b)
out.append("## Which outcomes and cut-offs does the analysis use?\n\n" + table + "\n\n" + cut + "\n")
tgt = notes_of("Does targeting reach more deficient people?").replace("24.8 per cent of the deficient people; the survey's regional averages 20.5; a random fifth 20 by definition; a perfect map 41.0",
                                                                    "25.1 per cent of the deficient people; the survey's regional averages 21.0; a random fifth 20 by definition; a perfect map 40.3 (post-fix)")
out.append(slide("Does targeting reach more deficient people?", tgt, F5 + "v4_targeting.png", "9.0in"))
out.append(keep("Accuracy pooled over all nutrients, in district pairs"))
out.append(keep("Every combination of training countries"))
civ = notes_of("The 2007 Côte d'Ivoire survey and the B12 prediction").replace("predicted = expit(logit(0.181) + 0.33 x 1.48 x z)", "predicted = expit(logit(0.181) + 0.32 x 1.48 x z)").replace(
    "0.33 how well that ranking held up", "0.32 how well that ranking held up").replace("9 to 31 per cent predicted against 0 to 49 measured, mean error 8.7 points", "9 to 31 per cent predicted against 0 to 49 measured, mean error 9.0 points")
out.append(slide("The 2007 Côte d'Ivoire survey and the B12 prediction", civ, F5 + "v4_civ_b12.png", "9.4in"))
out.append(keep("Does the model track the survey, district by district?"))
out.append(keep("The live dashboard"))

txt = "".join(out)
txt = re.sub(r"\n{3,}", "\n\n", txt)
open(V5, "w", encoding="utf-8").write(txt)
print("wrote", V5, "with", len(re.findall(r"^## ", txt, re.M)), "slides")
