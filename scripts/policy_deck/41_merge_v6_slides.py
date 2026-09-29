"""Copy into the rendered v6 MNF15 deck the slides that are not built from source (the
full-talk slides and Andrew's hand-edited slides from docs/slides/MNF15-talk-2026-09-v4-ANM.pptx,
read only), in the v6 order (see script 40). Uses the copy helpers of script 38.

    python scripts/policy_deck/41_merge_v6_slides.py

Against v5: "Both a big-data opportunity...", "How do we choose among algorithms?" and "How
much of what a survey can see..." go to the appendix; v4-ANM's old closing slide (59) is
not copied; the survey's regional average is quoted as scored fairly (script 39); the
worst-fifth line compares with a random fifth and a perfect map only (the regional-average
comparison there uses the older scoring, which under-rates it).
"""
import importlib.util
import os
import sys

from pptx import Presentation

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
spec = importlib.util.spec_from_file_location("m38", os.path.join(ROOT, "scripts/policy_deck/38_merge_v5_slides.py"))
m38 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m38)

ANM = m38.ANM
DST = os.path.join(ROOT, "docs/slides/MNF15-talk-2026-09-v6.pptx")
FIG = m38.FIG   # results/figures/mnf15_v5

PLAN = [
    (6, "Which of these causes can public data measure?"),
    (15, "But can proxy models recover national-level prevalence estimates?"),
    (22, "How well does it rank districts, nutrient by nutrient?"),
    # appendix
    (9, "How does the model work, and how was it tested?"), (12, None), (13, None),
    (18, "How well does it rank districts inside a surveyed country?"),
    (21, "What happens when a whole region is held out?"),
    (24, "Accuracy by nutrient, against the best possible"),
    (26, "How do proxy models compare with the survey's own information?"),
    (37, "Is it just geography? Comparison with DHS-style geostatistics"), (38, None),
    (42, "Which data sources contribute most to model prediction performance?"),
    (44, "What can we learn about unused predictor domains?"), (45, None), (46, None), (47, None),
    (49, "Which outcomes and cut-offs does the analysis use?"), (50, None),
    (52, "Does targeting reach more deficient people?"),
    (54, "Accuracy pooled over all nutrients, in district pairs"),
]
FIG6 = os.path.join(ROOT, "results/figures/mnf15_v6")
PICTURE = {k: (os.path.join(FIG6, v) if os.path.exists(os.path.join(FIG6, v)) else os.path.join(FIG, v)) for k, v in m38.PICTURE.items()}
PICTURE[38] = os.path.join(FIG6, "v6_anchor_designs.png")   # post-fix AR-01 with the flat national figure (LV-02)
TITLE_OF = {6: "575 proxy variables from 29 conceptual domains", 15: "How far off is the predicted prevalence, outcome by outcome?", 22: "What does a rank correlation of 0.40 mean on the ground?", 24: "How much of what a survey can see does the model recover?", 38: "Could a small national survey plus the model replace a district survey?", 52: "Does the model put districts in the right severity band?", 54: "What the models can and cannot do yet"}
TEXT = {k: dict(v) for k, v in m38.TEXT.items()}
TEXT[22].update({
    "In the fifth of districts it flags, deficiency runs at 32% against 24% nationally": "With a budget per person, rank by the model",
    "A map that knew every district would find 41% in its worst fifth.":
        "Covering a fifth of the people in its highest-rate districts reaches 31% of the deficient: random 20%, the survey's regional figures 28%, perfect knowledge 47%. With no survey, 25%.",
    "That fifth holds 22% of the deficient people": "With a budget in districts, start with the most populous",
    "Counted over every pair of districts. A random ordering gets 50% on average, the survey's own regional averages 56%; in a country the model has never seen, 55%.":
        "A coin toss gets 50%. The survey's own regional figures get 65%, but only where there is a survey; in a country with no survey the model still gets 62%.",
    "The survey's own regional averages capture 20%; a perfect map 48%. Real, modest, and better than today's number.":
        "The most populous fifth of districts holds 49% of the deficient people; ranking by the model's expected cases adds nothing to that.",
})
TEXT[24].update({
    "The survey itself is noisy at district level: a perfect predictor would score only 0.48": "The survey itself is noisy at district level: a perfect predictor would score only 0.47",
    "The model reaches 59% of that ceiling": "The model reaches about two-thirds of that ceiling",
    "0.28 on the prevalence target against a ceiling of 0.48. The gap is mostly the survey's noise, not the model's ignorance.": "0.30 on the prevalence target against a ceiling of 0.47. The gap is mostly the survey's noise, not the model's ignorance.",
    "And beats the best prediction you can make from survey data alone": "About as good as the survey's own regional figures",
    "The region's average, computed without the district itself: 0.31; the model 0.40, better in 13 of 18 outcome-country combinations.":
        "Scored as a planner would use them, the regional figures put 65 of 100 district pairs in order; the model 66, ahead in 9 of 14. Its value is where there is no survey.",
})
NOTES = dict(m38.NOTES)
NOTES[22] = ("ANDREW, about 40 seconds (100 words). What does that accuracy mean on the ground? Take any two districts: the model orders "
             "them as the survey does 66 times in 100, the average of the previous slide. The survey's own regional figures do about as well, "
             "65, but only where there is a survey; in a country with none, the model still gets 62. Inside a surveyed country, the fifth of "
             "districts it flags run at 31 per cent deficient against 22 nationally, and hold a quarter of the deficient people. In a country "
             "with no survey the flagged districts are still more deficient than average, but they are small, so they hold fewer of the "
             "deficient people than a random fifth.\n\n"
             "Backup, if asked: post-fix, each district held out (v3_percell_pairs.csv, v6_fair_regional_pairs.csv, nce_targeting_metrics.csv). "
             "Pairs inside a country: model 66 (14 combinations with an in-country test), the survey's regional figures 65 scored as a planner "
             "would use them (every district in a region given the same figure; model ahead in 10 of 14); held out 62 (16). The older scoring "
             "left each district out of its own region's figure and gave the regional figures 61, because it reverses two surveyed districts "
             "of the same region by construction. Worst fifth: deficiency 30.9 against 21.7 per cent nationally; a perfect map finds 47.1; "
             "share of deficient people reached 25.1 (a random fifth 20, a perfect map 40.3). Whole country held out (22 combinations): the flagged "
             "fifth runs at 22.0 per cent deficient against 19.1 nationally (higher in 13 of 22) but holds 16.8 per cent of the deficient people, "
             "under 20 in 14 of 22: the model ranks by rate, and the high-rate districts are small and rural. Ranking by predicted cases "
             "(rate times population) is the fix to test. Retitled on Andrew's comment ('why this "
             "number'): the earlier 0.40 was the pooled rank correlation.")
NOTES[6] = ("SONJA, about 25 seconds (56 words).\n\n"
            "From that exercise we assembled 575 public variables from 28 sources, in 29 conceptual domains, grouped here into eight. "
            "Every one is a district average, and all of it is free. Because there are far more variables than districts, each domain "
            "enters the model as a few summary scores, not as hundreds of raw columns.\n\n"
            "Backup, if asked: 575 columns against 14 to 87 districts per country cannot be fitted column by column, and the columns are "
            "heavily correlated within a source. Grouping them by what they measure and taking principal components inside each group keeps "
            "most of each group's variance, treats a source with 64 columns and a source with 3 as one construct each, and keeps the outcome "
            "out of the representation: the rotations are learned from the training districts' predictors alone, so nothing is selected on "
            "the biomarker. The grouping is written down: each source block brings a label, and 10 override rules "
            "(metadata/covariates/domain_overrides.csv) regroup columns by construct (infection and inflammation, immunisation, anaemia, "
            "child growth, infant feeding, supplementation, water and sanitation, education, child mortality, assets), so that Malaria Atlas, "
            "IHME, DHS and MICS versions of the same thing are weighed together. The headline model uses 383 of the 575 (no DHS microdata).")
NOTES[22] = NOTES[22].replace("ANDREW, about 40 seconds (100 words).", "ANDREW, about 50 seconds (121 words).")
assert "about 50 seconds" in NOTES[22]
NOTES_APPEND = {24: m38.NOTES_APPEND[24].replace(
    "Post-fix: model 0.40 against the regional average 0.31, better in 13 of 18.",
    "Post-fix: model 0.40 against the regional average 0.31 (Spearman, older scoring), better in 13 of 18. Scored as a planner would "
    "use them (every district in a region given the same figure), the regional figures put 65 of 100 pairs in order against the "
    "model's 66, model ahead in 9 of 14 (scripts/policy_deck/39_fair_regional_pairs.py).")}
assert NOTES_APPEND[24] != m38.NOTES_APPEND[24]
TEXT[38] = {"Four ways to get district numbers, costed against the size of survey each needs. A small national sample, 5% of a full survey, with the model's ranking spread across districts gives a median district error of 9.3 points. A survey designed to measure each district directly needs about 25% of its full size to do as well.":
            "Five ways to get district numbers, by the size of survey each needs. A small national sample, 5% of a full survey, gives a median district error of about 10 points, with or without the model's ranking spread across districts. A survey designed to measure each district directly needs about a quarter of its full size to do as well."}
TEXT[52] = {"Model: 64% of districts": "Model: 50% of districts", "Regional averages: 66%": "Regional figures: 53%",
            "Model: 91% of districts": "Model: 88% of districts", "Regional averages: 94%": "Regional figures: 90%",
            "Model: 9% of districts": "Model: 12% of districts", "Regional averages: 6%": "Regional figures: 10%"}
NOTES[22] = ("ANDREW, about 50 seconds (118 words). What does that accuracy mean on the ground? Take any two districts: the model orders "
    "them as the survey does 66 times in 100. The survey's own regional figures do about as well, 65, but only where there is a survey; in a "
    "country with none, the model still gets 62. For a programme, it depends on the budget. If it is per person, say supplements, rank by the "
    "model: covering a fifth of the people in its highest-rate districts reaches 31 per cent of the deficient, against 20 at random. If the "
    "budget is a number of districts, start with the most populous: that fifth holds almost half the deficient people, and the model adds nothing.\n\n"
    "Backup, if asked: pairs, post-fix, 14 in-country combinations (v3_percell_pairs.csv, v6_fair_regional_pairs.csv): model 66, the survey's "
    "regional figures 65 scored as a planner would use them (model ahead in 9 of 14); held out 62 (15). Budget rules (TC-01, "
    "tc01_targeting_by_cases_summary.csv; appendix figure): a fifth of the people in the model's highest-rate districts reaches 31 per cent of "
    "the deficient inside a surveyed country (regional figures 28, random 20, perfect 47) and 25 with no survey; a fifth of the districts: most "
    "populous 49 and 45 per cent, model expected cases 49 and 44, model rate 25 and 18, perfect 63 and 60. The earlier line (the worst fifth "
    "by rate holds 25 per cent) answered a question no programme asks. The next slide shows the budget rows as a figure; keep one of the two.")
TEXT[54] = {"Yes: 9.3 points at a 5% sample.": "About 10 points at a 5% sample; no closer than the national figure alone."}
NOTES[15] = m38.NOTES[15] + ("\n\nAdded 28 September: about half of these errors is the survey's own sampling noise (NZ-01: 47 per cent of the squared error over "
    "all 27 combinations, 60 over the 22 main ones; typical error against the true prevalence about 14 points rather than 19, root mean square). "
    "In a country with no survey, district levels anchored to the national figure are no closer than the national figure itself (LV-02: tie in 11 of 22 "
    "combinations at every sample size; external 11 of 29).")
NOTES[38] = ("Post-fix AR-01 (anchor_and_rank_summary.csv, climate-and-soil ranking, average status), with the design it lacked: every district given the "
    "national figure (LV-02, lv02_internal_summary.csv). Median district error at a 5% sample: national figure plus the model's ranking 9.6 points, national "
    "figure alone 10.1; per combination the two tie (11 of 22) at every sample size, so the ranking adds nothing to district levels. A survey measuring every "
    "district directly beats them from about a quarter of full size (9.2 at 25%) and is exact at full size by construction (it is scored against itself). "
    "Regional figures plus the ranking equal a regional survey. Use the model to order districts; take levels from the survey.")
NOTES[52] = ("Post-fix, WHO vitamin A bands (under 2, 2 to 10, 10 to 20, 20 per cent or more), each district held out, mean over the 8 vitamin A combinations "
    "(risk_category_accuracy.csv, scheme who_vitA): exact band, model 50 per cent, the survey's regional figures 53, a neighbour map 54; within one band 88, 90, "
    "92; two or more bands wrong 12, 10, 8. The survey's regional figures and the neighbour map do slightly better, because a band is mostly a level and they carry "
    "the country's own level. Down from 64 and 66 before the fixes: the vitamin A rule changed (RBP under 0.70), which moved many districts across band edges.")
NOTES_APPEND[54] = ("\n\nPost-fix (28 Sep): rows 1 and 2 hold (0.40 and 0.30). The anchor row is now about 10 points and no better than the national figure "
    "alone (LV-02). The last row, recomputed post-fix (block C of script 10, rerun by scripts/policy_deck/46): coverage 39 per cent (15 to 71 over 22 combinations); calibrated half-width 55 per cent of the list; worst-third calls right 52 per cent, chance 33.")


# run-level substring edits (keep bold labels such as "Who:") and notes substring edits, by v4-ANM slide number
RUN_SUB = {49: [(" children under five; non-pregnant women", " children 6-59 months; non-pregnant women")]}
# post-fix appendix redraws, 28 Sep evening (scripts/policy_deck/44, 45, 46)
for _n, _f in {21: "a6_where_rankings_work.png", 26: "a6_mean_concentration_error.png", 37: "a6_regional_level_maps.png",
               42: "a6_domain_ablation.png", 45: "a6_top20_child_iron.png", 46: "a6_malawi_b12_weights.png",
               47: "a6_malawi_b12_surrogate.png"}.items():
    PICTURE[_n] = os.path.join(FIG6, _f)
    assert os.path.exists(PICTURE[_n]), PICTURE[_n]
TEXT[21] = {
    "Where it works: iron, B12 and vitamin A where the deficiency is common, and everywhere the survey has several clusters per district.":
        "Where it works: B12, children's iron in Ghana, and vitamin A and women's iron in The Gambia, whose districts have more survey clusters than Ghana's or Malawi's.",
    "Why not elsewhere: an outcome under 1% prevalence has nothing to rank, zinc carries no district signal, and a single-cluster district gives the survey too little to be compared with.":
        "Why not elsewhere: the proxies carry no district signal for zinc, the survey cannot tell Malawi's districts apart on vitamin A, and a single-cluster district gives the survey too little to be compared with."}
TEXT[44] = {"The local staple price relative to the national level appears in 13 fits across 5 outcomes.":
            "The local staple price relative to the national level appears in 14 fits across 5 outcomes."}
TEXT[46] = {"What are the weights behind the best ranking? Malawi, women's B12": "What are the weights behind one of the best rankings? Malawi, women's B12"}
TEXT.setdefault(54, {})["Stable, yes; verified, no: coverage 38%, not 90%; a calibrated 90% interval is \u00b155% of the list. A worst-third call is right 52% of the time (chance 33%)."] = \
    "Stable, yes; verified, no: coverage 39%, not 90%; a calibrated 90% interval is \u00b155% of the list. A worst-third call is right 52% of the time (chance 33%)."
NOTES_SUB_EXTRA = {
    21: [("Above the line: iron, B12 and vitamin A where they are common enough to measure, and The Gambia everywhere. Below: the very rare outcomes, zinc, and the single-cluster countries.",
          "Above the line: B12, children's iron in Ghana, and vitamin A and women's iron in The Gambia. Below: zinc, Malawi's vitamin A, and the single-cluster districts.")],
    26: [("The model sits near one in most combinations even where it ranks well, which is the same message: order, not level.",
          "Calibrated for level, the model is below one in 16 of 18 combinations, every one except zinc, and beats the survey's regional average in 17 of 18."),
         ("The index is a ranker rescaled to the training spread, so its squared error is often at or above 1 even in combinations where it ranks well",
          "The figure shows the index calibrated for level (post-fix); the uncalibrated ranker, rescaled to the training spread, often has squared error at or above 1 even where it ranks well")],
    44: [("(34 of the 41 columns that recur in five or more fits are unanimous)", "(31 of the 37 columns that recur in five or more fits are unanimous)")],
    45: [("the starred ones (10 of twenty)", "the starred ones (9 of twenty)"),
         ("the rest is spread over the other 350 columns", "the rest is spread over the other 309 columns"),
         ("and the twenty-layer equal-weight composite scores about the same as the full index for a new country.",
          "and a twenty-layer equal-weight composite scores 0.24 for a new country, against 0.30 for the full index.")],
    46: [("the outcome-country combination with the best in-fill ranking of the 24 (Spearman 0.70; The Gambia child vitamin A is next at 0.69).",
          "the third-best in-fill ranking of the 24 after the survey fixes (Spearman 0.69; The Gambia children's vitamin A 0.71, women's vitamin A 0.70)."),
         ("districts that are poorer, grassier, further from permanent water and with more pigs have worse", "districts that are poorer, grassier and with more pigs have worse")],
}
NOTES_APPEND_EXTRA = {
    21: "\n\nPost-fix (28 Sep; benchmarks_v2_raw.csv, cell_master.csv): no in-fill combination is under 1 per cent any more (Malawi's vitamin A is 1.3 and 8.7 per cent), so rarity no longer explains the misses: women's vitamin A ranks 0.30 to 0.70 at 1.3 to 2.6 per cent. Country means 0.61 The Gambia, 0.44 Ghana, 0.27 Malawi; 6 of 18 at or above 0.5. The Gambia has 2.3 clusters per district against about 1.2 in Ghana and Malawi. Hollow points: not measurable.",
    26: "\n\nPost-fix (28 Sep; benchmarks_v2_raw.csv, rmse_sd, calibrated index): mean 0.90 against 1.01 for the survey's regional average. The Gambia women's vitamin A 0.71, Malawi B12 0.76, zinc 1.02 and 1.01.",
    37: "\n\nPost-fix (viz/oof_child_iron.csv, 28 Sep), children's iron: regions 0.77 (The Gambia), 0.76 (Ghana), 0.51 (Malawi, its 27 districts); at district level 0.37, 0.50 and 0.25.",
    42: "\n\nPost-fix (domain_ablation_loco.csv, 27 Sep): dropping climate costs 0.030, soil 0.029, anaemia 0.019; behind them come agriculture and livestock; only climate, soil and livestock have intervals clear of zero, and dropping household diet and consumption helps. A fixed climate-and-soil index transports at 0.37 at district level and 0.44 regionally, against 0.30 and 0.30 for the full set.",
    45: "\n\nPost-fix (index_importance_columns.csv, 28 Sep): pooled four-country index, children's iron, average status. The same twenty layers as before, lightly reordered; together they carry 21 per cent of the index.",
    46: "\n\nPost-fix (28 Sep): index refitted on 300 bootstrap resamples of Malawi's 87 surveyed districts (script 44); every weight shown keeps its sign in at least 99 per cent of refits. 18 of the earlier twenty remain: distance to permanent water and mild anaemia dropped out, child height-for-age and handwashing with soap came in.",
    47: "\n\nPost-fix (targets_v2.csv, 27 Sep): Spearman correlations across Malawi's 87 surveyed districts, -0.33 with modelled anaemia (was -0.40) and -0.45 with fish eaten (was -0.48). The anaemia link weakened after the survey fixes; the fish link holds. The model panel uses the headline columns without DHS (383 of 575), as the dashboard does.",
}

NOTES_SUB = {49: [("Children under five and non-pregnant women.", "Children aged 6 to 59 months and non-pregnant women.")]}
for _k, _v in NOTES_SUB_EXTRA.items():
    NOTES_SUB.setdefault(_k, []).extend(_v)
NOTES_APPEND.update(NOTES_APPEND_EXTRA)


def replace_runs(slide, pairs):
    hit = {a: 0 for a, _ in pairs}
    for sh in slide.shapes:
        if not sh.has_text_frame:
            continue
        for para in sh.text_frame.paragraphs:
            for r in para.runs:
                for a, b in pairs:
                    if a in r.text:
                        r.text = r.text.replace(a, b)
                        hit[a] += 1
    missing = [a for a, n in hit.items() if n == 0]
    if missing:
        sys.exit(f"run text not found: {missing}")


def main():
    anm = list(Presentation(ANM).slides)
    dst = Presentation(DST)
    last = None
    for num, after in PLAN:
        if num in TITLE_OF:
            assert m38.title_of(anm[num - 1]) == TITLE_OF[num], (num, m38.title_of(anm[num - 1]))
        new = m38.copy_slide(anm[num - 1], dst)
        if num in PICTURE:
            m38.replace_picture(new, PICTURE[num])
        if num in TEXT:
            m38.replace_text(new, TEXT[num])
        if num in NOTES:
            new.notes_slide.notes_text_frame.text = NOTES[num]
        if num in RUN_SUB:
            replace_runs(new, RUN_SUB[num])
        if num in NOTES_SUB:
            nt = new.notes_slide.notes_text_frame
            for x, y in NOTES_SUB[num]:
                assert x in nt.text, (num, x)
                nt.text = nt.text.replace(x, y)
        if num in NOTES_APPEND:
            nt = new.notes_slide.notes_text_frame
            nt.text = nt.text + NOTES_APPEND[num]
        idx = last if after is None else m38.index_of(dst, after)
        m38.move_after(dst, new, idx)
        last = idx + 1
        print(f"  v4-ANM {num:2d} -> after {m38.title_of(dst.slides[idx])[:60]!r}")
    dst.save(DST)
    print(f"v6: {len(dst.slides)} slides")
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import pptx_notes_size
    pptx_notes_size.main(DST, 18)


if __name__ == "__main__":
    main()
