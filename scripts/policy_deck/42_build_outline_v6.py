"""Write docs/slides/MNF15-talk-outline-2026-09-28.md (and .docx) for the v6 deck: the new
version of MNF15-talk-outline-2026-09-27.md. The script and backup for every main slide
are read from the v6 speaker notes, so the document and the deck cannot drift apart; the
visuals, sections and open items are written here.

    python scripts/policy_deck/42_build_outline_v6.py
"""
import re
import subprocess

from pptx import Presentation

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
DECK = ROOT + "docs/slides/MNF15-talk-2026-09-v6.pptx"
MD = ROOT + "docs/slides/MNF15-talk-outline-2026-09-28.md"
PANDOC = "C:/Users/andre/AppData/Local/Pandoc/pandoc.exe"

# keyed by slide title (apostrophes straightened), so inserting a slide cannot shift them
SECTION = {"Title": "Motivation", "Conceptual framework of iron deficiency": "Data assembly",
           "How do we get honest estimates of model performance?": "Methods",
           "But can proxy models recover national-level prevalence estimates?": "National prevalence",
           "How far off is the predicted prevalence, outcome by outcome?": "District prevalence",
           "How well does it rank districts, nutrient by nutrient?": "Ranking",
           "Can it rank Ghana's districts without Ghana's survey?": "Out of country",
           "What variables drive model predictions": "Variable importance",
           "Which deficiencies can it map?": "Policy relevance and next steps", "Conclusions": "Conclusions"}
VISUAL = {
    "Title": "Title, authors, the three logos and the Gates acknowledgement.",
    "Introduction and study objectives": "Sonja's bullets with the two specific objectives on the left; on the right the WHO VMNIS map of the latest national biomarker survey in each country, with counts by nutrient. Source as a footnote.",
    "Talk outline": "The five questions the talk answers, in talk order. The closing slide answers them in the same order.",
    "Conceptual framework of iron deficiency": "Sonja's published framework (Hess et al. 2023) with its footnotes.",
    "Which of these causes can public data measure?": "The same framework with every box coloured by public availability: dark green intrinsically ecological, light green individual data linkable as a district aggregate, amber partly, grey not public, teal the survey outcome.",
    "575 proxy variables from 29 conceptual domains": "Andrew's slide: predictors by source (bars) and the eight domain groups with their counts.",
    "Which surveys did we use to train the models?": "The four countries on a map and a table: survey, districts with blood samples as n of N, nutrients measured.",
    "How do we get honest estimates of model performance?": "Cross-validation schematic (five rounds) and three maps: one district in five hidden, a whole region hidden, a whole country hidden. Orange is the held-out test data.",
    "Why not a more complex machine-learning model?": "Three lines of text; below, six methods compared under the three hold-outs (the index, a SuperLearner of twelve, random forest, boosting, lasso, elastic net), post-fix.",
    "But can proxy models recover national-level prevalence estimates?": "Left: 15 country-outcomes, each country predicted from national indicators with its own survey left out, the gap labelled. Right: typical error in WHO's national survey panels, model against the average of the other countries.",
    "How far off is the predicted prevalence, outcome by outcome?": "Typical error in each district's predicted prevalence, each district hidden, for 18 combinations: the model calibrated for level, the survey's regional average, and the national estimate (X).",
    "How well does it rank districts, nutrient by nutrient?": "One panel, per nutrient and country: inside the country (circle, 95% interval), whole country held out (triangle), the survey's regional figures scored fairly (diamond). National prevalence in brackets. Zinc (Malawi women) has no held-out test.",
    "What does that accuracy mean on the ground?": "Three icon rows: pairs in the survey's order (model, regional figures, no survey); with a budget per person, rank by the model; with a budget in districts, start with the most populous. ALTERNATIVE to the next slide: keep one.",
    "Which districts should a programme go to first?": "The targeting figure: share of deficient people reached under each rule, inside a surveyed country and with no survey, for a budget of districts and a budget of people. ALTERNATIVE to the previous slide's budget rows: keep one.",
    "Where do the young children live?": "Ghana, children's iron: left, districts coloured by the model's predicted priority (highest-rate third: 33% of districts, 63% of the land); right, the same districts as circles sized by young children (that third holds 27%; the most populous fifth 41%).",
    "Is Côte d'Ivoire like the places we learned from?": "One row per country: each district's distance from the climate and soil conditions of the surveyed districts, with the edge of the training range. Côte d'Ivoire 33 of 33 inside; each surveyed country, held out, below it.",
    "Can it rank Ghana's districts without Ghana's survey?": "Ghana's children's iron: the survey's ranking and the model's with no Ghanaian data. Bars: each country held out in turn.",
    "What about countries with no survey?": "Left: Côte d'Ivoire's children's iron ranking from the climate-and-soil model. Right: the six countries whose WHO regional survey results checked the model.",
    "What variables drive model predictions": "Three findings on the left (climate and soil carry the cross-border ranking: dropping either makes it worse in all four countries; climate, soil and a satellite imagery summary are in the top five of all 22 cross-country and 21 of 24 single-country models; only B12 has a map of its own). Right, children's iron: top-ten domains in Ghana's, Malawi's and the four-country model, and the four-country model's top ten layers with their direction.",
    "Which deficiencies can it map?": "Traffic-light table by nutrient: inside a surveyed country, whole country held out, four African countries in WHO data, and whether to use it. One line underneath: only the B12 map is specific to its nutrient.",
    "What does this mean for policy? The vitamin B12 model": "Four cards: 70 to 75 of 100 pairs (65 to 71 held out, average B12 level); about 2x capture inside a surveyed country; about 9 points typical prevalence error; 4 of 4 WHO checks.",
    "Next steps": "Three asks as cards, the thesis banner and the QR code.",
    "Conclusions": "The five questions answered, each with its number, and the QR code.",
}


def mmss(s):
    return f"{s // 60}:{s % 60:02d}"


p = Presentation(DECK)
slides = list(p.slides)
rows, blocks, words_by = [], [], {"Sonja": 0, "Andrew": 0}
appendix = []
in_app = False
for i, s in enumerate(slides, 1):
    title = s.shapes.title.text_frame.text.strip().replace("\n", " ").replace("\u2019", "'") if s.shapes.title is not None else ""
    if title == "Appendix":
        in_app = True
        continue
    if in_app:
        appendix.append(title)
        continue
    n = s.notes_slide.notes_text_frame.text.strip()
    m = re.match(r"^(SONJA|ANDREW),? about (\d+) seconds(?: \(\d+ words\))?\.\s*", n)
    assert m, (i, n[:80])
    who, sec = m.group(1).title(), int(m.group(2))
    rest = n[m.end():]
    spoken, _, backup = rest.partition("Backup, if asked:")
    spoken = spoken.strip()
    w = len(spoken.split())
    words_by[who] += w
    if i == 1:
        title = "Title"
    assert title in VISUAL, ("no visual description for", title)
    rows.append((i, who, SECTION.get(title, ""), title, VISUAL[title], sec, w))
    sec_head = f"## {SECTION[title]}\n\n" if title in SECTION else ""
    quote = "\n>\n".join("> " + para.strip() for para in spoken.split("\n") if para.strip())
    blk = f"{sec_head}### Slide {i}. {title} ({who}, {mmss(sec)}, {w} words)\n\n**Visual.** {VISUAL[title]}\n\n{quote}\n"
    if backup.strip():
        blk += f"\n*Backup, if asked:* {backup.strip()}\n"
    blocks.append(blk)

tot_s = sum(r[5] for r in rows)
tot_w = sum(r[6] for r in rows)
spine = "| # | Speaker | Section | Title | Time | Words |\n|----|---------|-------------|----------------------------|-------|-------|\n" + "\n".join(
    f"| {i} | {who} | {sec_name} | {t} | {mmss(sec)} | {w} |" for i, who, sec_name, t, _, sec, w in rows)

md = f'''---
title: "MNF15 joint talk: outline and script for deck v6"
subtitle: "Micronutrient Forum 2026, Accra. Sonja Hess and Andrew Mertens. Session on national nutrition surveys, Thursday 1 October."
date: "2026-09-28"
---

# 0. What this version is

**Deck.** `docs/slides/MNF15-talk-2026-09-v6.pptx`: {len(rows)} main slides including the title, then the appendix ({len(appendix)} slides). It follows the order Andrew set on 28 September: motivation, data assembly, methods, national prevalence, district prevalence, ranking, out of country, variable importance, policy relevance and next steps, conclusions. It replaces the 27 September outline (an 11-slide spine built around three uses); the deck has since followed Andrew's edits.

**Time.** {tot_w:,} spoken words (Sonja {words_by["Sonja"]:,}, Andrew {words_by["Andrew"]:,}): about {tot_w / 140:.0f} minutes at 140 words a minute; the timings in the notes add up to {mmss(tot_s)}. The slot is 15 minutes, "plan for like 13" (Sonja, 18 September check-in), with questions pooled at the end of the session.

**Slides.** The Forum's guideline is at most one slide per two minutes of speaking, 7 or 8 for this slot. v6 has {len(rows)}: 22 once Andrew keeps one of the two targeting slides (13 and 14), including the two figures added on 29 September (the population cartogram and the Côte d'Ivoire applicability check). The pace (about 40 seconds a slide) is the main risk. Section 5 gives the cuts to about 10 minutes and 16 slides.

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

{spine}
| | | | **Total** | **{mmss(tot_s)}** | **{tot_w:,}** |

# 3. Slide by slide

The script is the speaker notes of v6, as they will appear in PowerPoint (18 pt). Sonja's lines are suggestions for her to rewrite. "Backup, if asked" carries the sources and the numbers for questions.

{chr(10).join(blocks)}
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

{chr(10).join(f"{k}. {t}" for k, t in enumerate(appendix, 1))}

# 7. Open items before 1 October

1. **Main talk:** every figure is now post-fix. The Côte d'Ivoire map was rebuilt on 28 September and barely moved (rank agreement 0.999).
2. **Appendix:** every figure was redrawn on post-fix data on 28 September except the geostatistics comparison (a 4 to 6 hour rerun) and the survey-planner line in the planning notes. The prevalence bands need script 66's national anchor file rebuilt on the fixed targets (a dashboard input too, so coordinate with the dashboard session).
3. **Slide 6** (Andrew's): the title overlaps the chart in the render; check in PowerPoint.
4. **Sonja's slides**: her introduction and framework-coverage notes were shortened (to about 40 and 45 seconds); she should read them.

# 8. Decisions

Made 28-29 September: keep 20 slides, plus the two figures added on 29 September; keep the five-question outline; show both targeting versions (13 and 14) for Andrew to choose; the nutrient finding stands (the session covers all micronutrients, not iron); slide 9 redrawn on the post-fix runs.

Still open:

1. Slide 13 (the budget rule in words) or slide 14 (the same result as a figure)?
'''
md = re.sub(r"\n{3,}", "\n\n", md)
open(MD, "w", encoding="utf-8").write(md)
print("wrote", MD, f"{len(md.split()):,} words; main slides {len(rows)}, spoken {tot_w}, {mmss(tot_s)}")
subprocess.run([PANDOC, MD, "-o", MD.replace(".md", ".docx")], check=True)
print("wrote", MD.replace(".md", ".docx"))
