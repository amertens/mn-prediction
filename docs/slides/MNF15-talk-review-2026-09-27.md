# MNF15 joint oral: CURRENT vs V2, an independent review (27 September 2026)

CURRENT = `docs/slides/MNF15-talk-2026-09.pptx` (body slides 1-15, with Sonja's two slides inserted after slide 3). V2 = `docs/slides/MNF15-talk-2026-09-v2.pptx` (body slides 1-12). Evidence: speaker notes, slide images, the 18 September check-in, and lookups in `results/tables/`.

## 1. Verdict

V2 is the better base for both audiences (about 75 per cent confident). For the policy majority it fits the slot (11.3 minutes of notes against CURRENT's 14.6 before Sonja's two slides), opens on a decision a programme officer recognises, and ends on uses, limits, a QR code and asks. For the experts it answers Sonja's 18 September objection ("1% better than looking at the neighbouring district?") by saying the neighbour map ties the model inside Ghana (0.53 each) and turning to countries without a survey, where public data earn their keep. But it adds overclaims experts will catch (the slide 7 title, "36 of 41 tests", "better than regional averages in countries it never saw", a clusters-per-district rule) and drops three things CURRENT does better. Neither deck meets the Forum guideline of 7-8 slides for a 15-minute slot (V2 11, CURRENT 17 with Sonja's two).

## 2. CURRENT

**Strengths**
- The three-map opener (slide 2) and the survey-team table (s3). The table is the only visible credit to the survey teams, and the Ghana team will be in the room.
- Accuracy for each nutrient with its own ceiling, women and children separate (s9), as Sonja asked.
- Domain weights (s10). This is the only evidence for the closing claim that most of the signal is in climate and soil. B12 leaning on infants fed meat or fish is the nutrition hook Sonja liked.
- The pooling curve behind the main ask (s14).

**Weaknesses**
- About 16 minutes long. Andrew first speaks about four minutes in, and five method and accuracy slides run back to back (s6-s10).
- The framework appears twice: Sonja's published figure, then the s4 redraw. The redraw has different tiers, and several boxes are missing (obesity, family planning, knowledge, built environment). Yet the script calls it "the same framework".
- The s8 "ruler" is the slide Sonja could not follow on 18 September. It pools 18 cells, is unclear on within-country versus across-border, and its 70 per cent "best nutrients" row is a post hoc pick.
- Five different number scales; no dashboard in the body; captions at 9-10 pt; title slide dated 18 September.

**Accuracy problems**
- **s5:** says 0.52 (`benchmarks_v2_cells.csv`: 0.529). The notes say the model "adds" public data to geography, but it uses no neighbour information and ties the neighbour map (0.529 against 0.526). Atiwa East "56th" comes from a single-draw file (`fig7_..._level.csv`); V2's ten-draw rebuild says 65th, so one district moves about ten places between draws. "Sends its first team somewhere else" is weak: 12th of 75 is inside any top-15 list.
- **s9:** "B12 and vitamin A map best, essentially at the ceiling." Vitamin A's 0.51 held out is The Gambia and Ghana only (0.61, 0.77, 0.37, 0.31); Malawi and Sierra Leone score 0.26, 0.24, 0.04, -0.02; externally in Africa vitamin A is weakest (0.23 over 14 cells) with both clear misses (-0.43); in Ghana children the regional average beats the model (0.41 against 0.28). "The immunoassay reads systematically lower" cannot explain a folate ranking collapse: a constant offset leaves ranks unchanged.
- **s11:** notes say dark means firm; the legend says dark means "could move 9 places". It shows stability, not uncertainty.
- **s13:** the tick on "where the next survey should sample" is contradicted by SP-01 (model-guided 0.362, random 0.357 at half the districts); "Today: an equal spread" ignores PPS sampling within strata; the tick on district prevalence ignores CP-01's ±25-point bands; "targeting likely better" hides that the worst fifth reaches 22 per cent of deficient people against 20 (perfect 48).
- **s14-15:** "every survey added improves every other country's map": the gain so far is in The Gambia and Ghana only. "Better than regional averages in a country never surveyed": held out the model averages 0.30 (22 cells); the in-country regional average averages 0.31.
- **s12:** 0.95 over nine zones on a steep north-south B12 gradient, with no latitude-only comparison.

## 3. V2

**Strengths**
- The order follows a programme officer's questions: the gap, Sonja's own framework with coverage badges (s3-4), one method slide, Ghana with its survey (s6) and without it (s7), no survey at all (s8), which nutrients, the next survey, uses, asks.
- s6's honest "so why bother?" turn; s7's held-out test where the answer is known; s8's external check on WHO VMNIS regional results (Sonja's own 18 September suggestion), the biggest gain in trust.
- s9's "Use it?" column turns accuracy into decisions; s10-11 speak to the session theme ("data collection, data targeting... that's what people want to hear").
- Numbers reproduce the benchmark cells; today's SP-01 and CP-01 are folded in; appendix s14 indexes expected questions.

**Weaknesses and losses**
- Loses the survey-team credit, the domain weights (so s12's "much of it in climate and soil" is asserted, never shown) and the pooling curve.
- s4 still carries DRAFT badge calls Sonja must confirm. s8 maps children's iron, while Côte d'Ivoire's own 2007 check is B12 and only spoken.
- s10's single "ceiling" of 0.54 is exceeded by cells the deck shows (B12 0.63; The Gambia 0.60 held out). "Put two or three clusters in a district, not one" is a projection presented as a rule to a session with survey designers; at fixed budget it means fewer districts, which are the model's effective sample, and Sierra Leone (4.3 clusters per district) is the weakest country. "Before a full survey" can be heard as "instead of".
- The QR code opens a dashboard whose start tile says 0.38 for a new country (climate-and-soil index, domains chosen on the same 22 cells); the talk says 0.30 (full index) and never says which model made the Côte d'Ivoire map.
- Appendix s43 looks blank on the contact sheet.

**Accuracy problems and overclaims**
- **s7 title:** "ranks Ghana just as well" rests on the single best cell (0.554). The slide's own bars show 0.35 for Ghana across nutrients, and 0.30 is the mean over all 22 cells. The "guessing" line is the 95th percentile of the null, not zero.
- **"36 of 41 tests":** it counts the average-status and prevalence versions of the same survey outcome as separate tests (12 + 21 in Africa), and mixes two soil sources for the 8 South Asian cells. The accurate version (`xv_transport_pooled.csv`):
  - Africa, 12 of 12 on average status (mean 0.40), 17 of 21 on prevalence (0.33).
  - South Asia, 7 of 8 (0.39).
  - Regions number 6 to 15 per country. Most single cells do not clear their own null; the pooled test does.
- **"Iron pointed the right way in every one of the six":** true of country means only (13 of 14 iron cells; Sudan women's prevalence is -0.06).
- **s9 header, "Six external countries":** it shows the four African countries' numbers. That hides Pakistan's B12 miss (-0.26): B12 is 4 of 5 across six countries. The held-out column is "measurable combinations only" (vitamin A 0.54/0.49 is The Gambia plus Ghana) and does not say so.
- **"No worse than inside our own four countries":** admin-1 regions against districts is not like for like.
- **s12:** repeats "better than regional averages ... in countries it never saw" (0.30 against 0.31).

## 4. Borrow list (into V2)

1. **Survey-team credit** (CURRENT s3) as a strip on the framework slide. Necessary in Accra.
2. **A domain-weight strip** (CURRENT s10): three bars (climate, soil, satellite: just under half the weight in every panel) plus the B12 and infant-feeding line. The V2 outline sends this to the short oral, but the talk's first conclusion needs its evidence on screen, and the audiences differ.
3. **Pooling curve** (CURRENT s14) as an inset on the closing slide, "on average", without the projected point.
4. **Per-nutrient ceilings** (CURRENT s9's grey bars) in place of the single 0.54.
5. **The 2007 B12 scatter** (CURRENT s12) as an inset on V2 s8, with the map's nutrient named.
6. **"Today: ... With the model: ..." contrasts** (CURRENT s13) inside V2 s11's three columns. It is the slide people will photograph.

## 5. Best version (plan 13 minutes; about 11.5 minutes scripted)

| # | Speaker | Title | Visual | s | Source |
|---|---|---|---|---|---|
| 1 | Sonja | Title | title, date fixed | 15 | V2 s1 |
| 2 | Sonja | A national survey measures regions. Programmes act on districts. | three maps, 185 of 260, 62 of 75, the question | 60 | V2 s2 |
| 3 | Sonja | We started from what causes deficiency (build) | Hess 2023 framework, then badges; survey-team strip | 75 | V2 s3-4 + CURRENT s3 |
| 4 | Andrew | Learn where blood was drawn, apply everywhere, test on what it never saw | four steps and testing banner; domain strip | 70 | V2 s5 + CURRENT s10 |
| 5 | Andrew | Ghana, with its survey and without it (build) | V2 s6 maps and dots; click: held-out map, 0.55, country bars | 130 | V2 s6-7, retitled |
| 6 | Andrew | No survey at all: a first map, checked against surveys it never saw | Côte d'Ivoire map (nutrient and model named), world map with corrected counts, B12 inset | 80 | V2 s8 + CURRENT s12 |
| 7 | Andrew | Some deficiencies can be mapped from public data; others cannot yet | corrected table | 60 | V2 s9 |
| 8 | Andrew | What the model can do for you, and for the next survey | three columns with "today" contrasts; anchor line; cluster trade-off spoken; QR | 110 | V2 s11 + s10 + CURRENT s13 |
| 9 | Andrew | What we found, and what we ask | V2 s12 with the pooling inset | 55 | V2 s12 + CURRENT s14 |

- Timing: 655 seconds plus about 30 for handovers, leaving about 1.5 minutes of buffer.
- V2 s10's chart moves to the appendix (it is already s19).
- The guideline allows 7 or 8 slides for the 15-minute slot. Folding slide 7 into slide 6, whose table already carries the external column, brings the deck to 8.

## 6. Top ten edits, in priority order

1. **Restructure to the nine slides above** (merge V2 s6-7 into a build; fold s10 into s11-12), then two timed run-throughs. *3-4 h.*
2. **Retitle V2 s7** "Ghana held out entirely: children's iron still ranked (0.55)"; say "one of our better cases" before the bars. *10 min.*
3. **Correct the external-check counts** (Africa 12 of 12 and 17 of 21, South Asia 7 of 8; iron "on average in all six"; fix the s9 header; add "measurable combinations only"). *30 min.*
4. **In both decks, replace "better than regional averages in countries it never saw"** with "with no survey at all, about as good as a survey's own regional averages (0.30 against 0.31)". *10 min.*
5. **Soften the survey-design claims:** per-nutrient ceilings (or label 0.54 an average); clusters per district as a projected trade-off against fewer districts; "alongside", not "before", a full survey. *45 min.*
6. **Align talk and dashboard** on one "new country" number and say which model made the Côte d'Ivoire map (hand the dashboard change to the session that owns it). *30 min.*
7. **One vitamin A message:** good in The Gambia and Ghana, weakest in the African check, not yet for supplementation targeting. *15 min.*
8. **Restore the survey-team strip, domain strip and pooling inset.** *2-3 h.*
9. **Sonja's sign-offs by Monday 28 September:** the DRAFT badges; retiring her bullet intro slide (if kept, "targeted machine learning" becomes "machine learning guided by conceptual frameworks": Berkeley statisticians will hear TMLE). *30 min of hers, 20 min edits.*
10. **Consistency sweep:** 0.53 everywhere; Atiwa East 65th; bands "about 25 points"; ceiling 0.54; "regional average", never "interpolation"; the blank s43. *1 h.*

## 7. Pooled Q&A: rehearse or have ready

Sonja takes biomarkers, framework and survey design; Andrew takes methods. Keep answers to 20 seconds and print the V2 s14 index.

- **"Isn't it just geography, or latitude?"** Inside a country, largely yes; across a border there are no neighbours. The DHS-style geostatistical model loses on order (0.27 against 0.38) but wins on levels (9.2 against 10.7 points). Before Thursday, compute a latitude-only ranking for the nine Côte d'Ivoire B12 zones and the held-out countries.
- **"Is your regional-average comparator handicapped?"** It leaves out the district's own data, which in regions of 3-5 districts hands the worst districts the lowest anchors (CF-01, `comparator_fairness.csv`). Repairs that never see the district (shrinkage, split halves) barely move it; only including the district's own answer does, and that leaks. Say "what a survey that missed the district would know"; never "no better than chance".
- **"Your truth is one cluster of about 12 children."** That is the ceiling on what anyone can show, and why the model gives ranks, not levels.
- **"Why not machine learning?"** Fourteen alternatives, 1,304 fits, none better by more than 0.03; 14-87 districts per country.
- **"Leakage?"** Hold-outs by district, region or country, selection inside the folds; anaemia layers external; Malawi's MNS clusters removed.
- **Biomarkers and timing.** Know each country's inflammation adjustment, cut-offs and RBP-to-retinol rule; the Côte d'Ivoire map uses current layers against a 2007 survey (Sonja raised this).
- **Prevalence and targeting.** About 40 draws gives a median district error near 9 points with ±25-point bands; have the no-model baseline (national figure everywhere) ready. The worst fifth reaches 22 per cent of deficient people against 20; burden follows population.
- **Survey planning (SP-01).** Model-guided choice about equals random; population-weighted is worse; design weights keep the national estimate unbiased.
- **Vitamin A, Africa against South Asia.** Unexplained; do not speculate.

*Checked afterwards against the V2 author's outline (`MNF15-talk-outline-2026-09-27.md`).* It caught, and I have added, the inverted s11 note, the PPS point, the folate offset logic, the redrawn framework, the Binduri sentence, the TMLE wording and the talk-dashboard mismatch. It still endorses "36 of 41", the slide 7 title, the single ceiling and the two-or-three-clusters rule, and V2 s12 turned its conclusion into "better than regional averages ... in countries it never saw". I disagree on all five.
