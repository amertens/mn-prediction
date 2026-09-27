# Policymaker dashboard revamp — brainstormed plan, 27 September 2026

Goal: the shinyapps.io dashboard, shown live to policymakers, should be
focused, readable to non-ML readers, current (17–22 Sep results and their
interpretation), optimistic-but-accurate about what more data buys,
uncertainty-honest everywhere, and should demonstrate survey planning
(where to sample, and sample-size reduction) using methods borrowed from the
geospatial / NTD-elimination survey-design literature.

Starting point (read 2026-09-27): the app is already strong. One estimator
(deployment index), 11 panels organised around the reader's questions, plain
prose, per-cell forest plots with CIs, the ceiling, the learning curve, the
worst-fifth probability layer, honest-confidence calibration, WHO bands, the
anchor-and-rank design curves, exact back-projected importance with sign
stability, a searchable predictor catalogue, CIV rank intervals (child iron
only). So this plan is **targeted additions, not a rebuild**.

Effort marks: **[DAY]** = achievable in a focused day. **[TODO]** = needs a
longer window (new script/compute/dependency); park on the roadmap.

---

## A. Currency: newest results and their interpretation

The bundles were built 2026-09-13; the interpretation source of truth is now
`docs/findings/TWO_READINGS_2026-09g.md` (17 Sep, revision g supersedes f)
plus XV-01/02 (22 Sep). Nothing after 13 Sep is in the app.

- **A1 [DAY] Rebuild bundles on the current tables.** Re-run
  `dashboard/data-raw/05_build_protocol_v2_bundles.R`, then `smoke_test.R`
  and `test_server.R`. Traps: the two protocol-v2 rebuild traps
  (memory `protocol_v2`), Sierra Leone spelling ("Sierra Leone" with space),
  lake polygons.
- **A2 [DAY] External validation panel (the headline addition).** New bundle
  from `results/tables/external_validation/xv_transport*.csv`; new panel
  under "Can we trust it" (suggested title: *"Tested in six more
  countries"*): dot plot of per-country mean ρ with the country-block null
  band; Africa (Zambia 0.65–0.80 cells, Ethiopia, Sudan, Nigeria; pooled
  level 0.402, 12/12 positive, p = 0.001) and off-continent (India 0.431 —
  the only cells that clear their own per-cell null — Pakistan 0.381; pooled
  prevalence 0.394). Plus a Start-here hero/checklist row: *"Does it work
  outside these four countries? Yes — checked against WHO-deposited surveys
  in six countries on two continents, none used in training."* Carry the
  caveats from XV-01 §"does not": admin-1 aggregates, not the pre-registered
  P1 scoring, Nigeria's 6 units, the vitamin A reversal (fails in Africa,
  strongest in South Asia — reported, unexplained). SoilGrids-ties-iSDA
  lifts the Africa-only constraint → belongs in the roadmap panel (C1).
- **A3 [DAY if targets built / TODO otherwise] Malawi selenium + iodine.**
  Configured 2026-09-15 (memory `malawi_selenium_iodine_outcomes`; selenium
  index 0.43–0.45 in-fill, the strongest Malawi cell) but absent from
  `outcome_short`, caveats, and the builder. Needs: pipeline targets built,
  outcome labels + a caveat ("Malawi only; in-country evidence only"),
  builder + module label additions. Decide first whether policymakers should
  see 8 or 10 outcomes (open question Q4).
- **A4 [DAY] Interpretation & terminology pass.** Reconcile every tab's prose
  with TWO_READINGS-g and the 18 Sep terminology decisions used in the MNF15
  decks, so the talk and the dashboard use the same words for the same
  things. Also fold in the survey-report reconciliation caveats where
  missing (SL child file = anaemic subset; iron binaries; RBP rule as
  sensitivity).
- **A5 [DAY] Independent corroboration into the CIV tab.** The 2007 CIV B12
  eco-region check (rank agreement 0.95, one outcome, nine zones) and the
  WFP/MIMI Nature Food 2026 ADM2 vitamin A corroboration
  (`docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md`) as a short "two
  independent checks" card, stated with their limits.

## B. Focus and readability for non-ML readers

- **B1 [DAY] "If you only have five minutes" path.** A short ordered strip on
  Start here: 1) the map for your country, 2) how well it works, 3) what a
  survey planner should do next — each a `go_to()` button with one sentence.
  (A full interactive guided tour/stepper is **[TODO]** and probably
  unnecessary.)
- **B2 [DAY, after user decision] Trim the nav for this audience.** 11 panels
  is a lot. Proposal: fold "Predictor catalogue" and "What tracks which
  nutrient" into a single "The data" panel (catalogue as a tab of it);
  demote "What else was tried" further down Methods. Keep Plan a survey
  top-level (it is the pitch). Decide with the user (open question Q2).
- **B3 [DAY] Message-title pass.** Every card header becomes the takeaway
  sentence, not the axis description (e.g. "The worst fifth holds twice the
  burden a random fifth would" instead of "Burden reached"). Keep chart
  titles ≤ 12 words; move mechanics to the methods_note. Audit hover text
  for jargon (rho → "level skill" everywhere a reader sees it).
- **B4 [DAY] Presentation hygiene.** Colour-blind-safe check on YlOrRd +
  the blues (they're near-safe; verify with a simulator), consistent
  legend formats, mobile/projector font sizes, and a browser-tab
  title/favicon. Boundary simplification if map loading is slow on
  shinyapps (`rmapshaper::ms_simplify` in the builder) — quick but verify
  the lake-polygon trap.
- **B5 [TODO] One-page country brief.** A downloadable per-country PDF/HTML
  brief (ranking map, worst fifth, confidence, next-survey advice) built
  from `dashboard/report/` machinery. High value for handing to a minister's
  office; roughly 1–2 days including layout polish.

## C. Optimistic-but-accurate: "what more data buys"

- **C1 [DAY] New panel: "The roadmap" (under Can we trust it, or a
  Start-here card).** Assembles, from existing tables and findings:
  1. the learning curve ("each survey added buys ~`lc_step` in a country
     never seen; the curve has not flattened") — exists;
  2. the ceiling story: the model is at the survey-noise ceiling in
     `vc_at`/`vc_n` cells, so **better surveys (more clusters per district)
     raise the ceiling itself** — exists (`ceiling`, `headroom_by_cell.csv`);
  3. the measurability screen: scored against what the survey can measure,
     the index reaches ~0.49 vs 0.40 (memory
     `measurability_screen_and_cross_cutting`) — needs a small bundle;
  4. XV-02: global SoilGrids ties Africa-only iSDA → **the recipe is not
     Africa-bound**;
  5. XO-01: a same-nutrient biomarker from another population beats the
     index (0.454 vs 0.403) → collecting *any* related biomarker helps;
  6. pre-registered candidates waiting on country five (climate+soil,
     cs_top5, soft-threshold/decorrelated weights) — framed as "the next
     country is a prediction, not a retrofit".
  Tone: "the ceiling is set by survey noise, not by the idea; every survey
  any country runs improves every other country's map."
- **C2 [DAY] Update the Start-here checklist** with the external-validation
  row and a "What would make it better?" row (more surveys, more clusters
  per district, MICS microdata, a fifth country).
- **C3 [DAY] Scenario text box.** Tanzania incorporation status and what the
  pre-registration commits to (PREREGISTRATION_NEW_COUNTRIES doc), stated as
  plans, clearly separated from results.

## D. Survey-planning demonstrations (the centrepiece)

Methods borrowed from the geospatial / NTD survey-design literature (all
verified 2026-09-27):

- Adaptive geostatistical design — sample next where prediction variance is
  high or prediction is near a decision threshold: Chipeta et al. 2016,
  *Spatial Statistics* (https://www.sciencedirect.com/science/article/abs/pii/S2211675315001153);
  rolling malaria surveys in Malawi: Kabaghe et al. 2017, *PLOS ONE*
  (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0172266).
- Adaptive hotspot-finding (Bayesian-optimisation batch sampling; "same
  hotspot classification with a fraction of the sample size"):
  Andrade-Pacheco et al. 2020, *Scientific Reports*
  (https://www.nature.com/articles/s41598-020-67666-3) — schisto/LF, with
  Sturrock; closest in spirit to what our index enables.
- NTD elimination-survey design and threshold exceedance: Fronterrè et al.
  2020, *J Infect Dis* (https://pubmed.ncbi.nlm.nih.gov/31930383/);
  the geospatial-paradigm argument for NTD surveys: Diggle et al. 2021,
  *Trans R Soc Trop Med Hyg* (https://academic.oup.com/trstmh/article/115/3/208/6136278);
  MBG for more precise trachoma prevalence in elimination settings: Amoah
  et al. 2022, *Int J Epidemiol* (https://pubmed.ncbi.nlm.nih.gov/34791259/);
  MBG for trachoma elimination assessment: Harding-Esch et al. 2023, *PLOS
  NTD* (https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0011476);
  WHO-adjacent technical consultation: *Am J Trop Med Hyg* 2025
  (https://www.ajtmh.org/view/journals/tpmd/113/4/article-p930.xml).

The honest framing for all of D: these are **design demonstrations on the
four surveys and screening tools**, not a validated sampling protocol —
exactly the "Yes, in design / not yet a pilot" line the checklist already
takes. The exceedance-probability idea transfers directly because we already
compute resampling-ensemble probabilities (worst-fifth chance; CIV rank
intervals; VZ-01 checked their calibration).

- **D1 [DAY–2 days] "Where to survey next" planner, v1.** New tab section
  (inside Plan a survey): choose country + outcome; every district scored by
  a transparent weighted combination of (i) *uncertainty* — worst-fifth
  chance in the uncertain middle (0.2–0.8), later rank-interval width from
  D3; (ii) *threshold proximity* — planning prevalence near the WHO severity
  cut (the exceedance/elimination-survey idea); (iii) *population*; (iv)
  *never-surveyed* bonus. Three sliders + presets ("Confirm the worst",
  "Shrink uncertainty", "Cover the unknown"); map highlights the chosen k
  districts; each carries a one-line "why" ("large, never surveyed, and the
  model cannot tell which side of the 20% line it sits on"). Cite the AGD
  papers in a "where these ideas come from" note (D6). Uses only existing
  bundle columns → genuinely day-sized for the v1.
- **D2 [TODO, ~2–4 days] Retrospective validation of the planner.** The
  claim that makes D1 more than a heuristic: with the four surveys, simulate
  surveying only k districts chosen by the planner vs at random vs
  population-proportional; fit the index anchor on the chosen set; score on
  the held-back surveyed districts (machinery parallel to
  `model_augmented_survey.csv` / `survey_size_symmetric.csv` scripts).
  Deliverable: "picking districts this way, a survey half the size loses
  only X of ranking accuracy / classifies Y% of worst-fifth districts
  correctly" — the Andrade-Pacheco-style result on our data. New protocol_v2
  script + a findings note; pre-register the objective before running to
  keep the audit discipline.
- **D3 [TODO, ~1–2 days] Rank + exceedance intervals for every district.**
  Extend the CIV resampling ensemble (currently child iron only) to all four
  countries × outcomes and CIV × both candidates: rank_lo/rank_hi, worst-
  fifth chance for **unsurveyed** districts (today it exists only for
  surveyed ones), and a new **P(planning prevalence ≥ WHO threshold)**
  exceedance layer — the NTD-elimination quantity policymakers act on.
  Computation lives in the builder (offline), so shinyapps stays light.
  VZ-01's calibration check re-run on the extended ensembles.
- **D4 [DAY] Re-frame "Plan a survey" for a budget-holder.** Translate the
  x-axis fractions into respondents and clusters for a selected country
  (from survey metadata + deff ≈ 2.4); "I can afford N clusters" slider
  reading off error/ranking accuracy for the four designs; keep the capture
  curve as "what the cheap design gives up". Include the model-augmented
  and symmetric-size results *honestly* (augmentation does not beat direct
  estimates where clusters exist — say so; it inoculates against
  overclaiming). Add the D6 citations note.
- **D5 [TODO, exploratory] Threshold-classification calculator.** "To
  classify district X against the 20% band with 90% confidence: n clusters
  without the model, m with the model's prior" — Fronterrè-style elimination
  survey logic adapted to severity bands. Needs simulation and careful
  framing (the methodology is the least mature here); defer unless D2 lands
  cleanly.
- **D6 [DAY, folded into D1/D4] "Where these ideas come from" note** with
  the citations above, one line each, in the planning tab.

## E. Variable importance: searchable, cross-model, with uncertainty

- **E1 [DAY] Master searchable importance explorer.** One reactable over
  outcome × target × scope (pooled / each country's own fit / LOCO fits):
  plain-language label, domain, source, signed weight, share, sign-stability
  ("same sign in a of b fits"), searchable and filterable by domain/source/
  outcome, CSV download. Data: `weight_sources_raw.csv`,
  `weight_sources_cells.csv`, existing `importance_top` — mostly builder
  work. Honour the sign convention (positive = more deficiency; memory
  `signal_scan_sign_convention`).
- **E2 interim [DAY] Cross-fit ranges as the uncertainty.** Show each
  weight's min–max across the four in-country fits and across LOCO fits as
  an interval bar behind the pooled weight (data exist per-cell). Uncertainty
  = "does this predictor's role replicate", which is the question a reader
  actually has.
- **E2 full [TODO, ~1–2 days] Resampling intervals on weights and domain
  shares.** Reuse the 40-refit ensemble (index weights are correlation
  weights — cheap) to give every weight and every domain share an interval;
  dot-and-interval plots replace bare bars on "Leading predictors" and the
  domain scatter. Same builder run as D3 — do them together.
- **E3 [DAY, after E2] Domain view with uncertainty** (share + transport
  cost with across-cell intervals).
- **E4 [TODO, external] Finish plain-language labels.** 79 of 373 variables
  undocumented in the annotation worksheet (RA task, memory
  `variable_sheet_and_ra_briefs`); the explorer should show coverage and
  fall back to the raw column name gracefully.

## F. Uncertainty everywhere (cross-cutting audit)

- **F1 [DAY] Uncertainty inventory.** Walk every displayed number: does it
  carry an interval, a band, a "how sure" companion, or an explicit caveat?
  Known gaps: planning prevalence (point + rho band today; gets a range from
  D3), people affected (propagate the range), Admin-1 aggregates, CIV
  outcomes other than child iron (D3 fixes), unsurveyed districts' "how
  sure" (D3). Produce the gap list, fix the text-level ones same day.
- **F2 [DAY, after D3] Uncertainty on the map itself.** Optional hatching or
  desaturation for districts whose 90% rank interval spans more than half
  the country (a VSUP-lite); legend sentence: "pale districts are ones the
  model cannot place firmly."
- **F3 [DAY] Evidence that the intervals are honest.** Surface VZ-01's
  calibration check beside the worst-fifth calibration panel, in one
  sentence.

## G. Plumbing

- **G1 [each rebuild]** Builder → `smoke_test.R` → `test_server.R` → deploy
  checklist (shinyapps.io via `dashboard/deploy.R`).
- **G2 [TODO]** Boundary simplification / lazy tabs if new layers slow the
  free shinyapps tier (roadmap already lists the tier question).
- **G3 [TODO]** URL state for sharing a specific country/outcome view
  (roadmap open item; nice for follow-up emails to policymakers).

---

## Suggested sequence

**Sprint (3 focused days, demo-ready):**
1. A1 rebuild + A2 external-validation panel + A4 terminology pass.
2. D4 plan-a-survey reframe + D1 planner v1 + D6 citations + C1 roadmap
   panel + C2 checklist.
3. E1 importance explorer + E2-interim ranges + F1 audit + B1/B3 readability
   + G1 deploy.

**Longer to-dos, in value order:** D3 ensemble intervals (unlocks F2 and the
exceedance layer) → D2 retrospective planner validation (the quantitative
survey-planning claim) → E2-full importance intervals (same run as D3) → B5
country briefs → A3 selenium/iodine → D5 threshold calculator → E4 labels
(RA), G2, G3.

## Open questions for Andrew

1. **When is the policymaker showing?** Determines sprint vs full list.
2. **Trim which tabs?** Proposal: merge catalogue + nutrient-signal under
   one "The data" panel; keep Plan a survey top-level.
3. **Planner default objective** (D1): confirm-the-worst, shrink-uncertainty,
   or cover-the-unsurveyed? Proposal: uncertainty × population default with
   visible sliders.
4. **Show Malawi selenium/iodine** to policymakers, or keep the established
   outcome families for focus? (Selenium is the strongest Malawi cell but
   single-country.)
5. **Cluster costs:** any defensible per-cluster / per-assay cost figures to
   put currency on D4, or stay with "share of a full survey"? (Staying
   abstract is safer.)

*Written by Claude from a read of `dashboard/` and
`results/tables/{protocol_v2,policy_deck,external_validation}` on
2026-09-27. Not yet reviewed.*
