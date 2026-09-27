# Dashboard roadmap — big-ticket items

## Policymaker revamp (2026-09-27)

Before showing the app to policymakers. Plan and effort marks:
`docs/superpowers/specs/2026-09-27-policymaker-dashboard-revamp-design.md`.

- **Headline-tier fit (consistency fix).** Builder 05's deployment ranking had
  silently been fitting on all 575 columns including the 105 DHS-tier ones,
  while every quoted accuracy and the app's own text said DHS is held out. The
  fit now uses `drop_near_outcome_v2` under `open,survey_public` (383 columns);
  the catalogue's full-column rank matrix is kept for the small maps.
- **External validation (XV-01/02)** as a "Tested in six more countries" panel
  and Start-here row; bundles `xv_cells/xv_summary/xv_pooled`.
- **Stability ensembles (UE-01,** `data-raw/06_build_uncertainty_ensembles.R`):
  per-cell bootstrap of the deployment fit (CV-01/VZ-01 design, 200 draws) →
  rank ranges, planning-prevalence bands, and WHO-threshold exceedance for
  EVERY district (surveyed or not) and for the Malawi selenium/iodine cells;
  pooled-fit weight and domain-share ranges (150 draws) for the importance tab.
  Always labelled stability, never coverage: `stability_note()` quotes the
  38% LOCO coverage (VZ-01) wherever these appear.
- **Plan a survey** rebuilt as three tabs: an interactive district-visit
  prioritiser (objectives = confirm worst / shrink uncertainty / cover the
  never-surveyed, borrowed from adaptive geostatistical design; citations in
  Methods), the AR-01 size story with respondents/clusters translation and a
  national-anchor precision line, and **SP-01**
  (`scripts/protocol_v2/64_survey_planner_validation.R`, pre-registered in
  `docs/findings/SP-01_SURVEY_PLANNER_DESIGN_2026-09-27.md`): model-stratified
  district selection beats random consistently but modestly, PPS is worse than
  random, and the amendment measures what k-district selection does to the
  NATIONAL estimate (bias under informative selection, precision loss) — the
  survey's primary product stays a probability sample's job.
- **Importance explorer**: searchable table over every scope × outcome ×
  target (bundle `importance_all`), with resampling ranges on pooled weights
  and domain shares.
- **What more data buys** panel: headroom-vs-achieved per cell, learning
  curve, XV, XO-01 same-nutrient borrowing, selenium as the new-outcome story,
  pre-registered candidates.
- **Malawi selenium + iodine** ranked via the targets adapter
  (`data-raw/00_read_targets.R`); labels + caveats added; no WHO bands.
- Map explorer: rank-stability and WHO-exceedance layers, fade-unstable
  toggle, stability lines in the click-through; district profiles carry rank
  ranges and prevalence bands; CIV uses CV-01 stability for all six outcomes
  plus the corroboration card (2007 B12, WFP/MIMI, candidate agreement).
- Nav: catalogue + nutrient-signal merged under "The data behind it";
  banner headline is the talks' thesis line; terminology per the 18 Sep
  decisions (geographic interpolation, infectious disease, anaemia (modelled)).
- One-page country briefs (`data-raw/07_build_country_briefs.R` →
  `dashboard/briefs/`, download button on the map tab); deploy list extended
  (uncertainty_ensembles.rds, briefs/).

Rebuild order: 64 → 05 → 06 → 07 → smoke_test → test_server → deploy.

**Same-day follow-up (2026-09-27, second deploy):**

- **Boundary simplification done** (`data-raw/08_simplify_boundaries.R`,
  mapshaper keep=0.10/0.20, topology-preserving, keys and validity asserted;
  originals kept as `*_full.rds`): admin2 10.9→1.1 MB, admin1 3.7→0.8 MB.
  rmapshaper needed Rcpp ≥1.1.0, which was installed to a scratch `R_LIBS`
  because another project's running Rscripts held the shared Rcpp DLL.
- **Every predictor now has a plain-language name**
  (`scripts/protocol_v2/65_generate_plain_names.R` → `R/predictor_plain_names_generated.R`,
  462 systematic per-family translations; curated map still wins; catalogue
  shows `name_source` = curated/generated and the review sheet for the RA is
  `results/tables/protocol_v2/plain_names_generated.csv`).
- **URL state done** (app.R): the address bar tracks tab + map country/outcome/
  layer/level, district profile, CIV outcome, planner country/outcome/preset/k;
  a pasted link restores them (applied three times ~1.4 s apart, because a
  country update repopulates outcome choices on the next client round-trip and
  would otherwise overwrite the deep-linked outcome and k).

**CP-01 (2026-09-27, third deploy): the prevalence bands are now calibrated.**
`scripts/protocol_v2/66_conformal_prevalence_bands.R` (design first in
docs/findings/CP-01_...md): split conformal on out-of-fold residuals of the
deployed quantity (calibrated index + national anchor, 5-fold x 10 draws),
target of coverage = the survey's own measured district value. All 27 cells
LOO-cover at 0.90-0.93; the stability bands they replace covered a mean of
8% (0-23%) and are retired from the prevalence displays (they remain the
right object for rank firmness). Median honest half-width 26 pp — much of it
the survey's own single-cluster noise. Exceedance is now a conformal
predictive distribution: Ghana child vitA ">=80% above the moderate line"
drops 204 -> 67 districts of 260, the honest middle grows to 193. Builder 05
attaches `prev_cal_lo/hi`, `p_modplus_cal`, `p_severe_cal`, `conformal_n`;
the map layer, district click-through and district profiles use them with
`calibrated_note()`; bundle `conformal_prev` carries the per-cell checks.
This clears the LQAS calculator's main blocker (the residual sets in
`conformal_prev_residuals.csv` are its calibrated prior); remaining blockers
are per-district design effects and missing WHO conventions for
selenium/iodine.

**LQAS-style per-district sample-size calculator — design sketch (not built).**
Deliverable: per district, "n samples classify it against the WHO band with
95% assurance, with vs without the model prior", using the UE-01 anchored-
prevalence draws as the prior in a Bayesian-assurance LQAS (n, d) rule with
deff-inflated effective n. Blocked, in order, on: (1) the prior is a STABILITY
distribution (38% LOCO coverage) — it must be conformally recalibrated (the
machinery is block C of scripts/policy_deck/10_viz_tables.R) and its coverage
re-verified per cell before any "95%" is printed; (2) district-level design
effects are not estimable where districts hold one PSU, and LQAS operating
characteristics are sensitive to exactly that; (3) threshold classification is
a LEVEL claim and the survey's own regional average currently beats the model
on WHO-band accuracy (66% vs 64%); (4) several outcomes have no per-district
WHO classification convention at all. Path: recalibrate intervals → verify
coverage → simulate the (n, d) rules retrospectively on the four surveys
(pre-registered, SP-01 style) → ship behind a "design study" label. ~2–3 days.

## Rebuild on the corrected protocol (2026-09-13)

Ahead of the Micronutrient Forum stakeholder lunch the dashboard was rebuilt
so that it says what the decks and the manuscript say. Before this, none of
its 28 data bundles read the protocol-v2 tables: every number came from the
pre-audit pipeline, and the entry page carried two withdrawn results (the
anchoring gain and the "no proxy survives correction" claim).

Removed: Model diagnostics, Methods comparison (P1 to P8), Benchmarks (the old
leaderboard), National burden, Resolution and anchoring, the GBD placeholder
module and bundle, the Sierra Leone Admin-3 layer, the five estimator layers
(person-level SL, area, Fay-Herriot, BYM2, recipe), `app_public.R` (it
referenced a module that no longer existed), the technical annex report, and
the builders for all of the above.

Replaced or added: one estimator (the zero-tuning index fitted on all surveyed
districts and applied to every district) with four views; the worst-fifth
probability as the "how sure" layer; a planning prevalence anchored to the
national survey; exact per-predictor decomposition of any district's score;
How well it works (three tests, ceiling, learning curve, geostatistical
comparator); What the ranking buys (burden, calibration, WHO bands); Plan a
survey on the anchor-and-rank design; Cote d'Ivoire from climate and soil;
What drives the estimate on back-projected weights; the Predictor catalogue;
Start here and Methods rewritten with the closing checklist. Builder:
`dashboard/data-raw/05_build_protocol_v2_bundles.R`.

Still open from the list below: shinyapps tier, boundary simplification, URL
state, accessibility. New: write per-variable descriptions (the annotation
sheet holds mechanism templates per sub-domain, not definitions; the plain-name
map in `R/predictor_plain_names.R` covers 81 of 454 and the catalogue is the
worklist); compute the worst-fifth probability for
unsurveyed districts (would need refits over draws in the deployment fit);
decide whether to add the geostatistical prevalence with intervals as a
"measured number" layer for surveyed districts.

---

Plan for the larger dashboard work flagged after the UC Davis / BMGF bi-weekly
call (June 2026), updated 2026-08-27. Items below that refer to the recipe,
Fay-Herriot, BYM2, GBD or Admin-3 layers are historical. Smaller call items (#1 Start-here guide,
#2 Ghana use case, #3 misclassification layer, #5 biomarker caveats) are **done
and deployed**.

Status key: ✅ done · 🟡 in progress · ⬜ planned · ⛔ blocked (decision/data)

Live app: <https://amertens.shinyapps.io/micronutrient-burden/>

---

## #7 Stability + public/analyst split  (do first — gates external sharing)

### A. Public vs analyst split — ✅ resolved, the other way
- ✅ **Decided against a second app (2026-08-27).** The split assumed `app.R`
  was internal. It is not: shinyapps serves it with no authentication, so
  anyone with the URL can already open every tab. A second app would have
  bought presentation, not privacy.
- ✅ Instead the *one* app was reorganised: seven policy tabs at the top level
  and the six analyst tabs (Model diagnostics, Benchmarks, Resolution,
  Transportability, Côte d'Ivoire OOS, Methods comparison) grouped under a
  **Technical appendix** menu. Same reduction in noise, no second URL to keep
  in sync.
- `app_public.R` is retained and still passes the smoke test, but is **not
  deployed** and has no rsconnect record. Delete it if the appendix grouping
  proves sufficient.
- *Side note:* the account is at its **5-application shinyapps limit** anyway
  (`Regression_Explorer`, `Sampling`, `imic_boxplots`,
  `imic_descriptive_stats`, `micronutrient-burden`, `postnatal-bw-imputation`;
  terminated apps do not count). A second app would have required freeing a
  slot or upgrading.

### B. Stability
- ⬜ **Upgrade shinyapps tier** (Starter → Basic): more RAM, longer idle
  timeout, multiple instances. ⛔ billing/account decision. Still the single
  biggest fix — a deploy failed outright on 2026-08-27 with
  `Timeout during request` while starting instances, and succeeded unchanged on
  retry.
- ⬜ **Reduce startup cost.** The bundle is ~16 MB, of which **14 MB is two
  boundary files**, all loaded eagerly in `global.R`. Geometry is simplified at
  `dTolerance = 0.001` (~110 m), far finer than a national choropleth needs;
  ~0.01 plus topojson should cut it by most of an order of magnitude. Cheapest
  available win and independent of the tier decision.
- ⬜ **Prune deps** (91 packages) and set explicit idle-timeout / instance count
  in `deployApp`.

---

## #6 GBD / cross-model comparison  (highest strategic value)

### A. Covariate comparison (no data dependency)
- 🟡 **Predictor-family comparison table** across this project, GBD (DisMod-MR)
  and the WP nutrient-inadequacy models, with the zinc example, on the
  **Methods** tab. Family-level for now.
- ⬜ **Refine to per-nutrient covariate lists** by compiling from GBD/WP methods
  docs and engaging those teams.

### B. Estimate comparison (blocked on real GBD data)
- ⛔ **Source GBD Results Tool exports** (RA task open) — prevalence by
  country/year with uncertainty bounds.
- ⚠️ **The placeholder still ships.** `gbd_estimates.rds` (April, marked
  PLACEHOLDER) is bundled to the public app and `mod_gbd_compare` is still in
  `R/`, with only the `nav_panel` commented out at `app.R`. Either land the
  export or remove the stub — a hidden placeholder is a credibility risk if
  someone finds it.

---

## #8 Finer resolution + SSA expansion

### A. Sierra Leone Admin-3 (contained win) — ✅ done
153 chiefdoms via a Fay-Herriot / empirical-Bayes area model on chiefdom GEE
covariates. Builder: `dashboard/data-raw/_build_sl_admin3.R`. The map explorer
offers "Admin 3 (chiefdom)" only when Sierra Leone is selected. No chiefdom
population denominator, so counts are not shown at Admin-3.

### B. Whole-SSA prediction — deferred per Andrew
- ⬜ Hold until the core models are stronger. When ready: automate GEE
  extraction for any GADM country plus a country picker, **with trust /
  out-of-support flags front and centre**.

---

## Done since the June plan (2026-08-27)

- ✅ **One estimator across the app.** The recipe layer was documented as "the
  recommended primary district estimator" but referenced by no module — it
  could not be selected. The map opened on area-level HAL, district profiles on
  Fay-Herriot, and the two disagreed by up to 4.9 pp on the same district set.
  `DEFAULT_PRED_MODEL` in `global.R` now feeds every tab, per
  `docs/AREA_LEVEL_RECIPE_SPEC.md`.
- ✅ **Admin-2 key hygiene reaching the map.** GADM ships Malawi's lakes as
  Admin-2 polygons (Lake Malawi as 8 separate features) and repeats four real
  Traditional Authority names across districts. The area and recipe layers
  carried both, so the default map painted deficiency on open water and gave
  both `TA Lundu` polygons whichever value sorted first. All five layers,
  boundaries and population now agree on 239 Malawi units.
- ✅ **Scope statement replaced the blanket disclaimer.** "Not for citation or
  external distribution" on a page served without a login was either untrue or
  unfollowable, and a caveat people scroll past takes the real warnings with it.
- ✅ **Scenario uncertainty.** Cases averted now carries a range from each
  district's 95% interval; the what-if mode is separated and permanently
  labelled illustrative. Scenarios had also been pinned to the person-level
  table — modelling interventions across 75 of Ghana's 260 districts and
  reporting the total as national.
- ✅ **Printable country briefs** — `dashboard/report/`, one per country to PDF,
  self-contained HTML, markdown and CSV, sourcing `global.R` so the printed
  page and the screen cannot disagree.
- ✅ **Regression tests** — `data-raw/smoke_test.R` (restored) and
  `data-raw/test_server.R` (new): every layer × country × outcome through the
  map helpers, both UIs constructed, water/duplicate-key assertions, and a
  layer-sanity gate comparing each interval layer's median district against the
  national survey figure.

---

## Open issues found while doing the above

- ⚠️ **Four cells have no district-level signal.** Ghana and Malawi women's
  vitamin A, Sierra Leone child vitamin A and child iron return an *identical*
  estimate for every district — the near-null pre-filter drops every predictor
  and the fit reduces to an intercept. Two more vary by under a percentage
  point. The app and briefs now say so and draw no map, but **if these cells
  appear in the manuscript as subnational estimates, they should not.** This is
  a model finding, not a display problem.
- ⚠️ **BYM2 and single-PSU districts.** `survey::svymean` returns SE ~1e-17
  where an area has one PSU. `_build_fh_layer.R` rejects those and is fine.
  Porting the same guard to BYM2 collapsed 10 of 24 fits toward zero (Ghana
  women's iron: 0.10% against a 7.6% survey figure), because D enters INLA as
  the precision multiplier `1/D`, not as a shrinkage weight. Reverted; BYM2 now
  sits at 2 of 24 cells more than 4× from the survey. **The underlying
  single-PSU problem is still unfixed in BYM2** and needs the likelihood
  weighting reworked, not D rescaled.
- ⬜ **Nothing in the app is linkable.** No `enableBookmarking`, no query-string
  state. An advisor cannot send a colleague "the Ghana women's-iron map". The
  cheapest large multiplier on how far this travels.
- ⬜ **Accessibility.** No `aria-label`, alt text, or keyboard path into the
  Leaflet map; no `@media` rules, and a 320 px sidebar plus a fixed 650 px map
  is unusable on a phone. Often a procurement precondition for government
  users, who are frequently phone-first.
- ⬜ **Name-keyed joins.** Admin-2 *name* is the join key throughout. Four real
  Malawi TAs are therefore modelled and reported as one unit each. Closing it
  means keying on `Admin1 + Admin2` from covariate extraction onward.
- ⬜ **`dashboard/README.md` is stale** — it lists test files that were archived.

---

## Recommended sequence

1. **Now:** decide the shinyapps tier (#7B); it gates reliability for every
   external viewer. Geometry simplification can proceed regardless.
2. **Next:** resolve the four no-signal cells before they reach the manuscript.
3. **Then:** URL state, then the GBD stub (land or delete).
4. **Deferred:** SSA expansion (#8B), Admin1+Admin2 keying.

## Open decisions

- shinyapps **tier upgrade** (billing)?
- **Real GBD data** access, or keep #6 as covariate-comparison only?
- Delete `app_public.R`, now that the appendix grouping replaces it?
- Confirm **SSA expansion stays deferred**.
