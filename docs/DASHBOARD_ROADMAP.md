# Dashboard roadmap — big-ticket items

Plan for the larger dashboard work flagged after the UC Davis / BMGF bi-weekly
call (June 2026), updated 2026-08-27. Smaller call items (#1 Start-here guide,
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
