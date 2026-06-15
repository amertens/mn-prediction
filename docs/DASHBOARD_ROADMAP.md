# Dashboard roadmap — big-ticket items

Plan for the larger dashboard work flagged after the UC Davis / BMGF bi-weekly
call (June 2026). Smaller call items (#1 Start-here guide, #2 Ghana use case,
#3 misclassification layer, #5 biomarker caveats) are **done and deployed**.
This roadmap covers the three big ones: GBD/cross-model comparison, stability +
public/analyst split, and finer-resolution / SSA expansion.

Status key: ✅ done · 🟡 in progress · ⬜ planned · ⛔ blocked (decision/data)

---

## #7 Stability + public/analyst split  (do first — gates external sharing)
The live demo kept disconnecting, and stakeholders want to start sharing.

### A. Public vs analyst split
- 🟡 **Two builds sharing the same modules** (this commit): `app.R` (full,
  internal) and **`app_public.R`** (lean, policymaker-facing: Start here, Map,
  District profiles, National burden, Decision value, Scenarios, Importance,
  Methods). Shared About/Glossary moved to `global.R`. Deploy script
  `deploy_public.R` (separate app name `micronutrient-burden-public`).
- ⬜ **Deploy the public URL** — held pending your go (creates a new public-facing
  app/URL). The internal app keeps its current URL.

### B. Stability
- ⬜ **Upgrade shinyapps tier** (Starter → Basic): more RAM, longer idle timeout,
  multiple instances — the biggest single fix for "reconnecting." ⛔ billing/account decision.
- ⬜ **Reduce startup cost:** lazy-load heavy bundles only when their tab opens;
  further simplify district geometries (or ship topojson); drop bundles the
  public app doesn't use.
- ⬜ **Prune deps/bundle** (~15 MB, 91 packages) and set explicit idle-timeout /
  instance count in `deployApp`.
- *Note:* the public build already helps — fewer simultaneous heavy renders.

---

## #6 GBD / cross-model comparison  (highest strategic value)
Tied to the IHME relationship; raised twice in the call (esp. zinc).

### A. Covariate comparison (no data dependency)
- 🟡 **Predictor-family comparison table** across this project, GBD (DisMod-MR),
  and the WP nutrient-inadequacy models, with the zinc example — added to the
  **Methods** tab (this commit). Family-level for now.
- ⬜ **Refine to per-nutrient covariate lists** by compiling from GBD/WP methods
  docs and engaging those teams (the call's intended deliverable).

### B. Estimate comparison (blocked on real GBD data)
- ⛔ **Source GBD Results Tool exports** (RA task already open) — prevalence by
  country/year with uncertainty bounds.
- ⬜ Once data lands: populate the (currently hidden) GBD tab — model-vs-GBD
  divergence + "population mischaracterized," with careful outcome/definition
  harmonisation (GBD anaemia vs biomarker iron; sparse GBD zinc/B12).

---

## #8 Finer resolution + SSA expansion
### A. Sierra Leone Admin-3 (contained win) — ✅ done
- ✅ **Built and wired** (153 chiefdoms vs 14 districts). Chiefdom prevalence is
  reconstructed from individuals (`sierraleone_cluster_to_admin3.rds` maps
  `gw_cnum` → `NAME_3`, survey-weighted binary deficiency), then a Fay-Herriot /
  empirical-Bayes area model on chiefdom GEE covariates (greedy-decorrelated,
  standardized) extends to all 153 chiefdoms with a 95% interval (tight where
  surveyed, wider elsewhere). Builder: `dashboard/data-raw/_build_sl_admin3.R`
  → `dashboard/data/admin3_predictions.rds` (918 rows = 153 × 6 outcomes) +
  `admin3_boundaries.rds` (GADM level-3).
- ✅ **Dashboard:** `global.R` loads `admin3_pred`/`admin3_bnds` (+
  `admin3_countries`); `get_country_admin3()` helper in `fct_helpers.R` mirrors
  the Admin-2 schema; `national_aggregate()` falls back to an unweighted mean
  when no population denominator exists. Map explorer shows "Admin 3 (chiefdom)"
  in the level toggle **only when Sierra Leone is selected** (dynamic
  `updateRadioButtons`); headline, caption and district panel handle the
  no-population case. All display layers (prevalence, CI width, WHO class,
  misclassification, modeled−survey) verified via `testServer`.
- *Note:* no chiefdom-level population denominator, so population-at-risk counts
  are not shown at Admin-3 (prevalence, intervals and WHO class only).

### B. Whole-SSA prediction (deferred per Andrew)
- ⬜ Pipeline can predict any SSA country (CIV is the example) but hold until the
  core models are stronger. When ready: automate GEE extraction for any GADM
  country + a country picker, **with trust / out-of-support flags front-and-
  centre** (transported predictions are no better than chance for some outcomes).

---

## Recommended sequence
1. **Now:** public/analyst split + stability code fixes (this commit) → reliable,
   shareable app. Fire off the GBD-data request so #6B isn't blocked later.
2. **Next:** GBD covariate table refinement (#6A); wire #6B when data arrives.
3. **Then:** Sierra Leone Admin-3 (#8A).
4. **Deferred:** SSA expansion (#8B).

## Open decisions
- shinyapps **tier upgrade** (billing)?
- **Real GBD data** access, or keep #6 as covariate-comparison only for now?
- Deploy the **public URL** now, or keep `app_public.R` staged?
- Confirm **SSA expansion stays deferred**.
