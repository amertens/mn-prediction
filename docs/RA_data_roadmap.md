# RA data roadmap — micronutrient proxy-prediction project

Consolidated from three landscape investigations (2026-06). Companion docs:
[`micronutrient_survey_candidates.md`](micronutrient_survey_candidates.md) ·
[`gee_landscape_new_domains.md`](gee_landscape_new_domains.md) ·
[`gee_admin3_audit_findings.md`](gee_admin3_audit_findings.md).

Cross-cutting reality check: the documented driver of weak cross-country
transport is the **cross-survey biomarker LEVEL offset**, not missing predictors
or geography. Prioritise data work that attacks *comparability*, then *country
range*, then *predictors*. Almost every candidate survey is **admin-1 only**, so
adding countries broadens range but does **not** fix the admin-2 resolution gap.

---

## Tier 1 — Harmonization & comparability (highest leverage, cheapest)

1. **Obtain the BRINDA pooled dataset / definitions.** BRINDA already standardized
   ferritin + RBP/retinol + CRP/AGP inflammation adjustment across 17+ surveys
   (incl. Kenya 2007/2010, Cameroon, Côte d'Ivoire, Liberia, Malawi). Two wins:
   (a) re-anchor our four countries' biomarkers to ONE adjustment standard —
   directly targeting the level-offset problem; (b) any BRINDA country we add
   comes pre-harmonized. **Action:** request BRINDA data/recipes (brinda-nutrition.org).
2. **Extend the VMNIS prevalence extraction** (the NA-discard bug is now fixed in
   `build_kenya_benchmark.R`; the `Prevalenceofdeficiency` columns are usable) to
   build national held-out benchmarks for ALL countries, not just Kenya.
3. **Compile published national deficiency prevalences** from each survey's final
   report (GMNS/GMS/SLMS/MNS) to validate our national estimates and confirm the
   biology-vs-assay calls (Gambia high iron = real; Malawi folate-replete = real).

## Tier 2 — Add micronutrient-BIOMARKER surveys (no anemia-only)

| Priority | Survey | Yr | Biomarkers | Subnational | Access |
|---|---|---|---|---|---|
| 1 | **Nigeria NFCMS** | 2021 | full panel + iodine + CRP/AGP | zone/state | GHDx request |
| 2 | **Ethiopia ENMS (EPHI)** | 2015 | ferritin, sTfR, RBP, zinc, folate, B12, iodine | 9 regions | EPHI / BRINDA |
| 3 | **Kenya KNMS** | 2011 | ferritin, RBP, zinc, folate | 8 provinces | KNBS NADA (national summary already extracted) |
| 4 | **Cameroon NMS** | 2009 | ferritin, sTfR, retinol, zinc, folate, B12 | 3 macro-regions | Engle-Stone (UC Davis) / BRINDA |
| 5 | **Tanzania DHS (NUT5)** | 2010 | RBP, ferritin, iodine | national subsample | DHS Program (free) |
| — | Pakistan NNS 2018, Bangladesh NMS 2011 | | full panel | Pakistan = **district** (rare) | non-SSA comparators |

- **DHS biomarker surveys are a short list** (beyond anemia): Tanzania 2010, Uganda (retinol/RBP several rounds), Zimbabwe (RBP), plus the full-panel DHS-linked MN surveys Cambodia 2014 and Malawi 2015-16 (have). Confirm Rwanda RDHS 2019-20 module.
- **MICS is a dead end for serum biomarkers** — no MICS round measures ferritin/RBP/zinc/folate/B12. Use MICS for **covariates only**, not outcomes.

## Tier 3 — Additional GEE predictor domains (window to survey fieldwork dates)

Add a *parsimonious* set (effective n = 14–87 areas — don't dump all bands), chosen by pathway:
- **Food availability:** Accessibility-to-Cities/travel-time (`Oxford/MAP/accessibility_to_cities_2015_v1_0`), JRC surface water, Gridded Livestock of the World, GFSAD irrigation.
- **Agro-climatic stress:** TerraClimate (`IDAHO_EPSCOR/TERRACLIMATE` — PDSI, water deficit, VPD), SMAP soil moisture, MODIS ET, derived SPI/SPEI.
- **Infection/inflammation (bears on the ferritin offset):** PM2.5 (ACAG), Sentinel-5P aerosols, FireCCI burned area.
- **Structural:** SRTM-derived terrain (slope/TPI).

GEE housekeeping (from the completeness audit): **drop `gee_ndvi_2022`** (all-NA everywhere); treat the **ESA WorldCereal tile columns** as country-specific (structural, exclude from the shared transport set). Overall GEE coverage is solid — not a transport bottleneck.

## Admin-3 resolution
Only **Sierra Leone** (153 chiefdoms — already implemented) and **Malawi** (3,126 units — exists in GADM, *not* implemented, likely too fine to help). Gambia & Ghana have no admin-3 in GADM. Low priority.

---

### Suggested sequence
1. BRINDA data + VMNIS/report validation (Tier 1) — firms up what we can already claim.
2. Add Nigeria NFCMS 2021 + Ethiopia ENMS 2015 as the first new BIOMARKER countries (full panel; BRINDA-harmonizable).
3. Trial a small set of new GEE domains (inflammation + market access + drought), windowed to survey dates.
