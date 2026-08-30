# Pipeline review — cleaning, documentation, feature engineering (2026-05-13)

Second-pass audit, focused on the **upstream cleaning pipeline**, **dormant
helpers**, **documentation gaps**, and **feature-engineering opportunities**.
Builds on [pipeline_review_2026-05-12.md](pipeline_review_2026-05-12.md) —
findings already documented there are not repeated.

All findings verified against current `main`.

---

## 🔴 Critical

### N-1. FSEC (Cadre Harmonisé / HFID) defined but never invoked

- **Files:**
  - [R/food_security.R:260](R/food_security.R) — `merge_food_security(merged_data, cc, hfid_path, ch_path)` is fully implemented.
  - Verified: zero callers across `src/`, `R/`, `scripts/` (only `R/food_security.R` references the function name).
  - Existing data: `data/CadreHarmonise/cadre_harmonise_caf_ipc_dec25.xlsx`, `data/HFID/hfid_hv1.csv` are present.
- **Consequence:** The `fsec_` domain in [R/config.R](R/config.R) is configured for all four countries and the manuscript appendix ([docs/data_sources_appendix.qmd](docs/data_sources_appendix.qmd)) explicitly documents Cadre Harmonisé + HFID as wired sources. **They are not — the merged datasets carry zero `fsec_*` columns.** This contradicts the appendix.
- **Fix:** Add a `merge_food_security(df, cc, here("data/HFID/hfid_hv1.csv"), here("data/CadreHarmonise/cadre_harmonise_caf_ipc_dec25.xlsx"))` call near the end of each Stage 2 merge script (after the IHME admin-1 join), and append `fsec_vars` to the metadata save block.

### N-2. Ghana ferritin not Brinda-adjusted — cross-country outcome-definition inconsistency

- **File:** [R/config.R](R/config.R) — `continuous` field per country/outcome:
  - **Gambia** (children + women iron): `gw_LogFerAdj` — Brinda-adjusted log ferritin
  - **Ghana** (children): `gw_cLogFerr` — **no "Adj" suffix; un-adjusted log ferritin**
  - **Ghana** (women): `gw_wLogFerr` — un-adjusted
  - **Sierra Leone**: `gw_cFerrAdj` / `gw_wFerAdjBR1` — Brinda-adjusted
  - **Malawi**: raw `sf_reg` (different schema again)
- **Consequence:** The Brinda inflammation adjustment (corrects ferritin for elevated CRP/AGP) materially shifts deficiency prevalence in malaria-endemic populations. Without it, Ghana's ferritin-based iron-deficiency outcome is systematically biased compared to Gambia and SL. Cross-country transfer experiments (`R/transportability.R`) will absorb this bias as country-specific noise.
- **Fix:** Verify whether Ghana has a Brinda-adjusted column (likely `gw_cFerAdjBR1` or `gw_cLogFerAdj` — check the survey codebook). If yes, switch config. If no, derive it in [src/Ghana/1_GW_Ghana_data_clean.R](src/Ghana/1_GW_Ghana_data_clean.R) from raw ferritin + CRP + AGP using `BRINDA::apply_adjustments()`. **Until corrected, document the inconsistency in the appendix as a known limitation**, and report a sensitivity model that uses BRINDA-adjusted ferritin for all four countries.

---

## 🟡 Important

### N-3. `extract_soilgrids()` defined but never called — entire SoilGrids domain absent

- **File:** [R/external_data.R:674](R/external_data.R) — fully implemented helper for SoilGrids pH, soil organic carbon, nitrogen, clay, sand, silt, CEC at Admin-2.
- Only references: the function itself, plus `scripts/download_external_predictors.R` (a one-shot orchestrator that *also* isn't called from any merge script).
- **Consequence:** SoilGrids macronutrient / pH layers are missing from the predictor matrix. The Soil domain currently consists only of iSDAsoil micronutrients (Fe/Zn/Ca/Mg/N/P/K/S/Al/TC/CEC × 4 bands = 44 cols of `gee_a2_Soil*` from the GEE folder). pH and SOC — which mediate plant micronutrient bioavailability — are absent.
- **Fix:** Pre-compute SoilGrids at Admin-2 via a standalone build script (analogous to `scripts/build_gee_admin2.R`), write to `data/external/<Country>_soilgrids_admin2.csv`, and left-join in each merge script as `soil_*` columns.

### N-4. `R/external_data.R::extract_all_external()` orchestrator is fully dormant

- **File:** [R/external_data.R:~2017](R/external_data.R) — wraps `extract_soilgrids()`, `extract_chirps()`, `extract_worldpop()`, `extract_ntl()` (Nighttime Lights), `extract_map_rasters()`, etc.
- Verified zero callers from any merge script or pipeline target.
- **Consequence:** A significant body of admin-2 extraction code is unused. WorldPop, CHIRPS rainfall, NTL — these may or may not duplicate things in the GEE buffer / admin-2 outputs, but the appendix doesn't clarify which is canonical.
- **Fix:** Audit which extracts are duplicated by the GEE admin-2 pipeline (CHIRPS ≈ TRMM, WorldPop ≈ GHSPOP, NTL ≈ CCNL), then either (a) wire the non-duplicates into the merge scripts, or (b) delete the dormant code from `R/external_data.R` and update the appendix to remove it.

### N-5. Malawi Stage 1 cleaning script uses non-standard naming and produces a non-standard output path

- **Files:**
  - Gambia: [src/Gambia/1_GW_Gambia_data_clean.R](src/Gambia/1_GW_Gambia_data_clean.R) → `data/IPD/Gambia/Gambia_GMS_cleaned.rds`
  - Ghana: [src/Ghana/1_GW_Ghana_data_clean.R](src/Ghana/1_GW_Ghana_data_clean.R) → `data/IPD/Ghana/Ghana_GMS_cleaned.rds`
  - SL: [src/SierraLeone/1_GW_SierraLeone_data_clean.R](src/SierraLeone/1_GW_SierraLeone_data_clean.R) → `.../SierraLeone_GMS_cleaned.rds`
  - **Malawi:** [src/malawi/1_DHS_mn_data.R](src/malawi/1_DHS_mn_data.R) → `data/IPD/Malawi/clean_malawi_mn_data.RDS`
- **Consequence:** Malawi's upstream script has a misleading name (says "DHS" — it processes the MNS biomarker survey, not DHS) and a different output filename pattern (`.RDS` capitalised, prefix `clean_`). A new collaborator running the pipeline by convention would not find Malawi's cleaning script. Not in the methodology appendix.
- **Fix:** Rename to `src/malawi/1_GW_Malawi_data_clean.R`, save as `Malawi_GMS_cleaned.rds`, and update the Stage 2 read path. Document the per-country Stage 1 inputs/outputs in a new section of the appendix.

### N-6. `read.csv()` used throughout merge scripts — type-coercion fragility

- **Files:** All four `2_GW_*_data_merge.R` scripts call `read.csv()` ~5-10 times each (e.g., [src/Gambia/2_GW_Gambia_data_merge.R:35,136,156,325](src/Gambia/2_GW_Gambia_data_merge.R)). Stage 1 scripts also.
- **Consequence:** Base `read.csv()` auto-infers types loosely. Joins fail silently if a numeric ID column is parsed as character on one side and integer on the other (left_join silently produces NA). The WFP food-price CSV in particular has a 2-row header (variable label + units) that base read.csv handles ungracefully.
- **Fix:** Replace with `readr::read_csv(file, show_col_types = FALSE)` and explicit `col_types = cols(cluster_id = col_integer(), ...)` for join-key columns. The WFP files specifically need `skip = 1` or `col_names = FALSE` handling.

### N-7. GEE buffer CSVs hard-coded with export date in filename

- **Files:** All four merge scripts —
  - `data/GEE/gambia2018_buffers_01.08.2026.csv`
  - `data/GEE/ghana2017_buffers_01.08.2026.csv`
  - `data/GEE/SL2012_buffers_01.08.2026.csv`
  - `data/GEE/malawi2015_buffers_01.08.2026.csv`
- **Consequence:** Next time the GEE export is run, the date changes (e.g., `06.15.2026`) and all four merge scripts break with cryptic "no such file" errors. Reproducibility risk.
- **Fix:** Glob the country folder for `<country><survey_year>_buffers_*.csv`, pick the latest by `file.info()$mtime`. Or just rename the files to drop the date (and keep date in metadata).

### N-8. The `0-functions.R` / `DHS_functions.R` / `DHS_variable_recode.R` helpers are sourced everywhere but documented nowhere

- **Files:**
  - [src/0-functions.R](src/0-functions.R) — 704 lines. Defines `fuzzy_match_admin2()`, `clean_DHS()`, `getDHSdata()`, `makeVlist()`, `evaluate_sl_performance()`, etc.
  - [src/0-SL-setup.R](src/0-SL-setup.R) — 463 lines, sl3 setup.
  - [src/DHS/DHS_functions.R](src/DHS/DHS_functions.R), [src/DHS/DHS_variable_recode.R](src/DHS/DHS_variable_recode.R) — sourced by every Stage 1 script.
- **Consequence:** New contributors have no map of the helper landscape. Some functions probably overlap with `R/data_prep.R::load_dhs_admin1()` etc. — duplicate logic with subtle divergences.
- **Fix:** Add a "Pipeline architecture" section to [docs/data_sources_appendix.qmd](docs/data_sources_appendix.qmd) listing each helper file, its exported functions, and intended scope. Audit for duplication with `R/*.R` and consolidate.

---

## 🟢 Minor / cleanup

### N-9. Stage 1 date parsing falls back silently

- **File:** [src/Ghana/1_GW_Ghana_data_clean.R:~114](src/Ghana/1_GW_Ghana_data_clean.R) — `parse_date_time(..., orders = c("a b d H:M:S z Y", "a b d H:M:S Y", "Y-m-d H:M:S"))`
- **Issue:** If a date is in none of the three formats, the result is NA without warning.
- **Fix:** Assert `mean(!is.na(parsed_dates)) > 0.95` after parsing.

### N-10. WFP nearest-market join has no distance cutoff

- **Files:** All merge scripts — every cluster gets its single nearest market's prices, regardless of distance.
- **Issue:** A cluster 250 km from any WFP market gets prices that are noise. The `nearest_market_distance_km` column is computed but never used as a quality filter.
- **Fix:** Set `wfp_*` to NA for clusters where `nearest_market_distance_km > 100` (or report this distance distribution and choose an empirical threshold).

### N-11. The Malawi `1_DHS_mn_data.R` script loads non-existent data paths and may be stale

- **File:** [src/malawi/1_DHS_mn_data.R:~10](src/malawi/1_DHS_mn_data.R) — `read_dta(here("data/Malawi/MW_PSC.DTA"))` etc. Path uses `data/Malawi/` (capital M, no IPD).
- **Issue:** Need to verify the script still runs end-to-end. The `data/Malawi/` path is distinct from `data/IPD/Malawi/` used by Stage 2.
- **Fix:** Run the script; document the exact .dta files it needs in the appendix.

---

## 💡 High-value feature engineering

The pipeline currently treats each predictor as raw. Below are derived features that are computable with **no new data** — all from columns already in the merged datasets.

### F-1. Seasonality indices from monthly GEE buffer columns

The GEE buffer CSVs include 12 monthly columns each for `trmm_*` (rainfall), `ndvi_*`, `temp_*` at 10/25/50 km buffers — ~108 monthly columns per cluster. Most modelers won't surface seasonality from raw monthly columns; SuperLearner's tree learners will pick a few months as splits, missing the seasonal *shape*.

For each (family × scale) — i.e., 9 combinations — derive 4 summary features:

| Feature | Formula | Why it matters |
|---|---|---|
| `<family>_<scale>_annual_mean` | mean across 12 months | year-level baseline (also new direct measure if not already present) |
| `<family>_<scale>_seasonal_cv` | sd / mean across months | seasonality intensity (drought vulnerability) |
| `<family>_<scale>_seasonal_range` | max − min across months | absolute spread |
| `<family>_<scale>_peak_month` | which.max | timing of peak (rainfall, vegetation) |
| `<family>_<scale>_concentration` | Herfindahl-Hirschman index on normalised values | how concentrated in one season |

That's ~36 new features replacing 108 raw monthly columns — much higher signal density.

### F-2. Year-over-year change from annual GEE rasters

The country GEE folders have annual rasters for 2010–2023 for several families (NDVI, TRMM, PopDensity). Extracting:
- `gee_<family>_yoy_change_t-1` = value at survey_year minus value at survey_year-1
- `gee_<family>_3yr_volatility` = sd over (survey_year-3) to survey_year

These capture *recent change* — a drought 2 years before survey, or rapid urbanisation. ~10 new features.

### F-3. Relative-position features (cluster prevalence vs. district / region / national)

For each outcome and each cluster:
- `<outcome>_above_district_avg` = cluster prevalence − district pop-weighted average
- `<outcome>_above_admin1_avg` = cluster − admin1 average
- `<outcome>_above_national_avg` = cluster − national average

These capture spatial autocorrelation (clusters in low-prevalence districts that are outliers reveal local context the model can use). Worth ~3-5 features per outcome.

### F-4. Distance-to-X features from existing rasters

Already on disk:
- **Distance to nearest large city** — derivable from MAP `Accessibility` rasters (`Accessibility_<Country>_2019.tif`) at admin-2.
- **Distance to coast / inland water** — derivable from existing land-cover rasters (mask water pixels, compute `terra::distance`).
- **Distance to nearest WFP market** — `nearest_market_distance_km` already exists; needs to be retained, not dropped by `check_categorical_levels` (since it's numeric, this should be fine — but verify).

### F-5. Interaction terms with strong biological motivation

A handful of pairwise interactions worth pre-computing (saves the model from having to discover them through tree-splits on millions of rows):

| Interaction | Reasoning |
|---|---|
| `malaria_prevalence × hbs_allele_freq` | HbS protects from malaria but confounds ferritin → joint effect on iron-deficiency assessment |
| `WASH_index × diarrhea_prevalence` | WASH access modulates diarrhea, which depletes nutrients |
| `wealth_quintile × food_diversity` | Wealth and dietary diversity are correlated but the *interaction* captures unequal access |
| `rural_indicator × accessibility_travel_time` | Rural + far from healthcare is a different exposure than rural alone |
| `cropland_fraction × precipitation_anomaly` | Cropland is vulnerable to drought; non-cropland less so |

Pre-computing these adds ~5 columns; the model would discover them eventually but with much weaker signal.

### F-6. Spatial-lag features (Tobler's first law)

For each Admin-2 polygon, compute:
- `<outcome>_neighbor_avg` = mean of outcome across queen-contiguous neighboring polygons (where survey data exists)
- `<predictor>_spatial_lag` = same for key continuous predictors (NDVI, malaria prevalence)

This is the standard Moran's-I style construction. Already partially implemented in `R/gwr_analysis.R` for residuals — generalise to a feature builder.

### F-7. Dietary diversity proxies from LSMS (Ghana)

[src/Ghana/2_GW_Ghana_data_merge.R](src/Ghana/2_GW_Ghana_data_merge.R) joins LSMS but only carries the raw 275 columns. Derive composite indices:
- **HDDS** (Household Dietary Diversity Score) = sum of food groups consumed in past 24h (0-12)
- **MDD-W** (Minimum Dietary Diversity for Women) = binary, ≥5 of 10 food groups
- **Animal-source food consumption frequency** — a single column distilled from ~30 individual food columns

These are validated nutrition indicators that map directly to micronutrient status. Currently each food column enters the predictor matrix in isolation; learners would have to discover the right combination from 275 columns.

---

## Priorities (after the 2026-05-12 fixes are in)

If you have a small budget:

1. **N-1** (FSEC wiring) — high impact, ~1 hour to add `merge_food_security()` calls + verify.
2. **N-2** (Ghana ferritin Brinda) — outcome correctness; first check the codebook for an existing adjusted column; otherwise derive it.
3. **N-3** (SoilGrids wiring) — half-day; high biological relevance (pH, organic carbon mediate iron/zinc bioavailability).
4. **F-1** (seasonality indices) — pure feature engineering, no new data, ~3 hours; likely largest model-AUC improvement of anything in this review.
5. **F-7** (Ghana LSMS dietary diversity) — high biological relevance for the iron and B12 outcomes specifically.

Everything else (N-4 through N-11, F-2 through F-6) is incremental.
