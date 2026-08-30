# Pipeline review — data aggregation & merging (2026-05-12)

Critical review of the per-country data merge scripts and supporting
machinery. Each finding cites file + line, what's wrong, the consequence,
and a one-line fix.

Verified against the actual code; the original audit had two false positives
which are noted at the bottom for record.

---

## 🔴 Critical

### C-1. Malawi IHME admin-2 merge joins on the wrong key

- **File:** [src/malawi/2_GW_Malawi_data_merge.R:423](src/malawi/2_GW_Malawi_data_merge.R)
- **Bug:** The fuzzy-matched IHME admin-2 lookup is joined as
  `left_join(good_matches, by = c("Admin1"="Admin2"))`. The right-hand
  side (`good_matches`) has an `Admin2` column (district from the survey)
  mapped to `ihme_adm2_name` (IHME's district name). The join key on the
  left side should be `Admin2`, not `Admin1`.
- **Consequence:** Every Malawi row gets NA for all `ihme_*` admin-2
  columns, because regions (Admin1 = "Northern", "Central", "Southern")
  never match district names. The IHME admin-2 domain has been silently
  empty for Malawi in every pipeline run.
- **Fix:** Change to `by = "Admin2"` (or `by = c("Admin2"="Admin2")`).

### C-2. Malawi WFP food prices use FUTURE prices (temporal leak)

- **File:** [src/malawi/2_GW_Malawi_data_merge.R:151](src/malawi/2_GW_Malawi_data_merge.R)
- **Bug:** Malawi survey fieldwork was Dec 2015 – Feb 2016, but the WFP
  food-price filter is `filter(year == 2017)` — using prices 1–2 years
  *after* the survey was conducted.
- **Consequence:** `wfp_*` predictors encode future market conditions and
  could artificially inflate predictive performance through lookahead.
- **Fix:** Change to `filter(year %in% c(2015, 2016))` (average across
  fieldwork window) or `filter(year == 2015)`.

---

## 🟡 Important

### I-1. WFP year filters in Gambia and Sierra Leone use lag, not synchronous

- **Files:**
  - [src/Gambia/2_GW_Gambia_data_merge.R:196](src/Gambia/2_GW_Gambia_data_merge.R) — `filter(year == 2017)` for a 2018 fieldwork (1-year lag)
  - [src/SierraLeone/2_GW_SierraLeone_data_merge.R:142](src/SierraLeone/2_GW_SierraLeone_data_merge.R) — `filter(year == 2011)` for a Nov-Dec 2013 fieldwork (2-year lag)
- **Bug:** Past-year prices are valid predictors (not a leak in the C-2
  sense), but the temporal misalignment is asymmetric across countries
  and undocumented. Sierra Leone's 2-year lag is conspicuously long.
- **Consequence:** Cross-country comparability of `wfp_*` is weakened;
  Sierra Leone's `wfp_*` columns reflect conditions two years before
  survey, mostly noise relative to the deficiency-at-survey outcome.
- **Fix:** Standardise on `year == survey_year` (preferred) or document
  why a lag is used; if past-year is intentional, apply the same offset
  uniformly across countries.

### I-2. Country asymmetry: predictor domains present in only 1–2 countries

- **Files:** [R/config.R](R/config.R) per-country `domains` lists; merge
  scripts that conditionally execute joins
- **Bug:** Predictors that exist for some countries but not others
  produce wide cross-country NA blocks rather than a clean missing-data
  pattern. Specifically:
  - **LSMS** — Ghana only ([src/Ghana/2_GW_Ghana_data_merge.R:~600](src/Ghana/2_GW_Ghana_data_merge.R)); Gambia stub commented out at line 532; absent for SL / Malawi.
  - **MICS** — Gambia and Ghana only.
  - **FluNet** — Ghana only.
  - **DHS** — Ghana carries three rounds (`dhs2014_`, `dhs2016_`, `dhs2017_`), others carry one.
- **Consequence:** In cross-country transfer experiments
  ([R/transportability.R](R/transportability.R)) and the joint-country
  outcome targets (`outcome_data_GW_Gambia_SL_*`), columns missing for
  some countries become all-NA blocks. SuperLearner's NA-handling treats
  this as "no information" rather than "domain not available", which
  inflates model complexity without adding signal.
- **Fix:** Either (a) impute a documented sentinel + add a "domain
  availability" indicator column per country, or (b) drop these
  unilateral domains from cross-country joint outcomes and only use them
  in country-specific models.

### I-3. Ghana DHS still pinned to 2014 (3-year gap)

- **File:** [R/config.R:115-122](R/config.R) — `dhs_year = 2014L` with TODO comment
- **Bug:** Ghana survey is 2017, but admin-1 DHS predictors come from the
  2014 round. The pre-computed admin-2 BYM2 file
  ([src/DHS/DHS_admin2_aggregation.R](src/DHS/DHS_admin2_aggregation.R))
  is already configured for Ghana 2017 but hasn't been run.
- **Consequence:** Ghana DHS predictors lag by 3 years; large secular
  changes (especially in vaccination, WASH, anemia) are unrepresented.
- **Fix:** Run `src/DHS/DHS_admin2_aggregation.R` for Ghana 2017 and
  update `dhs_year = 2017L` in config.

### I-4. Three DHS rounds for Ghana — likely correlation between columns

- **File:** [R/config.R:207-209](R/config.R) — `DHS2014`, `DHS2016`, `DHS2017` all wired
- **Bug:** Ghana's merge brings in admin-1 estimates from three DHS
  rounds (2014, 2016, 2017). The same indicator across rounds is highly
  correlated.
- **Consequence:** Domain ablation can't distinguish round effects;
  variable importance is split arbitrarily across 3× collinear columns;
  prescreening may drop the "wrong" round.
- **Fix:** Pick one round (the closest to survey year) as primary; keep
  the other two as a separate `DHS_TREND` domain or drop them.

### I-5. IHME admin-1 merge lacks duplicate-key guard

- **File:** [src/IHME/build_ihme_admin1.R](src/IHME/build_ihme_admin1.R) — `merge_ihme_admin1()`
- **Bug:** The admin-1 IHME wide CSV is assumed to have unique
  `ihme_adm1_name` keys, but the build script collapses duplicates with
  `summarise(value = mean(value))` *only after pivot*. If the same
  (Admin1, variable) appears twice with different values pre-pivot,
  pivot_wider raises a warning but returns list-columns. The merge then
  fans out rows silently.
- **Consequence:** Possible row multiplication for countries where
  IHME admin-1 has duplicates.
- **Fix:** Add `stopifnot(!any(duplicated(ihme1$ihme_adm1_name)))` in
  `merge_ihme_admin1()` after reading the CSV; or use `dplyr::distinct()`.

### I-6. Ghana, Gambia: cluster-buffer GEE join risks duplicate-key fanout

- **Files:**
  - [src/Gambia/2_GW_Gambia_data_merge.R:53](src/Gambia/2_GW_Gambia_data_merge.R) — `left_join(d, gee, by = c("gw_MICS_Cluster_Number" = "gee_MICS_Cluster_Number.x"))`
  - [src/Ghana/2_GW_Ghana_data_merge.R](src/Ghana/2_GW_Ghana_data_merge.R) — similar pattern
- **Bug:** The GEE buffer CSV is assumed to have unique cluster IDs.
  If duplicates exist (e.g. two rows for the same cluster from different
  EE export runs), the join fans out.
- **Consequence:** Potential silent row multiplication; affects all
  downstream merges (admin-2 joins on a duplicated row produce identical
  duplicate rows in admin-2 results).
- **Fix:** Add `dplyr::distinct(gw_MICS_Cluster_Number, .keep_all = TRUE)`
  on the GEE side before joining.

### I-7. The `gw_` cluster_id is dropped by `gw_exclude_patterns` regex if it matches

- **File:** [R/config.R](R/config.R) — `gw_exclude_patterns` per country
- **Bug:** Patterns like `"Fol"` (for folate) match `"Folder"` or other
  unrelated `gw_` column names. Patterns are case-insensitive substring
  match. Not currently an issue because biomarker abbreviations are
  short, but with country-specific `gw_` columns there's nontrivial
  collision risk (e.g. `gw_Fol_eligible`).
- **Consequence:** Unintended `gw_` columns dropped from the
  sensitivity-with-`gw_` model.
- **Fix:** Anchor patterns with word boundaries: `"\\bFol\\b"` instead of
  `"Fol"`, or use exact column-name allowlists per country.

### I-8. Spatial join in Gambia/Ghana/SL drops clusters silently when GPS falls outside polygons

- **Files:**
  - [src/Gambia/2_GW_Gambia_data_merge.R:115](src/Gambia/2_GW_Gambia_data_merge.R)
  - [src/Ghana/2_GW_Ghana_data_merge.R:88](src/Ghana/2_GW_Ghana_data_merge.R)
  - [src/SierraLeone/2_GW_SierraLeone_data_merge.R:65](src/SierraLeone/2_GW_SierraLeone_data_merge.R)
- **Bug:** `st_join(d_sf, poly.adm, join = st_within)` returns NA for
  points outside any polygon (coastal clusters with slight GPS offset,
  border clusters). These NA-Admin clusters keep their rows but lose all
  admin-keyed predictors silently.
- **Consequence:** Rows with NA Admin1/Admin2 carry incomplete predictor
  matrices; not dropped, just degraded. May not be obvious in summary
  statistics.
- **Fix:** Log the count of NA-Admin rows post-join; for coastal NA
  cases, retry with `st_nearest_feature` to attach the closest polygon
  within a tolerance.

### I-9. Malawi's WFP path missing for survey-year prices

- **File:** [src/malawi/2_GW_Malawi_data_merge.R:151](src/malawi/2_GW_Malawi_data_merge.R)
- **Bug:** (See C-2.) Malawi's WFP price file may simply not have 2015–16
  observations — that may be why year=2017 was hard-coded. Need to check
  `data/food_price/wfp_food_prices_mwi.csv` coverage before fixing C-2.
- **Consequence:** If 2015–16 prices are missing from the WFP file,
  fixing the year filter to 2015–16 will produce all-NA `wfp_*` columns.
- **Fix:** Inspect the WFP file. If 2015–16 unavailable, document the
  fallback (use nearest available year before survey, not after).

---

## 🟢 Minor / cleanup

### M-1. `rm(list = ls())` at the top of every merge script

- **Files:** All four `src/<Country>/2_GW_*_data_merge.R`, line 3
- **Issue:** Defeats modularity; if scripts are ever wrapped as targets
  or sourced from another, this wipes the caller's environment.
- **Fix:** Remove `rm(list=ls())`; wrap each merge in a function with
  explicit input/output paths (or run as `Rscript` in a fresh process).

### M-2. Inconsistent CRS handling

- **Files:** All merge scripts
- **Issue:** Some country scripts call `st_transform(poly.adm, crs = 4326)`
  explicitly after `load_gadm_cached`; others don't. GADM rasters are
  already 4326, so the call is a no-op, but the inconsistency makes the
  scripts hard to compare.
- **Fix:** Standardise: always assume `crs = 4326` after `load_gadm_cached`,
  drop the `st_transform` call (or always keep it for safety).

### M-3. Multiple debugging variants of merge scripts in `src/Ghana/`

- **Files:** `src/Ghana/2_GW_Ghana_data_merge_admin2_V2.R`,
  `2_GW_Ghana_data_merge_cluster_def_prev.R`
- **Issue:** Several Ghana merge-script variants exist alongside the
  canonical `2_GW_Ghana_data_merge.R`. Unclear which is canonical and
  what the experimental variants are for.
- **Fix:** Delete or move to `src/_archive/`; ensure
  `2_GW_Ghana_data_merge.R` is the only one referenced from any pipeline
  or documentation.

### M-4. `check_categorical_levels` drops columns with > 100 levels silently

- **Files:** All four merge scripts (function definition + invocation)
- **Issue:** A 101-level factor (cluster ID parsed as char, free-text
  observation field, market name) is dropped without logging which
  columns were removed.
- **Fix:** `print(very_high_fact$high_level_variables)` and save the list
  to the metadata RDS so we know what was dropped.

### M-5. Hard-coded year suffixes on raster filenames

- **Files:** [src/Gambia/2_GW_Gambia_data_merge.R:217-261](src/Gambia/2_GW_Gambia_data_merge.R) (and 3 siblings) — `rasters <- c("Malaria__202206_*", ...)`
- **Issue:** The Malaria Atlas `rasters` vector hard-codes specific
  release dates (202206, 202406, 202508 — for some countries only). When
  a new MAP release is published, all four scripts must be edited.
- **Fix:** Glob the country folder for `Malaria__*.tif` and pick the
  latest version per family, similar to what `select_rasters_for_year`
  does for GEE.

### M-6. No `gee_a2_vars` / `ihme_adm1_vars` separate metadata entries

- **Files:** All four merge scripts — metadata save block at bottom
- **Issue:** The metadata RDS lists `gee_vars` and `ihme_vars` but
  doesn't separate the new admin-2 zonal (`gee_a2_*`) and admin-1 IHME
  (`ihme_adm1_*`) sub-domains. They roll up into the parent domain,
  which is fine for domain ablation but loses granularity for diagnostics.
- **Fix:** Optionally add `gee_a2_vars` and `ihme_adm1_vars` to the
  saved metadata list (a quick win — just extra `grep("^gee_a2_", ...)`
  calls).

### M-7. ESA_WorldCereal columns have long, unwieldy names

- **Files:** Output of `scripts/build_gee_admin2.R`
- **Issue:** Column names like
  `gee_a2_ESA_WorldCereal_2021_32121_TC_MAIZE_MAIN_ACTIVECROPLAND_20210327_20211119_classification`
  are >80 chars and contain raster-tile / date strings that are not
  meaningful predictors.
- **Fix:** Add a post-processing step in `extract_gee_admin2.R` that
  shortens common-stem ESA column names to something like
  `gee_a2_ESA_WorldCereal_maize_main`.

---

## Original audit items I **rejected** after verification

These were flagged in an initial automated audit but turn out to be
non-issues after inspecting the code:

- **"Malawi missing `st_join` to assign Admin1/Admin2"** — false. Malawi's
  input data (`data/IPD/Malawi/clean_malawi_mn_data.RDS`) already
  contains `Admin1`, `Admin2`, `mregion`, `admin2.name.full` from
  upstream cleaning ([src/malawi/malawi_IHME_clean.R](src/malawi/malawi_IHME_clean.R)).
  The script intentionally skips spatial join.
- **"`gee_a2_*` and `ihme_adm1_*` invisible to domain ablation"** —
  false. [R/data_prep.R:370,377](R/data_prep.R) uses
  `startsWith(all_cols, dom$prefix)`, which catches both
  `gee_*`/`gee_a2_*` and `ihme_*`/`ihme_adm1_*` under their respective
  parent domains. They appear together in domain ablation, which is the
  intended behavior. (If you want to separately ablate the admin-2
  zonal vs cluster-buffer GEE, you'd need to add a new domain entry
  with prefix `gee_a2_` and pre-pend it to the prefix loop so it's
  matched before `gee_` — but that's a feature request, not a bug.)
- **"WFP year mismatch in Gambia / SL is a temporal leak"** — partial
  reject. Gambia (year==2017 for 2018 survey) and SL (year==2011 for
  2013 survey) use *past-year* prices, which is valid as a predictor.
  The misalignment is still worth standardising (see I-1) but it's not
  a leak. Only Malawi (year==2017 for 2015–16 survey) is a real
  future-prices leak (C-2).

---

## Priorities

If you have a small fix budget:

1. **C-1** (Malawi IHME join key) — 1-line fix, recovers ~120 IHME admin-2 predictors for Malawi.
2. **C-2** (Malawi WFP future prices) — 1-line fix, removes a real temporal leak.
3. **I-3** (Ghana DHS upgrade to 2017) — half-day run of `src/DHS/DHS_admin2_aggregation.R`, then config change.
4. **I-2** (cross-country NA blocks) — design decision; document then implement.
5. **I-1, I-4** (DHS / WFP alignment) — small fixes, document the policy.

Everything else is cleanup; can wait for a refactor sweep.
