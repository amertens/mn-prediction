# Adding a new country — checklist

A fill-in-the-blanks guide for incorporating a new country into the LOCO panel,
organized by how much manual work each dataset needs. **Tanzania (TDHS 2009-10,
vitamin A only)** is the worked example throughout; its scripts under
`src/Tanzania/` are the copy-from templates, and
`src/Tanzania/README_TANZANIA_TODO.md` is the per-domain status log to mirror.

The guiding fact: the covariate domains that carry the *transportable* signal
(earth-observation, soil, IHME) now auto-extract from global sources, so the real
per-country work is the **outcome-side survey cleaning + DHS proxies + config**.
Work top-down — Tier 3 (outcome + config) first, since nothing runs without it,
then let Tier 1/2 extractions fill in behind the `file.exists()`-guarded merge.

Replace `<Country>` / `<ISO3>` / `<year>` throughout. Commands run from the repo
root with `"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla`.

---

## Tier 3 — Manual, country-specific (do first; nothing runs without these)

### 3.1 Config block — `R/config.R`
Add an entry to `get_country_configs()` (copy the `Tanzania = list(...)` block at
the end). Set the identity and design fields, define only the outcomes with
usable biomarker data, and list the domain prefixes.

- [ ] `country`, `gadm_code` (`<ISO3>`), `survey_year`, `dhs_year`, `data_path`
      (`data/IPD/<Country> <year>/<Country>_merged_dataset.rds`), `raster_dir`.
- [ ] Design columns: `cluster_id`, `admin1_col`, `admin2_col`, `psu_col`,
      `weight_col`, `strata_col`, `child_flag`.
- [ ] `outcomes = list(...)` — one block per usable biomarker: `tag`, `label`,
      `population`, `child_flag_val`, `continuous`, `binary`, `cutoff`,
      `cutoff_dir`, `cutoff_scale`. **Exclude biomarkers that aren't
      cross-country comparable** (e.g. Tanzania omits sTfR-based iron — not
      ferritin-comparable — and single-country iodine).
- [ ] `domains = list(...)` — set each source's prefix; note the DHS prefix is
      **year-stamped** (`dhs2010_` for Tanzania, `dhs2019_` for Malawi).
- [ ] `gw_exclude_patterns` — regex list so raw biomarker columns never leak into
      the proxy predictors.

> The pipeline loops over `get_country_configs()`, so once the block exists and
> `data_path` resolves, `_targets` picks the country up automatically.

### 3.2 Outcome clean script — `src/<Country>/1_GW_<Country>_data_clean.R`
**The core effort; every country differs.** Copy
`src/Tanzania/1_GW_Tanzania_data_clean.R` and adapt. Produces
`<Country>_GMS_cleaned.rds` (one row per individual with the biomarker outcome,
GPS/cluster, weights, child flag).

- [ ] Locate the biomarker microdata (DHS OB/biomarker recode, MICS, or national
      survey). **Verify variable names against a live `makeVlist()` read**, not
      just the dictionary.
- [ ] Barcode → PR crosswalk if needed. **Known gotcha:** R 4.4.2 on this machine
      segfaults on some `haven::read_dta`/recode calls — do the recode in Python
      (`0_extract_pr_crosswalk.py`, pandas) as Tanzania did.
- [ ] **Unit conversion** — confirm the biomarker scale. Tanzania RBP is **mg/L,
      not µmol/L** (÷21.2 before the <0.70 VAD cutoff). Check against the survey's
      biomarker documentation.
- [ ] **Inflammation adjustment** — apply BRINDA where markers exist. Note what's
      available: Tanzania is **CRP-only** (no AGP), so it uses the OB `rbpadcrp`
      and can't do the full two-marker correction. Check
      `apply_brinda_vita_binary()` in `R/data_prep.R`.
- [ ] **Weighting** — use the micronutrient-subsample weight if present; if not
      (Tanzania, Gambia), fall back to the household weight and **document the
      subsample-weight caveat** for the manuscript.
- [ ] GPS → cluster join; emit the design columns named in the config block.

### 3.3 Merge script — `src/<Country>/2_GW_<Country>_data_merge.R`
Copy `src/Tanzania/2_GW_Tanzania_data_merge.R`. GADM admin-2 join + every
covariate-domain join, **each guarded with `file.exists()`** so it runs to
completion as extractions come online and can be re-run as they land.

- [ ] Confirm the GADM admin-2 join matches most clusters (`gadm41_<ISO3>_2`
      downloads on first run).
- [ ] Produces `<Country>_merged_dataset.rds` + `metadata/<XX>_variable_categories.rds`.

### 3.4 DHS proxy aggregates — `dhs<year>_` (account-gated)
- [ ] Obtain DHS recode files (IR/KR/PR/HR) — **DHS account/approval required**
      (`rdhs`).
- [ ] Run `src/DHS/DHS_admin2_aggregation.R` wired for `<Country> <year>`.
      Analog: `src/SierraLeone/SL_DHS_aggregation_admin2.R`.

### 3.5 Optional registration-gated surveys
- [ ] **LSMS** (`lsms_`, admin-1, coarse) — World Bank Microdata Library (free
      account + terms). Then `src/<Country>/<XX>_LSMS_clean.R`. Only some
      countries use it.
- [ ] **MICS** (`mics_`) — UNICEF registration; skip if no round in scope
      (Tanzania: N/A).

---

## Tier 1 — Auto-extraction (add the country, point at its boundaries)

These pull on-demand from global sources. One-time machine setup for GEE is
already done (see below); per country it's a small runner or an ISO code.

### 1.1 GEE earth-observation (`gee_`) — the transportable signal

Two separate products, and they are **not interchangeable**. Read this before
running anything.

**(a) Legacy-parity admin-2 covariates — REQUIRED for any pooled/LOCO analysis.**

- [ ] Run `Rscript scripts/build_gee_legacy_parity.R <Country>` →
      `data/GEE/<ISO3>_legacy_parity_admin2_gee.csv` (135 columns).
      `extract_gee_admin2()` / `extract_area_covariates()` pick this up
      automatically for any country with no `data/<Country>_GEE_rasters/`.
- [ ] Sanity-check the shared vocabulary before and after:
      `Rscript scripts/check_gee_legacy_parity.R` prints, per country, how many
      variables it costs the intersection.

Why this exists: the original four countries' admin-2 covariates are zonal means
over local `.tif` exports, and `.append_gee_zonal_cols()` names each column after
the **raster filename** (`gee_ndvi_2013`, `gee_soilzinc_mean_0_20`, ...). The
Earth Engine API names them after **EE bands** (`gee_a2_NDVI`,
`gee_a2_SoilZinc_mean_0_20`). Those two vocabularies overlap **0%**, and every
pooled / LOCO builder takes a strict name **intersection** — so a country added
via the API path contributes nothing *and deletes the whole GEE block for every
other country*, with no error. That is what
`scripts/build_gee_legacy_parity.R` fixes.

- [ ] **If you change the extractor, re-validate it.** Run it against a country
      that has both representations and compare:
      `Rscript scripts/build_gee_legacy_parity.R SierraLeone --validate`
      (14 districts, ~5 min). Every column must pass **r ≥ 0.9 against the
      raster-derived values AND a mean within 10%**. Rank agreement alone is not
      enough — a constant scale offset on one country is read by the pooled model
      as a real country difference. Results land in
      `results/sensitivity/gee_legacy_parity_validation_<country>.csv`.
      This check is not optional: the first version of the extractor passed only
      22/148 columns and the failures were silent scale errors (a year-over-year
      band averaged into the base value, a missing iSDAsoil back-transform, and
      mask-vs-zero semantics).
- [ ] **Do not raise `--scale-floor` above the 250 m default** to speed a large
      country up. Measured on Sierra Leone: at 1000 m every depth-*mean* column
      still passes, but the soil `stdev` bands and `gee_wsf_2015` come out 10-20%
      high (135/135 pass at 250 m; 117/135 at 1000 m) — variance and mask
      quantities do not survive coarse aggregation. Keep `--batch 40`; larger
      batches exceed EE's 10 MB request limit on detailed admin-2 geometry.
      Budget a few hours for a large country instead; the per-(family, year)
      cache means an interrupted run resumes where it stopped.

**(b) Full EE layer set (`gee_a2_`, ~162 cols) — within-country modelling only.**

- [ ] Copy `src/Tanzania/4_TZ_GEE_extract.R` → `src/<Country>/4_<XX>_GEE_extract.R`;
      point it at the country's GADM polygons + cluster buffers.
- [ ] Run it → writes `data/GEE/<Country>_<year>_admin2_gee.csv` +
      `data/GEE/<XX><year>_buffers_<date>.csv` (merge picks them up).
- [ ] **First-run validation:** watch per-layer `ok:` logs, note any deprecated
      asset substitutions, and spot-check `gee_a2_SoilZinc`.

These columns are richer than the parity set but are country-specific by name,
so they are dropped by every cross-country intersection. Useful for the
within-country model; irrelevant to LOCO.

One-time GEE setup (per `src/GEE/README_GEE_API_SETUP.md`, already configured on
this machine): EE account `amertens@berkeley.edu`, Cloud project
`mn-prediction-420517`, the `rgee` conda env. Only re-run `ee_Authenticate()` if
the token is stale, then `src/GEE/00_ee_connect_test.R` must print
"Connectivity test PASSED".

### 1.2 SoilGrids (`soil_`, 7 props) — *most* transportable covariate; run early
- [ ] Add `<Country>` (code `<ISO3>`, matching the merge's GADM) to the
      `countries` list in `scripts/build_soilgrids_admin2.R`.
- [ ] Run: `Rscript scripts/build_soilgrids_admin2.R <Country>` →
      `data/SoilGrids/<Country>_soilgrids_admin2.csv`. Fetched on-demand from the
      global ISRIC VRT; no download.

### 1.3 IHME / LBD (`ihme_`, ~200 admin-2 + admin-1 cols)
- [ ] Copy `src/Tanzania/3_TZ_IHME_clean.R`; it calls the generic
      `src/IHME/build_ihme_admin2.R` + `build_ihme_admin1.R` over the global LBD
      CSVs already in `data/IHME/<family>/`. Run it → `..._merged_IHME_data.csv`.
- [ ] Check the printed **PARITY** line vs Ghana for column alignment.
- [ ] **Known gap:** no admin-2 CSVs for **malaria** or **education** — those
      `ihme_` admin-2 cols will be absent (malaria is covered by `MAP_`,
      education by DHS proxies). To close it, pull the LMIC malaria/education
      *ADMIN_2* CSVs from GHDx.

### 1.4 GADM boundaries + global statics — no action
`gadm41_<ISO3>` auto-downloads on the first merge run. Global static layers
(SPAM crop, Köppen-Geiger, AEZ, FEWS NDVI anomaly, GLW4 livestock) are in-repo
and flow through the local-raster path.

---

## Tier 2 — Scripted, needs a per-country public fetch (semi-auto)

### 2.1 WFP food prices (`wfp_`)
- [ ] Download the country CSV from HDX → `data/food_price/wfp_food_prices_<iso3>.csv`.
      Merge block parses it (nearest-market join).

### 2.2 Malaria Atlas (`MAP_`)
- [ ] Copy an analog (`src/SierraLeone/SL_MAP_download.R`, or
      `src/Tanzania/TZ_MAP_download.R`) and run against the global MAP rasters.

### 2.3 Food security (`fsec_`)
- [ ] Add the country to the config maps in `R/food_security.R` (HFID is global;
      Cadre Harmonisé is West/Central Africa only).
- [ ] **Gotcha:** the HFID admin-2 names must crosswalk to GADM or the columns
      come out **all-NA** (Tanzania hit this — a name-crosswalk fix is pending).

---

## Validate

- [ ] Add `<Country>` to the `countries` list in
      `scripts/validate_merged_datasets.R`; confirm per-prefix column counts and
      NA rates match the other countries.
- [ ] Confirm `_targets` includes the country in the LOCO targets (it loops over
      `get_country_configs()`).

## Two caveats to carry every time

1. **LOCO column-name parity.** Covered in detail in 1.1(a): run
   `scripts/build_gee_legacy_parity.R`, and re-validate with `--validate` if you
   touch the extractor. Two known residual gaps for an API-only country:
   - 14 of the legacy 149 admin-2 variables cannot be reproduced on the same
     scale and are deliberately excluded (`GEE_PARITY_EXCLUDED` in
     `src/GEE/extract_gee_legacy_parity.R`): the 5 ESA-WorldCereal summaries,
     `gee_accessibility_2019` (asset corrupted server-side), and 7 of the 9
     `gee_soiltotalcarbon_*` columns. Including the country removes those from
     the shared vocabulary for everyone.
   - The **cluster-buffer** `gee_` columns (186 variables — monthly
     `ndvi`/`temp`/`trmm` at 10/25/50 km plus `fpp_`/`tpp_`/`ghm`/`grassland`/
     `smod`/`built_surface`) came from an analyst's EE Code Editor export whose
     source is not in the repo; `fpp_`/`tpp_` are not identifiable from the
     column names alone. A new country cannot reproduce them, so including it in
     the **individual-level** LOCO costs those 186 covariates. The area-level
     (Admin-2) LOCO is unaffected. Recovering the original export script would
     close this.
2. **Temporal spread.** Older survey years (Tanzania 2010 vs the 2013–2018 panel)
   are handled by year-stamped proxies, but widen the temporal spread — note it.

## The one account-gated dependency (besides the outcome data)

**DHS access** for the `dhs_` proxies. Everything else in Tier 1–2 is public or
already in-repo.
