# Tanzania (TDHS 2009-10) — incorporation status & remaining work

**Scope:** vitamin A only (RBP-based VAD, women + children), as an East-African
addition to the LOCO panel. Iron (sTfR, not ferritin-comparable) and iodine
(population-median metric, single country, no LOCO) are intentionally excluded —
see the rationale in the Tanzania block of `R/config.R`.

**Current state:** the *outcome half* is scaffolded. The DHS files under
`data/IPD/Tanzania 2010/` are sufficient to build the RBP outcome, link it to
admin-2, and build the DHS proxies — but **not** the earth-observation / modelled
covariate domains, which carry the transportable signal and must be extracted
separately. Tanzania is **not runnable end-to-end** until the items below exist.

---

## Done (scaffold)

- [x] `R/config.R` — `Tanzania` config block (vit A outcomes, dhs2010_ prefix).
- [x] `src/Tanzania/1_GW_Tanzania_data_clean.R` — OB bar-code → PR crosswalk →
      RBP outcome → GPS join → `Tanzania_GMS_cleaned.rds`.
      **Verify the dictionary-derived assumptions in the script header against a
      live `makeVlist()` read before trusting the output.**
- [x] `src/Tanzania/2_GW_Tanzania_data_merge.R` — GADM admin-2 join + all
      covariate-domain joins, each **guarded with `file.exists()`** so it runs
      to completion as extractions come online. Produces
      `Tanzania_merged_dataset.rds` + `metadata/TZ_variable_categories.rds`.

## ✅ MERGED DATASET BUILT (2026-07-01)

`data/IPD/Tanzania 2010/Tanzania_merged_dataset.rds` — 16,357 rows x 423 cols.
Domains in: outcome (RBP/VAD), GEE (162, ~3% NA), soil (7), IHME admin-2 (101,
157/165 districts matched) + admin-1 (117), WFP (3). Verified via pyreadr.
Ran IHME (`3_TZ_IHME_clean.R`, 89% col parity vs Ghana) + merge
(`2_GW_Tanzania_data_merge.R`, RUN_MAP=false). `_targets` picks Tanzania up
automatically (config wired; the `aggregation_level_national` target skips it
until now — it's included on next pipeline definition since data now exists).

Remaining/optional: **fsec** columns came out ~all-NA (Tanzania HFID Admin2
names don't join to GADM — needs a name-crosswalk fix in `R/food_security.R`);
**MAP** (malaria) toggled off (needs global rasters pointed at the TZ extraction);
**dhs2010_** DHS proxy aggregates (run `src/DHS/DHS_admin2_aggregation.R` for TZ);
**LSMS** (NPS download). None block modeling.

NOTE: this machine's R (4.4.2) has an intermittent native SEGFAULT on
`haven::read_dta` and some `readRDS` calls (non-deterministic). Reads succeed on
retry; pyreadr (Python) reads reliably. Background `Rscript` runs completed fine.

## 2026-08-18 — cross-country parity work (the real blocker)

The merged dataset had been built since July, but Tanzania was still absent from
**every** results table, and simply running the pipeline would have made things
worse, not better. Two structural problems, both now addressed:

### 1. Tanzania's GEE covariates shared ZERO column names with the other four

`src/Tanzania/4_TZ_GEE_extract.R` names columns after Earth Engine **bands**
(`gee_a2_NDVI`, `gee_a2_SoilZinc_mean_0_20`). The four legacy countries' admin-2
covariates are named after **raster filenames** by `.append_gee_zonal_cols()`
(`gee_ndvi_2013`, `gee_soilzinc_mean_0_20`). Overlap: **0 of 149**.

Because `build_pooled_dataset()` and `build_area_loco_dataset()` take a strict
name **intersection**, adding Tanzania would have silently:
- dropped the area-level LOCO GEE set from **149 to 0** covariates, and
- dropped the individual-level pooled proxy set from **332 to 139**,

for *every* country, with no error. That is why Tanzania must never be switched
on by config alone.

**Fix:** `scripts/build_gee_legacy_parity.R` +
`src/GEE/extract_gee_legacy_parity.R` pull the same EE assets and emit the
**legacy** column names. `extract_gee_admin2()` / `extract_area_covariates()`
now fall back to `data/GEE/<ISO3>_legacy_parity_admin2_gee.csv` whenever a
country has no `data/<Country>_GEE_rasters/`, and warn loudly if neither exists.

**Validated** against Sierra Leone, which has both representations
(`--validate`, 14 districts): **135/135 columns pass** r >= 0.9 against the
raster-derived values *and* mean within 10%. The first version passed only
22/148 and the failures were silent scale errors — a year-over-year delta band
being averaged into the base value (this halved TRMM), a missing iSDAsoil
`exp(x/10) - 1` back-transform, and `unmask` semantics on WSF. Report:
`results/sensitivity/gee_legacy_parity_validation_sierraleone.csv`.

**Deliberately excluded** (14 of the legacy 149; see `GEE_PARITY_EXCLUDED`) —
these leave the shared vocabulary for all countries:
`gee_accessibility_2019` (EE asset corrupted server-side), the 5
`gee_esa_worldcereal_2021_*` summaries, and 7 of 9 `gee_soiltotalcarbon_*`
(its stdev bands do not match under the back-transform that works for the other
ten soil properties). `gee_ndvi_2022` is beyond AVHRR CDR coverage.

**Still open:** the 186 **cluster-buffer** `gee_` columns (monthly
`ndvi`/`temp`/`trmm` at 10/25/50 km plus `fpp_`/`tpp_`/`ghm`/`grassland`/`smod`/
`built_surface`) come from an analyst EE Code Editor export that is not in the
repo, and `fpp_`/`tpp_` are not identifiable from the names. Tanzania cannot
reproduce them, so including it in the **individual-level** LOCO costs those
covariates. That analysis is the documented *sensitivity* analysis
(`sensitivity/README.md`); the **primary area-level (Admin-2) SAE and its LOCO
are unaffected**.

### 2. The vitamin A outcome used a different inflammation adjustment

`brinda_rbp_cols()` had no Tanzania entry, so `apply_brinda_vita_binary()` fell
through a warning path and kept the configured binary — i.e. Tanzania silently
used the DHS `rbpadcrp` while the other four used two-marker BRINDA (CRP+AGP).

**Fix:** TDHS 2010 has no AGP, and CRP was assayed on only ~27% of the RBP
sample, so `brinda_adjust_rbp()` gained a **CRP-only** mode (`agp = NULL`) and
`brinda_rbp_cols()` gained a `fallback` slot: rows without CRP take the survey
agency's own adjusted RBP rather than being left raw and mixed in with corrected
ones. `brinda_country_method()` labels what each country actually got, and the
method is printed on every run. `1_GW_Tanzania_data_clean.R` now keeps
`gw_rbp_raw_umol` (the config's `gw_exclude_patterns` keeps it out of the
predictor set).

Harmonized result: **child VAD 23.3%, women 7.2%** (configured binaries: 23.9% /
7.2%). This is *not* method-identical to the other four countries — say so in
the manuscript.

### Remaining steps (exact commands, in order)

**1. Earth Engine extraction** — started 2026-08-18, multi-hour. The 13 annual
AVHRR NDVI composites over a country ~13x Sierra Leone dominate the runtime.
Resumable: every (family, year) is cached under `data/GEE/.cache_parity_TZA/`,
so re-running only pulls what is missing.

```
Rscript scripts/build_gee_legacy_parity.R Tanzania --scale-floor 250 --batch 40
```

**Do not** raise `--scale-floor` to speed it up: at 1000 m the soil `stdev` bands
and `gee_wsf_2015` come out 10-20% high (117/135 pass vs 135/135 at 250 m).
Expect `135/136 expected columns present` (`gee_ndvi_2022` is beyond AVHRR CDR
coverage) and `data/GEE/TZA_legacy_parity_admin2_gee.csv` at the end.

**2. Confirm what admitting Tanzania costs the shared vocabulary**

```
Rscript scripts/check_gee_legacy_parity.R
```

Prints, per country, how many variables the intersection loses by including it,
and writes `metadata/gee_legacy_common_vars.txt`.

**3. Full pipeline rebuild.** ~14 h sequential, ~4-6 h at 4 workers. The fixes
touch `build_outcome_dataset()`, `compute_svy_admin2()` and
`extract_gee_admin2()`, so those targets and everything downstream re-run for all
five countries. The four raster countries' `gee_admin2_*` / `area_covariates_*`
(the slowest prerequisites, ~1.5 h) were pre-built on 2026-08-18 and should be
skipped.

```
PIPELINE_MODE=full TARGETS_WORKERS=4 Rscript -e "targets::tar_make_future(workers = 4)"
```

The `gee_parity_stamp_*` targets re-check the parity CSV's md5 on every run, so
step 1 finishing after a partial rebuild correctly re-triggers Tanzania's
covariates rather than serving a covariate-free cached snapshot.

**4. Verify Tanzania actually reached the results** (it never had before):

```
Rscript -e "for (f in c('area_loco_comparison.csv','transportability_area_loco_metrics.csv','benchmarks_all.csv')) { p <- file.path('results/tables', f); if (file.exists(p)) { d <- read.csv(p); cat(f, ':', sum(grepl('anzania', unlist(d), ignore.case = TRUE)), 'Tanzania cells
') } }"
```

Also check the run log for `[pool:domain]` lines: `gee_` should read
`common=135`, not `common=0`.

**5. Interpret.** Fill in `docs/transportability_loco_methods.md` §3 with the
5-country numbers and carry §4.2's caveats (CRP-only adjustment, household
weight, 2010 survey year, Tanzania's 6:1 sample dominance in the unweighted
individual-level pool) into the manuscript.

---

## Required to run (blocking)

- [x] **Clean script — DONE & VALIDATED (2026-07-01).** Ran on this machine.
      OB↔PR match 99.4%; VAD weighted child 24.6% / women 7.7% (plausible).
      `Tanzania_GMS_cleaned.rds` written (16,977 rows). Two bugs fixed: R
      segfaults on the PR recode (→ `0_extract_pr_crosswalk.py`, pandas); RBP is
      mg/L not µmol/L (→ ÷21.2 before the <0.70 cutoff — confirm vs biomarker DOC).
- [ ] **Run the merge script.** With no extractions yet it produces a minimal
      merged dataset (gw_ outcome + Admin1/2 from GADM only). Confirm the GADM
      admin-2 join matches most clusters (`gadm41_TZA_2` downloads on first run).
- [ ] **DHS admin-2 aggregation (`dhs2010_`)** — run
      `src/DHS/DHS_admin2_aggregation.R` wired for Tanzania 2010 (IR/KR/PR/HR
      are all present). Analog: `src/SierraLeone/SL_DHS_aggregation_admin2.R`.
- [x] **GEE admin-2 zonal means — DONE.** `data/GEE/Tanzania_2010_admin2_gee.csv`
      extracted via the EE API (186 districts x 162 `gee_a2_` cols, 0 all-NA;
      ~34/35 layers — only healthcare-Accessibility skipped, GEE-side corrupted
      asset). Cluster buffers still pending `Tanzania_GMS_cleaned.rds`.
      **Auto-extraction scaffold** (no manual Code-Editor export):
      - `src/GEE/gee_layer_manifest.R` — canonical layer→EE-asset manifest
        (from `data/GEE/GEE_export_metadata.xlsx`).
      - `src/GEE/extract_gee_ee_api.R` — rgee extractor (admin-2 polygons +
        cluster buffers, server-side `reduceRegions`).
      - `src/Tanzania/4_TZ_GEE_extract.R` — Tanzania 2010 runner →
        `data/GEE/Tanzania_2010_admin2_gee.csv` + `TZ2010_buffers_*.csv`.
      **BLOCKED ON: Earth Engine API access.** Follow
      `src/GEE/README_GEE_API_SETUP.md` (register a Cloud project, `rgee`,
      `ee_Initialize()`), then run the script. First run needs validation:
      watch per-layer `ok:` logs, the deprecated Oxford-MAP-LST substitute, and
      the column-name parity caveat (EE band names differ from the legacy
      terra-extracted names → matters for pooled/LOCO, not within-Tanzania).
- [x] **SoilGrids admin-2 — DONE (ran 2026-07-01).** 186 districts x 7 props
      (`Tanzania_soilgrids_admin2.csv`), no all-NA cols, sane values. No download needed: SoilGrids is
      fetched on-demand from the global ISRIC VRT (`/vsicurl/`) and cropped to
      Tanzania's bbox. Tanzania is now in the `countries` list of
      `scripts/build_soilgrids_admin2.R` (code `"TZA"`, matching the merge's
      GADM so the exact Admin1/Admin2 join lines up). Run:
      `Rscript scripts/build_soilgrids_admin2.R Tanzania` →
      `data/SoilGrids/Tanzania_soilgrids_admin2.csv`. (Soil micronutrients are
      the single most transportable covariate for this project — worth running
      early.)
- [ ] **Malaria Atlas (`MAP_`)** — analog `src/SierraLeone/SL_MAP_download.R`.
- [x] **IHME (`ihme_`) — DONE, run it.** No download needed: Tanzania is already
      in the global LBD CSVs in `data/IHME/<family>/`. Run
      `src/Tanzania/3_TZ_IHME_clean.R`, which calls the new generic builders
      `src/IHME/build_ihme_admin2.R` + `src/IHME/build_ihme_admin1.R` and writes
      `data/IHME/tanzania_2010_merged_IHME_data.csv` and
      `..._admin1_data.csv` (year 2010).
      - Check the printed **PARITY** line vs Ghana to confirm column alignment.
      - **Known gap:** the repo has no admin-2 CSVs for **Malaria** or
        **education**, so those `ihme_` admin-2 columns will be absent for
        Tanzania (malaria is covered by `MAP_`, education by DHS proxies). To
        close the gap, obtain the LMIC malaria + education *ADMIN_2* CSVs from
        GHDx, or extract from their GeoTIFFs.

## Optional domains — status

- [x] **WFP prices (`wfp_`) — DONE.** `data/food_price/wfp_food_prices_tza.csv`
      downloaded from HDX (64k rows, covers 2010). Merge block parses it
      (nearest-market join); domain added to config.
- [x] **Food security (`fsec_`) — WIRED.** HFID covers Tanzania (from 2011;
      nearest-year fallback for the 2010 survey). Tanzania added to the config
      maps in `R/food_security.R`; merge block calls `merge_food_security()`.
      Cadre Harmonisé N/A (West/Central Africa only).
- [~] **LSMS (`lsms_`) — SCAFFOLDED, needs a registered download.** Admin-1
      (coarse), optional — only Ghana/Gambia use it. Requires the Tanzania
      **National Panel Survey (NPS 2010-11, LSMS-ISA)** from the World Bank
      Microdata Library (free account + terms). Then run
      `src/Tanzania/TZ_LSMS_clean.R` (confirm the NPS file/variable names in its
      header) → `data/LSMS/Tanzania_LSMS_clean.RDS`. Merge block + config domain
      already wired.
- [ ] MICS (`mics_`) — **N/A**, Tanzania has no MICS round in scope.
- [ ] FluNet (`flunet_`) — optional; other countries mostly omit.

## Validation

- [ ] Add Tanzania to `scripts/validate_merged_datasets.R` (`countries` list)
      and confirm per-prefix column counts + NA rates match the other countries.
- [ ] Confirm `_targets` picks up Tanzania automatically (the pipeline loops over
      `get_country_configs()`), then check LOCO targets include it.

---

## Key caveats to carry into the manuscript

1. **Weighting.** No micronutrient-subsample weight exists in TDHS 2010; the
   clean script uses the household weight (`HV005/1e6`). The biomarkers were a
   subsample, so this is the known subsample-weight limitation (cf. Gambia).
2. **Inflammation adjustment is CRP-only** — TDHS 2010 has CRP but no AGP, so the
   full two-marker BRINDA correction isn't possible. The OB file's `rbpadcrp`
   (RBP adjusted for CRP) is used. Check interaction with
   `apply_brinda_vita_binary()` in `R/data_prep.R`.
3. **Survey year 2010** is older than the 2013–2018 panel; proxies are
   year-stamped so this is handled, but it widens the temporal spread.
4. **Excluded biomarkers:** sTfR (iron) and urinary/salt iodine are present in
   the raw data but excluded by design — revisit only if a separate sTfR-based
   iron track or a descriptive iodine product is wanted.
