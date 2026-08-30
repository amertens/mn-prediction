# Admin-3 boundaries & GEE covariate-completeness audit

_Read-only investigation, 2026-06-23. Targets store: `_targets_full`._

## 1. Admin-3 (GADM 4.1) availability for the four study countries

Queried via `geodata::gadm(iso, level=3)` against GADM 4.1.

| Country      | ISO3 | Admin-1 | Admin-2 | Admin-3 | Admin-3 available? |
|--------------|------|---------|---------|---------|--------------------|
| Gambia       | GMB  | 6       | 37      | —       | **No** (GADM has only levels 0-2) |
| Ghana        | GHA  | 16      | 260     | —       | **No** (GADM has only levels 0-2) |
| Sierra Leone | SLE  | 4       | 14      | 153     | **Yes** — 153 chiefdoms |
| Malawi       | MWI  | 28      | 256     | 3126    | **Yes** — 3126 units (Traditional Authorities / wards) |

Notes:
- Gambia and Ghana return "this file does not exist" for level 3 from the GADM server — admin-3 simply does not exist in GADM 4.1 for these countries.
- Sierra Leone admin-3 = chiefdoms (153). Already cached at `data/gadm/gadm/gadm41_SLE_3_pk.rds` and `data/external_cache/gadm/gadm/gadm41_SLE_3_pk.rds`.
- Malawi admin-3 = 3126 units — an order of magnitude finer than admin-2 (256). Very small areas; likely far below survey design resolution (effective n problem) if used as a modeling unit.

### Is admin-3 already wired into the repo?

**Sierra Leone: yes, fully wired** (and referenced in recent git history — commit `ca8dce9` "Add Sierra Leone Admin-3 (chiefdom) map layer").
- Prototype/analysis pipeline: `sandbox/21_sierraleone_admin3.R` — spatially joins DHS clusters to the 153 GADM chiefdoms, aggregates biomarkers, runs LOCO mixed models at admin3 vs admin2. Outputs:
  - `results/transportability/sierraleone_admin3_svy.rds`
  - `results/transportability/sierraleone_admin3_gee.rds`
  - `results/transportability/sierraleone_cluster_to_admin3.rds`
  - `sandbox/results/21_sl_admin3_vs_admin2.csv`
- Dashboard layer: `dashboard/data-raw/_build_sl_admin3.R` — reconstructs chiefdom prevalence from individuals + survey weights, fits a Fay-Herriot/empirical-Bayes area model on chiefdom GEE covariates to cover all 153 chiefdoms, fetches GADM level-3 boundaries, writes `dashboard/data/admin3_predictions.rds` and `dashboard/data/admin3_boundaries.rds`.
- Dashboard wiring: `admin3` referenced in `dashboard/global.R`, `dashboard/R/mod_map_explorer.R`, `dashboard/R/mod_resolution.R`, `dashboard/R/fct_helpers.R`, `dashboard/data-raw/01_prepare_dashboard_data.R`.

**Gambia, Ghana, Malawi: no admin-3 usage** anywhere. The core pipeline operates at admin-2 only (`R/external_data.R` loads `geodata::gadm(..., level = 2)`; `R/config.R` stores only `gadm_code` per country, no level field). Gambia/Ghana cannot have admin-3 (does not exist in GADM); Malawi *could* (3126 units exist) but is not implemented.

## 2. GEE covariate completeness audit

Read each country's `gee_admin2_<country>` target from `_targets_full`.

| Country      | Admin-2 rows | `gee_` columns | All-NA cols | >80%-NA cols |
|--------------|--------------|----------------|-------------|--------------|
| Gambia       | 37           | 556            | 1           | 0            |
| Ghana        | 260          | 559            | 1           | 0            |
| Malawi       | 243          | 539            | 1           | 0            |
| Sierra Leone | 14           | 543            | 1           | 0            |

- The single all-NA column in every country is **`gee_ndvi_2022`** — empty everywhere (a dead/failed extraction; harmless since it's all-NA universally, no transport asymmetry).
- No other columns exceed 80% NA in any country. Within-country coverage is good.

### Cross-country comparison (transport-coverage gaps)

The raw column-name universe across countries is 1446 distinct `gee_` names, with ~1297 differing — **but this is almost entirely a survey-year artifact**: each country's covariates are extracted for that country's survey year(s), so year-suffixed columns (e.g. `gee_..._2017` vs `_2018`) differ by name even though the underlying variable exists for all.

After collapsing 4-digit year tokens, there are **341 distinct variable bases**, of which:
- **323 (95%) are present and populated in all four countries** — the transport-safe covariate core.
- **18 bases are populated in some countries but absent/empty in others** — and 17 of these 18 are **ESA WorldCereal cropland tiles**, plus 1 landcover layer:

| Variable base (year-collapsed) | Populated in | Gap in |
|---|---|---|
| `gee_esa_worldcereal_*_32121_*` (maize_main, maize_second, wintercereals) | Gambia, Ghana, SLE | Malawi |
| `gee_esa_worldcereal_*_34120_*` (maize_main, wintercereals) | Gambia, Ghana | Malawi, SLE |
| `gee_esa_worldcereal_*_7091_*` (maize_main, maize_second, wintercereals) | Ghana, SLE | Gambia, Malawi |
| `gee_esa_worldcereal_*_10162/12048/30112/9030_*` (maize/wintercereals) | Malawi | Gambia, Ghana, SLE |
| `gee_landcoverlayers_*_change_confidence` | Gambia, Ghana | Malawi, SLE |

**Interpretation:** the ESA WorldCereal numeric codes (`32121`, `34120`, `7091`, `10162`, ...) are **agro-ecological / processing tile identifiers** — a country only receives the tiles that geographically overlap its territory. So these "gaps" are structural (different countries fall in different WorldCereal tiles) rather than failed extractions. They are *not* directly comparable cross-country and should be treated as country-specific, not used as shared transport predictors.

The only non-cropland gap is `gee_landcoverlayers_*_change_confidence`, populated (0% NA) for Gambia and Ghana but the column is entirely absent for Malawi and Sierra Leone.

### Bottom line for cross-country transport

- The shared, populated covariate core (~323 variable bases) is broad and consistent across all four countries — covariate coverage is **not** a meaningful source of transport weakness.
- The only genuine cross-country asymmetries are: (a) `gee_ndvi_2022` (all-NA everywhere — drop it), (b) ESA WorldCereal tile columns (structural, tile-geography driven — exclude from shared predictor set), and (c) `gee_landcoverlayers_*_change_confidence` (missing for Malawi & SLE).
- These are minor relative to the documented driver of transport failure (cross-survey biomarker *level* offsets), so closing GEE coverage gaps will not materially help transportability.
