# RA Ingestion Specs — Complementary Predictor Data Sources

**Date:** 2026-06-15 (rev. 2026-07-28) · **Companion to:** `docs/pipeline_improvements_todo.qmd`
§"New predictor data sources", `docs/feature_engineering_and_sae_notes.md`.

Sources **1–4** (ESPEN, GFDx, GDD, MapSPAM) fill the gaps the predictor audit identified — mechanistic
predictors for the **dietary micronutrients (folate, B12, zinc)** and for **helminth-driven iron
deficiency**, which the current sources (malaria/environment-heavy) serve worst. Sources **5–7** (VAS,
FAOSTAT, FPN CoAHD/CoNA) were added 2026-07 from the West Africa data-landscape review as open,
broad-coverage **national causal-driver** proxies (supplementation, food supply, diet affordability) that
strengthen the cross-country pooled / transfer model. Each spec is self-contained: access, resolution,
output schema, join keys, caveats, and acceptance criteria.

## Access mechanism (R) — at a glance

Whether each source can be pulled unattended from R (true query API + package) vs. a stable-URL file
download vs. a manual portal export. Prioritise the true-API sources — they are the cheapest RA wins.

| # | Source | Prefix | R access mechanism | Query API? | Key |
|---|--------|--------|--------------------|------------|-----|
| 1 | **ESPEN** (helminths) | `espen_` | **Manual portal download** (no key) → `build_espen_admin2.R` ingests `data/ESPEN/raw/`; ESPENAPI only if you already hold a key | API keys **PAUSED** by WHO | — |
| 2 | **GFDx** (fortification) | `gfdx_` | `download.file()` on the published dataset CSV | No (file download) | No |
| 3 | **GDD** (diet intake) | `gdd_` | `nutriR` pkg for inadequacy distributions; bulk file for full estimate set | Partial | No |
| 4 | **MapSPAM** (crop mix) | `spam_` | base `download.file()` on the Dataverse access API (file ids in `build_mapspam_admin2.R`) → CSV → zonal | **Yes** (no pkg needed) | No |
| 5 | **UNICEF VAS** (vit-A supp.) | `vas_` | `WDI::WDI("SN.ITK.VITA.ZS")` (WB mirror of the UNICEF estimate) | **Yes** (WB API) | No |
| 6 | **FAOSTAT** (food supply) | `fao_` | base `download.file()` of the public FBS bulk zip → `fread` (the `FAOSTAT` pkg forces a needless login) | **Bulk file** (no login) | No |
| 7 | **FPN** (CoAHD/CoNA) | `fpn_` | **Manual DataBank CSV export** → `build_fpn_affordability.R` ingests it (the `source=88` data API is dead; migrated to Data360) | API **migrated/closed** | No |

**Fully scriptable & verified running:** MapSPAM, VAS, FAOSTAT. **Manual one-time download, then scripted
ingest:** ESPEN (WHO paused API keys), FPN (WB retired the `source=88` data API; `CoNA_*` codes replaced by
`CoHD_*`/`CoCA_*`). **Still spec-only:** GFDx, GDD.

**Build status (2026-07-28) — built, wired, and RUN where possible:**
- **MapSPAM** (`spam_`) — `scripts/build_mapspam_admin2.R`, base-R Harvard-Dataverse download. **Ran for all
  5 countries** (Gambia 36 / Ghana 244 / SL 14 / Malawi 222 / Tanzania 181 Admin-2 rows). Join coverage vs
  survey Admin-2: Ghana 69/75, Malawi 83/89, Gambia 29/30 (misses = no-cropland/urban districts → NA).
- **VAS** (`vas_`) — `scripts/build_vas_national.R`, `WDI`. **Ran, all 5** (coverage at exact survey years:
  Gambia 30%, Ghana 50%, SL 99%, Malawi 91%, Tanzania 99%). Join coverage 100% (30/30, 75/75, 89/89).
- **FAOSTAT** (`fao_`) — `scripts/build_faostat_supply.R`, direct FBS bulk. **Ran, all 5**, 8 supply fields
  at exact survey years. Join coverage 100%.
- **ESPEN** (`espen_`) — `scripts/build_espen_admin2.R` ready; **awaits the one-time manual portal download**
  (keys paused). **FPN** (`fpn_`) — `scripts/build_fpn_affordability.R` ready; **awaits a one-time DataBank
  CSV export**.
- All five domains are registered in `R/config.R` (all 5 country blocks) and joined via `file.exists()`-guarded
  blocks in every `src/<Country>/2_GW_<Country>_data_merge.R` (+ metadata). Missing CSVs → domain simply
  skipped, pipeline unaffected.
- **Packages:** `WDI`/`wbstats`/`FAOSTAT` now installed (though FAOSTAT/FPN bypass their packages). `dataverse`
  not needed. GFDx/GDD (sources 2–3) remain spec-only.

## Conventions for every source (read first)

- **Output:** one CSV per country, `data/<SOURCE>/<Country>_<source>_admin2.csv` (and, where the raster is
  finer than Admin2, also a cluster-buffer CSV keyed on cluster id).
- **Key column:** `Admin2` matching **GADM `NAME_2`** (the pipeline's join key). Use the project's fuzzy
  matcher (Jaro-Winkler ≤ 0.15, mirrors `R/food_security.R::build_name_lookup()` / the IHME merge) — never
  assume exact name equality. Carry the raw source name in a `*_src_name` column for audit.
- **Column prefix:** pick a unique domain prefix (`espen_`, `gfdx_`, `gdd_`, `spam_`) so the columns enter a
  new domain. Wire by adding `DOMAIN = list(prefix = "<prefix>")` to each country block in `R/config.R` and a
  left-join on `Admin2` near the end of each `src/<Country>/2_GW_<Country>_data_merge.R`.
- **No imputation at ingest.** Unmatched Admin2 → `NA` (the SL pipeline imputes later).
- **Verify on first run.** Print the source's column names and confirm the field mapping before trusting it.
- **Countries / survey years:** Gambia 2018, Ghana 2017, Sierra Leone 2013, Malawi 2016. Pull the survey-year
  value (or nearest available; record the year used in a `*_year` column).
- **Acceptance criteria (all sources):** ≥ 80% of Admin2 rows non-NA per country (national-only sources
  exempt — they are constant within country by design); values in a plausible range; a quick
  `cor(value, area outcome prevalence)` is finite and not absurd.

---

## 1. ESPEN — soil-transmitted helminths & schistosomiasis  (P2, iron)

**Why.** Hookworm causes chronic intestinal blood loss → iron-deficiency anaemia; schistosomiasis (esp.
*S. haematobium*) likewise. Absent from the current predictors and **mechanistically the most relevant**
missing covariate for the iron outcomes. Subnational → helps within-country *and* transfer.

**Access — the API is currently CLOSED; use the public no-key download page.**
- **WHO/AFRO has PAUSED issuing ESPEN API keys** (confirmed 2026-07 via the ESPENAPI package README:
  *"ESPEN has paused issuing new API keys. Existing keys may also have expired."*). The keyed endpoint
  (`https://espenjapapi.afro.who.int/api/data?...`) and the `/espen-api-platform` page are gated. Requesting
  access, if it reopens: email **ntd.espen@who.int** or see <https://espen.afro.who.int/espen-api-platform>.
- **The route that works today (no key):** the portal's public **Download data** page,
  <https://espen.afro.who.int/maps-data/data-query-tools/download-data> — select Country, Disease = STH (then
  SCH), Level = IU, full year range, click **Download Data**, and save each CSV into `data/ESPEN/raw/`.
  (The public site's data lives behind an unauthenticated dashboard API, `/api/espen-dashboard/<disease>/…`,
  but it returns Drupal render-JSON, not clean CSV — the download page is the reliable extract.)
- WHO GHO (`ESPENAPI::gho_ntd_data()`) is the only no-key *API*, but it covers **oncho/trachoma at country
  level only** — useless for our STH/SCH-subnational need.
- **Runnable puller (rewritten for this reality):** `scripts/build_espen_admin2.R` now **ingests the manual
  raw downloads** from `data/ESPEN/raw/` (auto-detects country+disease from filename/contents), fuzzy-matches
  IU→GADM Admin2, and writes `data/ESPEN/<Country>_espen_admin2.csv`. It still uses `ESPENAPI` automatically
  if you happen to hold a valid key. **Outstanding step: the one-time manual portal download per country.**

**Resolution.** Implementation Unit (IU ≈ district ≈ Admin2). Years vary by survey campaign; aggregate
multiple surveys per IU to the mean (the script does this).

**Output schema.** `Admin2, espen_sth_prev, espen_sch_prev, helminth_sth_any_prev` (proportions 0–1).
Consider separate hookworm/ascaris/trichuris if the API exposes species-level fields (verify column names).

**Caveats.** Coverage is campaign-driven and uneven — some districts unmapped (→ NA, fine). IU names need
fuzzy matching. **Confirm the prevalence field name on first run** (`Prevalence`/`prev`/endemicity class).

**Acceptance.** STH prevalence present for a majority of districts in ≥3 of 4 countries; hookworm/STH
prevalence positively associated with area iron-deficiency prevalence in at least the higher-burden
countries (test via `sandbox_fe/21_espen_iron_test.R`, which auto-activates once the CSVs exist).

**Effort.** ~4 h (mostly key request + IU→GADM name reconciliation).

---

## 2. GFDx — food fortification  (P2, folate/B12/iron/vitamin A/zinc)

**Why.** Mandatory/voluntary fortification of flour (folic acid + iron + B12 + zinc), oil & sugar (vitamin A),
salt (iodine) is a *direct causal* determinant of population micronutrient status — especially **folate**,
which the current predictors barely address.

**Access (easy download, no key).**
- Datasets: <https://www.fortificationdata.org/download-datasets/> and full compendium
  <https://fortificationdata.org/full-gfdx-datasets/> — row-based CSVs, 196 countries.
- Pull: fortification *standards* (which nutrients, at what levels, mandatory vs voluntary) and *coverage*
  (% of food industrially processed / fortified) for the four countries.
- **R access:** no query API. Script the file fetch with `download.file()` on the published dataset CSV,
  then filter to the study countries — the download-page URLs change between releases, so **confirm the
  current CSV link on the download page first** (not a stable API contract).

**Resolution.** **National** (with year). Constant within country → enters as a country-level covariate.

**Output schema.** `Admin2` (every Admin2 in the country gets the national value),
`gfdx_flour_folic_acid, gfdx_flour_iron, gfdx_oil_vitA, gfdx_salt_iodine, gfdx_*_mandatory, gfdx_*_coverage_pct, gfdx_year`.

**Caveats.** National-only → **does not improve within-country area discrimination**; its value is for the
**multi-country pooled / transfer** model (explains why one country's folate burden differs from another's).
Standards ≠ actual intake; pair with GDD (#3) for intake.

**Acceptance.** Each country has a fortification profile for its survey year; the folate/folic-acid field
varies across the four countries (it should — fortification mandates differ).

**Effort.** ~3 h.

---

## 3. Global Dietary Database (GDD, Tufts) — modelled nutrient intake  (P3, B12/folate/zinc/vit A)

**Why.** Direct modelled dietary intake of specific micronutrients — a far better dietary signal than the
Ghana-only LSMS expenditure-share proxies, and available for all countries.

**Access (public).** <https://tuftsfoodismedicine.org/project/global-dietary-database/>; estimates and code on
the project GitHub; subnational distributions via the companion `nutriR` package / Passarelli et al. 2024
(Lancet GH) inadequate-intake estimates.

**R access:** partial. The Passarelli 2024 inadequate-intake distributions are packaged in `nutriR`
(`remotes::install_github(...)`) and usable directly in R; the full 15-nutrient GDD estimate set is a bulk
file (often gated behind a project request/email), so **script the `nutriR` pull and treat the full estimate
set as a one-time manual download** if the packaged distributions are insufficient.

**Resolution.** National × age-sex (some urbanicity/education strata); 185 countries, 15 micronutrients.
Largely **national** for our purposes.

**Output schema.** `Admin2, gdd_intake_vitA, gdd_intake_folate, gdd_intake_b12, gdd_intake_zinc, gdd_intake_iron, gdd_inadequacy_<nutrient>, gdd_year` (match age-sex group to the outcome population: women vs preschool children).

**Caveats.** National/age-sex → within-country constant (transfer value, like GFDx). Match the age-sex group
to the outcome (e.g. women_folate → women of reproductive age).

**Acceptance.** Intake/inadequacy values present for all four countries and the right age-sex groups; values
in plausible physiological ranges.

**Effort.** ~4 h.

---

## 4. MapSPAM / HarvestChoice — crop composition  (P3, vitamin A / dietary)

**Why.** You have cropland *class* (ESA WorldCereal) and productivity (WAPOR/GPP) but not crop *composition*.
The mix of staples vs pulses (iron/zinc) vs vegetables/orange-fleshed crops (provitamin A) vs roots maps to
local micronutrient availability. Subnational → helps within-country too.

**Access (download).** MapSPAM gridded production / harvested-area / yield by crop and technology
(IFPRI Dataverse / mapspam.info); HarvestChoice layers (already flagged for download in §Predictors). GeoTIFF
rasters → zonal means per Admin2 with the existing `extract_gee_admin2()` machinery.

**R access:** true API, no extra package. MapSPAM 2010 v2r0 is on Harvard Dataverse
(DOI `10.7910/DVN/PRFF8V`). `scripts/build_mapspam_admin2.R` downloads the all-technology production and
physical-area **CSV** zips with base `download.file()` from the Dataverse access API
(`https://dataverse.harvard.edu/api/access/datafile/<id>` — verified 2026-07: HTTP 206, `application/zip`;
prod-csv id `3984975`, phys-area-csv id `3984973`), filters pixels to the country ISO3, spatial-joins the
5-arcmin pixel centroids to GADM Admin-2 (`sf`), and aggregates crop groups (cereals / roots / pulses /
oilcrops / vegetables) to per-Admin2 production and shares. Runs off the installed `data.table`+`sf` stack;
the `dataverse` package is **not** required.

**Resolution.** ~10 km grid → Admin2 zonal mean (and cluster-buffer if desired).

**Output schema.** `Admin2, spam_maize_prod, spam_pulses_prod, spam_vegetables_prod, spam_roots_prod, spam_cereals_share, spam_pulses_share, ...` (production and/or share of total).

**Caveats.** SPAM reference years are discrete (2005/2010/2017/2020) — pick the nearest to the survey year and
record it. Production ≠ consumption (trade), but a reasonable local-availability proxy.

**Acceptance.** Crop layers non-NA for ≥80% of Admin2; crop shares sum sensibly; pulses/vegetable shares show
some spatial variation within each country.

**Effort.** ~3 h (reuses the GEE zonal-extraction code).

---

## 5. UNICEF Vitamin A supplementation coverage  (P2, vitamin A)

**Why.** Two-dose VAS coverage in children 6–59 mo is the proximate driver of population vitamin A status and
is the covariate **GBD leans on for its modelled VAD decline**. GFDx (#2) captures vitamin A *fortification*
of oil/sugar; **supplementation is a distinct — and for under-fives, larger — channel** that nothing in the
current predictor set represents. Directly relevant to the vitamin A outcome (incl. the new Tanzania panel).

**Access (public, no key).**
- Primary R route: the **World Bank WDI mirror** of the UNICEF estimate — indicator `SN.ITK.VITA.ZS`
  ("Vitamin A supplementation coverage rate, % of children 6–59 months"). Pull with
  `WDI::WDI(indicator = "SN.ITK.VITA.ZS", country = c("GH","TZ","GM","SL","MW"))` or `wbstats::wb_data()`.
- Native source / fallback: UNICEF Data Warehouse SDMX, dataflow `UNICEF:NUTRITION(1.0)`, via `rsdmx`
  (base `https://sdmx.data.unicef.org/ws/public/sdmxapi/rest/`). **Confirm the VAS indicator code on first
  run** (two-dose full coverage) — the WDI mirror avoids this and is preferred unless the mirror is stale.

**Resolution.** National × year, children under 5. Constant within country → transfer/pooled value.

**Output schema.** `Admin2` (national value broadcast to every Admin2), `vas_vita_coverage_pct`, `vas_year`.

**Caveats.** National-only → helps the pooled/transfer model, **not within-country Admin2 discrimination**.
The WDI series has **gaps in the most recent years** (UNICEF reporting changes) — take the value at (or
nearest to) each survey year and record it in `vas_year`. Match to the **child** vitamin A outcome, not women.

**Acceptance.** A coverage value present for each country at its survey year or nearest available; values in
0–100% and varying across countries.

**Effort.** ~1 h (single `WDI` pull).

---

## 6. FAOSTAT Food Balance Sheets / SUA — national food & nutrient supply  (P3, zinc/iron/vit A/folate)

**Why.** National per-capita supply of energy, protein, and food-group quantities is the upstream
*availability* signal for dietary micronutrients — the same backbone FAO/GBD use for the **zinc**
inadequate-intake model. Complements GDD (#3): supply (FAOSTAT) vs modelled intake (GDD).

**Access (public, no key, no login).** The installed `FAOSTAT` CRAN package unexpectedly forces an API login
(`get_faostat_bulk` → `faostat_login`), so the builder bypasses it and fetches the **public FBS bulk zip**
directly: `https://bulks-faostat.fao.org/production/FoodBalanceSheets_E_All_Data_(Normalized).zip` (~55 MB,
verified HTTP 200), unzips the normalized CSV, and `fread`s it. Filter to the study countries + survey years.
Runs fully unattended.

**Resolution.** National × year. Constant within country.

**Output schema.** `Admin2` (national value broadcast), `fao_supply_kcal`, `fao_supply_protein_g`,
`fao_supply_<foodgroup>_kg` for nutrient-relevant items (pulses, vegetables, animal-source foods, oils),
`fao_year`.

**Caveats.** National-only → transfer/pooled value. Supply ≠ intake (waste, distribution, trade). Item and
element codes must be selected deliberately — **verify the item/element mapping on first run.**

**Acceptance.** Supply values for all countries and survey years; per-capita quantities in plausible ranges.

**Effort.** ~2 h.

---

## 7. World Bank "Food Prices for Nutrition" — CoAHD / CoNA  (P3, all micronutrients, economic access)

**Why.** Cost and affordability of a **nutrient-adequate** (CoNA) and a healthy (CoAHD) diet — the *economic
access* dimension of micronutrient-dense foods, distinct from the WFP raw commodity prices already in the
pipeline. A plausible determinant of dietary-micronutrient status (folate/B12/zinc) that no current source
captures.

**Access — API migrated; use a manual DataBank export (no key).** Testing 2026-07 found the classic
`source=88` **data** endpoint (`/v2/country/GHA/indicator/<code>?source=88`) returns *"indicator not found"*
even for valid codes — FPN's data moved to the **Data360** platform, whose API returned generic retrieval
errors in testing. Also the **`CoNA_*` (nutrient-adequate) codes are RETIRED** in FPN 4.0; the live
`source=88` catalogue lists **`CoHD_*`** (healthy diet) and **`CoCA_*`** (calorie-adequate) only. So the
builder ingests a one-time **manual CSV export** (reliable today): from
<https://databank.worldbank.org/source/food-prices-for-nutrition> select the countries + the CoHD/CoCA series
+ all years → Download CSV → drop into `data/FPN/raw/`. `build_fpn_affordability.R` auto-detects the DataBank
wide layout (Country Code, Series Code, `YYYY [YRYYYY]` columns) and picks each country's survey-year value.

**Indicator codes (live 2026-07, `source=88` catalogue).** Diet cost/affordability: `CoHD_PPP`, `CoHD_LCU`,
`CoHD_fexp`, `CoHD_headcount` (% unable to afford a healthy diet), `CoHD_CoCA`; energy-sufficient:
`CoCA_PPP`, `CoCA_headcount`. Micronutrient-relevant cost shares: `CoHD_lns_prop` (legumes/nuts/seeds),
`CoHD_asf_prop` (animal-source foods), `CoHD_f_prop` (fruits). (`CoNA_*` no longer exist.)

**Resolution.** National × year. Constant within country.

**Output schema.** `Admin2` (national value broadcast), `fpn_cona_ppp`, `fpn_cona_headcount`,
`fpn_cohd_headcount`, `fpn_year`. Prefer **CoNA** for nutrient adequacy; `headcount` is the % unable to afford.

**Caveats.** National-only → transfer/pooled value. Series is short and anchored to ICP benchmark years; it
is a diet-cost index, **not micronutrient-specific**.

**Acceptance.** CoNA affordability present for all countries; `*_headcount` in 0–100%.

**Effort.** ~2 h.

---

## Sequencing recommendation

1. **ESPEN** first — the only subnational, mechanistically-iron-causal addition; `scripts/build_espen_admin2.R`
   is ready, just needs a key. It immediately answers "do helminths improve iron prediction?" via
   `sandbox_fe/21_espen_iron_test.R`.
2. **National causal-driver block — VAS, FAOSTAT (both DONE, unattended), then FPN (manual export).** VAS
   (`WDI`) and FAOSTAT (direct FBS bulk) already ran for all 5 countries; FPN needs the one-time DataBank CSV
   export. These populate the national covariates the pooled/transfer model most needs (VAS for vitamin A in
   particular).
3. **GFDx** — trivial file fetch, and with VAS it completes the vitamin A pathway (fortification +
   supplementation); GFDx is also the only source that plausibly explains **folate** burden.
4. **GDD** and **MapSPAM** as the dietary layer matures.

Remember the effective-n caveat: national sources (VAS, GFDx, GDD, FAOSTAT, FPN) help **cross-country transfer
and the pooled model**, not within-country Admin2 maps; subnational sources (ESPEN, MapSPAM) help both. Given
that transport fails mainly on absolute *level*, the national causal drivers are worth having despite being
constant within country.
