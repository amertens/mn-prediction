# RA Ingestion Specs — Complementary Predictor Data Sources

**Date:** 2026-06-15 · **Companion to:** `docs/pipeline_improvements_todo.qmd` §"New predictor data sources",
`docs/feature_engineering_and_sae_notes.md`.

These four sources fill the gaps the predictor audit identified — mechanistic predictors for the **dietary
micronutrients (folate, B12, zinc)** and for **helminth-driven iron deficiency**, which the current sources
(malaria/environment-heavy) serve worst. Each spec is self-contained: access, resolution, output schema,
join keys, caveats, and acceptance criteria.

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

**Access (API key required, free).**
- Request a key: <https://espen.afro.who.int> (or `ESPENAPI::espen_key_setup()`); set `ESPEN_API_KEY`.
- Data API endpoint: `https://espenjapapi.afro.who.int/api/data?iso2=<C>&disease=<sth|sch>&level=<iu|sitelevel>&start_year=&end_year=&api_key=<KEY>`
- R package: `remotes::install_github("olatunjijohnson/ESPENAPI")`, `ESPEN_API_data(iso2=, disease=, level=)`.
- **Runnable puller already drafted:** `scripts/build_espen_admin2.R` (prefers `ESPENAPI`, else direct httr;
  pulls STH+SCH at IU level, fuzzy-matches IU→GADM Admin2, writes `data/ESPEN/<Country>_espen_admin2.csv`).

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

**Resolution.** ~10 km grid → Admin2 zonal mean (and cluster-buffer if desired).

**Output schema.** `Admin2, spam_maize_prod, spam_pulses_prod, spam_vegetables_prod, spam_roots_prod, spam_cereals_share, spam_pulses_share, ...` (production and/or share of total).

**Caveats.** SPAM reference years are discrete (2005/2010/2017/2020) — pick the nearest to the survey year and
record it. Production ≠ consumption (trade), but a reasonable local-availability proxy.

**Acceptance.** Crop layers non-NA for ≥80% of Admin2; crop shares sum sensibly; pulses/vegetable shares show
some spatial variation within each country.

**Effort.** ~3 h (reuses the GEE zonal-extraction code).

---

## Sequencing recommendation

1. **ESPEN** first — the only subnational, mechanistically-iron-causal addition; `scripts/build_espen_admin2.R`
   is ready, just needs a key. It immediately answers "do helminths improve iron prediction?" via
   `sandbox_fe/21_espen_iron_test.R`.
2. **GFDx** second — trivial download, and the only source that plausibly explains **folate** burden, which
   nothing else captures; pairs with the multi-country pooled model (roadmap P1).
3. **GDD** and **MapSPAM** as the dietary layer matures.

Remember the effective-n caveat: national sources (GFDx, GDD) help **cross-country transfer and the pooled
model**, not within-country Admin2 maps; subnational sources (ESPEN, MapSPAM) help both.
