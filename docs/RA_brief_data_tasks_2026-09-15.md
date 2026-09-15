# RA brief — data downloads, research and cleaning after the 2026-09-15 rebuild

**Status:** ready to start. **Estimated effort:** 2–4 weeks depending on how many
tiers are taken. **Owner of the specs:** Andrew. Every task says what to
download or build, where it goes, and what "done" looks like, so it can be
picked up without the conversation that produced it.

Context in two sentences. The shared Admin-2 predictor set
(`data/covariates/harmonized/predictors_admin2_shared.csv`, 552 predictors
after the HCES and RTFP blocks of 2026-09-15, four countries) is built from
public sources by scripts; the sources that are
still missing are the ones that sit behind a registration, a manual portal
export, or a PDF table. The sandbox results say the missing information that
would matter most is **household diet** (no HCES block outside Ghana) and
**more biomarker surveys** (each training country buys ~0.05 transported rank
accuracy), so those come first. Background: `docs/findings/AB-01_SOURCE_ADDBACK_2026-09-15.md`.

Conventions for every task
- Raw downloads go under `data/<SOURCE>/<Country>/` with the provider's
  documentation kept next to the data (codebook, questionnaire, report PDF).
  `data/` is git-ignored, so nothing is committed by accident.
- Add one row per acquired file to `metadata/external_provenance.csv`
  (source, access method, URL, date, version, licence, who downloaded).
- Never rename provider files; put derived tables in a separate folder.
- Region names must be matched to GADM 4.1 names (`data/admin_boundaries/`);
  keep the provider's original label in a `*_src_name` column.
- If a download needs an account, register in your own name; do not share
  credentials in the repo.

---

## Status after the first folder review (2026-09-15, `data/RA_2026-09/`)

The folder delivered on 2026-09-15 (`Proxy VMD database`, now filed as
`data/RA_2026-09/` with the used archives unpacked under `extracted/`) closed
three Tier-1 items and part of a fourth. What was built from it:

| Delivered file | Used for | Result |
|---|---|---|
| MICS microdata for Sierra Leone 2017 and Malawi 2013-14 (added to `data/MICS/` 2026-09-15) | `scripts/covariates/build_mics_admin2_block.R` | 30 `mics_` columns in the shared set (salt iodine test, WASH, wealth, anthropometry, IYCF, illness, women's education, IPTp) for all four countries; Ghana at region level and Malawi at district level until their GPS arrive (T1.2); T3.2 superseded |
| `LSMS.zip` (IHS4 2016-17, IHS 2015/16, SLIHS 2018) + `MICS GPS.zip` (SL 2017, Gambia 2018) | `scripts/covariates/build_hces_diet_block.R` | 15 `hces_` columns in the shared set: food share and own-production share for all four countries, consumption per adult equivalent for three, HDDS and 11 food-group indicators for Malawi and The Gambia (T1.1 **done**, T3.1 superseded; the household tables are `data/covariates/harmonized/hces_household_<Country>.csv`) |
| `RTFP market food price.zip` | `scripts/covariates/build_rtfp_price_block.R` | 7 `rtfp_` columns (local price level, 12-month inflation, volatility, seasonal range, staple price, distance to market) for The Gambia (28 markets) and Malawi (129 markets); Ghana and Sierra Leone are not in the RTFP panel |
| `DataBank_Food Prices for Nutrition.csv` | `scripts/build_fpn_affordability.R` | `data/FPN/<Country>_fpn_admin2.csv` (T1.4 **done**; national constants, series years mixed 2017/2021 because the cost-share and headcount series start in 2021) |
| `GFDx.zip` (compendium + data dictionary) | `scripts/protocol_v2/59_build_addback_sources.R` | nutrient codes corrected (8 = niacin, 15 = zinc); T1.5 **done** |
| `ACLED data-Western_Africa-Ghana_2015_2017.csv` | `scripts/protocol_v2/60_build_acled_conflict.R` | Ghana 2015-17 only; the script waits for the four-country export (T1.3 **open**) |
| `_Cleaned datasets (1).zip` | - | the same `combined_dataset.dta` the national track already reads |
| `FPMA_3countries (2).csv` | - | 2025 only (24 months) and no Ghana; needs the full-history export (T1.7 below) |
| `GDD.zip`, `HANCI.zip`, `Cadre Harmonise.zip`, `HFID.zip`, `FAOSTAT.zip`, `FAO.zip`, `GHED_data.XLSX`, `MPI district xlsx`, `STATcompiler.zip`, `VMNIS.zip`, GAIN FSD csv, GDL csv | not yet wired | GDL and VMNIS duplicate what the pipeline already fetches; the rest are national or single-country (AB-01, "sources that don't make sense to harmonise") |

Still open in Tier 1, in priority order: **T1.2 (MICS GPS for Malawi 2013-14 and Ghana 2017-18)**,
**T1.3 (ACLED, all four countries)**, **T1.6 (SLIHS 2018 diary codes, GLSS7
food module, SLIHS 2011)**, **T1.7 (FPMA full history)**.

---

## Tier 1 — downloads that unblock a domain (registration, then a script exists)

### T1.1 HCES microdata - DONE 2026-09-15

The three surveys arrived in `LSMS.zip` and the block is built (see Status).
What the block still lacks is in T1.6.

### T1.2 MICS microdata and GPS - what is needed, survey by survey

MICS is the one survey programme besides DHS that measures, on the same
households, the things the biomarkers respond to: **household salt iodine
test**, **infant and young child feeding food groups** (eggs, flesh foods,
dairy, vitamin-A-rich fruit and vegetables; minimum dietary diversity),
vitamin A supplementation and deworming, iron-folic acid in pregnancy,
child anthropometry, water quality (E. coli, MICS6), WASH and wealth. It is
not nested in the micronutrient surveys (unlike the Malawi MDHS), so nothing
has to be excluded for leakage. The GPS files are requested separately on
mics.unicef.org ("GIS datasets", a short justification form), the route you
already used for Sierra Leone 2017 and The Gambia 2018.

| Survey | Have | Need | Why |
|---|---|---|---|
| The Gambia MICS6 2018 | microdata (`data/MICS/Gambia 2018/`) and GPS (`data/RA_2026-09/extracted/MICS GPS/GMB2018/`) | nothing | complete; Andrew runs the cluster model on it |
| Sierra Leone MICS6 2017 | microdata (`data/MICS/Sierra Leone MICS6 Datasets/`: hh with the salt test SA1-SA2 and the water-quality module WQ, hl, wm, ch, bh, fs, fg, mn, pn, tn) and GPS (`MICS GPS/SL2017/`, 600 of 600 clusters match `HH1`) | nothing | complete (checked 2026-09-15) |
| Malawi MICS5 2013-14 | microdata (`data/MICS/Malawi MICS 2013-14 SPSS Datasets/`: hh with salt test SI1, hl, wm, ch, bh, mm, mn, tn; `HH7` = 31 districts and cities, 1,139 clusters) | **GPS request** (mics.unicef.org, GIS datasets; say whether it exists for this round) | the only independent household survey near the MNS 2015-16 (the MDHS is nested in the MNS); without GPS the block is broadcast from the 31 districts to the TAs, with GPS it goes to TA level like the DHS block |
| Ghana MICS6 2017-18 | microdata (`data/MICS/Ghana 2017/`; `HH7` = 10 regions, 660 clusters) | **GPS request** | same year as the GMS 2017; without GPS it stays at 10 regions broadcast to 260 districts |
| Sierra Leone MICS4 2010 | - | optional: hh, hl, wm, ch (no GPS for MICS4) | the other side of the 2013 MNS in time; lower priority now that DHS 2013 and MICS 2017 bracket it |

Save any new files under `data/MICS/<Country> <year>/` with the survey's
questionnaire and report. **Done when:** the two GPS requests are filed
(request date and outcome in the provenance file) and each MICS folder has
a `README_download.md` listing the files and the region/district variable
(`HH7` / `HH7A`).

### T1.3 ACLED export (conflict domain; the script is ready) - PARTIAL

Only `ACLED data-Western_Africa-Ghana_2015_2017.csv` arrived (Ghana, three
years). The script needs all four countries from 2007-01-01, so the export
must be redone. Create an ACLED account (https://acleddata.com, research access), open the
Data Export Tool, select **Gambia, Ghana, Malawi, Sierra Leone** (add Nigeria,
Burkina Faso, Mali, Niger, Ethiopia if the account allows — those are the next
countries), all event types, dates 2007-01-01 to today, and export CSV. Save
the file(s) under `data/ACLED/` (any name, `*.csv`). **Done when:**
`Rscript -e "source('scripts/protocol_v2/60_build_acled_conflict.R')"` runs and
prints one line per country with the event counts; paste that output into the
provenance row. Licence: ACLED data must not be redistributed.

### T1.4 Food Prices for Nutrition - DONE 2026-09-15

`scripts/build_fpn_affordability.R` now reads the series-wide DataBank
layout you exported and wrote the five country files. Nothing more to do.

### T1.5 GFDx codebook check - DONE 2026-09-15

The compendium in `GFDx.zip` settled it: 8 = niacin, 15 = zinc; the `NUT`
map in `59_build_addback_sources.R` was corrected and the block rebuilt.

### T1.6 The pieces the HCES block still lacks

1. **SLIHS 2018 diary item codes.** The released `slihs2018_x.dta` diary has
   item codes without labels, so Sierra Leone has no food-group indicators.
   Get the item code list: the SLIHS 2018 questionnaire Book 4 (household
   consumption diary) or the codes annex from Statistics Sierra Leone / the
   World Bank Microdata Library entry (`SLE_2018_IHS`). Save it as
   `data/RA_2026-09/extracted/LSMS/SLE_2018/slihs2018_item_codes.csv`
   (`code, label`). **Done when:** every code in `slihs2018_x.dta` has a label.
2. **GLSS7 food consumption module.** Only the expenditure aggregates
   (`padq_hh_R` etc.) were used for Ghana. Download or locate the GLSS7
   Section 9 Part B food consumption file (item level, 7-day recall, with
   own-production quantities) from the GSS / World Bank entry and file it
   under `data/RA_2026-09/extracted/LSMS/GHA_2017/`. **Done when:** the file
   is there with its item code list.
3. **SLIHS 2011.** The 2013 Sierra Leone MNS sits between SLIHS 2011 and
   SLIHS 2018; 2011 is the closer round for diet. World Bank Microdata
   Library `SLE_2011_IHS` (public use after registration). Save under
   `data/RA_2026-09/extracted/LSMS/SLE_2011/`. **Done when:** the
   consumption module and the cluster/district file are identified by name.

### T1.7 FAO FPMA full-history export (extends the price block to four countries)

RTFP covers only The Gambia and Malawi. The FAO Food Price Monitoring and
Analysis tool has Sierra Leone (Bo, Freetown, Kailahun, Kenema, Koinadugu:
rice, imported rice, cassava, palm oil) and Ghana (Accra, Kumasi, Tamale,
Techiman, Bolgatanga, ...: maize, rice, cassava, yam, plantain, ...) retail
series back to the 2000s. The export in the folder covers 2025 only and
three countries. Re-export from https://fpma.fao.org/giews/fpmat4/ for
**all four countries, all markets, all commodities, monthly, 2007-01 to the
present**, in local currency as well as USD/kg, one CSV per country, into
`data/RA_2026-09/FPMA/`. **Done when:** each file covers the survey year for
its country (2013 Sierra Leone, 2016 Malawi, 2017 Ghana, 2018 The Gambia)
and the market names are listed with their districts in a `markets.csv`
(market, district, GADM Admin-2 name) so they can be geolocated.

---

## Tier 2 — research and transcription (no account needed; PDFs and web)

### T2.0 National estimates missing from VMNIS, from the four reports already in hand (do this first, ~1 day)

`data/national/vmnis_national_rep.rds` (WHO VMNIS) holds **no Gambia 2018
rows at all**, no vitamin A for Ghana 2017 or Sierra Leone 2013, and no iron
for any survey (VMNIS has no iron panel). The four survey reports and the
secondary publications collected under `Downloads/Micronutrient survey
reports/` contain those national estimates. Transcribe them into
`metadata/vmnis_supplement_reports.csv` with the VMNIS column names
(`iso3c, country, year, Beginyear, Endyear, Population, mn_group, Indicator,
Deficiencycutoff, Dataadjustedfor, Samplesize, Prevalenceofdeficiency, Mean,
L95CI, U95CI, source_document, source_table`) — one row per survey ×
population × nutrient (vitamin A, iron/ferritin, folate, B12, zinc, iodine,
selenium where measured), stating the cut-off and inflammation adjustment the
report used. **Done when:** GMNS 2018, GMS 2017, SLMS 2013 and MNS 2015-16
each have their published national rows and Andrew can merge the file into
`vmnis_national()`.

### T2.1 Registry of standalone national micronutrient surveys, 18 West African countries + panel candidates

This serves two things at once: the supplementary table the data-landscape
paper lacks, and the panel-expansion decision. Start from
`mn-proxies/Landscape data sources/Landscape_data source_Nutrition Surveys_2025-01-27.xlsx`
and `docs/micronutrient_survey_candidates.md`; extend by web search (survey
reports, VMNIS survey list, BRINDA membership, GAIN/UNICEF/CDC publications).
One row per survey with: country, survey name, fieldwork dates, populations
(PSC / SAC / WRA / pregnant / men), biomarkers measured (ferritin, sTfR, RBP,
retinol, zinc, folate, B12, urinary iodine, CRP/AGP), **stratification level
and number of strata** (national only / zones / regions / districts),
whether cluster GPS were collected, report URL, microdata access route
(public / request / BRINDA / none known), contact, and notes. Template:
`metadata/mn surveys/mn_survey_registry_TEMPLATE.csv` (create it with these
columns). **Done when:** every one of the 18 countries has at least a "none
found" row, and the ten candidates in `micronutrient_survey_candidates.md` are
included with their stratification level filled in.

### T2.2 Report-derived regional outcomes for candidate training countries

For surveys with a public report but no microdata in hand (Nigeria NFCMS 2021,
Cameroon 2009, Côte d'Ivoire 2007, Liberia 2011, Burkina Faso 2010 and 2014,
Senegal 2010, Kenya 2011, Ethiopia 2015), transcribe the **regional** tables:
one row per survey × region × population × nutrient with prevalence, n,
95% CI, mean or median biomarker where printed, the cut-off and the
inflammation adjustment the report used, and the table/page reference.
Template `data/outcomes_reports/TEMPLATE_regional_outcomes.csv` with columns:
`iso3, survey, year, region_label_report, gadm_admin1_match, population,
nutrient, indicator, cutoff, adjustment, n, prevalence, ci_low, ci_high, mean,
median, source_table, source_page`. Match `gadm_admin1_match` to GADM 4.1
level-1 names; where a report zone spans several GADM regions, list them
separated by `;`. **Done when:** each transcribed survey has a companion
`README` stating how the region labels were matched and any zone that could
not be matched. These rows let the regional transport test (pre-registration
P1) run for a new country without its microdata.

### T2.3 SMART survey district tables

Sierra Leone national SMART surveys (2010, 2014, 2017, 2021) and The Gambia
(2012, 2015; check for 2018) publish district-level wasting, stunting, MUAC and
retrospective mortality. Collect the report PDFs under `data/SMART/reports/`
and transcribe the district tables into
`data/SMART/smart_district_tables.csv` (`iso3, survey_year, district_report,
gadm_admin2_match, indicator, estimate, ci_low, ci_high, n, source_page`).
Check whether Malawi or Ghana have comparable national SMART rounds near
2016 / 2017 (Ghana had northern-region SMART surveys only). **Done when:** the
CSV covers the rounds nearest each survey year (SL 2014, Gambia 2015) and the
match rate to GADM districts is reported.

### T2.4 WHO Health Inequality Data Repository: other datasets with MICS rows

`scripts/protocol_v2/59a_fetch_who_heat.py` pulled the RMNCH, immunisation and
malaria datasets; only the immunisation one carried MICS rows for our
countries. Check the repository's other survey-based datasets (child
malnutrition, WASH, adolescent health, NCD risk factors, any "MICS" or
"UNICEF" labelled dataset) for subnational-region rows for GMB / GHA / MWI /
SLE and, if found, note the dataset ids so the fetch script can be extended.
**Done when:** a short table of dataset id → countries → rounds → number of
indicators with `Subnational region` disaggregation is added to the
provenance notes.

### T2.5 VMNIS subnational entries

In the WHO VMNIS extract (`mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta`,
fields `Representativeness`, `Areacovered`) list every subnational entry for
the 18 countries and the four panel countries: survey, year, area, population,
nutrient, prevalence, n. **Done when:** `results/tables/vmnis_subnational_entries.csv`
exists and a paragraph says whether VMNIS can supply regional outcomes for any
survey we lack a report for.

### T2.6 MAPS estimates for Malawi

Check whether the MAPS tool (Micronutrient Action Policy Support; LSHTM /
BMGF) exposes downloadable subnational nutrient-inadequacy estimates for
Malawi (IHS4-based) and for any other panel or candidate country, at what
level, and under what licence; download to `data/MAPS/` if open. **Done when:**
a note records coverage and access, and any downloaded file has a provenance row.

### T2.7 GAIN FACT Ghana

Write to GAIN (fortification assessment team) asking for the Ghana FACT survey
regional coverage tables (fortified wheat flour and oil, household and
individual coverage). Ghana-only, so an in-fill block; low priority but one
email. **Done when:** the reply and any tables are filed under `data/GAIN_FACT/`.

---

## Tier 3 — cleaning to an existing pattern (after the Tier 1 downloads land)

### T3.1 HCES first-pass clean - SUPERSEDED

Done in `scripts/covariates/build_hces_diet_block.R`; the household files
are `data/covariates/harmonized/hces_household_<Country>.csv` and the item
classifications `metadata/hces_food_groups_{Malawi,Gambia}.csv`. A useful
review task instead: read the two food-group files line by line and flag any
item whose HDDS group looks wrong (the classification is keyword-based).

### T3.2 MICS region summaries for Sierra Leone and Malawi - SUPERSEDED by the MICS block (2026-09-15)

Run the existing MICS clean pattern (`data/IPD/Malawi/malawi_MICS_clean.R`,
which currently processes the Gambia 2018 files; `metadata/MICS/` holds the
variable metadata) on the T1.2 downloads to produce
`data/MICS/mics_<country>_<year>_region_summary.csv` with the same layout as
the Gambia and Ghana files. **Done when:** the four countries' summaries share
a documented indicator list (`metadata/MICS/mics_common_indicators.csv`) and
Malawi 2013-14 is summarised at district level.

### T3.3 Annotate the new predictor columns

`results/tables/protocol_v2/variable_sheet.csv` is built by
`scripts/protocol_v2/06_build_variable_sheet.R` and, per
`docs/RA_brief_predictor_annotation.md`, still had 79 rows without a
definition before today; the rebuild added 77 columns (17 restored DHS, 60
add-back). Re-run the sheet builder and fill `definition`, `unit`, `direction`
(more = better / worse for nutrition) and `sub_domain` for every row marked
`MISSING_needs_RA`, using `predictors_admin2_shared_metadata.csv`,
`predictors_admin2_addback_metadata.csv` (its `assumption` column) and the
provider documentation. **Done when:** zero `MISSING_needs_RA` rows.

### T3.4 Region crosswalks for new surveys

For every survey transcribed in T2.2 and every HCES in T1.1, a crosswalk CSV
`metadata/crosswalks/<iso3>_<survey>_regions_to_gadm.csv` (report label →
GADM level and name, with a note when boundaries changed between the survey
and GADM 4.1). Ghana's `data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv`
is the model. **Done when:** every transcribed region resolves to a GADM unit
or is explicitly marked unmatched.

---

## Tier 4 — for the data-landscape paper revision

### T4.1 Catalogue rows for products the review does not list

Draft Table 2 rows (same columns: managing source, source function, DVC stage,
indicator domains, data level, temporal resolution, geographic coverage,
access status, URL, inputs, subnational availability, API/bulk access) for the
products the proxy work uses that the paper omits: AlphaEarth Foundations,
iSDAsoil, MapSPAM, Gridded Livestock of the World, Köppen-Geiger climate
zones, HarvestChoice AEZ, Meta Relative Wealth Index, JRC GHSL, JRC Global
Surface Water, World Settlement Footprint, CCNL night-time lights, Global
Human Modification, Global Pasture Watch, Global Data Lab SHDI. **Done when:**
a `docs/paper_supplement_extra_products.csv` in the paper's Table 2 format is
ready for Andrew to review.

### T4.2 The standalone-survey supplementary table

T2.1's registry, filtered to the 18 countries and formatted to the paper's
Table 1 style, with the stratification level and access route columns kept.

---

## Not RA tasks (need Andrew or a decision)

Re-running the models on the 552-column set; deciding whether Tang et al. may
be used as a predictor; the Gambia legacy DAG rebuild; anything that changes
`R/`, `_targets.R` or the protocol scripts.
