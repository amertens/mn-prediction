# LK-02 / AB-01 — leakage policy corrected, sources added back, what still needs a download

15 September 2026. Scope: the shared Admin-2 predictor set
(`data/covariates/harmonized/predictors_admin2_shared.csv`) and the national
VMNIS track. **No model was re-run**; the set is rebuilt to the point where a
predictor audit can sign it off (`results/tables/predictor_audit_2026-09-15*`).

## 1. The leakage rule, restated

A predictor is leakage only when it is measured on the same individuals as the
outcome, i.e. it comes from the same survey instance. Everything from an
independent source is admissible whatever it measures. Consequences:

| Layer | Before | After |
|---|---|---|
| `metadata/covariates/exclusions.csv` | 7 `leakage` rows excluding any `dhs_` column naming anaemia / Hb / iron / vitamin A / zinc / folate / B12 / iodine | rows removed (backup `exclusions.csv.pre_LK02`); the 6 `data_defect` rows stay |
| `dhs_leakage_review.csv` | 4 patterns "exclude" | "include", with the LK-02 reasoning |
| `drop_near_outcome_v2()` (fit time) | dropped the same columns | drops nothing by default; `V2_DROP_MODELLED=1` still removes the modelled outcome-adjacent surfaces as a sensitivity |
| `gw_` guard (`R/data_prep.R`) | blocks the biomarker survey's own blood-draw / Hb columns | **unchanged** — this is the real guard |
| Malawi DHS block | all 850 MDHS clusters, including the 105 the MNS re-sampled | `metadata/mns_dhs_overlap_clusters.csv` + `R/mns_dhs_overlap.R`; every DHS-derived Malawi predictor is built from the other 745 clusters (surveyPrev aggregation, custom indicators, cluster-model BYM2, builder direct block) |
| National panel (`R/national_covariates.R`) | all measured anaemia/Hb columns dropped | columns kept; cells of DHS-linked VMNIS country-years blanked before the nearest-year carry (`dhs_linked_vmnis_surveys()`: methodology text + `metadata/vmnis_dhs_linked_surveys.csv`); new cache `panel_national_cov50_lk02.rds` |

Back in the DHS block (per country, where the round carries the item): women's
and children's anaemia prevalence (any / moderate), mean child haemoglobin, iron
in pregnancy (any / days / 90+), postpartum vitamin A, child vitamin A
supplement / capsule, iron- and vitamin A-rich foods, vitamin A-rich fruit and
vegetables food group, zinc for diarrhoea, surveyPrev `AN_ANEM_W_ANY` and
`CN_ANMC_C_ANY`. Nothing outside the DHS block had been dropped by the rules.

Tests: `tests/testthat/test-mns-dhs-overlap.R` (18 assertions),
`tests/testthat/test-national-same-survey.R` (14).

## 2. Added to the harmonized set

Built by `scripts/protocol_v2/59_build_addback_sources.R` (order: builder →
07 → 08 → **59** → 53). Nothing is filtered for coverage; the consumer reads
`n_countries` and `subnational` in the metadata.

**Result of the rebuild (2026-09-15, 09:59; re-run of 59/53 at 10:20 after the
audit):** `predictors_admin2_shared.csv` = 554 Admin-2 rows × **530**
predictors (was 454): +17 restored DHS columns, +60 add-back columns, −1
all-NA column (`dhs_w_barrier_transport`, now a `data_defect` rule). 482
columns present in all four countries, 48 partial, 45 national constants, 0
value-identical duplicates, 161 DHS-derived, 21 modelled outcome-adjacent
surfaces. DS-01 cluster-model BYM2 estimates now stand in 130 DHS columns
(Malawi's fitted on the 745 non-MNS clusters: median Spearman 0.967 against
the old 850-cluster columns over 143 columns, 2 % below 0.9). The audit caught
one defect in the first pass of script 59 — `match_names()` kept the source
spelling on exact hits, which left the HEAT block populated for Sierra Leone
only — fixed and re-run; `results/tables/predictor_audit_2026-09-15*` and
`results/tables/source_inventory_2026-09-15*` describe the final set.

| Block | Source (landscape paper table) | Columns | Countries | Level | Alignment |
|---|---|---:|---|---|---|
| A | UNICEF vitamin A supplementation coverage (T2, UNICEF Global Databases) | 1 | 4 | national | nearest year |
| B | WHO FluNet (T1) | 3 | Ghana, Sierra Leone | national | survey year |
| C | GFDx **programme** fields: mandatory legislation in force, years in force, vehicle intake, share industrially processed, potential nutrient delivery for iron / vitamin A / folic acid / B12 / zinc / iodine (T2) | 21 | 4 | national | survey year |
| D | WHO/UNICEF Global Anaemia Estimates: children, non-pregnant, pregnant women + 5-year trend (T2) | 5 | 4 | national | survey year |
| E | Tang et al. 2026 (WFP / MIMI) HCES-based nutrient inadequacy (T2, MIMI) | 6 | Ghana | Admin-1 (10 old regions) | GLSS 2016/17 |
| F | Global Data Lab subnational HDI + 5-year trend (*not in the paper*) | 2 | 4 | Admin-1 (Sierra Leone: Admin-2) | nearest year |
| G | MICS immunisation indicators by subnational region, WHO Health Inequality Data Repository (T1 MICS via T2 WHO repository) | 20 (18 in all four; rotavirus in three) | 4 | Admin-1; Malawi districts | nearest MICS round (2018 / 2017 / 2010 / 2014) |

Block G is the first MICS block that covers Sierra Leone and Malawi: the HEAT
export (`scripts/protocol_v2/59a_fetch_who_heat.py`) carries every MICS round
at region level, but only its immunisation dataset holds MICS rows for these
countries, so the full MICS harmonisation still needs microdata (§4).

Ready but waiting on a download: **ACLED** conflict exposure
(`scripts/protocol_v2/60_build_acled_conflict.R`, 8 columns per district over
the 36 months before fieldwork; tested on a synthetic export).

Two defects found on the way, fixed here:
- the GFDx block that script 08 took was WHO 2011 anaemia and Wessells & Brown
  2012 zinc prevalence, not fortification programme data (those four columns
  stay, relabelled as national modelled surfaces; the programme block is new);
- the legacy Gambia merge (`src/Gambia/2_GW_Gambia_data_merge.R`, the
  `Admin1_old` block) swapped the two Central River LGAs (Niaminas / Fulladu
  West sent to Kuntaur, Saloums / Nianija / Sami to Janjanbureh) and folded
  Kanifing into Banjul. Measured in `Gambia_merged_dataset.rds`: the 20
  `dhs2019_` Admin-1 columns are NA for every division except Banjul anyway,
  because `load_dhs_admin1()` keys on GADM division names while `Admin1_old`
  holds LGA names, so that join never matched; the 109 `_adm2` columns (joined
  on Admin2) are fine. The mapping is corrected in the script (takes effect at
  the next merge re-run) and `gambia_lga()` in script 59 carries the same
  correction. The shared Admin-2 set is unaffected (GPS point-in-polygon).

## 3. Sources that should NOT be harmonised (and why)

| Source (paper) | Reason |
|---|---|
| FRAT, GAIN FACT (except Ghana), PMA, VACS, WEAI | none of the four countries at the survey years, or project areas only |
| WHO STEPS, WHO NCD microdata | adult NCD risk factors, national, no Admin-2 signal |
| Malaria Indicator Surveys | malaria already at 5 km from the Malaria Atlas at the survey year |
| JMP, UNAIDS, WUENIC, JME | national versions of what IHME surfaces and the DHS block already give subnationally |
| GBD, WHO GHE, WHO Mortality DB, IHME SDG, UN SDG DB, GHI, HANCI, GII, GDI, ILOSTAT, GHED, PIP, HDR | national indices or accounts: can move a country level, never a district ranking; the WDI block already carries the ones with a causal story (national track) |
| Climate Engine, Copernicus CDS, NASA EarthData, GLOSIS | the same rasters are already extracted directly (TerraClimate, CHIRPS/TRMM, MODIS, SoilGrids/iSDA) |
| CHIRPS national aggregates | tested as a time-matched gridded add-on (script 41): no gain |
| WFP HungerMap / VAM portal, IPC portal, Food Systems Dashboard, OWID, GHO, NLiS, WHO Nutrition Data Portal, GNR, HDX, GHDx | repositories / dashboards, not data sources; the underlying series are taken directly |
| Cadre Harmonisé standalone block | already harmonised through HFID (which merges IPC / CH / FEWS NET); the country-level `fsec_ch_*` block is a December-2025 snapshot, time-mismatched with 2013–2018 surveys |
| GEMS/Food | national contaminant monitoring, coverage uneven, no subnational detail |
| GDELT | open but the daily event files for a 3-year window are gigabytes per country and machine-coded; ACLED is the curated equivalent |
| Gallup GDQP, GDD | national diet quality by strata; low value at Admin-2 and both need a data request |
| FEWS NET beyond HFID | already in (IPC phase via HFID; NDVI anomaly via the GEE block) |
| LSMS Ghana-only block (134 region means) | in-fill for one country only and it lost to the base (AD-lsms); harmonise HCES across the four instead (§4) |

## 4. Tasks that need you (a download, a login, or a decision)

Ordered by what the sandbox says would move the results.

1. **HCES microdata for the three countries without it** — the dietary
   domain the ablations say is missing. Register / log in at the World Bank
   Microdata Library and download: Malawi **IHS4 2016-17** (LSMS-ISA,
   `MWI_2016_IHS-IV_v04_M`), Sierra Leone **SLIHS 2018** (Stats SL / WB),
   The Gambia **IHS 2015/16** (GBoS / WB). Drop the zips under
   `data/LSMS/<Country>/`. I then build a harmonised block at each survey's
   design level (Malawi & SL districts, Gambia LGAs, Ghana regions): food
   expenditure share, animal-source and pulse shares, dietary energy per adult
   equivalent, purchase vs own production, and MIMI-style nutrient adequacy.
   Alternative for two countries without microdata: the MAPS tool estimates
   for Malawi (IHS4-based) and Tang et al. for Ghana (already in, block E).
2. **MICS microdata for the full MICS block** — Sierra Leone MICS 2010 (and
   2017) and Malawi MICS 2013-14 from mics.unicef.org (registration). Then the
   existing MICS clean scripts (`data/IPD/Malawi/malawi_MICS_clean.R`,
   `src/Gambia/...`) are run for all four and a curated ~30-indicator set is
   harmonised at region / district level.
3. **ACLED export** — acleddata.com Data Export Tool, four countries, all
   event types, 2007 to date, saved under `data/ACLED/*.csv`; then
   `Rscript -e "source('scripts/protocol_v2/60_build_acled_conflict.R')"`.
   Little variance in these four; needed before Nigeria / Burkina / Mali.
4. **Food Prices for Nutrition (CoAHD / CoNA)** — one manual DataBank CSV
   export into `data/FPN/`; `scripts/build_fpn_affordability.R` is ready.
   National only.
5. **GFDx codebook check** — confirm the `standard_nutrient` codes I read
   as 7 iron / 12 vitamin A / 5 folic acid / 2 B12 / 8 zinc / 6 iodine
   (inferred from the target levels against the published standards).
6. **SMART survey reports** — Sierra Leone national SMART 2010 / 2014 / 2017
   and Gambia 2012 / 2015 give district wasting and stunting at or near the
   survey years; an RA can transcribe the district tables (PDFs are public,
   microdata is not centralised).
7. **GAIN FACT Ghana** — request the survey (regional coverage of fortified
   wheat flour and oil); Ghana-only, so an in-fill block.
8. **GDD / Gallup GDQP** — data requests; national only, lowest priority.
9. **Decision:** whether Tang et al. (block E) may be used as a predictor at
   all, since it was the Ghana external check
   (`docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md`); the Côte d'Ivoire check
   is unaffected either way. It is in the set, flagged MODELLED SURFACE.
10. **Decision:** fix the Central River LGA swap in the legacy Gambia merge
    (§2) and rebuild `Gambia_merged_dataset.rds`, or leave the legacy country
    dataset as it is.

## 5. Can published subnational estimates serve as outcomes for more countries?

For the tier where the signal is, yes. WS4a found Admin-1 the best resolution
(9 of 14 cells; beats Admin-2 in 13 of 14), transport works at the regional
tier (mean Spearman ~0.31–0.45) and is noise-limited at district level, and
each added training country buys ~0.05. A survey report's regional table
gives what the regional model needs: prevalence per region, n (for effective
n via the design effect, or from the published CI: n_eff = p(1−p)/SE²), and
sometimes a mean or median biomarker. VMNIS itself carries some subnational
rows. What you lose without microdata: the continuous target (+0.10 in every
estimand, unavailable where only prevalence is published), the uniform
BRINDA / cut-off re-derivation (rank-normalisation within country removes the
level offset, but not adjustment-induced re-ordering), district targets, the
cluster track, and the individual-level sensitivity arms. Two practical
limits: reports stratify by 3–6 zones in several countries (Ghana 3 belts,
Malawi 3 regions, Sierra Leone 4), below the ≥ 8 units the regional protocol
needs; and region definitions must be crosswalked to GADM Admin-1. So:
report-derived regional outcomes are enough to (a) score the pre-registered
P1 predictions for a new country and (b) add training regions, but not to
replace microdata for the district product. The 18-country landscape supplement
(§6) should record, per standalone micronutrient survey, the stratification
level and whether microdata are obtainable, because that decides which use a
survey can serve.

## 6. For the landscape paper

- The four micronutrient surveys this project depends on fall in the one
  category the paper did not catalogue (standalone national micronutrient
  surveys). A supplementary table of known standalone surveys in the 18
  countries, with stratification level and access route, fixes the paper's
  largest gap relative to the modelling it was written to support; the
  material exists in `mn-proxies/Landscape data sources/Landscape_data
  source_Nutrition Surveys_2025-01-27.xlsx` and
  `docs/micronutrient_survey_candidates.md`.
- Ten products the proxy set relies on are absent from the paper: AlphaEarth,
  iSDAsoil, MapSPAM, GLW4, Köppen-Geiger, HarvestChoice AEZ, Meta RWI, JRC GHSL,
  JRC GSW, LSIB — plus GDL SHDI and the WHO Health Inequality Data Repository's
  subnational DHS/MICS extracts, which turned out to be the fastest route to a
  cross-country MICS block.
