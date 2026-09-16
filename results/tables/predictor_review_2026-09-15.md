# Predictor set review, 2026-09-15

554 units x 570 predictors. Tiers: open 356, survey_dhs 150, survey_public 64.

## 1. Design matrix by arm (per-country prep: coverage >= 0.7 and sd > 0; transport = 4-country intersection)

| tiers | declared (after policy) | Gambia | Ghana | Malawi | Sierra Leone | transport common | domains in common |
|---|---:|---:|---:|---:|---:|---:|---:|
| open | 311 | 298 | 299 | 297 | 291 | 280 | 18 |
| open+survey_public | 375 | 360 | 350 | 359 | 340 | 326 | 20 |
| open+survey_public+survey_dhs | 525 | 509 | 490 | 508 | 488 | 464 | 24 |

Transport headline arm (open + survey_public): columns in the common matrix by domain

| domain | columns |
|---|---:|
| Satellite embedding | 64 |
| Climate and weather | 56 |
| Soil characteristics | 45 |
| Infection and inflammation burden | 25 |
| Immunisation | 18 |
| Water and sanitation | 15 |
| Ecosystem productivity/greenness | 14 |
| Education, employment, SES | 13 |
| Agricultural production, land use | 12 |
| Built environment | 10 |
| Child anthropometry | 8 |
| Infant and young child feeding | 8 |
| Livestock density | 8 |
| Household assets and characteristics | 7 |
| Ruralness, population density, built environment | 7 |
| Anaemia and haemoglobin | 4 |
| Food prices and supply | 4 |
| Water and coast proximity | 4 |
| Food fortification and supplementation | 2 |
| Household diet and consumption (HCES) | 2 |

Columns of the headline arm that do NOT reach the transport matrix (missing or constant in the country named)

| source | absent in | columns |
|---|---|---:|
| HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Ghana;SierraLeone | 12 |
| Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Gambia;Malawi;SierraLeone | 6 |
| World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Ghana;SierraLeone | 6 |
| WFP / HDX market prices | Malawi | 3 |
| WHO ESPEN implementation-unit database 2014-2025 (portal export, no key) | Malawi | 3 |
| HFID (FEWS NET / IPC / WFP mVAM) | Ghana;SierraLeone | 2 |
| MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | Malawi | 2 |
| WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region | SierraLeone | 2 |
| HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia | 1 |
| HFID (FEWS NET / IPC / WFP mVAM) | Gambia;Ghana | 1 |
| HFID (FEWS NET / IPC / WFP mVAM) | Ghana | 1 |
| IHME (modelled surfaces) | Gambia | 1 |
| JRC Global Surface Water 1.4; USDOS LSIB 2017 (Earth Engine) | Gambia | 1 |
| JRC Global Surface Water 1.4; USDOS LSIB 2017 (Earth Engine) | Gambia;SierraLeone | 1 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | Gambia;Ghana;SierraLeone | 1 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | Gambia;Malawi;SierraLeone | 1 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | Gambia;SierraLeone | 1 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | SierraLeone | 1 |
| MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | Gambia;Ghana;SierraLeone | 1 |
| Malaria Atlas Project (latest release, layer at each survey year) | SierraLeone | 1 |
| WFP / HDX market prices | Ghana;Malawi | 1 |

## 2. Value sanity

- share-like columns outside [0, 1]: 71 (dhs_c_mean_birthweight, dhs_c_mean_hemoglobin, dhs_c_mean_waz, dhs_c_mean_whz, dhs_hh_wealth_mean, dhs_w_anc_visits, dhs_w_mean_bmi, grassland_frac, lcover_bare_frac_t0, lcover_crops_frac_t0, lcover_grass_frac_t0, lcover_shrub_frac_t0)
- near-constant subnational columns (CV < 1%): 0 
- |skew| > 5 (rank-normalised before any fit, so cosmetic): 17
- columns with partial missingness inside a country (5-95% of units): 58 (dhs_c_diet_diversity_score, dhs_c_mean_birthweight, dhs_c_mean_haz, dhs_c_mean_hemoglobin, dhs_c_mean_waz, dhs_c_mean_whz, dhs_hh_mean_hhsize, dhs_hh_wealth_mean, dhs_w_anc_visits, dhs_w_birth_interval_short, dhs_w_edu_years_mean, dhs_w_mean_age_first_birth, dhs_w_mean_bmi, dhs_w_mean_parity, spam_parea_total)

## 3. Redundancy (pooled within-country rank-normalised |r| >= 0.98)

- 0 pairs, 0 across sources; by source pair:

- none


## 4. Cross-source agreement (within-country Spearman; NA = fewer than 8 units with both)

| a | b | Gambia | Ghana | Malawi | Sierra Leone |
|---|---|---:|---:|---:|---:|
| `mics_water_improved` | `dhs_hh_improved_water` | 0.06 | -0.07 | -0.11 | 0.59 |
| `mics_sanitation_improved` | `dhs_hh_improved_sanitation` | 0.54 | 0.26 | 0.12 | 0.43 |
| `mics_open_defecation` | `dhs_hh_open_defecation` | 0.46 | 0.75 | 0.36 | 0.52 |
| `mics_c_stunted` | `dhs_c_stunted` | 0.78 | 0.51 | 0.31 | 0.31 |
| `mics_c_stunted` | `ihme_stuntingprevalence` | 0.62 | 0.63 | 0.4 | 0.28 |
| `mics_wealth_score_mean` | `dhs_hh_wealth_mean` | 0.73 | 0.68 | 0.33 | 0.89 |
| `mics_wealth_score_mean` | `rwi_mean` | 0.62 | 0.57 | -0.27 | 0.59 |
| `mics_w_secondary_plus` | `dhs_w_secondary_plus` | 0.93 | 0.72 | 0.57 | 0.89 |
| `mics_w_literate` | `dhs_w_literate` | 0.55 | 0.47 | 0.65 | 0.35 |
| `mics_c_diarrhoea_2wk` | `dhs_c_diarrhea_2wk` | 0.08 | 0.43 | 0.36 | -0.3 |
| `mics_c_fever_2wk` | `dhs_c_fever_2wk` | 0.15 | 0.2 | 0.34 | -0.09 |
| `mics_c_mdd` | `dhs_c_mdd_4plus` | 0 | 0.35 | -0.09 | -0.36 |
| `mics_electricity` | `dhs_hh_electricity` | 0.48 | 0.68 | NA | 0.69 |
| `mics_heat_vmsl_sy` | `dhs_c_measles1` | 0.51 | 0.17 | -0.07 | 0.1 |
| `ndvi_modis_t0` | `evi_t0` | 0.88 | 0.91 | 0.77 | 0.8 |
| `ndvi_modis_win` | `ndvi_modis_t0` | 0.95 | 0.99 | 0.95 | 1 |
| `hces_food_share` | `dhs_hh_wealth_mean` | -0.68 | -0.35 | -0.27 | -0.45 |
| `hces_log_cons_pae_rel` | `rwi_mean` | NA | 0.6 | 0.07 | 0.72 |
| `hces_hdds` | `dhs_c_diet_diversity_score` | -0.02 | NA | 0.14 | NA |
| `map_sy_pf_parasite_rate` | `ihme_malaria_pfpr` | -0.72 | 0.78 | 0.85 | 0.64 |
| `wpop_log_density_survey_year` | `ghs_pop` | 0.97 | 0.95 | 1 | 1 |
| `rtfp_fpi_rel_national` | `fprice_staple_rel` | 0.69 | NA | 0.26 | NA |

## 5. Broadcast structure (columns where, in some country, fewer than half the units carry distinct values: parent-level estimates broadcast to districts)

| source | columns |
|---|---:|
| MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | 30 |
| GFDx (Global Fortification Data Exchange) programme fields | 23 |
| WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region | 19 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | 9 |
| FAOSTAT (national, broadcast) | 8 |
| Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | 6 |
| WHO ESPEN implementation-unit database 2014-2025 (portal export, no key) | 6 |
| HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | 5 |
| IHME (modelled surfaces) | 5 |
| WHO Global Anaemia Estimates (GHO, WHO/UNICEF joint estimates) | 5 |
| DHS | 4 |
| GFDx | 4 |
| HFID (FEWS NET / IPC / WFP mVAM) | 4 |
| MapSPAM | 3 |
| WHO FluNet (VIW_FNT export) | 3 |
| Global Data Lab Subnational HDI (SHDI) | 2 |
| JRC Global Surface Water 1.4; USDOS LSIB 2017 (Earth Engine) | 2 |
| Malaria Atlas Project (latest release, layer at each survey year) | 2 |
| WFP / HDX market prices | 2 |
| GEE | 1 |
| Meta/Data for Good RWI; WorldPop; CIESIN GPW v4.11 (Earth Engine) | 1 |
| UNICEF VAS (World Bank WDI mirror SN.ITK.VITA.ZS) | 1 |

