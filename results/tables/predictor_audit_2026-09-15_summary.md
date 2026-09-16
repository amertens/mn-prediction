# Predictor set audit, 2026-09-15

`predictors_admin2_shared.csv`: 554 Admin-2 rows x 542 predictors, 4 countries (Gambia, Ghana, Malawi, SierraLeone).

| flag | columns |
|---|---:|
| present in all 4 countries | 470 |
| partial coverage (< 4 countries) | 72 |
| national constants (no within-country variation) | 45 |
| completeness < 70% | 56 |
| value-identical duplicates | 0 |
| survey-derived (DHS) | 150 |
| modelled outcome-adjacent surfaces (V2_DROP_MODELLED sensitivity) | 48 |
| data year >= 3 years from the survey in some country (TA-01) | 306 |
| declared subnational flag disagrees with the data | 0 |

## By tier (TP-01: V2_PREDICTOR_TIERS; national constants dropped at fit time unless V2_KEEP_NATIONAL=1)

| tier | columns | sources | subnational | national const. |
|---|---:|---:|---:|---:|
| open | 328 | 24 | 283 | 45 |
| survey_dhs | 150 | 1 | 150 | 0 |
| survey_public | 64 | 3 | 64 | 0 |

## By source

| source | tier | columns | domains | all countries | subnational | national const. | <70% complete | duplicates | modelled | offset >= 3 y | max offset |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DHS | survey_dhs | 150 | 16 | 138 | 150 | 0 | 10 | 0 | 0 | 140 | 3 |
| GEE | open | 80 | 4 | 80 | 80 | 0 | 0 | 0 | 0 | 3 | 6 |
| AlphaEarth (GEE) | open | 64 | 1 | 64 | 64 | 0 | 0 | 0 | 0 | 64 | 4 |
| IHME (modelled surfaces) | open | 34 | 7 | 21 | 34 | 0 | 7 | 0 | 34 | 0 | 1 |
| MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | survey_public | 30 | 7 | 27 | 30 | 0 | 3 | 0 | 0 | 29 | 4 |
| GFDx (Global Fortification Data Exchange) programme fields | open | 23 | 1 | 21 | 0 | 23 | 1 | 0 | 0 | 0 | 0 |
| WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region | survey_public | 19 | 1 | 17 | 19 | 0 | 0 | 0 | 0 | 17 | 3 |
| SoilGrids / iSDA | open | 17 | 1 | 17 | 17 | 0 | 0 | 0 | 0 | 0 | - |
| HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | survey_public | 15 | 1 | 2 | 15 | 0 | 12 | 0 | 0 | 15 | 5 |
| MapSPAM | open | 12 | 1 | 12 | 12 | 0 | 0 | 0 | 0 | 12 | 8 |
| Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16 | open | 9 | 1 | 9 | 9 | 0 | 0 | 0 | 0 | 0 | - |
| WFP / HDX market prices | open | 9 | 1 | 5 | 8 | 1 | 4 | 0 | 0 | 0 | 0 |
| FAOSTAT (national, broadcast) | open | 8 | 1 | 8 | 0 | 8 | 0 | 0 | 0 | 0 | 0 |
| Gridded Livestock of the World 4 (2020, 5 arc-min, dasymetric; FAO catalog) | open | 8 | 1 | 8 | 8 | 0 | 0 | 0 | 0 | 8 | 3 |
| Malaria Atlas Project (latest release, layer at each survey year) | open | 8 | 1 | 8 | 8 | 0 | 0 | 0 | 0 | 0 | 0 |
| JRC Global Surface Water 1.4; USDOS LSIB 2017 (Earth Engine) | open | 6 | 1 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | - |
| Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | open | 6 | 1 | 0 | 6 | 0 | 6 | 0 | 6 | 0 | 1 |
| WHO ESPEN implementation-unit database 2014-2025 (portal export, no key) | open | 6 | 1 | 6 | 6 | 0 | 0 | 0 | 0 | 4 | 3 |
| World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | open | 6 | 1 | 0 | 6 | 0 | 6 | 0 | 0 | 0 | 0 |
| WorldPop age-sex / JRC GHS-SMOD (Earth Engine) | open | 6 | 2 | 6 | 6 | 0 | 0 | 0 | 0 | 5 | 7 |
| WHO Global Anaemia Estimates (GHO, WHO/UNICEF joint estimates) | open | 5 | 1 | 5 | 0 | 5 | 0 | 0 | 5 | 0 | 0 |
| GFDx | open | 4 | 3 | 0 | 0 | 4 | 0 | 0 | 3 | 3 | 6 |
| HFID (FEWS NET / IPC / WFP mVAM) | open | 4 | 1 | 0 | 4 | 0 | 4 | 0 | 0 | 0 | 0 |
| Meta/Data for Good RWI; WorldPop; CIESIN GPW v4.11 (Earth Engine) | open | 4 | 2 | 4 | 4 | 0 | 0 | 0 | 0 | 3 | 6 |
| Malaria Atlas | open | 3 | 1 | 3 | 3 | 0 | 0 | 0 | 0 | 3 | 6 |
| WHO FluNet (VIW_FNT export) | open | 3 | 1 | 0 | 0 | 3 | 3 | 0 | 0 | 0 | 0 |
| Global Data Lab Subnational HDI (SHDI) | open | 2 | 1 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| UNICEF VAS (World Bank WDI mirror SN.ITK.VITA.ZS) | open | 1 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |

## Temporal alignment (TA-01): columns whose data year is >= 3 years from the survey

| rule | columns | example | year used | max offset |
|---|---:|---|---|---:|
| spam_2010 | 12 | `spam_parea_total` | 2010 | 8 |
| wpop_2020 | 5 | `wpop_share_under5` | 2020 | 7 |
| access_2019 | 1 | `access_healthcare_min` | 2019 | 6 |
| gfdx_who_anaemia_2011 | 2 | `gfdx_anemia_ch_prevalence` | 2011 | 6 |
| map_blood_2012 | 3 | `map_blooddisorders201201africahbcallelefrequency` | 2012 | 6 |
| rwi_2019 | 3 | `rwi_mean` | 2019 | 6 |
| ccnl_2013 | 1 | `ntl_ccnl` | 2013 | 5 |
| gfdx_zinc_2012 | 1 | `gfdx_zinc_def_pop_prevalence` | 2012 | 5 |
| hces_year | 15 | `hces_food_share` | Gambia=2015;Ghana=2016;Malawi=2016;SierraLeone=2018 | 5 |
| alphaearth_2017 | 64 | `aef_A00` | 2017 | 4 |
| mics_microdata | 29 | `mics_salt_iodised_15ppm` | Gambia=2018;Ghana=2017;Malawi=2014;SierraLeone=2017 | 4 |
| dhs_round | 140 | `dhs_c_anemia_any` | Gambia=2019;Ghana=2014;Malawi=2015;SierraLeone=2013 | 3 |
| espen_mda | 4 | `espen_sth_mda_share` | 2016 | 3 |
| ghm_2016 | 1 | `human_modification` | 2016 | 3 |
| glw4 | 8 | `glw_cattle_km2` | 2015 | 3 |
| mics_heat_round | 17 | `mics_heat_vbcg_sy` | Gambia=2018;Ghana=2017;Malawi=2014;SierraLeone=2010 | 3 |

## Partial-coverage columns

| column | source | countries | completeness |
|---|---|---|---:|
| dhs_w_sib_n | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_dead_prop | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_adult_deaths | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_adult_death_any | DHS | Gambia|Malawi|SierraLeone | 0.53 |
| dhs_w_sib_mean_age_death | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_maternal_deaths | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_maternal_any | DHS | Gambia|Malawi|SierraLeone | 0.53 |
| dhs_w_sib_female_n | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_w_sib_female_dead | DHS | Gambia|Malawi|SierraLeone | 0.46 |
| dhs_hh_cattle | DHS | Gambia|Ghana|Malawi | 0.83 |
| dhs_hh_cattle_any | DHS | Gambia|Ghana|Malawi | 0.97 |
| dhs_hh_cook_indoors | DHS | SierraLeone | 0.02 |
| fprice_animal_rel | WFP / HDX market prices | Gambia;Ghana;SierraLeone | 0.56 |
| fprice_vegfruit_rel | WFP / HDX market prices | Gambia;Ghana;SierraLeone | 0.56 |
| fprice_oils_rel | WFP / HDX market prices | Gambia;SierraLeone | 0.09 |
| fprice_animal_to_staple | WFP / HDX market prices | Gambia;Ghana;SierraLeone | 0.56 |
| ihme_edu_0y_share | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.64 |
| ihme_edu_12plus_share | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.64 |
| ihme_edu_6_11y_share | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.64 |
| ihme_meanyearsofattainment | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.74 |
| ihme_wimpother | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.74 |
| ihme_oralrehydrationsolution | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.74 |
| ihme_orsorrhf | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.74 |
| ihme_recommendedhomefluids | IHME (modelled surfaces) | Gambia;Ghana;Malawi;SierraLeone | 0.74 |
| ihme_malecircumcisionprevalence | IHME (modelled surfaces) | Ghana;Malawi;SierraLeone | 0.80 |
| ihme_simp | IHME (modelled surfaces) | Ghana;Malawi;SierraLeone | 0.67 |
| ihme_simpother | IHME (modelled surfaces) | Ghana;Malawi;SierraLeone | 0.67 |
| ihme_spiped | IHME (modelled surfaces) | Ghana;Malawi;SierraLeone | 0.67 |
| ihme_sunimp | IHME (modelled surfaces) | Ghana;Malawi;SierraLeone | 0.67 |
| fsec_ipc_phase_fews | HFID (FEWS NET / IPC / WFP mVAM) | Malawi;SierraLeone | 0.45 |
| fsec_ipc_phase_ipcch | HFID (FEWS NET / IPC / WFP mVAM) | Gambia;Malawi | 0.50 |
| fsec_fcs_lit | HFID (FEWS NET / IPC / WFP mVAM) | Gambia;Malawi;SierraLeone | 0.53 |
| fsec_rcsi_lit | HFID (FEWS NET / IPC / WFP mVAM) | Gambia;Malawi | 0.50 |
| gfdx_anemia_ch_prevalence | GFDx | Ghana;Malawi;SierraLeone | 0.93 |
| gfdx_anemia_nonpreg_prevalence | GFDx | Ghana;Malawi;SierraLeone | 0.93 |
| gfdx_zinc_def_pop_prevalence | GFDx | Ghana;Malawi;SierraLeone | 0.93 |
| gfdx_ntd_per10k_prevalence | GFDx | Ghana;Malawi;SierraLeone | 0.93 |
| flunet_specimens_per_week_sy | WHO FluNet (VIW_FNT export) | Ghana;SierraLeone | 0.49 |
| flunet_share_positive_sy | WHO FluNet (VIW_FNT export) | Ghana;SierraLeone | 0.49 |
| flunet_influenza_a_share_sy | WHO FluNet (VIW_FNT export) | Ghana;SierraLeone | 0.49 |
| gfdx_maize_ip_pc_sy | GFDx (Global Fortification Data Exchange) programme fields | Gambia;Ghana;Malawi | 0.97 |
| gfdx_salt_ip_pc_sy | GFDx (Global Fortification Data Exchange) programme fields | Gambia;Ghana | 0.54 |
| mimi_mpi | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mimi_vita_inadequate_pct | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mimi_folate_inadequate_pct | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mimi_b12_inadequate_pct | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mimi_iron_inadequate_pct | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mimi_zinc_inadequate_pct | Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1 | Ghana | 0.47 |
| mics_heat_vrota_sy | WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region | Gambia;Ghana;Malawi | 0.97 |
| mics_heat_vrota24_35_sy | WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region | Gambia;Ghana;Malawi | 0.97 |
| hces_log_cons_pae_rel | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Ghana;Malawi;SierraLeone | 0.93 |
| hces_hdds | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_asf | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_fish | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_meat | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_eggs | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_dairy | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_pulses_nuts | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_fruit | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_veg | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_dgl | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_any_vita_fv | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| hces_asf_purchase_share | HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7) | Gambia;Malawi | 0.50 |
| rtfp_fpi_rel_national | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| rtfp_fpi_inflation_12m | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| rtfp_fpi_volatility | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| rtfp_fpi_seasonal_range | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| rtfp_staple_rel_national | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| rtfp_staple_inflation_12m | World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10 | Gambia;Malawi | 0.50 |
| mics_electricity | MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | Gambia;Ghana;SierraLeone | 0.56 |
| mics_c_cough_2wk | MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | Gambia;Ghana;SierraLeone | 0.56 |
| mics_c_vas_6mo | MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14) | Malawi | 0.44 |

## Duplicates

- none
