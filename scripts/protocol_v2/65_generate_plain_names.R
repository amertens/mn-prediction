# =============================================================================
# scripts/protocol_v2/65_generate_plain_names.R   [PN-01, 2026-09-27]
#
# Plain-language names for every predictor column that the curated map
# (R/predictor_plain_names.R) does not cover. Systematic translations of the
# column codes, one handler per source family, so the dashboard's catalogue and
# importance explorer never fall back to a raw code. Curated names always win;
# these are marked "generated" in the catalogue and the whole set is written to
# a review sheet for the RA pass (the annotation worksheet remains the
# definitive review path).
#
#   Rscript scripts/protocol_v2/65_generate_plain_names.R
# -> R/predictor_plain_names_generated.R            (PLAIN_GENERATED, sourced by builder 05)
#    results/tables/protocol_v2/plain_names_generated.csv   (RA review sheet)
# =============================================================================
md <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
source("R/predictor_plain_names.R")
cols <- setdiff(md$column, names(PLAIN))
lab <- stats::setNames(rep(NA_character_, length(cols)), cols)
put <- function(map) { k <- intersect(names(map), names(lab)); lab[k] <<- unname(map[k]); invisible(NULL) }
tc <- function(x) { s <- strsplit(x, "")[[1]]; paste0(toupper(s[1]), paste(s[-1], collapse = "")) }

# ── DHS (sensitivity tier; not in the headline model) ────────────────────────
put(c(
  dhs_c_anemia_any = "Children with any anaemia (DHS)", dhs_c_anemia_moderate = "Children with moderate anaemia (DHS)",
  dhs_c_ari_2wk = "Children with cough and fast breathing, last 2 weeks", dhs_c_bcg = "BCG vaccination",
  dhs_c_bottle_fed = "Bottle-fed children", dhs_c_diarrhea_2wk = "Children with diarrhoea, last 2 weeks",
  dhs_c_diet_diversity_score = "Child diet diversity score", dhs_c_dpt1 = "DPT first dose", dhs_c_dpt3 = "DPT third dose",
  dhs_c_early_initiation = "Early initiation of breastfeeding", dhs_c_fever_2wk = "Children with fever, last 2 weeks",
  dhs_c_fully_vaccinated = "Fully vaccinated children", dhs_c_infant_death = "Infant deaths",
  dhs_c_iron_rich_food = "Children eating iron-rich food", dhs_c_low_birthweight = "Low birthweight",
  dhs_c_mean_birthweight = "Mean birthweight", dhs_c_mean_haz = "Child height-for-age (mean)",
  dhs_c_mean_hemoglobin = "Child haemoglobin (mean)", dhs_c_mean_waz = "Child weight-for-age (mean)",
  dhs_c_mean_whz = "Child weight-for-height (mean)", dhs_c_measles1 = "Measles first dose",
  dhs_c_min_diet_diversity = "Minimum diet diversity met", dhs_c_neonatal_death = "Neonatal deaths",
  dhs_c_no_vaccination = "Children with no vaccinations", dhs_c_overweight = "Overweight children",
  dhs_c_polio3 = "Polio third dose", dhs_c_severely_stunted = "Severely stunted children",
  dhs_c_severely_wasted = "Severely wasted children", dhs_c_slept_itn = "Children who slept under a bed net",
  dhs_c_small_at_birth = "Reported small at birth", dhs_c_stunted = "Stunted children",
  dhs_c_under5_death = "Under-five deaths", dhs_c_underweight = "Underweight children",
  dhs_c_vita_supplement = "Child vitamin A supplementation", dhs_c_wasted = "Wasted children",
  dhs_c_zinc_diarrhea = "Zinc for diarrhoea", dhs_c_fg_vitA_fruitveg = "Children eating vitamin A-rich fruit and vegetables",
  dhs_c_diet_groups_n = "Child food groups eaten (count)", dhs_c_mdd_4plus = "Children eating 4+ food groups",
  dhs_c_vita_capsule = "Vitamin A capsule, last 6 months",
  dhs_CH_DIAT_C_ORT = "Diarrhoea treated with ORT", dhs_CH_VACC_C_BAS = "Basic vaccinations completed",
  dhs_CH_VACC_C_DP1 = "DPT first dose (API)", dhs_CH_VACC_C_DP3 = "DPT third dose (API)",
  dhs_CH_VACC_C_MSL = "Measles vaccination (API)", dhs_CH_VACC_C_NON = "No vaccinations (API)",
  dhs_CM_ECMR_C_NNR = "Neonatal mortality rate", dhs_RH_DELA_C_SKP = "Skilled birth attendance",
  dhs_WS_TLET_H_IMP = "Improved toilet (households)", dhs_WS_TLET_P_BAS = "Basic sanitation (population)",
  dhs_hh_electricity = "Households with electricity", dhs_hh_fridge = "Households with a fridge",
  dhs_hh_handwash_facility = "Handwashing facility at home", dhs_hh_improved_floor = "Improved floor material",
  dhs_hh_improved_sanitation = "Improved sanitation", dhs_hh_mean_hhsize = "Household size (mean)",
  dhs_hh_open_defecation = "Open defecation", dhs_hh_radio = "Households with a radio",
  dhs_hh_tv = "Households with a television", dhs_hh_urban = "Urban households",
  dhs_hh_water_30min = "Water within 30 minutes", dhs_hh_water_onpremise = "Water on premises",
  dhs_hh_wealth_mean = "Household wealth index (mean)", dhs_hh_wealth_poor = "Households in the poorest two fifths",
  dhs_hh_wealth_poorest = "Households in the poorest fifth",
  dhs_hh_horses = "Household horses", dhs_hh_horses_any = "Any household horses",
  dhs_hh_sheep_any = "Any household sheep", dhs_hh_chickens = "Household chickens",
  dhs_hh_chickens_any = "Any household chickens", dhs_hh_livestock_total = "Household livestock (total)",
  dhs_hh_livestock_any = "Any household livestock", dhs_hh_agland_ha = "Agricultural land owned (ha)",
  dhs_hh_agland_any = "Any agricultural land", dhs_hh_owns_land = "Households owning land",
  dhs_hh_clean_fuel = "Clean cooking fuel", dhs_hh_solid_fuel = "Solid cooking fuel",
  dhs_hh_cook_indoors = "Cooking indoors",
  dhs_w_anc_visits = "Antenatal visits (mean)", dhs_w_anc1 = "At least one antenatal visit",
  dhs_w_anc4plus = "Four or more antenatal visits", dhs_w_anemia_any = "Women with any anaemia (DHS)",
  dhs_w_anemia_moderate = "Women with moderate anaemia (DHS)", dhs_w_any_media = "Women reached by any media",
  dhs_w_bmi_overweight = "Overweight women (BMI)", dhs_w_delivery_csection = "Caesarean deliveries",
  dhs_w_delivery_facility = "Facility deliveries", dhs_w_edu_years_mean = "Women's schooling years (mean)",
  dhs_w_health_decision = "Women deciding on their own health care", dhs_w_high_parity = "High parity (5+ births)",
  dhs_w_listens_radio = "Women listening to radio", dhs_w_literate = "Literate women",
  dhs_w_mean_age_first_birth = "Age at first birth (mean)", dhs_w_mean_bmi = "Women's BMI (mean)",
  dhs_w_mean_parity = "Births per woman (mean)", dhs_w_reads_newspaper = "Women reading newspapers",
  dhs_w_secondary_plus = "Women with secondary schooling or more", dhs_w_skilled_attendant = "Skilled attendant at birth",
  dhs_w_teen_pregnancy = "Teenage pregnancy", dhs_w_unmet_fp = "Unmet family-planning need",
  dhs_w_watches_tv = "Women watching television", dhs_w_iron_pregnancy = "Iron in pregnancy",
  dhs_w_iron_days = "Days of iron supplementation", dhs_w_iron_90plus = "Iron for 90+ days in pregnancy",
  dhs_w_barrier_money = "Money as a barrier to care", dhs_w_barrier_distance = "Distance as a barrier to care",
  dhs_w_barriers_n = "Barriers to care (count)", dhs_w_sib_n = "Siblings reported (mean)",
  dhs_w_sib_dead_prop = "Share of siblings who died", dhs_w_sib_adult_deaths = "Adult sibling deaths",
  dhs_w_sib_adult_death_any = "Any adult sibling death", dhs_w_sib_mean_age_death = "Sibling age at death (mean)",
  dhs_w_sib_maternal_deaths = "Maternal sibling deaths", dhs_w_sib_female_n = "Sisters reported (mean)",
  dhs_w_sib_female_dead = "Sisters who died"
))

# ── AlphaEarth satellite embedding ──────────────────────────────────────────
ae <- grep("^aef_A[0-9]+$", cols, value = TRUE)
put(stats::setNames(sprintf("Satellite embedding, axis %d", as.integer(sub("^aef_A", "", ae))), ae))

# ── Soil (iSDA / SoilGrids) ─────────────────────────────────────────────────
soil_prop <- c(aluminium = "aluminium", calcium = "calcium", cec = "cation-exchange capacity",
               iron = "extractable iron", magnesium = "magnesium", nitrogen = "nitrogen",
               ph = "pH", phosphorous = "phosphorus", phosphorus = "phosphorus", potassium = "potassium",
               zinc = "extractable zinc", carbon = "organic carbon", oc = "organic carbon",
               sulfur = "extractable sulphur", totalcarbon = "total carbon",
               sand = "sand share", silt = "silt share", clay = "clay share", bulk = "bulk density",
               texture = "texture class", stone = "stone content")
for (cc in grep("^soil_", cols, value = TRUE)) {
  m <- regmatches(cc, regexec("^soil_([a-z]+)_(mean|stdev|sd)_([0-9]+)_([0-9]+)$", cc))[[1]]
  if (length(m) == 5 && m[2] %in% names(soil_prop))
    lab[cc] <- sprintf("Soil %s, %s-%s cm%s", soil_prop[[m[2]]], m[4], m[5], if (m[3] != "mean") " (variability)" else "")
}
put(c(soilgrids_cec = "Soil cation-exchange capacity (global)", soilgrids_clay = "Soil clay share (global)",
      soilgrids_nitrogen = "Soil nitrogen (global)", soilgrids_organic_carbon = "Soil organic carbon (global)",
      soilgrids_ph = "Soil pH (global)", soilgrids_sand = "Soil sand share (global)", soilgrids_silt = "Soil silt share (global)"))

# ── MICS via WHO HEAT (vaccination) + MICS microdata ────────────────────────
vax <- c(vbcg = "BCG vaccination", vdpt = "DPT third dose", vfull = "Fully vaccinated",
         vhib = "Hib vaccination", vmsl = "Measles vaccination", vpolio = "Polio third dose",
         vrota = "Rotavirus vaccination", vtet = "Tetanus-protected births", vtetprot = "Tetanus-protected births",
         vzdpt = "Zero-dose DPT", vzero = "Zero-dose children")
for (cc in grep("^mics_heat_", cols, value = TRUE)) {
  m <- regmatches(cc, regexec("^mics_heat_(v[a-z]+?)(24_35)?_sy$", cc))[[1]]
  if (length(m) >= 2 && m[2] %in% names(vax))
    lab[cc] <- paste0(vax[[m[2]]], if (identical(m[3], "24_35")) ", 24-35 months" else "", " (MICS)")
}
put(c(mics_salt_iodised_15ppm = "Adequately iodised household salt (15+ ppm)", mics_salt_any_iodine = "Any iodine in household salt",
      mics_water_improved = "Improved water source (MICS)", mics_water_piped = "Piped water (MICS)",
      mics_sanitation_improved = "Improved sanitation (MICS)", mics_handwash_soap = "Handwashing place with soap (MICS)",
      mics_head_no_education = "Household head with no schooling", mics_c_wasted = "Wasted children (MICS)",
      mics_c_underweight = "Underweight children (MICS)", mics_c_mean_haz = "Child height-for-age (MICS mean)",
      mics_c_fever_2wk = "Children with fever, last 2 weeks (MICS)", mics_c_cough_2wk = "Children with cough, last 2 weeks (MICS)",
      mics_c_fg_n = "Child food groups eaten (MICS)", mics_w_literate = "Literate women (MICS)"))

# ── GFDx fortification ──────────────────────────────────────────────────────
veh <- c(wheat = "wheat flour", maize = "maize flour", oil = "edible oil", salt = "salt", rice = "rice")
nut <- c(iron = "iron", vita = "vitamin A", folic = "folic acid", b12 = "vitamin B12", zinc = "zinc", iodine = "iodine")
for (cc in grep("^gfdx_", cols, value = TRUE)) {
  m <- regmatches(cc, regexec("^gfdx_([a-z0-9]+)_(mandatory|years_mandatory|intake_g|ip_pc|delivery_mg)_sy$", cc))[[1]]
  if (length(m) == 3) {
    v <- m[2]; f <- m[3]
    lab[cc] <- if (v %in% names(veh)) switch(f,
        mandatory = sprintf("Mandatory fortification of %s", veh[[v]]),
        years_mandatory = sprintf("Years of mandatory %s fortification", veh[[v]]),
        intake_g = sprintf("Daily intake of %s (g)", veh[[v]]),
        ip_pc = sprintf("Industrially processed share of %s", veh[[v]]))
      else if (v %in% names(nut) && f == "delivery_mg") sprintf("Fortification delivery of %s (mg/day)", nut[[v]])
      else NA_character_
  }
}
put(c(gfdx_anemia_ch_prevalence = "Anaemia in children (GFDx national)", gfdx_anemia_nonpreg_prevalence = "Anaemia in non-pregnant women (GFDx national)",
      gfdx_zinc_def_pop_prevalence = "Zinc deficiency (GFDx national)", gfdx_ntd_per10k_prevalence = "Neural-tube defects per 10,000 births",
      gfdx_n_vehicles_mandatory_sy = "Food vehicles under mandatory fortification (count)"))

# ── IHME modelled surfaces ──────────────────────────────────────────────────
put(c(ihme_edu_0y_share = "Women with no schooling (modelled)", ihme_edu_12plus_share = "Women with 12+ years of schooling (modelled)",
      ihme_edu_6_11y_share = "Women with 6-11 years of schooling (modelled)", ihme_hivprevalence = "HIV prevalence (modelled)",
      ihme_malaria_pfpr = "Malaria parasite rate (modelled)", ihme_mcvcoverage = "Measles vaccine coverage (modelled)",
      ihme_meanyearsofattainment = "Women's schooling years (modelled)", ihme_oncho_prevalence = "Onchocerciasis prevalence (modelled)",
      ihme_lf_prevalence = "Lymphatic filariasis prevalence (modelled)", ihme_wimp = "Improved water (modelled)",
      ihme_wimpother = "Other improved water (modelled)", ihme_wpiped = "Piped water (modelled)",
      ihme_wsurface = "Surface water use (modelled)", ihme_wunimp = "Unimproved water (modelled)",
      ihme_oralrehydrationsolution = "ORS for diarrhoea (modelled)", ihme_orsorrhf = "ORS or home fluids for diarrhoea (modelled)",
      ihme_recommendedhomefluids = "Recommended home fluids for diarrhoea (modelled)", ihme_u5_diarrhoea_prev = "Under-five diarrhoea prevalence (modelled)",
      ihme_malecircumcisionprevalence = "Male circumcision prevalence (modelled)", ihme_sod = "Open defecation (modelled)",
      ihme_spiped = "Piped sanitation (modelled)", ihme_sunimp = "Unimproved sanitation (modelled)"))

# ── Climate: TerraClimate normals / survey-year / night temperature ─────────
put(c(clim_pr_ann_sd = "Rainfall, year-to-year variability", clim_pr_cv = "Rainfall variability (CV)",
      clim_pr_top3_share = "Rainfall concentration (wettest 3 months)", clim_pr_win_anom_z = "Fieldwork-window rainfall anomaly",
      clim_tmax_range = "Daytime temperature, seasonal range", clim_tmax_sy_anom = "Survey-year heat anomaly",
      clim_tmin_ann = "Night temperature, annual mean", clim_pet_ann = "Potential evapotranspiration (annual)",
      clim_def_ann = "Climatic water deficit (annual)", clim_aet_ann = "Actual evapotranspiration (annual)",
      clim_soilm_ann = "Soil moisture (annual)", clim_vpd_ann = "Vapour-pressure deficit (annual)",
      clim_srad_ann = "Solar radiation (annual)", clim_lstd_ann = "Land surface temperature, day (annual)",
      clim_lstn_ann = "Land surface temperature, night (annual)"))
tstem <- c(aet = "Actual evapotranspiration", def = "Climatic water deficit", pet = "Potential evapotranspiration",
           pr = "Rainfall", ro = "Runoff", soil = "Soil moisture", srad = "Solar radiation",
           tmmn = "Minimum temperature", tmmx = "Maximum temperature", vap = "Vapour pressure",
           vpd = "Vapour-pressure deficit", vs = "Wind speed")
for (cc in grep("^tclim_", cols, value = TRUE)) {
  s <- sub("^tclim_([a-z]+)_t0$", "\\1", cc)
  if (s %in% names(tstem)) lab[cc] <- paste0(tstem[[s]], " (survey year)")
}
mon <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
for (cc in grep("^lst_night_", cols, value = TRUE)) {
  if (grepl("annual_(mean|max|min|range|sd)_t0$", cc)) {
    s <- sub(".*annual_([a-z]+)_t0$", "\\1", cc)
    lab[cc] <- paste0("Night land temperature, annual ", c(mean = "mean", max = "maximum", min = "minimum", range = "range", sd = "variability")[[s]])
  } else if (grepl("m[0-9]{2}_t0$", cc)) {
    lab[cc] <- paste0("Night land temperature, ", mon[as.integer(sub(".*m([0-9]{2})_t0$", "\\1", cc))])
  }
}

# ── Household budget surveys, crops, livestock, prices, food security ───────
put(c(hces_log_cons_pae_rel = "Household consumption per adult (relative)", hces_any_asf = "Households eating animal-source food",
      hces_any_eggs = "Households eating eggs", hces_any_pulses_nuts = "Households eating pulses or nuts",
      hces_any_fruit = "Households eating fruit", hces_any_veg = "Households eating vegetables",
      hces_any_dgl = "Households eating dark green leafy vegetables", hces_any_vita_fv = "Households eating vitamin A-rich fruit and vegetables",
      hces_asf_purchase_share = "Animal-source food bought rather than produced",
      spam_parea_total = "Total cropped area", spam_prod_cereals = "Cereal production", spam_prod_oilcrops = "Oil-crop production",
      spam_prod_pulses = "Pulse production", spam_prod_roots = "Root-crop production", spam_prod_total = "Total crop production",
      spam_prod_vegetables = "Vegetable production", spam_share_pulses = "Cropped-area share: pulses",
      spam_share_vegetables = "Cropped-area share: vegetables",
      glw_goats_km2 = "Goat density", glw_chickens_km2 = "Chicken density", glw_tlu_km2 = "Livestock density (tropical units)",
      fao_kcal_cap_day = "Calorie supply per person (national)", fao_protein_g = "Protein supply per person (national)",
      fao_fat_g = "Fat supply per person (national)", fao_supply_cereals_kg = "Cereal supply per person (national)",
      fao_supply_roots_kg = "Root-crop supply per person (national)", fao_supply_pulses_kg = "Pulse supply per person (national)",
      fao_supply_veg_kg = "Vegetable supply per person (national)", fao_supply_fruit_kg = "Fruit supply per person (national)",
      fprice_dist_nearest_market_km = "Distance to nearest food market", fprice_n_markets_100km = "Food markets within 100 km",
      fprice_animal_rel = "Animal-food price (relative)", fprice_pulses_rel = "Pulse price (relative)",
      fprice_vegfruit_rel = "Fruit and vegetable price (relative)", fprice_oils_rel = "Cooking-oil price (relative)",
      fprice_staple_volatility = "Staple price volatility", fprice_animal_to_staple = "Animal-to-staple price ratio",
      rtfp_fpi_rel_national = "Food-price level vs national (RTFP)", rtfp_staple_inflation_12m = "Staple inflation, 12 months (RTFP)",
      fsec_ipc_phase_fews = "Food-insecurity phase (FEWS NET)", fsec_ipc_phase_ipcch = "Food-insecurity phase (IPC)",
      fsec_fcs_lit = "Food consumption score (reported)", fsec_rcsi_lit = "Coping strategies index (reported)",
      mimi_mpi = "Micronutrient inadequacy index (modelled)", mimi_vita_inadequate_pct = "Vitamin A intake inadequacy (modelled)",
      mimi_folate_inadequate_pct = "Folate intake inadequacy (modelled)", mimi_b12_inadequate_pct = "Vitamin B12 intake inadequacy (modelled)",
      mimi_iron_inadequate_pct = "Iron intake inadequacy (modelled)", mimi_zinc_inadequate_pct = "Zinc intake inadequacy (modelled)"))

# ── Conflict, environment, access, remaining singletons ─────────────────────
put(c(acled_months_with_event_36m = "Months with conflict events, last 3 years", acled_any_event_36m = "Any conflict event, last 3 years",
      acled_dist_nearest_event_km = "Distance to nearest conflict event", acled_events_36m = "Conflict events, last 3 years",
      acled_fatalities_36m = "Conflict fatalities, last 3 years", acled_violent_events_36m = "Violent conflict events, last 3 years",
      acled_civilian_targeting_36m = "Civilian-targeting events, last 3 years", acled_protest_riot_36m = "Protests and riots, last 3 years",
      ndvi_anomaly_t0 = "Greenness anomaly (survey year)", ndvi_modis_t0 = "Greenness (survey year)",
      ndvi_modis_win = "Greenness in the fieldwork window", ndvi_modis_peak_t0 = "Peak greenness (survey year)",
      ndvi_modis_amp_t0 = "Greenness seasonal amplitude", ndvi_modis_clim = "Greenness (long-run average)",
      ndvi_modis_anom_t0 = "Greenness anomaly vs long run",
      espen_sth_prev_mid = "Soil-transmitted helminth prevalence", espen_sth_mda_share = "Deworming campaign share (STH)",
      espen_sch_prev_mid = "Schistosomiasis prevalence", espen_sch_mda_share = "Treatment campaign share (schistosomiasis)",
      espen_sch_cov_mean = "Schistosomiasis treatment coverage",
      map_blooddisorders201201globalg6pddallelefrequency = "G6PD deficiency allele frequency (modelled)",
      map_sy_itn_access = "Bed-net access (modelled, survey year)", map_sy_itn_use = "Bed-net use (modelled, survey year)",
      map_sy_irs_coverage = "Indoor spraying coverage (modelled, survey year)", map_sy_effective_treatment = "Effective malaria treatment (modelled, survey year)",
      who_anaemia_child_pct_sy = "Child anaemia (WHO estimate, survey year)", who_anaemia_child_trend5_pp = "Child anaemia trend, 5 years (WHO)",
      who_anaemia_nonpreg_pct_sy = "Anaemia in non-pregnant women (WHO estimate)", who_anaemia_nonpreg_trend5_pp = "Women's anaemia trend, 5 years (WHO)",
      who_anaemia_pregnant_pct_sy = "Anaemia in pregnant women (WHO estimate)",
      aez16_purity = "Agro-ecological uniformity", aez16_n_classes = "Agro-ecological classes (count)",
      koppen_purity = "Climate-zone uniformity", koppen_n_classes = "Climate zones (count)",
      lcover_bare_frac_t0 = "Bare-ground cover", lcover_shrub_frac_t0 = "Shrubland cover",
      lcover_tree_frac_t0 = "Tree cover", lcover_urban_frac_t0 = "Built-up cover",
      wpop_share_wra = "Share of women of reproductive age", wpop_share_over60 = "Share of people over 60",
      wpop_sex_ratio_wra = "Sex ratio, reproductive ages", wpop_log_density_survey_year = "Population density (log, survey year)",
      built_surface = "Built-up surface", built_surface_nres = "Non-residential built-up surface",
      gdl_shdi_sy = "Subnational Human Development Index", gdl_shdi_trend5 = "Human Development Index trend, 5 years",
      rwi_mean = "Relative wealth index (modelled)", rwi_n_points = "Wealth-index data points (count)",
      wdist_perm_km_min = "Distance to permanent water", wdist_any_km_min = "Distance to any surface water",
      access_healthcare_min = "Travel time to health care", aod_t0 = "Air pollution (aerosol optical depth)",
      evi_t0 = "Vegetation index (EVI)", ghs_pop = "Population (GHS)", ghsl_smod_mean = "Urbanisation grade (GHS-SMOD)",
      human_modification = "Human modification of land", lai_t0 = "Leaf-area index", ntl_ccnl = "Night-time lights",
      wapor_mean_t0 = "Vegetation productivity (survey year)"))
for (cc in grep("^(aez16|koppen)_is_[0-9]+$", cols, value = TRUE))
  lab[cc] <- sprintf("Share in %s class %s", if (startsWith(cc, "aez16")) "agro-ecological" else "climate", sub(".*_is_", "", cc))
for (cc in grep("^vas_", cols, value = TRUE)) lab[cc] <- "Vitamin A supplementation coverage (UNICEF)"

# ── write outputs ────────────────────────────────────────────────────────────
ok <- !is.na(lab)
cat(sprintf("generated %d of %d missing labels (%d still unnamed)\n", sum(ok), length(lab), sum(!ok)))
if (any(!ok)) cat("unnamed:", paste(names(lab)[!ok], collapse = ", "), "\n")
gen <- lab[ok]
rev <- data.frame(column = names(lab), generated_label = unname(lab),
                  source = md$source[match(names(lab), md$column)],
                  domain = md$domain[match(names(lab), md$column)], stringsAsFactors = FALSE)
write.csv(rev, "results/tables/protocol_v2/plain_names_generated.csv", row.names = FALSE)
f <- file("R/predictor_plain_names_generated.R", "w", encoding = "UTF-8")
writeLines(c(
  "# =============================================================================",
  "# R/predictor_plain_names_generated.R -- WRITTEN BY scripts/protocol_v2/65_generate_plain_names.R",
  "#",
  "# Systematic plain-language translations of predictor column codes, one",
  "# handler per source family. Curated names in R/predictor_plain_names.R",
  "# always take precedence; edit THERE (or rerun script 65), never here.",
  sprintf("# Generated %s: %d labels. Review sheet:", format(Sys.time(), "%Y-%m-%d"), length(gen)),
  "# results/tables/protocol_v2/plain_names_generated.csv",
  "# =============================================================================",
  "", "PLAIN_GENERATED <- c("), f)
writeLines(paste0("  ", sprintf('`%s` = "%s"', names(gen), gen), c(rep(",", length(gen) - 1), "")), f)
writeLines(")", f)
close(f)
cat("wrote R/predictor_plain_names_generated.R and the review CSV\n")
