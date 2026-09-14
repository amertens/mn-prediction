# =============================================================================
# R/predictor_plain_names.R
#
# Plain-language names for predictor columns, shared by the policy-deck figures
# (scripts/policy_deck/01_figures_main.R) and the dashboard's predictor
# catalogue (dashboard/data-raw/05_build_protocol_v2_bundles.R). These are the
# only per-variable descriptions the project has: the variable annotation sheet
# (results/tables/protocol_v2/variable_sheet.csv) carries a mechanism template
# per sub-domain, shared by many columns, and its RA-verified columns are empty.
# Add a name here when a column reaches a leading list or a catalogue user asks.
# =============================================================================

PLAIN <- c(
  dhs_CN_NUTS_C_HA2 = "Stunted children (survey)", dhs_CN_NUTS_C_WH2 = "Wasted children (survey)",
  dhs_w_height_low = "Short-stature women", dhs_w_birth_interval_short = "Short birth intervals",
  dhs_w_health_insurance = "Women with health insurance", dhs_FP_CUSA_W_MOD = "Modern contraceptive use",
  dhs_w_modern_fp = "Modern family planning use", dhs_hh_soap_available = "Households with soap",
  dhs_w_barrier_permission = "Need permission to seek care", dhs_w_barrier_alone = "Not wanting to go alone for care",
  dhs_w_bmi_low = "Thin women (BMI)", dhs_AN_NUTS_W_THN = "Thin women (survey)",
  dhs_c_fg_roots = "Children eating roots and tubers", dhs_c_fg_legumes = "Children eating legumes",
  dhs_c_fg_grains = "Children eating grains", dhs_c_fg_dairy = "Children eating dairy",
  dhs_c_fg_flesh = "Children eating meat or fish", dhs_c_fg_eggs = "Children eating eggs",
  dhs_c_fg_other_fruitveg = "Children eating other fruit and vegetables",
  dhs_CN_BRFS_C_EXB = "Exclusive breastfeeding", dhs_hh_cows = "Household cattle ownership",
  dhs_hh_cows_any = "Any household cattle", dhs_hh_cattle = "Household cattle", dhs_hh_cattle_any = "Any household cattle",
  dhs_hh_goats = "Household goats", dhs_hh_goats_any = "Any household goats", dhs_hh_sheep = "Household sheep",
  dhs_w_decides_earnings = "Women deciding on earnings", dhs_w_occ_agric = "Women in farm work",
  dhs_hh_improved_water = "Improved water source", dhs_hh_crowding = "Household crowding",
  dhs_w_deworm_pregnancy = "Deworming in pregnancy", dhs_hh_itn_any = "Bed net in household",
  dhs_ML_NETP_H_IT2 = "Bed net per two people", dhs_c_deworm_6mo = "Child deworming",
  dhs_w_sib_maternal_any = "Sibling maternal death", dhs_w_owns_house = "Women owning their home",
  dhs_w_working = "Women in paid work", dhs_w_primary_edu = "Women with primary schooling",
  dhs_w_no_education = "Women with no schooling",
  ihme_severeanemia = "Modelled severe anaemia", ihme_allanemia = "Modelled anaemia (any)",
  ihme_moderateanemia = "Modelled moderate anaemia", ihme_mildanemia = "Modelled mild anaemia",
  ihme_anemia = "Modelled anaemia", ihme_wastingprevalence = "Modelled child wasting",
  ihme_stuntingprevalence = "Modelled child stunting", ihme_underweightprevalence = "Modelled child underweight",
  espen_sth_cov_mean = "Deworming coverage", tclim_pdsi_t0 = "Drought index",
  lcover_crops_frac_t0 = "Cropland cover", lcover_grass_frac_t0 = "Grassland cover",
  lcover_water_seasonal_frac_t0 = "Seasonal water cover", wapor_sd_t0 = "Vegetation seasonality",
  grassland_frac = "Grassland share", npp_gpp_t0 = "Vegetation growth (gross)", npp_npp_t0 = "Vegetation growth (net)",
  wdist_coast_km_mean = "Distance to the coast", wdist_coast_km_min = "Distance to the coast (nearest)",
  wdist_perm_km_mean = "Distance to permanent water", wdist_any_km_mean = "Distance to any water",
  elevation = "Elevation", glw_ruminant_share = "Ruminant share of livestock", glw_cattle_km2 = "Cattle density",
  glw_pigs_km2 = "Pig density", glw_sheep_km2 = "Sheep density", glw_tlu_per_capita = "Livestock per person",
  spam_share_cereals = "Cereal share of cropland", spam_share_oilcrops = "Oil-crop share of cropland",
  spam_share_roots = "Root-crop share of cropland",
  soil_phosphorus_stdev_0_20 = "Soil phosphorus variability", soil_zinc_stdev_0_20 = "Soil zinc variability",
  map_blooddisorders201201africahbcallelefrequency = "Haemoglobin C gene frequency",
  map_blooddisorders201201globalsicklehaemoglobinhbsallelefrequency = "Sickle-cell gene frequency",
  map_sy_pf_mortality_rate = "Malaria mortality", map_sy_pf_incidence_rate = "Malaria incidence",
  map_sy_pf_parasite_rate = "Malaria parasite rate", map_sy_itn_use_rate = "Bed-net use rate",
  rwi_sd = "Wealth index spread", wpop_share_under5 = "Share of people under five",
  wpop_dependency_ratio = "Dependency ratio", fprice_staple_rel = "Relative staple food price"
)
PLAIN <- PLAIN[!duplicated(names(PLAIN))]

#' A readable code for columns without a plain name: strip the source prefix
#' and the time suffix, and replace underscores.
clean_code <- function(x) {
  gsub("_", " ", sub("_t0$", "", sub("^(dhs|glw|ihme|map|spam|lcover|wdist|tclim|soil|wpop|npp|wapor|fprice|espen|rwi)_", "", x)))
}

#' Plain name where one exists, cleaned code otherwise. `warn = TRUE` reports
#' the unnamed columns (the figure scripts use it so a new leading predictor
#' gets a name rather than a code in print).
plain_of <- function(x, warn = FALSE) {
  out <- unname(PLAIN[x])
  if (any(is.na(out))) {
    if (warn) warning("no plain-language name for: ", paste(x[is.na(out)], collapse = ", "), " (cleaned code used)")
    out[is.na(out)] <- clean_code(x[is.na(out)])
  }
  out
}
