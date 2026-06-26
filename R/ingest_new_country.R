# =============================================================================
# R/ingest_new_country.R
#
# Country-agnostic ingestion + harmonization helpers for ADDING a new country
# to the LOCO panel. Wires next-step candidates from scripts/ingest_new_country/.
#
# Status notes (2026-06-04):
#   - Kenya KNMS 2011        — waiting on MoH/Nutrition International DUA
#   - Tanzania NNS 2018      — waiting on TFNC DUA
#   - Ethiopia NMS 2016      — waiting on EPHI DUA
#   - KDHS 2022              — user can download via existing DHS account
#   - Tanzania DHS 2022      — user can download via existing DHS account
#
# Once a file lands in data/raw/<country>/, add a get_country_configs() entry
# in R/config.R that references the new harmonized .rds via data_path, and the
# entire pipeline (LOCO, calibration, SHAP, dashboard) picks it up.
#
# This module exports:
#   harmonize_dhs_hb_to_anaemia()  DHS Hb → WHO-cutoff anaemia binary
#   harmonize_rbp_thurnham()       BRINDA cRBP inflammation adjustment
#   cap_b12_outliers()             Drop B12 > 1500 pmol/L (per Ghana fix)
#   cap_folate_outliers()          Drop folate > 50 nmol/L (per Ghana fix)
#   apply_who_cutoff()             Generic binary outcome from continuous biomarker
#   ingest_dhs_anaemia_only()      Build minimal LOCO-ready dataset from DHS Hb
#   stub_new_country_config()      Template config block for a new country
# =============================================================================

#' DHS Hb (g/dL) → WHO anaemia binary (any severity)
#'
#' WHO cutoffs (2024 update):
#'   Children 6-59 months:        Hb < 11.0 g/dL
#'   Non-pregnant women 15-49:    Hb < 12.0 g/dL
#'   Pregnant women 15-49:        Hb < 11.0 g/dL
#'
#' DHS data already adjusts Hb for altitude (>1000m) and smoking; further
#' adjustment is not needed. We don't apply BRINDA inflammation adjustment
#' here since DHS doesn't release individual AGP/CRP values.
#'
#' @param hb_g_dL numeric vector, haemoglobin in g/dL (DHS standard)
#' @param population one of "children_6_59m", "women_nonpreg", "women_pregnant"
#' @return integer 0/1 vector (1 = anaemia)
harmonize_dhs_hb_to_anaemia <- function(hb_g_dL, population) {
  stopifnot(population %in% c("children_6_59m", "women_nonpreg", "women_pregnant"))
  cutoff <- switch(population,
                    children_6_59m  = 11.0,
                    women_nonpreg   = 12.0,
                    women_pregnant  = 11.0)
  as.integer(hb_g_dL < cutoff)
}

#' BRINDA Thurnham inflammation adjustment for serum/plasma RBP
#'
#' Reference: Thurnham et al. 2003, Am J Clin Nutr (re-derived in BRINDA
#' Suchdev et al. 2016). Multiplicative corrections:
#'   incubation:    RBP * 1.13
#'   early conval.: RBP * 1.24
#'   late conval.:  RBP * 1.11
#' "Incubation" = CRP > 5 mg/L AND AGP <= 1 g/L
#' "Early conv." = CRP > 5 mg/L AND AGP > 1 g/L
#' "Late conv."  = CRP <= 5 mg/L AND AGP > 1 g/L
#'
#' @param rbp serum/plasma RBP in umol/L
#' @param crp_mg_L CRP in mg/L
#' @param agp_g_L AGP in g/L
#' @return adjusted RBP vector
harmonize_rbp_thurnham <- function(rbp, crp_mg_L, agp_g_L) {
  factor <- rep(1.0, length(rbp))
  inc  <- !is.na(crp_mg_L) & !is.na(agp_g_L) & crp_mg_L > 5 & agp_g_L <= 1
  ec   <- !is.na(crp_mg_L) & !is.na(agp_g_L) & crp_mg_L > 5 & agp_g_L > 1
  lc   <- !is.na(crp_mg_L) & !is.na(agp_g_L) & crp_mg_L <= 5 & agp_g_L > 1
  factor[inc] <- 1.13
  factor[ec]  <- 1.24
  factor[lc]  <- 1.11
  rbp * factor
}

#' Outlier cap for serum B12: anything above 1500 pmol/L is biologically
#' implausible without supplementation. Matches the Ghana cleanup in
#' sandbox/20_loco_data_quality_fixes.R.
#' @param x numeric B12 in pmol/L
#' @return capped vector
cap_b12_outliers <- function(x) {
  out <- as.numeric(x)
  out[!is.finite(out)] <- NA_real_
  pmin(out, 1500, na.rm = FALSE)
}

#' Outlier cap for serum folate: 50 nmol/L. Matches Ghana cleanup.
cap_folate_outliers <- function(x) {
  out <- as.numeric(x)
  out[!is.finite(out)] <- NA_real_
  pmin(out, 50, na.rm = FALSE)
}

#' Generic biomarker → binary using a WHO/IOM cutoff.
#'
#' @param x numeric biomarker
#' @param cutoff numeric threshold
#' @param dir "less" if deficient when x<cutoff (default) else "greater"
#' @return integer 0/1
apply_who_cutoff <- function(x, cutoff, dir = c("less", "greater")) {
  dir <- match.arg(dir)
  x <- as.numeric(x); ok <- is.finite(x)
  out <- rep(NA_integer_, length(x))
  if (dir == "less")    out[ok] <- as.integer(x[ok] < cutoff)
  if (dir == "greater") out[ok] <- as.integer(x[ok] > cutoff)
  out
}

#' Build a minimal LOCO-ready dataset from a DHS individual-recode file with
#' Hb (variable name typically hw53 for children, hb55 for women). Adds anaemia
#' binary and standard survey design columns. Designed for KDHS 2022, Tanzania
#' DHS 2022, Burkina Faso 2021 etc. that lack a full BRINDA biomarker panel.
#'
#' NOTE: this assumes you have already aggregated DHS to an admin-2 polygon
#' membership via src/DHS/DHS_admin2_aggregation.R. The resulting dataset will
#' have ONE outcome (women_iron_anaemia) — vitA/B12/folate/zinc require a full
#' biomarker survey.
#'
#' @param dhs_rds path to a harmonized DHS rds with columns:
#'   hb_g_dL, age_months, is_pregnant, gw_cnum (cluster), Admin1, Admin2,
#'   gw_svy_weight, gw_strata
#' @param population "children_6_59m" | "women_nonpreg" | "women_pregnant"
#' @return data.frame ready to drop into the targets pipeline as merged_<cc>
ingest_dhs_anaemia_only <- function(dhs_rds, population) {
  d <- readRDS(dhs_rds)
  required <- c("hb_g_dL", "Admin1", "Admin2",
                "gw_cnum", "gw_svy_weight", "gw_strata")
  missing <- setdiff(required, colnames(d))
  if (length(missing)) stop("DHS file missing required cols: ",
                              paste(missing, collapse = ", "))
  d$anaemia <- harmonize_dhs_hb_to_anaemia(d$hb_g_dL, population)
  d$gw_child_flag <- if (population == "children_6_59m") 1L else 0L
  cat(sprintf("[ingest] %s: %d rows, %d anaemia cases (%.1f%%)\n",
              basename(dhs_rds), nrow(d), sum(d$anaemia, na.rm = TRUE),
              100 * mean(d$anaemia, na.rm = TRUE)))
  d
}

#' Return a placeholder country config block for a new country.
#'
#' Drop this in get_country_configs() in R/config.R, then point data_path at
#' a harmonized rds in data/IPD/<country>/. Adjust the outcomes list to match
#' the biomarkers actually available in the source survey.
#'
#' @param iso3 e.g. "KEN"
#' @param country_name e.g. "Kenya"
#' @param survey_year e.g. 2011L
#' @param dhs_year e.g. 2014L
#' @return named list usable as a get_country_configs() entry
stub_new_country_config <- function(iso3, country_name, survey_year, dhs_year) {
  list(
    country     = country_name,
    gadm_code   = iso3,
    survey_year = as.integer(survey_year),
    dhs_year    = as.integer(dhs_year),
    data_path   = here::here("data", "IPD", country_name,
                              sprintf("%s_merged_dataset.rds",
                                       gsub(" ", "_", country_name))),
    raster_dir  = here::here("data", sprintf("%s_GEE_rasters",
                                                gsub(" ", "_", country_name))),

    cluster_id    = "gw_cnum",
    admin1_col    = "Admin1",
    admin2_col    = "Admin2",
    psu_col       = "gw_cnum",
    weight_col    = "gw_svy_weight",
    strata_col    = "gw_strata",
    child_flag    = "gw_child_flag",

    outcomes = list(
      # Same outcome tags as existing countries; fill in survey-specific
      # variable names when the data lands.
      child_vitA = list(tag = "child_vitA", label = "Vitamin A (children)",
                         population = "children", child_flag_val = 1L,
                         continuous = NA_character_, binary = NA_character_,
                         cutoff = 0.70, cutoff_dir = "less",
                         cutoff_scale = "original"),
      women_vitA = list(tag = "women_vitA", label = "Vitamin A (women)",
                         population = "women", child_flag_val = 0L,
                         continuous = NA_character_, binary = NA_character_,
                         cutoff = 0.70, cutoff_dir = "less",
                         cutoff_scale = "original"),
      child_iron = list(tag = "child_iron", label = "Iron deficiency (children)",
                         population = "children", child_flag_val = 1L,
                         continuous = NA_character_, binary = NA_character_,
                         cutoff = log(12), cutoff_dir = "less",
                         cutoff_scale = "log"),
      women_iron = list(tag = "women_iron", label = "Iron deficiency (women)",
                         population = "women", child_flag_val = 0L,
                         continuous = NA_character_, binary = NA_character_,
                         cutoff = log(15), cutoff_dir = "less",
                         cutoff_scale = "log")
    ),

    domains = list(
      GW     = list(prefix = "gw_"),
      DHS    = list(prefix = sprintf("dhs%d_", dhs_year)),
      MICS   = list(prefix = "mics_"),
      IHME   = list(prefix = "ihme_"),
      LSMS   = list(prefix = "lsms_"),
      MAP    = list(prefix = "MAP_"),
      WFP    = list(prefix = "wfp_", extra = "nearest_market_id"),
      FLUNET = list(prefix = "flunet_"),
      GEE    = list(prefix = "gee_"),
      FSEC   = list(prefix = "fsec_"),
      SOIL   = list(prefix = "soil_")
    ),

    gw_exclude_patterns = c(
      "RBP", "rbp", "VAD", "LogFer", "logfer",
      "IDA", "Brinda", "Thurn",
      "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
    )
  )
}

# ── Pre-baked country stubs ────────────────────────────────────────────
# Uncomment in R/config.R::get_country_configs() once data lands.

CANDIDATE_STUBS <- list(
  Kenya       = stub_new_country_config("KEN", "Kenya",       2011L, 2014L),
  Tanzania    = stub_new_country_config("TZA", "Tanzania",    2018L, 2022L),
  Ethiopia    = stub_new_country_config("ETH", "Ethiopia",    2016L, 2019L),

  # DHS-only candidates (anaemia-only, women_iron outcome only):
  KenyaDHS2022       = stub_new_country_config("KEN", "Kenya DHS 2022",
                                                  2022L, 2022L),
  TanzaniaDHS2022    = stub_new_country_config("TZA", "Tanzania DHS 2022",
                                                  2022L, 2022L),
  BurkinaFasoDHS2021 = stub_new_country_config("BFA", "Burkina Faso DHS 2021",
                                                  2021L, 2021L)
)
