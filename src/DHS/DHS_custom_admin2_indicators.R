# =============================================================================
# src/DHS/DHS_custom_admin2_indicators.R
#
# Derive ~80-100 custom indicators from raw DHS microdata (IR, KR, PR, HR, BR)
# and compute direct survey-weighted admin-2 prevalences / means.
#
# This complements DHS_admin2_aggregation.R (which uses surveyPrev built-ins)
# by extracting many more risk factors for micronutrient deficiency from the
# raw recode variables.
#
# Outputs per country-year:
#   data/DHS/clean/{Country}_{Year}_dhs_custom_admin2_wide.rds
#
# Requirements:
#   - Raw DHS microdata downloaded (via getDHSdata) as .RDS in data/DHS/
#   - rdhs configured for credential access
#   - surveyPrev + SUMMER for directEST
#   - geodata for GADM boundaries
#
# Usage:
#   source("src/DHS/DHS_custom_admin2_indicators.R")
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(sf)
  library(geodata)
  library(surveyPrev)
  library(SUMMER)
  library(rdhs)
})

# =============================================================================
# 0. Configuration
# =============================================================================

COUNTRY_REGISTRY <- list(
  list(country = "Gambia",       iso3 = "GMB", years = 2019),
  list(country = "Ghana",        iso3 = "GHA", years = c(2014, 2017)),
  list(country = "Sierra Leone", iso3 = "SLE", years = 2013),
  list(country = "Malawi",       iso3 = "MWI", years = 2015)
)

# Set to a subset to only run specific countries, e.g. c("Gambia")
COUNTRIES_TO_RUN <- NULL

OUT_DIR <- here("data", "DHS", "clean")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Small helpers
# =============================================================================

#' Robust numeric coercion for haven_labelled / factor / character
as_num <- function(x) {
  if (inherits(x, "haven_labelled")) {
    suppressWarnings(as.numeric(as.character(x)))
  } else if (is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else {
    suppressWarnings(as.numeric(x))
  }
}

#' Check if all vars are present in a data frame
has_vars <- function(df, vars) all(vars %in% names(df))

#' Standardize a derived indicator into surveyPrev-compatible format
#' @param df  Raw DHS data frame (IR/KR/PR/HR/BR)
#' @param which  Which recode type: "IR", "KR", "PR", "HR", "BR"
#' @param value_var  Column name of the derived indicator
#' @return tibble with columns: cluster, householdID, v024, weight, strata, value
standardize_indicator <- function(df, which, value_var) {
  if (!value_var %in% names(df)) return(NULL)

  if (which %in% c("PR", "HR")) {
    out <- tibble(
      cluster     = as_num(df[["hv001"]]),
      householdID = as_num(df[["hv002"]]),
      v024        = df[["hv024"]],
      weight      = as_num(df[["hv005"]]) / 1e6,
      strata      = if ("hv023" %in% names(df)) df[["hv023"]] else df[["hv025"]],
      value       = as_num(df[[value_var]])
    )
  } else {
    # IR, KR, BR
    out <- tibble(
      cluster     = as_num(df[["v001"]]),
      householdID = as_num(df[["v002"]]),
      v024        = df[["v024"]],
      weight      = as_num(df[["v005"]]) / 1e6,
      strata      = if ("v023" %in% names(df)) df[["v023"]] else df[["v025"]],
      value       = as_num(df[[value_var]])
    )
  }
  out <- out %>% filter(!is.na(value) & !is.na(cluster) & !is.na(weight))
  if (nrow(out) == 0) return(NULL)
  out
}


# =============================================================================
# 2. Indicator derivation functions
#
# Each function takes a raw DHS data frame, derives indicators in-place,
# and returns a named list of indicator names that were successfully created.
# =============================================================================

# ---------------------------------------------------------------------------
# Domain 1: Maternal Nutrition & Anthropometry (IR or PR)
# ---------------------------------------------------------------------------
derive_maternal_nutrition <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  # Women's anemia (v457 = anemia level: 1=severe, 2=moderate, 3=mild, 4=not anemic)
  if ("v457" %in% names(IRdata)) {
    x <- as_num(IRdata$v457)
    IRdata$w_anemia_any      <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:3, 1L, 0L))
    IRdata$w_anemia_moderate <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "w_anemia_any", "w_anemia_moderate")
  }

  # BMI (v445 = BMI * 100, e.g., 2150 = 21.5; flagged values 9998/9999 are missing)
  if ("v445" %in% names(IRdata)) {
    bmi_raw <- as_num(IRdata$v445)
    bmi <- ifelse(bmi_raw >= 1200 & bmi_raw <= 6000, bmi_raw / 100, NA_real_)
    IRdata$w_bmi_low       <- ifelse(is.na(bmi), NA_integer_, ifelse(bmi < 18.5, 1L, 0L))
    IRdata$w_bmi_overweight <- ifelse(is.na(bmi), NA_integer_, ifelse(bmi >= 25, 1L, 0L))
    IRdata$w_mean_bmi      <- bmi  # continuous
    added <- c(added, "w_bmi_low", "w_bmi_overweight", "w_mean_bmi")
  }

  # Height (v438 = height in mm * 10; e.g., 15320 = 153.2 cm)
  if ("v438" %in% names(IRdata)) {
    ht <- as_num(IRdata$v438)
    ht_cm <- ifelse(ht >= 1000 & ht <= 2500, ht / 10, NA_real_)
    IRdata$w_height_low <- ifelse(is.na(ht_cm), NA_integer_, ifelse(ht_cm < 145, 1L, 0L))
    added <- c(added, "w_height_low")
  }

  list(df = IRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 2: Maternal Diet & Supplements (IR)
# ---------------------------------------------------------------------------
derive_maternal_diet <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  # Iron supplement during pregnancy (m45_1 = iron tablets/syrup)
  for (v in c("m45_1", "m45")) {
    if (v %in% names(IRdata) && is.null(IRdata$w_iron_supplement)) {
      x <- as_num(IRdata[[v]])
      IRdata$w_iron_supplement <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
      added <- c(added, "w_iron_supplement")
      break
    }
  }

  # Vitamin A postpartum (m54_1 = received vitamin A dose after delivery)
  for (v in c("m54_1", "m54")) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      IRdata$w_vitA_postpartum <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
      added <- c(added, "w_vitA_postpartum")
      break
    }
  }

  # Women's dietary diversity (v472a-v472g food groups consumed yesterday)
  food_vars <- c("v472a", "v472b", "v472c", "v472d", "v472e", "v472f", "v472g")
  avail_food <- intersect(food_vars, names(IRdata))
  if (length(avail_food) >= 5) {
    food_mat <- sapply(avail_food, function(v) {
      x <- as_num(IRdata[[v]])
      ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    })
    food_count <- rowSums(food_mat, na.rm = FALSE)
    IRdata$w_diet_diversity_score <- food_count  # continuous
    IRdata$w_mdd_w <- ifelse(is.na(food_count), NA_integer_, ifelse(food_count >= 5, 1L, 0L))
    added <- c(added, "w_diet_diversity_score", "w_mdd_w")

    # Individual food group indicators
    food_names <- c("w_ate_starchy", "w_ate_legumes", "w_ate_meat",
                    "w_ate_eggs", "w_ate_dairy", "w_ate_vita_fruit", "w_ate_other_fruit")
    for (i in seq_along(avail_food)) {
      if (i <= length(food_names)) {
        x <- as_num(IRdata[[avail_food[i]]])
        IRdata[[food_names[i]]] <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
        added <- c(added, food_names[i])
      }
    }
  }

  list(df = IRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 3: Reproductive Health & ANC (IR)
# ---------------------------------------------------------------------------
derive_reproductive_health <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  # ANC visits (m14_1 = number of antenatal visits for most recent birth)
  for (v in c("m14_1", "m14")) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      x <- ifelse(x >= 0 & x <= 20, x, NA_real_)  # flag impossibles
      IRdata$w_anc4plus   <- ifelse(is.na(x), NA_integer_, ifelse(x >= 4, 1L, 0L))
      IRdata$w_anc1       <- ifelse(is.na(x), NA_integer_, ifelse(x >= 1, 1L, 0L))
      IRdata$w_anc_visits <- x  # continuous
      added <- c(added, "w_anc4plus", "w_anc1", "w_anc_visits")
      break
    }
  }

  # Delivery at facility (m15_1 = place of delivery)
  for (v in c("m15_1", "m15")) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      # Facility = codes 20-36 (public/private health facilities)
      IRdata$w_delivery_facility <- ifelse(is.na(x), NA_integer_,
                                            ifelse(x >= 20 & x <= 36, 1L,
                                                   ifelse(x %in% c(11, 12, 96), 0L, NA_integer_)))
      added <- c(added, "w_delivery_facility")
      break
    }
  }

  # C-section (m17_1)
  for (v in c("m17_1", "m17")) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      IRdata$w_delivery_csection <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
      added <- c(added, "w_delivery_csection")
      break
    }
  }

  # Skilled birth attendant (m3a_1 through m3g_1: doctor, nurse/midwife, etc.)
  sba_vars <- paste0("m3a_1")  # doctor
  sba_nurse <- paste0("m3b_1") # nurse/midwife
  sba_avail <- intersect(c("m3a_1", "m3b_1", "m3c_1", "m3a", "m3b", "m3c"), names(IRdata))
  if (length(sba_avail) >= 1) {
    sba_mat <- sapply(sba_avail, function(v) as_num(IRdata[[v]]) == 1)
    IRdata$w_skilled_attendant <- ifelse(rowSums(sba_mat, na.rm = TRUE) > 0, 1L, 0L)
    IRdata$w_skilled_attendant[rowSums(!is.na(sba_mat)) == 0] <- NA_integer_
    added <- c(added, "w_skilled_attendant")
  }

  # Unmet need for family planning (v626a: 1-2 = unmet need)
  if ("v626a" %in% names(IRdata)) {
    x <- as_num(IRdata$v626a)
    IRdata$w_unmet_fp <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "w_unmet_fp")
  }

  # Modern contraception (v313: 3 = modern method)
  if ("v313" %in% names(IRdata)) {
    x <- as_num(IRdata$v313)
    IRdata$w_modern_fp <- ifelse(is.na(x), NA_integer_, ifelse(x == 3, 1L, 0L))
    added <- c(added, "w_modern_fp")
  }

  # High parity (v201 = total children ever born)
  if ("v201" %in% names(IRdata)) {
    x <- as_num(IRdata$v201)
    IRdata$w_high_parity <- ifelse(is.na(x), NA_integer_, ifelse(x >= 5, 1L, 0L))
    IRdata$w_mean_parity <- ifelse(x >= 0 & x <= 20, x, NA_real_)
    added <- c(added, "w_high_parity", "w_mean_parity")
  }

  # Teen pregnancy (v012 = age, v213 = currently pregnant)
  if (has_vars(IRdata, c("v012"))) {
    age <- as_num(IRdata$v012)
    # Teen = 15-19 who have ever had a birth or are currently pregnant
    preg <- if ("v213" %in% names(IRdata)) as_num(IRdata$v213) == 1 else FALSE
    parity <- if ("v201" %in% names(IRdata)) as_num(IRdata$v201) > 0 else FALSE
    IRdata$w_teen_pregnancy <- ifelse(is.na(age), NA_integer_,
                                      ifelse(age >= 15 & age <= 19 & (preg | parity), 1L,
                                             ifelse(age >= 15 & age <= 19, 0L, NA_integer_)))
    added <- c(added, "w_teen_pregnancy")
  }

  # Short birth interval (b11_01 = preceding birth interval in months)
  for (v in c("b11_01", "b11")) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      x <- ifelse(x > 0 & x <= 120, x, NA_real_)
      IRdata$w_birth_interval_short <- ifelse(is.na(x), NA_integer_, ifelse(x < 24, 1L, 0L))
      added <- c(added, "w_birth_interval_short")
      break
    }
  }

  list(df = IRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 4: Maternal Education & Empowerment (IR)
# ---------------------------------------------------------------------------
derive_maternal_education <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  # Education level (v106: 0=none, 1=primary, 2=secondary, 3=higher)
  if ("v106" %in% names(IRdata)) {
    x <- as_num(IRdata$v106)
    IRdata$w_no_education    <- ifelse(is.na(x), NA_integer_, ifelse(x == 0, 1L, 0L))
    IRdata$w_primary_edu     <- ifelse(is.na(x), NA_integer_, ifelse(x >= 1, 1L, 0L))
    IRdata$w_secondary_plus  <- ifelse(is.na(x), NA_integer_, ifelse(x >= 2, 1L, 0L))
    added <- c(added, "w_no_education", "w_primary_edu", "w_secondary_plus")
  }

  # Years of education (v133)
  if ("v133" %in% names(IRdata)) {
    x <- as_num(IRdata$v133)
    IRdata$w_edu_years_mean <- ifelse(x >= 0 & x <= 30, x, NA_real_)
    added <- c(added, "w_edu_years_mean")
  }

  # Literacy (v155: 0=cannot read, 1=parts, 2=whole sentence)
  if ("v155" %in% names(IRdata)) {
    x <- as_num(IRdata$v155)
    IRdata$w_literate <- ifelse(is.na(x), NA_integer_, ifelse(x == 2, 1L, 0L))
    added <- c(added, "w_literate")
  }

  # Health decision-making (v743a: 1=respondent alone, 2=respondent+husband)
  if ("v743a" %in% names(IRdata)) {
    x <- as_num(IRdata$v743a)
    IRdata$w_health_decision <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "w_health_decision")
  }

  # Media exposure: reads newspaper (v157), watches TV (v158), listens radio (v159)
  # 0=not at all, 1=less than once a week, 2=at least once a week
  media_vars <- c(v157 = "w_reads_newspaper", v158 = "w_watches_tv", v159 = "w_listens_radio")
  media_count <- rep(0L, nrow(IRdata))
  media_n <- 0
  for (v in names(media_vars)) {
    if (v %in% names(IRdata)) {
      x <- as_num(IRdata[[v]])
      IRdata[[media_vars[v]]] <- ifelse(is.na(x), NA_integer_, ifelse(x >= 1, 1L, 0L))
      added <- c(added, media_vars[v])
      media_count <- media_count + ifelse(!is.na(x) & x >= 1, 1L, 0L)
      media_n <- media_n + 1
    }
  }
  if (media_n > 0) {
    IRdata$w_any_media <- ifelse(media_count > 0, 1L, 0L)
    added <- c(added, "w_any_media")
  }

  # Age at first birth (v212)
  if ("v212" %in% names(IRdata)) {
    x <- as_num(IRdata$v212)
    IRdata$w_mean_age_first_birth <- ifelse(x >= 10 & x <= 50, x, NA_real_)
    added <- c(added, "w_mean_age_first_birth")
  }

  list(df = IRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 5: Child Nutrition & Feeding (KR / PR)
# ---------------------------------------------------------------------------
derive_child_nutrition_KR <- function(KRdata) {
  added <- character()
  if (is.null(KRdata) || nrow(KRdata) == 0) return(list(df = KRdata, added = added))

  # Exclusive breastfeeding (m4: 95 = still breastfeeding, others = not)
  # Only for children 0-5 months
  if (has_vars(KRdata, c("m4", "b19"))) {
    age_m <- as_num(KRdata$b19)  # current age in months
    bf_status <- as_num(KRdata$m4)
    # Currently breastfeeding = 95; exclusively BF means 95 AND no other foods
    # Simplified: m4==95 as proxy for still breastfeeding among 0-5 month olds
    IRdata_excl <- ifelse(is.na(age_m) | is.na(bf_status), NA_integer_,
                          ifelse(age_m <= 5, ifelse(bf_status == 95, 1L, 0L), NA_integer_))
    KRdata$c_excl_breastfed <- IRdata_excl
    added <- c(added, "c_excl_breastfed")
  }

  # Early initiation of breastfeeding (m34: time to breastfeeding after birth)
  # m34: 0=immediately, 100-123=hours, 200-223=days; <1 hour = immediately or within 1 hr
  if ("m34" %in% names(KRdata)) {
    x <- as_num(KRdata$m34)
    # immediately (0) or within 1 hour (100) = early initiation
    KRdata$c_early_initiation <- ifelse(is.na(x), NA_integer_,
                                        ifelse(x == 0 | x == 100, 1L,
                                               ifelse(x %in% c(999, 998, 997), NA_integer_, 0L)))
    added <- c(added, "c_early_initiation")
  }

  # Bottle feeding (m38: used bottle with nipple)
  if ("m38" %in% names(KRdata)) {
    x <- as_num(KRdata$m38)
    KRdata$c_bottle_fed <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "c_bottle_fed")
  }

  # Low birth weight (m19: birth weight in grams; 9996+ = not weighed)
  if ("m19" %in% names(KRdata)) {
    x <- as_num(KRdata$m19)
    x <- ifelse(x >= 500 & x < 9996, x, NA_real_)
    KRdata$c_low_birthweight <- ifelse(is.na(x), NA_integer_, ifelse(x < 2500, 1L, 0L))
    KRdata$c_mean_birthweight <- x  # continuous, in grams
    added <- c(added, "c_low_birthweight", "c_mean_birthweight")
  }

  # Small at birth (m18: 1=very large, 2=larger, 3=average, 4=smaller, 5=very small)
  if ("m18" %in% names(KRdata)) {
    x <- as_num(KRdata$m18)
    KRdata$c_small_at_birth <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 4:5, 1L, 0L))
    added <- c(added, "c_small_at_birth")
  }

  # Child diet indicators (v414 series: food groups consumed in last 24h)
  # These are for children 6-23 months
  iron_vars <- intersect(c("v414g", "v414h", "v414m"), names(KRdata))
  if (length(iron_vars) >= 1) {
    iron_mat <- sapply(iron_vars, function(v) as_num(KRdata[[v]]) == 1)
    KRdata$c_iron_rich_food <- ifelse(rowSums(iron_mat, na.rm = TRUE) > 0, 1L, 0L)
    KRdata$c_iron_rich_food[rowSums(!is.na(iron_mat)) == 0] <- NA_integer_
    added <- c(added, "c_iron_rich_food")
  }

  vita_vars <- intersect(c("v414i", "v414j", "v414k"), names(KRdata))
  if (length(vita_vars) >= 1) {
    vita_mat <- sapply(vita_vars, function(v) as_num(KRdata[[v]]) == 1)
    KRdata$c_vita_rich_food <- ifelse(rowSums(vita_mat, na.rm = TRUE) > 0, 1L, 0L)
    KRdata$c_vita_rich_food[rowSums(!is.na(vita_mat)) == 0] <- NA_integer_
    added <- c(added, "c_vita_rich_food")
  }

  # Minimum dietary diversity (4+ food groups out of 7, for children 6-23 months)
  all_food_vars <- intersect(paste0("v414", letters[1:23]), names(KRdata))
  if (length(all_food_vars) >= 4) {
    food_mat <- sapply(all_food_vars, function(v) as_num(KRdata[[v]]) == 1)
    food_count <- rowSums(food_mat, na.rm = TRUE)
    food_any <- rowSums(!is.na(food_mat)) > 0
    KRdata$c_diet_diversity_score <- ifelse(food_any, food_count, NA_real_)
    KRdata$c_min_diet_diversity   <- ifelse(food_any, ifelse(food_count >= 4, 1L, 0L), NA_integer_)
    added <- c(added, "c_diet_diversity_score", "c_min_diet_diversity")
  }

  list(df = KRdata, added = added)
}


derive_child_nutrition_PR <- function(PRdata) {
  added <- character()
  if (is.null(PRdata) || nrow(PRdata) == 0) return(list(df = PRdata, added = added))

  # Child anthropometry z-scores (hc70=HAZ, hc71=WAZ, hc72=WHZ; stored as z*100)
  zscore_defs <- list(
    c_stunted          = list(var = "hc70", threshold = -200, dir = "lt"),
    c_severely_stunted = list(var = "hc70", threshold = -300, dir = "lt"),
    c_wasted           = list(var = "hc72", threshold = -200, dir = "lt"),
    c_severely_wasted  = list(var = "hc72", threshold = -300, dir = "lt"),
    c_underweight      = list(var = "hc71", threshold = -200, dir = "lt"),
    c_overweight       = list(var = "hc72", threshold =  200, dir = "gt")
  )

  for (ind_name in names(zscore_defs)) {
    def <- zscore_defs[[ind_name]]
    if (def$var %in% names(PRdata)) {
      z <- as_num(PRdata[[def$var]])
      # DHS flags: 9996-9999 are missing/flagged
      z <- ifelse(z >= -600 & z <= 600, z, NA_real_)
      if (def$dir == "lt") {
        PRdata[[ind_name]] <- ifelse(is.na(z), NA_integer_, ifelse(z < def$threshold, 1L, 0L))
      } else {
        PRdata[[ind_name]] <- ifelse(is.na(z), NA_integer_, ifelse(z > def$threshold, 1L, 0L))
      }
      added <- c(added, ind_name)
    }
  }

  # Continuous z-score means
  zscore_means <- list(c_mean_haz = "hc70", c_mean_waz = "hc71", c_mean_whz = "hc72")
  for (nm in names(zscore_means)) {
    v <- zscore_means[[nm]]
    if (v %in% names(PRdata)) {
      z <- as_num(PRdata[[v]])
      PRdata[[nm]] <- ifelse(z >= -600 & z <= 600, z / 100, NA_real_)  # back to z-score units
      added <- c(added, nm)
    }
  }

  # Child anemia (hc56 = hemoglobin altitude-adjusted, g/dL * 10)
  if ("hc56" %in% names(PRdata)) {
    hb <- as_num(PRdata$hc56)
    hb_gdl <- ifelse(hb >= 10 & hb < 999, hb / 10, NA_real_)
    PRdata$c_anemia_any      <- ifelse(is.na(hb_gdl), NA_integer_, ifelse(hb_gdl < 11, 1L, 0L))
    PRdata$c_anemia_moderate <- ifelse(is.na(hb_gdl), NA_integer_, ifelse(hb_gdl < 10, 1L, 0L))
    PRdata$c_mean_hemoglobin <- hb_gdl  # continuous
    added <- c(added, "c_anemia_any", "c_anemia_moderate", "c_mean_hemoglobin")
  }

  list(df = PRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 6: Child Vaccination & Health (KR)
# ---------------------------------------------------------------------------
derive_child_health <- function(KRdata) {
  added <- character()
  if (is.null(KRdata) || nrow(KRdata) == 0) return(list(df = KRdata, added = added))

  # Vaccination variables: h2=BCG, h3=DPT1, h5=DPT2, h7=DPT3, h9=measles
  # Codes: 0=no, 1=vaccination date on card, 2=reported by mother, 3=vacc marked on card
  vacc_defs <- list(
    c_bcg     = "h2",
    c_dpt1    = "h3",
    c_dpt3    = "h7",
    c_polio3  = "h8",
    c_measles1 = "h9"
  )

  for (ind_name in names(vacc_defs)) {
    v <- vacc_defs[[ind_name]]
    if (v %in% names(KRdata)) {
      x <- as_num(KRdata[[v]])
      KRdata[[ind_name]] <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:3, 1L, 0L))
      added <- c(added, ind_name)
    }
  }

  # Fully vaccinated (BCG + DPT3 + Polio3 + Measles)
  basic_vacc <- intersect(c("c_bcg", "c_dpt3", "c_polio3", "c_measles1"), added)
  if (length(basic_vacc) >= 3) {
    vacc_mat <- as.matrix(KRdata[, basic_vacc, drop = FALSE])
    KRdata$c_fully_vaccinated <- ifelse(rowSums(!is.na(vacc_mat)) == length(basic_vacc),
                                        ifelse(rowSums(vacc_mat, na.rm = FALSE) == length(basic_vacc), 1L, 0L),
                                        NA_integer_)
    KRdata$c_no_vaccination   <- ifelse(rowSums(!is.na(vacc_mat)) == length(basic_vacc),
                                        ifelse(rowSums(vacc_mat, na.rm = FALSE) == 0, 1L, 0L),
                                        NA_integer_)
    added <- c(added, "c_fully_vaccinated", "c_no_vaccination")
  }

  # Vitamin A supplement (h33 or h34)
  for (v in c("h33", "h34")) {
    if (v %in% names(KRdata)) {
      x <- as_num(KRdata[[v]])
      KRdata$c_vita_supplement <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
      added <- c(added, "c_vita_supplement")
      break
    }
  }

  # Diarrhea in last 2 weeks (h11: 0=no, 1=yes in last 24h, 2=yes in last 2 weeks)
  if ("h11" %in% names(KRdata)) {
    x <- as_num(KRdata$h11)
    KRdata$c_diarrhea_2wk <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "c_diarrhea_2wk")
  }

  # Fever in last 2 weeks (h22)
  if ("h22" %in% names(KRdata)) {
    x <- as_num(KRdata$h22)
    KRdata$c_fever_2wk <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "c_fever_2wk")
  }

  # ARI symptoms (h31: cough + short rapid breaths)
  if ("h31" %in% names(KRdata)) {
    x <- as_num(KRdata$h31)
    KRdata$c_ari_2wk <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "c_ari_2wk")
  }

  # ORT for diarrhea (h13: ORS, h13b: RHS, h14: given increased fluids)
  ort_vars <- intersect(c("h13", "h13b", "h14"), names(KRdata))
  if (length(ort_vars) >= 1 && "h11" %in% names(KRdata)) {
    # Only among children with diarrhea
    had_diar <- as_num(KRdata$h11) %in% 1:2
    ort_mat <- sapply(ort_vars, function(v) as_num(KRdata[[v]]) == 1)
    got_ort <- rowSums(ort_mat, na.rm = TRUE) > 0
    KRdata$c_ort_diarrhea <- ifelse(!had_diar, NA_integer_,
                                    ifelse(got_ort, 1L, 0L))
    added <- c(added, "c_ort_diarrhea")
  }

  # Zinc for diarrhea (h15e)
  if ("h15e" %in% names(KRdata) && "h11" %in% names(KRdata)) {
    had_diar <- as_num(KRdata$h11) %in% 1:2
    x <- as_num(KRdata$h15e)
    KRdata$c_zinc_diarrhea <- ifelse(!had_diar, NA_integer_,
                                     ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L)))
    added <- c(added, "c_zinc_diarrhea")
  }

  # Care seeking for fever (h32: sought advice/treatment)
  if ("h32" %in% names(KRdata) && "h22" %in% names(KRdata)) {
    had_fever <- as_num(KRdata$h22) == 1
    x <- as_num(KRdata$h32)
    KRdata$c_care_seeking_fever <- ifelse(!had_fever, NA_integer_,
                                          ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L)))
    added <- c(added, "c_care_seeking_fever")
  }

  list(df = KRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 7: WASH (HR)
# ---------------------------------------------------------------------------
derive_wash <- function(HRdata) {
  added <- character()
  if (is.null(HRdata) || nrow(HRdata) == 0) return(list(df = HRdata, added = added))

  # Improved water source (hv201)
  if ("hv201" %in% names(HRdata)) {
    x <- as_num(HRdata$hv201)
    improved_water <- c(11, 12, 13, 14, 21, 31, 32, 41, 51, 61, 71)
    # More conservative improved set:
    improved_water_strict <- c(11, 12, 13, 14, 31, 32, 41)
    HRdata$hh_improved_water <- ifelse(is.na(x), NA_integer_,
                                       ifelse(x %in% improved_water_strict, 1L, 0L))
    added <- c(added, "hh_improved_water")
  }

  # Improved sanitation (hv205 + hv225)
  if ("hv205" %in% names(HRdata)) {
    tlet <- as_num(HRdata$hv205)
    improved_san <- c(11, 12, 13, 21, 22, 41)
    imp <- ifelse(is.na(tlet), NA_integer_, ifelse(tlet %in% improved_san, 1L, 0L))

    if ("hv225" %in% names(HRdata)) {
      shared <- as_num(HRdata$hv225)
      HRdata$hh_improved_sanitation <- ifelse(is.na(imp) | is.na(shared), NA_integer_,
                                              ifelse(imp == 1L & shared == 0, 1L, 0L))
    } else {
      HRdata$hh_improved_sanitation <- imp
    }
    added <- c(added, "hh_improved_sanitation")

    # Open defecation
    HRdata$hh_open_defecation <- ifelse(is.na(tlet), NA_integer_,
                                        ifelse(tlet %in% c(31, 61, 96), 1L, 0L))
    added <- c(added, "hh_open_defecation")
  }

  # Water on premises / <=30 min (hv204 = time to water source in minutes)
  if ("hv204" %in% names(HRdata)) {
    t <- as_num(HRdata$hv204)
    t <- ifelse(t >= 0 & t < 996, t, NA_real_)
    HRdata$hh_water_onpremise <- ifelse(is.na(t), NA_integer_, ifelse(t == 0, 1L, 0L))
    HRdata$hh_water_30min     <- ifelse(is.na(t), NA_integer_, ifelse(t <= 30, 1L, 0L))
    added <- c(added, "hh_water_onpremise", "hh_water_30min")
  }

  # Handwashing facility (hv230a: 1=observed fixed, 2=observed mobile, 3=not observed)
  if ("hv230a" %in% names(HRdata)) {
    x <- as_num(HRdata$hv230a)
    HRdata$hh_handwash_facility <- ifelse(is.na(x), NA_integer_, ifelse(x %in% 1:2, 1L, 0L))
    added <- c(added, "hh_handwash_facility")
  }

  # Soap available (hv232: 1=soap/detergent, 0=none)
  if ("hv232" %in% names(HRdata)) {
    x <- as_num(HRdata$hv232)
    HRdata$hh_soap_available <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "hh_soap_available")
  }

  list(df = HRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 8: Household Socioeconomic (HR)
# ---------------------------------------------------------------------------
derive_household_ses <- function(HRdata) {
  added <- character()
  if (is.null(HRdata) || nrow(HRdata) == 0) return(list(df = HRdata, added = added))

  # Wealth quintile (hv270: 1=poorest ... 5=richest)
  if ("hv270" %in% names(HRdata)) {
    x <- as_num(HRdata$hv270)
    HRdata$hh_wealth_poorest <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    HRdata$hh_wealth_poor    <- ifelse(is.na(x), NA_integer_, ifelse(x <= 2, 1L, 0L))
    added <- c(added, "hh_wealth_poorest", "hh_wealth_poor")
  }

  # Wealth index score (hv271: continuous, factor score * 100000)
  if ("hv271" %in% names(HRdata)) {
    x <- as_num(HRdata$hv271)
    HRdata$hh_wealth_mean <- x / 100000  # rescale to original factor score
    added <- c(added, "hh_wealth_mean")
  }

  # Electricity (hv206: 0=no, 1=yes)
  if ("hv206" %in% names(HRdata)) {
    x <- as_num(HRdata$hv206)
    HRdata$hh_electricity <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "hh_electricity")
  }

  # Floor material (hv213: improved = not earth/sand/dung; codes vary but 30+ typically finished)
  if ("hv213" %in% names(HRdata)) {
    x <- as_num(HRdata$hv213)
    # Natural/rudimentary floors: 11-19 (earth, sand, dung), 21-29 (rudimentary)
    # Finished floors: 31+ (parquet, vinyl, ceramic, cement, carpet, etc.)
    HRdata$hh_improved_floor <- ifelse(is.na(x), NA_integer_, ifelse(x >= 30, 1L, 0L))
    added <- c(added, "hh_improved_floor")
  }

  # Urban residence (hv025: 1=urban, 2=rural)
  if ("hv025" %in% names(HRdata)) {
    x <- as_num(HRdata$hv025)
    HRdata$hh_urban <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "hh_urban")
  }

  # Household assets
  asset_vars <- list(
    hh_radio  = "hv207",
    hh_tv     = "hv208",
    hh_fridge = "hv209"
  )
  for (nm in names(asset_vars)) {
    v <- asset_vars[[nm]]
    if (v %in% names(HRdata)) {
      x <- as_num(HRdata[[v]])
      HRdata[[nm]] <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
      added <- c(added, nm)
    }
  }

  # Household size (hv009 = number of household members)
  if ("hv009" %in% names(HRdata)) {
    x <- as_num(HRdata$hv009)
    HRdata$hh_mean_hhsize <- ifelse(x >= 1 & x <= 50, x, NA_real_)
    added <- c(added, "hh_mean_hhsize")
  }

  # Crowding (hv009 / hv216 = sleeping rooms)
  if (has_vars(HRdata, c("hv009", "hv216"))) {
    members <- as_num(HRdata$hv009)
    rooms <- as_num(HRdata$hv216)
    rooms <- ifelse(rooms >= 1 & rooms <= 20, rooms, NA_real_)
    ratio <- members / rooms
    HRdata$hh_crowding <- ifelse(is.na(ratio), NA_integer_, ifelse(ratio > 3, 1L, 0L))
    added <- c(added, "hh_crowding")
  }

  list(df = HRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 9: Malaria & Vector Control (HR/PR)
# ---------------------------------------------------------------------------
derive_malaria_HR <- function(HRdata) {
  added <- character()
  if (is.null(HRdata) || nrow(HRdata) == 0) return(list(df = HRdata, added = added))

  # Any ITN in household (hv227)
  if ("hv227" %in% names(HRdata)) {
    x <- as_num(HRdata$hv227)
    HRdata$hh_itn_any <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "hh_itn_any")
  }

  # IRS in last 12 months (hv253)
  if ("hv253" %in% names(HRdata)) {
    x <- as_num(HRdata$hv253)
    HRdata$hh_irs <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "hh_irs")
  }

  list(df = HRdata, added = added)
}


derive_malaria_PR <- function(PRdata) {
  added <- character()
  if (is.null(PRdata) || nrow(PRdata) == 0) return(list(df = PRdata, added = added))

  # Child <5 slept under ITN (hml12)
  if (has_vars(PRdata, c("hv105", "hml12"))) {
    age <- as_num(PRdata$hv105)
    used <- as_num(PRdata$hml12)
    PRdata$c_slept_itn <- ifelse(is.na(age), NA_integer_,
                                 ifelse(age < 5,
                                        ifelse(is.na(used), NA_integer_, ifelse(used == 1, 1L, 0L)),
                                        NA_integer_))
    added <- c(added, "c_slept_itn")
  }

  # Malaria RDT positive (hml35: 0=neg, 1=pos)
  if ("hml35" %in% names(PRdata)) {
    x <- as_num(PRdata$hml35)
    PRdata$c_malaria_rdt <- ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L))
    added <- c(added, "c_malaria_rdt")
  }

  list(df = PRdata, added = added)
}


# ---------------------------------------------------------------------------
# Domain 10: Mortality & Birth Outcomes (BR)
# ---------------------------------------------------------------------------
derive_mortality <- function(BRdata) {
  added <- character()
  if (is.null(BRdata) || nrow(BRdata) == 0) return(list(df = BRdata, added = added))

  # Neonatal death (b7 = age at death in months; b5 = child alive)
  if (has_vars(BRdata, c("b5", "b7"))) {
    alive <- as_num(BRdata$b5)
    age_death <- as_num(BRdata$b7)

    # b5: 0 = dead, 1 = alive
    BRdata$c_neonatal_death <- ifelse(is.na(alive), NA_integer_,
                                      ifelse(alive == 0 & !is.na(age_death) & age_death == 0, 1L, 0L))
    BRdata$c_infant_death   <- ifelse(is.na(alive), NA_integer_,
                                      ifelse(alive == 0 & !is.na(age_death) & age_death < 12, 1L, 0L))
    BRdata$c_under5_death   <- ifelse(is.na(alive), NA_integer_,
                                      ifelse(alive == 0 & !is.na(age_death) & age_death < 60, 1L, 0L))
    added <- c(added, "c_neonatal_death", "c_infant_death", "c_under5_death")
  }

  list(df = BRdata, added = added)
}


# =============================================================================
# 3. Master function: derive all indicators from one country-year's microdata
# =============================================================================

derive_all_indicators <- function(dhs_data_list) {
  # dhs_data_list should have named elements: IRdata, KRdata, PRdata, HRdata, BRdata

  all_standardized <- list()  # named list of indicator_name -> standardized data frame
  log <- tibble(indicator = character(), recode = character(), n = integer(), ok = logical())

  # --- IR (Individual Recode / Women's data) ---
  IRdata <- dhs_data_list$IRdata
  if (!is.null(IRdata) && nrow(IRdata) > 0) {
    res1 <- derive_maternal_nutrition(IRdata)
    res2 <- derive_maternal_diet(res1$df)
    res3 <- derive_reproductive_health(res2$df)
    res4 <- derive_maternal_education(res3$df)
    IRdata <- res4$df
    ir_added <- c(res1$added, res2$added, res3$added, res4$added)

    for (ind in ir_added) {
      std <- standardize_indicator(IRdata, "IR", ind)
      if (!is.null(std) && nrow(std) > 0) {
        all_standardized[[ind]] <- std
        log <- bind_rows(log, tibble(indicator = ind, recode = "IR", n = nrow(std), ok = TRUE))
      } else {
        log <- bind_rows(log, tibble(indicator = ind, recode = "IR", n = 0L, ok = FALSE))
      }
    }
    cat(sprintf("    IR: %d indicators derived (%d standardized)\n",
                length(ir_added), sum(ir_added %in% names(all_standardized))))
  }

  # --- KR (Children's Recode) ---
  KRdata <- dhs_data_list$KRdata
  if (!is.null(KRdata) && nrow(KRdata) > 0) {
    res5 <- derive_child_nutrition_KR(KRdata)
    res6 <- derive_child_health(res5$df)
    KRdata <- res6$df
    kr_added <- c(res5$added, res6$added)

    for (ind in kr_added) {
      std <- standardize_indicator(KRdata, "KR", ind)
      if (!is.null(std) && nrow(std) > 0) {
        all_standardized[[ind]] <- std
        log <- bind_rows(log, tibble(indicator = ind, recode = "KR", n = nrow(std), ok = TRUE))
      } else {
        log <- bind_rows(log, tibble(indicator = ind, recode = "KR", n = 0L, ok = FALSE))
      }
    }
    cat(sprintf("    KR: %d indicators derived (%d standardized)\n",
                length(kr_added), sum(kr_added %in% names(all_standardized))))
  }

  # --- PR (Household Member Recode) ---
  PRdata <- dhs_data_list$PRdata
  if (!is.null(PRdata) && nrow(PRdata) > 0) {
    res7 <- derive_child_nutrition_PR(PRdata)
    res8 <- derive_malaria_PR(res7$df)
    PRdata <- res8$df
    pr_added <- c(res7$added, res8$added)

    for (ind in pr_added) {
      std <- standardize_indicator(PRdata, "PR", ind)
      if (!is.null(std) && nrow(std) > 0) {
        all_standardized[[ind]] <- std
        log <- bind_rows(log, tibble(indicator = ind, recode = "PR", n = nrow(std), ok = TRUE))
      } else {
        log <- bind_rows(log, tibble(indicator = ind, recode = "PR", n = 0L, ok = FALSE))
      }
    }
    cat(sprintf("    PR: %d indicators derived (%d standardized)\n",
                length(pr_added), sum(pr_added %in% names(all_standardized))))
  }

  # --- HR (Household Recode) ---
  HRdata <- dhs_data_list$HRdata
  if (!is.null(HRdata) && nrow(HRdata) > 0) {
    res9  <- derive_wash(HRdata)
    res10 <- derive_household_ses(res9$df)
    res11 <- derive_malaria_HR(res10$df)
    HRdata <- res11$df
    hr_added <- c(res9$added, res10$added, res11$added)

    for (ind in hr_added) {
      std <- standardize_indicator(HRdata, "HR", ind)
      if (!is.null(std) && nrow(std) > 0) {
        all_standardized[[ind]] <- std
        log <- bind_rows(log, tibble(indicator = ind, recode = "HR", n = nrow(std), ok = TRUE))
      } else {
        log <- bind_rows(log, tibble(indicator = ind, recode = "HR", n = 0L, ok = FALSE))
      }
    }
    cat(sprintf("    HR: %d indicators derived (%d standardized)\n",
                length(hr_added), sum(hr_added %in% names(all_standardized))))
  }

  # --- BR (Births Recode) ---
  BRdata <- dhs_data_list$BRdata
  if (!is.null(BRdata) && nrow(BRdata) > 0) {
    res12 <- derive_mortality(BRdata)
    BRdata <- res12$df
    br_added <- res12$added

    for (ind in br_added) {
      std <- standardize_indicator(BRdata, "BR", ind)
      if (!is.null(std) && nrow(std) > 0) {
        all_standardized[[ind]] <- std
        log <- bind_rows(log, tibble(indicator = ind, recode = "BR", n = nrow(std), ok = TRUE))
      } else {
        log <- bind_rows(log, tibble(indicator = ind, recode = "BR", n = 0L, ok = FALSE))
      }
    }
    cat(sprintf("    BR: %d indicators derived (%d standardized)\n",
                length(br_added), sum(br_added %in% names(all_standardized))))
  }

  cat(sprintf("  Total: %d custom indicators standardized\n", length(all_standardized)))

  list(standardized = all_standardized, log = log)
}


# =============================================================================
# 4. Compute admin-2 direct estimates for all custom indicators
# =============================================================================

compute_custom_admin2 <- function(standardized_list, cluster.info) {

  results <- list()

  for (ind_name in names(standardized_list)) {
    dat <- standardized_list[[ind_name]]

    tryCatch({
      est <- directEST(data = dat, cluster.info = cluster.info, admin = 2, aggregation = FALSE)
      n_est <- sum(!is.na(est$res.admin2$direct.est))
      cat(sprintf("    [%s] %d areas\n", ind_name, n_est))
      results[[ind_name]] <- est$res.admin2
    }, error = function(e) {
      cat(sprintf("    [%s] FAILED: %s\n", ind_name, e$message))
    })
  }

  results
}


# =============================================================================
# 5. Build wide output
# =============================================================================

build_custom_wide <- function(admin2_results, prefix) {
  if (length(admin2_results) == 0) return(tibble())

  rows <- list()
  for (ind_name in names(admin2_results)) {
    df <- admin2_results[[ind_name]]
    if (!is.null(df) && "admin2.name.full" %in% names(df) && "direct.est" %in% names(df)) {
      rows[[ind_name]] <- tibble(
        admin2.name.full = df$admin2.name.full,
        indicator = ind_name,
        est = df$direct.est
      )
    }
  }

  if (length(rows) == 0) return(tibble())

  long <- bind_rows(rows)

  wide <- long %>%
    pivot_wider(
      names_from   = indicator,
      values_from  = est,
      names_prefix = prefix
    )

  # Add _adm2 suffix to all indicator columns
  wide <- wide %>%
    rename_with(~ paste0(.x, "_adm2"), -admin2.name.full)

  # Extract Admin2 name
  wide$Admin2 <- sub("^[^_]+_", "", wide$admin2.name.full)
  wide$admin2.name.full <- NULL
  wide <- wide %>% select(Admin2, everything())

  wide
}


# =============================================================================
# 6. Main loop
# =============================================================================

cat("=== DHS Custom Admin-2 Indicator Aggregation ===\n\n")

for (entry in COUNTRY_REGISTRY) {

  if (!is.null(COUNTRIES_TO_RUN) && !(entry$country %in% COUNTRIES_TO_RUN)) next

  for (yr in entry$years) {

    cat(sprintf("\n%s\n  %s %d\n%s\n\n",
                strrep("=", 70), entry$country, yr, strrep("=", 70)))

    prefix <- paste0("dhs", yr, "_")

    # --- Step 1: Load or download raw microdata ---
    cat("  [1/4] Loading raw DHS microdata...\n")

    # Try loading cached data first
    cached_file <- list.files(
      here("data", "DHS"),
      pattern = paste0("dhs_", entry$country, "_", yr, "\\.RDS$"),
      full.names = TRUE
    )

    dhs_raw <- NULL
    if (length(cached_file) > 0) {
      cat(sprintf("    Loading cached: %s\n", basename(cached_file[1])))
      dhs_loaded <- tryCatch(readRDS(cached_file[1]), error = function(e) NULL)
      if (!is.null(dhs_loaded)) {
        # Extract the inner list (files are saved as list(country_year = list(...)))
        inner_name <- paste0(entry$country, "_", yr)
        if (inner_name %in% names(dhs_loaded)) {
          dhs_raw <- dhs_loaded[[inner_name]]
        } else {
          dhs_raw <- dhs_loaded
        }
      }
    }

    # If not cached, download
    if (is.null(dhs_raw)) {
      cat("    Downloading from DHS...\n")
      dhs_raw <- tryCatch(
        getDHSdata(country = entry$country, indicator = NULL, year = yr),
        error = function(e) {
          cat(sprintf("    FATAL: Download failed: %s\n", e$message))
          NULL
        }
      )
    }

    if (is.null(dhs_raw)) {
      cat(sprintf("  Skipping %s %d — no data available\n", entry$country, yr))
      next
    }

    # Identify recode datasets
    recode_names <- names(dhs_raw)
    cat(sprintf("    Available recodes: %s\n", paste(recode_names, collapse = ", ")))

    # Map to standard names
    dhs_data_list <- list(
      IRdata = dhs_raw[["IRdata"]],
      KRdata = dhs_raw[["KRdata"]],
      PRdata = dhs_raw[["PRdata"]],
      HRdata = dhs_raw[["HRdata"]],
      BRdata = dhs_raw[["BRdata"]]
    )

    # --- Step 2: Derive custom indicators ---
    cat("  [2/4] Deriving custom indicators from raw microdata...\n")
    custom_res <- tryCatch(
      derive_all_indicators(dhs_data_list),
      error = function(e) {
        cat(sprintf("    FATAL: Indicator derivation failed: %s\n", e$message))
        NULL
      }
    )

    if (is.null(custom_res) || length(custom_res$standardized) == 0) {
      cat(sprintf("  Skipping %s %d — no custom indicators derived\n", entry$country, yr))
      next
    }

    # --- Step 3: Build spatial info ---
    cat("  [3/4] Building spatial info (GPS + GADM boundaries)...\n")
    spatial <- tryCatch({
      geo <- getDHSgeo(country = entry$country, year = yr)
      cat(sprintf("    GPS clusters: %d\n", nrow(geo)))

      source(here::here("R", "data_prep.R"))
      poly.adm1 <- load_gadm_cached(entry$iso3, level = 1)
      poly.adm2 <- load_gadm_cached(entry$iso3, level = 2)

      cat(sprintf("    Admin1: %d | Admin2: %d polygons\n", nrow(poly.adm1), nrow(poly.adm2)))

      cluster.info <- clusterInfo(
        geo = geo, poly.adm1 = poly.adm1, poly.adm2 = poly.adm2,
        by.adm1 = "NAME_1", by.adm2 = "NAME_2"
      )

      n_wrong <- length(cluster.info$wrong.points)
      if (n_wrong > 0) {
        cat(sprintf("    WARNING: %d clusters with invalid GPS excluded\n", n_wrong))
      }

      list(cluster.info = cluster.info, poly.adm1 = poly.adm1, poly.adm2 = poly.adm2)
    }, error = function(e) {
      cat(sprintf("    FATAL: Spatial info failed: %s\n", e$message))
      NULL
    })

    if (is.null(spatial)) next

    # --- Step 4: Compute direct admin-2 estimates ---
    cat("  [4/4] Computing direct admin-2 estimates...\n")

    admin2_results <- tryCatch(
      compute_custom_admin2(custom_res$standardized, spatial$cluster.info),
      error = function(e) {
        cat(sprintf("    FATAL: Admin-2 estimation failed: %s\n", e$message))
        list()
      }
    )

    # --- Step 5: Compile and save ---
    tryCatch({
      wide <- build_custom_wide(admin2_results, prefix)

      file_prefix <- paste0(entry$country, "_", yr)

      if (nrow(wide) > 0) {
        f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_custom_admin2_wide.rds"))
        saveRDS(wide, f)
        cat(sprintf("    Saved: %s (%d admin2 areas, %d indicator columns)\n",
                    basename(f), nrow(wide), ncol(wide) - 1L))
      }

      # Save the derivation log
      f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_custom_admin2_log.rds"))
      saveRDS(custom_res$log, f)
      cat(sprintf("    Saved: %s\n", basename(f)))

      cat(sprintf("\n  %s %d complete: %d custom indicators at admin-2\n",
                  entry$country, yr, length(admin2_results)))

    }, error = function(e) {
      cat(sprintf("  ERROR in compile/save for %s %d: %s\n", entry$country, yr, e$message))
    })
  }
}

cat("\n=== Custom admin-2 aggregation complete ===\n")
