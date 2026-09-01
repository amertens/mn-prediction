# =============================================================================
# src/DHS/DHS_extra_indicators_2026-09.R
#
# Additional DHS-derived Admin-2 indicators, covering biological, dietary and
# health-access pathways the existing 97 columns do not reach.
#
# WHY THESE, AND WHY NOW
# ----------------------
# The harmonized predictor set carries 97 DHS columns: 20 surveyPrev built-ins
# and 77 custom derivations covering vaccination, WASH, assets, education and
# fertility. Screening 294 predictors found two survivors, and a targeted scan
# of 31 nutrition-proximal predictors found none. Both of those scans were run
# on covariates that are DISTAL to intake: soil, climate, crop area, programme
# coverage.
#
# The DHS carries something none of them do -- what people actually ate in the
# last 24 hours, and what animals the household owns. Those are the closest
# thing to dietary intake obtainable without a survey of the target country,
# because they can be aggregated from a DHS the country already ran.
#
# Eight blocks, chosen for a mechanism rather than for availability, and each
# verified present in all four countries' cached recodes before being written:
#
#   child dietary diversity     23 food-group variables, identical in all four
#   maternal dietary recall     the same instrument for women
#   livestock ownership         14 variables; animal-source foods are the main
#                               dietary route for iron, zinc and B12
#   supplementation in pregnancy iron duration, IPTp, deworming
#   child preventive contact    vitamin A capsule, deworming
#   cooking fuel                 inflammation pathway via household air
#   health-access barriers      distance, money, permission, insurance
#   women's work and earnings   control over household resources
#
# THE LEAK GUARD APPLIES HERE AND IS NOT OPTIONAL.
# Anaemia and haemoglobin (hw53, v453, ha53, hc53) and malaria test results
# (hml32-hml36) are DELIBERATELY EXCLUDED. They are measurements from a blood
# draw, several are the outcome itself under another name, and DHS names them
# by questionnaire position so no analyte-name guard catches them. See
# dhs_measurement_class() in R/data_prep.R.
#
# Sourced alongside src/DHS/DHS_custom_admin2_indicators.R; the derive_*
# contract is the same, returning list(df, added).
# =============================================================================

# `as_num` and `has_vars` come from DHS_custom_admin2_indicators.R.

# ---------------------------------------------------------------------------
# Child dietary diversity, from the 24-hour recall food groups (KR).
#
# v414* are the individual food groups. The WHO minimum-dietary-diversity
# indicator counts groups out of eight; the count itself is more informative
# than the binary, so both are emitted.
# ---------------------------------------------------------------------------
derive_child_diet_groups <- function(KRdata) {
  added <- character()
  if (is.null(KRdata) || nrow(KRdata) == 0) return(list(df = KRdata, added = added))

  fg <- list(
    grains      = c("v412a", "v414e", "v414f"),
    roots       = c("v414f"),
    legumes     = c("v414o"),
    dairy       = c("v411", "v411a", "v414p"),
    flesh       = c("v414h", "v414m", "v414n"),
    eggs        = c("v414g"),
    vitA_fruitveg = c("v414i", "v414j", "v414k"),
    other_fruitveg = c("v414l")
  )
  present <- 0L
  cnt <- rep(0, nrow(KRdata)); denom <- rep(0, nrow(KRdata))
  for (g in names(fg)) {
    vv <- intersect(fg[[g]], names(KRdata))
    if (!length(vv)) next
    present <- present + 1L
    m <- do.call(cbind, lapply(vv, function(v) {
      x <- as_num(KRdata[[v]]); ifelse(is.na(x), NA, ifelse(x == 1, 1, 0))
    }))
    any1 <- apply(m, 1, function(z) if (all(is.na(z))) NA else as.integer(any(z == 1, na.rm = TRUE)))
    nm <- paste0("c_fg_", g)
    KRdata[[nm]] <- any1
    added <- c(added, nm)
    cnt <- cnt + ifelse(is.na(any1), 0, any1)
    denom <- denom + ifelse(is.na(any1), 0, 1)
  }
  if (present >= 4L) {
    KRdata$c_diet_groups_n <- ifelse(denom >= 4, cnt, NA_real_)
    KRdata$c_mdd_4plus <- ifelse(denom >= 4, as.integer(cnt >= 4), NA_integer_)
    added <- c(added, "c_diet_groups_n", "c_mdd_4plus")
  }
  list(df = KRdata, added = added)
}

# ---------------------------------------------------------------------------
# Livestock and land, from the household recode.
#
# hv246a-h are counts of cattle, cows, horses, goats, sheep, chickens. Animal
# ownership is the household's own route to animal-source foods, which carry
# haem iron, zinc and B12 in forms plant sources do not.
# ---------------------------------------------------------------------------
derive_livestock <- function(HRdata) {
  added <- character()
  if (is.null(HRdata) || nrow(HRdata) == 0) return(list(df = HRdata, added = added))

  spec <- c(hh_cattle = "hv246a", hh_cows = "hv246b", hh_horses = "hv246c",
            hh_goats = "hv246d", hh_sheep = "hv246e", hh_chickens = "hv246f")
  owned <- NULL
  for (nm in names(spec)) {
    v <- spec[[nm]]
    if (!v %in% names(HRdata)) next
    x <- as_num(HRdata[[v]])
    # DHS codes 95+ as "95 or more" and 98/99 as unknown.
    x[x >= 96] <- NA
    HRdata[[nm]] <- x
    HRdata[[paste0(nm, "_any")]] <- ifelse(is.na(x), NA_integer_, as.integer(x > 0))
    added <- c(added, nm, paste0(nm, "_any"))
    owned <- if (is.null(owned)) ifelse(is.na(x), 0, x) else owned + ifelse(is.na(x), 0, x)
  }
  if (!is.null(owned)) {
    HRdata$hh_livestock_total <- owned
    HRdata$hh_livestock_any <- as.integer(owned > 0)
    added <- c(added, "hh_livestock_total", "hh_livestock_any")
  }
  if ("hv245" %in% names(HRdata)) {
    x <- as_num(HRdata$hv245); x[x >= 95] <- NA
    HRdata$hh_agland_ha <- x
    HRdata$hh_agland_any <- ifelse(is.na(x), NA_integer_, as.integer(x > 0))
    added <- c(added, "hh_agland_ha", "hh_agland_any")
  }
  if ("hv244" %in% names(HRdata)) {
    x <- as_num(HRdata$hv244)
    HRdata$hh_owns_land <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "hh_owns_land")
  }
  list(df = HRdata, added = added)
}

# ---------------------------------------------------------------------------
# Micronutrient supplementation and preventive contact during pregnancy (IR).
# ---------------------------------------------------------------------------
derive_pregnancy_supplementation <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  # m45_1: took iron tablets/syrup during last pregnancy.
  if ("m45_1" %in% names(IRdata)) {
    x <- as_num(IRdata$m45_1)
    IRdata$w_iron_pregnancy <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "w_iron_pregnancy")
  }
  # m46_1: number of days iron was taken. 998 is "don't know".
  if ("m46_1" %in% names(IRdata)) {
    x <- as_num(IRdata$m46_1); x[x >= 998] <- NA
    IRdata$w_iron_days <- x
    IRdata$w_iron_90plus <- ifelse(is.na(x), NA_integer_, as.integer(x >= 90))
    added <- c(added, "w_iron_days", "w_iron_90plus")
  }
  # m49a: took intermittent preventive treatment for malaria in pregnancy.
  if ("m49a" %in% names(IRdata)) {
    x <- as_num(IRdata$m49a)
    IRdata$w_iptp_any <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "w_iptp_any")
  }
  # m60_1: took deworming medication during pregnancy.
  if ("m60_1" %in% names(IRdata)) {
    x <- as_num(IRdata$m60_1)
    IRdata$w_deworm_pregnancy <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "w_deworm_pregnancy")
  }
  list(df = IRdata, added = added)
}

# ---------------------------------------------------------------------------
# Child preventive contact (KR): vitamin A capsule and deworming.
# ---------------------------------------------------------------------------
derive_child_preventive <- function(KRdata) {
  added <- character()
  if (is.null(KRdata) || nrow(KRdata) == 0) return(list(df = KRdata, added = added))

  if ("h33" %in% names(KRdata)) {
    x <- as_num(KRdata$h33)
    # 0 = no, 1/2/3 = yes by various sources, 8 = don't know
    KRdata$c_vita_capsule <- ifelse(is.na(x) | x == 8, NA_integer_,
                                    as.integer(x %in% c(1, 2, 3)))
    added <- c(added, "c_vita_capsule")
  }
  if ("h43" %in% names(KRdata)) {
    x <- as_num(KRdata$h43)
    KRdata$c_deworm_6mo <- ifelse(is.na(x) | x == 8, NA_integer_, as.integer(x == 1))
    added <- c(added, "c_deworm_6mo")
  }
  list(df = KRdata, added = added)
}

# ---------------------------------------------------------------------------
# Cooking fuel (HR). Solid fuel use is the household-air-pollution pathway,
# which acts on micronutrient status through chronic inflammation rather than
# through intake.
# ---------------------------------------------------------------------------
derive_cooking_fuel <- function(HRdata) {
  added <- character()
  if (is.null(HRdata) || nrow(HRdata) == 0) return(list(df = HRdata, added = added))
  if ("hv226" %in% names(HRdata)) {
    x <- as_num(HRdata$hv226)
    # 1-4 clean (electricity, LPG, natural gas, biogas); 5+ solid/kerosene;
    # 95/96/99 other or missing.
    HRdata$hh_clean_fuel <- ifelse(is.na(x) | x >= 95, NA_integer_, as.integer(x <= 4))
    HRdata$hh_solid_fuel <- ifelse(is.na(x) | x >= 95, NA_integer_,
                                   as.integer(x >= 6 & x <= 11))
    added <- c(added, "hh_clean_fuel", "hh_solid_fuel")
  }
  if ("hv239" %in% names(HRdata)) {
    x <- as_num(HRdata$hv239)
    HRdata$hh_cook_indoors <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "hh_cook_indoors")
  }
  list(df = HRdata, added = added)
}

# ---------------------------------------------------------------------------
# Barriers to health care (IR). v467b-f are "big problem" flags.
# ---------------------------------------------------------------------------
derive_health_access <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  spec <- c(w_barrier_money = "v467c", w_barrier_distance = "v467d",
            w_barrier_transport = "v467e", w_barrier_alone = "v467f",
            w_barrier_permission = "v467b")
  n_bar <- NULL; d_bar <- NULL
  for (nm in names(spec)) {
    v <- spec[[nm]]
    if (!v %in% names(IRdata)) next
    x <- as_num(IRdata[[v]])
    b <- ifelse(is.na(x) | x >= 8, NA_integer_, as.integer(x == 1))
    IRdata[[nm]] <- b
    added <- c(added, nm)
    n_bar <- if (is.null(n_bar)) ifelse(is.na(b), 0, b) else n_bar + ifelse(is.na(b), 0, b)
    d_bar <- if (is.null(d_bar)) ifelse(is.na(b), 0, 1) else d_bar + ifelse(is.na(b), 0, 1)
  }
  if (!is.null(n_bar)) {
    IRdata$w_barriers_n <- ifelse(d_bar >= 3, n_bar, NA_real_)
    added <- c(added, "w_barriers_n")
  }
  if ("v481" %in% names(IRdata)) {
    x <- as_num(IRdata$v481)
    IRdata$w_health_insurance <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "w_health_insurance")
  }
  list(df = IRdata, added = added)
}

# ---------------------------------------------------------------------------
# Women's work and control over earnings (IR).
# ---------------------------------------------------------------------------
derive_womens_work <- function(IRdata) {
  added <- character()
  if (is.null(IRdata) || nrow(IRdata) == 0) return(list(df = IRdata, added = added))

  if ("v714" %in% names(IRdata)) {
    x <- as_num(IRdata$v714)
    IRdata$w_working <- ifelse(is.na(x), NA_integer_, as.integer(x == 1))
    added <- c(added, "w_working")
  }
  if ("v717" %in% names(IRdata)) {
    x <- as_num(IRdata$v717)
    IRdata$w_occ_agric <- ifelse(is.na(x), NA_integer_, as.integer(x %in% c(4, 5)))
    added <- c(added, "w_occ_agric")
  }
  if ("v739" %in% names(IRdata)) {
    x <- as_num(IRdata$v739)
    IRdata$w_decides_earnings <- ifelse(is.na(x) | x >= 8, NA_integer_,
                                        as.integer(x %in% c(1, 2)))
    added <- c(added, "w_decides_earnings")
  }
  if ("v745a" %in% names(IRdata)) {
    x <- as_num(IRdata$v745a)
    IRdata$w_owns_house <- ifelse(is.na(x), NA_integer_, as.integer(x > 0))
    added <- c(added, "w_owns_house")
  }
  list(df = IRdata, added = added)
}

# ---------------------------------------------------------------------------
# The registry, mirroring derive_all_indicators()'s contract so the caller can
# fold these in without knowing which file they came from.
# ---------------------------------------------------------------------------
EXTRA_DERIVERS <- list(
  IR = list(derive_pregnancy_supplementation, derive_health_access,
            derive_womens_work),
  KR = list(derive_child_diet_groups, derive_child_preventive),
  HR = list(derive_livestock, derive_cooking_fuel)
)
