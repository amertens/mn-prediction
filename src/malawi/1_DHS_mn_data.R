

library(tidyverse)
library(haven)
library(here)
source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

mw_psc <- read_dta(here("data/Malawi/MW_PSC.DTA")) %>% rename(svy_weight=mweight) #preschool children
mw_sac <- read_dta(here("data/Malawi/MW_SAC.DTA")) %>% rename(svy_weight=mweight) #school-age children
mw_wra <- read_dta(here("data/Malawi/MW_WRA.DTA")) %>% rename(svy_weight=mweight) #women
mw_men <- read_dta(here("data/Malawi/MW_MEN.DTA")) %>% rename(svy_weight=mweight) #men


gambia_variables <- makeVlist(mw_wra)
gambia_variables[grep("weight", gambia_variables$label),]

gambia_variables <- makeVlist(mw_psc)
gambia_variables[grep("weight", gambia_variables$label),]


head(mw_psc)
head(mw_wra)

#define deficiency levels
define_mn_deficiency <- function(
    df,
    population = c("young_child", "school_child", "woman"),
    rbp,
    fer,
    zn,
    crp,
    agp,
    time_blood_draw,
    fast
) {

  population <- match.arg(population)

  # --- helper: BRINDA regression adjustment ---
  # 2026-09-15 FIX. The previous version took exp(residual + intercept) for
  # every observation, i.e. it re-centred EVERYONE to log(CRP) = log(AGP) = 0
  # (CRP = 1 mg/L, AGP = 1 g/L). For the majority of children with CRP well
  # below 1 mg/L that RAISED ferritin (and lowered RBP), so iron deficiency
  # came out at 9.7% against the survey's own BRINDA flag of 20.1% (sf_c1;
  # report 21.7%). BRINDA (Namaste 2017; Larson 2017) instead uses the 10th
  # percentile of CRP and AGP as the reference and adjusts only observations
  # ABOVE it, leaving the uninflamed untouched; direction is fixed by `sign`
  # (+1: inflammation raises the marker, ferritin; -1: depresses it, RBP),
  # with a coefficient of the wrong sign clamped to zero.
  brinda_adjust <- function(y, crp, agp, sign = +1) {
    ok <- is.finite(y) & is.finite(crp) & is.finite(agp) &
      y > 0 & crp > 0 & agp > 0
    y_adj <- rep(NA_real_, length(y))
    if (sum(ok) > 20) {
      ly <- log(y[ok]); lc <- log(crp[ok]); la <- log(agp[ok])
      c_ref <- as.numeric(quantile(lc, 0.10)); a_ref <- as.numeric(quantile(la, 0.10))
      b  <- coef(lm(ly ~ lc + la))
      bC <- if (sign > 0) max(b[["lc"]], 0) else min(b[["lc"]], 0)
      bA <- if (sign > 0) max(b[["la"]], 0) else min(b[["la"]], 0)
      corr <- bC * pmax(lc - c_ref, 0) + bA * pmax(la - a_ref, 0)
      y_adj[ok] <- exp(ly - corr)
    }
    y_adj
  }

  # --- pull vectors ---
  rbp  <- df[[rbp]]
  fer  <- df[[fer]]
  zn   <- df[[zn]]
  crp  <- df[[crp]]
  agp  <- df[[agp]]
  time <- df[[time_blood_draw]]
  fast <- df[[fast]]

  # treat missing fast as non-fasting
  fast[is.na(fast)] <- 0

  # --- Vitamin A ---
  # Inflammation depresses RBP: sign = -1. (The pipeline re-derives the VAD
  # binary at run time with R/brinda_adjustment.R; this column is a fallback.)
  if (population == "young_child") {
    rbp_adj <- brinda_adjust(rbp, crp, agp, sign = -1)
    vad <- as.integer(rbp_adj < 0.70)
  } else {
    vad <- as.integer(rbp < 0.70)
  }

  # --- Iron ---
  # Prefer the survey's own BRINDA internal-regression ferritin (`sf_reg`,
  # report section 2.9) when the file carries it; its cut-off flag `sf_c1` is
  # exactly sf_reg < 12 (PSC) / < 15 (others). Fall back to the local BRINDA
  # (inflammation RAISES ferritin: sign = +1) only when sf_reg is absent.
  fer_adj <- if ("sf_reg" %in% names(df)) df[["sf_reg"]] else brinda_adjust(fer, crp, agp, sign = +1)

  # Report Table 2.3: ferritin < 12 ug/L for preschool children, < 15 ug/L for
  # school-aged children, women and men (the old code used 12 for SAC too).
  if (population == "young_child") {
    id <- as.integer(fer_adj < 12)
  } else {
    id <- as.integer(fer_adj < 15)
  }

  # --- Zinc ---
  # zn_gdl carries sentinel values (-100, 0) that are not measurements.
  zn[!is.na(zn) & zn <= 0] <- NA_real_
  zn_def <- rep(NA_integer_, length(zn))

  for (i in seq_along(zn)) {
    if (is.na(zn[i])) next

    # assume non-fasting
    if (population %in% c("young_child", "school_child")) {
      if (time[i] == 1) {
        zn_def[i] <- as.integer(zn[i] < 65)
      } else if (time[i] == 2) {
        zn_def[i] <- as.integer(zn[i] < 57)
      } else {
        zn_def[i] <- as.integer(zn[i] < 65)  # conservative fallback
      }
    }

    if (population == "woman") {
      if (time[i] == 1) {
        zn_def[i] <- as.integer(zn[i] < 66)
      } else if (time[i] == 2) {
        zn_def[i] <- as.integer(zn[i] < 59)
      } else {
        zn_def[i] <- as.integer(zn[i] < 66)  # conservative fallback
      }
    }
  }

  # --- return ---
  df_out <- df
  df_out$vitA_def <- vad
  df_out$iron_def <- id
  df_out$zinc_def <- zn_def

  return(df_out)
}
#Women
mw_wra <- define_mn_deficiency(
  df = mw_wra,
  population = "woman",
  rbp = "rbp",
  fer = "fer",
  zn  = "zn_gdl",
  crp = "crp",
  agp = "agp",
  time_blood_draw = "time_blood_draw",
  fast = "fast"
)

table(mw_wra$vitA_def)
table(mw_wra$iron_def)
table(mw_wra$zinc_def)

#Young children
mw_psc <- define_mn_deficiency(
  df = mw_psc,
  population = "young_child",
  rbp = "rbp",
  fer = "fer",
  zn  = "zn_gdl",
  crp = "crp",
  agp = "agp",
  time_blood_draw = "time_blood_draw",
  fast = "fast"
)
#School-age children
mw_sac <- define_mn_deficiency(
  df = mw_sac,
  population = "school_child",
  rbp = "rbp",
  fer = "fer",
  zn  = "zn_gdl",
  crp = "crp",
  agp = "agp",
  time_blood_draw = "time_blood_draw",
  fast = "fast"
)



#load the gps data
gps <- read.csv(here("data/DHS/dhs_Malawi_2015_cluster_gps.csv"))


df <- bind_rows(mw_psc %>% mutate(population="preschool children"),
                mw_sac %>% mutate(population="school-age children"),
                mw_wra %>% mutate(population="women"),
                mw_men %>% mutate(population="men")) %>%
  rename(cluster = mcluster) %>% left_join(gps, by="cluster") %>%
  rename(longitude = LONGNUM, latitude = LATNUM, Admin1=admin1.name, Admin2= admin2.name)
head(df)
#once the dataset is created, still merge in aggregated DHS as if it's from a proxy source

saveRDS(df, file=here("data/Malawi/clean_malawi_mn_data.RDS"))

table(df$population)
table(df$population, df$Admin1)
#outcomes
#fer  stfr   rbp   crp   agp zn_gdl incap_dr incap_retinol mrdr_ratio

#check for
