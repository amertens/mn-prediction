

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
  brinda_adjust <- function(y, crp, agp) {
    ok <- is.finite(y) & is.finite(crp) & is.finite(agp) &
      y > 0 & crp > 0 & agp > 0

    y_adj <- rep(NA_real_, length(y))

    if (sum(ok) > 20) {
      fit <- lm(log(y[ok]) ~ log(crp[ok]) + log(agp[ok]))
      y_adj[ok] <- exp(resid(fit) + coef(fit)[1])
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
  if (population == "young_child") {
    rbp_adj <- brinda_adjust(rbp, crp, agp)
    vad <- as.integer(rbp_adj < 0.70)
  } else {
    vad <- as.integer(rbp < 0.70)
  }

  # --- Iron ---
  fer_adj <- brinda_adjust(fer, crp, agp)

  if (population %in% c("young_child", "school_child")) {
    id <- as.integer(fer_adj < 12)
  } else {
    id <- as.integer(fer_adj < 15)
  }

  # --- Zinc ---
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
