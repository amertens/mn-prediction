
rm(list=ls())
library(dplyr)
library(tidyverse)
library(haven)
library(here)
library(purrr)
library(labelled)
library(sl3)
library(origami)
library(tlverse)
library(caret)
library(data.table)
library(ck37r)
library(rdhs)
library(maps)
library(sf)
library(terra)
library(caret)
library(readxl)
library(survey)
library(srvyr)     # tidy interface to 'survey'
library(stringr)
library(purrr)
library(viridis)
library(scales)


d <- readRDS(here("data", "IPD", "Ghana", "Ghana_GMS_cleaned.rds"))


summary(d$gw_cRBPAdjBrinda)
table(is.na(d$gw_cRBPAdjBrinda))





#-------------------------------------------------------------------------------
# GEE data
#-------------------------------------------------------------------------------


gee <- read.csv(here("data/gee/GMS_gee_merge_GH_2017.csv")) %>% select(latnum, longnum, trmm_mean_10km:temp_mean_50km)
colnames(gee)= paste0("gee_",colnames(gee))
head(gee)

d <- left_join(d, gee, by = c("longitude"="gee_longnum", "latitude"="gee_latnum"))


#merged_proxies <- readRDS(here::here("data/DHS/clean/dhs_Ghana_2019_gee_fp_map_merge.RDS"))

#-------------------------------------------------------------------------------
#  Admin 1 and 2 membership
#-------------------------------------------------------------------------------

summary(d$latitude)
summary(d$longitude)
d$lat= as.numeric(d$latitude)
d$lon= as.numeric(d$longitude)

poly.adm <- geodata::gadm(country="GH", level=2, path=tempdir())
poly.adm <- sf::st_as_sf(poly.adm) %>% select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
d_sf <- st_as_sf(d, coords = c("longitude","latitude"), crs = 4326)
poly.adm <- st_transform(poly.adm, crs = 4326)
#df <- as.data.frame(st_join(d_sf, poly.adm, join = st_within))
df <- (st_join(d_sf, poly.adm, join = st_within))

#get old borders for some merges
poly.adm1_old <- st_read(here("data/old_ghana_admin_boundaries/gadm36_GHA_1.shp")) %>%
  select(NAME_1) %>%
  rename(Admin1_old = NAME_1)
poly.adm2_old <- st_read(here("data/old_ghana_admin_boundaries/gadm36_GHA_2.shp")) %>%
  select(NAME_2) %>%
  rename(Admin2_old = NAME_2)
df <- (st_join(df, poly.adm1_old, join = st_within))
df <- (st_join(df, poly.adm2_old, join = st_within))

# --- ensure the six outcomes are coded 0/1 numeric (already in your script) ---
df <- df %>%
  mutate(
    gw_cIDAdjBrinda = as.numeric(gw_cIDAdjBrinda),
    gw_wIDAdjBrinda = as.numeric(gw_wIDAdjBrinda),
    gw_cVAD         = as.numeric(gw_cVAD),
    gw_wFolateDef   = as.numeric(gw_wFolateDef),
    gw_wB12Def      = as.numeric(gw_wB12Def),
    gw_wVADAdjThurn = as.numeric(gw_wVADAdjThurn)
  )


# Drop geometry
d_nog <- sf::st_drop_geometry(df)

# --- Set PSU and strata identifiers manually ---
# Replace with the correct names in your dataset
psu_var    <- "gw_EACode"   # e.g., cluster / EA / PSU column
strata_var <- "gw_Strata"      # e.g., survey strata column

# Make sure they're factors/numeric as needed
d_nog[[psu_var]]    <- as.factor(d_nog[[psu_var]])
d_nog[[strata_var]] <- as.factor(d_nog[[strata_var]])
d_nog$gw_sWeight       <- as.numeric(d_nog$gw_sWeight)

#set options for survey design
options(survey.lonely.psu = "adjust")

# --- Helper: compute admin-2 prevalence for outcomes ---
compute_adm2_prev <- function(data, outcomes) {
  dat <- data[rowSums(sapply(outcomes, function(v) !is.na(data[[v]]))) > 0, , drop = FALSE]

  des <- srvyr::as_survey_design(
    dat,
    ids     = !!rlang::sym(psu_var),
    strata  = !!rlang::sym(strata_var),
    weights = !!rlang::sym("gw_sWeight"),
    nest    = TRUE
  )

  purrr::map_dfr(outcomes, function(v) {
    des %>%
      dplyr::group_by(Admin2) %>%
      dplyr::summarise(prev = survey_mean(!!rlang::sym(v), vartype = c("se","ci"), na.rm = TRUE)) %>%
      dplyr::mutate(outcome = v)
  })
}

# --- All outcomes (children + women use same sWeight now) ---
outcomes <- c("gw_cIDAdjBrinda","gw_wIDAdjBrinda","gw_cVAD",
              "gw_wFolateDef","gw_wB12Def","gw_wVADAdjThurn")

adm2_prev <- compute_adm2_prev(d_nog, outcomes)

# --- Relabel outcomes for clarity ---
adm2_prev <- adm2_prev %>%
  dplyr::mutate(outcome = dplyr::recode(outcome,
                                        gw_cIDAdjBrinda = "Child iron deficiency (adj)",
                                        gw_wIDAdjBrinda = "Women iron deficiency (adj)",
                                        gw_cVAD         = "Child VAD",
                                        gw_wFolateDef   = "Women folate deficiency",
                                        gw_wB12Def      = "Women B12 deficiency",
                                        gw_wVADAdjThurn = "Women VAD (adj)"
  ))


#save outcomes for merge with predictors
saveRDS(adm2_prev, here("data", "IPD", "Ghana", "Ghana_GMS_admin2_prevalence.rds"))

head(adm2_prev)

#------------------------------------------------------------------------------
# Optional mapping
#------------------------------------------------------------------------------

# Join survey estimates to polygons
poly_adm_joined <- poly.adm %>%
  left_join(adm2_prev, by = "Admin2") %>%
  mutate(
    cv = if_else(!is.na(prev) & prev > 0, prev_se / prev, NA_real_)
  )

# Quick diagnostic: which Admin2 have NA SE?
na_se <- poly_adm_joined %>%
  st_drop_geometry() %>%
  filter(is.na(prev_se)) %>%
  distinct(Admin2, outcome)

print(na_se)

# 1) Faceted prevalence map (0–100%)
p_prev <- ggplot(poly_adm_joined) +
  geom_sf(aes(fill = prev), color = "white", size = 0.15) +
  facet_wrap(~ outcome, ncol = 2) +
  scale_fill_viridis_c(
    name   = "Prevalence",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    na.value = "grey90"
  ) +
  labs(
    title    = "GMS 2017: Admin-2 deficiency prevalences (design-based)",
    subtitle = "Weighted with sWeight; CIs account for clustering & stratification",
    caption  = "Ghana Micronutrient Survey 2017"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",
    strip.text       = element_text(face = "bold")
  )

print(p_prev)

# 2) Faceted uncertainty map: Coefficient of Variation (SE / prevalence)
poly_adm_joined <- poly_adm_joined %>%
  mutate(cv = if_else(!is.na(prev) & prev > 0, prev_se / prev, NA_real_))

poly_adm_joined$prev
poly_adm_joined$se

poly_adm_joined$cv

p_cv <- ggplot(poly_adm_joined) +
  geom_sf(aes(fill = cv), color = "white", size = 0.15) +
  facet_wrap(~ outcome, ncol = 2) +
  scale_fill_viridis_c(
    name   = "CV (SE / prev)",
    labels = number_format(accuracy = 0.01),
    na.value = "grey90"
  ) +
  labs(
    title    = "Uncertainty in admin-2 estimates",
    subtitle = "Higher CV = lower precision; grey = no data",
    caption  = "Ghana Micronutrient Survey 2017"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",
    strip.text       = element_text(face = "bold")
  )

print(p_cv)

# # Optional: save to files
# ggsave(filename = here::here("figures", "GMS2017_admin2_prevalence_faceted.png"),
#        plot = p_prev, width = 10, height = 8, dpi = 300)
# ggsave(filename = here::here("figures", "GMS2017_admin2_uncertainty_faceted.png"),
#        plot = p_cv, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------
# (Optional) Single-outcome helper if you want one map at a time
# ------------------------------------------------------------
map_outcome <- function(outcome_label, value = c("prev","cv")) {
  value <- match.arg(value)
  dat <- poly_adm_joined %>% filter(outcome == outcome_label)
  ggplot(dat) +
    geom_sf(aes(fill = .data[[value]]), color = "white", size = 0.15) +
    scale_fill_viridis_c(
      name = if (value == "prev") "Prevalence" else "CV",
      #labels = if (value == "prev") percent_format(accuracy = 1) else number_format(accuracy = 0.01),
      limits = if (value == "prev") c(0,1) else NULL,
      na.value = "grey90"
    ) +
    labs(
      title = paste0("GMS 2017: ", outcome_label),
      subtitle = if (value == "prev") "Design-based prevalence by district"
      else "Uncertainty (SE / prevalence) by district",
      caption = "Ghana Micronutrient Survey 2017"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      legend.position  = "right"
    )
}
