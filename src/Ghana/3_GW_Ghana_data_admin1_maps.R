# GMS 2017 – Admin-1 aggregation + maps (design-based)
rm(list = ls())

# ---- Packages ----
library(dplyr)
library(tidyverse)
library(haven)
library(here)
library(purrr)
library(labelled)
library(sf)
library(terra)
library(readxl)
library(survey)
library(srvyr)        # tidy interface to 'survey'
library(stringr)
library(viridis)
library(scales)       # percent_format(), number_format()
library(ggplot2)
library(geodata)      # for GADM boundaries

# ---- Data ----
d <- readRDS(here("data", "IPD", "Ghana", "Ghana_GMS_cleaned.rds"))

# (Optional quick check)
# summary(d$gw_cRBPAdjBrinda); table(is.na(d$gw_cRBPAdjBrinda))

# ---- Merge GEE data (same as your script) ----
gee <- read.csv(here("data/gee/GMS_gee_merge_GH_2017.csv")) %>%
  select(latnum, longnum, trmm_mean_10km:temp_mean_50km)
colnames(gee) <- paste0("gee_", colnames(gee))
d <- left_join(d, gee, by = c("longitude" = "gee_longnum", "latitude" = "gee_latnum"))

# ---- Admin membership (get Admin-1 and Admin-2) ----
d$lat <- as.numeric(d$latitude)
d$lon <- as.numeric(d$longitude)

# polygons for membership assignment (level 2, contains level 1 names)
poly.adm2 <- geodata::gadm(country = "GH", level = 2, path = tempdir()) |>
  sf::st_as_sf() |>
  dplyr::select(NAME_1, NAME_2) |>
  dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)

# polygons for mapping at Admin-1
poly.adm1 <- geodata::gadm(country = "GH", level = 1, path = tempdir()) |>
  sf::st_as_sf() |>
  dplyr::select(NAME_1) |>
  dplyr::rename(Admin1 = NAME_1)

d_sf      <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
poly.adm2 <- st_transform(poly.adm2, crs = 4326)
poly.adm1 <- st_transform(poly.adm1, crs = 4326)

# attach Admin1/Admin2 to points
df <- st_join(d_sf, poly.adm2, join = st_within)

# (Old borders you used earlier, not needed for Admin-1 maps—keep if you still need them for merges)
# poly.adm1_old <- st_read(here("data/old_ghana_admin_boundaries/gadm36_GHA_1.shp")) %>%
#   select(NAME_1) %>% rename(Admin1_old = NAME_1)
# poly.adm2_old <- st_read(here("data/old_ghana_admin_boundaries/gadm36_GHA_2.shp")) %>%
#   select(NAME_2) %>% rename(Admin2_old = NAME_2)
# df <- st_join(df, poly.adm1_old, join = st_within)
# df <- st_join(df, poly.adm2_old, join = st_within)

# ---- Outcomes to numeric 0/1 ----
df <- df %>%
  mutate(
    gw_cIDAdjBrinda = as.numeric(gw_cIDAdjBrinda),
    gw_wIDAdjBrinda = as.numeric(gw_wIDAdjBrinda),
    gw_cVAD         = as.numeric(gw_cVAD),
    gw_wFolateDef   = as.numeric(gw_wFolateDef),
    gw_wB12Def      = as.numeric(gw_wB12Def),
    gw_wVADAdjThurn = as.numeric(gw_wVADAdjThurn)
  )

# ---- Survey design setup ----
# Drop geometry
d_nog <- sf::st_drop_geometry(df)

# Set your survey design identifiers (EDIT if your names differ)
psu_var    <- "gw_EACode"   # cluster / EA / PSU column
strata_var <- "gw_Strata"   # strata column
wt_var     <- "gw_sWeight"  # individual analysis weight

# Coerce types
d_nog[[psu_var]]    <- as.factor(d_nog[[psu_var]])
d_nog[[strata_var]] <- as.factor(d_nog[[strata_var]])
d_nog[[wt_var]]     <- as.numeric(d_nog[[wt_var]])

# Handle lonely PSUs for SE computation
options(survey.lonely.psu = "adjust")

# ---- Helper: Admin-1 prevalence (design-based) ----
compute_admin1_prev <- function(data, outcomes) {
  dat <- data[rowSums(sapply(outcomes, function(v) !is.na(data[[v]]))) > 0, , drop = FALSE]

  des <- srvyr::as_survey_design(
    dat,
    ids     = !!rlang::sym(psu_var),
    strata  = !!rlang::sym(strata_var),
    weights = !!rlang::sym(wt_var),
    nest    = TRUE
  )

  purrr::map_dfr(outcomes, function(v) {
    des %>%
      group_by(Admin1) %>%
      summarise(
        prev = survey_mean(!!rlang::sym(v), vartype = c("se", "ci"), na.rm = TRUE),
        # helpful diagnostics
        n_unw = unweighted(sum(!is.na(!!rlang::sym(v)))),
        n_psu = dplyr::n_distinct(.data[[psu_var]][!is.na(.data[[v]])])
      ) %>%
      mutate(outcome = v)
  })
}

# ---- Compute Admin-1 estimates for all outcomes ----
outcomes <- c("gw_cIDAdjBrinda", "gw_wIDAdjBrinda", "gw_cVAD",
              "gw_wFolateDef",  "gw_wB12Def",      "gw_wVADAdjThurn")

adm1_prev <- compute_admin1_prev(d_nog, outcomes) %>%
  mutate(outcome = recode(outcome,
                          gw_cIDAdjBrinda = "Child iron deficiency (adj)",
                          gw_wIDAdjBrinda = "Women iron deficiency (adj)",
                          gw_cVAD         = "Child VAD",
                          gw_wFolateDef   = "Women folate deficiency",
                          gw_wB12Def      = "Women B12 deficiency",
                          gw_wVADAdjThurn = "Women VAD (adj)"
  ))

# ---- Join to Admin-1 polygons + CV ----
poly_adm1_joined <- poly.adm1 %>%
  left_join(adm1_prev, by = "Admin1") %>%
  mutate(cv = if_else(!is.na(prev) & prev > 0, prev_se / prev, NA_real_))

# Quick diagnostic: which Admin1 have NA SE?
na_se <- poly_adm1_joined %>%
  st_drop_geometry() %>%
  filter(is.na(prev_se)) %>%
  distinct(Admin1, outcome)
if (nrow(na_se) > 0) {
  message("Admin-1 with NA SE (likely sparse / one-PSU strata):")
  print(na_se)
}

# ---- Maps: prevalence + uncertainty (CV) ----
# 1) Prevalence map (0–100%)
p_prev <- ggplot(poly_adm1_joined) +
  geom_sf(aes(fill = prev), color = "white", size = 0.15) +
  facet_wrap(~ outcome, ncol = 2) +
  scale_fill_viridis_c(
    name   = "Prevalence",
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    na.value = "grey90"
  ) +
  labs(
    title    = "GMS 2017: Admin-1 deficiency prevalences (design-based)",
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

# 2) Uncertainty map: Coefficient of Variation (SE / prevalence)
p_cv <- ggplot(poly_adm1_joined) +
  geom_sf(aes(fill = cv), color = "white", size = 0.15) +
  facet_wrap(~ outcome, ncol = 2) +
  scale_fill_viridis_c(
    name   = "CV (SE / prev)",
    labels = scales::number_format(accuracy = 0.01),
    na.value = "grey90"
  ) +
  labs(
    title    = "Uncertainty in Admin-1 estimates",
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

# ---- Save outputs (optional) ----
# dir.create(here("figures"), showWarnings = FALSE)
# ggsave(here("figures", "GMS2017_admin1_prevalence_faceted.png"),  p_prev, width = 10, height = 8, dpi = 300)
# ggsave(here("figures", "GMS2017_admin1_uncertainty_faceted.png"), p_cv,   width = 10, height = 8, dpi = 300)

# ---- Helper to map a single outcome at Admin-1 (optional) ----
map_outcome_admin1 <- function(outcome_label, value = c("prev", "cv")) {
  value <- match.arg(value)
  dat <- poly_adm1_joined %>% filter(outcome == outcome_label)
  ggplot(dat) +
    geom_sf(aes(fill = .data[[value]]), color = "white", size = 0.15) +
    scale_fill_viridis_c(
      name   = if (value == "prev") "Prevalence" else "CV",
      labels = if (value == "prev") scales::percent_format(accuracy = 1) else scales::number_format(accuracy = 0.01),
      limits = if (value == "prev") c(0,1) else NULL,
      na.value = "grey90"
    ) +
    labs(
      title    = paste0("GMS 2017: ", outcome_label),
      subtitle = if (value == "prev") "Design-based prevalence by region" else "Uncertainty (SE / prevalence) by region",
      caption  = "Ghana Micronutrient Survey 2017"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      legend.position  = "right"
    )
}

# Example:
 print(map_outcome_admin1("Women iron deficiency (adj)", value = "prev"))
 print(map_outcome_admin1("Women iron deficiency (adj)", value = "cv"))
