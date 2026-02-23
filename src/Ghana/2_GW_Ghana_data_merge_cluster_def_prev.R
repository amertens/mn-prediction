
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


df <- readRDS(here("data", "IPD", "Ghana", "Ghana_GMS_cleaned.rds"))

# Drop geometry
d_nog <- sf::st_drop_geometry(df)
# --- Set PSU and strata identifiers manually ---
psu_var <- "gw_EACode"  # e.g., cluster / EA / PSU column
strata_var <- "gw_Strata" # e.g., survey strata column
# Make sure they're factors/numeric as needed
d_nog[[psu_var]] <- as.factor(d_nog[[psu_var]])
d_nog[[strata_var]] <- as.factor(d_nog[[strata_var]])
d_nog$gw_sWeight <- as.numeric(d_nog$gw_sWeight)


# 0) Ensure cnum exists; if not, default to PSU id
if (!("gw_cnum" %in% names(d_nog))) {
  message("`cnum` not found; using PSU variable as cluster id.")
  d_nog$gw_cnum <- d_nog[[psu_var]]
}
d_nog$gw_cnum <- as.factor(d_nog$gw_cnum)

# 1) Outcomes
outcomes <- c("gw_cIDAdjBrinda","gw_wIDAdjBrinda","gw_cVAD",
              "gw_wFolateDef","gw_wB12Def","gw_wVADAdjThurn")

# 2) Survey design (same as you used before)
options(survey.lonely.psu = "adjust")   # handle 1-PSU strata gracefully

des_all <- srvyr::as_survey_design(
  d_nog,
  ids     = !!rlang::sym(psu_var),
  strata  = !!rlang::sym(strata_var),
  weights = !!rlang::sym("gw_sWeight"),
  nest    = TRUE
)

# 3) Cluster-level prevalences (domain = cnum)
compute_cluster_prev <- function(des, outcomes) {
  purrr::map_dfr(outcomes, function(v) {
    des %>%
      group_by(gw_cnum) %>%
      summarise(
        prev = survey_mean(!!rlang::sym(v), vartype = c("se","ci"), na.rm = TRUE),
        # diagnostics (unweighted count with non-missing outcome, and # of PSUs contributing)
        n_unw = unweighted(sum(!is.na(!!rlang::sym(v)))),
        n_psu = dplyr::n_distinct(.data[[psu_var]][!is.na(.data[[v]])])
      ) %>%
      mutate(outcome = v)
  }) %>%
    mutate(outcome = dplyr::recode(outcome,
                                   gw_cIDAdjBrinda = "Child iron deficiency (adj)",
                                   gw_wIDAdjBrinda = "Women iron deficiency (adj)",
                                   gw_cVAD         = "Child VAD",
                                   gw_wFolateDef   = "Women folate deficiency",
                                   gw_wB12Def      = "Women B12 deficiency",
                                   gw_wVADAdjThurn = "Women VAD (adj)"
    ))
}

cluster_prev <- compute_cluster_prev(des_all, outcomes)

# 4) (Optional) attach cluster point geometry (one row per cnum)
#    Here we take the mean lon/lat within cluster for plotting later.
cluster_coords <- d_nog %>%
  group_by(gw_cnum) %>%
  summarise(
    lon = mean(as.numeric(longitude), na.rm = TRUE),
    lat = mean(as.numeric(latitude), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(lon), is.finite(lat))

cluster_prev <- cluster_prev %>%
  left_join(cluster_coords, by = "gw_cnum")

cluster_prev_sf <- st_as_sf(cluster_prev, coords = c("lon","lat"), crs = 4326)

# 5) Save results
saveRDS(cluster_prev,   here("data","IPD","Ghana","GMS_cluster_prevalence.rds"))
saveRDS(cluster_prev_sf,here("data","IPD","Ghana","GMS_cluster_prevalence_points.rds"))

# 6) Quick peek
cluster_prev %>% arrange(outcome, gw_cnum) %>% print(n = 20)

#7) (Optional) simple scatter of cluster prevalences for one outcome
ggplot(cluster_prev %>% filter(outcome == "Women iron deficiency (adj)"),
       aes(x = reorder(gw_cnum, prev), y = prev)) +
  geom_point() + coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster (cnum)", y = "Prevalence",
       title = "Cluster-level prevalence: Women iron deficiency (adj)") +
  theme_minimal()



df <- cluster_prev_sf %>% select(gw_cnum, outcome, prev, geometry) %>% mutate(sex=1, age=55)

head(df)
