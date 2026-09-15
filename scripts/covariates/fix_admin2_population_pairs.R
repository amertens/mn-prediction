# =============================================================================
# scripts/covariates/fix_admin2_population_pairs.R   [JK-02, 2026-09-15]
#
# THE POPULATION TABLE WAS MISSING ONE TA OF EACH SAME-NAME PAIR IN MALAWI
#
# dashboard/data/admin2_population.rds is built by
# dashboard/data-raw/01_prepare_dashboard_data.R from the legacy external-
# predictor cache, which was deduplicated BY NAME (`!duplicated(ext$Admin2)`)
# - the defect class WS8g removed everywhere else. Malawi has four district
# names that occur in two regions (TA Lundu, TA Ngabu, TA Pemba, TA Malemia),
# so four spine units had no population: Chikwawa/TA Lundu, Nsanje/TA Ngabu,
# Salima/TA Pemba and Zomba/TA Malemia, two of them surveyed. Every
# burden-weighted protocol arm (scripts 12, 16, 21, 22, 27, 28, 32, 35)
# filters `is.finite(pop) & pop > 0`, so those two districts dropped out of
# the burden figures silently.
#
# Fix: fill the missing spine units from the survey-year WorldPop density the
# shared set already carries (`wpop_log_density_survey_year`, log1p persons
# per km2, GEE) times the polygon area, then apply the country's own
# 2023 projection factor and child / women shares as implied by its existing
# rows. The four units' populations are a few tens of thousands each, in line
# with their neighbours. Backup: admin2_population.rds.pre_pairfix.
#
#   Rscript -e "source('scripts/covariates/fix_admin2_population_pairs.R')"
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/admin2_keys.R")
P_PATH <- "dashboard/data/admin2_population.rds"
if (!file.exists(paste0(P_PATH, ".pre_pairfix"))) file.copy(P_PATH, paste0(P_PATH, ".pre_pairfix"))
P <- readRDS(P_PATH); P <- as.data.frame(unclass(P), stringsAsFactors = FALSE)
P$country_key <- admin2_country_label(P$country)
sp <- admin2_spine()
G <- read.csv("data/covariates/harmonized/gee_rwi_density_admin2.csv", check.names = FALSE, stringsAsFactors = FALSE)
BND <- readRDS("dashboard/data/admin2_boundaries.rds"); names(BND) <- c(gambia = "Gambia", ghana = "Ghana", sierraleone = "SierraLeone", malawi = "Malawi")[names(BND)]

key <- function(d, cn = d$country) paste(cn, d$Admin1, d$Admin2)
missing <- sp[!key(sp) %in% key(P, P$country_key), ]
cat(sprintf("population table: %d rows; spine units without population: %d\n", nrow(P), nrow(missing)))
if (nrow(missing)) {
  print(missing[, c("country", "Admin1", "Admin2", "engtype_2")], row.names = FALSE)
  add <- lapply(seq_len(nrow(missing)), function(i) {
    cn <- missing$country[i]; a1 <- missing$Admin1[i]; a2 <- missing$Admin2[i]
    b <- BND[[cn]]; poly <- b[as.character(b$Admin1) == a1 & as.character(b$Admin2) == a2, ]
    area_km2 <- as.numeric(sf::st_area(sf::st_transform(poly, 4326))) / 1e6
    g <- G[G$country == cn & G$Admin1 == a1 & G$Admin2 == a2, ]
    stopifnot(nrow(poly) == 1L, nrow(g) == 1L)
    ref <- P[P$country_key == cn & is.finite(P$population) & P$population > 0, ]
    # calibrate density x area to the table's own scale (the existing rows sit a
    # constant 1.16x above density x polygon area in Malawi; rank agreement 1.00)
    gg <- merge(merge(ref[, c("Admin1", "Admin2", "population")], G[G$country == cn, c("Admin1", "Admin2", "wpop_log_density_survey_year")], by = c("Admin1", "Admin2")),
                data.frame(Admin1 = as.character(b$Admin1), Admin2 = as.character(b$Admin2), area = as.numeric(sf::st_area(sf::st_transform(b, 4326))) / 1e6), by = c("Admin1", "Admin2"))
    calib <- median(gg$population / (expm1(gg$wpop_log_density_survey_year) * gg$area))
    pop <- expm1(g$wpop_log_density_survey_year) * area_km2 * calib
    proj <- median(ref$population_2023 / ref$population); sh_c <- median(ref$pop_child / ref$population); sh_w <- median(ref$pop_women / ref$population)
    data.frame(country = ref$country[1], Admin2 = a2, Admin1 = a1, population = pop, pop_year = as.integer(g$wpop_density_year),
               population_2023 = pop * proj, pop_child = pop * sh_c, pop_women = pop * sh_w, pop_child_2023 = pop * proj * sh_c, pop_women_2023 = pop * proj * sh_w,
               stringsAsFactors = FALSE)
  })
  ADD <- bind_rows(add)
  cat("\nfilled from WorldPop density x polygon area:\n"); print(ADD[, c("Admin1", "Admin2", "population", "population_2023", "pop_child", "pop_women")], row.names = FALSE)
  P2 <- bind_rows(P[, setdiff(names(P), "country_key")], ADD[, setdiff(names(P), "country_key")])
  stopifnot(!anyDuplicated(paste(admin2_country_label(P2$country), P2$Admin1, P2$Admin2)))
  saveRDS(P2, P_PATH)
  cat(sprintf("\nwritten: %d rows (was %d); every spine unit now has a population row: %s\n", nrow(P2), nrow(P), all(key(sp) %in% key(P2, admin2_country_label(P2$country)))))
} else cat("nothing to fix\n")
cat("DONE\n")
