# =============================================================================
# scripts/covariates/build_admin2_spine.R   [JK-01, 2026-09-15]
#
# THE CANONICAL ADMIN-2 SPINE WITH STABLE CODES
#
# Every table in the project is keyed by (country, Admin1, Admin2) names taken
# from GADM 4.1 (the polygons in dashboard/data/admin2_boundaries.rds are GADM
# 4.1 level 2 minus Malawi's 13 water bodies). Names are a fragile key: a
# name-only join fanned rows in Malawi (six district names occur in two
# regions), fuzzy matching has mis-linked a district once already, and no
# provider uses the same spelling. This table gives every unit its GADM code
# (GID_2, stable within GADM 4.1), so joins can be made and checked on codes,
# and carries a `pcode` column for the OCHA COD-AB P-code once the crosswalk
# is built (scripts/covariates/build_pcode_crosswalk.R; empty until then).
#
# Checks: one row per polygon; (country, Admin1, Admin2) unique; every unit of
# the shared predictor set is in the spine and vice versa.
#
#   Rscript -e "source('scripts/covariates/build_admin2_spine.R')"
# -> metadata/admin2_spine.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
names(BND) <- c(gambia = "Gambia", ghana = "Ghana", sierraleone = "SierraLeone", malawi = "Malawi")[names(BND)]

rows <- lapply(names(ISO), function(cn) {
  b <- sf::st_drop_geometry(BND[[cn]]); b$Admin1 <- as.character(b$Admin1); b$Admin2 <- as.character(b$Admin2)
  g <- sf::st_drop_geometry(readRDS(sprintf("data/admin_boundaries/gadm41_%s_2.rds", ISO[[cn]])))
  j <- match(paste(b$Admin1, b$Admin2), paste(g$NAME_1, g$NAME_2))
  if (anyNA(j)) stop(cn, ": polygons without a GADM 4.1 match: ", paste(paste(b$Admin1, b$Admin2)[is.na(j)], collapse = "; "))
  data.frame(country = cn, iso3 = ISO[[cn]], Admin1 = b$Admin1, Admin2 = b$Admin2,
             gid_1 = g$GID_1[j], gid_2 = g$GID_2[j], hasc_2 = g$HASC_2[j], engtype_2 = g$ENGTYPE_2[j],
             varname_2 = g$VARNAME_2[j], pcode = NA_character_, stringsAsFactors = FALSE)
})
SP <- bind_rows(rows)
stopifnot(!anyDuplicated(paste(SP$country, SP$Admin1, SP$Admin2)), !anyDuplicated(SP$gid_2))

# keep the P-codes of an earlier build, if any
old <- "metadata/admin2_spine.csv"
if (file.exists(old)) { o <- read.csv(old, stringsAsFactors = FALSE); if ("pcode" %in% names(o)) SP$pcode <- o$pcode[match(SP$gid_2, o$gid_2)] }

S <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE, stringsAsFactors = FALSE)[, c("country", "Admin1", "Admin2")]
k <- function(d) paste(d$country, d$Admin1, d$Admin2)
stopifnot(setequal(k(S), k(SP)), nrow(S) == nrow(SP))
write.csv(SP, old, row.names = FALSE)
cat(sprintf("admin2 spine: %d units (%s); P-codes filled: %d\n", nrow(SP), paste(sprintf("%s %d", names(table(SP$country)), table(SP$country)), collapse = ", "), sum(!is.na(SP$pcode))))
print(table(SP$country, SP$engtype_2))
cat("DONE\n")
