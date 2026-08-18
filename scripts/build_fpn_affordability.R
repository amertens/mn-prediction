#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_fpn_affordability.R
#
# Cost & affordability of a healthy diet (World Bank "Food Prices for Nutrition")
# as a NATIONAL predictor, broadcast to every GADM Admin-2 ->
# data/FPN/<Country>_fpn_admin2.csv with `fpn_` columns.
#
# WHY: the economic-access dimension of micronutrient-dense foods, distinct from
# the WFP raw commodity prices already in the pipeline. National-only -> helps
# the cross-country pooled / transfer model. See docs/ra_new_data_sources_specs.md §7.
#
# -----------------------------------------------------------------------------
# ACCESS NOTE (verified 2026-07): FPN has MIGRATED off the classic World Bank
# indicator API. The old `source=88` data endpoint returns "indicator not found"
# even for valid codes, and the Data360 replacement API returned generic errors
# in testing. The `CoNA_*` (nutrient-adequate) codes are RETIRED; FPN 4.0 uses
# `CoHD_*` (healthy diet) and `CoCA_*` (calorie-adequate). So this builder uses a
# one-time MANUAL CSV export (no key), which is reliable today:
#
#   1. Go to  https://databank.worldbank.org/source/food-prices-for-nutrition
#   2. Select our countries, the CoHD/CoCA series below, and all years.
#   3. Download as CSV (default wide layout) into  data/FPN/raw/  (any filename).
#
# The script auto-detects the DataBank wide layout (Country Code, Series Code,
# "YYYY [YRYYYY]" year columns), picks each country's survey-year value (nearest
# available), and broadcasts it to every Admin-2. Prints what it matched.
#
# USAGE:  Rscript scripts/build_fpn_affordability.R        # all countries
# =============================================================================
suppressMessages({library(here); library(dplyr)})

ISO3 <- c(Gambia="GMB",Ghana="GHA",SierraLeone="SLE",Malawi="MWI",Tanzania="TZA")
SURVEY_YEAR <- c(Gambia=2018, Ghana=2017, SierraLeone=2013, Malawi=2016, Tanzania=2010)
# FPN 4.0 series code -> output column (CoHD = healthy diet; CoCA = calorie-adequate)
SERIES <- c(fpn_cohd_ppp        = "CoHD_PPP",
            fpn_cohd_headcount  = "CoHD_headcount",
            fpn_cohd_fexp       = "CoHD_fexp",
            fpn_coca_headcount  = "CoCA_headcount",
            fpn_cohd_lns_prop   = "CoHD_lns_prop",   # legumes/nuts/seeds cost share
            fpn_cohd_asf_prop   = "CoHD_asf_prop")   # animal-source-food cost share
RAW_DIR <- here::here("data","FPN","raw"); dir.create(RAW_DIR, showWarnings=FALSE, recursive=TRUE)
OUT_DIR <- here::here("data","FPN");        dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

args <- commandArgs(trailingOnly=TRUE)
countries <- if (length(args)) args else names(ISO3)
bad <- setdiff(countries, names(ISO3)); if (length(bad)) stop("Unknown country: ", paste(bad, collapse=", "))

raw_files <- list.files(RAW_DIR, pattern="\\.csv$", full.names=TRUE, ignore.case=TRUE)
if (!length(raw_files))
  stop("No FPN CSV in data/FPN/raw/. Export from ",
       "https://databank.worldbank.org/source/food-prices-for-nutrition (see script header).")

# read + stack all raw CSVs; detect DataBank wide layout
read_databank <- function(f) {
  d <- utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  nm <- names(d)
  cc <- nm[grep("Country Code|Country.Code|iso3|REF_AREA", nm, ignore.case=TRUE)][1]
  sc <- nm[grep("Series Code|Series.Code|INDICATOR|Indicator Code", nm, ignore.case=TRUE)][1]
  yr <- nm[grepl("^\\s*\\d{4}", nm) | grepl("YR\\d{4}", nm)]
  if (is.na(cc) || is.na(sc) || !length(yr)) return(NULL)
  long <- do.call(rbind, lapply(yr, function(y){
    year <- as.integer(sub(".*?(\\d{4}).*", "\\1", y))
    data.frame(iso3=toupper(trimws(d[[cc]])), series=trimws(d[[sc]]),
               year=year, value=suppressWarnings(as.numeric(d[[y]])),
               stringsAsFactors=FALSE)
  }))
  long[is.finite(long$value), ]
}
fpn <- do.call(rbind, Filter(Negate(is.null), lapply(raw_files, read_databank)))
if (is.null(fpn) || !nrow(fpn))
  stop("Could not parse any DataBank CSV in data/FPN/raw/ (expected Country Code + Series Code + year cols).")
cat("## FPN affordability — parsed", nrow(fpn), "obs;",
    length(unique(fpn$series)), "series,", length(unique(fpn$iso3)), "countries\n")

source(here::here("R","data_prep.R"))  # load_gadm_cached()

for (cn in countries) {
  cat("\n#####", cn, "#####\n")
  yr <- SURVEY_YEAR[[cn]]; vals <- list()
  for (out in names(SERIES)) {
    d <- fpn[fpn$iso3 == ISO3[[cn]] & fpn$series == SERIES[[out]], ]
    if (!nrow(d)) next
    r <- d[which.min(abs(d$year - yr)), ]
    vals[[out]] <- r$value
    vals[["fpn_year"]] <- if (is.null(vals[["fpn_year"]])) r$year else vals[["fpn_year"]]
  }
  if (!length(vals)) { cat("  no matching series (check the exported Series Codes)\n"); next }
  cat(sprintf("  CoHD headcount=%s ; year used=%s (survey %d)\n",
              if(!is.null(vals$fpn_cohd_headcount)) round(vals$fpn_cohd_headcount,1) else "NA",
              if(!is.null(vals$fpn_year)) vals$fpn_year else "NA", yr))
  gadm <- tryCatch(load_gadm_cached(ISO3[[cn]], level=2), error=function(e) NULL)
  if (is.null(gadm)) { cat("  GADM load failed; skipping\n"); next }
  g <- unique(sf::st_drop_geometry(gadm)[, intersect(c("NAME_1","NAME_2"), names(gadm)), drop=FALSE])
  names(g)[names(g)=="NAME_1"] <- "Admin1"; names(g)[names(g)=="NAME_2"] <- "Admin2"
  g$Admin1 <- trimws(g$Admin1); g$Admin2 <- trimws(g$Admin2)
  for (k in names(vals)) g[[k]] <- vals[[k]]
  fn <- file.path(OUT_DIR, sprintf("%s_fpn_admin2.csv", cn))
  write.csv(g, fn, row.names=FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows, %d fpn_ cols)\n",
              basename(fn), nrow(g), sum(grepl("^fpn_", names(g)))))
}
cat("\nDONE. Domain FPN=list(prefix=\"fpn_\") + guarded merge join are wired in config/merges.\n")
