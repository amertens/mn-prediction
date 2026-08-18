#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_faostat_supply.R
#
# National food & nutrient SUPPLY from FAOSTAT Food Balance Sheets (FBS),
# broadcast to every GADM Admin-2 -> data/FAOSTAT/<Country>_fao_admin2.csv
# with `fao_` columns (constant within country).
#
# WHY: national per-capita supply of energy/protein and food-group quantities is
# the upstream availability signal for dietary micronutrients — the backbone FAO
# and GBD use for the zinc inadequate-intake model. National-only -> helps the
# cross-country pooled / transfer model. See docs/ra_new_data_sources_specs.md §6.
#
# ACCESS (public, no key, no login): FAOSTAT public BULK download of Food Balance
# Sheets, normalized CSV (the installed FAOSTAT CRAN package forces an API login
# it should not need, so we fetch the bulk zip directly with download.file()).
# ~55 MB, cached under data/FAOSTAT/raw/. Covers our survey years 2010-2018.
# Written defensively — prints the element/item labels it matched.
#
# USAGE:  Rscript scripts/build_faostat_supply.R           # all countries
#         Rscript scripts/build_faostat_supply.R Ghana
# =============================================================================
suppressMessages({library(here); library(dplyr); library(data.table)})
FBS_URL <- "https://bulks-faostat.fao.org/production/FoodBalanceSheets_E_All_Data_(Normalized).zip"

ISO3 <- c(Gambia="GMB",Ghana="GHA",SierraLeone="SLE",Malawi="MWI",Tanzania="TZA")
FAO_AREA <- c(Gambia="Gambia", Ghana="Ghana", SierraLeone="Sierra Leone",
              Malawi="Malawi", Tanzania="United Republic of Tanzania")
SURVEY_YEAR <- c(Gambia=2018, Ghana=2017, SierraLeone=2013, Malawi=2016, Tanzania=2010)
RAW_DIR <- here::here("data","FAOSTAT","raw"); dir.create(RAW_DIR, showWarnings=FALSE, recursive=TRUE)
OUT_DIR <- here::here("data","FAOSTAT");       dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

args <- commandArgs(trailingOnly=TRUE)
countries <- if (length(args)) args else names(ISO3)
bad <- setdiff(countries, names(ISO3)); if (length(bad)) stop("Unknown country: ", paste(bad, collapse=", "))

# element label -> output column (nutrient supply, whole-diet)
ELEM_NUTRIENT <- c(
  "fao_kcal_cap_day" = "Food supply \\(kcal/capita/day\\)",
  "fao_protein_g"    = "Protein supply quantity \\(g/capita/day\\)",
  "fao_fat_g"        = "Fat supply quantity \\(g/capita/day\\)")
# item groups for food-supply quantity (kg/capita/yr)
ITEM_GROUP <- c(
  "fao_supply_cereals_kg" = "^Cereals - Excluding Beer$",
  "fao_supply_roots_kg"   = "^Starchy Roots$",
  "fao_supply_pulses_kg"  = "^Pulses$",
  "fao_supply_veg_kg"     = "^Vegetables$",
  "fao_supply_fruit_kg"   = "^Fruits - Excluding Wine$",
  "fao_supply_animal_kg"  = "^Animal Products$")
ELEM_QTY <- "Food supply quantity \\(kg/capita/yr\\)"

cat("## FAOSTAT Food Balance Sheets (FBS)\n")
zip_path <- file.path(RAW_DIR, "FBS_normalized.zip")
if (!file.exists(zip_path)) {
  cat("  downloading FBS bulk (~55 MB) ...\n")
  old <- options(timeout = 3600); on.exit(options(old), add = TRUE)
  utils::download.file(FBS_URL, zip_path, mode = "wb", quiet = TRUE)
}
members <- utils::unzip(zip_path, list = TRUE)$Name
csv <- grep("Normalized.*\\.csv$", members, value = TRUE, ignore.case = TRUE)[1]
if (is.na(csv)) csv <- grep("\\.csv$", members, value = TRUE)[1]
csv_path <- file.path(RAW_DIR, basename(csv))
if (!file.exists(csv_path)) utils::unzip(zip_path, files = csv, exdir = RAW_DIR, junkpaths = TRUE)
cat("  reading", basename(csv), "...\n")
fbs <- as.data.frame(data.table::fread(csv_path, encoding = "Latin-1", showProgress = FALSE))
# normalise column names across FAOSTAT versions
nm <- tolower(names(fbs))
col <- function(...) { for (k in c(...)) { i <- match(k, nm); if (!is.na(i)) return(names(fbs)[i]) }; NA }
c_area <- col("area"); c_item <- col("item"); c_elem <- col("element")
c_year <- col("year"); c_val <- col("value")
if (any(is.na(c(c_area,c_item,c_elem,c_year,c_val))))
  stop("Unexpected FBS columns: ", paste(names(fbs), collapse=", "))

source(here::here("R","data_prep.R"))  # load_gadm_cached()

pick_year <- function(d, yr) { d <- d[is.finite(d[[c_val]]), ]; if(!nrow(d)) return(NULL)
  d[which.min(abs(d[[c_year]] - yr)), , drop=FALSE] }

for (cn in countries) {
  cat("\n#####", cn, "#####\n")
  sub <- fbs[fbs[[c_area]] == FAO_AREA[[cn]], ]
  if (!nrow(sub)) { cat("  no FBS rows for", FAO_AREA[[cn]], "— check area name\n"); next }
  yr <- SURVEY_YEAR[[cn]]; vals <- list()

  for (out in names(ELEM_NUTRIENT)) {
    d <- sub[grepl(ELEM_NUTRIENT[[out]], sub[[c_elem]]) &
             grepl("Grand Total|Total", sub[[c_item]], ignore.case=TRUE), ]
    r <- pick_year(d, yr); if (!is.null(r)) { vals[[out]] <- r[[c_val]]; vals[["fao_year"]] <- r[[c_year]] }
  }
  qty <- sub[grepl(ELEM_QTY, sub[[c_elem]]), ]
  for (out in names(ITEM_GROUP)) {
    d <- qty[grepl(ITEM_GROUP[[out]], qty[[c_item]]), ]
    r <- pick_year(d, yr); if (!is.null(r)) vals[[out]] <- r[[c_val]]
  }
  if (!length(vals)) { cat("  no matching elements/items; available elements:\n    ",
      paste(head(unique(sub[[c_elem]]),8), collapse=" | "), "\n"); next }
  cat(sprintf("  matched %d fields (year used: %s)\n", sum(grepl("^fao_(kcal|protein|fat|supply)", names(vals))),
              if (!is.null(vals$fao_year)) vals$fao_year else "NA"))

  gadm <- tryCatch(load_gadm_cached(ISO3[[cn]], level=2), error=function(e) NULL)
  if (is.null(gadm)) { cat("  GADM load failed; skipping\n"); next }
  g <- unique(sf::st_drop_geometry(gadm)[, intersect(c("NAME_1","NAME_2"), names(gadm)), drop=FALSE])
  names(g)[names(g)=="NAME_1"] <- "Admin1"; names(g)[names(g)=="NAME_2"] <- "Admin2"
  g$Admin1 <- trimws(g$Admin1); g$Admin2 <- trimws(g$Admin2)
  for (k in names(vals)) g[[k]] <- vals[[k]]
  fn <- file.path(OUT_DIR, sprintf("%s_fao_admin2.csv", cn))
  write.csv(g, fn, row.names=FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows, %d fao_ cols)\n",
              basename(fn), nrow(g), sum(grepl("^fao_", names(g)))))
}
cat("\nDONE. Domain FAO=list(prefix=\"fao_\") + guarded merge join are wired in config/merges.\n")
