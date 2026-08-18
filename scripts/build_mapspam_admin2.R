#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_mapspam_admin2.R
#
# Build MapSPAM crop-composition predictors aggregated to GADM Admin-2, writing
# data/MapSPAM/<Country>_spam_admin2.csv with pipeline-ready `spam_` columns.
#
# WHY: the pipeline has cropland *class* (land cover) and productivity (GEE GPP)
# but not crop *composition*. The mix of staples vs pulses (iron/zinc) vs
# vegetables/roots (provitamin A) maps to local micronutrient availability, and
# it is SUBNATIONAL — so it can help within-country ranking AND transfer.
# See docs/ra_new_data_sources_specs.md §4.
#
# SOURCE / ACCESS (public, no key, no extra package):
#   MapSPAM 2010 v2r0 "Global Spatially-Disaggregated Crop Production Statistics"
#   Harvard Dataverse, DOI 10.7910/DVN/PRFF8V (IFPRI, 2019). 5-arcmin (~10 km)
#   grid, 42 crops, CSV/GeoTIFF. We use the all-technology (_TA) CSVs for
#   production (P) and physical area (A). Downloaded with base download.file()
#   from the Dataverse access API (verified 2026-07: HTTP 206, application/zip);
#   ~150 MB per zip, cached under data/MapSPAM/raw/. No API key or package needed.
#
# USAGE:
#   Rscript scripts/build_mapspam_admin2.R                # all countries
#   Rscript scripts/build_mapspam_admin2.R Ghana Tanzania # a subset
#
# The heavy download + unzip runs once; per-country processing reuses the cache.
# Written defensively: on first run it PRINTS the CSV filenames and detected crop
# columns — confirm them before trusting the output.
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(data.table); library(sf)
})

# ---- config ----------------------------------------------------------------
ISO3 <- c(Gambia = "GMB", Ghana = "GHA", SierraLeone = "SLE",
          Malawi = "MWI", Tanzania = "TZA")

# Harvard Dataverse access API (download by numeric file id; verified 2026-07).
DVERSE   <- "https://dataverse.harvard.edu/api/access/datafile/%d"
RAW_DIR  <- here::here("data", "MapSPAM", "raw")
OUT_DIR  <- here::here("data", "MapSPAM")
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Dataverse file names + numeric ids (dataset doi:10.7910/DVN/PRFF8V, v2r0):
FILES <- list(
  prod = list(zip = "spam2010v2r0_global_prod.csv.zip",      id = 3984975L),
  area = list(zip = "spam2010v2r0_global_phys_area.csv.zip", id = 3984973L)
)

# SPAM 42-crop codes grouped by micronutrient relevance -----------------------
CROP_GROUPS <- list(
  cereals    = c("whea","rice","maiz","barl","pmil","smil","sorg","ocer"),
  roots      = c("pota","swpo","yams","cass","orts"),
  pulses     = c("bean","chic","cowp","pige","lent","opul","soyb"),
  oilcrops   = c("grou","cnut","oilp","sunf","rape","sesa","ooil"),
  vegetables = c("vege")
)

args <- commandArgs(trailingOnly = TRUE)
countries <- if (length(args)) args else names(ISO3)
bad <- setdiff(countries, names(ISO3))
if (length(bad)) stop("Unknown country: ", paste(bad, collapse = ", "),
                      " (known: ", paste(names(ISO3), collapse = ", "), ")")

# ---- download + unzip the two zips once (cached) ---------------------------
fetch_and_unzip <- function(spec) {
  zip_path <- file.path(RAW_DIR, spec$zip)
  if (!file.exists(zip_path)) {
    url <- sprintf(DVERSE, spec$id)
    cat("  downloading", spec$zip, "(~150 MB) ...\n")
    old <- options(timeout = 3600); on.exit(options(old), add = TRUE)
    utils::download.file(url, zip_path, mode = "wb", quiet = TRUE)
  }
  # unzip any CSVs not yet extracted
  members <- utils::unzip(zip_path, list = TRUE)$Name
  csvs <- grep("\\.csv$", members, value = TRUE, ignore.case = TRUE)
  need <- csvs[!file.exists(file.path(RAW_DIR, basename(csvs)))]
  if (length(need)) utils::unzip(zip_path, files = need, exdir = RAW_DIR, junkpaths = TRUE)
  file.path(RAW_DIR, basename(csvs))
}

# Pick the all-technology (_TA) CSV for a given variable letter (P or A).
pick_ta_csv <- function(csv_paths, var_letter) {
  ta <- grep("_TA\\.csv$", csv_paths, value = TRUE, ignore.case = TRUE)
  # prefer the file carrying the variable letter, e.g. _P_TA / _A_TA
  vv <- grep(sprintf("_%s_TA\\.csv$", var_letter), ta, value = TRUE, ignore.case = TRUE)
  out <- if (length(vv)) vv[1] else if (length(ta)) ta[1] else csv_paths[1]
  out
}

cat("## MapSPAM 2010 v2r0 — Admin-2 crop composition\n")
prod_csvs <- fetch_and_unzip(FILES[["prod"]])
area_csvs <- fetch_and_unzip(FILES[["area"]])
prod_csv  <- pick_ta_csv(prod_csvs, "P")
area_csv  <- pick_ta_csv(area_csvs, "A")
cat(sprintf("  production CSV : %s\n  phys-area CSV  : %s\n",
            basename(prod_csv), basename(area_csv)))

# ---- read once, keep only cols we need, filter later per country -----------
# SPAM CSVs carry: iso3, x (lon), y (lat), alloc_key, plus one column per crop.
all_crops <- unique(unlist(CROP_GROUPS))
read_spam <- function(path) {
  hdr <- names(data.table::fread(path, nrows = 0))
  lc  <- tolower(hdr)
  iso_col <- hdr[match(TRUE, lc %in% c("iso3","cntry_code","cntrycode"))]
  x_col   <- hdr[match(TRUE, lc %in% c("x","cell_x","lon","longitude"))]
  y_col   <- hdr[match(TRUE, lc %in% c("y","cell_y","lat","latitude"))]
  if (any(is.na(c(iso_col, x_col, y_col))))
    stop("Could not find iso3/x/y columns in ", basename(path),
         " — headers: ", paste(head(hdr, 20), collapse = ", "))
  # SPAM v2r0 crop columns are "<code>_<tech>" (e.g. whea_a for all-technology);
  # strip the trailing "_<letter>" to match the bare 4-char crop codes, then
  # rename the kept columns to bare codes so downstream grouping is simple.
  bare_of  <- sub("_[a-z]$", "", lc)
  crop_raw <- hdr[bare_of %in% all_crops]
  keep <- c(iso_col, x_col, y_col, crop_raw)
  dt <- data.table::fread(path, select = keep)
  data.table::setnames(dt, c(iso_col, x_col, y_col), c("iso3", "x", "y"))
  bare_crops <- sub("_[a-z]$", "", tolower(crop_raw))
  data.table::setnames(dt, crop_raw, bare_crops)
  dt[, iso3 := toupper(trimws(iso3))]
  list(dt = dt, crop_cols = bare_crops)
}

cat("  reading production grid ...\n"); P <- read_spam(prod_csv)
cat("  reading physical-area grid ...\n"); A <- read_spam(area_csv)
cat(sprintf("  detected %d crop columns (e.g. %s)\n",
            length(P$crop_cols), paste(head(P$crop_cols, 8), collapse = ", ")))

source(here::here("R", "data_prep.R"))  # load_gadm_cached()

group_sum <- function(dt, crop_cols) {
  present <- intersect(unlist(CROP_GROUPS), tolower(crop_cols))
  # map actual (possibly mixed-case) column names by lowercase
  colmap <- setNames(crop_cols, tolower(crop_cols))
  for (g in names(CROP_GROUPS)) {
    cc <- colmap[intersect(CROP_GROUPS[[g]], names(colmap))]
    dt[[paste0("g_", g)]] <- if (length(cc)) rowSums(dt[, ..cc], na.rm = TRUE) else 0
  }
  dt
}

for (cn in countries) {
  cat("\n#####", cn, "(", ISO3[[cn]], ") #####\n")
  gadm <- tryCatch(load_gadm_cached(ISO3[[cn]], level = 2), error = function(e) NULL)
  if (is.null(gadm)) { cat("  !! GADM load failed; skipping\n"); next }
  gadm <- sf::st_as_sf(gadm)
  gkeep <- intersect(c("NAME_1", "NAME_2"), names(gadm))
  gadm  <- gadm[, gkeep]

  pc <- P$dt[iso3 == ISO3[[cn]]]
  ac <- A$dt[iso3 == ISO3[[cn]]]
  if (!nrow(pc)) { cat("  no SPAM pixels for", ISO3[[cn]], "\n"); next }
  pc <- group_sum(data.table::copy(pc), P$crop_cols)
  ac <- group_sum(data.table::copy(ac), A$crop_cols)

  # spatial join pixel centroids -> Admin-2 polygon
  assign_admin2 <- function(dt) {
    pts <- sf::st_as_sf(dt, coords = c("x", "y"), crs = 4326, remove = FALSE)
    idx <- sf::st_join(pts, gadm, join = sf::st_intersects)
    sf::st_drop_geometry(idx)
  }
  pj <- assign_admin2(pc); aj <- assign_admin2(ac)

  gcols <- paste0("g_", names(CROP_GROUPS))
  prod_a2 <- as.data.table(pj)[!is.na(NAME_2),
              lapply(.SD, sum, na.rm = TRUE), by = .(NAME_1, NAME_2), .SDcols = gcols]
  area_a2 <- as.data.table(aj)[!is.na(NAME_2),
              .(spam_parea_total = sum(rowSums(.SD, na.rm = TRUE))),
              by = .(NAME_1, NAME_2), .SDcols = gcols]

  data.table::setnames(prod_a2, gcols, paste0("spam_prod_", names(CROP_GROUPS)))
  prod_a2[, spam_prod_total := rowSums(.SD, na.rm = TRUE),
          .SDcols = paste0("spam_prod_", names(CROP_GROUPS))]
  for (g in names(CROP_GROUPS)) {
    prod_a2[[paste0("spam_share_", g)]] <-
      ifelse(prod_a2$spam_prod_total > 0,
             prod_a2[[paste0("spam_prod_", g)]] / prod_a2$spam_prod_total, NA_real_)
  }

  out <- merge(prod_a2, area_a2, by = c("NAME_1", "NAME_2"), all.x = TRUE)
  data.table::setnames(out, c("NAME_1", "NAME_2"), c("Admin1", "Admin2"))
  out <- as.data.frame(out)
  out$Admin1 <- trimws(out$Admin1); out$Admin2 <- trimws(out$Admin2)

  fn <- file.path(OUT_DIR, sprintf("%s_spam_admin2.csv", cn))
  write.csv(out, fn, row.names = FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows, %d spam_ cols)\n",
              basename(fn), nrow(out), sum(grepl("^spam_", names(out)))))
}

cat("\nDONE. Wire in R/config.R (SPAM = list(prefix = \"spam_\")) — already added —\n",
    "and the guarded join in each src/<Country>/2_GW_<Country>_data_merge.R (already added).\n")
