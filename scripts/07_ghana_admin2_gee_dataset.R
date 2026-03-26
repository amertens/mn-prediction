# =============================================================================
# scripts/07_ghana_admin2_gee_dataset.R
#
# Build Admin1- and Admin2-level analytic datasets for Ghana by:
#   1. Computing survey-weighted micronutrient deficiency prevalence
#      (Vitamin A, children) at Admin1 and Admin2 using the stratified,
#      clustered survey design (PSU = gw_cnum, strata = gw_Strata,
#      weight = gw_sWeight).
#   2. Extracting ALL GEE rasters from data/Ghana_GEE rasters/,
#      handling multi-layer .tif files by extracting each layer separately,
#      clipped to Admin1/Admin2 boundaries and averaged over each unit.
#   3. Merging survey prevalence with raster summaries into analytic
#      datasets ready for ML prediction models.
#
# Inputs:
#   - data/IPD/Ghana/Ghana_merged_dataset.rds  (individual-level survey)
#   - data/Ghana_GEE rasters/*.tif  (all rasters in folder)
#
# Outputs:
#   - data/Ghana_admin1_analytic.rds / .csv / _sf.rds
#   - data/Ghana_admin2_analytic.rds / .csv / _sf.rds
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(sf)
  library(terra)
  library(geodata)
  library(exactextractr)
  library(readr)
  library(srvyr)
  library(survey)
})

options(survey.lonely.psu = "adjust")

cat("=== Ghana Admin1 + Admin2 GEE Dataset Builder ===\n\n")

# =============================================================================
# §01  LOAD SURVEY DATA & COMPUTE SURVEY-WEIGHTED PREVALENCE
# =============================================================================
cat("[01] Loading survey data...\n")

d <- readRDS(here::here("data", "IPD", "Ghana", "Ghana_merged_dataset.rds"))
cat(sprintf("  %d individuals × %d columns\n", nrow(d), ncol(d)))
cat(sprintf("  Survey year: %s\n", unique(d$gw_year)))

# Strip haven labels to plain numeric for key columns
d$gw_cVAD          <- as.numeric(unclass(d$gw_cVAD))
d$gw_cRBPAdjBrinda <- as.numeric(unclass(d$gw_cRBPAdjBrinda))
d$gw_sWeight       <- as.numeric(unclass(d$gw_sWeight))
d$gw_cnum          <- as.factor(d$gw_cnum)
d$gw_Strata        <- as.factor(d$gw_Strata)
d$Admin1           <- as.character(d$Admin1)
d$Admin2           <- as.character(d$Admin2)

# Children with Vitamin A data (gw_cVAD is non-NA for children)
d_children <- d %>%
  filter(!is.na(gw_cVAD), !is.na(Admin2), !is.na(Admin1),
         !is.na(gw_sWeight), !is.na(gw_cnum), !is.na(gw_Strata))

cat(sprintf("  Children with complete VAD + survey design data: %d\n",
            nrow(d_children)))

# ── Build survey design ──────────────────────────────────────────────────────
cat("  Building survey design (stratified, clustered, weighted)...\n")
svy_des <- d_children %>%
  srvyr::as_survey_design(
    ids     = gw_cnum,
    strata  = gw_Strata,
    weights = gw_sWeight,
    nest    = TRUE
  )

# ── National-level weighted prevalence ────────────────────────────────────────
svy_national <- svy_des %>%
  summarise(
    svy_prev = survey_mean(gw_cVAD, vartype = c("se", "ci"), na.rm = TRUE),
    n_children = n()
  )
cat(sprintf("  National weighted VAD prevalence: %.1f%% [%.1f%%, %.1f%%] (n=%d)\n",
            svy_national$svy_prev * 100,
            svy_national$svy_prev_low * 100,
            svy_national$svy_prev_upp * 100,
            svy_national$n_children))

# ── Admin1-level weighted prevalence ──────────────────────────────────────────
admin1_prev <- svy_des %>%
  group_by(Admin1) %>%
  summarise(
    svy_prev_vad = survey_mean(gw_cVAD, vartype = c("se", "ci"), na.rm = TRUE),
    raw_prev_vad = pick(gw_cVAD) %>% pull(gw_cVAD) %>% mean(na.rm = TRUE),
    mean_rbp     = pick(gw_cRBPAdjBrinda) %>% pull(gw_cRBPAdjBrinda) %>%
                     mean(na.rm = TRUE),
    n_children   = n(),
    .groups = "drop"
  ) %>%
  mutate(n_vad = round(svy_prev_vad * n_children))

cat(sprintf("  Admin1 regions with data: %d\n", nrow(admin1_prev)))

# ── Admin2-level weighted prevalence ──────────────────────────────────────────
admin2_prev <- svy_des %>%
  group_by(Admin1, Admin2) %>%
  summarise(
    svy_prev_vad = survey_mean(gw_cVAD, vartype = c("se", "ci"), na.rm = TRUE),
    raw_prev_vad = pick(gw_cVAD) %>% pull(gw_cVAD) %>% mean(na.rm = TRUE),
    mean_rbp     = pick(gw_cRBPAdjBrinda) %>% pull(gw_cRBPAdjBrinda) %>%
                     mean(na.rm = TRUE),
    n_children   = n(),
    .groups = "drop"
  ) %>%
  mutate(n_vad = round(svy_prev_vad * n_children))

cat(sprintf("  Admin2 districts with data: %d\n", nrow(admin2_prev)))


# =============================================================================
# §02  DOWNLOAD GADM BOUNDARIES (Admin1 + Admin2)
# =============================================================================
cat("\n[02] Downloading GADM boundaries...\n")

gadm_path <- tempdir()

poly_a1 <- geodata::gadm(country = "GHA", level = 1, path = gadm_path)
poly_a1 <- sf::st_as_sf(poly_a1) %>%
  select(NAME_1) %>%
  rename(Admin1_gadm = NAME_1) %>%
  sf::st_transform(crs = 4326)
cat(sprintf("  %d Admin1 polygons\n", nrow(poly_a1)))

poly_a2 <- geodata::gadm(country = "GHA", level = 2, path = gadm_path)
poly_a2 <- sf::st_as_sf(poly_a2) %>%
  select(NAME_1, NAME_2) %>%
  rename(Admin1_gadm = NAME_1, Admin2_gadm = NAME_2) %>%
  sf::st_transform(crs = 4326)
cat(sprintf("  %d Admin2 polygons\n", nrow(poly_a2)))


# =============================================================================
# §03  EXTRACT ALL GEE RASTERS AT ADMIN1 AND ADMIN2 LEVEL
# =============================================================================
cat("\n[03] Extracting all GEE rasters...\n")

raster_dir <- here::here("data", "Ghana_GEE rasters")
tif_files  <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
cat(sprintf("  Found %d .tif files in %s\n", length(tif_files), basename(raster_dir)))

# Helper: clean a filename + layer name into a valid R variable name
make_varname <- function(filename, layer_name = NULL, layer_idx = NULL, nlayers = 1) {
  # Strip path and extension
  base <- tools::file_path_sans_ext(basename(filename))
  # Remove "_Ghana" suffix for brevity
  base <- sub("_Ghana$", "", base)
  # Convert to lowercase and replace non-alnum with underscore
  base <- tolower(gsub("[^A-Za-z0-9]+", "_", base))
  base <- sub("^_|_$", "", base)


  if (nlayers == 1) return(paste0("gee_", base))

  # Multi-layer: append cleaned layer name
  if (!is.null(layer_name) && nzchar(layer_name) && layer_name != "b1") {
    lyr <- tolower(gsub("[^A-Za-z0-9]+", "_", layer_name))
    lyr <- sub("^_|_$", "", lyr)
    # Truncate very long layer names (GPW demographic has ~100-char names)
    if (nchar(lyr) > 60) lyr <- substr(lyr, 1, 60)
    return(paste0("gee_", base, "__", lyr))
  }

  # Fallback: use numeric index
  paste0("gee_", base, "__lyr", layer_idx)
}

# Helper: extract zonal means from a single .tif, handling multiple layers.
# Returns a named list of vectors (one per layer), each length = nrow(polygons).
extract_tif_zonal <- function(tif_path, polygons) {
  r  <- tryCatch(terra::rast(tif_path), error = function(e) {
    warning(sprintf("  Could not read %s: %s", basename(tif_path), e$message))
    return(NULL)
  })
  if (is.null(r)) return(NULL)

  nl <- terra::nlyr(r)
  layer_names <- names(r)

  cat(sprintf("    %s | %d layer(s) | res: %.5f\n",
              basename(tif_path), nl, terra::res(r)[1]))

  result <- list()
  for (i in seq_len(nl)) {
    varname <- make_varname(tif_path, layer_names[i], i, nl)
    vals <- tryCatch(
      exactextractr::exact_extract(r[[i]], polygons, fun = "mean"),
      error = function(e) {
        warning(sprintf("    Extract failed layer %d (%s): %s",
                        i, layer_names[i], e$message))
        rep(NA_real_, nrow(polygons))
      }
    )
    result[[varname]] <- vals
  }
  result
}

# ── Extract for Admin2 and Admin1 ────────────────────────────────────────────
raster_admin2 <- data.frame(Admin2_gadm = poly_a2$Admin2_gadm)
raster_admin1 <- data.frame(Admin1_gadm = poly_a1$Admin1_gadm)

for (idx in seq_along(tif_files)) {
  tif <- tif_files[idx]
  cat(sprintf("  [%d/%d] %s\n", idx, length(tif_files), basename(tif)))

  # Admin2
  a2_vals <- extract_tif_zonal(tif, poly_a2)
  if (!is.null(a2_vals)) {
    for (nm in names(a2_vals)) raster_admin2[[nm]] <- a2_vals[[nm]]
  }

  # Admin1
  a1_vals <- extract_tif_zonal(tif, poly_a1)
  if (!is.null(a1_vals)) {
    for (nm in names(a1_vals)) raster_admin1[[nm]] <- a1_vals[[nm]]
  }

  # Free memory after large rasters
  gc(verbose = FALSE)
}

n_gee_vars <- ncol(raster_admin2) - 1
cat(sprintf("\n  Raster extraction complete: %d GEE variables extracted\n", n_gee_vars))
cat(sprintf("    Admin2: %d polygons × %d vars\n", nrow(raster_admin2), n_gee_vars))
cat(sprintf("    Admin1: %d polygons × %d vars\n", nrow(raster_admin1), n_gee_vars))


# =============================================================================
# §04  FUZZY-MATCH SURVEY NAMES TO GADM NAMES & MERGE
# =============================================================================
cat("\n[04] Matching survey names to GADM polygons...\n")

match_names <- function(survey_nms, gadm_nms) {
  do.call(rbind, lapply(survey_nms, function(nm) {
    dists <- adist(nm, gadm_nms, ignore.case = TRUE)[1, ]
    best  <- which.min(dists)
    data.frame(survey_name   = nm,
               gadm_name     = gadm_nms[best],
               edit_distance = dists[best],
               stringsAsFactors = FALSE)
  }))
}

# ── Admin2 matching ──────────────────────────────────────────────────────────
name_map_a2 <- match_names(unique(admin2_prev$Admin2), poly_a2$Admin2_gadm)
imperfect_a2 <- name_map_a2[name_map_a2$edit_distance > 0, ]
if (nrow(imperfect_a2) > 0) {
  cat("  Admin2 imperfect matches:\n"); print(imperfect_a2, row.names = FALSE)
} else {
  cat("  Admin2: all names matched exactly.\n")
}

admin2_merged <- admin2_prev %>%
  left_join(name_map_a2, by = c("Admin2" = "survey_name")) %>%
  rename(Admin2_gadm = gadm_name) %>%
  select(-edit_distance) %>%
  left_join(raster_admin2, by = "Admin2_gadm")

admin2_sf <- poly_a2 %>%
  left_join(admin2_merged, by = "Admin2_gadm")

# ── Admin1 matching ──────────────────────────────────────────────────────────
name_map_a1 <- match_names(unique(admin1_prev$Admin1), poly_a1$Admin1_gadm)
imperfect_a1 <- name_map_a1[name_map_a1$edit_distance > 0, ]
if (nrow(imperfect_a1) > 0) {
  cat("  Admin1 imperfect matches:\n"); print(imperfect_a1, row.names = FALSE)
} else {
  cat("  Admin1: all names matched exactly.\n")
}

admin1_merged <- admin1_prev %>%
  left_join(name_map_a1, by = c("Admin1" = "survey_name")) %>%
  rename(Admin1_gadm = gadm_name) %>%
  select(-edit_distance) %>%
  left_join(raster_admin1, by = "Admin1_gadm")

admin1_sf <- poly_a1 %>%
  left_join(admin1_merged, by = "Admin1_gadm")


# =============================================================================
# §05  SAVE ANALYTIC DATASETS
# =============================================================================
cat("\n[05] Saving analytic datasets...\n")

# Identify GEE columns (all columns starting with gee_)
gee_cols <- grep("^gee_", colnames(raster_admin2), value = TRUE)

# ── Admin2 flat ──────────────────────────────────────────────────────────────
admin2_flat <- admin2_merged %>%
  select(Admin1, Admin2, Admin2_gadm, n_children, n_vad,
         svy_prev_vad, svy_prev_vad_se, svy_prev_vad_low, svy_prev_vad_upp,
         raw_prev_vad, mean_rbp,
         all_of(gee_cols))

saveRDS(admin2_flat, here::here("data", "Ghana_admin2_analytic.rds"))
readr::write_csv(admin2_flat, here::here("data", "Ghana_admin2_analytic.csv"))
saveRDS(admin2_sf,   here::here("data", "Ghana_admin2_analytic_sf.rds"))

cat(sprintf("  Admin2: %d rows × %d cols (%d survey + %d GEE)\n",
            nrow(admin2_flat), ncol(admin2_flat),
            ncol(admin2_flat) - length(gee_cols), length(gee_cols)))

# ── Admin1 flat ──────────────────────────────────────────────────────────────
admin1_flat <- admin1_merged %>%
  select(Admin1, Admin1_gadm, n_children, n_vad,
         svy_prev_vad, svy_prev_vad_se, svy_prev_vad_low, svy_prev_vad_upp,
         raw_prev_vad, mean_rbp,
         all_of(gee_cols))

saveRDS(admin1_flat, here::here("data", "Ghana_admin1_analytic.rds"))
readr::write_csv(admin1_flat, here::here("data", "Ghana_admin1_analytic.csv"))
saveRDS(admin1_sf,   here::here("data", "Ghana_admin1_analytic_sf.rds"))

cat(sprintf("  Admin1: %d rows × %d cols (%d survey + %d GEE)\n",
            nrow(admin1_flat), ncol(admin1_flat),
            ncol(admin1_flat) - length(gee_cols), length(gee_cols)))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat(sprintf("  National weighted VAD prevalence: %.1f%%\n",
            svy_national$svy_prev * 100))
cat(sprintf("  GEE variables extracted: %d\n", length(gee_cols)))
cat(sprintf("  Admin1 regions: %d / %d GADM | VAD range: %.1f%% – %.1f%%\n",
            sum(!is.na(admin1_sf$svy_prev_vad)), nrow(poly_a1),
            min(admin1_flat$svy_prev_vad, na.rm = TRUE) * 100,
            max(admin1_flat$svy_prev_vad, na.rm = TRUE) * 100))
cat(sprintf("  Admin2 districts: %d / %d GADM | VAD range: %.1f%% – %.1f%%\n",
            sum(!is.na(admin2_sf$svy_prev_vad)), nrow(poly_a2),
            min(admin2_flat$svy_prev_vad, na.rm = TRUE) * 100,
            max(admin2_flat$svy_prev_vad, na.rm = TRUE) * 100))

cat("\nDone.\n")
