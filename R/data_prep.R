# =============================================================================
# R/data_prep.R
#
# Functions for loading and constructing analysis-ready datasets.
# Wrapped as pure functions so targets can track inputs/outputs.
# =============================================================================


#' Load GADM admin boundaries with local caching
#'
#' Downloads once from geodata::gadm(), saves as RDS in data/admin_boundaries/,
#' then loads from cache on subsequent runs.
#'
#' @param country_code ISO country code (e.g., "GM", "GH", "SLE")
#' @param level Admin level (1 or 2)
#' @param cache_dir Directory to cache boundaries (default: data/admin_boundaries/)
#' @return sf object with admin boundaries
load_gadm_cached <- function(country_code, level = 2,
                             cache_dir = here::here("data", "admin_boundaries")) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, sprintf("gadm41_%s_%d.rds", country_code, level))

  if (file.exists(cache_file)) {
    cat(sprintf("[load_gadm_cached] Loading cached: %s\n", basename(cache_file)))
    return(readRDS(cache_file))
  }

  cat(sprintf("[load_gadm_cached] Downloading GADM %s level %d...\n", country_code, level))
  poly <- geodata::gadm(country = country_code, level = level, path = tempdir())
  poly_sf <- sf::st_as_sf(poly)
  saveRDS(poly_sf, cache_file)
  cat(sprintf("[load_gadm_cached] Cached to: %s\n", basename(cache_file)))
  poly_sf
}

#' Load the merged country dataset from disk and harmonize columns
#'
#' Handles cross-country differences:
#' - Derives gw_child_flag if missing (Ghana) from child biomarker presence
#' - Derives binary VAD columns if missing (Sierra Leone) from continuous RBP
#' - Malawi uses raw biomarker names (rbp, fer, vitA_def) and a text
#'   `population` column — no harmonization needed, just load.
#'
#' @param data_path Path to the .rds file
#' @return data.frame with harmonized columns
load_merged_data <- function(data_path) {
  stopifnot(file.exists(data_path))
  d <- readRDS(data_path)
  cat(sprintf("[data_prep] Loaded %s: %d rows x %d cols\n",
              basename(data_path), nrow(d), ncol(d)))

  # ── Detect Malawi by its column naming pattern ─────────────────────────
  # Malawi has a text `population` column and no gw_-prefixed biomarkers.
  # The config handles population filtering via child_flag_val text matching,
  # so no harmonization is needed here.
  is_malawi <- "population" %in% colnames(d) &&
    is.character(d[["population"]]) &&
    any(d[["population"]] == "preschool children", na.rm = TRUE)

  if (is_malawi) {
    cat("  Malawi MNS data detected — using raw biomarker columns\n")
    # Strip haven labels even for Malawi
    for (col in colnames(d)) {
      if (inherits(d[[col]], "haven_labelled")) {
        d[[col]] <- as.double(unclass(haven::zap_labels(d[[col]])))
      }
    }
    return(d)
  }

  # ── Derive gw_child_flag if missing ────────────────────────────────────
  # Ghana doesn't have this column. Derive from child-specific biomarkers:
  # rows with non-NA child RBP (gw_cRBPAdjThurn or gw_cRBPAdj) are children.
  if (!"gw_child_flag" %in% colnames(d)) {
    child_rbp <- NULL
    for (col in c("gw_cRBPAdjThurn", "gw_cRBPAdjBrinda", "gw_cRBPAdj",
                  "gw_cRBP", "gw_childid")) {
      if (col %in% colnames(d) && any(!is.na(d[[col]]))) {
        child_rbp <- col
        break
      }
    }
    if (!is.null(child_rbp)) {
      d$gw_child_flag <- ifelse(!is.na(d[[child_rbp]]), 1L, 0L)
      cat(sprintf("  Derived gw_child_flag from %s: %d children, %d women\n",
                  child_rbp, sum(d$gw_child_flag == 1), sum(d$gw_child_flag == 0)))
    } else {
      warning("Cannot derive gw_child_flag — no child biomarker column found")
    }
  }

  # ── Derive binary VAD columns if missing ────────────────────────────────
  # Sierra Leone has continuous RBP but no pre-computed binary VAD column.
  # Derive as RBP < 0.70 µmol/L (standard WHO threshold).
  derive_binary <- function(d, cont_col, bin_col, threshold = 0.70) {
    if (cont_col %in% colnames(d) && !(bin_col %in% colnames(d))) {
      d[[bin_col]] <- ifelse(d[[cont_col]] < threshold, 1L, 0L)
      d[[bin_col]][is.na(d[[cont_col]])] <- NA_integer_
      n_def <- sum(d[[bin_col]] == 1, na.rm = TRUE)
      n_tot <- sum(!is.na(d[[bin_col]]))
      cat(sprintf("  Derived %s from %s < %.2f: %d/%d deficient\n",
                  bin_col, cont_col, threshold, n_def, n_tot))
    }
    d
  }

  d <- derive_binary(d, "gw_cRBPAdj", "gw_cVAD", 0.70)
  d <- derive_binary(d, "gw_wRBPAdj", "gw_wVAD", 0.70)

  # ── Strip haven labels from all columns ──────────────────────────────
  # haven_labelled columns from .dta files cause arithmetic errors downstream
  # (e.g., <double> * <haven_labelled> is not permitted). Also convert
  # user-defined NA values (.a, .b etc.) to real NA.
  n_haven <- 0L
  for (col in colnames(d)) {
    if (inherits(d[[col]], "haven_labelled")) {
      d[[col]] <- as.double(unclass(haven::zap_labels(d[[col]])))
      n_haven <- n_haven + 1L
    }
  }
  if (n_haven > 0) {
    cat(sprintf("  Stripped haven labels from %d columns\n", n_haven))
  }

  d
}


#' Load DHS admin-1 estimates and pivot to wide format
#'
#' Converts the long-format `*_dhs_admin1_direct.rds` file produced by
#' DHS_admin2_aggregation.R into the legacy wide format that merge scripts
#' expect: one row per admin-1 region, with columns named
#' `dhs{YEAR}_{indicator}` and a key column `dhs{YEAR}_DHSREGEN`.
#'
#' Falls back to the old `*_dhs_aggregation.rds` file if it still exists.
#'
#' @param dhs_dir  Path to data/DHS/clean/ directory
#' @param country  Country name (e.g., "Gambia")
#' @param year     DHS survey year (e.g., 2019)
#' @return data.frame in wide format, or NULL if no file found
load_dhs_admin1 <- function(dhs_dir, country, year) {

  prefix <- paste0("dhs", year, "_")


  # --- Try new long-format file first ---
  long_path <- file.path(dhs_dir, paste0(country, "_", year, "_dhs_admin1_direct.rds"))
  if (file.exists(long_path)) {
    long_df <- readRDS(long_path)
    cat(sprintf("[load_dhs_admin1] Loaded long format: %s (%d rows, %d indicators)\n",
                basename(long_path), nrow(long_df), length(unique(long_df$indicator))))

    # Pivot to wide: one row per admin1, columns = dhs{YEAR}_{indicator}
    wide <- long_df %>%
      dplyr::select(admin1.name, indicator, direct.est) %>%
      tidyr::pivot_wider(
        names_from  = indicator,
        values_from = direct.est,
        names_prefix = prefix
      ) %>%
      dplyr::rename(!!paste0(prefix, "DHSREGEN") := admin1.name)

    cat(sprintf("[load_dhs_admin1] Pivoted to wide: %d regions x %d columns\n",
                nrow(wide), ncol(wide)))
    return(as.data.frame(wide))
  }

  # --- Fallback: old legacy wide file ---
  # Try multiple naming conventions used across countries
  legacy_candidates <- c(
    file.path(dhs_dir, paste0(country, "_", year, "_dhs_aggregation.rds")),
    file.path(dhs_dir, paste0("SL_", year, "_dhs_aggregation.rds")),     # Sierra Leone
    file.path(dhs_dir, paste0(country, "_", year, "_dhs_aggregation.rds"))
  )
  for (lp in unique(legacy_candidates)) {
    if (file.exists(lp)) {
      wide <- readRDS(lp)
      cat(sprintf("[load_dhs_admin1] Loaded legacy file: %s (%d regions x %d cols)\n",
                  basename(lp), nrow(wide), ncol(wide)))
      # Add prefix if not already present
      if (!any(grepl(paste0("^", prefix), names(wide)))) {
        names(wide) <- paste0(prefix, names(wide))
      }
      return(as.data.frame(wide))
    }
  }

  warning(sprintf("[load_dhs_admin1] No admin-1 DHS file found for %s %d", country, year))
  NULL
}


#' Load DHS admin-2 estimates for merging into the pipeline dataset
#'
#' Loads the pre-built wide .rds file produced by DHS_admin2_aggregation.R.
#' Column names follow the convention {prefix}{indicator}_adm2, e.g.,
#' "dhs2019_womananemia_adm2". The Admin2 column matches GADM NAME_2.
#'
#' @param dhs_dir Path to data/DHS/clean/ directory
#' @param country Country name (e.g., "Gambia")
#' @param year DHS survey year (e.g., 2019)
#' @param merge_col Name of the Admin2 column in the target dataset (default "Admin2")
#' @return data.frame with merge_col + indicator columns, or NULL if file missing
load_dhs_admin2 <- function(dhs_dir, country, year, merge_col = "Admin2") {

  wide_path   <- file.path(dhs_dir, paste0(country, "_", year, "_dhs_admin2_wide.rds"))
  custom_path <- file.path(dhs_dir, paste0(country, "_", year, "_dhs_custom_admin2_wide.rds"))

  wide <- NULL
  custom <- NULL

  # Load surveyPrev-based admin-2 estimates
  if (file.exists(wide_path)) {
    wide <- readRDS(wide_path)
    cat(sprintf("[load_dhs_admin2] Loaded surveyPrev: %s (%d areas x %d cols)\n",
                basename(wide_path), nrow(wide), ncol(wide) - 1L))
  }


  # Load custom indicator admin-2 estimates
  if (file.exists(custom_path)) {
    custom <- readRDS(custom_path)
    cat(sprintf("[load_dhs_admin2] Loaded custom: %s (%d areas x %d cols)\n",
                basename(custom_path), nrow(custom), ncol(custom) - 1L))
  }

  # Merge both if available
  if (!is.null(wide) && !is.null(custom)) {
    # Remove duplicate indicator columns (custom takes precedence if overlap)
    overlap_cols <- setdiff(intersect(names(wide), names(custom)), "Admin2")
    if (length(overlap_cols) > 0) {
      cat(sprintf("[load_dhs_admin2] %d overlapping columns — keeping custom version\n",
                  length(overlap_cols)))
      wide <- wide[, !(names(wide) %in% overlap_cols), drop = FALSE]
    }
    combined <- merge(wide, custom, by = "Admin2", all = TRUE)
    cat(sprintf("[load_dhs_admin2] Combined: %d Admin2 areas x %d indicator columns\n",
                nrow(combined), ncol(combined) - 1L))
    return(combined)
  }

  if (!is.null(wide)) return(wide)
  if (!is.null(custom)) return(custom)

  warning(sprintf("[load_dhs_admin2] No admin-2 DHS files found for %s %d.\n  Run src/DHS/DHS_admin2_aggregation.R and/or DHS_custom_admin2_indicators.R first.",
                  country, year))
  NULL
}


#' Build an outcome-specific dataset (one population x one micronutrient)
#'
#' Filters to the correct population, selects predictors, removes leakage
#' columns, and drops rows with missing outcome.
#'
#' @param merged_data The full merged dataset
#' @param cc Country config (from get_country_configs())
#' @param oc Outcome config (one element of cc$outcomes)
#' @return list with components: data, Xvars_full, Xvars, domain_vars
build_outcome_dataset <- function(merged_data, cc, oc) {

  d <- merged_data

  # Filter to population
  pop_col <- cc$child_flag
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    d <- d[d[[pop_col]] == oc$child_flag_val, ]
  }

  # Require non-missing outcome
  if (!is.null(oc$continuous) && oc$continuous %in% colnames(d))
    d <- d[!is.na(d[[oc$continuous]]), ]
  if (!is.null(oc$binary) && oc$binary %in% colnames(d))
    d <- d[!is.na(d[[oc$binary]]), ]

  # Recode binary outcome if coded as factor-like 1/2 instead of 0/1.
  # Some surveys (e.g., Sierra Leone MNS) code deficiency as 1=yes, 2=no.
  # Convert to standard 0/1 where 1=deficient, 0=not deficient.
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    yb <- d[[oc$binary]]
    yb_vals <- sort(unique(as.numeric(yb[!is.na(yb)])))
    if (length(yb_vals) == 2 && all(yb_vals == c(1, 2))) {
      cat(sprintf("  Recoding %s from {1,2} to {1,0} (1=deficient, 2=not → 0)\n",
                  oc$binary))
      d[[oc$binary]] <- ifelse(d[[oc$binary]] == 1, 1, 0)
    }
  }

  cat(sprintf("[build_outcome] %s — %s: %d observations\n",
              cc$country, oc$tag, nrow(d)))

  # Identify predictor columns by domain prefix
  all_cols <- colnames(d)
  domain_vars <- list()
  Xvars_all <- character()

  # Add external predictor domains if the columns exist in the data
  # (merged via merge_external_predictors)
  external_domains <- list(
    CHIRPS    = list(prefix = "chirps_"),
    WorldPop  = list(prefix = "worldpop_"),
    MAP2      = list(prefix = "map2_"),
    Soil      = list(prefix = "soil_"),
    Nightlights = list(prefix = "ntl_"),
    GDL       = list(prefix = "gdl_"),
    Crop      = list(prefix = "crop_"),
    IPC       = list(prefix = "ipc_"),
    ACLED     = list(prefix = "acled_"),
    WFP_ext   = list(prefix = "wfp_")
  )

  # Merge external domains into config domains (only if columns exist)
  all_domains <- cc$domains
  for (nm in names(external_domains)) {
    pfx <- external_domains[[nm]]$prefix
    if (any(startsWith(all_cols, pfx))) {
      all_domains[[nm]] <- external_domains[[nm]]
    }
  }

  for (nm in names(all_domains)) {
    dom <- all_domains[[nm]]
    prefix_cols <- all_cols[startsWith(all_cols, dom$prefix)]
    if (!is.null(dom$extra)) prefix_cols <- c(prefix_cols, intersect(dom$extra, all_cols))

    # Remove leakage columns from GW domain
    if (nm == "GW" && !is.null(cc$gw_exclude_patterns)) {
      leak <- Reduce("|", lapply(cc$gw_exclude_patterns,
                                  function(p) grepl(p, prefix_cols, ignore.case = TRUE)))
      prefix_cols <- prefix_cols[!leak]

      # Also exclude the outcome columns themselves
      exclude_exact <- c(oc$continuous, oc$binary, cc$admin1_col, cc$admin2_col,
                         cc$cluster_id, cc$weight_col, cc$strata_col, cc$child_flag)
      prefix_cols <- setdiff(prefix_cols, exclude_exact)
    }

    domain_vars[[nm]] <- prefix_cols
    Xvars_all <- c(Xvars_all, prefix_cols)
  }

  Xvars_full <- unique(Xvars_all)
  Xvars_no_gw <- unique(setdiff(Xvars_all, domain_vars[["GW"]]))

  # Keep only columns we need
  keep_cols <- unique(c(
    oc$continuous, oc$binary,
    cc$cluster_id, cc$admin1_col, cc$admin2_col,
    cc$weight_col, cc$strata_col, cc$psu_col,
    Xvars_full
  ))
  keep_cols <- intersect(keep_cols, colnames(d))
  d <- d[, keep_cols, drop = FALSE]

  list(
    data        = d,
    Xvars_full  = Xvars_full,
    Xvars       = Xvars_no_gw,
    domain_vars = domain_vars
  )
}
