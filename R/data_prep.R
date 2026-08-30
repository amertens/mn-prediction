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

    # ── Derive folate deficiency (Malawi) ──────────────────────────────────
    # Malawi `fol` is serum folate ALREADY IN nmol/L (median ~17.2). Verified
    # against Qi et al. 2024 (same survey): 8.5% deficient at <7 nmol/L vs 9.3%
    # computed here from `fol` directly. Apply the WHO serum-folate deficiency
    # cutoff <10 nmol/L — the SAME definition Ghana and Sierra Leone use
    # (gw_wFolate < 10 nmol/L).
    # 2026-06-16 fixes: (a) the earlier <3 ng/mL survey definition was far
    # stricter than WHO; (b) a `fol * 2.266` ng/mL->nmol/L conversion was ALSO
    # wrong because `fol` is already nmol/L. Together these reported Malawi
    # folate deficiency as ~0% (an artifact) instead of the correct ~21% at
    # <10 nmol/L. The residual gap vs Ghana/SL (Malawi women have higher folate)
    # is real — consistent with Malawi's 2015 folic-acid flour fortification.
    if ("fol" %in% colnames(d)) {
      d$fol_nmol <- d$fol  # already nmol/L — do NOT multiply by 2.266
      if (!"folate_def" %in% colnames(d)) {
        d$folate_def <- ifelse(d$fol_nmol < 10, 1L, 0L)
        d$folate_def[is.na(d$fol_nmol)] <- NA_integer_
        n_def <- sum(d$folate_def == 1, na.rm = TRUE)
        n_tot <- sum(!is.na(d$folate_def))
        cat(sprintf("  Derived folate_def from fol_nmol < 10 nmol/L (WHO): %d/%d deficient\n", n_def, n_tot))
      }
    }

    # ── Derive vitamin B12 deficiency (Malawi) ─────────────────────────────
    # Malawi `vitb12` is serum B12 in pmol/L. WHO cutoff: <148 pmol/L.
    if ("vitb12" %in% colnames(d)) {
      if (!"b12_def" %in% colnames(d)) {
        d$b12_def <- ifelse(d$vitb12 < 148, 1L, 0L)
        d$b12_def[is.na(d$vitb12)] <- NA_integer_
        n_def <- sum(d$b12_def == 1, na.rm = TRUE)
        n_tot <- sum(!is.na(d$b12_def))
        cat(sprintf("  Derived b12_def from vitb12 < 148 pmol/L: %d/%d deficient\n", n_def, n_tot))
      }
    }

    # ── Zinc: zinc_def already exists in raw data (IZiNCG time-of-day cutoffs)
    # No derivation needed, but log for completeness
    if ("zinc_def" %in% colnames(d) && "zn_gdl" %in% colnames(d)) {
      n_def <- sum(d$zinc_def == 1, na.rm = TRUE)
      n_tot <- sum(!is.na(d$zinc_def))
      cat(sprintf("  Zinc: zinc_def present, %d/%d deficient\n", n_def, n_tot))
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

  # ── Recode Sierra Leone 1/2 → 1/0 for folate and B12 ──────────────────
  # Sierra Leone codes deficiency as 1=deficient, 2=not deficient.
  recode_12 <- function(d, col) {
    if (col %in% colnames(d)) {
      vals <- d[[col]]
      if (is.numeric(vals) && all(vals[!is.na(vals)] %in% c(1, 2))) {
        d[[col]] <- ifelse(vals == 1, 1L, 0L)
        n_def <- sum(d[[col]] == 1, na.rm = TRUE)
        n_tot <- sum(!is.na(d[[col]]))
        cat(sprintf("  Recoded %s from {1,2} to {1,0}: %d/%d deficient\n", col, n_def, n_tot))
      }
    }
    d
  }
  d <- recode_12(d, "gw_wFolDef")
  d <- recode_12(d, "gw_wB12DefWHO")
  d <- recode_12(d, "gw_wB12DefHerb")

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


#' Resolve the DHS survey year a country's merge will use.
#'
#' Mirrors the resolution in merge_external_predictors() exactly (config first,
#' then the legacy hardcoded map) so the stamp below covers precisely the files
#' that will actually be read.
dhs_year_for <- function(cc) {
  y <- cc$dhs_year
  if (!is.null(y) && !is.na(y)) return(as.integer(y))
  fallback <- c("Gambia" = 2019L, "Ghana" = 2014L,
                "Sierra Leone" = 2013L, "Malawi" = 2015L)
  unname(fallback[cc$country])
}

#' DHS clean files that merge_external_predictors() reads for one country.
#'
#' @return character vector of EXISTING paths (possibly empty)
dhs_clean_inputs <- function(cc, dhs_dir = here::here("data", "DHS", "clean")) {
  yr <- dhs_year_for(cc)
  if (is.null(yr) || is.na(yr) || !dir.exists(dhs_dir)) return(character(0))
  p <- file.path(dhs_dir, paste0(cc$country, "_", yr,
                                 c("_dhs_admin2_wide.rds",
                                   "_dhs_custom_admin2_wide.rds")))
  p[file.exists(p)]
}

#' Content stamp for a country's DHS clean inputs.
#'
#' load_dhs_admin2() is called from INSIDE merge_external_predictors(), so
#' {targets} cannot see those .rds files: regenerating them (e.g. re-running
#' src/DHS/DHS_admin2_aggregation.R) would leave merged_ext_* serving a stale
#' cached result forever. The dhs_stamp_* targets in _targets.R depend on this,
#' so a changed file invalidates the merge. Same hazard and same remedy as
#' path_merged_* and gee_parity_stamp_* -- see docs/pipeline_audit_2026-08.md #2.
#'
#' @return named character vector of md5 sums, or NA_character_ when absent
dhs_clean_stamp <- function(cc, dhs_dir = here::here("data", "DHS", "clean")) {
  f <- dhs_clean_inputs(cc, dhs_dir)
  if (!length(f)) return(NA_character_)
  stats::setNames(unname(tools::md5sum(f)), basename(f))
}

#' The external-predictor cache file merge_external_predictors() reads.
#'
#' Path construction mirrors merge_external_predictors() exactly, so the stamp
#' covers precisely the file that will be read.
ext_cache_input <- function(cc, cache_dir = here::here("data", "external_cache")) {
  p <- file.path(cache_dir,
                 paste0(tolower(gsub(" ", "_", cc$country)),
                        "_external_predictors.rds"))
  p[file.exists(p)]
}

#' Content stamp for a country's external-predictor cache.
#'
#' SAME HAZARD AS dhs_clean_stamp(), SAME FUNCTION, AND IT BIT US.
#' merge_external_predictors() reads
#' data/external_cache/<country>_external_predictors.rds directly, so {targets}
#' never sees it. The caches were refreshed on 2026-08-27 (note the
#' .bak_pre_wfp_refresh siblings) but every merged_ext_* target had been built
#' 2026-08-18..23 and was never invalidated -- so the refresh reached NO result.
#' Measured cost when it was found (2026-08-30): Sierra Leone was serving 0
#' map2_ and 0 wfp_ columns while its cache held 57 and 13; rebuilding recovers
#' +196 columns for Sierra Leone, +9 Ghana, +3 Malawi. Sierra Leone's absence
#' also zeroed the shared map2_accessibility_ domain in every pooled/LOCO model,
#' which is what surfaced as the "[pool] contributes 0 pooled predictors"
#' warning.
#'
#' @return named character vector of md5 sums, or NA_character_ when absent
ext_cache_stamp <- function(cc, cache_dir = here::here("data", "external_cache")) {
  f <- ext_cache_input(cc, cache_dir)
  if (!length(f)) return(NA_character_)
  stats::setNames(unname(tools::md5sum(f)), basename(f))
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


#' Identify population-scaled count columns and near-duplicate MAP snapshots to drop
#'
#' Pure on column names (no data needed), so it is leakage-free and cheap.
#'
#' counts -> rates:
#'   - MAP: drop `MAP_<snapshot>_<ind>_Count` when the population-free
#'     `..._Rate` sibling exists (the rate IS the converted version).
#'   - IHME: drop population-scaled counts (`*_count(s)`, `*number_of_*`). IHME
#'     ships rate/prevalence versions of these indicators; arithmetic conversion
#'     is avoided because the matching age/sex/year denominator is not reliably
#'     identifiable from column names. Counts encode population size, which is
#'     epidemiologically uninformative and degrades cross-country transfer.
#'
#' MAP snapshot dedup:
#'   - `MAP_<YYYYMM>_<indicator>` families carry multiple release dates
#'     (e.g. 2022.06 + 2024.06). Keep only the snapshot whose year is nearest the
#'     survey year; drop the rest (they are near-collinear duplicates).
#'
#' @param cols character vector of (proxy) column names
#' @param survey_year integer survey year, for snapshot selection
#' @return character vector of column names to drop
prune_predictor_cols <- function(cols, survey_year = NA_integer_) {
  drop <- character(0)

  # --- counts -> rates ---
  map_cnt <- grep("^MAP_.*_Count$", cols, value = TRUE)
  for (cc_col in map_cnt) {
    rate <- sub("_Count$", "_Rate", cc_col)
    if (rate %in% cols) drop <- c(drop, cc_col)
  }
  ihme_cnt <- grep("^ihme_.*(_counts?$|number_of_)", cols, value = TRUE, ignore.case = TRUE)
  drop <- c(drop, ihme_cnt)

  # --- MAP snapshot dedup: keep the date nearest survey_year per indicator ---
  map <- grep("^MAP_[0-9]{6}_", cols, value = TRUE)
  if (length(map) && !is.na(survey_year)) {
    yr   <- as.integer(substr(sub("^MAP_", "", map), 1, 4))   # YYYY from YYYYMM
    stem <- sub("^MAP_[0-9]{6}_", "", map)                    # indicator stem
    for (s in unique(stem)) {
      idx <- which(stem == s)
      if (length(idx) < 2) next
      keep <- idx[which.min(abs(yr[idx] - survey_year))]
      drop <- c(drop, map[setdiff(idx, keep)])
    }
  }

  # --- GEE covariate hygiene (opt-in, GEE_COVARIATE_HYGIENE=true) ---
  # Cross-band _annual_* summaries over non-commensurable bands, plus static
  # layers' identical per-year copies. See R/gee_band_semantics.R.
  drop <- c(drop, prune_gee_covariates(cols, verbose = FALSE))

  unique(drop)
}


#' Build an outcome-specific dataset (one population x one micronutrient)
#'
#' Filters to the correct population, selects predictors, removes leakage
#' columns, and drops rows with missing outcome.
#'
#' @param merged_data The full merged dataset
#' @param cc Country config (from get_country_configs())
#' @param oc Outcome config (one element of cc$outcomes)
#' @param clean_predictors Drop population-scaled count columns (keeping rate
#'   siblings) and collapse near-duplicate MAP temporal snapshots to the one
#'   nearest the survey year. Defaults to TRUE; disable with
#'   FE_CLEAN_PREDICTORS=false. See prune_predictor_cols().
#' @return list with components: data, Xvars_full, Xvars, Xvars_bundle, domain_vars
build_outcome_dataset <- function(merged_data, cc, oc,
                                  clean_predictors =
                                    !identical(tolower(Sys.getenv("FE_CLEAN_PREDICTORS", "true")), "false")) {

  d <- merged_data

  # Filter to population. The !is.na() guard matters: a bare `d[[col]] == val`
  # comparison yields NA for missing flags, and `d[NA, ]` injects all-NA phantom
  # rows rather than dropping them. Matches cluster_aggregation.R's filter.
  pop_col <- cc$child_flag
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    d <- d[!is.na(d[[pop_col]]) & d[[pop_col]] == oc$child_flag_val, , drop = FALSE]
  }

  # 2026-06-24: Gambia biomarker outcomes — the configured `gw_svy_weight` is
  # NA/zero across the blood (biomarker) subsample, so it mis-weights every
  # assayed estimate. Swap in the population-specific blood-draw weight
  # (`gw_c_blood_weight` for children, `gw_w_blood_weight` for women), which fully
  # covers the assayed sample. d is already population-filtered here, so all rows
  # share one population. All downstream consumers read cc$weight_col unchanged.
  # See docs/survey_weighting_critique.md.
  if (identical(cc$country, "Gambia") && !is.null(cc$weight_col) &&
      cc$weight_col %in% colnames(d)) {
    bw <- if (!is.null(oc$population) && oc$population == "children")
            "gw_c_blood_weight" else "gw_w_blood_weight"
    if (bw %in% colnames(d)) {
      bw_vals <- suppressWarnings(as.numeric(d[[bw]]))
      old <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
      n_bad <- sum(is.na(old) | old <= 0)
      if (sum(is.finite(bw_vals) & bw_vals > 0) >= 0.9 * nrow(d)) {
        d[[cc$weight_col]] <- bw_vals
        cat(sprintf("  [weights] Gambia %s: using blood-subsample weight %s (replaced %d NA/zero of %d svy weights)\n",
                    oc$tag, bw, n_bad, nrow(d)))
      }
    }
  }

  # 2026-06-24 (DC-H2 #1): make the uniform BRINDA-adjusted VAD the VitA binary
  # here (no-op for non-VitA), so national + within-country + area-transport all
  # use one definition. Done before the non-missing filter so it operates on the
  # BRINDA binary. CRP/AGP are still present (pre keep_cols). See
  # docs/dc_h2_brinda_validation.md.
  d <- apply_brinda_vita_binary(d, cc, oc, "[build]")

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
    if (nm == "GW") {
      # (a) per-country configured patterns
      cfg_pat <- cc$gw_exclude_patterns
      leak <- if (!is.null(cfg_pat) && length(cfg_pat))
        Reduce("|", lapply(cfg_pat, function(p) grepl(p, prefix_cols, ignore.case = TRUE)))
      else rep(FALSE, length(prefix_cols))
      # (b) hardcoded outcome-biomarker guard. 2026-06-23 fix for outcome leakage:
      # the per-country lists missed one-"r" "FerAdj" (e.g. gw_LNwFerAdjBR1),
      # "VitADef"/"VitAInsuff" (pattern was "VAD"), and hemoglobin. These are
      # outcome-DERIVED biomarkers, not legitimate diet/behaviour covariates
      # (gw_wVitASuppl, gw_wVitARichFood, gw_wIodSalt are deliberately kept).
      guard_ci <- paste(c("FerAdj", "RBPAdj", "Retinol",
                          "VitADef", "VitAInsuff", "VitADefic",
                          "FeDef", "FolDef", "B12Def", "ZincDef"),
                        collapse = "|")
      leak <- leak |
        grepl(guard_ci, prefix_cols, ignore.case = TRUE) |
        grepl("Hb", prefix_cols)   # case-SENSITIVE: matches gw_wHb/gw_cHb/gw_HbCat,
                                    # not lower-case household vars (gw_hBuy*, gw_hBirds*)
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

  # Predictor hygiene: counts->rates + dedup near-duplicate MAP snapshots.
  if (isTRUE(clean_predictors)) {
    drop_hy <- prune_predictor_cols(Xvars_no_gw, survey_year = cc$survey_year)
    if (length(drop_hy) > 0) {
      cat(sprintf("  [hygiene] dropping %d predictor cols (population-scaled counts + duplicate MAP snapshots)\n",
                  length(drop_hy)))
      Xvars_no_gw <- setdiff(Xvars_no_gw, drop_hy)
      Xvars_full  <- setdiff(Xvars_full,  drop_hy)
    }
  }

  # Outcome-specific, biology-driven bundle (proxy-only). NULL prefixes -> use
  # the full proxy set. Always computed so callers can opt in via config
  # (params$use_outcome_bundles); see bundle_prefixes_for_outcome() in config.R.
  bundle_pfx <- tryCatch(bundle_prefixes_for_outcome(oc$tag), error = function(e) NULL)
  Xvars_bundle <- if (is.null(bundle_pfx)) Xvars_no_gw else {
    keep <- Reduce("|", lapply(bundle_pfx, function(p) startsWith(Xvars_no_gw, p)))
    sub  <- Xvars_no_gw[keep]
    if (length(sub) >= 2) sub else Xvars_no_gw  # guard: never empty the model
  }

  # Retain raw RBP/CRP/AGP so compute_svy_admin2() can build the uniform BRINDA
  # VitA transport outcome (DC-H2). These go into $data only, NOT into the
  # predictor set (Xvars_full), so they cannot leak into models. Scoped to VitA
  # outcomes so non-VitA outcome_data hashes (and their downstream) are unchanged.
  brinda_keep <- if (!is.null(oc$tag) && grepl("vitA", oc$tag, ignore.case = TRUE))
    tryCatch(unlist(brinda_rbp_cols(cc$country), use.names = FALSE), error = function(e) NULL)
  else NULL

  # Keep only columns we need
  keep_cols <- unique(c(
    oc$continuous, oc$binary,
    cc$cluster_id, cc$admin1_col, cc$admin2_col,
    cc$weight_col, cc$strata_col, cc$psu_col,
    Xvars_full, brinda_keep
  ))
  keep_cols <- intersect(keep_cols, colnames(d))
  d <- d[, keep_cols, drop = FALSE]

  list(
    data         = d,
    Xvars_full   = Xvars_full,
    Xvars        = Xvars_no_gw,
    Xvars_bundle = Xvars_bundle,
    domain_vars  = domain_vars
  )
}
