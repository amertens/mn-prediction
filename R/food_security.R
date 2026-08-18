# =============================================================================
# R/food_security.R
#
# Merge food security indicators from HFID and Cadre Harmonisé into the
# per-country merged datasets. Creates "fsec_" prefixed predictors.
#
# HFID: Harmonized Food Insecurity Dataset (Admin2-level, monthly, 2007-2024)
#   - FCS (food consumption score) prevalence
#   - rCSI (reduced coping strategy index) prevalence
#   - IPC/CH phase classifications
#
# Cadre Harmonisé: IPC-compatible food security classification for West/Central
# Africa (Admin1 or Admin2, biannual exercises, 2014-2025)
#   - IPC phase population counts (phase 1-5, phase 3+5)
#   - Component scores: FCG, HDDS, HHS, LhHCSCat, rCSI
# =============================================================================

# Survey year for each country (used to match food security data temporally)
# Fieldwork dates:
#   Gambia GMNS:         13 Mar – 4 May 2018
#   Ghana GMS:           27 Apr – 9 Jun 2017
#   Sierra Leone SLMS:   11 Nov – 2 Dec 2013
#   Malawi MNS:          8 Dec 2015 – 16 Feb 2016
# These are fallbacks — prefer cc$survey_year from R/config.R when available.
COUNTRY_SURVEY_YEARS <- list(
  Gambia         = 2018L,
  Ghana          = 2017L,
  SierraLeone    = 2013L,
  Malawi         = 2016L,
  Tanzania       = 2010L   # TDHS 2009-10
)

# Mapping from our config country names to HFID ADMIN0 names
HFID_COUNTRY_MAP <- c(
  Gambia      = "Gambia",
  Ghana       = "Ghana",
  SierraLeone = "Sierra Leone",
  Malawi      = "Malawi",
  # HFID Tanzania coverage starts 2011 (no 2010) — load_hfid_country's
  # nearest-year fallback uses 2011, a 1-year gap from the 2010 survey.
  Tanzania    = "Tanzania"
)

# Mapping from config country names to CH adm0_name
CH_COUNTRY_MAP <- c(
  Gambia      = "Gambia",
  Ghana       = "Ghana",
  SierraLeone = "Sierra Leone"
  # Malawi and Tanzania not in Cadre Harmonisé (West/Central Africa only)
)


#' Load and prepare HFID data for one country
#'
#' Reads the HFID CSV, filters to the country and a time window around the
#' survey year, and aggregates monthly data to Admin2 annual averages.
#'
#' @param hfid_path Path to hfid_hv1.csv
#' @param country_key Config key (e.g. "Gambia", "SierraLeone")
#' @param window_months Months around survey midpoint to average (default 12)
#' @return data.frame with Admin1, Admin2, and fsec_hfid_* columns
load_hfid_country <- function(hfid_path, country_key, window_months = 12L) {

  hfid_country_name <- HFID_COUNTRY_MAP[[country_key]]
  survey_year <- COUNTRY_SURVEY_YEARS[[country_key]]
  if (is.null(hfid_country_name) || is.null(survey_year)) return(NULL)

  hfid <- data.table::fread(hfid_path, showProgress = FALSE)

  # Filter to country
  hfid <- hfid[hfid$ADMIN0 == hfid_country_name, ]
  if (nrow(hfid) == 0) return(NULL)

  # Parse year_month to Date
  hfid$date <- as.Date(paste0(hfid$year_month, "-01"))
  hfid$year <- as.integer(format(hfid$date, "%Y"))

  # Time window: use survey year and the preceding window_months
  # (prefer data at or before survey, not after)
  end   <- as.Date(paste0(survey_year, "-12-31"))
  start <- end - window_months * 30
  hfid <- hfid[hfid$date >= start & hfid$date <= end, ]

  if (nrow(hfid) == 0) {
    # Fallback: use closest year at or before survey; then nearest overall
    all_hfid <- data.table::fread(hfid_path, showProgress = FALSE)
    all_hfid <- all_hfid[all_hfid$ADMIN0 == hfid_country_name, ]
    all_hfid$year <- as.integer(substr(all_hfid$year_month, 1, 4))
    available_years <- sort(unique(all_hfid$year))
    prior_years <- available_years[available_years <= survey_year]
    if (length(prior_years) > 0) {
      nearest <- max(prior_years)
    } else {
      nearest <- available_years[which.min(abs(available_years - survey_year))]
    }
    cat(sprintf("  [hfid] %s: no data in survey window, using year %d\n",
                country_key, nearest))
    hfid <- all_hfid[all_hfid$year == nearest, ]
  }

  cat(sprintf("  [hfid] %s: %d rows in time window around %d\n",
              country_key, nrow(hfid), survey_year))

  # Clean column names for aggregation
  # Food security indicators to aggregate
  fs_cols <- c("ipc_phase_fews", "ipc_phase_ipcch",
               "fcs_lit", "rcsi_lit",
               "fcs_rt mean", "fcs_rt max", "fcs_rt min",
               "rcsi_rt mean", "rcsi_rt max", "rcsi_rt min")
  fs_cols <- intersect(fs_cols, colnames(hfid))

  # Aggregate to Admin2 means (across months in the window)
  agg <- hfid %>%
    dplyr::group_by(ADMIN1, ADMIN2) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(fs_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  # Rename with fsec_hfid_ prefix
  new_names <- gsub(" ", "_", fs_cols)
  new_names <- paste0("fsec_hfid_", new_names)
  colnames(agg)[match(fs_cols, colnames(agg))] <- new_names

  colnames(agg)[colnames(agg) == "ADMIN1"] <- "hfid_admin1"
  colnames(agg)[colnames(agg) == "ADMIN2"] <- "hfid_admin2"

  agg
}


#' Load and prepare Cadre Harmonisé data for one country
#'
#' Reads the CH Excel file, filters to the country and the exercise year
#' closest to the survey, and computes IPC phase proportions and component
#' indicators at the finest available admin level.
#'
#' @param ch_path Path to cadre_harmonise_caf_ipc_dec25.xlsx
#' @param country_key Config key
#' @return data.frame with admin columns and fsec_ch_* columns, or NULL
load_ch_country <- function(ch_path, country_key) {

  ch_country_name <- CH_COUNTRY_MAP[[country_key]]
  survey_year <- COUNTRY_SURVEY_YEARS[[country_key]]
  if (is.null(ch_country_name) || is.null(survey_year)) return(NULL)

  ch <- readxl::read_xlsx(ch_path)
  ch <- ch[ch$adm0_name == ch_country_name, ]
  if (nrow(ch) == 0) return(NULL)

  # Use "current" analysis only
  ch <- ch[ch$chtype == "current", ]

  # Find closest year at or before the survey; fall back to nearest after
  ch$exercise_year <- as.integer(ch$exercise_year)
  available_years <- sort(unique(ch$exercise_year))
  prior_years <- available_years[available_years <= survey_year]
  if (length(prior_years) > 0) {
    nearest_year <- max(prior_years)
  } else {
    nearest_year <- available_years[which.min(abs(available_years - survey_year))]
    cat(sprintf("  [ch] %s: no data at or before survey year %d, using %d\n",
                country_key, survey_year, nearest_year))
  }
  ch <- ch[ch$exercise_year == nearest_year, ]

  cat(sprintf("  [ch] %s: %d rows from exercise year %d (survey: %d)\n",
              country_key, nrow(ch), nearest_year, survey_year))

  # Determine admin level: use adm2_name if available, else adm1_name
  has_adm2 <- any(!is.na(ch$adm2_name) & ch$adm2_name != "")

  # Compute IPC phase proportions from population counts
  phase_cols <- c("phase1", "phase2", "phase3", "phase4", "phase5")
  phase_cols <- intersect(phase_cols, colnames(ch))

  for (pc in phase_cols) ch[[pc]] <- as.numeric(ch[[pc]])
  ch$population <- as.numeric(ch$population)

  ch$total_classified <- rowSums(ch[, phase_cols, drop = FALSE], na.rm = TRUE)
  for (pc in phase_cols) {
    prop_name <- paste0("prop_", pc)
    ch[[prop_name]] <- ifelse(ch$total_classified > 0,
                              ch[[pc]] / ch$total_classified, NA_real_)
  }
  # Compute combined phase 3+4+5 proportion (food crisis + emergency + famine)
  # phase35 may not exist as a pre-computed column in all CH data versions
  if ("phase35" %in% colnames(ch)) {
    ch$prop_phase35 <- ifelse(ch$total_classified > 0,
                              as.numeric(ch$phase35) / ch$total_classified, NA_real_)
  } else {
    # Sum phases 3-5 manually
    p345_cols <- intersect(c("phase3", "phase4", "phase5"), colnames(ch))
    if (length(p345_cols) > 0) {
      ch$phase35_sum <- rowSums(ch[, p345_cols, drop = FALSE], na.rm = TRUE)
      ch$prop_phase35 <- ifelse(ch$total_classified > 0,
                                ch$phase35_sum / ch$total_classified, NA_real_)
    } else {
      ch$prop_phase35 <- NA_real_
    }
  }

  # Component indicators (already present as phase classifications)
  component_cols <- c("phase_class",
                      "FCG_Poor", "FCG_Borderline", "FCG_Acceptable",
                      "rCSI_Phase1", "rCSI_Phase2", "rCSI_Phase3",
                      "HHS_Phase1", "HHS_Phase2", "HHS_Phase3", "HHS_Phase4", "HHS_Phase5",
                      "HDDS_Phase1", "HDDS_Phase2", "HDDS_Phase3", "HDDS_Phase4", "HDDS_Phase5")
  component_cols <- intersect(component_cols, colnames(ch))

  # Convert to numeric
  for (cc_col in component_cols) {
    ch[[cc_col]] <- suppressWarnings(as.numeric(ch[[cc_col]]))
  }

  # Select columns to keep
  keep_cols <- c("prop_phase1", "prop_phase2", "prop_phase3",
                 "prop_phase4", "prop_phase5", "prop_phase35",
                 component_cols)
  keep_cols <- intersect(keep_cols, colnames(ch))

  # Aggregate if multiple exercises/rows per area
  if (has_adm2) {
    grp <- c("adm1_name", "adm2_name")
  } else {
    grp <- "adm1_name"
  }

  agg <- ch %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(keep_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  # Rename with fsec_ch_ prefix
  for (col in keep_cols) {
    new_name <- paste0("fsec_ch_", col)
    colnames(agg)[colnames(agg) == col] <- new_name
  }

  # Rename admin columns for merging
  colnames(agg)[colnames(agg) == "adm1_name"] <- "ch_admin1"
  if ("adm2_name" %in% colnames(agg))
    colnames(agg)[colnames(agg) == "adm2_name"] <- "ch_admin2"

  agg
}


#' Merge food security data into a country's merged dataset
#'
#' Performs fuzzy name matching on Admin1/Admin2 names between the food security
#' sources and the MN survey data. Falls back to Admin1-level merge when Admin2
#' matching is poor.
#'
#' @param merged_data The country's merged dataset (from load_merged_data)
#' @param cc Country config
#' @param hfid_path Path to HFID CSV
#' @param ch_path Path to CH Excel
#' @return merged_data with additional fsec_* columns
merge_food_security <- function(merged_data, cc, hfid_path, ch_path) {

  country_key <- gsub(" ", "", cc$country)  # "Sierra Leone" -> "SierraLeone"
  cat(sprintf("[fsec] Merging food security data for %s\n", cc$country))

  n_before <- ncol(merged_data)

  # ── HFID merge ──────────────────────────────────────────────────────────
  hfid <- tryCatch(load_hfid_country(hfid_path, country_key),
                   error = function(e) { warning(e$message); NULL })

  if (!is.null(hfid) && nrow(hfid) > 0) {
    # Try Admin2 merge first
    merged_data <- merge_by_admin(merged_data, hfid,
                                  admin1_col = cc$admin1_col,
                                  admin2_col = cc$admin2_col,
                                  src_admin1 = "hfid_admin1",
                                  src_admin2 = "hfid_admin2",
                                  source_name = "HFID")
  }

  # ── Cadre Harmonisé merge ───────────────────────────────────────────────
  ch <- tryCatch(load_ch_country(ch_path, country_key),
                 error = function(e) { warning(e$message); NULL })

  if (!is.null(ch) && nrow(ch) > 0) {
    has_ch_admin2 <- "ch_admin2" %in% colnames(ch)
    merged_data <- merge_by_admin(merged_data, ch,
                                  admin1_col = cc$admin1_col,
                                  admin2_col = if (has_ch_admin2) cc$admin2_col else NULL,
                                  src_admin1 = "ch_admin1",
                                  src_admin2 = if (has_ch_admin2) "ch_admin2" else NULL,
                                  source_name = "CH")
  }

  # Drop fsec_ columns that came out entirely NA — indicators absent from this
  # country/year's source coverage (e.g. HFID's early years carry only FCS, not
  # IPC/rCSI or the rolling-window stats). An all-NA column carries no signal and
  # only adds imputation burden, so keep the merged dataset honest.
  fsec_cols <- grep("^fsec_", colnames(merged_data), value = TRUE)
  if (length(fsec_cols) > 0) {
    all_na <- fsec_cols[vapply(merged_data[fsec_cols],
                               function(x) all(is.na(x)), logical(1))]
    if (length(all_na) > 0) {
      merged_data <- merged_data[, setdiff(colnames(merged_data), all_na), drop = FALSE]
      cat(sprintf("  [fsec] Dropped %d all-NA fsec_ column(s): %s\n",
                  length(all_na),
                  paste(sub("^fsec_(hfid|ch)_", "", all_na), collapse = ", ")))
    }
  }

  n_after <- ncol(merged_data)
  n_fsec <- sum(startsWith(colnames(merged_data), "fsec_"))
  cat(sprintf("  [fsec] Added %d food security columns (%d total fsec_ cols)\n",
              n_after - n_before, n_fsec))

  merged_data
}


#' Fuzzy Admin name matching merge helper
#'
#' Tries exact match on Admin2, then falls back to fuzzy string matching,
#' then falls back to Admin1-level merge.
#'
#' @param target Target data.frame (merged MN data)
#' @param source Source data.frame (food security data)
#' @param admin1_col Admin1 column name in target
#' @param admin2_col Admin2 column name in target (or NULL for admin1-only)
#' @param src_admin1 Admin1 column name in source
#' @param src_admin2 Admin2 column name in source (or NULL)
#' @param source_name Label for logging
#' @return target with source columns merged in
merge_by_admin <- function(target, source, admin1_col, admin2_col,
                           src_admin1, src_admin2, source_name) {

  # Get the data columns (not admin columns)
  data_cols <- setdiff(colnames(source), c(src_admin1, src_admin2))

  # ── Try Admin2 merge ────────────────────────────────────────────────────
  if (!is.null(admin2_col) && !is.null(src_admin2) && src_admin2 %in% colnames(source)) {

    target_a2 <- unique(target[[admin2_col]])
    target_a2 <- target_a2[!is.na(target_a2)]
    source_a2 <- unique(source[[src_admin2]])

    # Build fuzzy match lookup
    match_lookup <- build_name_lookup(target_a2, source_a2)
    match_rate <- mean(!is.na(match_lookup))

    if (match_rate >= 0.3) {
      # Good enough Admin2 match — use it
      source$merge_key <- match_lookup[source[[src_admin2]]]
      # For unmatched source rows, try Admin1 fallback
      source_matched <- source[!is.na(source$merge_key), ]

      if (nrow(source_matched) > 0) {
        merge_df <- source_matched[, c("merge_key", data_cols), drop = FALSE]
        # Average if multiple source rows map to same target Admin2
        merge_df <- merge_df %>%
          dplyr::group_by(merge_key) %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)),
                           .groups = "drop") %>%
          as.data.frame()

        colnames(merge_df)[1] <- admin2_col
        target <- dplyr::left_join(target, merge_df, by = admin2_col)
        # Count rows that got ANY non-NA indicator, not just data_cols[1] — the
        # first data column can itself be all-NA (e.g. HFID ipc_phase_fews),
        # which would falsely report "0 rows got data" even when others merged.
        got_any <- rowSums(!is.na(as.data.frame(target)[, data_cols, drop = FALSE])) > 0
        n_matched <- sum(got_any)
        cat(sprintf("  [%s] Admin2 merge: %.0f%% name match, %d/%d rows got data\n",
                    source_name, match_rate * 100, n_matched, nrow(target)))
        return(target)
      }
    }
  }

  # ── Fallback: Admin1 merge ──────────────────────────────────────────────
  if (!is.null(src_admin1) && src_admin1 %in% colnames(source)) {

    # Aggregate source to Admin1
    source_a1 <- source %>%
      dplyr::group_by(.data[[src_admin1]]) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(data_cols),
                                     ~ mean(.x, na.rm = TRUE)),
                       .groups = "drop") %>%
      as.data.frame()

    target_a1 <- unique(target[[admin1_col]])
    target_a1 <- target_a1[!is.na(target_a1)]
    source_a1_names <- source_a1[[src_admin1]]

    match_lookup <- build_name_lookup(target_a1, source_a1_names)

    source_a1$merge_key <- match_lookup[source_a1[[src_admin1]]]
    source_a1 <- source_a1[!is.na(source_a1$merge_key), ]

    if (nrow(source_a1) > 0) {
      merge_df <- source_a1[, c("merge_key", data_cols), drop = FALSE]
      colnames(merge_df)[1] <- admin1_col
      target <- dplyr::left_join(target, merge_df, by = admin1_col)
      got_any <- rowSums(!is.na(as.data.frame(target)[, data_cols, drop = FALSE])) > 0
      n_matched <- sum(got_any)
      cat(sprintf("  [%s] Admin1 merge: %d/%d target areas matched, %d/%d rows got data\n",
                  source_name, nrow(source_a1), length(target_a1),
                  n_matched, nrow(target)))
    } else {
      cat(sprintf("  [%s] No admin names matched — skipping\n", source_name))
    }
  }

  target
}


#' Build a name matching lookup from target names to source names
#'
#' Uses exact match first, then case-insensitive, then fuzzy string distance.
#'
#' @param target_names Character vector of target admin names
#' @param source_names Character vector of source admin names
#' @return Named character vector: source_name -> matched target_name (or NA)
build_name_lookup <- function(target_names, source_names) {

  target_clean <- tolower(trimws(target_names))
  source_clean <- tolower(trimws(source_names))

  # Result: for each source name, which target name it maps to
  lookup <- stats::setNames(rep(NA_character_, length(source_names)), source_names)

  # Exact match (case insensitive)
  for (i in seq_along(source_names)) {
    exact <- which(target_clean == source_clean[i])
    if (length(exact) == 1) {
      lookup[source_names[i]] <- target_names[exact]
    }
  }

  # Fuzzy match for unmatched (using adist)
  unmatched <- is.na(lookup)
  if (any(unmatched) && length(target_names) > 0) {
    um_source <- source_names[unmatched]
    um_clean  <- source_clean[unmatched]

    # Use full string matching (not partial) to avoid false positives
    # like "Bo" matching "Bong County" or "Bole" matching "Bolgatanga".
    # Tighter threshold (0.3 instead of 0.4) reduces false matches.
    dist_mat <- utils::adist(um_clean, target_clean, ignore.case = TRUE, partial = FALSE)
    # Normalize by max string length
    len_mat <- outer(nchar(um_clean), nchar(target_clean), pmax)
    norm_dist <- dist_mat / pmax(len_mat, 1)

    for (i in seq_along(um_source)) {
      best <- which.min(norm_dist[i, ])
      if (length(best) == 1 && norm_dist[i, best] < 0.3) {
        lookup[um_source[i]] <- target_names[best]
      }
    }
  }

  lookup
}
