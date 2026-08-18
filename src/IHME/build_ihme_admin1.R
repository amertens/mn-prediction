# =============================================================================
# src/IHME/build_ihme_admin1.R
#
# Build per-country admin1-level IHME predictor tables by reading every
# `*_ADMIN_1_*` / `*_AD1_*` / `*_ADM1_*` CSV under data/IHME/<indicator>/, filtering
# to the requested country + year, and pivoting wide on adm1_name. Output is
# written to data/IHME/<country>_<year>_merged_IHME_admin1_data.csv with all
# columns prefixed `ihme_adm1_` so they are picked up by the IHME domain in the
# proxy predictor set (R/config.R domains$IHME prefix = "ihme_").
#
# Deduplication rules (same logic the admin2 clean scripts now use):
#   - HIV 2017 admin1 files are not shipped under data/IHME/HIV — only the 2018
#     SSA release is, so no extra exclusion is required.
#   - DBM wasting is dropped because CGF wasting covers it (CGF is the focused
#     dataset); DBM overweight is kept (it is not in any other indicator).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(here)
})

#' Build admin1 IHME wide table for one country.
#'
#' @param country_name  IHME `ADM0_NAME` value (e.g. "Ghana"). Pass a vector for
#'   countries with multiple spellings (e.g. c("Gambia", "The Gambia")).
#' @param year_filter   Year to extract (single integer).
#' @param out_path      Output CSV path. If NULL, returns the data frame.
#' @return invisibly returns the wide data frame.
build_ihme_admin1 <- function(country_name, year_filter, out_path = NULL) {

  ihme_root <- here("data", "IHME")

  # Discover all admin1 CSVs (handles ADMIN_1 / AD1 / ADM1 naming variations)
  all_csv <- list.files(ihme_root, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
  is_adm1 <- grepl("ADMIN_?1|_AD1[^0-9]|_ADM1[^0-9]", basename(all_csv),
                   ignore.case = TRUE)
  files <- all_csv[is_adm1]

  # Drop codebooks, input-source manifests, disputed-territory files
  files <- files[!grepl("codebook|disputed|sources|info_sheet",
                        basename(files), ignore.case = TRUE)]

  # DBM wasting duplicates CGF wasting — keep DBM overweight only.
  files <- files[!grepl("DBM.*WASTING", basename(files), ignore.case = TRUE)]

  cat(sprintf("[ihme-admin1] %s %d: scanning %d files\n",
              paste(country_name, collapse = "/"), year_filter, length(files)))

  long_list <- list()

  for (f in files) {
    df <- tryCatch(
      readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                      name_repair = "minimal"),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) next
    names(df) <- tolower(names(df))

    # Need adm0_name, adm1_name, year, and a numeric value column
    if (!all(c("adm0_name", "adm1_name", "year") %in% names(df))) next
    value_col <- if ("mean"  %in% names(df)) "mean"
                 else if ("value" %in% names(df)) "value"
                 else NA_character_
    if (is.na(value_col)) next

    df_f <- df[df$adm0_name %in% country_name & df$year == year_filter, ,
               drop = FALSE]
    if (nrow(df_f) == 0) next

    # File-level tag (strips date stamps, ADMIN markers, IHME_ prefix).
    file_tag <- basename(f)
    file_tag <- sub("\\.csv$", "", file_tag, ignore.case = TRUE)
    file_tag <- gsub("^IHME_", "", file_tag, ignore.case = TRUE)
    file_tag <- gsub("_Y\\d{4}M\\d{2}D\\d{2}.*$", "", file_tag,
                     ignore.case = TRUE)
    file_tag <- gsub("_ADMIN_?1|_AD1|_ADM1", "", file_tag, ignore.case = TRUE)
    file_tag <- gsub("\\s+", "_", file_tag)

    # Build a stratifier from whatever optional cols exist.
    strat_cols <- intersect(
      c("age_group_name", "age_group", "age_group_id_label",
        "sex", "sex_id_label", "measure", "metric"),
      names(df_f)
    )

    if (length(strat_cols) > 0) {
      strat <- apply(df_f[, strat_cols, drop = FALSE], 1L,
                     function(r) paste(r, collapse = "_"))
    } else {
      strat <- ""
    }
    variable <- paste(file_tag, strat, sep = "__")

    long_list[[length(long_list) + 1L]] <- data.frame(
      adm1_name = df_f$adm1_name,
      variable  = variable,
      value     = suppressWarnings(as.numeric(df_f[[value_col]])),
      stringsAsFactors = FALSE
    )
  }

  if (length(long_list) == 0) {
    warning(sprintf("No admin1 IHME rows found for %s %d.",
                    paste(country_name, collapse = "/"), year_filter))
    return(invisible(NULL))
  }

  long <- dplyr::bind_rows(long_list) |>
    dplyr::filter(!is.na(value), !is.na(adm1_name), adm1_name != "")

  # Collapse any duplicate (adm1, variable) keys by mean.
  long <- long |>
    dplyr::group_by(adm1_name, variable) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  wide <- tidyr::pivot_wider(long, id_cols = adm1_name,
                             names_from = variable, values_from = value) |>
    janitor::clean_names()

  # Prefix every column with ihme_adm1_ (including the id column → ihme_adm1_name)
  names(wide) <- ifelse(names(wide) == "adm1_name",
                        "ihme_adm1_name",
                        paste0("ihme_adm1_", names(wide)))

  cat(sprintf("[ihme-admin1] %s %d: %d admin1 units x %d variables\n",
              paste(country_name, collapse = "/"), year_filter,
              nrow(wide), ncol(wide) - 1L))

  if (!is.null(out_path)) {
    readr::write_csv(wide, out_path)
    cat(sprintf("[ihme-admin1] wrote %s\n", out_path))
  }
  invisible(wide)
}


#' Merge admin1 IHME wide table into a survey data frame.
#'
#' Joins on the survey's Admin1 column (or a user-supplied column) using a light
#' fuzzy match (case-insensitive, stripped of administrative suffixes) so that
#' minor GADM-vs-IHME name discrepancies don't drop rows.
#'
#' @param df              Survey data frame (must contain `admin1_col`).
#' @param admin1_csv      Path to the wide admin1 IHME CSV written by
#'   build_ihme_admin1().
#' @param admin1_col      Name of the admin1 column in `df` (default "Admin1").
#' @param min_similarity  Jaro-Winkler similarity threshold for fuzzy match.
#' @return df with new `ihme_adm1_*` columns left-joined.
merge_ihme_admin1 <- function(df, admin1_csv, admin1_col = "Admin1",
                              min_similarity = 0.85) {
  if (!file.exists(admin1_csv)) {
    warning(sprintf("Admin1 IHME CSV not found: %s — skipping merge.",
                    admin1_csv))
    return(df)
  }
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("Package 'stringdist' is required for fuzzy admin1 matching.")
  }
  ihme1 <- readr::read_csv(admin1_csv, show_col_types = FALSE, progress = FALSE)

  clean <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("\\s+region$",   "", x)
    x <- gsub("\\s+province$", "", x)
    x <- gsub("\\s+district$", "", x)
    x <- gsub("[/.,'-]",      " ", x)
    x <- gsub("\\s+",         " ", x)
    trimws(x)
  }

  df_a1     <- unique(df[[admin1_col]])
  ihme_a1   <- ihme1$ihme_adm1_name

  sim <- 1 - stringdist::stringdistmatrix(clean(df_a1), clean(ihme_a1),
                                          method = "jw")
  # Robust to all-NA rows (e.g. an NA admin1 from GPS-unmatched clusters):
  # max.col always returns a clean integer vector, unlike apply(which.max)
  # which yields a list when a row is empty.
  sim[is.na(sim)] <- -Inf
  best_idx <- max.col(sim, ties.method = "first")
  best_sim <- sim[cbind(seq_len(nrow(sim)), best_idx)]
  best_sim[!is.finite(best_sim)] <- NA_real_
  lookup <- data.frame(
    .a1_key       = df_a1,
    ihme_adm1_name = ifelse(!is.na(best_sim) & best_sim >= min_similarity,
                             ihme_a1[best_idx], NA_character_),
    a1_similarity = best_sim,
    stringsAsFactors = FALSE
  )
  names(lookup)[1L] <- admin1_col

  n_matched <- sum(!is.na(lookup$ihme_adm1_name))
  cat(sprintf("[ihme-admin1] Admin1 matched %d/%d units (min sim = %.2f)\n",
              n_matched, nrow(lookup), min_similarity))

  df <- df |>
    dplyr::left_join(lookup[, c(admin1_col, "ihme_adm1_name")],
                     by = admin1_col) |>
    dplyr::left_join(ihme1, by = "ihme_adm1_name")

  df
}

