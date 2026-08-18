# =============================================================================
# src/IHME/build_ihme_admin2.R
#
# Self-contained R replacement for the original Stata rollup (data/IHME/
# data_merge.do) + SL_IHME_clean.R, so a NEW country's admin-2 IHME predictor
# table can be built without Stata and without a hand-made per-country rollup.
#
# It reads every `*_ADMIN_2_*` / `*_ADM2_*` CSV under data/IHME/<family>/,
# filters to the requested country + year, harmonises the (slightly varying)
# per-family column layouts, applies the same dedup + measure-rename rules the
# original pipeline used, and pivots wide on
#   (age_group_name, age_group_id, sex, measure, metric)
# with every column prefixed `ihme_` — matching the column naming of the four
# existing countries' *_merged_IHME_data.csv files so the pooled / LOCO
# common-vocabulary set lines up.
#
# IMPORTANT COVERAGE NOTE: the repo ships admin-2 CSVs for 11 of the 13 IHME
# families. **Malaria and education have NO admin-2 CSVs here** (the existing
# countries' malaria/education ihme_ columns came from the original Box drive).
# So a country built by this script will lack `ihme_*malaria*` and
# `ihme_*education*` admin-2 columns. Malaria is largely covered by the MAP_
# domain (Malaria Atlas) in the merge step, and education by DHS proxies, so the
# practical loss is small — but it is a real difference from Ghana/SL/Malawi/
# Gambia and is reported by the parity diagnostic below.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(here)
})

#' Build admin-2 IHME wide table for one country.
#'
#' @param country_name IHME `ADM0_NAME` value(s) (e.g. "Tanzania"; pass a vector
#'   for multiple spellings, e.g. c("Gambia", "The Gambia")).
#' @param year_filter  Survey year to extract (single integer). IHME geospatial
#'   files run through 2017/2018/2019 depending on family; rows are filtered to
#'   exactly this year (use the survey year, matching the existing countries).
#' @param out_path     Output CSV path; if NULL, returns the data frame only.
#' @param reference_csv Optional path to an existing country's
#'   *_merged_IHME_data.csv (e.g. Ghana) — if given, prints a column-name parity
#'   report so you can confirm cross-country alignment.
#' @return invisibly, the wide data frame.
build_ihme_admin2 <- function(country_name, year_filter, out_path = NULL,
                              reference_csv = NULL) {

  ihme_root <- here("data", "IHME")

  all_csv <- list.files(ihme_root, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
  is_adm2 <- grepl("ADMIN_?2|_AD2[^0-9]|_ADM2[^0-9]", basename(all_csv),
                   ignore.case = TRUE)
  files <- all_csv[is_adm2]
  files <- files[!grepl("codebook|disputed|sources|info_sheet",
                        basename(files), ignore.case = TRUE)]
  # DBM wasting duplicates CGF wasting — keep DBM overweight only.
  files <- files[!grepl("DBM.*WASTING", basename(files), ignore.case = TRUE)]

  # Family = the data/IHME/<family>/ folder a file lives under. Match against the
  # known family names by path substring (robust to OS path separators).
  known_families <- c("Anemia", "CGF", "DBM", "EBF", "HIV", "LF", "Malaria",
                      "mcv1", "oral rehydration", "u5 diarrhea", "u5 mort",
                      "WASH access", "education")
  family_of <- function(path) {
    norm <- gsub("\\\\", "/", path)
    hit <- known_families[vapply(known_families,
                                 function(fam) grepl(paste0("/", fam, "/"), norm,
                                                     fixed = TRUE),
                                 logical(1))]
    if (length(hit)) hit[1] else NA_character_
  }

  cat(sprintf("[ihme-admin2] %s %d: scanning %d admin-2 files\n",
              paste(country_name, collapse = "/"), year_filter, length(files)))

  long_list <- list()
  fams_seen <- character(0)

  for (f in files) {
    df <- tryCatch(
      readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                      name_repair = "minimal"),
      error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    names(df) <- tolower(names(df))

    # Harmonise the few layout variants across families.
    if (!"age_group_name" %in% names(df) && "age_group" %in% names(df))
      names(df)[names(df) == "age_group"] <- "age_group_name"
    if (!"mean" %in% names(df) && "value" %in% names(df))
      names(df)[names(df) == "value"] <- "mean"

    if (!all(c("adm0_name", "adm2_name", "year", "mean") %in% names(df))) next

    fam <- family_of(f)
    df_f <- df[df$adm0_name %in% country_name & df$year == year_filter, ,
               drop = FALSE]
    df_f <- df_f[!is.na(df_f$adm2_name) & df_f$adm2_name != "", , drop = FALSE]
    if (nrow(df_f) == 0) next
    fams_seen <- union(fams_seen, fam)

    # Stratifier columns (fill missing with "" so the unite key is stable).
    get <- function(col) if (col %in% names(df_f)) as.character(df_f[[col]]) else ""
    measure <- get("measure")

    # Measure renames (mirror SL_IHME_clean.R) so generic Incidence/Prevalence
    # from different families don't collide in the wide pivot and the column
    # names match the existing countries.
    if (fam == "u5 diarrhea") {
      measure[measure == "Incidence"]  <- "u5 diarrhea incidence"
      measure[measure == "Mortality"]  <- "u5 diarrhea mortality"
      measure[measure == "Prevalence"] <- "u5 diarrhea prevalence"
    }
    if (fam == "LF") {
      measure[measure == "Prevalence"] <- "LF prevalence"
    }
    # (Malaria / onchocerciasis renames are no-ops here — no admin-2 CSVs.)

    long_list[[length(long_list) + 1L]] <- data.frame(
      adm2_name      = as.character(df_f$adm2_name),
      age_group_name = get("age_group_name"),
      age_group_id   = get("age_group_id"),
      sex            = get("sex"),
      measure        = measure,
      metric         = get("metric"),
      value          = suppressWarnings(as.numeric(df_f$mean)),
      family         = fam,
      stringsAsFactors = FALSE)
  }

  if (length(long_list) == 0) {
    warning(sprintf("No admin-2 IHME rows for %s %d.",
                    paste(country_name, collapse = "/"), year_filter))
    return(invisible(NULL))
  }

  long <- dplyr::bind_rows(long_list) |>
    dplyr::filter(!is.na(value), !is.na(adm2_name), adm2_name != "")

  cat(sprintf("[ihme-admin2] families with data: %s\n",
              paste(sort(fams_seen), collapse = ", ")))
  missing_fams <- setdiff(c("Anemia","CGF","DBM","EBF","HIV","LF","Malaria",
                            "mcv1","oral rehydration","u5 diarrhea","u5 mort",
                            "WASH access","education"), fams_seen)
  if (length(missing_fams))
    cat(sprintf("[ihme-admin2] NOTE: no admin-2 rows for: %s\n",
                paste(missing_fams, collapse = ", ")))

  long <- long |>
    tidyr::unite("variable", age_group_name, age_group_id, sex, measure, metric,
                 sep = "_", remove = TRUE) |>
    dplyr::group_by(adm2_name, variable) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  wide <- tidyr::pivot_wider(long, id_cols = adm2_name,
                             names_from = variable, values_from = value) |>
    janitor::clean_names()

  # Prefix ihme_, then strip the clean_names numeric-guard "x" (mirrors
  # SL_IHME_clean.R: colnames <- gsub("_x","_", colnames)) so columns read
  # ihme_15_49_years_... exactly like the existing countries.
  names(wide) <- paste0("ihme_", names(wide))
  names(wide) <- gsub("_x", "_", names(wide))
  names(wide)[names(wide) == "ihme_adm2_name"] <- "ihme_adm2_name"  # keep id col

  cat(sprintf("[ihme-admin2] %s %d: %d admin-2 units x %d variables\n",
              paste(country_name, collapse = "/"), year_filter,
              nrow(wide), ncol(wide) - 1L))

  # ── Column-parity diagnostic vs an existing country ──────────────────────
  if (!is.null(reference_csv) && file.exists(reference_csv)) {
    ref_cols <- names(readr::read_csv(reference_csv, n_max = 0,
                                      show_col_types = FALSE))
    new_cols <- names(wide)
    shared   <- intersect(new_cols, ref_cols)
    cat(sprintf(
      "[ihme-admin2] PARITY vs %s: %d/%d new cols match (%.0f%%); %d ref-only, %d new-only\n",
      basename(reference_csv), length(shared), length(new_cols),
      100 * length(shared) / max(1, length(new_cols)),
      length(setdiff(ref_cols, new_cols)), length(setdiff(new_cols, ref_cols))))
    no_match <- setdiff(new_cols, ref_cols)
    if (length(no_match))
      cat("  new-only (first 15): ",
          paste(head(no_match, 15), collapse = " | "), "\n", sep = "")
  }

  if (!is.null(out_path)) {
    readr::write_csv(wide, out_path)
    cat(sprintf("[ihme-admin2] wrote %s\n", out_path))
  }
  invisible(wide)
}
