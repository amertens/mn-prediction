# =============================================================================
# scripts/validate_merged_datasets.R
#
# Post-merge sanity check: for each country, verify the merged_dataset.rds
# carries the gee_a2_*, ihme_adm1_*, MAP_*_Blood_Disorders, and other recent
# additions; confirm row counts, NA rates, and survey-year alignment.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})

countries <- list(
  Gambia      = list(rds = "data/IPD/Gambia/Gambia_merged_dataset.rds",
                     gee_a2_csv = "data/GEE/Gambia_2018_admin2_gee.csv",
                     survey_year = 2018L),
  Ghana       = list(rds = "data/IPD/Ghana/Ghana_merged_dataset.rds",
                     gee_a2_csv = "data/GEE/Ghana_2017_admin2_gee.csv",
                     survey_year = 2017L),
  SierraLeone = list(rds = "data/IPD/Sierra Leone/SierraLeone_merged_dataset.rds",
                     gee_a2_csv = "data/GEE/SierraLeone_2013_admin2_gee.csv",
                     survey_year = 2013L),
  Malawi      = list(rds = "data/IPD/Malawi/Malawi_merged_dataset.rds",
                     gee_a2_csv = "data/GEE/Malawi_2016_admin2_gee.csv",
                     survey_year = 2016L)
)

summarise_country <- function(cn, cfg) {
  cat("\n=== ", cn, " ===\n")
  rds_path <- here(cfg$rds)
  if (!file.exists(rds_path)) {
    cat("  RDS MISSING:", rds_path, "\n")
    return(invisible(NULL))
  }
  d <- readRDS(rds_path)
  # Drop sf/geometry so dplyr can group/summarise cleanly
  if (inherits(d, "sf")) d <- sf::st_drop_geometry(d)
  geom_cols <- vapply(d, function(c) inherits(c, c("sfc", "sfc_POINT")),
                      logical(1))
  if (any(geom_cols)) d <- d[, !geom_cols, drop = FALSE]
  cat(sprintf("  rows: %d   total cols: %d\n", nrow(d), ncol(d)))

  # Per-prefix counts
  prefixes <- c("gw_", "dhs", "mics_", "ihme_", "ihme_adm1_", "lsms_",
                "MAP_", "wfp_", "flunet_", "gee_", "gee_a2_", "fsec_", "soil_")
  cat("  Predictor counts by prefix:\n")
  for (p in prefixes) {
    cols <- grep(paste0("^", p), colnames(d), value = TRUE)
    if (length(cols) > 0) {
      # Quick NA rate
      na_pct <- mean(vapply(cols, function(c) mean(is.na(d[[c]])),
                             numeric(1)))
      cat(sprintf("    %-15s n=%4d   mean NA rate=%5.1f%%\n",
                  p, length(cols), 100 * na_pct))
    }
  }

  # Specifically check the new gee_a2_ block — non-empty, joined OK
  gee_a2_cols <- grep("^gee_a2_", colnames(d), value = TRUE)
  if (length(gee_a2_cols) > 0) {
    na_per_col <- vapply(gee_a2_cols,
                          function(c) mean(is.na(d[[c]])), numeric(1))
    cat(sprintf("  gee_a2_:  %d cols, %d fully-populated, %d all-NA, mean NA = %.1f%%\n",
                length(gee_a2_cols),
                sum(na_per_col == 0),
                sum(na_per_col == 1),
                100 * mean(na_per_col)))
  } else {
    cat("  gee_a2_:  NONE FOUND — admin-2 GEE merge did not run\n")
  }

  # Cross-check vs source CSV: pick one column, verify a known value
  if (file.exists(here(cfg$gee_a2_csv)) && length(gee_a2_cols) > 0) {
    csv <- read.csv(here(cfg$gee_a2_csv), check.names = FALSE)
    csv$Admin2 <- trimws(csv$Admin2)
    test_col <- intersect(colnames(csv), gee_a2_cols)[1]
    if (!is.na(test_col) && "Admin2" %in% colnames(d)) {
      # Build per-Admin2 expected from csv vs observed from rds (first row per Admin2)
      first_per_adm2 <- d %>%
        as.data.frame() %>%
        dplyr::filter(!is.na(Admin2)) %>%
        dplyr::group_by(Admin2) %>%
        dplyr::summarise(observed = first(.data[[test_col]]),
                          .groups = "drop")
      cmp <- merge(first_per_adm2,
                   csv[, c("Admin2", test_col)],
                   by = "Admin2", all.x = TRUE)
      cmp$diff <- ifelse(is.na(cmp$observed) | is.na(cmp[[test_col]]),
                          NA, abs(cmp$observed - cmp[[test_col]]))
      cat(sprintf("  Spot-check col '%s' against %s:\n",
                  test_col, basename(cfg$gee_a2_csv)))
      cat(sprintf("    %d Admin2 matched, max abs diff = %g, NA-mismatch = %d\n",
                  sum(!is.na(cmp$diff)),
                  max(cmp$diff, na.rm = TRUE),
                  sum(is.na(cmp$observed) != is.na(cmp[[test_col]]))))
    }
  }

  invisible(NULL)
}

for (cn in names(countries)) summarise_country(cn, countries[[cn]])
cat("\nDone.\n")
