# =============================================================================
# scripts/pull_stevens2022.R
#
# Inspect & summarise the publicly-released Stevens et al. 2022 Lancet GH
# "hidden hunger" data files (downloaded from
# https://github.com/GAINAlliance/hiddenhunger).
#
# The repo includes:
#   - NPW_input_data.dta        Non-pregnant women individual-survey rows
#   - PSC_input_data.dta        Preschool-aged children individual-survey rows
#   - country_covariates.dta    Country-level covariates used in the model
# Model posterior outputs are NOT in the repo (Stata .do scripts only).
#
# What this script does:
#   1. Read all three .dta files
#   2. Print schema/summary for each
#   3. Filter to target countries (Gambia/Ghana/SLE/Malawi/CIV/Kenya/TZA/ETH)
#   4. Save filtered subsets to results/external/stevens2022_*.csv as the
#      gold-standard pooled comparison dataset
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(haven); library(readr)
})

CACHE  <- "C:/Users/andre/OneDrive/Documents/mn-prediction/data/external_cache/stevens2022"
OUT    <- here::here("results", "external")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

# A loose target-country regex (Stevens uses ISO names — Cote d'Ivoire with accent etc.)
TARGETS_RE <- paste(c(
  "gambia", "ghana", "sierra leone", "malawi",
  "cote d'ivoire", "côte d'ivoire", "ivoire",
  "kenya", "tanzania", "ethiopia"
), collapse = "|")

read_and_describe <- function(label, path) {
  cat(sprintf("\n=========== %s ===========\n", label))
  cat(sprintf("file: %s\n", path))
  if (!file.exists(path)) { cat("  [missing]\n"); return(NULL) }
  d <- haven::read_dta(path)
  cat(sprintf("rows: %d, cols: %d\n", nrow(d), ncol(d)))
  cat("columns:\n")
  print(colnames(d))

  # If country column present (most likely as country / Country / iso / iso3),
  # report unique countries and how many target rows.
  country_col <- intersect(c("country", "Country", "iso", "iso3", "ISO3", "ISO"),
                            colnames(d))
  if (length(country_col)) {
    cc <- country_col[1]
    cat(sprintf("country column: %s\n", cc))
    vals <- unique(d[[cc]])
    cat(sprintf("unique countries: %d\n", length(vals)))
    mask <- grepl(TARGETS_RE, tolower(d[[cc]]))
    cat(sprintf("rows in target countries: %d\n", sum(mask)))
    if (any(mask)) {
      cat("target-country rows preview:\n")
      print(head(d[mask, ], 3))
    }
  }
  d
}

npw <- read_and_describe("NPW_input_data.dta (non-pregnant women)",
                          file.path(CACHE, "NPW_input_data.dta"))
psc <- read_and_describe("PSC_input_data.dta (preschool children)",
                          file.path(CACHE, "PSC_input_data.dta"))
cov <- read_and_describe("country_covariates.dta",
                          file.path(CACHE, "country_covariates.dta"))

# Save filtered subsets (if a country column exists)
maybe_save <- function(d, label) {
  if (is.null(d)) return(invisible(NULL))
  country_col <- intersect(c("country", "Country", "iso", "iso3", "ISO3", "ISO"),
                            colnames(d))
  if (!length(country_col)) {
    readr::write_csv(d, file.path(OUT, sprintf("stevens2022_%s.csv", label)))
    cat(sprintf("[save] full file -> stevens2022_%s.csv (no country col)\n", label))
    return(invisible(NULL))
  }
  cc   <- country_col[1]
  mask <- grepl(TARGETS_RE, tolower(d[[cc]]))
  if (!any(mask)) {
    cat(sprintf("[save] %s: 0 target rows; saving full file\n", label))
    out <- d
  } else {
    out <- d[mask, ]
  }
  readr::write_csv(out, file.path(OUT, sprintf("stevens2022_%s.csv", label)))
  cat(sprintf("[save] stevens2022_%s.csv (%d rows)\n", label, nrow(out)))
}
maybe_save(npw, "npw")
maybe_save(psc, "psc")
maybe_save(cov, "country_covariates")
