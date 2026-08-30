# =============================================================================
# R/covariates/load_harmonized.R
#
# Readers for the stage-3 harmonised covariate dataset, for use by the modelling
# pipeline.
#
# These are deliberately additive: nothing here changes the existing analysis of
# record, which continues to read the raster-derived `gee_*` vocabulary through
# extract_area_covariates(). Switching the models over is a separate, explicit
# decision that requires a full pipeline re-run and a re-reported comparison --
# it should not happen as a side effect of this file existing.
#
# To use the harmonised set in the area-level model, pass
# `load_harmonized_admin2(country)` where `gee_admin2` is expected. Column names
# carry no `gee_` prefix, so the two vocabularies cannot be silently mixed.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

.harm_dir <- function() here::here("data", "covariates", "harmonized")

#' Pooled harmonised admin-2 predictors: the canonical variables shared by
#' every country, one row per admin-2 area, with a `country` column.
load_harmonized_pooled <- function() {
  f <- file.path(.harm_dir(), "predictors_admin2_harmonized.csv")
  if (!file.exists(f))
    stop("harmonised covariates not built. Run: Rscript scripts/covariates/run_all.R")
  suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>% as.data.frame()
}

#' One country's slice of the harmonised (shared-vocabulary) table.
load_harmonized_admin2 <- function(country) {
  d <- load_harmonized_pooled()
  out <- d[d$country == country, setdiff(names(d), "country"), drop = FALSE]
  if (!nrow(out))
    stop("no harmonised rows for '", country, "'. Countries present: ",
         paste(sort(unique(d$country)), collapse = ", "))
  out
}

#' One country's FULL canonical table, including variables that are not present
#' in every country. Use for within-country models, where the cross-country name
#' intersection is not a constraint and discarding country-specific covariates
#' costs information for no benefit.
load_canonical_admin2 <- function(country, keep_excluded = FALSE) {
  source(here::here("R", "covariates", "canonicalize.R"), local = TRUE)
  f <- here::here("data", "covariates", "country",
                  sprintf("%s_predictors_admin2_raw.csv", country))
  if (!file.exists(f))
    stop("no stage-2 table for '", country, "'. Run: ",
         "Rscript scripts/covariates/02_build_country_predictors.R ", country)
  raw <- suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>% as.data.frame()
  cov_harmonize_country(raw, country, keep_excluded = keep_excluded)$data
}

#' The data dictionary for the harmonised set.
load_covariate_dictionary <- function() {
  f <- file.path(.harm_dir(), "data_dictionary.csv")
  if (!file.exists(f))
    stop("dictionary not built. Run: Rscript scripts/covariates/04_document_and_qc.R")
  suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>% as.data.frame()
}

#' Canonical variables belonging to one conceptual domain, for domain ablation.
#'
#' @param domain a value of the `domain` column in source_registry.csv, e.g.
#'   "Vegetation", "Soil chemistry", "Land surface (learned)"
harmonized_domain_vars <- function(domain) {
  d <- load_covariate_dictionary()
  d$canonical[!is.na(d$domain) & d$domain == domain]
}
