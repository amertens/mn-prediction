# =============================================================================
# R/covariates/covariate_selection.R
#
# ONE place that answers "which columns of this admin-2 table are covariates?"
#
# WHY THIS EXISTS
# ---------------
# Nine call sites across seven files selected area covariates with a hardcoded
# `grep("^gee_", ...)`. That is correct for the legacy raster vocabulary and
# silently returns ZERO columns for the harmonised one, whose names carry no
# prefix. On the 2026-08-29 harmonised launch that took out the corrected-
# methods analysis of record (p2_p6, p7, p9), the area recipe, transportability
# and feature engineering -- some erroring, others liable to run on an empty
# covariate matrix and produce output that looks fine.
#
# A prefix test is the wrong abstraction: it asks "what is this column called"
# when the question is "is this column a predictor". These helpers ask the
# latter, so a third vocabulary would not break them either.
#
# THE TWO REGIMES
# ---------------
# The area recipe distinguishes a "universal" covariate set from an "enriched"
# one. Under the legacy vocabulary that was simply gee_ (Earth observation)
# versus gee_ + MAP/IHME (modelled epidemiological surfaces).
#
# The harmonised vocabulary folds SoilGrids and MapSPAM into the same table, so
# the boundary has to be drawn on MEANING rather than on prefix. The rule kept
# here is the one the distinction was always about -- availability:
#
#   universal  covariates derivable ANYWHERE on earth from open global products
#              (Earth observation, soil, crop allocation)
#   enriched   universal + sources that require a survey or a modelled
#              epidemiological surface for that specific country (DHS, ESPEN,
#              MAP, IHME)
#
# So DHS admin-2 indicators are enrichment, not universal: a country with no DHS
# has none of them, which is exactly the case the universal regime exists to
# represent.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

# Keys, survey outcome columns, and geometry helpers -- never predictors.
AREA_NON_COVARIATE_COLS <- c(
  "Admin1", "Admin2", "country", "region", "geometry",
  "svy_prev", "svy_prev_se", "svy_prev_low", "svy_prev_upp", "svy_cv",
  "svy_ci_width", "n_svy", "sl_prev", "n_sl", "area_pred", "area_pred_oof",
  "area_pred_prev", "area_pred_prev_unbenchmarked", "mrp_prev",
  "downscale_prev", "lon", "lat")

# Canonical prefixes that are NOT universally available (see header).
AREA_ENRICHMENT_PATTERNS <- c("^dhs_", "^espen_", "^MAP_", "^map2?_", "^ihme_")

#' Which vocabulary is this table written in?
#' @return "legacy" if any gee_ columns are present, "harmonized" otherwise
area_vocab_of <- function(df) {
  if (any(grepl("^gee_", names(df)))) "legacy" else "harmonized"
}

#' Covariate columns of an admin-2 table, in either vocabulary.
#'
#' @param df admin-2 table (keys + covariates, possibly + survey outcomes)
#' @param regime "all" (default) or "universal" -- see the header for what the
#'   distinction means and why it is drawn on availability rather than prefix
#' @param numeric_only drop non-numeric columns (default TRUE; a categorical
#'   covariate is not usable by these estimators without encoding)
#' @param min_unique require at least this many distinct non-missing values, so
#'   constants and all-NA columns never reach a model
#' @return character vector of column names
area_covariate_cols <- function(df, regime = c("all", "universal"),
                                numeric_only = TRUE, min_unique = 2L) {
  regime <- match.arg(regime)
  if (is.null(df) || !ncol(df)) return(character(0))
  cols <- setdiff(names(df), AREA_NON_COVARIATE_COLS)

  if (identical(area_vocab_of(df), "legacy")) {
    # Legacy: the EO block is exactly the gee_ columns; anything else present is
    # enrichment (MAP/IHME merged in by the caller).
    cols <- if (regime == "universal") grep("^gee_", cols, value = TRUE) else cols
  } else if (regime == "universal") {
    drop <- Reduce(`|`, lapply(AREA_ENRICHMENT_PATTERNS,
                               function(p) grepl(p, cols, perl = TRUE)),
                   init = rep(FALSE, length(cols)))
    cols <- cols[!drop]
  }

  if (numeric_only)
    cols <- cols[vapply(cols, function(v) {
      x <- df[[v]]
      is.numeric(x) || is.integer(x) || inherits(x, "haven_labelled")
    }, logical(1))]

  if (min_unique > 1L)
    cols <- cols[vapply(cols, function(v) {
      x <- suppressWarnings(as.numeric(df[[v]]))
      length(unique(x[is.finite(x)])) >= min_unique
    }, logical(1))]

  cols
}

#' Soil-chemistry depth-mean columns, in either vocabulary.
#' legacy  gee_soilzinc_mean_0_20      harmonised  soil_zinc_mean_0_20
area_soil_cols <- function(df) {
  cols <- setdiff(names(df), AREA_NON_COVARIATE_COLS)
  grep("^(gee_)?soil_?[a-z]+_mean_[0-9]+_[0-9]+$", cols, value = TRUE, perl = TRUE)
}

#' Cropland / crop-mix columns, in either vocabulary.
#' legacy  gee_esa_worldcereal_*_annual_mean   harmonised  spam_share_* / spam_prod_*
area_crop_cols <- function(df) {
  cols <- setdiff(names(df), AREA_NON_COVARIATE_COLS)
  c(grep("worldcereal.*(annual_mean|annual_max)$", cols, value = TRUE),
    grep("^spam_(share|prod|parea)_", cols, value = TRUE))
}

#' Is this NAME an area-level covariate, in either vocabulary?
#'
#' Needed where only a character vector of predictor names is available (the
#' pooled individual-level design matrix), so the data-frame helpers above
#' cannot be used. Legacy names are recognised by prefix; harmonised names by
#' membership in the canonical vocabulary that stage 3 wrote out.
#'
#' Falls back to prefix-only when coverage.csv is absent, which keeps a clean
#' checkout working before the covariate pipeline has ever been run.
is_area_covariate_name <- function(x) {
  legacy <- grepl("^gee_", x)
  f <- here::here("data", "covariates", "harmonized", "coverage.csv")
  if (!file.exists(f)) return(legacy)
  canon <- tryCatch(
    suppressMessages(readr::read_csv(f, show_col_types = FALSE))$canonical,
    error = function(e) character(0))
  legacy | x %in% canon
}
