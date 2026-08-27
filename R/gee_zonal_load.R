# =============================================================================
# R/gee_zonal_load.R
#
# Read Admin-2 GEE covariates from the ZONAL CSVs written by
# src/GEE/extract_gee_ee_api.R, rather than re-deriving them from GeoTIFFs.
#
# WHY
# ---
# The pipeline had two parallel covariate paths:
#
#   rasters  data/<Country>_GEE_rasters/*.tif  -> canonicalize_gee_varname()
#            used by predict_oos_pooled()
#   zonal    data/GEE/<Country>_<year>_admin2_gee.csv (gee_a2_ columns)
#            used by the area-level and LOCO models
#
# Two paths means two chances for a country to end up in a private vocabulary,
# and the pooled model intersects covariates by EXACT NAME, so one odd country
# zeroes the shared block for everyone. That is what happened with Tanzania.
# Collapsing onto the zonal CSVs removes the second path entirely.
#
# THE CATCH, MEASURED 2026-08-24
# ------------------------------
# The zonal CSVs on disk today were written by three different generations of
# the extractor and DO NOT share a vocabulary:
#
#   Ghana / Malawi / SierraLeone   gee_a2_AerosolOptical_2017   family + year
#   Gambia                         gee_a2_NDVI_NDVI             family + band
#   Tanzania                       gee_a2_NDVI                  family only
#
# Their column intersection is 0 -- worse than the rasters' 56 across four
# countries. So reading them is only correct AFTER a consistent re-extraction:
#
#   Rscript scripts/refresh_gee_coverage.R --run --rebuild
#
# Rather than paper over three vintages with renaming heuristics -- which would
# silently pair columns that are not the same measurement -- these helpers check
# whether the vocabularies actually agree and report honestly when they do not.
# The caller then keeps using rasters until the rebuild has been run.
# =============================================================================

#' Path to a country's Admin-2 zonal GEE CSV.
#' @param country_name display name as used in get_country_configs()
#' @param cc that country's config entry (for survey_year)
gee_zonal_csv_path <- function(country_name, cc) {
  here::here("data", "GEE",
             sprintf("%s_%s_admin2_gee.csv", country_name, cc$survey_year))
}

#' Load one country's Admin-2 zonal covariates.
#'
#' Column names are used EXACTLY as written by the extractor. No renaming, no
#' year-stripping, no fuzzy matching: after a `--rebuild` every country is
#' written by the same code path and the names line up by construction. Any
#' normalisation applied here would be a heuristic guess about whether two
#' differently-named columns are the same measurement, which is precisely the
#' class of silent error this consolidation exists to remove.
#'
#' @return data.frame with Admin2 + gee_a2_* columns, or NULL if absent.
load_gee_zonal_admin2 <- function(country_name, cc, quiet = FALSE) {
  path <- gee_zonal_csv_path(country_name, cc)
  if (!file.exists(path)) {
    if (!quiet) cat(sprintf("[gee_zonal] %s: no zonal CSV at %s\n",
                            country_name, basename(path)))
    return(NULL)
  }
  d <- tryCatch(
    suppressMessages(readr::read_csv(path, show_col_types = FALSE)),
    error = function(e) {
      if (!quiet) cat(sprintf("[gee_zonal] %s: unreadable (%s)\n",
                              country_name, conditionMessage(e)))
      NULL
    })
  if (is.null(d)) return(NULL)
  if (!"Admin2" %in% names(d)) {
    if (!quiet) cat(sprintf("[gee_zonal] %s: no Admin2 column, skipping\n", country_name))
    return(NULL)
  }
  keep <- c("Admin2", grep("^gee_a2_", names(d), value = TRUE))
  d <- d[, keep, drop = FALSE]

  # Same key hygiene as everywhere else, so a later join cannot fan out on
  # duplicate district names or carry lake polygons in as districts.
  d <- tryCatch(clean_admin2_keys(d, what = paste(country_name, "zonal GEE")),
                error = function(e) d)
  if (!quiet) cat(sprintf("[gee_zonal] %s: %d areas x %d covariates\n",
                          country_name, nrow(d), ncol(d) - 1L))
  d
}

#' Can the zonal CSVs serve this set of countries?
#'
#' All-or-nothing by design. Serving some countries from the zonal CSVs and
#' others from rasters would mix two vocabularies in one model matrix, which
#' produces an empty or -- worse -- a partially-populated intersection that
#' still looks plausible. Either every country resolves, or none do.
#'
#' @param countries named list: country name -> config entry
#' @param min_shared minimum shared covariates to consider the switch usable
#' @return list(ok, data, shared, reason)
gee_zonal_available <- function(countries, min_shared = 20L, quiet = FALSE) {
  dfs <- list()
  for (cn in names(countries)) {
    d <- load_gee_zonal_admin2(cn, countries[[cn]], quiet = quiet)
    if (is.null(d)) {
      return(list(ok = FALSE, data = NULL, shared = character(0),
                  reason = sprintf("%s has no usable zonal CSV", cn)))
    }
    dfs[[cn]] <- d
  }
  shared <- Reduce(intersect, lapply(dfs, function(d) setdiff(names(d), "Admin2")))
  if (!quiet) {
    cat(sprintf("[gee_zonal] %d covariates shared across %d countries\n",
                length(shared), length(dfs)))
  }
  if (length(shared) < min_shared) {
    return(list(
      ok = FALSE, data = dfs, shared = shared,
      reason = sprintf(
        paste("only %d covariate(s) shared across %d countries (need >= %d).",
              "The zonal CSVs were written by different generations of the",
              "extractor and use incompatible column names. Run",
              "`Rscript scripts/refresh_gee_coverage.R --run --rebuild` to",
              "regenerate them all in one vocabulary."),
        length(shared), length(dfs), min_shared)))
  }
  list(ok = TRUE, data = dfs, shared = shared, reason = NA_character_)
}
