# =============================================================================
# R/covariates/area_vocabulary.R
#
# Which covariate vocabulary the AREA-LEVEL model runs on.
#
#   COVARIATE_VOCAB=legacy      (default) the raster-derived gee_* columns --
#                               the analysis of record as published
#   COVARIATE_VOCAB=harmonized  the canonical set built by
#                               scripts/covariates/run_all.R
#
# The switch lives here rather than in each model so that flipping it changes
# every downstream consumer at once (area_model_*, run_area_comparison,
# build_area_loco_dataset, the benchmark suite) and cannot leave two
# vocabularies mixed in one design matrix -- which is the failure mode that
# produced a zero-column intersection in the first place.
#
# WHAT CHANGES UNDER `harmonized`
# -------------------------------
#   * cross-country shared predictors: 140 -> 208, from 22 source families
#     across 19 providers in 11 conceptual domains
#   * unit defects fixed (Tanzania night LST was Kelvin against Celsius
#     elsewhere; NPP and LAI differed by 1e4 and 10)
#   * categorical class-code means, QC bands and non-commensurable cross-band
#     summaries removed
#   * soil_nitrogen dropped: irreconcilable across countries, and it was inside
#     the previous intersection
#   * DHS admin-2 indicators, SoilGrids, MapSPAM and 64 AlphaEarth dimensions
#     become available to the area-level model for the first time
#
# Because this changes the covariate set, it changes the results. It is opt-in
# on purpose: the published analysis of record must stay reproducible from a
# clean checkout with no environment set.
#
# After flipping it:
#   targets::tar_invalidate(matches("area_|loco|benchmark"))
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

#' Active covariate vocabulary: "legacy" (default) or "harmonized".
covariate_vocabulary <- function() {
  v <- tolower(Sys.getenv("COVARIATE_VOCAB", "legacy"))
  if (!v %in% c("legacy", "harmonized")) {
    warning("COVARIATE_VOCAB='", v, "' not recognised; using 'legacy'")
    return("legacy")
  }
  v
}

#' Is the harmonised vocabulary active AND actually built?
#'
#' Returns FALSE with a loud message when the switch is on but stage 3 has not
#' been run, so a mis-set environment variable degrades to the legacy set with
#' an explanation rather than to an empty covariate matrix.
harmonized_vocab_active <- function(country = NULL) {
  if (!identical(covariate_vocabulary(), "harmonized")) return(FALSE)
  f <- here::here("data", "covariates", "country",
                  sprintf("%s_predictors_admin2_raw.csv", cov_country_key(country %||% "")))
  if (!is.null(country) && !file.exists(f)) {
    # STOP, do not degrade. A per-country fallback is far more dangerous than a
    # failure: the countries that resolve get canonical names, the one that does
    # not gets legacy gee_* names, and the pooled/LOCO intersection of the two
    # is ~zero -- the models then run on almost no covariates and look fine.
    #
    # That is exactly what happened on 2026-08-29: Gambia, Ghana and Malawi
    # picked up 220-221 harmonised covariates while Sierra Leone silently fell
    # back, because get_country_configs() spells it "Sierra Leone" and the
    # stage-2 artefact is keyed "SierraLeone". The run had to be killed.
    # cov_country_key() fixes that spelling gap; this stop() makes sure any
    # future gap is loud instead of silent.
    stop(sprintf(paste0("COVARIATE_VOCAB=harmonized but no stage-2 table for '%s' ",
                        "(looked for %s). Run scripts/covariates/run_all.R. ",
                        "Refusing to fall back for one country only -- that would ",
                        "mix two covariate vocabularies in the pooled models."),
                 country, basename(f)))
  }
  TRUE
}

#' Attach harmonised canonical covariates to an admin-2 base table.
#'
#' @param base data.frame carrying Admin2 (and ideally Admin1) in the order the
#'   caller needs preserved -- for the area model this is GADM polygon order,
#'   which the prediction-to-map step depends on.
#' @param country country name as used in get_country_configs()
#' @param shared_only TRUE keeps only canonical variables present in every
#'   country (what the pooled/LOCO models can use); FALSE keeps the country's
#'   full canonical set (right for within-country models, where the
#'   cross-country intersection is not a constraint).
#' @return `base` with canonical covariate columns appended, same rows, same order.
append_harmonized_admin2 <- function(base, country, shared_only = FALSE) {
  src <- here::here("R", "covariates")
  source(file.path(src, "canonicalize.R"), local = TRUE)

  f <- here::here("data", "covariates", "country",
                  sprintf("%s_predictors_admin2_raw.csv", cov_country_key(country)))
  if (!file.exists(f)) {
    warning("[vocab] no stage-2 table for ", country, "; returning base unchanged")
    return(base)
  }
  raw <- suppressMessages(readr::read_csv(f, show_col_types = FALSE)) |> as.data.frame()
  can <- cov_harmonize_country(raw, country)$data

  if (shared_only) {
    cf <- here::here("data", "covariates", "harmonized", "coverage.csv")
    if (file.exists(cf)) {
      cvg <- suppressMessages(readr::read_csv(cf, show_col_types = FALSE))
      shared <- cvg$canonical[cvg$in_all]
      keep <- c(intersect(c("Admin1", "Admin2"), names(can)), intersect(shared, names(can)))
      can <- can[, keep, drop = FALSE]
    } else {
      warning("[vocab] coverage.csv missing; keeping the full canonical set")
    }
  }

  # Join on the (Admin1, Admin2) PAIR when both sides carry Admin1. GADM
  # admin-2 names repeat within a country (Malawi: 256 polygons, 243 names), so
  # a name-only join silently multiplies rows. match() rather than a join keeps
  # `base` in its original order and its original row count, both of which the
  # area model relies on.
  # The separator is load-bearing: pasting with sep="" would make ("AB","C")
  # and ("A","BC") the same key and silently cross-join two different
  # districts. ASCII unit separator cannot occur in a GADM name.
  KEY_SEP <- "\u001f"
  pair <- all(c("Admin1", "Admin2") %in% names(base)) &&
          all(c("Admin1", "Admin2") %in% names(can))
  key <- function(d) if (pair) paste(trimws(d$Admin1), trimws(d$Admin2), sep = KEY_SEP)
                     else trimws(as.character(d$Admin2))
  idx <- match(key(base), key(can))

  covs <- setdiff(names(can), c("Admin1", "Admin2"))
  for (v in covs) base[[v]] <- can[[v]][idx]

  cat(sprintf("[vocab] %s: %d harmonised covariates joined by %s (%d/%d areas matched)\n",
              country, length(covs), if (pair) "(Admin1, Admin2)" else "Admin2",
              sum(!is.na(idx)), nrow(base)))
  base
}

`%||%` <- function(a, b) if (is.null(a)) b else a
