# =============================================================================
# R/gee_band_semantics.R
#
# Covariate hygiene for the GEE Admin-2 predictors.
#
# THE PROBLEM
# -----------
# .append_gee_zonal_cols() (R/admin2_analysis.R) emits, for every raster with
# <= 24 bands, one column per band PLUS five cross-band summaries:
#
#   <var>_annual_mean  _annual_sd  _annual_min  _annual_max  _annual_range
#
# That is correct when the bands are a temporal series of ONE variable — for
# LST_Night_Monthly (bands 2013_01_Mean ... 2013_12_Mean) "annual mean night
# temperature" is exactly right. But the same code path runs on rasters whose
# bands are DIFFERENT quantities, and then it averages things that do not belong
# together. Verified band lists from data/Sierra_Leone_GEE_rasters/:
#
#   TerraClimate   14 bands: aet, def, pdsi, pet, pr, ro, ...   distinct variables
#   FLDAS          28 bands: Evap_tavg, Psurf_f_tavg, ...       distinct variables
#   iSDAsoil        4 bands: mean_0_20, mean_20_50,             location AND
#                            stdev_0_20, stdev_20_50            dispersion mixed
#   LST_Night_
#     Annual_Mean   2 bands: Mean, FilledProportion             Kelvin (~290) mixed
#                                                               with a 0-1 proportion
#   Productivity    3 bands: Gpp, Npp, Npp_QC                   averages a QUALITY
#                                                               FLAG into the value
#   LandCoverType  13 bands: LC_Type1 ... LC_Prop1_Assessment   categorical CLASS CODES
#   GPW_Demographic 77 bands: age-sex population bins           a scaled total count
#
# Consequences that reach the reported results:
#   * The summary is dominated by whichever band has the largest units. On Sierra
#     Leone, gee_fldas_..._annual_mean = 3640.9, which is surface pressure
#     (~98,000 Pa) / 27 — i.e. an elevation proxy wearing a climate label.
#   * Some summaries are EXACT duplicates of a real band: gee_soilzinc_annual_max
#     is bit-for-bit identical to gee_soilzinc_mean_0_20 (same for iron,
#     potassium, ...), and gee_terraclimate_2012_annual_min is identical to
#     gee_terraclimate_2012_pdsi.
#   * Duplicated columns split variable importance and make the lasso's choice
#     among identical columns arbitrary, which destabilises the selected-variable
#     lists reported per LOCO fold.
#
# On Sierra Leone: 543 Admin-2 gee_ columns, of which 248 (46%) are cross-band
# summaries, and 243 columns are an exact copy of another column.
#
# WHAT THIS FILE DOES
# -------------------
# Declares, per raster family, whether its bands are commensurable, and provides
# a NAME-BASED pruning function. Name-based (not value-based) is deliberate: a
# value-based duplicate filter would drop different columns in different
# countries, which would make the cross-country vocabulary country-dependent and
# reintroduce the silent-intersection-collapse failure this pipeline already had.
#
# Gated behind GEE_COVARIATE_HYGIENE (default OFF) so the change is opt-in and
# can be evaluated against the current behaviour — see
# scripts/compare_covariate_hygiene.R.
# =============================================================================

#' Band semantics per raster family.
#'
#' Keys are regular expressions matched against the family stem (the column name
#' with the `gee_` prefix, any year tokens and any `_annual_*` suffix removed).
#' Values:
#'   "temporal"     bands are a time series of ONE variable -> cross-band
#'                  summaries are meaningful, keep them
#'   "multivariate" bands are different physical quantities -> summaries are not
#'                  meaningful
#'   "categorical"  bands are class codes -> neither summaries nor arithmetic on
#'                  the band values themselves is meaningful
#'
#' Classified from the actual band names in data/<Country>_GEE_rasters/, not from
#' the family name. Add new families here; anything unmatched is reported by
#' gee_unclassified_families() so it cannot be silently forgotten.
GEE_BAND_SEMANTICS <- list(
  # --- temporal: keep the summaries -----------------------------------------
  "^lst_night_monthly$"        = "temporal",      # 12 monthly Mean bands
  "^wapor$"                    = "temporal",      # 36 dekadal NPP slices
  "^dailyevi$"                 = "temporal",
  "^lai8days$"                 = "temporal",

  # --- multivariate: bands are different quantities --------------------------
  "^soil[a-z]+$"               = "multivariate",  # depth means + their stdevs
  "^terraclimate$"             = "multivariate",  # aet/def/pdsi/pet/pr/ro/...
  "^fldas(_annual_mean|_monthly)?$" = "multivariate",  # Evap/Psurf/Qair/Wind/...
  "^lst_night_annual_mean$"    = "multivariate",  # Kelvin + FilledProportion
  "^productivity$"             = "multivariate",  # Gpp + Npp + Npp_QC (a QC flag)
  "^ghsbuilts$"                = "multivariate",  # built_surface + _nres
  "^atmosphere$"               = "multivariate",  # 628 mixed optical quantities
  "^gpw_demographic$"          = "multivariate",  # 77 age-sex count bins
  "^esa_worldcereal$"          = "multivariate",  # per-tile crop markers

  # --- categorical: class codes ---------------------------------------------
  "^landcovertype$"            = "categorical",   # LC_Type1..5 class codes
  "^landcoverlayers$"          = "categorical"    # class code + cover fractions
)

#' Families whose underlying layer is STATIC but which the export re-ran once per
#' survey year, producing identical year-stamped copies (e.g.
#' gee_soilzinc_mean_0_20 == gee_soilzinc_2012_mean_0_20 == ..._2013_...).
#'
#' Matched as a PREFIX, not a whole stem: these columns keep their band suffix
#' after gee_family_stem() (`soilzinc_mean_0_20`), unlike the summary columns
#' whose suffix is stripped (`soilzinc`).
GEE_STATIC_FAMILIES <- "^soil[a-z]+(_|$)"

#' Is covariate hygiene enabled? Off by default; set GEE_COVARIATE_HYGIENE=true.
gee_hygiene_enabled <- function() {
  isTRUE(tolower(Sys.getenv("GEE_COVARIATE_HYGIENE", "false")) %in%
           c("1", "true", "yes"))
}

#' Reduce a gee_ column name to its raster-family stem.
#' Strips the prefix, any 4-digit year tokens, and the cross-band summary suffix.
gee_family_stem <- function(cols) {
  x <- sub("^gee_", "", cols)
  x <- sub("_annual_(mean|sd|min|max|range)$", "", x)
  x <- gsub("_(19|20)[0-9]{2}", "", x)
  x <- gsub("_+", "_", x)
  sub("_$", "", x)
}

#' Look up the declared semantics for family stems. NA when unclassified.
gee_band_semantics_of <- function(stems) {
  out <- rep(NA_character_, length(stems))
  for (pat in names(GEE_BAND_SEMANTICS)) {
    hit <- is.na(out) & grepl(pat, stems)
    out[hit] <- GEE_BAND_SEMANTICS[[pat]]
  }
  out
}

#' Family stems that emit cross-band summaries but are not classified above.
#'
#' Kept visible on purpose: an unclassified family keeps its summaries (so this
#' file can never silently delete a covariate), but it should be classified.
gee_unclassified_families <- function(cols) {
  s <- grep("_annual_(mean|sd|min|max|range)$", cols, value = TRUE)
  if (!length(s)) return(character(0))
  stems <- unique(gee_family_stem(s))
  sort(stems[is.na(gee_band_semantics_of(stems))])
}

#' Cross-band summary columns whose bands are not commensurable.
gee_noncommensurable_summaries <- function(cols) {
  s <- grep("_annual_(mean|sd|min|max|range)$", cols, value = TRUE)
  if (!length(s)) return(character(0))
  sem <- gee_band_semantics_of(gee_family_stem(s))
  s[!is.na(sem) & sem != "temporal"]
}

#' Year-stamped copies of a static layer, when an unstamped copy is also present.
#'
#' Keeps the unstamped column (the one that survives the cross-country
#' intersection) and drops the identical per-year duplicates.
gee_static_year_duplicates <- function(cols) {
  gee <- grep("^gee_", cols, value = TRUE)
  if (!length(gee)) return(character(0))
  stems <- gee_family_stem(gee)
  is_static <- Reduce(`|`, lapply(GEE_STATIC_FAMILIES, function(p) grepl(p, stems)),
                      init = rep(FALSE, length(gee)))
  cand <- gee[is_static & grepl("_(19|20)[0-9]{2}_", gee)]
  if (!length(cand)) return(character(0))
  # The unstamped equivalent is the same name with the year token removed.
  unstamped <- gsub("_(19|20)[0-9]{2}", "", cand)
  cand[unstamped %in% cols]
}

#' Columns to drop from the GEE predictor set under covariate hygiene.
#'
#' @param cols character vector of candidate predictor names
#' @param verbose print a one-line summary of what was dropped
#' @return character vector of column names to drop (empty when hygiene is off)
prune_gee_covariates <- function(cols, verbose = TRUE) {
  if (!gee_hygiene_enabled()) return(character(0))
  drop <- unique(c(gee_noncommensurable_summaries(cols),
                   gee_static_year_duplicates(cols)))
  drop <- intersect(drop, cols)
  if (verbose && length(drop)) {
    unclassified <- gee_unclassified_families(cols)
    cat(sprintf(paste0("[gee_hygiene] dropping %d of %d gee_ columns ",
                       "(%d non-commensurable cross-band summaries, ",
                       "%d static per-year duplicates)%s\n"),
                length(drop), sum(startsWith(cols, "gee_")),
                length(gee_noncommensurable_summaries(cols)),
                length(gee_static_year_duplicates(cols)),
                if (length(unclassified))
                  sprintf("; %d unclassified family/families kept: %s",
                          length(unclassified),
                          paste(unclassified, collapse = ", ")) else ""))
  }
  drop
}
