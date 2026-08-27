# =============================================================================
# R/dhs_harmonization.R
#
# Harmonizes each country's DHS admin2-level survey indicators into a common
# cross-country naming scheme so they can enter the pooled / LOCO
# transportability analysis (see COMMON_DOMAIN_PREFIXES in R/transportability.R).
#
# --- The problem -------------------------------------------------------------
# src/DHS/DHS_admin2_aggregation.R and src/DHS/DHS_custom_admin2_indicators.R
# both write wide-format Admin2 tables with columns named
#     dhs<year>_<indicator>_adm2
# where <year> is each country's OWN DHS survey year (Gambia 2019, Ghana 2014,
# Sierra Leone 2013, Malawi 2015). merge_external.R's load_dhs_admin2() /
# merge_external_predictors() then merge these into the individual-level
# dataset, so every row already carries its area's dhs<year>_..._adm2 value.
#
# Because the year is baked into the column name, `Reduce(intersect, ...)`
# across countries in build_pooled_dataset() finds ZERO dhs*-prefixed columns
# in common -- the entire DHS admin2 domain is invisible to the pooled/LOCO
# model, even though these household-survey indicators (vaccination coverage,
# ANC visits, sanitation, wealth quintile, stunting, ...) are plausibly very
# relevant to micronutrient deficiency.
#
# --- Why this is fixable with simple prefix-stripping (no crosswalk needed) --
# Investigated whether metadata/dhs_indicator_lookup.csv / dhs_prefix_lookup.csv
# (a StatCompiler indicator-code crosswalk) would be needed to match different
# per-round column names to the same concept. It is NOT needed here, because
# the <indicator> part of the column name is ALREADY identical across
# countries/rounds for both indicator families this pipeline produces:
#
#   1. surveyPrev built-in indicators (DHS_admin2_aggregation.R): <indicator>
#      is the surveyPrev/DHS StatCompiler indicator ID itself (e.g.
#      "CH_VACC_C_MSL", "CN_NUTS_C_HA2") -- these are DHS's own standardized
#      indicator codes, defined with a fixed numerator/denominator in
#      metadata/dhs_indicators.csv, and used identically across surveys BY
#      DESIGN (that is the entire point of a StatCompiler code).
#
#   2. Custom indicators (DHS_custom_admin2_indicators.R): <indicator> is a
#      hand-rolled name (e.g. "c_stunted", "hh_improved_sanitation",
#      "w_anc4plus") computed by the SAME derive_*() function, reading the
#      SAME DHS *standard recode* variable (hc70, hv205, m14_1, ...) for every
#      country. DHS's standard recode variables are kept consistent across
#      phases/rounds by DHS itself, so these names/definitions are also
#      identical across countries.
#
# Confirmed empirically (2026, this session): loading the 4 countries' merged
# individual-level datasets and stripping the "dhs<year>_" prefix, 107 of the
# resulting indicator codes are present in ALL FOUR countries, with 0-17.7%
# missingness (see harmonize_dhs_admin2() below + verification script output
# in the session notes). No indicator crosswalk table was needed.
#
# --- A pipeline quirk this function has to route around ----------------------
# In the merged individual-level data, every genuine dhs<year>_<ind>_adm2
# column is present TWICE, suffixed ".x" and ".y" (an upstream duplicate-merge
# artifact -- not something this file fixes, since the task here is narrowly
# harmonization). Also present: a handful of dhs<year>_<ind> columns with NO
# "_adm2" suffix that are NOT admin2-resolved (constant, or resolved at a much
# coarser level than Admin2) -- these are excluded. harmonize_dhs_admin2()
# only looks at "..._adm2[.x/.y]" columns and coalesces the .x/.y duplicates
# (preferring .x, falling back to .y) rather than assuming they're identical
# -- for Gambia/Sierra Leone they are bit-identical, but for Ghana/Malawi a
# handful of areas differ slightly between .x and .y (likely direct-vs-smoothed
# estimate provenance from the same upstream duplication), so coalescing is
# safer than picking one side and silently discarding the other's information.
#
# --- Deliberately NOT harmonized ---------------------------------------------
# See DHS_HARMONIZE_EXCLUDE below for indicators present in all 4 countries
# that were held back on purpose, with reasons.
# =============================================================================

#' Indicator codes to withhold from harmonization even though the raw
#' `dhs<year>_<code>_adm2` column exists in all 4 countries.
#'
#' Being present in all 4 countries under the same code is necessary but not
#' sufficient for being *poolable* -- a couple of indicators are known, by
#' construction, to not be on a comparable scale/denominator across countries.
#' Per the task's conservative instruction ("a wrong harmonization is worse
#' than leaving one out"), these are excluded rather than guessed at:
#'
#'   - "hh_wealth_mean": DHS's household wealth index factor score (hv271).
#'     DHS constructs this via country-specific PCA over each survey's own
#'     asset list, standardized to mean ~0 within that SURVEY. DHS's own
#'     documentation is explicit that the raw wealth index score is NOT
#'     comparable in absolute terms across countries or over time (only
#'     within-survey rank/quintile is intended for comparison). This is the
#'     same class of problem as raw WFP commodity prices (R/external_data.R) --
#'     see extract_wfp_price_deviation()'s docs for the analogous fix there.
#'     `hh_wealth_poor` / `hh_wealth_poorest` (quintile-membership binaries,
#'     each ~20%/40% of a country's OWN population by construction) are kept:
#'     a relative-rank concept is at least consistently defined across
#'     countries, unlike the raw factor score's scale.
#'
#'   - "c_diet_diversity_score": a raw count of food groups consumed (from
#'     DHS_custom_admin2_indicators.R's `derive_child_nutrition_KR()`), built
#'     from `intersect(paste0("v414", letters[1:23]), names(KRdata))` -- i.e.
#'     the denominator (how many v414* food-group questions were actually
#'     asked in a given country's questionnaire) is not verified identical
#'     across the 4 surveys from the merged data alone (would need the raw
#'     KR recodes to confirm). A raw count is directly sensitive to that
#'     denominator, so it's excluded. `c_min_diet_diversity` (the binarized
#'     ">=4 food groups" MDD-C threshold) is kept: it is DHS/WHO's own
#'     standardized indicator definition and far less sensitive to a couple
#'     of missing food-group items than the raw count is.
#'
#' Everything else in the common set was judged safe: DHS StatCompiler
#' indicator IDs (fixed definitions in metadata/dhs_indicators.csv), z-scores
#' against the WHO growth standard (haz/waz/whz -- internationally
#' standardized by construction), and physical/absolute measures (BMI,
#' hemoglobin g/dL, birthweight grams, years of education, household size,
#' age, parity) that don't have a country-specific scale problem.
DHS_HARMONIZE_EXCLUDE <- c("hh_wealth_mean", "c_diet_diversity_score")


#' Harmonize one country's dhs<year>_<indicator>_adm2 columns into a common
#' cross-country vocabulary.
#'
#' Adds new `dhsharm_<indicator>_adm2` columns (originals are left in place,
#' so within-country / single-survey-year uses of the raw dhs<year>_ columns
#' are unaffected). Safe to call on data that has no dhs*_adm2 columns at all
#' (returns `d` unchanged).
#'
#' @param d One country's merged individual-level data.frame (as produced by
#'   load_merged_data() / used as an element of `all_merged` in
#'   build_pooled_dataset()).
#' @param exclude_indicators Character vector of indicator codes (the part of
#'   the column name between "dhs<year>_" and "_adm2") to skip. Defaults to
#'   DHS_HARMONIZE_EXCLUDE; pass character(0) to harmonize everything found.
#' @return `d` with additional `dhsharm_<indicator>_adm2` columns.
harmonize_dhs_admin2 <- function(d, exclude_indicators = DHS_HARMONIZE_EXCLUDE) {

  cols <- colnames(d)

  # Only true admin2-resolved indicator columns: dhs<4-digit-year>_<code>_adm2,
  # optionally with the upstream ".x"/".y" duplicate-merge suffix. Excludes
  # the handful of dhs<year>_<code> columns with no "_adm2" suffix, which are
  # NOT admin2-resolved (see file header) and would silently corrupt the
  # pooled model if treated as area-level covariates.
  dhs_adm2 <- grep("^dhs[0-9]{4}_.+_adm2([.][xy])?$", cols, value = TRUE)
  if (length(dhs_adm2) == 0) return(d)

  # Canonical indicator name: strip the country-specific year prefix and the
  # .x/.y suffix, keep "_adm2". This is identical across countries by
  # construction (see file header) -- no lookup table needed.
  canon <- sub("^dhs[0-9]{4}_(.+?)_adm2([.][xy])?$", "\\1_adm2", dhs_adm2)

  by_canon <- split(dhs_adm2, canon)

  n_added <- 0L
  for (ind_adm2 in names(by_canon)) {
    ind_code <- sub("_adm2$", "", ind_adm2)
    if (ind_code %in% exclude_indicators) next

    raw_cols <- by_canon[[ind_adm2]]
    vals <- d[[raw_cols[1]]]
    if (length(raw_cols) > 1) {
      # Coalesce rather than assume identical -- see file header: Ghana/Malawi
      # .x vs .y differ for a handful of areas.
      for (rc in raw_cols[-1]) {
        vals <- ifelse(is.na(vals), d[[rc]], vals)
      }
    }

    harm_col <- paste0("dhsharm_", ind_adm2)
    d[[harm_col]] <- vals
    n_added <- n_added + 1L
  }

  cat(sprintf("[harmonize_dhs_admin2] %d raw dhs*_adm2 columns -> %d dhsharm_ columns (%d excluded)\n",
              length(dhs_adm2), n_added, sum(names(by_canon) %in% paste0(exclude_indicators, "_adm2"))))

  d
}


#' Diagnostic: apply harmonize_dhs_admin2() to a named list of per-country
#' merged data.frames and report the actual common dhsharm_ column count and
#' per-country missingness. Mirrors the per-domain audit already printed by
#' build_pooled_dataset() in R/transportability.R.
#'
#' @param all_merged Named list of merged data.frames (one per country), as
#'   used elsewhere in this pipeline (e.g. build_pooled_dataset()).
#' @return invisibly, a list(common = character vector of common dhsharm_
#'   columns, missingness = data.frame of per-country % missing for each).
audit_dhs_harmonization <- function(all_merged) {
  harmonized <- lapply(all_merged, harmonize_dhs_admin2)

  per_country_cols <- lapply(harmonized, function(d) {
    grep("^dhsharm_", colnames(d), value = TRUE)
  })

  common <- Reduce(intersect, per_country_cols)
  cat(sprintf("[audit_dhs_harmonization] %d common dhsharm_ columns across %d countries\n",
              length(common), length(all_merged)))

  miss <- data.frame(indicator = sort(common))
  for (cn in names(harmonized)) {
    d <- harmonized[[cn]]
    miss[[cn]] <- vapply(miss$indicator, function(col) mean(is.na(d[[col]])) * 100, numeric(1))
  }

  invisible(list(common = sort(common), missingness = miss))
}
