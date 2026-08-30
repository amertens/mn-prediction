# =============================================================================
# R/covariates/canonicalize.R
#
# The covariate harmonisation engine. Turns each country's raw, export-specific
# covariate vocabulary into ONE canonical vocabulary shared by every country.
#
# THE PROBLEM IT SOLVES
# ---------------------
# The pooled / LOCO models match covariates by EXACT NAME. Measured 2026-08 on
# the built store: 1,515 distinct admin-2 covariate names exist across the five
# countries but only 140 are shared, and 99 of those 140 are the iSDAsoil block.
# Yet 33 of 35 raster families are present in ALL FIVE countries -- the loss is
# almost entirely a naming artefact, because dynamic layers carry the survey
# year in the column name (Gambia 2018 -> _2017/_2018, Sierra Leone 2013 ->
# _2012/_2013) and export vintages differ in how bands were suffixed.
#
# HOW IT WORKS
# ------------
# Four data files drive everything; no rule lives in code:
#   metadata/covariates/source_registry.csv      what each family IS
#   metadata/covariates/harmonization_rules.csv  raw name pattern -> canonical
#   metadata/covariates/unit_conversions.csv     values -> canonical units
#   metadata/covariates/exclusions.csv           canonical vars that are unsafe
#
# A rule is (order, regex, action, canonical template, collapse policy, reason).
# First match by `order` wins. `canonical` may use backreferences. `collapse`
# says what to do when several raw columns in ONE country map to the same
# canonical name:
#   none    1:1, nothing to do (calendar-year series kept as distinct columns)
#   latest  keep the column with the highest embedded year (2-year exports)
#   mean    average them (month climatologies spanning several years)
#
# EVERY dropped column keeps its reason, so the audit trail is complete and a
# reviewer can see what was removed and why.
#
# ADDING A SOURCE: add a row to source_registry.csv and one rule to
# harmonization_rules.csv. No code change.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

`%||%` <- function(a, b) if (is.null(a)) b else a

.cov_meta_dir <- function() here::here("metadata", "covariates")

cov_load_rules <- function()
  readr::read_csv(file.path(.cov_meta_dir(), "harmonization_rules.csv"),
                  show_col_types = FALSE) %>% dplyr::arrange(order)
cov_load_units <- function()
  readr::read_csv(file.path(.cov_meta_dir(), "unit_conversions.csv"), show_col_types = FALSE)
cov_load_exclusions <- function()
  readr::read_csv(file.path(.cov_meta_dir(), "exclusions.csv"), show_col_types = FALSE)
cov_load_registry <- function()
  readr::read_csv(file.path(.cov_meta_dir(), "source_registry.csv"), show_col_types = FALSE)

#' Canonical file/scope key for a country name.
#'
#' The pipeline carries TWO spellings of the same country: get_country_configs()
#' uses the display name ("Sierra Leone") while every covariate artefact is keyed
#' on the compact form ("SierraLeone"). On the first harmonised launch that
#' mismatch made Sierra Leone alone fall back to the legacy vocabulary while the
#' other three used the canonical one -- a pooled model whose name intersection
#' is ~zero, which is the exact failure this whole exercise exists to prevent.
#' Every path lookup and every country-scoped unit rule goes through here.
cov_country_key <- function(country) gsub("[^A-Za-z0-9]", "", as.character(country))

#' Strip the pipeline's covariate prefixes, leaving the export-native stem.
cov_base_name <- function(v) sub("^gee_(a2_)?", "", v)

#' Highest 4-digit year embedded in a name (NA when there is none).
cov_year_of <- function(v) {
  m <- regmatches(v, gregexpr("(?<![0-9])(19|20)[0-9]{2}(?![0-9])", v, perl = TRUE))
  vapply(m, function(x) if (!length(x)) NA_integer_ else max(as.integer(x)), integer(1))
}

#' Map raw covariate names onto the canonical vocabulary.
#'
#' @param vars character vector of raw column names (with or without gee_ prefix)
#' @param rules rules table; defaults to the on-disk file
#' @return data.frame: variable, base, year, rule_order, action, canonical,
#'   collapse, reason. Names matching no rule get action "unmatched" so a new
#'   export vintage is reported rather than silently dropped.
cov_map_names <- function(vars, rules = cov_load_rules()) {
  base <- cov_base_name(vars)
  out <- data.frame(variable = vars, base = base, year = cov_year_of(base),
                    rule_order = NA_integer_, action = "unmatched",
                    canonical = NA_character_, collapse = NA_character_,
                    reason = "no rule matched this name shape",
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(rules))) {
    hit <- out$action == "unmatched" & grepl(rules$match_regex[i], out$base, perl = TRUE)
    if (!any(hit)) next
    out$rule_order[hit] <- rules$order[i]
    out$action[hit]     <- rules$action[i]
    out$reason[hit]     <- rules$reason[i]
    out$collapse[hit]   <- rules$collapse[i]
    if (identical(rules$action[i], "keep"))
      out$canonical[hit] <- sub(rules$match_regex[i], rules$canonical[i],
                                out$base[hit], perl = TRUE)
  }
  out
}

#' Unit conversion factors for a set of canonical names in one country.
#'
#' First matching row wins, so a country-scoped row must be ordered before a
#' broader one for the same pattern in unit_conversions.csv.
cov_unit_factors <- function(canon, country, units = cov_load_units()) {
  n <- length(canon)
  mult <- rep(1, n); add <- rep(0, n)
  unit <- rep(NA_character_, n); why <- rep(NA_character_, n)
  applied <- rep(FALSE, n)
  for (i in seq_len(nrow(units))) {
    scope <- cov_country_key(trimws(strsplit(units$country[i], ";")[[1]]))
    if (!("ALL" %in% scope || cov_country_key(country) %in% scope)) next
    hit <- !applied & grepl(units$canonical_regex[i], canon, perl = TRUE)
    if (!any(hit)) next
    mult[hit] <- units$multiply[i]; add[hit] <- units$add[i]
    unit[hit] <- units$canonical_unit[i]; why[hit] <- units$reason[i]
    applied[hit] <- TRUE
  }
  data.frame(canonical = canon, multiply = mult, add = add,
             canonical_unit = unit, unit_reason = why, stringsAsFactors = FALSE)
}

#' Harmonise ONE country's admin-2 covariate table.
#'
#' @param df data.frame with Admin1/Admin2 keys plus raw covariate columns
#' @param country country name (selects country-scoped unit conversions)
#' @param keep_excluded keep canonical variables listed in exclusions.csv
#'   (they are always recorded in the map; this only controls the data columns)
#' @return list(data = harmonised wide table, map = full per-column audit trail)
cov_harmonize_country <- function(df, country, keep_excluded = FALSE,
                                  rules = cov_load_rules(),
                                  units = cov_load_units(),
                                  exclusions = cov_load_exclusions()) {
  keys <- intersect(c("Admin1", "Admin2"), names(df))
  vars <- setdiff(names(df), keys)
  map  <- cov_map_names(vars, rules)
  map$country <- country

  canon_chr <- ifelse(is.na(map$canonical), "", map$canonical)
  map$excluded <- FALSE
  for (i in seq_len(nrow(exclusions))) {
    h <- !is.na(map$canonical) & grepl(exclusions$canonical_regex[i], canon_chr, perl = TRUE)
    map$excluded[h] <- TRUE
    map$reason[h]   <- exclusions$reason[i]
  }

  keepers <- map[map$action == "keep" & (keep_excluded | !map$excluded), , drop = FALSE]
  out <- df[, keys, drop = FALSE]
  if (!nrow(keepers)) return(list(data = out, map = map))

  for (cn in unique(keepers$canonical)) {
    grp <- keepers[keepers$canonical == cn, , drop = FALSE]
    policy <- grp$collapse[1]
    if (nrow(grp) == 1L) {
      v <- suppressWarnings(as.numeric(df[[grp$variable[1]]]))
    } else if (identical(policy, "none")) {
      warning(sprintf("[cov] %s: %d raw columns map to '%s' under collapse='none'",
                      country, nrow(grp), cn))
      v <- suppressWarnings(as.numeric(df[[grp$variable[1]]]))
    } else if (identical(policy, "mean")) {
      M <- vapply(grp$variable, function(v) suppressWarnings(as.numeric(df[[v]])),
                  numeric(nrow(df)))
      v <- rowMeans(M, na.rm = TRUE); v[is.nan(v)] <- NA_real_
    } else {  # latest
      pick <- grp$variable[order(grp$year, decreasing = TRUE, na.last = TRUE)][1]
      v <- suppressWarnings(as.numeric(df[[pick]]))
    }
    out[[cn]] <- v
  }

  # Physically impossible values, clamped at source. GHSL population comes back
  # with small negatives from raster resampling (Ghana min -7.29). Left alone
  # they reach models as real data, and benchmark_area_predictions() silently
  # zeroes negative weights, dropping that district from the aggregate entirely.
  for (cn in intersect(c("ghs_pop", "built_surface", "built_surface_nres",
                         "grassland_frac"), names(out))) {
    n_neg <- sum(is.finite(out[[cn]]) & out[[cn]] < 0)
    if (n_neg) {
      out[[cn]][is.finite(out[[cn]]) & out[[cn]] < 0] <- 0
      message(sprintf("[cov] %s: clamped %d negative value(s) in %s to 0",
                      country, n_neg, cn))
    }
  }

  uf <- cov_unit_factors(setdiff(names(out), keys), country, units)
  for (i in seq_len(nrow(uf))) {
    cn <- uf$canonical[i]
    if (uf$multiply[i] != 1 || uf$add[i] != 0)
      out[[cn]] <- out[[cn]] * uf$multiply[i] + uf$add[i]
  }
  map <- dplyr::left_join(map, uf, by = "canonical")
  list(data = out, map = map)
}
