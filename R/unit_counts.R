# =============================================================================
# R/unit_counts.R
#
# HOW MANY ADMIN-2 UNITS SHOULD A TABLE HAVE, AND WHAT DOES IT MEAN WHEN IT HAS
# MORE?
#
# WHY THIS EXISTS
# ---------------
# The pair-key defect (Section 8 of docs/SESSION_FINDINGS_FOR_REVIEW.md) is
# silent in the code and loud in the row count. Joining Malawi's Admin-2 tables
# on the district name fans rows, and the tenth instance was caught only because
# a deployed dashboard listed 90 districts where the prediction table held 87.
# A row count is therefore the cheapest reliable detector the project has.
#
# THE ASYMMETRY THAT MAKES THIS A USABLE TEST
# -------------------------------------------
# A result table may legitimately hold FEWER units than the survey supports: a
# cell drops districts with no outcome data, aggregation drops units under a
# minimum sample size, and a model may fail to predict a region. It may never
# hold MORE. Every mechanism that raises the count above the analytic number is
# a defect, and a name-keyed join is the one this project keeps meeting. The
# assertion is therefore one-sided: observed <= reference.
#
# THE REFERENCE IS DERIVED FROM DATA, NOT DECLARED
# ------------------------------------------------
# Hard-coding "Malawi 87" would make the test agree with whatever the author
# believed. The analytic count is read from the survey Admin-2 targets in the
# built store, and the GADM counts from the boundary layer, so the test tracks
# the data and fails when a country's units genuinely change.
# =============================================================================

#' Analytic Admin-2 unit counts per country, read from the built targets store.
#'
#' The analytic count is the largest number of Admin-2 rows any of a country's
#' `svy_admin2_*` targets holds. Individual outcomes drop districts where the
#' outcome was not measured, so the maximum over outcomes is the count a table
#' covering that country may not exceed.
#'
#' @param store path to the targets store
#' @param countries named list of country configs; defaults to get_country_configs()
#' @return data.frame(country, analytic_units, n_outcomes, per_outcome)
admin2_analytic_counts <- function(store = here::here("_targets_full"),
                                   countries = NULL) {
  cfgs <- countries %||% get_country_configs()
  rows <- list()
  for (cn in names(cfgs)) {
    cc <- cfgs[[cn]]
    per <- integer(0)
    for (ocn in names(cc$outcomes)) {
      nm <- paste0("svy_admin2_", tolower(cn), "_", ocn)
      sv <- tryCatch(targets::tar_read_raw(nm, store = store), error = function(e) NULL)
      if (is.null(sv) || !nrow(sv)) next
      per[ocn] <- nrow(sv)
    }
    if (!length(per)) next
    rows[[cn]] <- data.frame(
      country        = cc$country,
      analytic_units = max(per),
      n_outcomes     = length(per),
      per_outcome    = paste(sprintf("%s=%d", names(per), per), collapse = ";"),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' GADM level-2 polygon and distinct-name counts per country.
#'
#' The gap between these two numbers is the size of the fan a name-only join can
#' produce. Malawi is the country where it is not zero.
#'
#' @param polys named list or data.frame supplying Admin1/Admin2 per country
#' @return data.frame(country, gadm_polygons, gadm_unique_names, dup_names)
admin2_gadm_counts <- function(polys) {
  if (is.null(polys)) return(NULL)
  rows <- list()
  for (cn in names(polys)) {
    p <- polys[[cn]]
    if (is.null(p) || !nrow(p) || !"Admin2" %in% names(p)) next
    a2 <- as.character(p$Admin2)
    tb <- table(a2)
    rows[[cn]] <- data.frame(
      country           = cn,
      gadm_polygons     = length(a2),
      gadm_unique_names = length(unique(a2)),
      dup_names         = sum(tb > 1),
      dup_name_list     = paste(names(tb)[tb > 1], collapse = ";"),
      stringsAsFactors  = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Check a result table's per-country unit counts against the reference.
#'
#' @param tbl a result table
#' @param count_col name of its unit-count column (e.g. "n_admin2", "n_units")
#' @param reference data.frame(country, analytic_units)
#' @param country_col name of its country column
#' @param label a name for the table, used in the returned rows
#' @return data.frame(label, country, observed_max, reference, status) where
#'   status is "ok" when observed_max <= reference and "OVER" otherwise.
check_unit_counts <- function(tbl, count_col, reference,
                              country_col = "country", label = "table") {
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  if (!all(c(count_col, country_col) %in% names(tbl))) return(NULL)
  # Country labels differ across tables ("SierraLeone", "Sierra Leone",
  # "sierraleone"). Compare on a squashed, case-folded key so a naming
  # difference does not read as a missing country.
  norm <- function(x) tolower(gsub("[^a-z]", "", tolower(as.character(x))))
  ref <- reference; ref$.k <- norm(ref$country)
  obs <- tapply(suppressWarnings(as.numeric(tbl[[count_col]])),
                norm(tbl[[country_col]]), max, na.rm = TRUE)
  rows <- list()
  for (k in names(obs)) {
    r <- ref$analytic_units[match(k, ref$.k)]
    rows[[k]] <- data.frame(
      label = label, country = k,
      observed_max = as.integer(obs[[k]]),
      reference = if (is.na(r)) NA_integer_ else as.integer(r),
      status = if (is.na(r)) "no reference"
               else if (obs[[k]] > r) "OVER" else "ok",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
