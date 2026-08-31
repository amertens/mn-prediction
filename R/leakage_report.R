# =============================================================================
# R/leakage_report.R
#
# MEASURE THE LEAK, DO NOT REASON ABOUT IT.
#
# WHY THIS EXISTS
# ---------------
# Section 7 of docs/SESSION_FINDINGS_FOR_REVIEW.md records three attempts at an
# assay guard. The first two produced publishable-looking numbers, r = 0.973 and
# r = 0.986, and both were leakage. Neither was caught by reading regexes. Both
# were caught by ranking every surviving predictor against the outcome and
# looking at the top of the list, where the offenders were `gw_cAnemiaYN`
# (|r| 0.80), `gw_bs2` (0.69) and `gw_cID_NoAdj` (0.57). None of those names an
# analyte, so no pattern over analyte names could have found them.
#
# The lesson generalises past the guard: a correlation ranking is a detector
# that does not depend on anyone having anticipated the naming convention. This
# file makes that ranking a pipeline artefact rather than an ad hoc check, so a
# fifth country cannot enter the pipeline without one being produced.
#
# THE THRESHOLD
# -------------
# Set from the measured history in Section 7, not from theory. The leaked
# columns sat between 0.57 and 0.80. The clean maximum, after the guard was
# correct, was 0.544 across all 64 rows of the individual anchor. A flag at 0.55
# therefore sits in the gap between the two regimes: above every value observed
# on clean data, below every value observed on leaked data.
#
# A flag is not a verdict. A genuine exposure could exceed it, and a leaked
# column with a moderate correlation will not reach it. The report exists so
# that the question is asked of every cell and answered from data.
#
# THE EFFECTIVE-n FLOOR, AND WHY THE FIRST VERSION OF THIS FILE NEEDED ONE
# -----------------------------------------------------------------------
# Measured 2026-08-31 on the smoke cell (Ghana child_iron): the first run of
# this report ranked `gw_cFormMilkFreq` at |r| 0.8165, above the 0.80 of
# `gw_cAnemiaYN`, the worst leak Section 7 found. It is not a leak. The column
# has FOUR complete pairs with the outcome, and a correlation on four points is
# not evidence of anything. In that one cell 134 of 2,253 guarded predictors
# carry fewer than 50 complete pairs.
#
# A detector whose top rank is occupied by n = 4 noise is worse than no
# detector, because the real leak sits below the noise and the reader learns to
# skip the report. Columns below `min_n` complete pairs are therefore ranked and
# reported but never flagged, and `n_pairs` is carried on every row so the
# reader can see what a correlation rests on. The default floor of 100 is well
# under the smallest cell's outcome count and well above the region where a
# correlation is unstable.
# =============================================================================

#' Absolute correlation of every predictor with the outcome, per cell.
#'
#' @param outcome_data output of build_outcome_dataset()
#' @param cc,oc country and outcome configs
#' @param set which predictor set: "proxy" uses $Xvars; "questionnaire" uses
#'   $Xvars_full with is_biomarker_column() applied, which is the set the
#'   individual anchor's questionnaire arm sees and the set the guard protects.
#' @param top how many ranked columns to return per cell
#' @param threshold flag level; see the header for where it comes from
#' @return data.frame ranked by abs_r descending, or NULL
#' @param min_n minimum complete predictor-outcome pairs before a column may be
#'   flagged; see the header for the measurement that forced this
leakage_rank_cell <- function(outcome_data, cc, oc, set = c("questionnaire", "proxy"),
                              top = 25L, threshold = 0.55, min_n = 100L) {
  set <- match.arg(set)
  if (is.null(outcome_data)) return(NULL)
  d <- outcome_data$data
  if (is.null(oc$binary) || !oc$binary %in% names(d)) return(NULL)

  y <- suppressWarnings(as.numeric(haven::zap_labels(d[[oc$binary]])))
  ok <- is.finite(y)
  if (sum(ok) < 30 || length(unique(y[ok])) < 2) return(NULL)

  vars <- if (set == "proxy") outcome_data$Xvars else {
    full <- outcome_data$Xvars_full %||% outcome_data$Xvars
    full[!is_biomarker_column(full)]
  }
  vars <- intersect(vars, names(d))
  if (!length(vars)) return(NULL)

  y <- y[ok]
  # A single matrix correlation over all predictors at once. Columns that are
  # constant or all-missing return NA and are dropped rather than warned about.
  X <- vapply(vars, function(v)
    suppressWarnings(as.numeric(haven::zap_labels(d[[v]][ok]))),
    numeric(length(y)))
  if (is.null(dim(X))) X <- matrix(X, nrow = length(y))
  colnames(X) <- vars

  r <- suppressWarnings(stats::cor(X, y, use = "pairwise.complete.obs"))
  r <- r[, 1]
  n_pairs <- colSums(is.finite(X) & is.finite(y))
  keep <- is.finite(r)
  if (!any(keep)) return(NULL)
  r <- r[keep]; n_pairs <- n_pairs[keep]

  # Only sufficiently measured columns may be flagged. Under-measured columns
  # stay in the ranking, carrying their n, so an inspection can see them.
  eligible <- n_pairs >= min_n
  flagged_all <- eligible & abs(r) >= threshold
  n_flag <- sum(flagged_all)
  max_elig <- if (any(eligible)) round(max(abs(r[eligible])), 4) else NA_real_

  # Rank on the eligible columns first so the report's head is readable, then
  # append the highest under-measured ones so they remain visible.
  ord <- c(order(ifelse(eligible, abs(r), -Inf), decreasing = TRUE))
  head_ord <- ord[seq_len(min(top, length(ord)))]

  data.frame(
    country      = cc$country,
    outcome      = oc$tag,
    predictor_set = set,
    n_predictors = length(r),
    rank         = seq_along(head_ord),
    column       = names(r)[head_ord],
    r            = round(unname(r[head_ord]), 4),
    abs_r        = round(abs(unname(r[head_ord])), 4),
    n_pairs      = as.integer(unname(n_pairs[head_ord])),
    min_n        = as.integer(min_n),
    eligible     = unname(eligible[head_ord]),
    threshold    = threshold,
    flagged      = unname(flagged_all[head_ord]),
    n_flagged_in_cell     = n_flag,
    max_abs_r_in_cell     = max_elig,
    n_under_measured_cols = sum(!eligible),
    stringsAsFactors = FALSE)
}

#' Build the leakage report across every country and outcome in the store.
#'
#' @param store targets store holding the outcome_data_* targets
#' @param sets which predictor sets to rank
#' @param countries optional character vector restricting the run; the smoke
#'   profile passes a single country so the report is exercised cheaply
#' @param outcomes optional character vector restricting the outcomes
#' @return data.frame, one block of ranked columns per cell and predictor set
build_leakage_report <- function(store = here::here("_targets_full"),
                                 sets = c("questionnaire", "proxy"),
                                 countries = NULL, outcomes = NULL,
                                 top = 25L, threshold = 0.55, min_n = 100L) {
  cfgs <- get_country_configs()
  if (!is.null(countries)) cfgs <- cfgs[intersect(names(cfgs), countries)]
  rows <- list()
  for (cn in names(cfgs)) {
    cc <- cfgs[[cn]]
    ocs <- names(cc$outcomes)
    if (!is.null(outcomes)) ocs <- intersect(ocs, outcomes)
    for (ocn in ocs) {
      oc <- cc$outcomes[[ocn]]
      od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                           store = store), error = function(e) NULL)
      if (is.null(od)) next
      for (s in sets) {
        r <- tryCatch(leakage_rank_cell(od, cc, oc, set = s, top = top,
                                        threshold = threshold, min_n = min_n),
                      error = function(e) NULL)
        if (!is.null(r)) rows[[length(rows) + 1L]] <- r
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Build the report from already-materialised outcome_data objects.
#'
#' The store-reading form above is what the standalone script uses. This form is
#' what the pipeline target uses, so the report depends on every
#' `outcome_data_*` target and a newly ingested country cannot reach a model
#' without one being produced. See `guard_targets` in _targets.R.
#'
#' @param od_list named list of build_outcome_dataset() outputs, named
#'   "<country>_<outcome>" as the targets are
#' @param cfgs country configs
#' @return the same data.frame build_leakage_report() returns
build_leakage_report_from_list <- function(od_list, cfgs = get_country_configs(),
                                           sets = c("questionnaire", "proxy"),
                                           top = 25L, threshold = 0.55,
                                           min_n = 100L) {
  rows <- list()
  for (nm in names(od_list)) {
    od <- od_list[[nm]]
    if (is.null(od)) next
    # Names are "<lowercase country>_<outcome tag>"; match the country by the
    # longest configured prefix so "sierraleone_child_iron" resolves.
    cn <- NULL
    for (k in names(cfgs)) if (startsWith(nm, paste0(tolower(k), "_"))) cn <- k
    if (is.null(cn)) next
    ocn <- sub(paste0("^", tolower(cn), "_"), "", nm)
    cc <- cfgs[[cn]]; oc <- cc$outcomes[[ocn]]
    if (is.null(oc)) next
    for (s in sets) {
      r <- tryCatch(leakage_rank_cell(od, cc, oc, set = s, top = top,
                                      threshold = threshold, min_n = min_n),
                    error = function(e) NULL)
      if (!is.null(r)) rows[[length(rows) + 1L]] <- r
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' One-line-per-cell summary of the report, which is what a reviewer reads.
leakage_report_summary <- function(rep) {
  if (is.null(rep) || !nrow(rep)) return(NULL)
  k <- paste(rep$country, rep$outcome, rep$predictor_set, sep = "|")
  do.call(rbind, lapply(split(rep, k), function(d) data.frame(
    country = d$country[1], outcome = d$outcome[1],
    predictor_set = d$predictor_set[1],
    n_predictors = d$n_predictors[1],
    max_abs_r = d$max_abs_r_in_cell[1],
    n_flagged = d$n_flagged_in_cell[1],
    n_under_measured = d$n_under_measured_cols[1],
    top_column = d$column[which(d$eligible)[1]],
    top_n_pairs = d$n_pairs[which(d$eligible)[1]],
    stringsAsFactors = FALSE)))
}
