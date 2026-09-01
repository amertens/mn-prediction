# =============================================================================
# R/calibration_gate.R
#
# WS-D. A gate that refuses to ship a fitted cell whose national aggregate does
# not reproduce the survey's own national estimate.
#
# WHY THIS EXISTS
# ---------------
# WS2c measured the un-anchored population-weighted national aggregate against
# the design-based survey national estimate across 24 cells and found a mean
# absolute gap of 9.60 pp, reaching 77.57 pp for Sierra Leone child vitamin A,
# where the leave-one-region-out ridge predicts 89.6 percent against a measured
# 12.0 percent. A model that wrong at the national level should never reach a
# district map, a dashboard or a partner export, and nothing currently stops it.
#
# The gate is deliberately crude. It does not ask whether a model is good. It
# asks whether it reproduces a quantity the survey measures directly and well,
# and it is the national aggregate precisely because that is the one number this
# project can check without a ceiling argument.
#
# THE THRESHOLD
# -------------
# The larger of 10 percentage points and twice the survey's own CI half-width.
# The absolute floor catches catastrophic extrapolation; the CI-relative term
# stops the gate from failing a cell whose survey estimate is itself imprecise,
# which would punish a model for the survey's noise rather than its own error.
# Both are configurable and both are recorded on every row.
# =============================================================================

#' Evaluate the calibration gate for one cell.
#'
#' @param pred numeric district predictions
#' @param pop population weights aligned to `pred`; NULL for equal weights
#' @param national_prev design-based national prevalence for the cell
#' @param national_se design-based standard error of that estimate; NA is
#'   allowed and falls back to the absolute floor alone
#' @param abs_floor_pp absolute threshold in percentage points
#' @param ci_multiple multiple of the CI half-width
#' @return one-row data.frame with the aggregate, the gap and the verdict
calibration_gate <- function(pred, pop = NULL, national_prev, national_se = NA_real_,
                             abs_floor_pp = 10, ci_multiple = 2) {
  ok <- is.finite(pred)
  w <- if (is.null(pop)) rep(1, length(pred)) else as.numeric(pop)
  w[!is.finite(w) | w <= 0] <- NA_real_
  ok <- ok & is.finite(w)
  agg <- if (any(ok)) stats::weighted.mean(pred[ok], w[ok]) else NA_real_
  half <- if (is.finite(national_se)) 1.96 * national_se else NA_real_
  thr_pp <- max(abs_floor_pp,
                if (is.finite(half)) ci_multiple * 100 * half else -Inf)
  gap_pp <- 100 * (agg - national_prev)
  data.frame(
    national_pred_pp = round(100 * agg, 3),
    national_survey_pp = round(100 * national_prev, 3),
    gap_pp = round(gap_pp, 3),
    abs_gap_pp = round(abs(gap_pp), 3),
    ci_half_width_pp = if (is.finite(half)) round(100 * half, 3) else NA_real_,
    threshold_pp = round(thr_pp, 3),
    status = if (!is.finite(gap_pp)) "not_computed"
             else if (abs(gap_pp) > thr_pp) "calibration_failed" else "passed",
    n_districts = sum(ok),
    stringsAsFactors = FALSE)
}

#' Filter a deliverable table to the cells that passed.
#'
#' Deliverables and harness exports call this rather than reimplementing the
#' rule, so a cell cannot be excluded from one artefact and present in another.
#'
#' @param tbl a table with country and outcome columns
#' @param gate the gate report, with country, outcome and status
#' @param arm optional arm to read the verdict from, when the report holds several
#' @return `tbl` with failed cells removed, carrying an attribute naming them
apply_calibration_gate <- function(tbl, gate, arm = NULL) {
  if (is.null(tbl) || !nrow(tbl) || is.null(gate) || !nrow(gate)) return(tbl)
  g <- gate
  if (!is.null(arm) && "arm" %in% names(g)) g <- g[g$arm == arm, , drop = FALSE]
  failed <- g[g$status == "calibration_failed", , drop = FALSE]
  if (!nrow(failed)) { attr(tbl, "gate_excluded") <- character(0); return(tbl) }
  k <- function(d) paste(tolower(gsub("[^a-z]", "", tolower(d$country))), d$outcome)
  drop <- k(failed)
  keep <- !(k(tbl) %in% drop)
  out <- tbl[keep, , drop = FALSE]
  attr(out, "gate_excluded") <- unique(paste(failed$country, failed$outcome))
  out
}
