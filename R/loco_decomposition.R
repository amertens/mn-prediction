# =============================================================================
# R/loco_decomposition.R
#
# The cross-country estimator: an anchored level plus a covariate rank tilt.
#
# WHY THIS SHAPE, AND NOT AN ENSEMBLE
# -----------------------------------
# Leave-one-country-out does not fail the way leave-one-region-out fails.
# Measured over 16 cells in `relative_target_loco.csv`:
#
#   arm                                    mean r    MAE pp   |bias| pp
#   1 absolute (covariates, no anchor)      0.201     13.50       10.67
#   2 relative + TRUE national anchor       0.140     10.81        7.00
#   3 relative + BLIND (predicted) anchor   0.141     14.23       11.01
#   4 anchor only (flat national)              --      9.85        5.68
#
# Two things follow, and neither is about the learner.
#
#   (a) Taking the held-out country's national prevalence and applying it FLAT
#       to every district beats every covariate model on error. The covariates
#       buy about 0.20 of correlation and cost 3.7 pp of accuracy.
#   (b) Arms 2 and 3 differ ONLY in whether the anchor is known or predicted,
#       and they differ by 3.4 pp. That gap is the whole LOCO problem.
#
# So the estimator separates the two quantities that fail separately: a LEVEL,
# which must come from an anchor, and a RANKING, which covariates supply
# partially and unevenly.
#
# WHY THE WEIGHTS ARE FIXED AND NOT LEARNED
# -----------------------------------------
# Within a country there are 14 to 87 districts to fit a stack on. Across
# countries there are FOUR, so a leave-one-country-out stack fits its weights on
# three points. Transportability is not estimable at that n, and a stack that
# pretended otherwise would be fitting noise and calling it adaptivity. The tilt
# weight is therefore a declared parameter, its default is conservative, and the
# sweep over it is reported rather than optimised away.
#
# Per-country transport, from the same table, is why no single weight is right:
#
#   Gambia 0.432   Ghana 0.374   Malawi 0.095   Sierra Leone -0.098
# =============================================================================

#' Within-country standardisation of a covariate prediction.
#'
#' Strips the LEVEL out of a cross-country prediction and keeps only the shape.
#' This is the step that makes a transported model usable at all: the raw
#' prediction carries the training countries' level, which is wrong by up to
#' 42 pp for the held-out country.
#'
#' @param pred numeric predictions for one country's districts
#' @return z-scores, or zeros where the prediction has no usable spread
loco_rank_signal <- function(pred) {
  ok <- is.finite(pred)
  if (sum(ok) < 3L) return(rep(0, length(pred)))
  s <- stats::sd(pred[ok])
  if (!is.finite(s) || s <= 0) return(rep(0, length(pred)))
  z <- rep(0, length(pred))
  z[ok] <- (pred[ok] - mean(pred[ok])) / s
  z
}

#' Anchor plus covariate tilt.
#'
#' @param anchor  the level. Either one number for the whole country, or a
#'   vector giving each district's anchor unit (region-level anchoring).
#' @param pred    the transported covariate prediction, any level
#' @param spread  assumed within-country SD of true district prevalence, on the
#'   probability scale. Taken from the training countries when not supplied.
#' @param tilt    how much of the covariate ranking to apply, 0 to 1. The
#'   DEFAULT IS DELIBERATELY LOW. At tilt = 0 this is the `anchor only` arm,
#'   which is the best-scoring arm in the measured table.
#' @return numeric vector of district predictions
loco_anchor_tilt <- function(anchor, pred, spread, tilt = 0.35) {
  z <- loco_rank_signal(pred)
  a <- if (length(anchor) == 1L) rep(anchor, length(pred)) else anchor
  out <- a + tilt * spread * z
  pmin(pmax(out, 0), 1)
}

#' The within-country spread of district prevalence, learned from the countries
#' that DO have labels.
#'
#' A transported ranking has no scale of its own. Borrowing the spread from the
#' training countries is an assumption, and it is the assumption most likely to
#' be wrong for a country whose districts are more or less unequal than the
#' training set's. It is reported alongside the estimate for that reason.
loco_training_spread <- function(p, country, n = NULL) {
  sp <- vapply(split(p, country), function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 4L) return(NA_real_)
    stats::sd(v)
  }, numeric(1))
  sp <- sp[is.finite(sp)]
  if (!length(sp)) return(NA_real_)
  stats::median(sp)
}

#' Score one held-out country under the decomposition, across a tilt sweep.
#'
#' Returns one row per tilt value so the trade-off is visible rather than
#' collapsed to a single chosen setting: the tilt buys ranking and costs error,
#' and which side of that a user wants depends on whether they are targeting
#' districts or estimating burden.
#'
#' @param truth   observed district prevalence in the held-out country
#' @param pred    transported covariate prediction for those districts
#' @param anchor_true  the held-out country's actual level (the ceiling case)
#' @param anchor_blind a predicted level, e.g. from VMNIS (the realistic case)
#' @param spread  within-country SD from the training countries
#' @param tilts   grid to sweep
loco_score_decomposition <- function(truth, pred, anchor_true, anchor_blind,
                                     spread, tilts = c(0, 0.2, 0.35, 0.5, 0.75, 1)) {
  ok <- is.finite(truth)
  rows <- list()
  for (nm in c("true", "blind")) {
    a <- if (nm == "true") anchor_true else anchor_blind
    if (!is.finite(a[1]) && length(a) == 1L) next
    for (tl in tilts) {
      est <- loco_anchor_tilt(a, pred, spread, tilt = tl)
      o <- ok & is.finite(est)
      if (sum(o) < 4L) next
      r <- if (stats::sd(est[o]) > 0 && stats::sd(truth[o]) > 0)
        suppressWarnings(stats::cor(est[o], truth[o])) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        anchor = nm, tilt = tl, n = sum(o),
        r = r,
        mae_pp = 100 * mean(abs(est[o] - truth[o])),
        bias_pp = 100 * mean(est[o] - truth[o]),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
