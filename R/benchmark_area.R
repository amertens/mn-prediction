# =============================================================================
# R/benchmark_area.R
#
# Benchmark Admin-2 predictions to a design-based national (or Admin-1) total.
#
# WHY
# ---
# These surveys are designed and weighted to deliver NATIONAL estimates. Admin-2
# is an unplanned domain: whatever sample lands in a district is incidental, the
# weights are not calibrated at that level, and the direct estimates are both
# noisy and potentially biased within-domain.
#
# The practical consequence is that the model should DISAGGREGATE a national
# figure the survey already delivers well, not re-estimate it. Left unbenchmarked,
# the area model's own population-weighted aggregate drifts away from the
# design-based national estimate -- and in the cross-country (LOCO) setting it
# misses by ~5 pp on average and by a factor of 2-3 on iron
# (sandbox_parsimony/FINDINGS.md section 19).
#
# HOW
# ---
# A single shift on the LOGIT scale, chosen so the population-weighted aggregate
# of the adjusted predictions equals the target:
#
#   find delta such that  sum_i w_i * plogis(logit(p_i) + delta) / sum_i w_i = P
#
# Chosen over the two textbook alternatives because:
#   * ratio benchmarking (p_i * P / Pbar) can push predictions above 1;
#   * difference benchmarking (p_i + P - Pbar) can push them below 0, and
#     distorts low-prevalence districts hardest;
#   * the logit shift is strictly monotone, so the DISTRICT RANKING is
#     preserved exactly. That matters here: the ranking is what the model is
#     comparatively good at, and the level is what it is bad at, so the
#     correction should touch only the level.
#
# Both alternatives are available via `method` for comparison.
# =============================================================================

#' Population-weighted aggregate of a set of area predictions
#'
#' @param pred numeric predictions in [0, 1]
#' @param w population weights; equal weights if NULL (rarely what you want --
#'   an unweighted mean over districts is a mean over POLYGONS, not people)
area_aggregate <- function(pred, w = NULL) {
  ok <- is.finite(pred)
  if (!any(ok)) return(NA_real_)
  if (is.null(w)) w <- rep(1, length(pred))
  w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
  if (sum(w[ok]) <= 0) w[ok] <- 1
  stats::weighted.mean(pred[ok], w[ok])
}

#' Benchmark area predictions to a known aggregate.
#'
#' @param pred numeric vector of predicted prevalences in [0, 1]
#' @param target the design-based aggregate to reproduce (a proportion)
#' @param w population weights, same length as `pred`
#' @param method "logit_shift" (default), "ratio" or "difference"
#' @param eps clamp applied before the logit so 0 and 1 stay finite
#' @return list(pred, delta, before, after, method); `pred` is the adjusted
#'   vector, `delta` the shift applied (NA for non-shift methods)
benchmark_area_predictions <- function(pred, target, w = NULL,
                                       method = c("logit_shift", "ratio", "difference"),
                                       eps = 1e-4) {
  method <- match.arg(method)
  before <- area_aggregate(pred, w)
  if (!is.finite(target) || !is.finite(before)) {
    return(list(pred = pred, delta = NA_real_, before = before,
                after = before, method = method))
  }

  if (method == "ratio") {
    out <- pmin(pmax(pred * (target / max(before, 1e-8)), 0), 1)
  } else if (method == "difference") {
    out <- pmin(pmax(pred + (target - before), 0), 1)
  } else {
    p <- pmin(pmax(pred, eps), 1 - eps)
    lg <- stats::qlogis(p)
    f <- function(d) area_aggregate(stats::plogis(lg + d), w) - target
    # target is a proportion, so a shift of +/- 20 on the logit scale brackets
    # anything achievable; fall back to the unadjusted vector if it does not.
    rt <- tryCatch(stats::uniroot(f, c(-20, 20), tol = 1e-10),
                   error = function(e) NULL)
    if (is.null(rt))
      return(list(pred = pred, delta = NA_real_, before = before,
                  after = before, method = method))
    out <- pmin(pmax(stats::plogis(lg + rt$root), 0), 1)
    return(list(pred = out, delta = rt$root, before = before,
                after = area_aggregate(out, w), method = method))
  }
  list(pred = out, delta = NA_real_, before = before,
       after = area_aggregate(out, w), method = method)
}

#' Design-based national prevalence from the individual-level survey data.
#'
#' Deliberately ignores district entirely: this is the survey-weighted mean over
#' all respondents, which is the estimate the sampling design supports.
#'
#' @param outcome_data output of build_outcome_dataset()
#' @param cc country config
#' @param oc outcome config
#' @return list(prev, n) or NULL
national_design_based <- function(outcome_data, cc, oc) {
  d <- outcome_data$data
  if (!all(c(cc$weight_col, oc$binary) %in% colnames(d))) return(NULL)
  y <- suppressWarnings(as.numeric(d[[oc$binary]]))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  ok <- is.finite(y)
  if (sum(ok) < 30) return(NULL)
  list(prev = stats::weighted.mean(y[ok], w[ok]), n = sum(ok))
}

#' Benchmark an Admin-2 prediction table in place.
#'
#' @param preds data.frame with Admin2 and a prediction column
#' @param pred_col name of the prediction column
#' @param target design-based aggregate to hit
#' @param pop optional data.frame(Admin2, pop) of population weights. When
#'   absent, districts are weighted equally and the result is flagged, because
#'   an unweighted aggregate answers a different question.
#' @return `preds` with the prediction column replaced and attributes recording
#'   what was done
benchmark_admin2_table <- function(preds, pred_col, target, pop = NULL,
                                   method = "logit_shift") {
  if (is.null(preds) || !nrow(preds) || !pred_col %in% names(preds)) return(preds)
  w <- NULL
  if (!is.null(pop) && all(c("Admin2", "pop") %in% names(pop))) {
    w <- pop$pop[match(preds$Admin2, pop$Admin2)]
    if (all(!is.finite(w))) w <- NULL
  }
  if (is.null(w))
    warning("benchmark_admin2_table: no population weights supplied; ",
            "benchmarking to an UNWEIGHTED district mean, which is not the ",
            "same quantity as the design-based national prevalence")
  b <- benchmark_area_predictions(preds[[pred_col]], target, w, method = method)
  preds[[pred_col]] <- b$pred
  attr(preds, "benchmark") <- b[c("delta", "before", "after", "method")]
  cat(sprintf("[benchmark] %s: aggregate %.1f%% -> %.1f%% (target %.1f%%, delta %.3f)\n",
              pred_col, 100 * b$before, 100 * b$after, 100 * target,
              b$delta %||% NA_real_))
  preds
}
