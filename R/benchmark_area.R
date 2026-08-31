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
    # Match on the (Admin1, Admin2) PAIR when both sides carry Admin1. Matching
    # on the district NAME alone hands a duplicate-named district whichever of
    # the pair happens to come first -- six Malawi names occur in more than one
    # region -- so the benchmark reweights to the wrong population. Falls back
    # to the name for callers whose tables have no Admin1.
    if (all(c("Admin1", "Admin2") %in% names(pop)) &&
        "Admin1" %in% names(preds)) {
      w <- pop$pop[match(paste0(preds$Admin1, "\r", preds$Admin2),
                         paste0(pop$Admin1, "\r", pop$Admin2))]
    } else {
      w <- pop$pop[match(preds$Admin2, pop$Admin2)]
    }
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


# ═══════════════════════════════════════════════════════════════════════════
# ADMIN-1 BENCHMARKING
#
# The national benchmark above fixes one number for the whole country. But the
# surveys support REGIONAL estimates too, and far better than district ones:
# measured across the 24 country x outcome cells, the noise ceiling r_max is
# 0.664 at admin-1 against 0.132 at admin-2 (deff = 1.5; the ordering holds
# across the whole 1.0-2.0 deff band). Admin-1 direct estimates are therefore
# trustworthy anchors, and anchoring each region separately corrects regional
# level error that a single national shift cannot touch.
#
# WHY THE TARGETS ARE SHRUNK
# --------------------------
# A region is not automatically reliable just because it is bigger than a
# district. Ghana's regions carry a median of 22-64 biomarker reads for some
# outcomes, and treating such an estimate as an exact anchor imports its
# sampling error straight into the predictions. Each region's target is
# therefore shrunk toward the national estimate by an empirical-Bayes weight
#
#     lambda_r = v_between / (v_between + v_r)
#
# where v_r is that region's sampling variance and v_between is the estimated
# true between-region variance. A precisely measured region keeps its own
# value; a thin one is pulled toward the national figure. Setting shrink = FALSE
# reproduces hard benchmarking.
#
# Regions with too little data to estimate at all fall back to the national
# shift rather than being left unbenchmarked.
# ═══════════════════════════════════════════════════════════════════════════

#' Design-based prevalence per Admin-1 region.
#'
#' @param outcome_data output of build_outcome_dataset()
#' @param cc,oc country and outcome configs
#' @param admin1_col column holding the region; defaults to the config's, else "Admin1"
#' @return data.frame(Admin1, prev, n, samp_var) or NULL
admin1_design_based <- function(outcome_data, cc, oc, admin1_col = NULL) {
  d <- outcome_data$data
  a1 <- admin1_col %||% cc$admin1_col %||% "Admin1"
  if (!a1 %in% names(d)) return(NULL)
  if (!all(c(cc$weight_col, oc$binary) %in% names(d))) return(NULL)

  # Same uniform-outcome resolution the district target uses, so the benchmark
  # and the thing being benchmarked are the same quantity.
  derived <- resolve_uniform_outcome(d, cc, oc, label = "[benchmark-a1]")
  if (!is.null(derived)) d[[oc$binary]] <- derived

  y <- suppressWarnings(as.numeric(d[[oc$binary]]))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  g <- trimws(as.character(d[[a1]]))
  ok <- is.finite(y) & !is.na(g)
  if (sum(ok) < 30) return(NULL)
  y <- y[ok]; w <- w[ok]; g <- g[ok]

  regions <- sort(unique(g))
  prev <- vapply(regions, function(r) stats::weighted.mean(y[g == r], w[g == r]), numeric(1))
  n    <- vapply(regions, function(r) sum(g == r), numeric(1))
  data.frame(Admin1 = regions, prev = as.numeric(prev), n = as.numeric(n),
             samp_var = area_sampling_var(as.numeric(prev), as.numeric(n)),
             stringsAsFactors = FALSE)
}

#' Shrink regional targets toward the national estimate (empirical Bayes).
#'
#' @param a1 data.frame from admin1_design_based()
#' @param national national prevalence to shrink toward
#' @return `a1` with `target` and `lambda` columns added
shrink_admin1_targets <- function(a1, national) {
  v_obs <- stats::var(a1$prev)
  v_between <- max(0, v_obs - mean(a1$samp_var, na.rm = TRUE))
  lam <- if (v_between <= 0) rep(0, nrow(a1)) else
    v_between / (v_between + pmax(a1$samp_var, .Machine$double.eps))
  lam[!is.finite(lam)] <- 0
  a1$lambda <- lam
  a1$target <- lam * a1$prev + (1 - lam) * national
  a1
}

#' Benchmark an Admin-2 prediction table region by region.
#'
#' @param preds data.frame with Admin2 (and Admin1, or supply `admin1_map`)
#' @param pred_col prediction column to adjust
#' @param a1_targets data.frame(Admin1, target) -- typically from
#'   shrink_admin1_targets()
#' @param national fallback target for regions with no usable estimate
#' @param admin1_map optional data.frame(Admin2, Admin1) when `preds` lacks Admin1
#' @param pop optional data.frame(Admin2, pop) population weights
#' @param min_n minimum region sample size to benchmark on its own target
#' @return `preds` with the prediction column replaced; attribute "benchmark_a1"
#'   records the per-region shift actually applied
benchmark_admin2_to_admin1 <- function(preds, pred_col, a1_targets, national = NULL,
                                       admin1_map = NULL, pop = NULL,
                                       min_n = 25, method = "logit_shift") {
  if (is.null(preds) || !nrow(preds) || !pred_col %in% names(preds)) return(preds)
  if (is.null(a1_targets) || !nrow(a1_targets)) return(preds)

  if (!"Admin1" %in% names(preds)) {
    if (is.null(admin1_map) || !all(c("Admin2", "Admin1") %in% names(admin1_map))) {
      warning("benchmark_admin2_to_admin1: no Admin1 on preds and no admin1_map; ",
              "falling back to the national benchmark")
      return(if (is.null(national)) preds else
               benchmark_admin2_table(preds, pred_col, national, pop, method))
    }
    preds$Admin1 <- admin1_map$Admin1[match(preds$Admin2, admin1_map$Admin2)]
  }

  w_all <- NULL
  if (!is.null(pop) && all(c("Admin2", "pop") %in% names(pop)))
    # Pair match where both sides have Admin1 -- see benchmark_admin2_table().
    w_all <- if (all(c("Admin1", "Admin2") %in% names(pop)) &&
                 "Admin1" %in% names(preds))
      pop$pop[match(paste0(preds$Admin1, "\r", preds$Admin2),
                    paste0(pop$Admin1, "\r", pop$Admin2))]
    else pop$pop[match(preds$Admin2, pop$Admin2)]

  tgt <- a1_targets$target %||% a1_targets$prev
  names(tgt) <- a1_targets$Admin1
  nn <- a1_targets$n %||% rep(Inf, nrow(a1_targets)); names(nn) <- a1_targets$Admin1

  log <- list()
  for (r in unique(stats::na.omit(preds$Admin1))) {
    i <- which(preds$Admin1 == r)
    if (!length(i)) next
    t_r <- if (r %in% names(tgt) && is.finite(tgt[[r]]) && (nn[[r]] %||% 0) >= min_n)
      tgt[[r]] else national
    if (is.null(t_r) || !is.finite(t_r)) next
    b <- benchmark_area_predictions(preds[[pred_col]][i], t_r,
                                    if (is.null(w_all)) NULL else w_all[i],
                                    method = method)
    preds[[pred_col]][i] <- b$pred
    log[[length(log) + 1L]] <- data.frame(
      Admin1 = r, n_areas = length(i), target = t_r, delta = b$delta,
      before = b$before, after = b$after,
      used = if (r %in% names(tgt) && (nn[[r]] %||% 0) >= min_n) "region" else "national",
      stringsAsFactors = FALSE)
  }
  lg <- dplyr::bind_rows(log)
  attr(preds, "benchmark_a1") <- lg
  if (nrow(lg))
    cat(sprintf("[benchmark-a1] %s: %d region(s) anchored (%d on their own target), median |delta| = %.3f\n",
                pred_col, nrow(lg), sum(lg$used == "region"),
                stats::median(abs(lg$delta), na.rm = TRUE)))
  preds
}
