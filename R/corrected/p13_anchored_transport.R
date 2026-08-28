# =============================================================================
# R/corrected/p13_anchored_transport.R  —  P13: anchored transport, honest
# new-country intervals, and exceedance probabilities
#
# WHY
# ---
# WS3 established that the transported map's remaining error is dominated by a
# country-level height the covariates cannot supply: on the log biomarker scale
# the pure country intercept is the larger term for three of four outcomes, and
# for child iron the training pool sits between 0.526 and 2.679 times the
# held-out country's level. A monotone shift that forces the map's aggregate to
# a known national value is the direct remedy for exactly that error, and it
# leaves district RANKING untouched, so it cannot manufacture spatial signal.
#
# THREE PIECES
# ------------
# 1. Anchored LOCO. Fit the area transport recipe on the training countries,
#    predict the held-out country, then shift on the logit scale so the
#    aggregate matches an anchor. Two anchor sources:
#
#      own_survey  the held-out country's own design-based national value. Not
#                  available in deployment; it measures the VALUE of anchoring,
#                  the best a perfect anchor could do.
#      external    a published national estimate (metadata/anchors/). Measures
#                  anchor QUALITY, which is the deployment-relevant question.
#
#    Coverage of the external arm is the binding constraint and is reported, not
#    worked around: see scripts/build_anchors.R.
#
# 2. New-country prediction interval. A country random-intercept model on the
#    continuous scale, giving the interval a genuinely new country's national
#    level would fall in. With four countries this interval is wide. That width
#    is the honest product, not a defect to be tuned away.
#
# 3. Exceedance probabilities. P(district prevalence exceeds a public-health
#    threshold), from the anchored prediction and its interval.
#
# WEIGHTS
# -------
# Population weights are not available: WorldPop ingestion is a separate
# workstream that was not run. Aggregation uses survey n as the weight, which is
# NOT the same quantity as a population-weighted national prevalence. Every row
# records this in `aggregation_weight` so the figures cannot be read as
# population-weighted.
# =============================================================================

ANCHOR_SOURCES <- c("unanchored", "own_survey", "external")

#' Anchor one country's predicted Admin-2 map to a national value.
#'
#' Delegates the shift to benchmark_area_predictions() (R/benchmark_area.R:62),
#' which is the same monotone logit-scale shift the project already uses, so
#' this is a new application of an existing estimator rather than a second
#' implementation of it.
#'
#' @param pred numeric predicted prevalences
#' @param target the national value to hit
#' @param w weights for the aggregate (survey n here; see file header)
#' @return list(pred, delta, before, after)
anchor_predictions <- function(pred, target, w = NULL) {
  if (!is.finite(target)) return(list(pred = pred, delta = NA_real_,
                                      before = NA_real_, after = NA_real_))
  b <- benchmark_area_predictions(pred, target, w, method = "logit_shift")
  list(pred = b$pred, delta = b$delta %||% NA_real_,
       before = b$before, after = b$after)
}

#' Metrics for one anchored or unanchored map.
.at_metrics <- function(pred, truth, w) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 4) return(NULL)
  data.frame(
    n_areas     = sum(ok),
    pearson_r   = suppressWarnings(stats::cor(pred[ok], truth[ok])),
    spearman_r  = suppressWarnings(stats::cor(pred[ok], truth[ok], method = "spearman")),
    mae_pp      = 100 * mean(abs(pred[ok] - truth[ok])),
    rmse_pp     = 100 * sqrt(mean((pred[ok] - truth[ok])^2)),
    nat_bias_pp = 100 * (stats::weighted.mean(pred[ok], w[ok]) -
                         stats::weighted.mean(truth[ok], w[ok])),
    stringsAsFactors = FALSE
  )
}

#' Run the anchored-transport LOCO comparison for one outcome.
#'
#' @param pool assemble/build_area_loco_dataset()-shaped list(data, predictors, countries)
#' @param outcome outcome tag
#' @param anchors metadata/anchors/national_anchors.csv, already read
#' @param recipe area transport recipe
#' @return data.frame, one row per held-out country per anchor source
run_anchored_loco <- function(pool, outcome, anchors, recipe = AREA_TRANSPORT_RECIPE) {
  res <- run_area_transport_loco(pool, recipe)
  if (is.null(res) || is.null(res$predictions) || !nrow(res$predictions)) return(NULL)
  p <- res$predictions

  ext <- anchors[anchors$outcome == outcome &
                 anchors$usable_as_anchor %in% TRUE, , drop = FALSE]

  rows <- list()
  for (ho in unique(p$country)) {
    d <- p[p$country == ho, , drop = FALSE]
    w <- d$n_svy; w[!is.finite(w) | w <= 0] <- 1
    truth <- d$survey_prev

    # own_survey: the held-out country's design-based national value. Computed
    # from its own survey rows, weighted the same way the aggregate is, so the
    # anchor and the aggregate are the same functional.
    own <- stats::weighted.mean(truth, w)

    e <- ext[ext$country == ho, , drop = FALSE]
    external <- if (nrow(e)) e$anchor_prevalence[which.min(e$year_gap)] else NA_real_
    ext_year_gap <- if (nrow(e)) e$year_gap[which.min(e$year_gap)] else NA_integer_

    for (src in ANCHOR_SOURCES) {
      target <- switch(src, unanchored = NA_real_, own_survey = own, external = external)
      if (src != "unanchored" && !is.finite(target)) next
      a <- if (src == "unanchored") list(pred = d$modeled_prev, delta = NA_real_,
                                         before = NA_real_, after = NA_real_)
           else anchor_predictions(d$modeled_prev, target, w)
      m <- .at_metrics(a$pred, truth, w)
      if (is.null(m)) next
      m$outcome <- outcome; m$held_out <- ho; m$anchor_source <- src
      m$anchor_target <- target
      m$anchor_delta_logit <- a$delta
      m$aggregate_before <- a$before; m$aggregate_after <- a$after
      m$external_year_gap <- if (src == "external") ext_year_gap else NA_integer_
      m$aggregation_weight <- "survey n (not population)"
      rows[[length(rows) + 1L]] <- m
    }
  }
  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)

  # A monotone shift cannot change rank order. If Spearman moves, the shift was
  # not monotone and the result is wrong; this asserts rather than assumes.
  chk <- out |>
    dplyr::group_by(.data$outcome, .data$held_out) |>
    dplyr::summarise(spread = diff(range(.data$spearman_r, na.rm = TRUE)),
                     .groups = "drop")
  bad <- chk[is.finite(chk$spread) & chk$spread > 1e-8, , drop = FALSE]
  if (nrow(bad))
    warning("anchoring changed the district ranking for ", nrow(bad),
            " fold(s); the shift is supposed to be monotone", call. = FALSE)
  out$ranking_preserved <- nrow(bad) == 0
  out
}

# ── New-country prediction interval ──────────────────────────────────────────

#' Prediction interval for the national level of a country not in the sample.
#'
#' Two estimators, reported side by side because with four countries they can
#' disagree and neither is obviously right:
#'
#'   `t_interval`  the textbook prediction interval for a new observation from a
#'                 sample of country means: mean +/- t(n-1) * s * sqrt(1 + 1/n).
#'                 Makes no distributional assumption beyond normality of the
#'                 country means and does not try to separate within- from
#'                 between-country variance.
#'   `lmer`        a country random-intercept model on the AREA-level response,
#'                 which does separate them: the interval uses the estimated
#'                 between-country SD plus the uncertainty in the grand mean.
#'
#' With four countries the between-country variance is estimated from three
#' degrees of freedom, so both intervals are wide and the lmer one can collapse
#' to a boundary estimate of zero. That is reported, not hidden.
#'
#' @param d data.frame with country and a numeric response column
#' @param response column name
#' @return data.frame with one row per estimator
new_country_interval <- function(d, response = "area_log_biomarker") {
  d <- d[is.finite(d[[response]]), , drop = FALSE]
  if (!nrow(d) || dplyr::n_distinct(d$country) < 3) return(NULL)
  cm <- tapply(d[[response]], d$country, mean, na.rm = TRUE)
  n  <- length(cm); m <- mean(cm); s <- stats::sd(cm)
  tq <- stats::qt(0.975, df = n - 1)

  rows <- list(data.frame(
    estimator = "t_interval", n_countries = n,
    centre = m, sd_between = s,
    lo = m - tq * s * sqrt(1 + 1 / n), hi = m + tq * s * sqrt(1 + 1 / n),
    boundary_fit = FALSE, stringsAsFactors = FALSE))

  if (requireNamespace("lme4", quietly = TRUE)) {
    fit <- tryCatch(
      lme4::lmer(stats::as.formula(paste(response, "~ 1 + (1 | country)")),
                 data = d, REML = TRUE),
      error = function(e) NULL)
    if (!is.null(fit)) {
      vc  <- as.data.frame(lme4::VarCorr(fit))
      tau <- sqrt(vc$vcov[vc$grp == "country"][1])
      fe  <- unname(lme4::fixef(fit)[1])
      se  <- sqrt(as.numeric(stats::vcov(fit))[1])
      wid <- stats::qt(0.975, df = n - 1) * sqrt(tau^2 + se^2)
      rows[[2]] <- data.frame(
        estimator = "lmer_random_intercept", n_countries = n,
        centre = fe, sd_between = tau,
        lo = fe - wid, hi = fe + wid,
        # A zero between-country SD is a boundary fit, not evidence that
        # countries are identical.
        boundary_fit = isTRUE(tau < 1e-6), stringsAsFactors = FALSE)
    }
  }
  out <- dplyr::bind_rows(rows)
  out$response <- response
  # exp() because the response is a log biomarker: the interval is then a
  # multiplicative range for a new country's level.
  out$centre_ratio_scale <- exp(out$centre)
  out$lo_ratio_scale <- exp(out$lo)
  out$hi_ratio_scale <- exp(out$hi)
  out
}

# ── Exceedance probabilities ─────────────────────────────────────────────────

#' P(district prevalence exceeds a threshold), on the logit scale.
#'
#' @param pred predicted prevalence
#' @param se standard error of the prediction on the PREVALENCE scale
#' @param threshold prevalence threshold
#' @return numeric probability, NA where inputs are unusable
exceedance_probability <- function(pred, se, threshold) {
  ok <- is.finite(pred) & is.finite(se) & se > 0 & is.finite(threshold) &
        pred > 0 & pred < 1
  out <- rep(NA_real_, length(pred))
  # Delta-method transfer of the SE to the logit scale, so the normal
  # approximation is made where the quantity is unbounded rather than on a
  # proportion that would put mass outside [0, 1].
  z <- (stats::qlogis(pmin(pmax(threshold, 1e-6), 1 - 1e-6)) -
        stats::qlogis(pmin(pmax(pred[ok], 1e-6), 1 - 1e-6))) /
       (se[ok] / (pred[ok] * (1 - pred[ok])))
  out[ok] <- 1 - stats::pnorm(z)
  out
}
