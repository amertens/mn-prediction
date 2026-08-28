# =============================================================================
# R/corrected/p11_level_decomposition.R  —  P11: what is the LOCO level bias
# made of?
#
# The dominant leave-one-country-out failure mode in this project is a national
# LEVEL offset: the transported map has roughly the right shape and the wrong
# height. WS3's first deliverable established that harmonizing the inflammation
# adjustment moves only a small part of it (uniform_brinda, see
# scripts/run_uniform_brinda.R). This file asks the complementary question on
# the continuous biomarker scale: of each fold's national offset, how much is a
# pure country intercept that no covariate could have supplied?
#
# THE DECOMPOSITION
# -----------------
# Work at Admin-2 level on the log biomarker. For country c and area a let
# y_ca be the area mean. Standardize WITHIN survey:
#
#   z_ca = (y_ca - mu_c) / sd_c
#
# mu_c and sd_c are that survey's own location and scale. Standardizing within
# survey is what removes the level and spread differences from the modelling
# problem, so the fitted model is about spatial PATTERN only.
#
# Under LOCO, a model fitted on the training countries predicts z for the
# held-out country. At deployment the held-out country's own mu and sd are
# unknown, so the transported prediction has to be reconstituted with the
# training-pooled values:
#
#   yhat_ha = mu_train + sd_train * zhat_ha
#
# The national bias on the log scale is then mean(yhat) - mean(y_true), and it
# splits exactly into two terms:
#
#   location_offset = mu_train - mu_holdout
#   pattern_term    = sd_train * mean(zhat) - sd_holdout * mean(z_true)
#
# mean(z_true) is 0 by construction, so pattern_term reduces to
# sd_train * mean(zhat): whatever net level the covariates themselves assert.
# The two terms sum to the total by algebra, not by fitting, so the shares are
# not an approximation.
#
# `location_share` is the fraction of the absolute bias attributable to the
# pure intercept. A share near 1 means the covariates are not the problem and
# no amount of covariate work will fix the height; it is an argument for
# anchoring (WS7) rather than for better predictors.
#
# All of this is additive: nothing here is called from the production path.
# =============================================================================

#' Area-level table of the log biomarker plus covariates, for one outcome.
#'
#' @param svy_list named list of compute_svy_admin2()-shaped tables, one per
#'   country, carrying Admin2 and a continuous area mean in `area_mean`
#' @param gee_list named list of gee_admin2_* tables
#' @return list(data, predictors, countries) or NULL
build_level_decomp_dataset <- function(svy_list, gee_list) {
  p <- build_area_loco_dataset(svy_list, gee_list)
  d <- p$pooled_data
  if (is.null(d) || !nrow(d)) return(NULL)
  if (!"n_svy" %in% names(d)) d$n_svy <- 1
  list(data = d, predictors = p$common_gee_vars, countries = p$country_names)
}

#' Within-survey standardization of an area-level response.
#'
#' @param y numeric area-level response
#' @param country character vector of country labels
#' @return list(z, mu, sd) where mu and sd are named by country
.ld_standardize <- function(y, country) {
  mu <- tapply(y, country, mean, na.rm = TRUE)
  sd <- tapply(y, country, stats::sd,  na.rm = TRUE)
  sd[!is.finite(sd) | sd <= 0] <- 1
  # as.numeric() strips the names tapply() attaches through mu[country].
  # A named response reaches glmnet::cv.glmnet() as a named vector and makes it
  # fail with "non-conformable arrays"; .tr_fit_predict() catches that error and
  # silently returns the training mean, which here would have looked exactly
  # like "the covariates contribute nothing" (n_selected 0, zhat_sd 0) rather
  # than like a failed fit.
  list(z = as.numeric((y - mu[country]) / sd[country]), mu = mu, sd = sd)
}

#' Decompose each LOCO national level bias into intercept and pattern terms.
#'
#' The model is deliberately the same elastic net the area transport recipe
#' uses (.tr_prep_X / .tr_screen / .tr_fit_predict in R/transportability_area.R),
#' so the decomposition describes the estimator the project actually reports
#' rather than a separate model invented for the purpose. The only change is the
#' response: the within-survey standardized log biomarker instead of the
#' logit prevalence.
#'
#' @param pool build_level_decomp_dataset() output
#' @param outcome outcome tag, for labelling
#' @param recipe list like AREA_TRANSPORT_RECIPE
#' @param response_col column holding the area-level continuous response
#' @return data.frame, one row per held-out country
decompose_level_bias <- function(pool, outcome, recipe = AREA_TRANSPORT_RECIPE,
                                 response_col = "area_log_biomarker") {
  d <- pool$data
  if (!response_col %in% names(d))
    stop("decompose_level_bias(): missing response column ", response_col, call. = FALSE)
  d <- d[is.finite(d[[response_col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  st <- .ld_standardize(d[[response_col]], d$country)
  d$.z <- st$z

  rows <- list()
  for (ho in unique(d$country)) {
    tr <- d[d$country != ho, , drop = FALSE]
    te <- d[d$country == ho, , drop = FALSE]
    if (nrow(tr) < 10 || nrow(te) < 4) next

    pp  <- .tr_prep_X(tr, te, pool$predictors)
    Xtr <- pp$Xtr; Xte <- pp$Xte
    w   <- if (isTRUE(recipe$weight)) pmax(tr$n_svy, 1, na.rm = TRUE) else NULL
    if (!is.null(recipe$screen_K) && is.finite(recipe$screen_K)) {
      sel <- .tr_screen(Xtr, tr$.z, recipe$screen_K)
      Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
    }
    fp   <- .tr_fit_predict(Xtr, tr$.z, Xte, recipe, w)
    zhat <- fp$pred

    # Training-pooled location and scale: what a deployment in an unsurveyed
    # country would actually have to use.
    mu_train <- mean(tr[[response_col]], na.rm = TRUE)
    sd_train <- stats::sd(tr[[response_col]], na.rm = TRUE)
    if (!is.finite(sd_train) || sd_train <= 0) sd_train <- 1
    mu_ho <- unname(st$mu[ho]); sd_ho <- unname(st$sd[ho])

    yhat  <- mu_train + sd_train * zhat
    ytrue <- te[[response_col]]

    total    <- mean(yhat, na.rm = TRUE) - mean(ytrue, na.rm = TRUE)
    location <- mu_train - mu_ho
    pattern  <- sd_train * mean(zhat, na.rm = TRUE) -
                sd_ho    * mean(te$.z, na.rm = TRUE)

    rows[[ho]] <- data.frame(
      outcome        = outcome,
      held_out       = ho,
      n_test         = nrow(te),
      n_train        = nrow(tr),
      n_predictors   = ncol(Xtr),
      # Diagnostics that distinguish "the covariates contribute nothing" from
      # "the fit failed". When n_selected is 0 the elastic net has shrunk to the
      # intercept, and because every country's z is centred within survey the
      # pooled training mean is exactly 0, so zhat is exactly 0 and the pattern
      # term vanishes by construction rather than by accident.
      n_selected     = length(fp$vars),
      zhat_sd        = stats::sd(zhat),
      response       = response_col,
      mu_train       = mu_train,
      mu_holdout     = mu_ho,
      sd_train       = sd_train,
      sd_holdout     = sd_ho,
      # exp() of a log-scale mean difference is the multiplicative level ratio
      # between the training pool and the held-out survey, which is the units a
      # reader can check against the raw biomarker medians.
      level_ratio_train_over_holdout = exp(mu_train - mu_ho),
      total_bias_log     = total,
      location_offset_log = location,
      pattern_term_log    = pattern,
      abs_location_log = abs(location),
      abs_pattern_log  = abs(pattern),
      # Share of the realised bias. UNBOUNDED on purpose: when the covariate
      # term opposes the intercept the two partly cancel, |total| < |location|,
      # and the share exceeds 1. A fold where they nearly cancel produces a very
      # large value, which is informative about that fold but useless as a
      # summary statistic, so it is reported per fold and never averaged.
      location_share = if (abs(total) < 1e-8) NA_real_ else abs(location) / abs(total),
      pattern_share  = if (abs(total) < 1e-8) NA_real_ else abs(pattern)  / abs(total),
      # Bounded companion in [0, 1]: of the total absolute movement the two
      # terms generate, how much of it is the intercept. This is the statistic
      # to summarise across folds.
      location_share_bounded = if ((abs(location) + abs(pattern)) < 1e-12) NA_real_
                               else abs(location) / (abs(location) + abs(pattern)),
      # Does the covariate term push in the same direction as the offset?
      pattern_opposes_location = sign(pattern) != sign(location),
      pearson_z = suppressWarnings(stats::cor(te$.z, zhat, use = "complete.obs")),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
