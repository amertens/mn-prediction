# =============================================================================
# Helper functions used across modules
# =============================================================================

#' Format a proportion as a percentage string
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "—", sprintf(paste0("%.", digits, "f%%"), x * 100))
}

#' Format an integer count with comma separators
fmt_count <- function(x) {
  ifelse(is.na(x), "—", formatC(round(x), format = "d", big.mark = ","))
}

#' Methods note styling — used at the bottom of each table
methods_note <- function(...) {
  htmltools::div(
    class = "methods-note",
    style = paste(
      "font-size: 0.85em;",
      "color: #555;",
      "background: #f8f9fa;",
      "border-left: 3px solid #2c7bb6;",
      "padding: 0.75em 1em;",
      "margin-top: 0.5em;",
      "margin-bottom: 1em;",
      "line-height: 1.5;"
    ),
    htmltools::tags$strong("Methods: "),
    ...
  )
}

#' Get the joined Admin-2 sf for a country with predictions merged
#'
#' Returns sf object with one row per Admin-2 polygon, with prediction columns
#' merged in (NA for districts not covered by the model).
#'
#' @param pop_year Either "survey" (default — uses survey-year population) or
#'   "2023" (uses projected 2023 population for current burden estimates).
get_country_admin2 <- function(ctry_key, oc_key,
                                admin2_bnds, admin2_pred, admin2_pop,
                                pop_year = "survey") {
  bnd <- admin2_bnds[[ctry_key]]
  if (is.null(bnd)) return(NULL)

  ctry_label <- meta$countries[ctry_key]

  pred_subset <- admin2_pred[admin2_pred$country == ctry_label &
                              admin2_pred$outcome == oc_key, ,
                              drop = FALSE]

  pop_subset <- admin2_pop[admin2_pop$country == ctry_label, , drop = FALSE]

  # Standardize keys for joining (trim whitespace, drop NAs)
  bnd$Admin2 <- trimws(bnd$Admin2)
  if ("Admin1" %in% names(bnd)) bnd$Admin1 <- trimws(as.character(bnd$Admin1))

  # JOIN ON THE (Admin1, Admin2) PAIR WHERE BOTH SIDES CARRY IT.
  #
  # GADM admin-2 names are not unique within a country: Malawi has 243
  # polygons under 239 names -- TA Lundu, TA Malemia, TA Ngabu and TA Pemba
  # each name two genuinely different districts in different regions. The old
  # code joined on the name alone and, worse, called
  #   pred_subset[!duplicated(pred_subset$Admin2), ]
  # "to be safe", which DISCARDED the second district of each pair and painted
  # both polygons with the first one's prevalence and population. That is a
  # silent wrong answer on the map, not a safety measure.
  #
  # Deduplication is still needed, but on the KEY that actually identifies a
  # district. Where a side lacks Admin1 we fall back to the name and warn, so a
  # regenerated data file that drops the column is visible rather than silent.
  key_of <- function(d) {
    if ("Admin1" %in% names(d) &&
        any(!is.na(d$Admin1) & nzchar(as.character(d$Admin1))))
      paste(trimws(as.character(d$Admin1)), trimws(as.character(d$Admin2)),
            sep = "")
    else trimws(as.character(d$Admin2))
  }
  pair_ok <- function(d) "Admin1" %in% names(d) &&
    any(!is.na(d$Admin1) & nzchar(as.character(d$Admin1)))

  by_pred <- if (pair_ok(bnd) && pair_ok(pred_subset)) c("Admin1", "Admin2") else "Admin2"
  by_pop  <- if (pair_ok(bnd) && pair_ok(pop_subset))  c("Admin1", "Admin2") else "Admin2"
  if (!identical(by_pred, c("Admin1", "Admin2")) && anyDuplicated(bnd$Admin2))
    warning(sprintf(paste0("[dashboard] %s: predictions lack Admin1 and this ",
                           "country has duplicate Admin-2 names; the map will ",
                           "paint same-named districts identically."), ctry_label))

  pred_subset <- pred_subset[!duplicated(key_of(pred_subset)), , drop = FALSE]
  pop_subset  <- pop_subset[!duplicated(key_of(pop_subset)), , drop = FALSE]

  joined <- merge(bnd, pred_subset[, c(by_pred, "pred_prev", "obs_prev",
                                        "ci_lo", "ci_hi", "ci_width",
                                        "n_survey", "who_class")],
                  by = by_pred, all.x = TRUE, sort = FALSE)

  pop_cols <- intersect(c(by_pop, "population", "population_2023",
                          "pop_child", "pop_women",
                          "pop_child_2023", "pop_women_2023"),
                        colnames(pop_subset))
  joined <- merge(joined, pop_subset[, pop_cols],
                  by = by_pop, all.x = TRUE, sort = FALSE)

  # Pick the subgroup-specific population matching the outcome.
  # Outcomes starting with "child_" use child population; "women_" use women.
  is_child_outcome <- startsWith(oc_key, "child_")
  if (pop_year == "2023") {
    joined$pop_total <- joined$population_2023
    joined$pop_subgroup <- if (is_child_outcome) joined$pop_child_2023
                            else joined$pop_women_2023
    joined$pop_year_label <- "2023 projection"
  } else {
    joined$pop_total <- joined$population
    joined$pop_subgroup <- if (is_child_outcome) joined$pop_child
                            else joined$pop_women
    joined$pop_year_label <- "Survey year"
  }
  joined$population <- joined$pop_subgroup  # the relevant denominator
  joined$subgroup_label <- if (is_child_outcome) "children 6–59 months"
                            else "women 15–49 years"

  # Compute population at risk (subgroup × prevalence)
  joined$pop_at_risk <- joined$pred_prev * joined$population

  # Difference between the modeled (within-country SuperLearner) prevalence and
  # the survey-observed prevalence (NA where there is no survey value).
  joined$diff_prev <- joined$pred_prev - joined$obs_prev

  # Transportability (leave-one-country-out) modeled prevalence and its
  # difference from the local survey. The transportability model uses a
  # HARMONIZED, cross-country comparable outcome (uniform WHO cutoffs:
  # ferritin<12/15 = ID, RBP<0.70 = VAD), so the difference is taken against
  # that file's OWN survey prevalence (loco_survey_prev) rather than the
  # dashboard's survey-native obs_prev. Available only for the country ×
  # outcome combinations covered by the transportability model; NA otherwise.
  joined$loco_pred_prev <- NA_real_
  joined$loco_survey_prev <- NA_real_
  joined$loco_diff <- NA_real_
  if (exists("loco_pred") && !is.null(loco_pred) && nrow(loco_pred) > 0) {
    lp <- loco_pred[loco_pred$country == ctry_key &
                      loco_pred$outcome == oc_key, , drop = FALSE]
    if (nrow(lp) > 0) {
      # Same pair-key reasoning as the prediction join above: matching on the
      # name alone assigns one district's transported prevalence to both
      # members of a duplicate-named pair.
      lp <- lp[!duplicated(key_of(lp)), , drop = FALSE]
      idx <- match(key_of(joined), key_of(lp))
      joined$loco_pred_prev <- lp$loco_modeled_prev[idx]
      joined$loco_survey_prev <- lp$loco_survey_prev[idx]
      joined$loco_diff <- joined$loco_pred_prev - joined$loco_survey_prev
    }
  }

  # WHO class for coloring (NA → "No data")
  joined$who_class[is.na(joined$who_class)] <- "No data"

  joined
}

#' Get the joined Admin-3 (chiefdom) sf for a country with predictions merged.
#'
#' Currently available for Sierra Leone only. Mirrors get_country_admin2()'s
#' output schema (one row per Admin-3 polygon, prediction columns merged, NA
#' where not covered) so the rest of the map logic can stay level-agnostic.
#' Population denominators are not available at Admin-3, so pop_* and
#' pop_at_risk are NA; the headline then falls back to a no-population display.
#' Transportability (LOCO) is not computed at Admin-3 either.
get_country_admin3 <- function(ctry_key, oc_key, admin3_bnds, admin3_pred) {
  if (is.null(admin3_bnds) || is.null(admin3_pred)) return(NULL)
  bnd <- admin3_bnds[[ctry_key]]
  if (is.null(bnd)) return(NULL)

  ctry_label <- meta$countries[ctry_key]
  pred_subset <- admin3_pred[admin3_pred$country == ctry_label &
                               admin3_pred$outcome == oc_key, , drop = FALSE]

  bnd$Admin3 <- trimws(bnd$Admin3)
  pred_subset <- pred_subset[!duplicated(pred_subset$Admin3), , drop = FALSE]

  joined <- merge(bnd, pred_subset[, c("Admin3", "pred_prev", "obs_prev",
                                       "ci_lo", "ci_hi", "ci_width",
                                       "n_survey", "who_class")],
                  by = "Admin3", all.x = TRUE, sort = FALSE)

  # No Admin-3 population denominator — keep schema parity but leave NA.
  is_child_outcome <- startsWith(oc_key, "child_")
  joined$pop_total      <- NA_real_
  joined$pop_subgroup   <- NA_real_
  joined$population     <- NA_real_
  joined$pop_at_risk    <- NA_real_
  joined$pop_year_label <- "Survey year"
  joined$subgroup_label <- if (is_child_outcome) "children 6–59 months"
                            else "women 15–49 years"

  joined$diff_prev        <- joined$pred_prev - joined$obs_prev
  joined$loco_pred_prev   <- NA_real_
  joined$loco_survey_prev <- NA_real_
  joined$loco_diff        <- NA_real_

  joined$who_class[is.na(joined$who_class)] <- "No data"
  joined
}

#' Aggregate Admin-2 predictions to Admin-1 polygons.
#'
#' Population-weighted aggregation of prevalence, CIs, and survey n. Conformal
#' CIs were already computed at Admin-1 in the pipeline and broadcast to
#' Admin-2, so all Admin-2 rows within an Admin-1 region share the same
#' ci_lo/ci_hi; we recover that value as the population-weighted mean (= the
#' shared value when populations are non-negative). WHO classification at
#' Admin-1 is the population-weighted modal class within the region — ties
#' broken by max severity.
#'
#' Returns an sf object keyed on Admin1, with the same column schema as
#' get_country_admin2() so the rest of the map logic can be level-agnostic.
get_country_admin1 <- function(ctry_key, oc_key,
                                admin1_bnds, admin2_bnds,
                                admin2_pred, admin2_pop,
                                pop_year = "survey") {
  a2 <- get_country_admin2(ctry_key, oc_key,
                            admin2_bnds, admin2_pred, admin2_pop,
                            pop_year = pop_year)
  if (is.null(a2) || nrow(a2) == 0) return(NULL)
  bnd1 <- admin1_bnds[[ctry_key]]
  if (is.null(bnd1)) return(NULL)

  a2_df <- as.data.frame(a2)
  a2_df$geometry <- NULL

  # Drop rows with no Admin1 assignment or no population — can't aggregate
  a2_df <- a2_df[!is.na(a2_df$Admin1) & !is.na(a2_df$population) &
                  a2_df$population > 0, , drop = FALSE]
  if (nrow(a2_df) == 0) return(NULL)

  pop_weighted <- function(x, w) {
    keep <- !is.na(x) & !is.na(w) & w > 0
    if (!any(keep)) return(NA_real_)
    sum(x[keep] * w[keep]) / sum(w[keep])
  }

  agg <- do.call(rbind, lapply(split(a2_df, a2_df$Admin1), function(g) {
    w <- g$population

    # Modal WHO class, ties broken by severity rank
    sev_rank <- c("Low" = 1, "Mild" = 2, "Moderate" = 3,
                  "Severe" = 4, "No data" = 0)
    who_tab <- tapply(w, g$who_class, sum, na.rm = TRUE)
    if (length(who_tab) == 0 || all(is.na(who_tab))) {
      who_a1 <- "No data"
    } else {
      max_w <- max(who_tab, na.rm = TRUE)
      cands <- names(who_tab)[!is.na(who_tab) & who_tab == max_w]
      who_a1 <- cands[which.max(sev_rank[cands])]
    }

    data.frame(
      Admin1        = g$Admin1[1],
      pred_prev     = pop_weighted(g$pred_prev,  w),
      obs_prev      = pop_weighted(g$obs_prev,   w),
      loco_pred_prev = pop_weighted(g$loco_pred_prev, w),
      loco_survey_prev = pop_weighted(g$loco_survey_prev, w),
      ci_lo         = pop_weighted(g$ci_lo,      w),
      ci_hi         = pop_weighted(g$ci_hi,      w),
      ci_width      = pop_weighted(g$ci_width,   w),
      n_survey      = sum(g$n_survey, na.rm = TRUE),
      population    = sum(g$population, na.rm = TRUE),
      who_class     = who_a1,
      stringsAsFactors = FALSE
    )
  }))
  agg$pop_at_risk <- agg$pred_prev * agg$population
  agg$diff_prev   <- agg$pred_prev - agg$obs_prev
  agg$loco_diff   <- agg$loco_pred_prev - agg$loco_survey_prev

  bnd1$Admin1 <- trimws(bnd1$Admin1)
  joined <- merge(bnd1, agg, by = "Admin1", all.x = TRUE, sort = FALSE)
  joined$who_class[is.na(joined$who_class)] <- "No data"

  # Add Admin2 column = NA so downstream code can reuse the same field names
  # safely (we expose "area" for label/click logic instead).
  joined$Admin2 <- NA_character_
  joined$pop_year_label <- attr(a2, "pop_year_label") %||%
                            (if (pop_year == "2023") "2023 projection"
                             else "Survey year")
  joined
}

#' Confidence badge — categorize CI width into low/moderate/high confidence
#'
#' Used as a visual companion to the numeric interval. CI widths are in
#' proportion units (0–1).
confidence_badge <- function(ci_width) {
  if (is.na(ci_width)) return("Unknown")
  if (ci_width < 0.10) return("High confidence")
  if (ci_width < 0.20) return("Moderate confidence")
  return("Low confidence")
}

#' Aggregate Admin-2 predictions to national level with population weighting
#'
#' Point estimate is a population-weighted average of district predictions.
#' Confidence intervals are propagated using a partial-correlation assumption:
#' districts within the same Admin-1 region share their CI (because the
#' pipeline computes conformal CIs at Admin-1 level), so for variance
#' propagation we treat each Admin-1 as one independent observation.
#'
#' Returns a one-row data.frame with point estimate and approximate 95% CI,
#' weighted by population.
national_aggregate <- function(df) {
  df <- as.data.frame(df)

  # No-population fallback (e.g. the Admin-3 layer, which has no chiefdom
  # denominators): report an UNWEIGHTED mean prevalence + interval across areas
  # and leave population counts NA. The headline then hides the population line.
  if (!("population" %in% colnames(df)) || all(is.na(df$population))) {
    dd <- df[!is.na(df$pred_prev), , drop = FALSE]
    if (nrow(dd) == 0) {
      return(data.frame(
        pop_total = NA_real_, pred_prev_natl = NA_real_,
        pop_at_risk_natl = NA_real_, ci_lo_natl = NA_real_,
        ci_hi_natl = NA_real_, ci_method = NA_character_,
        stringsAsFactors = FALSE))
    }
    return(data.frame(
      pop_total = NA_real_,
      pred_prev_natl = mean(dd$pred_prev, na.rm = TRUE),
      pop_at_risk_natl = NA_real_,
      ci_lo_natl = if (!all(is.na(dd$ci_lo))) mean(dd$ci_lo, na.rm = TRUE) else NA_real_,
      ci_hi_natl = if (!all(is.na(dd$ci_hi))) mean(dd$ci_hi, na.rm = TRUE) else NA_real_,
      ci_method = "unweighted mean (no population denominator)",
      stringsAsFactors = FALSE))
  }

  df <- df[!is.na(df$pred_prev) & !is.na(df$population), , drop = FALSE]
  if (nrow(df) == 0) {
    return(data.frame(
      pop_total = NA_real_,
      pred_prev_natl = NA_real_,
      pop_at_risk_natl = NA_real_,
      ci_lo_natl = NA_real_,
      ci_hi_natl = NA_real_,
      ci_method = NA_character_
    ))
  }

  w <- df$population / sum(df$population, na.rm = TRUE)
  pred_natl <- sum(w * df$pred_prev, na.rm = TRUE)

  # CI propagation: aggregate CIs at the Admin-1 level (where they are
  # actually computed by the pipeline), then population-weight across
  # Admin-1 regions. This avoids overstating uncertainty by treating
  # each district as independent.
  ci_method <- "unavailable"
  ci_lo_natl <- NA_real_
  ci_hi_natl <- NA_real_

  if ("Admin1" %in% colnames(df) &&
      !all(is.na(df$ci_lo)) && !all(is.na(df$ci_hi))) {
    # Group by Admin1, pop-weight predictions and CIs within each region
    a1 <- aggregate(
      list(
        pop = df$population,
        pop_x_pred = df$population * df$pred_prev,
        pop_x_lo = df$population * df$ci_lo,
        pop_x_hi = df$population * df$ci_hi
      ),
      by = list(Admin1 = df$Admin1),
      FUN = sum, na.rm = TRUE
    )
    a1$pred_a1 <- a1$pop_x_pred / a1$pop
    a1$lo_a1 <- a1$pop_x_lo / a1$pop
    a1$hi_a1 <- a1$pop_x_hi / a1$pop

    # National point estimate (matches earlier — sanity check)
    w_a1 <- a1$pop / sum(a1$pop)

    # Variance propagation (conservative): each Admin-1 region's
    # conformal half-width treated as 1 SD; combine via population
    # weights with assumed independence ACROSS regions.
    half_width_a1 <- (a1$hi_a1 - a1$lo_a1) / 2
    se_a1 <- half_width_a1 / qnorm(0.975)  # ≈ /1.96
    se_natl <- sqrt(sum(w_a1^2 * se_a1^2))

    ci_lo_natl <- max(0, pred_natl - qnorm(0.975) * se_natl)
    ci_hi_natl <- min(1, pred_natl + qnorm(0.975) * se_natl)
    ci_method <- "Admin-1 weighted variance propagation"
  }

  pop_total <- sum(df$population, na.rm = TRUE)

  data.frame(
    pop_total        = pop_total,
    pred_prev_natl   = pred_natl,
    pop_at_risk_natl = pred_natl * pop_total,
    ci_lo_natl       = ci_lo_natl,
    ci_hi_natl       = ci_hi_natl,
    ci_method        = ci_method,
    stringsAsFactors = FALSE
  )
}
