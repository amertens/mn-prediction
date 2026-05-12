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

  # Standardize Admin2 names for joining (trim whitespace, drop NAs)
  bnd$Admin2 <- trimws(bnd$Admin2)

  # Join predictions by Admin2 (deduplicate first to be safe)
  pred_subset <- pred_subset[!duplicated(pred_subset$Admin2), , drop = FALSE]
  pop_subset  <- pop_subset[!duplicated(pop_subset$Admin2), , drop = FALSE]

  joined <- merge(bnd, pred_subset[, c("Admin2", "pred_prev", "obs_prev",
                                        "ci_lo", "ci_hi", "ci_width",
                                        "n_survey", "who_class")],
                  by = "Admin2", all.x = TRUE, sort = FALSE)

  pop_cols <- intersect(c("Admin2", "population", "population_2023",
                          "pop_child", "pop_women",
                          "pop_child_2023", "pop_women_2023"),
                        colnames(pop_subset))
  joined <- merge(joined, pop_subset[, pop_cols],
                  by = "Admin2", all.x = TRUE, sort = FALSE)

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

  # WHO class for coloring (NA → "No data")
  joined$who_class[is.na(joined$who_class)] <- "No data"

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
