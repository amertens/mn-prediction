# =============================================================================
# R/aggregation_uncertainty.R
#
# Bootstrap-style uncertainty for population-weighted aggregated quantities
# such as national prevalence, "hidden burden" populations, and "cases averted"
# under coverage scenarios.
#
# We sample district-level prevalences from their conformal prediction
# intervals (treating the interval bounds as approximating a normal SD with
# scale half_width / qnorm(0.975)), recompute the aggregate quantity, and
# return its empirical distribution.
#
# This is a parametric bootstrap that propagates *prediction* uncertainty
# (not population uncertainty — population is treated as fixed).
# =============================================================================


#' Bootstrap a population-weighted aggregate from district-level predictions
#'
#' @param admin2_df data.frame with columns: pred_prev, ci_lo, ci_hi, population.
#' @param fun A function taking a numeric vector of bootstrap-sampled
#'   prevalences and a numeric vector of populations, returning a scalar.
#'   Examples below.
#' @param B Number of bootstrap replicates (default 1000)
#' @param level Confidence level (default 0.95)
#' @return list with point, ci_lo, ci_hi, samples (length B)
bootstrap_aggregate <- function(admin2_df, fun, B = 1000L, level = 0.95) {

  df <- admin2_df[!is.na(admin2_df$pred_prev) &
                    !is.na(admin2_df$population), , drop = FALSE]
  if (nrow(df) == 0) {
    return(list(point = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                samples = numeric(0)))
  }

  # Convert each district's CI to a normal SD (assumes interval is symmetric
  # around the point estimate with normal tails — an approximation).
  half_w <- (df$ci_hi - df$ci_lo) / 2
  half_w[is.na(half_w)] <- 0
  sigma <- half_w / stats::qnorm(0.975)

  # Bootstrap
  samples <- numeric(B)
  for (b in seq_len(B)) {
    perturbed <- df$pred_prev + stats::rnorm(nrow(df), 0, sigma)
    perturbed <- pmin(pmax(perturbed, 0), 1)  # clip to [0,1]
    samples[b] <- fun(perturbed, df$population)
  }

  alpha <- 1 - level
  list(
    point   = fun(df$pred_prev, df$population),
    ci_lo   = quantile(samples, alpha / 2, names = FALSE, na.rm = TRUE),
    ci_hi   = quantile(samples, 1 - alpha / 2, names = FALSE, na.rm = TRUE),
    samples = samples
  )
}


# ── Standard aggregate functions ───────────────────────────────────────────

#' Population-weighted national prevalence
agg_natl_prev <- function(prev, pop) {
  stats::weighted.mean(prev, pop, na.rm = TRUE)
}

#' Total "people affected" count (prevalence × population, summed)
agg_pop_at_risk <- function(prev, pop) {
  sum(prev * pop, na.rm = TRUE)
}

#' "Hidden burden": sum of populations in districts whose prevalence
#' deviates from the population-weighted mean by more than a threshold.
agg_hidden_burden <- function(prev, pop, threshold = 0.05) {
  natl <- stats::weighted.mean(prev, pop, na.rm = TRUE)
  diverging <- abs(prev - natl) > threshold
  sum(pop[diverging], na.rm = TRUE)
}

#' Cases averted under a coverage × effect-size intervention applied to a
#' targeted subset of districts.
#'
#' @param prev Vector of district baseline prevalences.
#' @param pop Vector of district populations.
#' @param target Logical vector indicating which districts are targeted.
#' @param coverage Proportion of population in targeted districts reached.
#' @param effect Proportional reduction in prevalence among reached individuals.
agg_cases_averted <- function(prev, pop, target, coverage, effect) {
  prev_after <- prev
  prev_after[target] <- prev[target] * (1 - coverage * effect)
  sum((prev - prev_after) * pop, na.rm = TRUE)
}
