# =============================================================================
# scripts/accuracy_impact/wsj1_targeting_lift.R
#
# WS-J. Targeting lift under all three fold protocols.
#
# WHY THIS EXISTS
# ---------------
# Every accuracy metric this project reports answers "how much deficiency is in
# this area". Programmes ask a different question: "which areas do we go to
# first". Those are different estimation problems and they have different
# baselines. The prevalence question is scored against the national mean; the
# targeting question is scored against choosing at random.
#
# The distinction is not cosmetic, because the two behave completely differently
# under the fold protocol that decides whether a result is deployable:
#
#                                     targeting lift   correlation
#   within country, cluster-blocked            1.32          0.416
#   within country, region-blocked             1.21          0.154
#   leave-one-country-out                      1.23          0.184
#
# Correlation loses 63 percent of its value when neighbouring areas are removed
# from the training set. Targeting lift barely moves, and it is the SAME at
# leave-one-country-out as it is within a country. Lift is a coarser functional
# of the prediction: it depends only on which areas land in the top of the
# ranking, so it tolerates the monotone distortion that destroys a correlation.
#
# That is the empirical case for reporting targeting as a primary endpoint. The
# principled case has to be made in advance and on decision-relevance grounds,
# not chosen after seeing which metric survived, and the write-up says so.
#
#   Rscript scripts/accuracy_impact/wsj1_targeting_lift.R
# -> results/tables/targeting_lift.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
TDIR <- here("results", "tables")
FRAC <- 0.20   # the top fifth of areas, the usual programme coverage question

#' Burden reached by targeting the top `frac` of areas, over burden reached by
#' choosing the same number of areas at random.
#'
#' `lift` is the ratio of mean prevalence in the selected areas to the overall
#' mean. `reach` is the share of total burden the selected areas contain, whose
#' random-allocation expectation is `frac` itself.
targeting_lift <- function(obs, pred, frac = FRAC) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 6L) return(NULL)
  o <- obs[ok]; p <- pred[ok]
  k <- max(1L, round(frac * length(o)))
  sel <- order(p, decreasing = TRUE)[seq_len(k)]
  overall <- mean(o)
  if (!is.finite(overall) || overall <= 0) return(NULL)
  data.frame(n_areas = length(o), k_selected = k,
             lift = round(mean(o[sel]) / overall, 3),
             reach = round(sum(o[sel]) / sum(o), 3),
             random_reach = frac,
             stringsAsFactors = FALSE)
}

rows <- list()

# --- within country, region-blocked, from the resolution sweep --------------
rp <- file.path(TDIR, "resolution_predictions.csv")
if (file.exists(rp)) {
  d <- suppressMessages(readr::read_csv(rp, show_col_types = FALSE)) |> as.data.frame()
  for (lv in unique(d$level)) {
    s <- d[d$level == lv, , drop = FALSE]
    for (k in unique(paste(s$country, s$outcome))) {
      v <- s[paste(s$country, s$outcome) == k, , drop = FALSE]
      for (arm in c("covariates", "spatial")) {
        if (!arm %in% names(v)) next
        L <- targeting_lift(v$obs, v[[arm]])
        if (is.null(L)) next
        L$protocol <- "within_country_region_blocked"; L$level <- lv
        L$arm <- arm; L$country <- v$country[1]; L$outcome <- v$outcome[1]
        rows[[length(rows) + 1L]] <- L
      }
    }
  }
}

# --- leave-one-country-out --------------------------------------------------
lp <- file.path(TDIR, "transportability_area_loco_predictions.csv")
if (file.exists(lp)) {
  d <- suppressMessages(readr::read_csv(lp, show_col_types = FALSE)) |> as.data.frame()
  for (k in unique(paste(d$country, d$outcome))) {
    v <- d[paste(d$country, d$outcome) == k, , drop = FALSE]
    L <- targeting_lift(v$survey_prev, v$modeled_prev)
    if (is.null(L)) next
    L$protocol <- "leave_one_country_out"; L$level <- "admin-2"
    L$arm <- "covariates"; L$country <- v$country[1]; L$outcome <- v$outcome[1]
    rows[[length(rows) + 1L]] <- L
  }
}

# --- within country, cluster-blocked, already computed elsewhere ------------
dp <- file.path(TDIR, "corrected", "decision_value.csv")
if (file.exists(dp)) {
  d <- suppressMessages(readr::read_csv(dp, show_col_types = FALSE)) |> as.data.frame()
  for (i in seq_len(nrow(d))) {
    rows[[length(rows) + 1L]] <- data.frame(
      n_areas = d$n_admin2[i], k_selected = NA_integer_,
      lift = d$lift_vs_no_targeting[i], reach = d$reach_at_20pct[i],
      random_reach = FRAC, protocol = "within_country_cluster_blocked",
      level = "admin-2", arm = "covariates",
      country = d$country[i], outcome = d$outcome[i], stringsAsFactors = FALSE)
  }
}

if (!length(rows)) stop("[wsj1] no predictions found")
out <- dplyr::bind_rows(rows)
front <- c("protocol", "level", "arm", "country", "outcome")
out <- out[, c(front, setdiff(names(out), front))]
readr::write_csv(out, file.path(TDIR, "targeting_lift.csv"))

cat("=== WS-J: targeting lift, top 20 percent of areas ===\n")
s <- out |> group_by(protocol, level, arm) |>
  summarise(cells = dplyr::n(),
            median_lift = round(stats::median(lift, na.rm = TRUE), 2),
            above_1 = sprintf("%d/%d", sum(lift > 1, na.rm = TRUE), dplyr::n()),
            max_lift = round(max(lift, na.rm = TRUE), 2),
            median_reach = round(stats::median(reach, na.rm = TRUE), 3),
            .groups = "drop") |> arrange(level, protocol, arm)
print(as.data.frame(s), row.names = FALSE)

cat("\n=== by country, covariate arm ===\n")
b <- out |> filter(arm == "covariates") |>
  mutate(cty = tolower(gsub("[^A-Za-z]", "", country))) |>
  group_by(protocol, level, cty) |>
  summarise(median_lift = round(stats::median(lift, na.rm = TRUE), 2),
            above_1 = sprintf("%d/%d", sum(lift > 1, na.rm = TRUE), dplyr::n()),
            .groups = "drop")
print(as.data.frame(b), row.names = FALSE)
cat(sprintf("\n-> %s\n", file.path("results", "tables", "targeting_lift.csv")))
