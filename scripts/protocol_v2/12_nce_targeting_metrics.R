# =============================================================================
# scripts/protocol_v2/12_nce_targeting_metrics.R
#
# The numbers the NCE actually needs, under the corrected protocol.
#
# WHY THIS IS A SEPARATE, FOCUSED RUN
# -----------------------------------
# The main leaderboard reports rank correlation and top-quartile OVERLAP. The
# NCE makes a different, more programmatic claim: "visiting the top fifth of
# districts reaches X percent of the country's deficiency burden". That is
# BURDEN capture, weighted by population, and it is not the same quantity as
# overlap. This script computes it directly, for only the arms and estimands
# the NCE cites, so it is minutes rather than the ~40 of a full re-run.
#
# BURDEN, DEFINED PROPERLY
#   burden_d = prevalence_d x population_d      (cases, not rates)
# using the child or women population as the outcome requires, from
# dashboard/data/admin2_population.rds. A district that is 40 percent deficient
# with 5,000 children carries less burden than one 20 percent deficient with
# 80,000, and a targeting claim that ignores that is not a programme claim.
#
# WHAT IS REPORTED PER CELL
#   capture_top20   share of national burden inside the model's worst-ranked 20%
#   lift            capture_top20 / 0.20  (1.0 = no better than picking at random)
#   prev_top20      population-weighted prevalence of the selected districts
#   prev_national   population-weighted national prevalence
#   ratio_prev      prev_top20 / prev_national - "how much worse are the
#                   districts we send you to than the country as a whole"
#
# ARMS: the model, and the two things a programme would otherwise do -
#   null_train_mean  apply one national number to every district (the status quo
#                    where no sub-national data exists)
#   region_mean_jk   use the survey's own regional averages, jackknifed so a
#                    district never sees its own respondents (in-fill only)
#   oracle           rank by the TRUE district prevalence - the ceiling, so a
#                    reader can see how much of the attainable gain is captured
#
#   Rscript scripts/protocol_v2/12_nce_targeting_metrics.R
# -> results/tables/protocol_v2/nce_targeting_metrics.csv
# -> results/tables/protocol_v2/nce_targeting_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
TOPFRAC <- 0.20
REPS <- as.integer(Sys.getenv("NCE_REPS", "10"))
set.seed(20260991L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
POP <- readRDS("dashboard/data/admin2_population.rds")
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

#' population relevant to the outcome's target group
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"

build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")],
               by = c("Admin1", "Admin2")) |>
    left_join(POP[POP$country == cn, c("Admin2", pop_for(on))],
              by = "Admin2")
  names(m)[names(m) == pop_for(on)] <- "pop"
  m <- m[is.finite(m$pop) & m$pop > 0, ]
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  list(country = cn, outcome = on, n = nrow(m), y = m$y_prev, pop = m$pop,
       X = Xr, D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = m$y_prev),
       Admin1 = m$Admin1, w = m$n_eff)
}

#' burden captured by the worst-ranked TOPFRAC of districts
capture <- function(y, pop, score) {
  ok <- is.finite(y) & is.finite(pop) & is.finite(score)
  if (sum(ok) < 5) return(c(NA, NA, NA, NA))
  y <- y[ok]; pop <- pop[ok]; score <- score[ok]
  burden <- y * pop
  k <- max(1, round(TOPFRAC * length(y)))
  sel <- order(score, decreasing = TRUE)[seq_len(k)]
  cap <- sum(burden[sel]) / sum(burden)
  prev_sel <- sum(burden[sel]) / sum(pop[sel])
  prev_nat <- sum(burden) / sum(pop)
  c(cap, cap / TOPFRAC, prev_sel, prev_nat)
}

rows <- list()
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  ymod <- .v2_logit(cl$y)

  # ---- A. in-fill, replicated ----
  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    for (a in c("null_train_mean", "region_mean_jk", "spatial",
                "domain_index", "spatial_plus_domain")) {
      fn <- ARMS_V2[[a]]
      pred <- rep(NA_real_, cl$n)
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12) next
        p <- tryCatch(fn(tr, te, ymod, cl$X, cl$D, cl$aux),
                      error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p
      }
      cp <- capture(cl$y, cl$pop, pred)
      rows[[paste(cn, on, "infill", a, r)]] <- data.frame(
        country = cn, outcome = on, estimand = "infill", arm = a, rep = r,
        n_areas = cl$n, capture_top20 = cp[1], lift = cp[2],
        prev_top20 = cp[3], prev_national = cp[4])
    }
  }
  # oracle ceiling and the random floor, once
  co <- capture(cl$y, cl$pop, cl$y)
  rows[[paste(cn, on, "oracle")]] <- data.frame(
    country = cn, outcome = on, estimand = "infill", arm = "oracle_ceiling",
    rep = 1L, n_areas = cl$n, capture_top20 = co[1], lift = co[2],
    prev_top20 = co[3], prev_national = co[4])
  cat("done", cn, on, "\n")
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "nce_targeting_metrics.csv"), row.names = FALSE)

CELL <- R |> group_by(country, outcome, estimand, arm) |>
  summarise(across(c(capture_top20, lift, prev_top20, prev_national),
                   ~ mean(.x, na.rm = TRUE)), n_areas = max(n_areas),
            .groups = "drop")
SUMM <- CELL |> group_by(estimand, arm) |>
  summarise(cells = n(),
            mean_capture = round(mean(capture_top20, na.rm = TRUE), 3),
            median_capture = round(median(capture_top20, na.rm = TRUE), 3),
            mean_lift = round(mean(lift, na.rm = TRUE), 3),
            cells_lift_gt1 = sum(lift > 1, na.rm = TRUE),
            mean_prev_top20 = round(100 * mean(prev_top20, na.rm = TRUE), 1),
            mean_prev_national = round(100 * mean(prev_national, na.rm = TRUE), 1),
            .groups = "drop") |> arrange(estimand, desc(mean_capture))
write.csv(SUMM, file.path(OUTDIR, "nce_targeting_summary.csv"), row.names = FALSE)

cat("\n===== BURDEN CAPTURED BY THE WORST-RANKED 20% OF DISTRICTS =====\n")
cat("(random selection captures 0.20 by definition; lift = capture / 0.20)\n\n")
print(as.data.frame(SUMM), row.names = FALSE)
cat("\n--- per-cell, model vs the two real alternatives ---\n")
wide <- CELL |> filter(estimand == "infill") |>
  select(country, outcome, arm, capture_top20) |>
  tidyr::pivot_wider(names_from = arm, values_from = capture_top20)
print(as.data.frame(wide |> mutate(across(where(is.numeric), ~ round(.x, 3)))),
      row.names = FALSE)
cat("\nDONE\n")
