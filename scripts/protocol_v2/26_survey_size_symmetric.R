# =============================================================================
# scripts/protocol_v2/26_survey_size_symmetric.R   [G4-02]
#
# DOES THE MODEL LET A COUNTRY FIELD A SMALLER SURVEY? (done symmetrically)
#
# Script 17 answered half of this. Its survey arm was fitted to a simulated
# smaller survey, but its model arm was trained on the FULL survey's district
# outcomes, so the comparison across sample fractions was not fair to the
# survey: the model got information the survey arm was denied. This script
# subsamples BOTH. For every fraction f and replicate:
#
#   y_sub    the full-survey district estimate plus the incremental sampling
#            noise a fraction-f survey would add: var = (1/f - 1) p(1-p)/n_eff
#            (zero at f = 1, so the curve anchors on the full survey; validated
#            against the empirical ws5 curve in script 17)
#   survey   jackknifed regional mean of y_sub  (flat within region)
#   model    the in-fill model FITTED ON y_sub, predicted out of fold (5-fold
#            district CV) -- the model sees only the smaller survey
#   aug      survey level + the model's within-region pattern (both from y_sub)
#
# All scored against the FULL survey's district estimates, the best available
# stand-in for truth. If the model's out-of-fold accuracy degrades more slowly
# than the survey's as f falls, the model substitutes for sample; if it
# degrades in step, it does not. Both the zero-tuning index and the spatial +
# covariate arm are run (G4_SPATIAL=0 skips the slower one).
#
#   Rscript scripts/protocol_v2/26_survey_size_symmetric.R
# -> results/tables/protocol_v2/survey_size_symmetric.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"
FRACTIONS <- c(0.15, 0.25, 0.40, 0.60, 0.80, 1.00)
REPS <- as.integer(Sys.getenv("G4_REPS", "8")); DO_SP <- Sys.getenv("G4_SPATIAL", "1") == "1"
set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2],
             stringsAsFactors = FALSE) }))
build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]; if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 2) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(n = nrow(m), y = m$y_prev, n_eff = m$n_eff, Admin1 = m$Admin1, X = Xr,
       D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = m$y_prev))
}
oof <- function(cl, ysub_logit, arm, rep_id) {
  folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = rep_id); pred <- rep(NA_real_, cl$n)
  aux <- cl$aux; aux$y_nat <- ysub_logit
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
    p <- tryCatch(ARMS_V2[[arm]](tr, te, ysub_logit, cl$X, cl$D, aux), error = function(e) rep(NA_real_, length(te)))
    if (length(p) == length(te)) pred[te] <- p
  }
  pred
}

rows <- list()
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL); if (is.null(cl)) next
  reg <- factor(cl$Admin1); pv <- pmin(pmax(cl$y, 1e-4), 1 - 1e-4)
  for (f in FRACTIONS) for (r in seq_len(REPS)) {
    nn <- pmax(round(f * cl$n_eff), 1)
    var_extra <- pmax((1 / f - 1) * pv * (1 - pv) / cl$n_eff, 0)
    y_sub <- pmin(pmax(cl$y + stats::rnorm(cl$n, 0, sqrt(var_extra)), 0), 1)
    ys_logit <- .v2_logit(y_sub)
    rs <- tapply(y_sub * nn, reg, sum)[as.character(reg)]; rn <- tapply(nn, reg, sum)[as.character(reg)]
    denom <- rn - nn
    rmean <- ifelse(denom > 0, (rs - y_sub * nn) / denom, stats::weighted.mean(y_sub, nn))
    p <- list(survey = rmean)
    mh <- oof(cl, ys_logit, "domain_index", 100 * r + round(100 * f))
    p$model_idx <- .v2_expit(mh)
    mreg <- tapply(mh, reg, mean)[as.character(reg)]
    p$aug_idx <- .v2_expit(.v2_logit(pmin(pmax(rmean, 1e-4), 1 - 1e-4)) + (mh - mreg))
    if (DO_SP) p$model_sp <- .v2_expit(oof(cl, ys_logit, "spatial_plus_domain", 100 * r + round(100 * f)))
    for (nm in names(p)) {
      ok <- is.finite(p[[nm]]) & is.finite(cl$y); if (!any(ok)) next
      rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, fraction = f, rep = r, arm = nm,
        n_areas = cl$n, mae = 100 * mean(abs(cl$y[ok] - p[[nm]][ok])),
        rho = suppressWarnings(stats::cor(cl$y[ok], p[[nm]][ok], method = "spearman")), stringsAsFactors = FALSE)
    }
  }
  cat("done", cn, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "survey_size_symmetric.csv"), row.names = FALSE)
CELL <- R |> group_by(country, outcome, fraction, arm) |> summarise(mae = mean(mae), rho = mean(rho, na.rm = TRUE), .groups = "drop")
SUMM <- CELL |> group_by(fraction, arm) |> summarise(cells = dplyr::n(), mae = round(median(mae), 2), rho = round(median(rho, na.rm = TRUE), 3), .groups = "drop")
write.csv(SUMM, file.path(OUTDIR, "survey_size_symmetric_summary.csv"), row.names = FALSE)

cat("\n===== G4-02: BOTH arms see only the smaller survey; scored vs the full survey =====\n")
cat("median district MAE (pp):\n"); print(as.data.frame(pivot_wider(SUMM[, c("fraction", "arm", "mae")], names_from = arm, values_from = mae)), row.names = FALSE)
cat("\nmedian Spearman:\n"); print(as.data.frame(pivot_wider(SUMM[, c("fraction", "arm", "rho")], names_from = arm, values_from = rho)), row.names = FALSE)
cat("\n--- head-to-head vs survey-only, per cell x fraction (identical draws) ---\n")
W <- pivot_wider(CELL[, c("country", "outcome", "fraction", "arm", "mae")], names_from = arm, values_from = mae)
for (a in setdiff(unique(CELL$arm), "survey")) { d <- W[[a]] - W$survey
  cat(sprintf("  %-10s lower MAE in %3d of %3d | median %+.2f pp\n", a, sum(d < 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE))) }
cat("\n--- degradation from f = 1.00 to f = 0.15 (median MAE, pp) ---\n")
for (a in unique(SUMM$arm)) { m1 <- SUMM$mae[SUMM$arm == a & SUMM$fraction == 1]; m0 <- SUMM$mae[SUMM$arm == a & SUMM$fraction == 0.15]
  cat(sprintf("  %-10s %.2f -> %.2f  (+%.2f)\n", a, m1, m0, m0 - m1)) }
cat("\nDONE\n")
