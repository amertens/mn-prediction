# =============================================================================
# scripts/protocol_v2/44_risk_category_accuracy.R   [RC-01]
#
# RISK-CATEGORY ACCURACY, REDONE UNDER THE PROTOCOL
#
# The January 2026 Ghana deck closed with a table of how often the model put a
# district in the right WHO public-health-significance class (none / mild /
# moderate / severe), exact and within one class, for child vitamin A at
# Admin-2 and Admin-1, full model vs proxies only (44/31, 75/49, 56/56,
# 100/69 percent). That came from the old leaderboard. This is the same
# question under protocol v2:
#   predictions   out-of-fold in-fill (5-fold by district, 10 draws) from the
#                 zero-tuning index, the spatial smoother, smoother + covariates
#                 and the jackknifed regional mean; plus leave-one-country-out
#                 transport with the held-out country's OWN national prevalence
#                 as the anchor (AR-01's design A1), which is the only honest
#                 way to put a transported ranking on the prevalence scale
#   classes       vitamin A: WHO 2009 bands (< 2 none, 2-10 mild, 10-20
#                 moderate, >= 20 severe) from metadata/who_severity_thresholds
#                 iron and the others: the WHO 20% threshold is the only one
#                 with a source, so a 2-class (above / below) accuracy is the
#                 primary, and anaemia-style 5 / 20 / 40 bands a labelled
#                 sensitivity with no WHO standing
#   scoring       exact-match and within-one-class share of districts, at
#                 Admin-2, and at Admin-1 after n_eff-weighted aggregation;
#                 the truth is the survey's own district / regional class
#
#   Rscript scripts/protocol_v2/44_risk_category_accuracy.R
# -> results/tables/protocol_v2/risk_category_accuracy.csv, _summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; REPS <- as.integer(Sys.getenv("RC_REPS", "10")); MIN_TRAIN <- 20L; set.seed(20260904L)
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) { b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
domains <- sort(unique(stats::na.omit(domain_of[PREDS]))); prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)]); CS <- c("Climate and weather", "Soil characteristics")
clamp <- function(p, eps = 0.005) pmin(pmax(p, eps), 1 - eps)

# ── classes ──────────────────────────────────────────────────────────────────
SCHEMES <- list(
  who_vitA   = list(outcomes = c("child_vitA", "women_vitA"), cuts = c(0.02, 0.10, 0.20), labels = c("none", "mild", "moderate", "severe"), source = "WHO 2009 vitamin A bands"),
  who20      = list(outcomes = c("child_iron", "women_iron", "women_folate", "women_b12", "child_zinc", "women_zinc"), cuts = 0.20, labels = c("below 20%", "20% or more"), source = "WHO / IZiNCG 20% threshold, 2 classes"),
  bands_5_20_40 = list(outcomes = c("child_iron", "women_iron", "women_folate", "women_b12", "child_zinc", "women_zinc"), cuts = c(0.05, 0.20, 0.40), labels = c("none", "mild", "moderate", "severe"), source = "anaemia-style bands, no WHO standing (sensitivity)"))
classify <- function(p, cuts) findInterval(p, cuts) + 1L
acc <- function(obs, pred) { ok <- is.finite(obs) & is.finite(pred); if (sum(ok) < 3) return(c(exact = NA_real_, within1 = NA_real_, n = sum(ok)))
  c(exact = mean(obs[ok] == pred[ok]), within1 = mean(abs(obs[ok] - pred[ok]) <= 1), n = sum(ok)) }

build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]; if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |> inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  m <- m[is.finite(m$lon), ]; if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 2) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, outcome = on, n = nrow(m), y = m$y_prev, w = m$n_eff, Admin1 = m$Admin1, Admin2 = m$Admin2, X = Xr, lon = m$lon, lat = m$lat)
}
agg1 <- function(x, w, g) { s <- tapply(x * w, g, sum) / tapply(w, g, sum); s[sort(names(s))] }
score_cell <- function(cl, pred_nat, scheme, arm, estimand, extra = list()) {
  cuts <- SCHEMES[[scheme]]$cuts
  o2 <- classify(cl$y, cuts); p2 <- classify(pred_nat, cuts); a2 <- acc(o2, p2)
  o1 <- classify(agg1(cl$y, cl$w, cl$Admin1), cuts); p1 <- classify(agg1(pred_nat, cl$w, cl$Admin1), cuts); a1 <- acc(o1, p1)
  out <- data.frame(country = cl$country, outcome = cl$outcome, estimand = estimand, arm = arm, scheme = scheme, rho_train = NA_real_,
                    n_admin2 = a2[["n"]], exact_admin2 = a2[["exact"]], within1_admin2 = a2[["within1"]], n_admin1 = a1[["n"]], exact_admin1 = a1[["exact"]], within1_admin1 = a1[["within1"]],
                    share_obs_top_admin2 = mean(o2 == max(seq_along(cuts) + 1L)), stringsAsFactors = FALSE)
  if (length(extra)) for (k in names(extra)) out[[k]] <- extra[[k]]
  out
}

rows <- list()
# ── in-fill ──────────────────────────────────────────────────────────────────
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells))) { cl <- tryCatch(build(cells$country[i], cells$outcome[i]), error = function(e) NULL); if (is.null(cl)) next
  Y <- .v2_logit(clamp(cl$y)); aux <- list(lon = cl$lon, lat = cl$lat, Admin1 = cl$Admin1, y_nat = Y)
  schemes <- names(SCHEMES)[vapply(SCHEMES, function(s) cl$outcome %in% s$outcomes, TRUE)]
  for (arm in c("domain_index", "spatial", "spatial_plus_domain", "region_mean_jk")) {
    P <- matrix(NA_real_, cl$n, REPS)
    for (r in seq_len(REPS)) { folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
      for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
        D <- domain_representation_v2(cl$X, domain_of, sign_rows = tr)
        p <- tryCatch(ARMS_V2[[arm]](tr, te, Y, cl$X, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) P[te, r] <- p } }
    pred <- .v2_expit(rowMeans(P, na.rm = TRUE))   # the draw-averaged out-of-fold prediction
    for (sc in schemes) rows[[length(rows) + 1L]] <- score_cell(cl, pred, sc, arm, "infill") }
  cat("infill done", cl$country, cl$outcome, "\n") }

# ── transport with the held-out country's own national anchor (A1) ──────────
for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(.v2_logit(clamp(z$y)))))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  for (h in names(cl)) { te <- which(ctry == h); tr <- which(ctry != h); if (length(tr) < MIN_TRAIN) next
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr); dom <- col_domain(colnames(Dm))
    Z <- cl[[h]]; p_nat <- sum(Z$y * Z$w) / sum(Z$w)
    sd_tr <- mean(vapply(setdiff(names(cl), h), function(t) stats::sd(.v2_logit(clamp(cl[[t]]$y))), numeric(1)))
    for (set in c("full", "climate_soil")) { keep <- if (set == "full") seq_len(ncol(Dm)) else which(dom %in% CS); if (!length(keep)) next
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm[, keep, drop = FALSE], aux), error = function(e) NULL); if (is.null(p) || length(p) != length(te) || stats::sd(p) == 0) next
      # nested rho among training countries for the BLP shrinkage, as in AR-01
      rt <- vapply(setdiff(names(cl), h), function(t) { te2 <- which(ctry == t); tr2 <- which(!ctry %in% c(h, t)); if (length(tr2) < MIN_TRAIN) return(NA_real_)
        D2 <- domain_representation_v2(Xm, domain_of, sign_rows = tr2); k2 <- if (set == "full") seq_len(ncol(D2)) else which(col_domain(colnames(D2)) %in% CS)
        p2 <- tryCatch(ARMS_V2[["domain_index"]](tr2, te2, Y, NULL, D2[, k2, drop = FALSE], aux), error = function(e) NULL); if (is.null(p2)) NA_real_ else suppressWarnings(stats::cor(cl[[t]]$y, p2, method = "spearman")) }, numeric(1))
      rho_tr <- max(0, mean(rt, na.rm = TRUE)); if (!is.finite(rho_tr)) rho_tr <- 0
      z <- as.numeric(scale(p)); pred <- .v2_expit(.v2_logit(clamp(p_nat)) + rho_tr * sd_tr * z)
      schemes <- names(SCHEMES)[vapply(SCHEMES, function(s) on %in% s$outcomes, TRUE)]
      for (sc in schemes) rows[[length(rows) + 1L]] <- score_cell(Z, pred, sc, paste0("transport_anchored_", set), "country", list(rho_train = rho_tr)) } }
  cat("transport done", on, "\n") }

R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "risk_category_accuracy.csv"), row.names = FALSE)
SUMM <- R |> group_by(estimand, scheme, arm) |> summarise(cells = dplyr::n(), exact_admin2 = mean(exact_admin2, na.rm = TRUE), within1_admin2 = mean(within1_admin2, na.rm = TRUE),
  exact_admin1 = mean(exact_admin1, na.rm = TRUE), within1_admin1 = mean(within1_admin1, na.rm = TRUE), .groups = "drop")
write.csv(SUMM, file.path(OUTDIR, "risk_category_accuracy_summary.csv"), row.names = FALSE)
cat("\n===== RC-01: WHO class accuracy (share of units), mean over cells =====\n")
for (sc in names(SCHEMES)) { cat(sprintf("\n-- %s [%s] --\n", sc, SCHEMES[[sc]]$source))
  print(as.data.frame(SUMM[SUMM$scheme == sc, ] |> mutate(across(where(is.numeric), ~ round(.x, 2))) |> arrange(estimand, desc(exact_admin2))), row.names = FALSE) }
cat("\n-- child vitamin A, in-fill, per country (the January table's cell) --\n")
print(as.data.frame(R |> filter(scheme == "who_vitA", outcome == "child_vitA", estimand == "infill") |> select(country, arm, exact_admin2, within1_admin2, exact_admin1, within1_admin1) |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
cat("\n-- prevalence of the top class among districts (how often 'severe' is the truth) --\n")
print(as.data.frame(R |> filter(estimand == "infill", arm == "domain_index") |> group_by(scheme, outcome) |> summarise(share_top = round(mean(share_obs_top_admin2), 2), .groups = "drop")), row.names = FALSE)
cat("\nDONE\n")
