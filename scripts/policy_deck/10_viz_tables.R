# =============================================================================
# scripts/policy_deck/10_viz_tables.R   [VZ-01, 2026-09-17]
#
# TABLES BEHIND THE NEW VISUALISATION SLIDES OF THE FULL TALK
#
# Four small computations the deck's figures read, kept out of the qmd so a
# render stays quick and every number has a file behind it:
#
#   A. oof_child_iron.csv       in-fill (5-fold by district, 10 draws) index
#                               predictions of child iron deficiency for every
#                               surveyed district of The Gambia, Ghana, Malawi:
#                               mean, min, max over the draws, survey value,
#                               clusters, an urbanicity composite, Admin1.
#                               Feeds: predicted-vs-observed scatter, Admin-1
#                               paired maps, accuracy-by-urbanicity.
#   B. exceedance_ghana_vitA.csv deployment fit for Ghana child vitamin A on
#                               all 75 surveyed districts, 200 stratified
#                               bootstrap refits, applied to all 260 districts:
#                               P(predicted prevalence >= 20%) and >= 10%, the
#                               WHO vitamin A bands. Feeds: exceedance map.
#   C. rank_interval_coverage.csv  does the 90% rank interval of script 06's
#                               bootstrap cover the truth? For each held-out
#                               training country and outcome, 100 resamples of
#                               the other three (climate + soil, level target),
#                               the interval per district, and the share of
#                               districts whose survey rank falls inside it.
#                               Feeds: calibration slide.
#   D. contributions_top12.csv  per-district contributions beta_j * x_ij of the
#                               twelve largest back-projected weights of the
#                               pooled index (level target), by outcome.
#                               Feeds: the SHAP-style beeswarm.
#
#   Rscript scripts/policy_deck/10_viz_tables.R [A,B,C,D]
# -> results/tables/policy_deck/viz/*.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R")
P2 <- "results/tables/protocol_v2"; OUT <- "results/tables/policy_deck/viz"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BLOCKS <- if (length(commandArgs(TRUE))) strsplit(commandArgs(TRUE)[1], ",")[[1]] else c("A", "B", "C", "D")
set.seed(20260917L)
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi")

# an urbanicity composite: mean of within-country ranks of night lights, population density, built surface, urban land cover
urb_cols <- intersect(c("ntl_ccnl", "wpop_log_density_survey_year", "built_surface", "lcover_urban_frac_t0"), names(S))
S$urbanicity <- ave(seq_len(nrow(S)), S$country, FUN = function(i) { m <- sapply(urb_cols, function(cc) rank(S[[cc]][i], na.last = "keep") / sum(is.finite(S[[cc]][i]))); rowMeans(m, na.rm = TRUE) })

# ── A. in-fill out-of-fold predictions, child iron ────────────────────────────
if ("A" %in% BLOCKS) {
  cat("A. out-of-fold child iron\n"); rows <- list()
  for (cn in COUNTRIES) {
    t <- TG[TG$country == cn & TG$outcome == "child_iron" & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
    m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", "urbanicity", PREDS)], by = c("Admin1", "Admin2"))
    Xr <- prep_predictors_v2(as.matrix(m[, PREDS])); Y <- .v2_logit(m$y_prev); aux <- list(Admin1 = m$Admin1, y_nat = Y)
    pred <- matrix(NA_real_, nrow(m), 10)
    for (r in 1:10) { folds <- make_folds_v2("kfold_district", nrow(m), k = 5, rep_id = r)
      for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
        pred[te, r] <- ARMS_V2[["domain_index"]](tr, te, Y, NULL, D, aux) } }
    pp <- .v2_expit(pred)
    rows[[cn]] <- data.frame(country = cn, Admin1 = m$Admin1, Admin2 = m$Admin2, y_prev = m$y_prev, n_psu = m$n_psu, n_raw = m$n_raw, n_eff = m$n_eff, urbanicity = m$urbanicity,
                             pred = rowMeans(pp), pred_lo = apply(pp, 1, min), pred_hi = apply(pp, 1, max), stringsAsFactors = FALSE)
    cat(sprintf("  %-8s %d districts, Spearman %.2f\n", cn, nrow(m), cor(m$y_prev, rowMeans(pp), method = "spearman"))) }
  write.csv(bind_rows(rows), file.path(OUT, "oof_child_iron.csv"), row.names = FALSE)
}

# ── B. exceedance probabilities, Ghana child vitamin A ───────────────────────
if ("B" %in% BLOCKS) {
  cat("B. exceedance, Ghana child vitamin A\n"); B <- 200L
  t <- TG[TG$country == "Ghana" & TG$outcome == "child_vitA" & is.finite(TG$y_prev), c("Admin1", "Admin2", "y_prev", "n_psu")]
  m <- left_join(S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)], t, by = c("Admin1", "Admin2"))
  tr0 <- which(is.finite(m$y_prev)); Xr <- prep_predictors_v2(as.matrix(m[, PREDS])); Y <- .v2_logit(m$y_prev); aux <- list(Admin1 = m$Admin1, y_nat = Y)
  P <- matrix(NA_real_, nrow(m), B)
  for (b in seq_len(B)) { tr <- unlist(lapply(split(tr0, m$Admin1[tr0]), function(i) i[sample.int(length(i), length(i), replace = TRUE)]))   # not sample(i): a one-district region would draw from 1:i
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
    P[, b] <- tryCatch(.v2_expit(ARMS_V2[["domain_index"]](tr, seq_len(nrow(m)), Y, NULL, D, aux)), error = function(e) NA_real_)
    if (b %% 50 == 0) cat("  ", b, "\n") }
  ok <- colSums(is.finite(P)) == nrow(m); P <- P[, ok, drop = FALSE]
  E <- data.frame(country = "Ghana", outcome = "child_vitA", Admin1 = m$Admin1, Admin2 = m$Admin2, surveyed = is.finite(m$y_prev), y_prev = m$y_prev, n_psu = m$n_psu,
                  pred_med = apply(P, 1, median), pred_lo = apply(P, 1, quantile, 0.05), pred_hi = apply(P, 1, quantile, 0.95),
                  p_ge20 = rowMeans(P >= 0.20), p_ge10 = rowMeans(P >= 0.10), refits = ncol(P), stringsAsFactors = FALSE)
  write.csv(E, file.path(OUT, "exceedance_ghana_vitA.csv"), row.names = FALSE)
  cat(sprintf("  %d districts, %d refits; P(>=20%%) >= 0.8 in %d districts, <= 0.2 in %d\n", nrow(E), ncol(P), sum(E$p_ge20 >= 0.8), sum(E$p_ge20 <= 0.2)))
}

# ── C. coverage of the bootstrap rank interval under leave-one-country-out ───
if ("C" %in% BLOCKS) {
  cat("C. rank-interval coverage under LOCO (climate + soil, level)\n"); NB <- 100L
  CS <- drop_near_outcome_v2(intersect(MD$column[MD$domain %in% c("Climate and weather", "Soil characteristics")], names(S)), MD)
  ALL4 <- c("Gambia", "Ghana", "Malawi", "SierraLeone"); rows <- list()
  for (on in c("child_iron", "child_vitA", "women_iron", "women_vitA", "women_folate", "women_b12")) {
    cl <- list()
    for (cn in ALL4) { t <- TG[TG$country == cn & TG$outcome == on & is.finite(TG$y_level), ]
      m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", CS)], by = c("Admin1", "Admin2")); if (nrow(m) < 12) next
      cl[[cn]] <- list(n = nrow(m), y = m$y_level, X = prep_predictors_v2(as.matrix(m[, CS]))) }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE])); ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y)))); yobs <- unlist(lapply(cl, function(z) z$y))
    for (h in names(cl)) { te <- which(ctry == h); pool <- setdiff(names(cl), h); truth <- rank(-yobs[te])
      R <- matrix(NA_real_, length(te), NB)
      for (b in seq_len(NB)) { tr <- unlist(lapply(pool, function(g) { i <- which(ctry == g); i[sample.int(length(i), length(i), replace = TRUE)] }))
        D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        p <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, D, NULL), error = function(e) rep(NA_real_, length(te)))
        if (all(is.finite(p))) R[, b] <- rank(-p) }
      R <- R[, colSums(is.finite(R)) == length(te), drop = FALSE]; if (!ncol(R)) next
      lo <- apply(R, 1, quantile, 0.05); hi <- apply(R, 1, quantile, 0.95)
      rows[[paste(on, h)]] <- data.frame(outcome = on, heldout = h, n_districts = length(te), refits = ncol(R), coverage_90 = mean(truth >= lo & truth <= hi),
        median_width = median(hi - lo), width_share = median(hi - lo) / length(te), spearman_med = cor(truth, apply(R, 1, median), method = "spearman"), stringsAsFactors = FALSE)
      cat(sprintf("  %-13s %-12s coverage %.2f width %.0f of %d\n", on, h, mean(truth >= lo & truth <= hi), median(hi - lo), length(te))) } }
  write.csv(bind_rows(rows), file.path(OUT, "rank_interval_coverage.csv"), row.names = FALSE)
}

# ── D. per-district contributions of the top twelve predictors, pooled index ──
if ("D" %in% BLOCKS) {
  cat("D. contributions of the top twelve predictors\n")
  IT <- read.csv(file.path(P2, "index_importance_top.csv"), stringsAsFactors = FALSE); rows <- list()
  Xall <- S[, c("country", "Admin1", "Admin2")]
  for (cn in unique(S$country)) { i <- which(S$country == cn); Xr <- prep_predictors_v2(as.matrix(S[i, PREDS])); for (cc in colnames(Xr)) Xall[i, cc] <- Xr[, cc] }
  for (on in unique(IT$outcome)) {
    top <- IT |> filter(outcome == on, target == "level") |> arrange(desc(abs(beta_std))) |> head(12)
    tt <- TG[TG$outcome == on, c("country", "Admin1", "Admin2")]
    for (k in seq_len(nrow(top))) { cc <- top$column[k]; if (!cc %in% names(Xall)) next
      j <- inner_join(tt, Xall[, c("country", "Admin1", "Admin2", cc)], by = c("country", "Admin1", "Admin2"))
      rows[[paste(on, cc)]] <- data.frame(outcome = on, column = cc, rank = k, beta = top$beta[k], beta_std = top$beta_std[k], domain = top$domain[k],
        country = j$country, Admin2 = j$Admin2, x = j[[cc]], contribution = top$beta[k] * j[[cc]], stringsAsFactors = FALSE) } }
  write.csv(bind_rows(rows), file.path(OUT, "contributions_top12.csv"), row.names = FALSE)
}
cat("DONE\n")
