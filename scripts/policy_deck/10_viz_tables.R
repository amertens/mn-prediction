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
#   Rscript scripts/policy_deck/10_viz_tables.R [A,B,C,D,E]   (E: Ghana source-map PNGs for the deck)
# -> results/tables/policy_deck/viz/*.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R")
P2 <- "results/tables/protocol_v2"; OUT <- "results/tables/policy_deck/viz"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BLOCKS <- if (length(commandArgs(TRUE))) strsplit(commandArgs(TRUE)[1], ",")[[1]] else c("A", "B", "C", "D", "E")
set.seed(20260917L)
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi")
# optional overrides (v4 MNF15 deck, 27 Sep): VZ_A_OUTCOMES / VZ_A_COUNTRIES for block A and
# VZ_B2_PAIRS = "Country:outcome:file.csv;..." for block B2. Unset, the script does what it always did.
env_list <- function(v, sep = ",") { x <- Sys.getenv(v); if (nzchar(x)) strsplit(x, sep)[[1]] else NULL }
if (!is.null(env_list("VZ_A_COUNTRIES"))) COUNTRIES <- env_list("VZ_A_COUNTRIES")

# an urbanicity composite: mean of within-country ranks of night lights, population density, built surface, urban land cover
urb_cols <- intersect(c("ntl_ccnl", "wpop_log_density_survey_year", "built_surface", "lcover_urban_frac_t0"), names(S))
S$urbanicity <- ave(seq_len(nrow(S)), S$country, FUN = function(i) { m <- sapply(urb_cols, function(cc) rank(S[[cc]][i], na.last = "keep") / sum(is.finite(S[[cc]][i]))); rowMeans(m, na.rm = TRUE) })

# ── A. in-fill out-of-fold predictions, child iron ────────────────────────────
if ("A" %in% BLOCKS) for (OUTC in (if (!is.null(env_list("VZ_A_OUTCOMES"))) env_list("VZ_A_OUTCOMES") else c("child_iron", "women_iron"))) {
  cat("A. out-of-fold", OUTC, "\n"); rows <- list()
  for (cn in COUNTRIES) {
    t <- TG[TG$country == cn & TG$outcome == OUTC & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
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
  write.csv(bind_rows(rows), file.path(OUT, paste0("oof_", OUTC, ".csv")), row.names = FALSE)
}

# ── B. deployment fit with bootstrap refits: exceedance and prediction intervals ─
# The index fitted on all surveyed districts of a country and applied to every district (as on the dashboard),
# refitted on 200 stratified resamples of the surveyed districts: median prediction, 90% band, P(prevalence >= 20%).
if ("B" %in% BLOCKS) for (pair in list(c("Ghana", "child_vitA"), c("Ghana", "women_iron"), c("Malawi", "women_iron"))) {
  CN <- pair[1]; ON <- pair[2]; cat("B. deployment bootstrap", CN, ON, "\n"); B <- 200L
  t <- TG[TG$country == CN & TG$outcome == ON & is.finite(TG$y_prev), c("Admin1", "Admin2", "y_prev", "n_psu")]
  m <- left_join(S[S$country == CN, c("Admin1", "Admin2", PREDS)], t, by = c("Admin1", "Admin2"))
  tr0 <- which(is.finite(m$y_prev)); Xr <- prep_predictors_v2(as.matrix(m[, PREDS])); Y <- .v2_logit(m$y_prev); aux <- list(Admin1 = m$Admin1, y_nat = Y)
  P <- matrix(NA_real_, nrow(m), B)
  for (b in seq_len(B)) { tr <- unlist(lapply(split(tr0, m$Admin1[tr0]), function(i) i[sample.int(length(i), length(i), replace = TRUE)]))   # not sample(i): a one-district region would draw from 1:i
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
    P[, b] <- tryCatch(.v2_expit(ARMS_V2[["domain_index"]](tr, seq_len(nrow(m)), Y, NULL, D, aux)), error = function(e) NA_real_)
    if (b %% 50 == 0) cat("  ", b, "\n") }
  ok <- colSums(is.finite(P)) == nrow(m); P <- P[, ok, drop = FALSE]
  D0 <- domain_representation_v2(Xr, domain_of, sign_rows = tr0); fit <- .v2_expit(ARMS_V2[["domain_index"]](tr0, seq_len(nrow(m)), Y, NULL, D0, aux))
  E <- data.frame(country = CN, outcome = ON, Admin1 = m$Admin1, Admin2 = m$Admin2, surveyed = is.finite(m$y_prev), y_prev = m$y_prev, n_psu = m$n_psu, pred_fit = fit,
                  pred_med = apply(P, 1, median), pred_lo = apply(P, 1, quantile, 0.05), pred_hi = apply(P, 1, quantile, 0.95),
                  p_ge20 = rowMeans(P >= 0.20), p_ge10 = rowMeans(P >= 0.10), refits = ncol(P), stringsAsFactors = FALSE)
  f <- if (CN == "Ghana" && ON == "child_vitA") "exceedance_ghana_vitA.csv" else sprintf("deploy_%s_%s.csv", tolower(CN), ON)
  write.csv(E, file.path(OUT, f), row.names = FALSE)
  cat(sprintf("  %d districts, %d refits; P(>=20%%) >= 0.8 in %d districts, <= 0.2 in %d\n", nrow(E), ncol(P), sum(E$p_ge20 >= 0.8), sum(E$p_ge20 <= 0.2)))
}

# ── B2. the CALIBRATED deployment tables (CP-01, 2026-09-27) ─────────────────
# The deck's exceedance and four-panel slides switched from block B's stability
# quantities to the calibrated ones: the deployed map (calibrated index +
# population anchor, exactly builder 05's fit) with the conformal 90% band and
# the conformal-predictive-distribution threshold chances from script 66's
# out-of-fold residuals. One committed csv per slide cell.
if ("B2" %in% BLOCKS) {
  source("R/admin2_key_hygiene.R")
  CFR  <- read.csv(file.path(P2, "conformal_prev_residuals.csv"), stringsAsFactors = FALSE)
  NE   <- read.csv("results/tables/national_estimates_all.csv", stringsAsFactors = FALSE)
  POP  <- readRDS("dashboard/data/admin2_population.rds")
  metaD <- readRDS("dashboard/data/metadata.rds")
  LBL2 <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
  B2_PAIRS <- if (!is.null(env_list("VZ_B2_PAIRS", ";"))) lapply(env_list("VZ_B2_PAIRS", ";"), function(z) strsplit(z, ":")[[1]]) else
    list(c("Ghana", "child_vitA", "exceedance_ghana_vitA_cal.csv"),
         c("Ghana", "women_iron", "deploy_ghana_women_iron_cal.csv"),
         c("Malawi", "women_iron", "deploy_malawi_women_iron_cal.csv"))
  for (pair in B2_PAIRS) {
    CN <- pair[1]; ON <- pair[2]; cat("B2. calibrated deployment", CN, ON, "\n")
    all_s <- S[S$country == CN & !is_water_admin2(S$Admin2), ]
    t <- TG[TG$country == CN & TG$outcome == ON & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
    key_all <- paste(all_s$Admin1, all_s$Admin2, sep = "|")
    tr <- match(paste(t$Admin1, t$Admin2, sep = "|"), key_all); keep <- is.finite(tr); tr <- tr[keep]; t <- t[keep, ]
    Xr <- prep_predictors_v2(as.matrix(all_s[, PREDS]))
    Y <- rep(NA_real_, nrow(all_s)); Y[tr] <- .v2_logit(t$y_prev)
    y_nat <- rep(NA_real_, nrow(all_s)); y_nat[tr] <- t$y_prev
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
    pred <- ARMS_V2[["domain_index_cal"]](tr, seq_len(nrow(all_s)), Y, NULL, D,
                                          list(Admin1 = all_s$Admin1, y_nat = y_nat, target = "prev"))
    pop <- POP[POP$country == LBL2[[CN]], ]
    popn <- (if (startsWith(ON, "child_")) pop$pop_child else pop$pop_women)[match(key_all, paste(pop$Admin1, pop$Admin2, sep = "|"))]
    p_nat <- NE$obs_prev[NE$country == LBL2[[CN]] & NE$outcome == ON][1]
    if (!is.finite(p_nat)) p_nat <- stats::weighted.mean(t$y_prev, t$n_raw)
    ok <- is.finite(popn) & popn > 0
    f <- function(c) sum(popn[ok] * .v2_expit(pred[ok] + c)) / sum(popn[ok]) - p_nat
    shift <- tryCatch(stats::uniroot(f, c(-12, 12))$root, error = function(e) 0)
    pa <- .v2_expit(pred + shift)
    e <- CFR$e[CFR$country == CN & CFR$outcome == ON]
    half <- sort(abs(e))[min(length(e), ceiling(0.9 * (length(e) + 1)))]
    th <- metaD$who_thresholds[[ON]]
    cps <- function(thv) vapply(pa, function(p) (sum(p + e >= thv) + 0.5) / (length(e) + 1), 0)
    E <- data.frame(country = CN, outcome = ON, Admin1 = all_s$Admin1, Admin2 = all_s$Admin2,
                    surveyed = seq_len(nrow(all_s)) %in% tr, y_prev = y_nat,
                    prev_anchored = pa,
                    prev_cal_lo = pmax(0, pa - half), prev_cal_hi = pmin(1, pa + half),
                    p_modplus_cal = if (!is.null(th)) cps(th[["mild"]]) else NA_real_,
                    p_sev_cal = if (!is.null(th)) cps(th[["moderate"]]) else NA_real_,
                    th_modplus = if (!is.null(th)) th[["mild"]] else NA_real_,
                    th_sev = if (!is.null(th)) th[["moderate"]] else NA_real_,
                    half_width = half, n_resid = length(e), stringsAsFactors = FALSE)
    write.csv(E, file.path(OUT, pair[3]), row.names = FALSE)
    cat(sprintf("  %d districts | half-width %.1f pp | P(>=%.0f%%)>=0.8 in %d, <=0.2 in %d\n",
                nrow(E), 100 * half, 100 * E$th_sev[1], sum(E$p_sev_cal >= 0.8), sum(E$p_sev_cal <= 0.2)))
  }
}

# ── C. coverage of the bootstrap rank interval under leave-one-country-out ───
if ("C" %in% BLOCKS) {
  cat("C. rank-interval coverage under LOCO (climate + soil, level)\n"); NB <- 100L
  CS <- drop_near_outcome_v2(intersect(MD$column[MD$domain %in% c("Climate and weather", "Soil characteristics")], names(S)), MD)
  ALL4 <- c("Gambia", "Ghana", "Malawi", "SierraLeone"); rows <- list(); drows <- list()
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
      lo <- apply(R, 1, quantile, 0.05); hi <- apply(R, 1, quantile, 0.95); med <- apply(R, 1, median)
      drows[[paste(on, h)]] <- data.frame(outcome = on, heldout = h, n_districts = length(te), truth_rank = truth, rank_med = med, rank_lo = lo, rank_hi = hi,
        score = abs(truth - med) / length(te), stringsAsFactors = FALSE)   # nonconformity: how far the truth sits from the point prediction, as a share of the list
      rows[[paste(on, h)]] <- data.frame(outcome = on, heldout = h, n_districts = length(te), refits = ncol(R), coverage_90 = mean(truth >= lo & truth <= hi),
        median_width = median(hi - lo), width_share = median(hi - lo) / length(te), spearman_med = cor(truth, apply(R, 1, median), method = "spearman"), stringsAsFactors = FALSE)
      cat(sprintf("  %-13s %-12s coverage %.2f width %.0f of %d\n", on, h, mean(truth >= lo & truth <= hi), median(hi - lo), length(te))) } }
  write.csv(bind_rows(rows), file.path(OUT, "rank_interval_coverage.csv"), row.names = FALSE)
  DR <- bind_rows(drows); write.csv(DR, file.path(OUT, "rank_interval_districts.csv"), row.names = FALSE)
  # conformal calibration: the 90th percentile of the nonconformity score over every held-out district is the
  # half-width (share of the list) a 90% interval needs to cover held-out survey ranks 90% of the time
  q90 <- quantile(DR$score, 0.9); cat(sprintf("  calibrated 90%% half-width: %.0f%% of the list (bootstrap median half-width %.0f%%)
", 100 * q90, 100 * median((DR$rank_hi - DR$rank_lo) / 2 / DR$n_districts)))
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
# ── E. Ghana maps of one representative layer per source group, as PNGs for the two "what each source is" slides ──
if ("E" %in% BLOCKS) {
  suppressPackageStartupMessages({library(ggplot2); library(sf); library(patchwork)})
  B <- sf::st_as_sf(readRDS("dashboard/data/admin2_boundaries.rds")[["ghana"]]); Sg <- S[S$country == "Ghana", ]
  sets <- list(land = c(clim_pr_ann_mean = "Rainfall, 30-year mean", soil_zinc_mean_0_20 = "Soil zinc", spam_share_cereals = "Cereal share of crops", glw_cattle_km2 = "Cattle per km2"),
               health = c(ihme_allanemia = "Modeled anemia (women)", map_sy_pf_parasite_rate = "Malaria parasite rate", mics_wealth_score_mean = "Household wealth (MICS)", ntl_ccnl = "Night-time lights"))
  for (nm in names(sets)) { vars <- sets[[nm]]; vars <- vars[names(vars) %in% names(Sg)]
    g <- dplyr::left_join(B, Sg[, c("Admin1", "Admin2", names(vars))], by = c("Admin1", "Admin2"))
    maps <- lapply(names(vars), function(v) ggplot(g) + geom_sf(aes(fill = .data[[v]]), colour = "white", linewidth = 0.08) + scale_fill_distiller(palette = "YlGnBu", direction = 1, na.value = "grey90", guide = "none") +
      labs(title = vars[[v]]) + theme_void(base_size = 11) + theme(plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5)))
    ggsave(sprintf("docs/slides/img/ghana_sources_%s.png", nm), patchwork::wrap_plots(maps, ncol = 2), width = 4.6, height = 5.2, dpi = 150, bg = "white")
    cat("  wrote ghana_sources_", nm, ".png\n", sep = "") }
}
cat("DONE\n")
