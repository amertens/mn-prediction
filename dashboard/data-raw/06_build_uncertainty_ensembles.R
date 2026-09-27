# =============================================================================
# dashboard/data-raw/06_build_uncertainty_ensembles.R   [UE-01, 2026-09-27]
#
# STABILITY ENSEMBLES FOR EVERY DISTRICT AND EVERY WEIGHT ON THE DASHBOARD
#
# Extends the two designs the deck already checked to the whole app:
#   - scripts/policy_deck/06_civ_rank_uncertainty.R (CV-01): resample the
#     training districts, refit, re-rank; and
#   - scripts/policy_deck/10_viz_tables.R block B: the same for the deployment
#     fit, with exceedance probabilities against the WHO bands.
# Per country x outcome cell: B stratified resamples (within Admin1, replace =
# TRUE) of the surveyed districts; per draw the domain-PC basis is rebuilt on
# the resampled rows (as CV-01 does) and the zero-tuning index refitted; every
# district of the country is re-ranked and re-mapped to a calibrated, anchored
# planning prevalence (IS-01 map with the cell's full-fit nested rho, held
# fixed across draws; the anchor shift is re-solved per draw).
#
# WHAT THESE RANGES ARE (the deck's VZ-01 language): stability under a
# different training draw, NOT calibrated coverage of the survey's own rank —
# under leave-one-country-out the 90% band covers the held-out survey rank 38%
# of the time (results/tables/policy_deck/viz/rank_interval_coverage.csv). The
# app says this wherever they appear. The out-of-fold worst-fifth probability
# (scripts/policy_deck/07, calibration checked) remains the honest "how sure"
# for surveyed districts; these ensembles add the unsurveyed districts, the
# planning-prevalence bands and the WHO-threshold exceedance, all labelled as
# stability.
#
# Also refits the POOLED four-country index (the importance tab's object) on
# the same resampling design, so every back-projected weight carries a range.
#
#   Rscript dashboard/data-raw/06_build_uncertainty_ensembles.R
#   UE_DRAWS (default 200), UE_DRAWS_POOLED (default 150) trim or extend.
# -> dashboard/data/uncertainty_ensembles.rds
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
setwd(here::here())
source("R/protocol_v2.R")
source("R/protocol_v2_weights.R")
source("R/protocol_v2_importance.R")
source("dashboard/data-raw/00_read_targets.R")
set.seed(20260927L)

P2  <- "results/tables/protocol_v2"
OUT <- "dashboard/data"
B_CELL   <- as.integer(Sys.getenv("UE_DRAWS", "200"))
B_POOLED <- as.integer(Sys.getenv("UE_DRAWS_POOLED", "150"))
TIERS <- "open,survey_public"                       # the headline set (TP-01)
Sys.setenv(V2_PREDICTOR_TIERS = TIERS)

TG  <- read_targets_with_extras(P2)
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
NE  <- { f <- "results/tables/national_estimates_all.csv"; if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
POP <- readRDS(file.path(OUT, "admin2_population.rds"))
meta <- readRDS(file.path(OUT, "metadata.rds"))
source("dashboard/R/fct_helpers.R")                 # is_water()
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
S <- S[!is_water(S$Admin2), ]
LABEL <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
KEY   <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
.k <- function(a, b) paste(trimws(a), trimws(b), sep = "|")
qs <- function(M, p) apply(M, 1, stats::quantile, p, na.rm = TRUE)
cat(sprintf("tiers %s: %d predictors; %d cell draws, %d pooled draws\n", TIERS, length(PREDS), B_CELL, B_POOLED))

# ── A. per-cell deployment ensembles ─────────────────────────────────────────
cells <- list(); cell_w <- list(); cell_dom <- list()
for (ctry in unique(TG$country)) {
  all_s <- S[S$country == ctry, ]
  Xr <- prep_predictors_v2(as.matrix(all_s[, PREDS]))
  cn_cols <- colnames(Xr)                    # prep drops constant / empty columns per country
  key_all <- .k(all_s$Admin1, all_s$Admin2)
  pop <- POP[POP$country == LABEL[[ctry]], ]
  for (oc in unique(TG$outcome[TG$country == ctry])) {
    t <- TG[TG$country == ctry & TG$outcome == oc & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
    tr0 <- match(.k(t$Admin1, t$Admin2), key_all); keep <- is.finite(tr0); tr0 <- tr0[keep]; t <- t[keep, ]
    n <- nrow(all_s); if (length(tr0) < 8) next
    Y <- rep(NA_real_, n); Y[tr0] <- .v2_logit(t$y_prev)
    popn <- (if (startsWith(oc, "child_")) pop$pop_child else pop$pop_women)[match(key_all, .k(pop$Admin1, pop$Admin2))]
    p_nat <- if (!is.null(NE)) NE$obs_prev[NE$country == LABEL[[ctry]] & NE$outcome == oc][1] else NA_real_
    if (!is.finite(p_nat)) p_nat <- stats::weighted.mean(t$y_prev, t$n_raw)
    D0 <- domain_representation_v2(Xr, domain_of, sign_rows = tr0)
    rho0 <- .index_rho_v2(tr0, Y, D0, how = "nested")
    th <- meta$who_thresholds[[oc]]
    anchor <- function(pred) {                       # shift so pop-weighted mean prevalence = the national figure
      ok <- is.finite(popn) & popn > 0 & is.finite(pred)
      if (!any(ok)) return(pred)
      f <- function(c) sum(popn[ok] * .v2_expit(pred[ok] + c)) / sum(popn[ok]) - p_nat
      sh <- tryCatch(stats::uniroot(f, c(-12, 12))$root, error = function(e) 0)
      pred + sh
    }
    RK <- matrix(NA_real_, n, B_CELL); PV <- matrix(NA_real_, n, B_CELL)
    BT <- matrix(NA_real_, length(cn_cols), B_CELL, dimnames = list(cn_cols, NULL))
    strata <- split(tr0, all_s$Admin1[tr0])
    for (b in seq_len(B_CELL)) {
      tr <- unlist(lapply(strata, function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
      D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
      w <- tryCatch(.ws_z_pooled(tr, Y, D), error = function(e) NULL)
      if (is.null(w)) next
      idx <- as.numeric(D %*% w)
      if (!is.finite(stats::sd(idx[tr])) || stats::sd(idx[tr]) == 0) next
      RK[, b] <- rank(-idx, ties.method = "average")
      z <- (idx - mean(idx[tr])) / stats::sd(idx[tr])
      PV[, b] <- .v2_expit(anchor(mean(Y[tr]) + max(rho0, 0.001) * stats::sd(Y[tr]) * z))
      BT[, b] <- index_backproject_v2(w, attr(D, "basis"), cn_cols)
    }
    ok <- colSums(is.finite(RK)) == n; RK <- RK[, ok, drop = FALSE]; PV <- PV[, ok, drop = FALSE]; BT <- BT[, ok, drop = FALSE]
    if (ncol(RK) < 20) { cat(sprintf("   %-12s %-14s: only %d usable draws, skipped\n", ctry, oc, ncol(RK))); next }
    kw <- ceiling(n / 5)
    cells[[length(cells) + 1]] <- data.frame(
      country = LABEL[[ctry]], country_key = KEY[[ctry]], outcome = oc,
      Admin1 = all_s$Admin1, Admin2 = all_s$Admin2, n_districts = n,
      rank_med = apply(RK, 1, stats::median), rank_lo = qs(RK, 0.05), rank_hi = qs(RK, 0.95),
      p_worst_fifth_boot = rowMeans(RK <= kw),
      prev_med = apply(PV, 1, stats::median), prev_lo = qs(PV, 0.05), prev_hi = qs(PV, 0.95),
      p_moderate_plus = if (!is.null(th)) rowMeans(PV >= th[["mild"]]) else NA_real_,
      p_severe = if (!is.null(th)) rowMeans(PV >= th[["moderate"]]) else NA_real_,
      th_moderate_plus = if (!is.null(th)) th[["mild"]] else NA_real_,
      th_severe = if (!is.null(th)) th[["moderate"]] else NA_real_,
      refits = ncol(RK), stringsAsFactors = FALSE)
    cells[[length(cells)]]$rank_width <- cells[[length(cells)]]$rank_hi - cells[[length(cells)]]$rank_lo
    cells[[length(cells)]]$width_share <- cells[[length(cells)]]$rank_width / n
    sgn <- pmax(rowMeans(BT > 0, na.rm = TRUE), rowMeans(BT < 0, na.rm = TRUE))
    cell_w[[length(cell_w) + 1]] <- data.frame(
      country_key = KEY[[ctry]], outcome = oc, column = cn_cols,
      beta_med = apply(BT, 1, stats::median, na.rm = TRUE), beta_lo = qs(BT, 0.05), beta_hi = qs(BT, 0.95),
      sign_stab = sgn, refits = ncol(BT), stringsAsFactors = FALSE)
    shr <- apply(BT, 2, function(bb) { a <- tapply(abs(bb), domain_of[cn_cols], sum, na.rm = TRUE); a / sum(a, na.rm = TRUE) })
    cell_dom[[length(cell_dom) + 1]] <- data.frame(
      country_key = KEY[[ctry]], outcome = oc, domain = rownames(shr),
      share_med = apply(shr, 1, stats::median, na.rm = TRUE), share_lo = qs(shr, 0.05), share_hi = qs(shr, 0.95),
      stringsAsFactors = FALSE)
    cat(sprintf("   %-12s %-14s %3d districts, %3d draws, rho %.2f, median rank width %.0f\n",
                LABEL[[ctry]], oc, n, ncol(RK), rho0, stats::median(cells[[length(cells)]]$rank_width)))
  }
}
cells <- bind_rows(cells); cell_w <- bind_rows(cell_w); cell_dom <- bind_rows(cell_dom)

# ── B. pooled-fit weight ensembles (the importance tab's object) ─────────────
pooled_w <- list(); pooled_dom <- list()
main_oc <- intersect(c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12", "child_zinc", "women_zinc"),
                     unique(TG$outcome))
for (tg_kind in c("level", "prev")) for (oc in main_oc) {
  cl <- list()
  for (cn in c("Gambia", "Ghana", "Malawi", "SierraLeone")) {
    ycol <- if (tg_kind == "level") "y_level" else "y_prev"
    t <- TG[TG$country == cn & TG$outcome == oc & is.finite(TG[[ycol]]), ]
    m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
    if (nrow(m) < 12) next
    y <- if (tg_kind == "level") as.numeric(scale(m$y_level)) else as.numeric(scale(.v2_logit(m$y_prev)))
    cl[[cn]] <- list(n = nrow(m), y = y, X = prep_predictors_v2(as.matrix(m[, PREDS])))
  }
  if (length(cl) < 2) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  Y <- unlist(lapply(cl, function(z) z$y))
  BT <- matrix(NA_real_, length(common), B_POOLED, dimnames = list(common, NULL))
  for (b in seq_len(B_POOLED)) {
    tr <- unlist(lapply(unique(ctry), function(g) { i <- which(ctry == g); i[sample.int(length(i), length(i), replace = TRUE)] }))
    D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    w <- tryCatch(.ws_z_pooled(tr, Y, D), error = function(e) NULL)
    if (is.null(w)) next
    BT[, b] <- index_backproject_v2(w, attr(D, "basis"), common)
  }
  BT <- BT[, colSums(is.finite(BT)) > 0, drop = FALSE]
  if (ncol(BT) < 20) next
  sgn <- pmax(rowMeans(BT > 0, na.rm = TRUE), rowMeans(BT < 0, na.rm = TRUE))
  pooled_w[[length(pooled_w) + 1]] <- data.frame(
    outcome = oc, target = tg_kind, column = common,
    beta_med = apply(BT, 1, stats::median, na.rm = TRUE), beta_lo = qs(BT, 0.05), beta_hi = qs(BT, 0.95),
    sign_stab = sgn, refits = ncol(BT), stringsAsFactors = FALSE)
  shr <- apply(BT, 2, function(bb) { a <- tapply(abs(bb), domain_of[common], sum, na.rm = TRUE); a / sum(a, na.rm = TRUE) })
  pooled_dom[[length(pooled_dom) + 1]] <- data.frame(
    outcome = oc, target = tg_kind, domain = rownames(shr),
    share_med = apply(shr, 1, stats::median, na.rm = TRUE), share_lo = qs(shr, 0.05), share_hi = qs(shr, 0.95),
    stringsAsFactors = FALSE)
  cat(sprintf("   pooled %-14s %-5s: %d columns, %d draws\n", oc, tg_kind, length(common), ncol(BT)))
}
pooled_w <- bind_rows(pooled_w); pooled_dom <- bind_rows(pooled_dom)

saveRDS(list(cells = cells, cell_weights = cell_w, cell_domains = cell_dom,
             pooled_weights = pooled_w, pooled_domains = pooled_dom,
             meta = list(draws_cell = B_CELL, draws_pooled = B_POOLED, tiers = TIERS,
                         note = paste("Stability under stratified resampling of the training districts (CV-01/VZ-01 design);",
                                      "basis and weights refitted per draw, IS-01 rho held at the full-fit value, anchor re-solved per draw.",
                                      "Not calibrated coverage: see rank_interval_coverage.csv (38% under LOCO).")),
             build_time = Sys.time()),
        file.path(OUT, "uncertainty_ensembles.rds"))
cat(sprintf("wrote uncertainty_ensembles.rds: %d district rows (%d cells), %d cell-weight rows, %d pooled-weight rows\n",
            nrow(cells), nrow(distinct(cells, country_key, outcome)), nrow(cell_w), nrow(pooled_w)))
