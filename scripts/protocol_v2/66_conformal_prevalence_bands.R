# =============================================================================
# scripts/protocol_v2/66_conformal_prevalence_bands.R   [CP-01, 2026-09-27]
#
# CALIBRATED 90% BANDS FOR THE PLANNING PREVALENCE, AND CALIBRATED WHO
# EXCEEDANCE, BY SPLIT CONFORMAL ON OUT-OF-FOLD RESIDUALS.
#
# Design pre-registered in docs/findings/CP-01_CONFORMAL_PREVALENCE_BANDS_
# 2026-09-27.md. Per cell: the deployed quantity (calibrated index, IS-01,
# anchored on population to the national figure) is predicted out-of-fold for
# every surveyed district (5-fold by district, 10 draws, the protocol's
# in-fill folds), residuals against the survey's own district value give a
# conformal half-width and a conformal predictive distribution for the
# threshold chances. Checks: leave-one-out coverage of the 90% band, and the
# empirical coverage of the UE-01 STABILITY band on the same districts (the
# number that motivated this).
#
#   Rscript scripts/protocol_v2/66_conformal_prevalence_bands.R
# -> results/tables/protocol_v2/conformal_prev_cells.csv       per-cell summary and checks
#    results/tables/protocol_v2/conformal_prev_districts.csv   per-district audit trail
#    results/tables/protocol_v2/conformal_prev_residuals.csv   residual sets (read by builder 05)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
source("dashboard/data-raw/00_read_targets.R")
source("dashboard/R/fct_helpers.R")            # is_water()
set.seed(20260927L)

P2 <- "results/tables/protocol_v2"
DRAWS <- as.integer(Sys.getenv("CP_DRAWS", "10"))
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")   # the headline set, as deployed

TG  <- read_targets_with_extras(P2)
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
NE  <- { f <- "results/tables/national_estimates_all.csv"; if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
POP <- readRDS("dashboard/data/admin2_population.rds")
meta <- readRDS("dashboard/data/metadata.rds")
UE  <- { f <- "dashboard/data/uncertainty_ensembles.rds"; if (file.exists(f)) readRDS(f) else NULL }
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
S <- S[!is_water(S$Admin2), ]
LBL <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
KEY <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
.kk <- function(a, b) paste(trimws(a), trimws(b), sep = "|")
conf_q <- function(ae, level = 0.9) sort(ae)[min(length(ae), ceiling(level * (length(ae) + 1)))]

cells <- list(); dists <- list(); resids <- list()
for (ctry in unique(TG$country)) {
  all_s <- S[S$country == ctry, ]
  Xr <- prep_predictors_v2(as.matrix(all_s[, PREDS]))
  key_all <- .kk(all_s$Admin1, all_s$Admin2)
  pop <- POP[POP$country == LBL[[ctry]], ]
  for (oc in unique(TG$outcome[TG$country == ctry])) {
    t <- TG[TG$country == ctry & TG$outcome == oc & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
    tr0 <- match(.kk(t$Admin1, t$Admin2), key_all); keep <- is.finite(tr0); tr0 <- tr0[keep]; t <- t[keep, ]
    nS <- length(tr0); nA <- nrow(all_s)
    if (nS < 12) { cat(sprintf("   skip %s / %s: %d surveyed\n", ctry, oc, nS)); next }
    Y <- rep(NA_real_, nA); Y[tr0] <- .v2_logit(t$y_prev)
    y_nat <- rep(NA_real_, nA); y_nat[tr0] <- t$y_prev
    popn <- (if (startsWith(oc, "child_")) pop$pop_child else pop$pop_women)[match(key_all, .kk(pop$Admin1, pop$Admin2))]
    p_nat <- if (!is.null(NE)) NE$obs_prev[NE$country == LBL[[ctry]] & NE$outcome == oc][1] else NA_real_
    if (!is.finite(p_nat)) p_nat <- stats::weighted.mean(t$y_prev, t$n_raw)
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr0)   # basis is X-only; every district has X
    aux <- list(Admin1 = all_s$Admin1, y_nat = y_nat, target = "prev")
    anchor <- function(pred) {
      ok <- is.finite(popn) & popn > 0 & is.finite(pred)
      if (!any(ok)) return(pred)
      f <- function(c) sum(popn[ok] * .v2_expit(pred[ok] + c)) / sum(popn[ok]) - p_nat
      sh <- tryCatch(stats::uniroot(f, c(-12, 12))$root, error = function(e) 0)
      pred + sh
    }
    # out-of-fold prediction of the DEPLOYED quantity for every surveyed district
    acc <- cnt <- rep(0, nA)
    for (r in seq_len(DRAWS)) {
      folds <- make_folds_v2("kfold_district", nS, k = 5, rep_id = r)
      for (f in unique(folds)) {
        te <- tr0[folds == f]; tr <- tr0[folds != f]
        if (length(tr) < 8) next
        pred <- tryCatch(ARMS_V2[["domain_index_cal"]](tr, seq_len(nA), Y, NULL, D, aux), error = function(e) NULL)
        if (is.null(pred) || !any(is.finite(pred))) next
        pv <- .v2_expit(anchor(pred))
        acc[te] <- acc[te] + pv[te]; cnt[te] <- cnt[te] + 1
      }
    }
    ok <- tr0[cnt[tr0] > 0]
    p_hat <- acc[ok] / cnt[ok]
    p_svy <- t$y_prev[match(ok, tr0)]
    e <- p_svy - p_hat
    if (length(e) < 12) { cat(sprintf("   skip %s / %s: only %d usable residuals\n", ctry, oc, length(e))); next }
    half <- conf_q(abs(e))
    loo <- mean(vapply(seq_along(e), function(i) abs(e[i]) <= conf_q(abs(e[-i])), TRUE))
    # what the STABILITY band covered on the same districts (the premise)
    stab_cov <- NA_real_
    if (!is.null(UE) && !is.null(UE$cells)) {
      u <- UE$cells[UE$cells$country_key == KEY[[ctry]] & UE$cells$outcome == oc, ]
      j <- match(key_all[ok], .kk(u$Admin1, u$Admin2))
      if (any(is.finite(j))) stab_cov <- mean(p_svy >= u$prev_lo[j] & p_svy <= u$prev_hi[j], na.rm = TRUE)
    }
    th <- meta$who_thresholds[[oc]]
    cells[[length(cells) + 1]] <- data.frame(
      country = LBL[[ctry]], outcome = oc, n = length(e),
      half_width_pp = 100 * half, median_abs_err_pp = 100 * stats::median(abs(e)),
      loo_coverage_90 = loo, stability_band_coverage = stab_cov,
      th_moderate_plus = if (!is.null(th)) th[["mild"]] else NA_real_,
      th_severe = if (!is.null(th)) th[["moderate"]] else NA_real_, stringsAsFactors = FALSE)
    dists[[length(dists) + 1]] <- data.frame(
      country = LBL[[ctry]], outcome = oc, Admin1 = all_s$Admin1[ok], Admin2 = all_s$Admin2[ok],
      p_survey = p_svy, p_oof = p_hat, resid = e, stringsAsFactors = FALSE)
    resids[[length(resids) + 1]] <- data.frame(country = ctry, outcome = oc, e = e, stringsAsFactors = FALSE)
    cat(sprintf("   %-12s %-14s n %3d | half-width %5.1f pp | LOO cov %.2f | stability band covered %.2f\n",
                LBL[[ctry]], oc, length(e), 100 * half, loo, stab_cov))
  }
}
C <- bind_rows(cells); Dd <- bind_rows(dists); R <- bind_rows(resids)
write.csv(C, file.path(P2, "conformal_prev_cells.csv"), row.names = FALSE)
write.csv(Dd, file.path(P2, "conformal_prev_districts.csv"), row.names = FALSE)
write.csv(R, file.path(P2, "conformal_prev_residuals.csv"), row.names = FALSE)
cat(sprintf("\n%d cells | median half-width %.1f pp | mean LOO coverage %.2f | mean stability-band coverage %.2f\n",
            nrow(C), stats::median(C$half_width_pp), mean(C$loo_coverage_90), mean(C$stability_band_coverage, na.rm = TRUE)))
cat("written conformal_prev_{cells,districts,residuals}.csv\n")
