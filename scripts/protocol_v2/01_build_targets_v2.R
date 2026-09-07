# =============================================================================
# scripts/protocol_v2/01_build_targets_v2.R
#
# FIX 1 (part) and FIX 2: build the corrected district outcome table.
#
# WHAT THIS PRODUCES THAT svy_admin2_* DOES NOT
# ---------------------------------------------
#   n_eff        effective sample size = n_raw / deff, with ONE deff estimated
#                per country x outcome at the national level. The audit
#                confirmed WS1b defect 2 is still live: svy_prev is a
#                design-based weighted estimate but n_svy is dplyr::n(), the
#                raw unweighted count, and downstream code uses it as if it
#                were the estimator's sample size. It also confirmed that a
#                district-level deff is NOT estimable here - Ghana 64 of 75,
#                Malawi 74 of 87 and Gambia 17 of 30 districts contain a single
#                PSU, so svy_prev_se is degenerate (< 1e-10) for them.
#                Estimating deff nationally, where every country has 60-103
#                PSUs, is the available honest route.
#   y_level      FIX 2: the survey-weighted district mean of the NEGATED log
#                biomarker concentration, so that higher = worse status on both
#                scales and signs are comparable with the prevalence target.
#                Measured cross-validated skill on this target is +0.24 against
#                -0.03 for the dichotomised one under identical strict folds.
#   clamp flags  how many districts sit exactly on the logit clamp, which for
#                women's vitamin A is most of them (the target is exactly zero
#                in 47-87 percent of districts).
#
# WEIGHT DIAGNOSTIC (an audit item this script resolves rather than assumes)
# -------------------------------------------------------------------------
# Two weight problems were flagged and are measured here rather than silently
# "fixed", because changing the weight changes every published estimate and
# that is the PI's call, not this script's:
#   Ghana   the configured gw_sWeight takes only 3 distinct values, matching
#           gw_strata exactly - it is a stratum constant, not a survey weight,
#           while gw_PSU_weight / gw_PSUStrat_weight carry 90 distinct values
#           (one per cluster).
#   Gambia  gw_c_blood_weight exists and is unused; biomarker outcomes are
#           weighted by the general gw_svy_weight.
# The table reports district prevalences under the configured weight and under
# each alternative, with the correlation and mean absolute difference between
# them, so the size of the decision is visible.
#
#   Rscript scripts/protocol_v2/01_build_targets_v2.R
# -> results/tables/protocol_v2/targets_v2.csv
# -> results/tables/protocol_v2/deff_v2.csv
# -> results/tables/protocol_v2/weight_diagnostic_v2.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
# The whole R/ tree is sourced because the binary outcome must be resolved with
# the pipeline's OWN resolve_uniform_outcome(): compute_svy_admin2() overwrites
# the configured binary with one re-derived from the adjusted continuous
# biomarker under a uniform cross-country cutoff (BRINDA for vitamin A, a WHO
# threshold otherwise). Using the raw configured binary instead reproduces the
# pipeline exactly for 19 of 24 cells and diverges by up to 50 pp on the iron
# cells - caught by the reconciliation check at the foot of this script, which
# is why that check exists.
targets::tar_source("R")

STORE  <- "_targets_full"
OUTDIR <- "results/tables/protocol_v2"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")
# alternative weights to measure against the configured one
ALT_WEIGHTS <- list(
  gambia      = c("gw_c_blood_weight", "gw_hhweight", "gw_cWeight"),
  ghana       = c("gw_PSU_weight", "gw_PSUStrat_weight", "gw_cwt"),
  malawi      = character(0),
  sierraleone = c("gw_wStatWt")
)

cfgs <- get_country_configs()
tgt_rows <- list(); deff_rows <- list(); wt_rows <- list()

for (lc in names(COUNTRIES)) {
  cn <- COUNTRIES[[lc]]
  cc <- cfgs[[cn]]
  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on),
                                         store = STORE), error = function(e) NULL)
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on),
                                         store = STORE), error = function(e) NULL)
    if (is.null(sv) || is.null(od)) next
    d <- od$data
    need <- c("Admin1", "Admin2", oc$binary)
    if (!all(need %in% names(d))) next

    w  <- .v2_num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- NA
    psu <- if (!is.null(cc$cluster_id) && cc$cluster_id %in% names(d))
      as.character(d[[cc$cluster_id]]) else NA_character_

    # Use the pipeline's own uniform outcome, not the configured binary.
    derived <- tryCatch(resolve_uniform_outcome(d, cc, oc, label = "[v2]"),
                        error = function(e) NULL)
    ybin <- if (!is.null(derived)) .v2_num(derived) else .v2_num(d[[oc$binary]])
    outcome_source <- if (!is.null(derived)) "uniform_derived" else "configured_binary"

    # continuous, on the modelling scale, negated so higher = worse.
    # Vitamin A: the SAME uniformly BRINDA-adjusted RBP the prevalence is cut
    # from (AU-01 finding 2, 2026-09-07); before this the level used each
    # survey's own adjusted column, raw RBP in Malawi, so the two targets
    # carried different inflammation adjustments.
    ycont <- rep(NA_real_, nrow(d)); level_source <- "configured_continuous"
    adj_vita <- if (!is.null(oc$tag) && grepl("vitA", oc$tag, ignore.case = TRUE))
      tryCatch(brinda_vad_adjusted(d, cc, oc, label = "[v2 level]"), error = function(e) NULL) else NULL
    if (!is.null(adj_vita)) {
      v <- as.numeric(adj_vita); v[!is.finite(v) | v <= 0] <- NA
      ycont <- -log(v); level_source <- "brinda_adjusted_rbp"
    } else if (!is.null(oc$continuous) && oc$continuous %in% names(d)) {
      v <- .v2_num(d[[oc$continuous]])
      t <- if (identical(oc$cutoff_scale, "log")) v else {
        v[!is.finite(v) | v <= 0] <- NA; log(v)
      }
      ycont <- -t
    }

    # ── deff, estimated nationally where PSUs are plentiful ────────────────
    db <- deff_national_v2(ybin, w, psu)
    dc <- deff_national_v2(ycont, w, psu)
    deff_rows[[paste(lc, on)]] <- data.frame(
      country = cn, outcome = on,
      n_raw = db$n_raw, n_psu = db$n_psu, n_kish = round(db$n_kish, 1),
      deff_binary = round(db$deff, 3), deff_cont = round(dc$deff, 3),
      method = db$method, outcome_source = outcome_source, level_source = level_source,
      weight_col = cc$weight_col,
      weight_ndistinct = dplyr::n_distinct(round(w, 6)))

    # ── district aggregation, both targets ────────────────────────────────
    dd <- data.frame(Admin1 = as.character(d$Admin1),
                     Admin2 = as.character(d$Admin2),
                     ybin = ybin, ycont = ycont, w = w,
                     stringsAsFactors = FALSE)
    agg <- dd |>
      group_by(Admin1, Admin2) |>
      summarise(
        n_raw      = sum(is.finite(ybin)),
        n_raw_cont = sum(is.finite(ycont)),
        n_kish     = kish_n_v2(w[is.finite(ybin)]),
        y_prev     = .v2_wmean(ybin, w),
        y_level    = .v2_wmean(ycont, w),
        sd_level   = sqrt(.v2_wvar(ycont, w)),
        .groups = "drop")

    agg$n_eff      <- effective_n_v2(agg$n_raw,      db$deff)
    agg$n_eff_cont <- effective_n_v2(agg$n_raw_cont, dc$deff)
    agg$deff_binary <- db$deff
    agg$deff_cont   <- dc$deff
    agg$country <- cn; agg$outcome <- on
    agg$at_clamp <- is.finite(agg$y_prev) &
      (agg$y_prev <= 0.005 | agg$y_prev >= 0.995)
    # reconciliation against the pipeline's own table
    key <- paste(sv$Admin1, sv$Admin2)
    agg$svy_prev_pipeline <- sv$svy_prev[match(paste(agg$Admin1, agg$Admin2), key)]
    agg$n_svy_pipeline    <- sv$n_svy[match(paste(agg$Admin1, agg$Admin2), key)]
    tgt_rows[[paste(lc, on)]] <- agg

    # ── weight diagnostic ─────────────────────────────────────────────────
    for (aw in ALT_WEIGHTS[[lc]]) {
      if (!aw %in% names(d)) next
      w2 <- .v2_num(d[[aw]]); w2[!is.finite(w2) | w2 <= 0] <- NA
      if (all(is.na(w2))) next
      a2 <- dd |> mutate(w2 = w2) |> group_by(Admin1, Admin2) |>
        summarise(p2 = .v2_wmean(ybin, w2), .groups = "drop")
      m <- inner_join(agg[, c("Admin1", "Admin2", "y_prev")], a2,
                      by = c("Admin1", "Admin2"))
      ok <- is.finite(m$y_prev) & is.finite(m$p2)
      wt_rows[[paste(lc, on, aw)]] <- data.frame(
        country = cn, outcome = on,
        configured = cc$weight_col, alternative = aw,
        alt_ndistinct = dplyr::n_distinct(round(w2, 6)),
        n_areas = sum(ok),
        pearson = if (sum(ok) > 3 && sd(m$p2[ok]) > 0)
          round(cor(m$y_prev[ok], m$p2[ok]), 4) else NA_real_,
        spearman = if (sum(ok) > 3 && sd(m$p2[ok]) > 0)
          round(cor(m$y_prev[ok], m$p2[ok], method = "spearman"), 4) else NA_real_,
        mean_abs_diff_pp = round(mean(abs(m$y_prev[ok] - m$p2[ok])) * 100, 3),
        max_abs_diff_pp = round(max(abs(m$y_prev[ok] - m$p2[ok])) * 100, 3))
    }
    cat("built", lc, on, "\n")
  }
}

TG <- bind_rows(tgt_rows)
DF <- bind_rows(deff_rows)
WT <- bind_rows(wt_rows)
write.csv(TG, file.path(OUTDIR, "targets_v2.csv"), row.names = FALSE)
write.csv(DF, file.path(OUTDIR, "deff_v2.csv"), row.names = FALSE)
write.csv(WT, file.path(OUTDIR, "weight_diagnostic_v2.csv"), row.names = FALSE)

cat("\n=== design effects (national) and effective n ===\n")
print(as.data.frame(DF), row.names = FALSE)

cat("\n=== reconciliation with the pipeline's svy_prev ===\n")
rec <- TG |> filter(is.finite(y_prev), is.finite(svy_prev_pipeline)) |>
  group_by(country, outcome) |>
  summarise(n = n(), r = round(cor(y_prev, svy_prev_pipeline), 4),
            max_abs_diff_pp = round(max(abs(y_prev - svy_prev_pipeline)) * 100, 3),
            .groups = "drop")
print(as.data.frame(rec), row.names = FALSE)

cat("\n=== effective n vs raw n, and clamp incidence ===\n")
sm <- TG |> group_by(country, outcome) |>
  summarise(areas = n(),
            median_n_raw = median(n_raw),
            median_n_eff = round(median(n_eff), 1),
            areas_n_raw_lt10 = sum(n_raw < 10),
            areas_n_eff_lt10 = sum(n_eff < 10),
            at_clamp = sum(at_clamp), .groups = "drop")
print(as.data.frame(sm), row.names = FALSE)

cat("\n=== weight diagnostic (configured vs alternatives) ===\n")
print(as.data.frame(WT), row.names = FALSE)
cat("\nDONE\n")
