# =============================================================================
# explore/scripts/12_level_transport.R   [probe LT-01]
#
# QUESTION. The project scores transport on RANKS only, because biomarker
# levels carry large cross-survey offsets. That is a global rule, and RV-01
# found the problem is not global: the between-country share of level variance
# is 0.80 for child iron and 0.71 for women iron, but 0.14 for women B12 and
# 0.28-0.35 for vitamin A. Is a transported LEVEL - not just a ranking - usable
# for B12 and vitamin A?
#
# DESIGN. Leave-one-country-out on the RAW level, with the outcome NOT
# standardised within country (that standardisation is exactly what the rank-
# only protocol does, and it is what this probe is testing the need for).
# Predictors are still rank-normalised within country (fix 3), which is
# outcome-free and necessary for pooling.
#
# FIVE ARMS, TWO OF THEM ORACLES USED AS MEASURING DEVICES, NOT METHODS
#   null_train        predict the training countries' mean level. The honest
#                     no-information transported level.
#   index_cs          climate+soil index, centred and scaled on the training
#                     countries - so it inherits their offset, like any real
#                     transported prediction
#   blup_cs           the same information through a REML-shrunk kernel
#   null_true_mean    ORACLE: predict every district at the held-out country's
#                     TRUE mean. No ranking at all, perfect anchor. This is what
#                     a small anchoring survey would buy you.
#   index_cs_anchored ORACLE: index_cs shifted so its mean matches the held-out
#                     country's true mean. Perfect anchor AND ranking.
#
# The gap between index_cs and index_cs_anchored is the cost of the offset.
# The gap between null_true_mean and index_cs_anchored is what ranking adds once
# the anchor is right. Together they say whether the rank-only rule is the
# binding constraint for a given nutrient.
#
# METRICS. MAE in the outcome's own log units is not comparable across
# biomarkers, so everything is also reported relative to the held-out country's
# own within-country sd. skill = 1 - MAE/MAE(null_train): positive means the
# transported level beats predicting the training countries' mean.
#
#   Rscript explore/scripts/12_level_transport.R
# -> explore/out/12_level_transport.csv, 12_level_transport_decomp.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

E <- exp_load()
DOM <- E$domain_of
CS <- c("Climate and weather", "Soil characteristics")
ix <- exp_cell_index(E)

rows <- list()
for (on in unique(ix$outcome)) {
  cl <- exp_all_cells(E, "level", outcomes = on)
  if (length(cl) < 3) { message("skip ", on, " (", length(cl), " countries)"); next }

  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
  if (length(common) < 20) { message("skip ", on, " (", length(common), " common cols)"); next }

  # pooled, outcome on its RAW scale - the whole point of this probe
  Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  y    <- unlist(lapply(cl, function(z) z$y_nat))
  ctry <- rep(vapply(cl, function(z) z$country, ""), vapply(cl, function(z) z$n, 0L))
  wv   <- unlist(lapply(cl, function(z) z$w))
  lon  <- unlist(lapply(cl, function(z) z$aux$lon))
  lat  <- unlist(lapply(cl, function(z) z$aux$lat))
  csc  <- common[which(DOM[common] %in% CS)]
  if (length(csc) < 5) { message("skip ", on, " (no climate+soil)"); next }

  for (held in unique(ctry)) {
    te <- which(ctry == held); tr <- which(ctry != held)
    if (length(tr) < 20 || length(te) < 8) next

    aux <- list(lon = lon, lat = lat, Admin1 = paste(ctry, seq_along(ctry)),
                y_nat = y, target = "level", w = wv, country = ctry)
    # domain axes built on the pooled matrix with training-country orientation
    Dm <- domain_representation_v2(Xm[, csc, drop = FALSE], DOM, sign_rows = tr)

    p_null  <- rep(mean(y[tr]), length(te))
    p_index <- tryCatch(arm_domain_index_v2(tr, te, y, Xm, Dm, aux),
                        error = function(e) p_null)
    p_blup  <- tryCatch(reml_blup(y, list(cs = k_linear(Xm[, csc, drop = FALSE])),
                                  tr, te)$pred, error = function(e) p_null)
    truemu  <- mean(y[te])
    p_anch  <- p_index - mean(p_index) + truemu
    p_tmean <- rep(truemu, length(te))

    # SPREAD CALIBRATION. index_cs is scaled to the TRAINING countries' sd,
    # i.e. as if the ranking were perfect. A prediction whose correlation with
    # truth is rho should have spread rho x sd, not sd - otherwise its own
    # over-dispersion adds error. rho is estimated by leave-one-TRAINING-
    # country-out inside tr, so the held-out country never informs it.
    rho <- local({
      trc <- unique(ctry[tr]); oof <- rep(NA_real_, length(tr))
      for (c0 in trc) {
        i2 <- tr[ctry[tr] != c0]; o2 <- which(ctry[tr] == c0)
        if (length(i2) < 20 || !length(o2)) next
        D2 <- domain_representation_v2(Xm[, csc, drop = FALSE], DOM, sign_rows = i2)
        pp <- tryCatch(arm_domain_index_v2(i2, tr[o2], y, Xm, D2, aux),
                       error = function(e) rep(NA_real_, length(o2)))
        if (length(pp) == length(o2)) oof[o2] <- pp
      }
      r <- suppressWarnings(stats::cor(oof, y[tr], use = "complete.obs"))
      if (!is.finite(r)) 0 else max(min(r, 1), 0)
    })
    sd_te <- stats::sd(y[te])
    z_te <- (p_index - mean(p_index)) / max(stats::sd(p_index), 1e-9)
    # anchor from the oracle, spread shrunk by rho: training scale, then
    # held-out scale (the fuller oracle, isolating the ranking's own worth)
    p_cal_tr <- truemu + rho * stats::sd(y[tr]) * z_te
    p_cal_te <- truemu + rho * sd_te * z_te

    mae_null <- mean(abs(y[te] - p_null))

    for (a in c("null_train", "index_cs", "blup_cs", "null_true_mean",
                "index_cs_anchored", "anch_cal_trainsd", "anch_cal_truesd")) {
      p <- switch(a, null_train = p_null, index_cs = p_index, blup_cs = p_blup,
                  null_true_mean = p_tmean, index_cs_anchored = p_anch,
                  anch_cal_trainsd = p_cal_tr, anch_cal_truesd = p_cal_te)
      ok <- is.finite(p) & is.finite(y[te])
      if (sum(ok) < 5) next
      err <- y[te][ok] - p[ok]
      rows[[paste(on, held, a)]] <- data.frame(
        outcome = on, held_out = held, arm = a,
        n_areas = sum(ok), countries_trained = length(unique(ctry[tr])),
        mae = mean(abs(err)),
        mae_over_sd = mean(abs(err)) / sd_te,
        skill_vs_null = 1 - mean(abs(err)) / mae_null,
        bias = -mean(err),                       # predicted minus observed
        bias_over_sd = -mean(err) / sd_te,
        rmse = sqrt(mean(err^2)),
        offset_share_of_mse = (mean(err)^2) / mean(err^2),
        spearman = suppressWarnings(stats::cor(y[te][ok], p[ok], method = "spearman")),
        sd_heldout = sd_te, rho_nested = rho,
        stringsAsFactors = FALSE)
    }
    message(sprintf("  %-13s held out %-12s n=%3d  trained on %d countries",
                    on, held, length(te), length(unique(ctry[tr]))))
  }
}

R <- dplyr::bind_rows(rows)
exp_write(R, "12_level_transport")

# ── how much of the level error is the country offset? ──────────────────────
cat("\n== TRANSPORTED LEVEL: error relative to the held-out country's own sd ==\n")
cat("   (1.0 means the transported level is no better than a value one sd off)\n\n")
w <- reshape(R[R$arm %in% c("index_cs", "index_cs_anchored", "null_train",
                            "null_true_mean", "anch_cal_trainsd", "anch_cal_truesd"),
               c("outcome", "held_out", "arm", "mae_over_sd")],
             idvar = c("outcome", "held_out"), timevar = "arm", direction = "wide")
names(w) <- sub("^mae_over_sd[.]", "", names(w))
w$offset_cost <- round(w$index_cs - w$index_cs_anchored, 3)
w <- w[order(w$outcome, w$held_out), ]
print(w, row.names = FALSE, digits = 3)

cat("\n== by outcome: mean over held-out countries ==\n")
agg <- aggregate(cbind(index_cs, index_cs_anchored, anch_cal_trainsd,
                       anch_cal_truesd, null_train, null_true_mean,
                       offset_cost) ~ outcome, data = w,
                 FUN = function(z) round(mean(z, na.rm = TRUE), 3))
agg$countries <- aggregate(held_out ~ outcome, data = w, FUN = length)$held_out
print(agg[order(agg$index_cs), ], row.names = FALSE)

cat("\n== share of squared error that is pure country offset (index_cs arm) ==\n")
o <- R[R$arm == "index_cs", ]
s <- aggregate(offset_share_of_mse ~ outcome, data = o,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
s$mean_bias_over_sd <- aggregate(bias_over_sd ~ outcome, data = o,
                                 FUN = function(z) round(mean(abs(z), na.rm = TRUE), 3))$bias_over_sd
print(s[order(s$offset_share_of_mse), ], row.names = FALSE)
exp_write(s, "12_level_transport_decomp")

cat("\n== does the transported level beat the training-countries' mean? ==\n")
k <- R[R$arm %in% c("index_cs", "blup_cs"), ]
sk <- aggregate(skill_vs_null ~ outcome + arm, data = k,
                FUN = function(z) round(mean(z, na.rm = TRUE), 3))
sk$positive <- aggregate(skill_vs_null ~ outcome + arm, data = k,
                         FUN = function(z) sum(z > 0, na.rm = TRUE))$skill_vs_null
sk$n <- aggregate(skill_vs_null ~ outcome + arm, data = k, FUN = length)$skill_vs_null
print(sk[order(sk$arm, -sk$skill_vs_null), ], row.names = FALSE)

cat("
== THE CALIBRATION TEST ==
")
cat("   A constant at the TRUE mean gives MAE/sd = sqrt(2/pi) = 0.798 for a
")
cat("   normal outcome. Any arm BELOW that is adding real district information.

")
cc <- R[R$arm %in% c("null_true_mean", "index_cs_anchored", "anch_cal_trainsd",
                     "anch_cal_truesd"), ]
tab <- aggregate(mae_over_sd ~ outcome + arm, data = cc,
                 FUN = function(z) round(mean(z, na.rm = TRUE), 3))
tab2 <- reshape(tab, idvar = "outcome", timevar = "arm", direction = "wide")
names(tab2) <- sub("^mae_over_sd[.]", "", names(tab2))
tab2$beats_constant <- ifelse(tab2$anch_cal_truesd < tab2$null_true_mean, "YES", "no")
print(tab2[order(tab2$anch_cal_truesd), ], row.names = FALSE)
cat("
mean nested rho used for the spread shrinkage:
")
print(aggregate(rho_nested ~ outcome, data = R[R$arm == "index_cs", ],
                FUN = function(z) round(mean(z), 3)), row.names = FALSE)
