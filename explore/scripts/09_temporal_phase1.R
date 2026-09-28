# =============================================================================
# explore/scripts/09_temporal_phase1.R   [probe TM-01, phase 1]
#
# QUESTION. Fieldwork windows are only 2-3 months per country, so calendar
# spread is NOT the temporal axis available here. What varies across clusters
# is SEASONAL PHASE: two clusters visited in the same week sit at different
# points in their own agricultural year, because their rainy seasons peak in
# different months. Ferritin and retinol integrate months of intake, so where
# a household sits relative to its own last harvest should matter more than the
# long-run climatology of its pixel.
#
# The cluster table already carries fieldwork-window values (_fw, _fw_anom,
# _prev3) for 7 dynamic layers and a peak month (_peak) for each. It does NOT
# carry phase. This probe constructs it and asks whether it adds anything:
#
#   phase_<layer> = circular months from the layer's PEAK month to the
#                   cluster's fieldwork month, plus its sin/cos encoding
#
# ARMS (all ridge on rank-normalised blocks, folds cut by DISTRICT so a
# cluster's own district never trains on it)
#   clim          climatology only (_mean, _amp, _peak and the static layers)
#   clim_fw       + the fieldwork-window block (FW-01's +0.02)
#   clim_fw_phase + seasonal phase
#   phase_only    phase alone, to see whether it carries anything on its own
#
# WHAT THIS CANNOT TEST. A real lag stack (0-24 months before the draw, so the
# PREVIOUS growing season is represented) needs monthly rasters covering each
# survey's window. On disk that exists only for parts of 2014-15, so the lag
# hypothesis is deferred to the gated Earth Engine step. This probe's job is to
# say whether that extraction is worth the wall-clock.
#
#   Rscript explore/scripts/09_temporal_phase1.R
# -> explore/out/09_temporal_cluster.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))

TC <- read.csv(file.path(EXP_ROOT, "results/tables/cluster_level/targets_cluster.csv"),
               stringsAsFactors = FALSE)
PC <- read.csv(file.path(EXP_ROOT, "data/covariates/cluster/predictors_cluster.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
MD <- read.csv(file.path(EXP_ROOT, "data/covariates/cluster/predictors_cluster_metadata.csv"),
               stringsAsFactors = FALSE)

role_of <- stats::setNames(MD$role, MD$column)
FW_COLS   <- MD$column[MD$role == "fieldwork"]
DYN_COLS  <- MD$column[MD$role == "dynamic"]
BASE_COLS <- MD$column[MD$role %in% c("static", "slow", "dynamic")]

# ── construct seasonal phase ────────────────────────────────────────────────
# _peak columns hold the peak MONTH of a dynamic layer; month_med is the
# cluster's fieldwork month. The circular difference is months since that
# layer's seasonal peak at the moment of the blood draw.
peak_cols <- grep("_peak$", DYN_COLS, value = TRUE)
PH <- PC[, c("country", "cluster"), drop = FALSE]
for (pk in peak_cols) {
  d <- (as.numeric(PC$month_med) - as.numeric(PC[[pk]])) %% 12
  lay <- sub("_peak$", "", pk)
  PH[[paste0("phase_", lay)]]     <- d
  PH[[paste0("phase_sin_", lay)]] <- sin(2 * pi * d / 12)
  PH[[paste0("phase_cos_", lay)]] <- cos(2 * pi * d / 12)
}
PHASE_COLS <- setdiff(names(PH), c("country", "cluster"))
message("phase columns built: ", length(PHASE_COLS), " from ",
        length(peak_cols), " peak-month layers")

M0 <- merge(PC, PH, by = c("country", "cluster"))

# ── score one cell ──────────────────────────────────────────────────────────
ridge_pred <- function(tr, te, y, X) {
  if (ncol(X) < 2) return(rep(mean(y[tr]), length(te)))
  .v2_enet(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], alpha = 0)
}

BLOCKS <- list(
  clim          = BASE_COLS,
  clim_fw       = c(BASE_COLS, FW_COLS),
  clim_fw_phase = c(BASE_COLS, FW_COLS, PHASE_COLS),
  phase_only    = PHASE_COLS
)

rows <- list()
cells <- unique(TC[, c("country", "outcome")])
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  t <- TC[TC$country == cn & TC$outcome == on, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  m <- merge(t, M0[M0$country == cn, ], by = c("country", "cluster"))
  if (nrow(m) < 25 || dplyr::n_distinct(m$Admin2) < 5) next

  y <- m$y_level
  dist <- m$Admin2
  ud <- unique(dist)

  for (arm in names(BLOCKS)) {
    cc <- intersect(BLOCKS[[arm]], names(m))
    Xr <- prep_predictors_v2(as.matrix(m[, cc, drop = FALSE]))
    if (ncol(Xr) < 2) next
    for (r in seq_len(REPS)) {
      set.seed(20260951L + r)
      fold_of_dist <- stats::setNames(
        sample(rep(seq_len(min(5, length(ud))), length.out = length(ud))), ud)
      folds <- fold_of_dist[dist]
      pred <- rep(NA_real_, nrow(m))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 20 || !length(te)) next
        pred[te] <- ridge_pred(tr, te, y, Xr)
      }
      s_cl <- score_v2(y, pred, m$n_eff_cont, scale = "level")
      # and the same predictions aggregated to districts
      ok <- is.finite(pred)
      ad <- stats::aggregate(cbind(y, pred, w = m$n_eff_cont)[ok, ],
                             by = list(d = dist[ok]), FUN = mean)
      s_d <- score_v2(ad$y, ad$pred, ad$w, scale = "level")
      rows[[paste(cn, on, arm, r)]] <- rbind(
        cbind(data.frame(country = cn, outcome = on, arm = arm, rep = r,
                         unit = "cluster", n_clusters = nrow(m),
                         n_districts = length(ud), n_cols = ncol(Xr),
                         stringsAsFactors = FALSE), s_cl),
        cbind(data.frame(country = cn, outcome = on, arm = arm, rep = r,
                         unit = "district", n_clusters = nrow(m),
                         n_districts = length(ud), n_cols = ncol(Xr),
                         stringsAsFactors = FALSE), s_d))
    }
  }
  message("  ", cn, " ", on, "  clusters ", nrow(m), "  districts ", length(ud))
}

raw <- dplyr::bind_rows(rows)
sm <- raw |>
  dplyr::group_by(country, outcome, arm, unit) |>
  dplyr::summarise(reps = dplyr::n(), n_clusters = dplyr::first(n_clusters),
                   n_districts = dplyr::first(n_districts),
                   n_cols = dplyr::first(n_cols),
                   spearman = mean(spearman, na.rm = TRUE),
                   pearson = mean(pearson, na.rm = TRUE),
                   rmse_sd = mean(rmse_sd, na.rm = TRUE), .groups = "drop") |>
  as.data.frame()
exp_write(sm, "09_temporal_cluster")

cat("\n== median Spearman over cells, by arm and scoring unit ==\n")
a <- aggregate(spearman ~ unit + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$unit, -a$spearman), ], row.names = FALSE)

cat("\n== paired gain from each temporal block (district scoring) ==\n")
d <- sm[sm$unit == "district", ]
w <- reshape(d[, c("country", "outcome", "arm", "spearman")],
             idvar = c("country", "outcome"), timevar = "arm", direction = "wide")
names(w) <- sub("^spearman\\.", "", names(w))
w$fw_gain    <- round(w$clim_fw - w$clim, 3)
w$phase_gain <- round(w$clim_fw_phase - w$clim_fw, 3)
print(w[order(-w$phase_gain), ], row.names = FALSE, digits = 3)
cat("\nmean fieldwork-window gain:", round(mean(w$fw_gain, na.rm = TRUE), 4),
    " cells improved:", sum(w$fw_gain > 0, na.rm = TRUE), "of", sum(is.finite(w$fw_gain)), "\n")
cat("mean phase gain           :", round(mean(w$phase_gain, na.rm = TRUE), 4),
    " cells improved:", sum(w$phase_gain > 0, na.rm = TRUE), "of", sum(is.finite(w$phase_gain)), "\n")
