# =============================================================================
# explore/scripts/11_lag_stack.R   [probe TM-01, phase 2 scoring]
#
# QUESTION. Micronutrient stores integrate months of intake, so the quality of
# the PREVIOUS GROWING SEASON should predict status better than the long-run
# climatology of the pixel and better than the 3-month window the cluster table
# already carries. Phase 1 showed seasonal PHASE is a coin flip; this is the
# mechanistically stronger version of the same idea and the reason the Earth
# Engine extraction was run.
#
# INPUT. explore/out/10_gee_cluster_monthly.csv - 36 monthly values per cluster
# per layer (CHIRPS rainfall, MODIS NDVI and LST, FLDAS soil moisture and
# evaporation), aligned so lag 0 is each cluster's own fieldwork month.
#
# ANOMALIES, NOT LEVELS. Raw monthly values are dominated by the season the
# survey happened to run in, which is near-constant within a country. Each
# value is therefore converted to an anomaly against that cluster's OWN mean
# for the same calendar month over the three years available, so a lag column
# means "how good was that month, here, relative to normal here".
#
# FOUR REPRESENTATIONS OF THE SAME 36 LAGS
#   lag_raw      all 36 lags per layer (180 columns) - the unregularised
#                comparator, expected to lose
#   lag_spline   the lag profile projected onto a 5-df natural spline in lag,
#                so a layer costs 5 coefficients not 36: a penalised
#                distributed lag, the standard answer in environmental
#                epidemiology to exactly this shape of problem
#   lag_window   8 interpretable windows (0-2, 3-5, 6-8, 9-11, 12-14, 15-17,
#                18-23, 24-35 months before the draw)
#   prev_season  lags 9-15 only, one column per layer - the hypothesis stated
#                as narrowly as it can be
#
#   Rscript explore/scripts/11_lag_stack.R
# -> explore/out/11_lag_stack.csv, 11_lag_profile.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
suppressPackageStartupMessages({library(splines)})

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))

G <- read.csv(file.path(EXP_ROOT, "explore/out/10_gee_cluster_monthly.csv"),
              stringsAsFactors = FALSE)
TC <- read.csv(file.path(EXP_ROOT, "results/tables/cluster_level/targets_cluster.csv"),
               stringsAsFactors = FALSE)
PC <- read.csv(file.path(EXP_ROOT, "data/covariates/cluster/predictors_cluster.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
MD <- read.csv(file.path(EXP_ROOT, "data/covariates/cluster/predictors_cluster_metadata.csv"),
               stringsAsFactors = FALSE)
BASE_COLS <- MD$column[MD$role %in% c("static", "slow", "dynamic")]

G$cluster <- as.character(G$cluster)
PC$cluster <- as.character(PC$cluster)
TC$cluster <- as.character(TC$cluster)
G$cal_month <- as.integer(substr(G$year_month, 6, 7))
message("extraction: ", nrow(G), " rows, ", dplyr::n_distinct(G$layer), " layers, ",
        dplyr::n_distinct(paste(G$country, G$cluster)), " clusters")

# ── anomaly against each cluster's own same-calendar-month mean ──────────────
G <- G |>
  dplyr::group_by(country, cluster, layer, cal_month) |>
  dplyr::mutate(anom = value - mean(value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  as.data.frame()

LAYERS <- sort(unique(G$layer))
LAGS <- 0:35

#' Wide matrix: one column per (layer, lag), rows = clusters
wide_lags <- function(g) {
  key <- unique(g[, c("country", "cluster")])
  key <- key[order(key$country, key$cluster), ]
  id <- paste(key$country, key$cluster)
  out <- matrix(NA_real_, nrow(key), length(LAYERS) * length(LAGS))
  cn <- character(0)
  j <- 0L
  for (ly in LAYERS) for (lg in LAGS) {
    j <- j + 1L
    s <- g[g$layer == ly & g$lag_months == lg, ]
    out[, j] <- s$anom[match(id, paste(s$country, s$cluster))]
    cn <- c(cn, sprintf("%s_L%02d", ly, lg))
  }
  colnames(out) <- cn
  list(key = key, X = out)
}
W <- wide_lags(G)
message("lag matrix: ", nrow(W$X), " clusters x ", ncol(W$X), " columns")

# ── the four representations ────────────────────────────────────────────────
Bspl <- splines::ns(LAGS, df = 5)            # 36 x 5 basis in lag
WINDOWS <- list(w0_2 = 0:2, w3_5 = 3:5, w6_8 = 6:8, w9_11 = 9:11,
                w12_14 = 12:14, w15_17 = 15:17, w18_23 = 18:23, w24_35 = 24:35)

block_for <- function(kind) {
  cols <- list()
  for (ly in LAYERS) {
    idx <- match(sprintf("%s_L%02d", ly, LAGS), colnames(W$X))
    L <- W$X[, idx, drop = FALSE]
    L[!is.finite(L)] <- 0
    if (kind == "raw") {
      m <- L; colnames(m) <- sprintf("%s_L%02d", ly, LAGS)
    } else if (kind == "spline") {
      m <- L %*% Bspl; colnames(m) <- sprintf("%s_S%d", ly, seq_len(ncol(Bspl)))
    } else if (kind == "window") {
      m <- do.call(cbind, lapply(WINDOWS, function(k) rowMeans(L[, k + 1, drop = FALSE])))
      colnames(m) <- sprintf("%s_%s", ly, names(WINDOWS))
    } else {
      m <- matrix(rowMeans(L[, 9:15 + 1, drop = FALSE]), ncol = 1)
      colnames(m) <- sprintf("%s_prev_season", ly)
    }
    cols[[ly]] <- m
  }
  do.call(cbind, cols)
}
BLK <- lapply(c(raw = "raw", spline = "spline", window = "window",
                prev = "prev"), block_for)
for (nm in names(BLK)) message("  block ", nm, ": ", ncol(BLK[[nm]]), " columns")

KEYDF <- cbind(W$key, as.data.frame(do.call(cbind, BLK)))
names(KEYDF) <- make.unique(names(KEYDF))

# ── score, folds cut by district ────────────────────────────────────────────
M0 <- merge(PC, KEYDF, by = c("country", "cluster"))
ARMSETS <- list(
  clim            = BASE_COLS,
  clim_lag_raw    = c(BASE_COLS, colnames(BLK$raw)),
  clim_lag_spline = c(BASE_COLS, colnames(BLK$spline)),
  clim_lag_window = c(BASE_COLS, colnames(BLK$window)),
  clim_prev_seas  = c(BASE_COLS, colnames(BLK$prev)),
  lag_spline_only = colnames(BLK$spline),
  prev_seas_only  = colnames(BLK$prev)
)

rows <- list()
cells <- unique(TC[, c("country", "outcome")])
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  t <- TC[TC$country == cn & TC$outcome == on, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  m <- merge(t, M0[M0$country == cn, ], by = c("country", "cluster"))
  if (nrow(m) < 25 || dplyr::n_distinct(m$Admin2) < 5) next
  y <- m$y_level; dist <- m$Admin2; ud <- unique(dist)

  for (arm in names(ARMSETS)) {
    cc <- intersect(ARMSETS[[arm]], names(m))
    if (length(cc) < 2) next
    Xr <- prep_predictors_v2(as.matrix(m[, cc, drop = FALSE]))
    if (ncol(Xr) < 2) next
    for (r in seq_len(REPS)) {
      set.seed(20260951L + r)
      fod <- stats::setNames(
        sample(rep(seq_len(min(5, length(ud))), length.out = length(ud))), ud)
      folds <- fod[dist]
      pred <- rep(NA_real_, nrow(m))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 20 || !length(te)) next
        pred[te] <- .v2_enet(Xr[tr, , drop = FALSE], y[tr], Xr[te, , drop = FALSE], alpha = 0)
      }
      ok <- is.finite(pred)
      ad <- stats::aggregate(cbind(y, pred, w = m$n_eff_cont)[ok, ],
                             by = list(d = dist[ok]), FUN = mean)
      rows[[paste(cn, on, arm, r)]] <- rbind(
        cbind(data.frame(country = cn, outcome = on, arm = arm, rep = r,
                         unit = "cluster", n_cols = ncol(Xr), stringsAsFactors = FALSE),
              score_v2(y, pred, m$n_eff_cont, scale = "level")),
        cbind(data.frame(country = cn, outcome = on, arm = arm, rep = r,
                         unit = "district", n_cols = ncol(Xr), stringsAsFactors = FALSE),
              score_v2(ad$y, ad$pred, ad$w, scale = "level")))
    }
  }
  message("  ", cn, " ", on, "  clusters ", nrow(m), "  districts ", length(ud))
}

raw <- dplyr::bind_rows(rows)
sm <- raw |>
  dplyr::group_by(country, outcome, arm, unit) |>
  dplyr::summarise(reps = dplyr::n(), n_cols = dplyr::first(n_cols),
                   spearman = mean(spearman, na.rm = TRUE),
                   pearson = mean(pearson, na.rm = TRUE), .groups = "drop") |>
  as.data.frame()
exp_write(sm, "11_lag_stack")

# ── descriptive: where in the lag profile does the association sit? ─────────
# Marginal correlation of the district-mean outcome with each (layer, lag).
prof <- list()
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  t <- TC[TC$country == cn & TC$outcome == on, ]
  m <- merge(t[is.finite(t$y_level), ], KEYDF[KEYDF$country == cn, ],
             by = c("country", "cluster"))
  if (nrow(m) < 25) next
  for (ly in LAYERS) for (lg in LAGS) {
    v <- m[[sprintf("%s_L%02d", ly, lg)]]
    if (is.null(v) || !is.finite(stats::sd(v, na.rm = TRUE)) || stats::sd(v, na.rm = TRUE) == 0) next
    prof[[length(prof) + 1L]] <- data.frame(
      country = cn, outcome = on, layer = ly, lag = lg,
      r = suppressWarnings(stats::cor(v, m$y_level, method = "spearman",
                                      use = "complete.obs")),
      stringsAsFactors = FALSE)
  }
}
PR <- dplyr::bind_rows(prof); exp_write(PR, "11_lag_profile")

cat("\n== median Spearman over cells, by arm and unit ==\n")
a <- aggregate(spearman ~ unit + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$unit, -a$spearman), ], row.names = FALSE)

cat("\n== paired gain over climatology (district scoring) ==\n")
d <- sm[sm$unit == "district", ]
w <- reshape(d[, c("country", "outcome", "arm", "spearman")],
             idvar = c("country", "outcome"), timevar = "arm", direction = "wide")
names(w) <- sub("^spearman\\.", "", names(w))
for (a in c("clim_lag_raw", "clim_lag_spline", "clim_lag_window", "clim_prev_seas")) {
  g <- w[[a]] - w$clim
  cat(sprintf("  %-16s mean gain %+.4f   improved %2d of %2d cells\n", a,
              mean(g, na.rm = TRUE), sum(g > 0, na.rm = TRUE), sum(is.finite(g))))
}

cat("\n== mean |marginal correlation| by lag, over all cells and layers ==\n")
mp <- aggregate(abs(r) ~ lag, data = PR, FUN = function(z) round(mean(z, na.rm = TRUE), 3))
names(mp)[2] <- "mean_abs_r"
print(mp, row.names = FALSE)
