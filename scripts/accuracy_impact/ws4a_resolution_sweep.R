# =============================================================================
# scripts/accuracy_impact/ws4a_resolution_sweep.R
#
# WS4a. At what spatial resolution does skill against a validated ceiling peak?
#
# THE QUESTION
# ------------
# Open question 2 of docs/SESSION_FINDINGS_FOR_REVIEW.md: is there a resolution
# BETWEEN Admin-1 and Admin-2 at which covariates remain informative and the
# target is adequately measured? WS5 bounds the answer without locating it:
# full-survey Admin-1 error is 0.297 pp and Admin-2 error is 5.09 to 13.66 pp.
# This script fills the interval.
#
# THE FOUR LEVELS
#   admin1          the region as the survey designed it
#   admin1_split2   each region cut into two spatially compact parts
#   admin1_split3   each region cut into three
#   admin2          the district
#
# THE SPLIT RULE, STATED BECAUSE IT IS A CHOICE
# ---------------------------------------------
# k-means on the district centroids WITHIN each region, k = 2 or 3, seeded.
# k-means on coordinates yields spatially compact parts, which is what makes the
# intermediate levels a plausible reporting geography rather than an arbitrary
# regrouping. It does not guarantee contiguity in the graph-theoretic sense: a
# region shaped around a bay can produce a part whose members are compact in
# Euclidean distance and not adjacent. The parts are therefore described as
# compact, not contiguous. A region with fewer than k districts yields fewer
# parts, and the realised count is recorded per cell.
#
# Centroids come from data/admin2_centroids.rds, which is keyed on country and
# district NAME with no Admin1. Malawi has six names that occur in more than one
# region, so centroids are averaged within (country, name) before use and the
# affected districts share a centroid. That is a known imprecision confined to
# Malawi and it is reported rather than hidden.
#
# WHAT IS MEASURED AT EACH LEVEL
#   reliability  the WS1a split-half estimator applied to the unit, so the
#                ceiling is re-measured at every resolution rather than assumed
#   skill        the same area-level ridge with an in-fold top-20 prescreen and
#                leave-one-region-out blocking, so the learner is held constant
#   n per unit   the median number of biomarker measurements per unit, which is
#                how the crossover is expressed in survey-design terms
#
# FOLDS. Blocking is always on Admin1, at every level, so a held-out region is
# held out whole. At admin1 that makes leave-one-region-out identical to
# leave-one-unit-out, which is thin: Sierra Leone has four regions, so a model
# trains on three units and is skipped. That is recorded, not worked around.
#
#   PROFILE=smoke   Ghana only
#   Rscript scripts/accuracy_impact/ws4a_resolution_sweep.R
# -> results/tables/resolution_sweep.csv
# -> results/figures/resolution_sweep.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
SEED <- 20260911L
B_SPLIT <- as.integer(Sys.getenv("WS4A_B", if (PROFILE == "smoke") "50" else "200"))
K_SCREEN <- 20L
MIN_UNITS_MODEL <- 8L      # units needed to fit and score at all
MIN_TRAIN_UNITS <- 5L      # training units needed inside a fold
TDIR <- here("results", "tables"); FDIR <- here("results", "figures")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
set.seed(SEED)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))

CENT <- readRDS(here("data", "admin2_centroids.rds")) |> as.data.frame()
CENT <- CENT |> group_by(country, Admin2) |>
  summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE),
            n_polygons = dplyr::n(), .groups = "drop") |> as.data.frame()

#' Assign each district to a unit at a given resolution.
#' @return character vector of unit labels, one per row of `key`
assign_units <- function(key, level, cent) {
  a1 <- key$Admin1; a2 <- key$Admin2
  if (level == "admin1") return(a1)
  if (level == "admin2") return(paste(a1, a2, sep = " | "))
  k <- if (level == "admin1_split2") 2L else 3L
  xy <- cent[match(a2, cent$Admin2), c("lon", "lat")]
  out <- a1
  for (r in unique(a1)) {
    i <- which(a1 == r)
    ok <- i[is.finite(xy$lon[i]) & is.finite(xy$lat[i])]
    kk <- min(k, length(unique(paste(xy$lon[ok], xy$lat[ok]))))
    if (length(ok) < 2L || kk < 2L) { out[i] <- paste(r, "p1"); next }
    m <- as.matrix(xy[ok, , drop = FALSE])
    m <- scale(m)
    m[!is.finite(m)] <- 0
    cl <- tryCatch(stats::kmeans(m, centers = kk, nstart = 25L)$cluster,
                   error = function(e) rep(1L, length(ok)))
    out[i] <- paste(r, "p1")                       # districts with no centroid
    out[ok] <- paste(r, paste0("p", cl))
  }
  out
}

LEVELS <- c("admin1", "admin1_split2", "admin1_split3", "admin2")
cfgs <- get_country_configs()
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]
rows <- list()

for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  hc <- H[H$country == cn, , drop = FALSE]
  ce <- CENT[CENT$country == cn, , drop = FALSE]
  if (!nrow(hc)) { cat("[skip]", cn, "no covariates\n"); next }
  n_shared_centroid <- sum(ce$n_polygons > 1)

  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od) || is.null(oc$binary)) next
    d <- od$data
    need <- c(cc$admin1_col, cc$admin2_col, cc$cluster_id, cc$weight_col, oc$binary)
    if (!all(need %in% names(d))) next
    y <- num(d[[oc$binary]])
    keep <- is.finite(y); d <- d[keep, , drop = FALSE]; y <- y[keep]
    if (length(unique(y)) < 2 || nrow(d) < 100) next
    d$.a1 <- trimws(as.character(d[[cc$admin1_col]]))
    d$.a2 <- trimws(as.character(d[[cc$admin2_col]]))
    d$.cl <- as.character(d[[cc$cluster_id]])
    d$.w  <- num(d[[cc$weight_col]]); d$.w[!is.finite(d$.w) | d$.w <= 0] <- 1
    d$.y  <- y
    ok <- !is.na(d$.a1) & !is.na(d$.a2)
    d <- d[ok, , drop = FALSE]
    if (nrow(d) < 100) next

    # District key: only districts the survey actually measured AND the
    # covariate table carries, so every level is built on the same districts.
    key <- unique(d[, c(".a1", ".a2")]); names(key) <- c("Admin1", "Admin2")
    key <- dplyr::inner_join(key, hc, by = admin2_join_by(key, hc))
    if (nrow(key) < MIN_UNITS_MODEL) next
    dk <- paste(d$.a1, d$.a2)
    d <- d[dk %in% paste(key$Admin1, key$Admin2), , drop = FALSE]
    if (nrow(d) < 100) next

    pop <- if ("ghs_pop" %in% names(key)) key$ghs_pop else rep(1, nrow(key))
    pop[!is.finite(pop) | pop <= 0] <- 1

    for (lv in LEVELS) {
      key$unit <- assign_units(key, lv, ce)
      d$.unit <- key$unit[match(paste(d$.a1, d$.a2), paste(key$Admin1, key$Admin2))]
      if (any(is.na(d$.unit))) d <- d[!is.na(d$.unit), , drop = FALSE]
      n_units <- dplyr::n_distinct(key$unit)

      # ---- the survey yardstick at this unit -------------------------------
      uw <- d |> group_by(.unit) |>
        summarise(prev = stats::weighted.mean(.y, .w), n = dplyr::n(),
                  .groups = "drop") |> as.data.frame()

      # ---- the ceiling, re-measured at this unit ---------------------------
      rel <- tryCatch(split_half_reliability(
        d, ".unit", ".cl", ".w", ".y", scheme = "within",
        B = B_SPLIT, seed = SEED, min_units = 4L), error = function(e) NULL)

      # ---- covariates aggregated to the unit, population weighted ----------
      agg <- do.call(rbind, lapply(split(seq_len(nrow(key)), key$unit), function(i)
        data.frame(unit = key$unit[i[1]], Admin1 = key$Admin1[i[1]],
                   as.list(vapply(COVS, function(v)
                     stats::weighted.mean(key[[v]][i], pop[i], na.rm = TRUE),
                     numeric(1))), stringsAsFactors = FALSE, check.names = FALSE)))
      m <- dplyr::inner_join(uw, agg, by = c(".unit" = "unit"))
      m <- m[is.finite(m$prev), , drop = FALSE]

      r <- NA_real_; mae <- NA_real_; n_scored <- 0L
      if (nrow(m) >= MIN_UNITS_MODEL && dplyr::n_distinct(m$Admin1) >= 3) {
        X <- as.matrix(m[, COVS, drop = FALSE]); X[!is.finite(X)] <- NA
        oof <- rep(NA_real_, nrow(m))
        for (rg in unique(m$Admin1)) {
          i <- which(m$Admin1 == rg); tr <- setdiff(seq_len(nrow(m)), i)
          if (length(tr) < MIN_TRAIN_UNITS) next
          fit <- .ds_fit(X[tr, , drop = FALSE], m$prev[tr], k_screen = K_SCREEN)
          pr <- .ds_predict(fit, X[i, , drop = FALSE])
          if (!is.null(pr)) oof[i] <- pr
        }
        fin <- is.finite(oof)
        if (sum(fin) >= MIN_UNITS_MODEL && stats::sd(m$prev[fin]) > 0) {
          r <- suppressWarnings(stats::cor(m$prev[fin], oof[fin]))
          mae <- 100 * mean(abs(oof[fin] - m$prev[fin]))
          n_scored <- sum(fin)
        }
      }

      rmx <- if (!is.null(rel)) rel$r_max_emp else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        country = cc$country, outcome = ocn, level = lv,
        n_units = n_units, n_scored = n_scored,
        median_n_per_unit = stats::median(uw$n),
        min_n_per_unit = min(uw$n),
        r_max_emp = if (is.finite(rmx)) round(rmx, 4) else NA_real_,
        r_oof = if (is.finite(r)) round(r, 4) else NA_real_,
        mae_pp = if (is.finite(mae)) round(mae, 2) else NA_real_,
        r_share = if (is.finite(rmx) && rmx > 0.05 && is.finite(r))
          round(r / rmx, 3) else NA_real_,
        n_shared_centroid_districts = n_shared_centroid,
        stringsAsFactors = FALSE)
      cat(sprintf("  %-13s %-13s %-14s units=%3d n/unit=%4.0f r_max=%5s r=%6s\n",
                  cn, ocn, lv, n_units, stats::median(uw$n),
                  ifelse(is.finite(rmx), sprintf("%.3f", rmx), "NA"),
                  ifelse(is.finite(r), sprintf("%+.3f", r), "NA")))
    }
  }
}

res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
res$level <- factor(res$level, levels = LEVELS)
readr::write_csv(res, file.path(TDIR, sprintf("resolution_sweep%s.csv", SUF)))

cat("\n=== WS4a: skill and ceiling by resolution ===\n")
s <- res |> group_by(level) |>
  summarise(cells = sum(is.finite(r_oof)),
            med_units = round(stats::median(n_units)),
            med_n_per_unit = round(stats::median(median_n_per_unit)),
            med_r_max = round(stats::median(r_max_emp, na.rm = TRUE), 3),
            med_r = round(stats::median(r_oof, na.rm = TRUE), 3),
            med_mae = round(stats::median(mae_pp, na.rm = TRUE), 2),
            med_r_share = round(stats::median(r_share, na.rm = TRUE), 3),
            .groups = "drop")
print(as.data.frame(s), row.names = FALSE)

cat("\n=== paired within cell: r_share by level, cells present at every level ===\n")
w <- res |> select(country, outcome, level, r_share) |>
  tidyr::pivot_wider(names_from = level, values_from = r_share)
cmp <- w[stats::complete.cases(w[, LEVELS]), , drop = FALSE]
if (nrow(cmp)) {
  print(as.data.frame(cmp), row.names = FALSE)
  best <- LEVELS[apply(cmp[, LEVELS], 1, which.max)]
  cat("\nlevel maximising r_share, per cell:\n"); print(table(best))
} else cat("no cell has r_share at all four levels\n")

cat("\n=== paired within cell: raw r by level ===\n")
wr <- res |> select(country, outcome, level, r_oof) |>
  tidyr::pivot_wider(names_from = level, values_from = r_oof)
cr <- wr[stats::complete.cases(wr[, LEVELS]), , drop = FALSE]
if (nrow(cr)) {
  print(as.data.frame(cr), row.names = FALSE)
  bestr <- LEVELS[apply(cr[, LEVELS], 1, which.max)]
  cat("\nlevel maximising raw r, per cell:\n"); print(table(bestr))
  cat(sprintf("\nmean r by level over the %d cells present at all four:\n", nrow(cr)))
  print(round(colMeans(cr[, LEVELS]), 3))
}

png(file.path(FDIR, sprintf("resolution_sweep%s.png", SUF)),
    width = 1200, height = 850, res = 130)
op <- par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
bx <- function(col, lab) {
  vals <- lapply(LEVELS, function(l) res[[col]][res$level == l])
  boxplot(vals, names = LEVELS, las = 2, ylab = lab, main = lab,
          col = "grey90", border = "grey30")
  pts <- vapply(vals, function(v) stats::median(v, na.rm = TRUE), numeric(1))
  lines(seq_along(LEVELS), pts, col = "firebrick", lwd = 2, type = "b", pch = 19)
}
bx("r_oof", "Out-of-fold r")
bx("r_share", "r / empirical ceiling")
par(op); dev.off()
cat(sprintf("\n-> %s\n", file.path("results", "figures",
                                   sprintf("resolution_sweep%s.png", SUF))))
