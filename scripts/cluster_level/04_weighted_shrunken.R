# =============================================================================
# scripts/cluster_level/04_weighted_shrunken.R   [CL-04]
#
# TWO REPAIRS FOR THE NOISY CLUSTER OUTCOME
#
# CL-03 found that fitting at the cluster does not beat fitting at the
# district, and blamed the outcome: a cluster of 8-19 respondents is a noisy
# unit, and the zero-tuning index weights each domain by an UNWEIGHTED
# Spearman correlation with it. Two repairs, crossed:
#
#   weighted index   the domain weights are precision-weighted Spearman
#                    correlations (weights = the cluster's Kish n), so a
#                    ten-respondent cluster counts less than a thirty-one;
#                    everything else as arm_domain_index_v2
#   shrunken outcome empirical-Bayes shrinkage of each cluster's outcome on
#                    the modelling scale toward its district mean (or its
#                    region mean where the district is a single cluster),
#                    with the between-cluster variance estimated by the method
#                    of moments and the sampling variance from the cluster's
#                    own n: y_s = B y + (1 - B) m,  B = tau2 / (tau2 + v)
#
# Shrinkage uses only training-fold clusters (folds are cut by district and
# region, so a held-out district's clusters are never in the pool), and every
# score is still against the UNSHRUNK aggregated district outcome, exactly as
# in CL-03. Predictor sets: transportable and climate + soil.
#
#   Rscript scripts/cluster_level/04_weighted_shrunken.R      CL_REPS=10
# -> results/tables/cluster_level/benchmarks_cluster_ws_{cells,loco}.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- Sys.getenv("CL_OUT_DIR", "results/tables/cluster_level"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE); CDIR <- "data/covariates/cluster"
PT <- Sys.getenv("CL_PRED_TAG", "")   # which predictor table: "" = 2/5 km buffers, "_r10" = 10 km rural
REPS <- as.integer(Sys.getenv("CL_REPS", "10")); MIN_TRAIN <- 20L; set.seed(20260904L)
TC <- read.csv("results/tables/cluster_level/targets_cluster.csv", stringsAsFactors = FALSE)
P  <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
MD <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, "_metadata.csv")), stringsAsFactors = FALSE)
V2 <- tryCatch(read.csv("results/tables/protocol_v2/benchmarks_v2_cells.csv", stringsAsFactors = FALSE), error = function(e) NULL)
A1T <- tryCatch(read.csv("results/tables/protocol_v2/admin1_transport.csv", stringsAsFactors = FALSE), error = function(e) NULL)
TC$cluster <- as.character(TC$cluster); P$cluster <- as.character(P$cluster); MD <- MD[MD$column %in% names(P), ]
CS <- c("Climate and weather", "Soil characteristics")
SETS <- list(transportable = MD$column[MD$role != "fieldwork"], climate_soil = MD$column[MD$role != "fieldwork" & MD$domain %in% CS])
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
emp_logit <- function(y, n) log((y * n + 0.5) / ((1 - y) * n + 0.5))
group_folds <- function(groups, k, rep_id) { set.seed(20260951L + rep_id); g <- unique(groups); f <- sample(rep(seq_len(min(k, length(g))), length.out = length(g))); f[match(groups, g)] }
agg <- function(x, w, unit) { s <- tapply(x * w, unit, sum) / tapply(w, unit, sum); s[sort(names(s))] }

# ── the weighted index ───────────────────────────────────────────────────────
wspearman <- function(x, y, w) { rx <- rank(x); ry <- rank(y); w <- w / sum(w); mx <- sum(w * rx); my <- sum(w * ry)
  sx <- sqrt(sum(w * (rx - mx)^2)); sy <- sqrt(sum(w * (ry - my)^2)); if (sx == 0 || sy == 0) 0 else sum(w * (rx - mx) * (ry - my)) / (sx * sy) }
arm_domain_index_w <- function(tr, te, y, X, D, aux) {
  ytr <- y[tr]; w <- aux$w[tr]; w[!is.finite(w) | w <= 0] <- 1
  z <- apply(D[tr, , drop = FALSE], 2, function(x) { if (stats::sd(x) == 0) return(0); r <- wspearman(x, ytr, w); r <- max(min(r, 0.999), -0.999); 0.5 * log((1 + r) / (1 - r)) })
  z[!is.finite(z)] <- 0; idx_tr <- as.numeric(D[tr, , drop = FALSE] %*% z); idx_te <- as.numeric(D[te, , drop = FALSE] %*% z)
  if (stats::sd(idx_tr) == 0) return(rep(stats::weighted.mean(ytr, w), length(te)))
  ((idx_te - mean(idx_tr)) / stats::sd(idx_tr)) * stats::sd(ytr) + mean(ytr) }
ARMS <- list(index = ARMS_V2[["domain_index"]], index_w = arm_domain_index_w)

# ── the shrunken outcome ─────────────────────────────────────────────────────
# ymod, v: value and sampling variance on the modelling scale; pool: rows to
# learn tau2 and the group means from (the training fold); returns ymod_s for
# every row in `rows` (test rows are never shrunk because they are never used)
shrink <- function(ymod, v, district, region, pool) {
  y <- ymod[pool]; vv <- v[pool]; d <- district[pool]; r <- region[pool]
  md <- tapply(y, d, mean)[d]; nd <- tapply(y, d, length)[d]
  multi <- nd >= 2
  tau2 <- if (sum(multi) > 5) max(0, sum((y[multi] - md[multi])^2) / sum(1 - 1 / nd[multi]) - mean(vv[multi])) else NA_real_
  mr <- tapply(y, r, mean)[r]; nr <- tapply(y, r, length)[r]
  tau2r <- if (sum(nr >= 2) > 5) max(0, sum((y[nr >= 2] - mr[nr >= 2])^2) / sum(1 - 1 / nr[nr >= 2]) - mean(vv[nr >= 2])) else NA_real_
  ys <- y
  if (is.finite(tau2)) { B <- tau2 / (tau2 + vv); ys[multi] <- B[multi] * y[multi] + (1 - B[multi]) * md[multi] }
  if (is.finite(tau2r)) { B <- tau2r / (tau2r + vv); ys[!multi] <- B[!multi] * y[!multi] + (1 - B[!multi]) * mr[!multi] }
  out <- ymod; out[pool] <- ys; attr(out, "tau2") <- c(district = tau2, region = tau2r); out }

build <- function(cn, on, cols) {
  t <- TC[TC$country == cn & TC$outcome == on & is.finite(TC$lat), ]
  m <- inner_join(t, P[P$country == cn, setdiff(names(P), c("lat", "lon", "month_med"))], by = c("country", "cluster"))
  if (nrow(m) < 15) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, intersect(cols, names(m)), drop = FALSE])); if (ncol(Xr) < 10) return(NULL)
  list(country = cn, outcome = on, n = nrow(m), cluster = m$cluster, Admin1 = m$Admin1, Admin2 = paste(m$Admin1, m$Admin2, sep = "|"), lon = m$lon, lat = m$lat,
       y_prev = m$y_prev, y_level = m$y_level, sd_level = m$sd_level, n_raw = m$n_raw, w_prev = pmax(1, m$n_eff), w_level = pmax(1, m$n_eff_cont), X = Xr)
}
targets_of <- function(cl, target) {
  if (target == "prev") { ok <- is.finite(cl$y_prev) & is.finite(cl$w_prev) & cl$n_raw > 0; n <- cl$w_prev
    list(ok = ok, ymod = emp_logit(cl$y_prev, cl$n_raw), v = 1 / (cl$y_prev * n + 0.5) + 1 / ((1 - cl$y_prev) * n + 0.5), ynat = cl$y_prev, w = cl$w_prev, back = .v2_expit, scale = "prev") }
  else { ok <- is.finite(cl$y_level) & is.finite(cl$w_level); s2 <- cl$sd_level^2; s2[!is.finite(s2)] <- stats::median(s2[is.finite(s2)], na.rm = TRUE)
    list(ok = ok, ymod = cl$y_level, v = s2 / cl$w_level, ynat = cl$y_level, w = cl$w_level, back = identity, scale = "level") } }
score_both <- function(ynat, pred_nat, w, unit, scale) {
  s1 <- score_v2(ynat, pred_nat, w, scale = scale)
  ok <- is.finite(pred_nat); o <- agg(ynat[ok], w[ok], unit[ok]); p <- agg(pred_nat[ok], w[ok], unit[ok]); ww <- tapply(w[ok], unit[ok], sum)[names(o)]
  s2 <- score_v2(o, p, ww, scale = scale)
  data.frame(n_clusters = s1$n, rho_cluster = s1$spearman, n_units = s2$n, rho_agg = s2$spearman)
}
VARIANTS <- expand.grid(outcome_v = c("raw", "shrunk"), arm = names(ARMS), stringsAsFactors = FALSE)

# ── estimands A and B ────────────────────────────────────────────────────────
rows <- list(); cells <- TC |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells))) for (set_name in names(SETS)) {
  cl <- tryCatch(build(cells$country[i], cells$outcome[i], SETS[[set_name]]), error = function(e) NULL); if (is.null(cl)) next
  for (target in c("prev", "level")) { tv <- targets_of(cl, target); idx <- which(tv$ok); if (length(idx) < 15) next
    Y <- tv$ymod[idx]; V <- tv$v[idx]; X <- cl$X[idx, , drop = FALSE]; unit <- cl$Admin2[idx]; reg <- cl$Admin1[idx]; w <- tv$w[idx]
    aux <- list(lon = cl$lon[idx], lat = cl$lat[idx], Admin1 = reg, y_nat = Y, w = w)
    for (est in c("infill", "region")) { if (est == "region" && dplyr::n_distinct(reg) < 2) next
      for (r in seq_len(if (est == "infill") REPS else 1L)) {
        folds <- if (est == "infill") group_folds(unit, 5, r) else as.integer(factor(reg))
        preds <- lapply(seq_len(nrow(VARIANTS)), function(k) rep(NA_real_, length(Y))); tau <- c()
        for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
          D <- domain_representation_v2(X, domain_of, sign_rows = tr)
          Ys <- shrink(Y, V, unit, reg, tr); tau <- rbind(tau, attr(Ys, "tau2"))
          for (k in seq_len(nrow(VARIANTS))) { yy <- if (VARIANTS$outcome_v[k] == "shrunk") Ys else Y
            p <- tryCatch(ARMS[[VARIANTS$arm[k]]](tr, te, yy, X, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) preds[[k]][te] <- p } }
        for (k in seq_len(nrow(VARIANTS))) { s <- score_both(tv$ynat[idx], tv$back(preds[[k]]), w, unit, tv$scale)
          rows[[length(rows) + 1L]] <- cbind(data.frame(country = cl$country, outcome = cl$outcome, set = set_name, target = target, estimand = est, rep = r, arm = VARIANTS$arm[k], outcome_v = VARIANTS$outcome_v[k],
            tau2_district = if (length(tau)) mean(tau[, "district"], na.rm = TRUE) else NA_real_, mean_v = mean(V), stringsAsFactors = FALSE), s) }
      } }
  }
  cat("done", cells$country[i], cells$outcome[i], set_name, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "benchmarks_cluster_ws_raw.csv"), row.names = FALSE)
CELL <- R |> group_by(country, outcome, set, target, estimand, arm, outcome_v) |> summarise(rho_cluster = median(rho_cluster, na.rm = TRUE), rho_agg = median(rho_agg, na.rm = TRUE), shrink_B = median(1 - mean_v / (tau2_district + mean_v), na.rm = TRUE), .groups = "drop")
write.csv(CELL, file.path(OUTDIR, "benchmarks_cluster_ws_cells.csv"), row.names = FALSE)

# ── estimand C ───────────────────────────────────────────────────────────────
lrows <- list()
for (set_name in names(SETS)) for (on in unique(TC$outcome)) for (target in c("prev", "level")) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on, SETS[[set_name]]), error = function(e) NULL); if (is.null(z)) next
    tv <- targets_of(z, target); if (sum(tv$ok) >= 15) { z$tv <- tv; cl[[cn]] <- z } }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 10) next
  parts <- lapply(cl, function(z) { idx <- which(z$tv$ok); ys <- shrink(z$tv$ymod, z$tv$v, z$Admin2, z$Admin1, idx)   # within-country shrinkage, then standardise
    list(Y = as.numeric(scale(z$tv$ymod[idx])), Ys = as.numeric(scale(ys[idx])), X = z$X[idx, common, drop = FALSE], ynat = z$tv$ynat[idx], w = z$tv$w[idx],
         a2 = z$Admin2[idx], a1 = z$Admin1[idx], ctry = rep(z$country, length(idx)), lon = z$lon[idx], lat = z$lat[idx]) })
  g <- function(k) unlist(lapply(parts, `[[`, k))
  Y <- g("Y"); Ys <- g("Ys"); Xm <- do.call(rbind, lapply(parts, `[[`, "X")); ynat <- g("ynat"); wv <- g("w"); a2 <- g("a2"); a1 <- g("a1"); ctry <- g("ctry")
  aux <- list(lon = g("lon"), lat = g("lat"), Admin1 = paste(ctry, a1), y_nat = Y, w = wv); scale <- if (target == "prev") "prev" else "level"
  for (h in unique(ctry)) { te <- which(ctry == h); tr <- which(ctry != h); if (length(tr) < MIN_TRAIN) next
    D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    for (k in seq_len(nrow(VARIANTS))) { yy <- if (VARIANTS$outcome_v[k] == "shrunk") Ys else Y
      p <- tryCatch(ARMS[[VARIANTS$arm[k]]](tr, te, yy, Xm, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) != length(te)) next
      s <- score_both(ynat[te], p, wv[te], a2[te], scale)
      o1 <- agg(ynat[te], wv[te], a1[te]); p1 <- agg(p, wv[te], a1[te]); s1 <- score_v2(o1, p1, tapply(wv[te], a1[te], sum)[names(o1)], scale = scale)
      lrows[[length(lrows) + 1L]] <- data.frame(set = set_name, outcome = on, target = target, heldout = h, arm = VARIANTS$arm[k], outcome_v = VARIANTS$outcome_v[k],
        rho_cluster = s$rho_cluster, rho_admin2 = s$rho_agg, rho_admin1 = s1$spearman, stringsAsFactors = FALSE) } }
  cat("loco done", set_name, on, target, "\n")
}
LO <- bind_rows(lrows); write.csv(LO, file.path(OUTDIR, "benchmarks_cluster_ws_loco.csv"), row.names = FALSE)

# ── report ───────────────────────────────────────────────────────────────────
cat("\n===== CL-04: weighted index x shrunken outcome, cluster-fitted, aggregated to districts =====\n")
CELL$variant <- paste(CELL$outcome_v, CELL$arm, sep = "/")
cat("\nshrinkage factor B (median over cells; 1 = no shrinkage):\n"); print(as.data.frame(CELL |> filter(arm == "index", outcome_v == "shrunk", estimand == "infill", set == "transportable") |> group_by(target) |> summarise(B = round(median(shrink_B, na.rm = TRUE), 2), .groups = "drop")), row.names = FALSE)
v2c <- if (!is.null(V2)) V2 |> filter(estimand %in% c("infill", "region"), arm == "domain_index") |> group_by(country, outcome, target, estimand) |> summarise(rho_district_v2 = mean(spearman, na.rm = TRUE), .groups = "drop") else NULL
CMP <- if (!is.null(v2c)) left_join(CELL, v2c, by = c("country", "outcome", "target", "estimand")) else mutate(CELL, rho_district_v2 = NA_real_)
base <- CELL |> filter(variant == "raw/index") |> select(country, outcome, set, target, estimand, rho_base = rho_agg)
CMP <- left_join(CMP, base, by = c("country", "outcome", "set", "target", "estimand")) |>
  mutate(beats_base = is.finite(rho_agg) & is.finite(rho_base) & rho_agg > rho_base, beats_v2 = is.finite(rho_agg) & is.finite(rho_district_v2) & rho_agg > rho_district_v2, has_v2 = is.finite(rho_agg) & is.finite(rho_district_v2))
cat("\n-- in-fill and region: mean Spearman over cells (aggregated to districts); base = raw outcome, unweighted index (CL-03) --\n")
print(as.data.frame(CMP |> group_by(set, estimand, target, variant) |> summarise(cells = dplyr::n(), rho_agg_m = round(mean(rho_agg, na.rm = TRUE), 3), rho_cluster_m = round(mean(rho_cluster, na.rm = TRUE), 3),
  vs_base = sum(beats_base), district_v2 = round(mean(rho_district_v2, na.rm = TRUE), 3), beats_v2 = sum(beats_v2), paired_v2 = sum(has_v2), .groups = "drop") |> arrange(set, estimand, target, desc(rho_agg_m))), row.names = FALSE)
if (nrow(LO)) { LO$variant <- paste(LO$outcome_v, LO$arm, sep = "/")
  cat("\n-- transport (leave-one-country-out), cluster-fitted, scored at Admin-2 and Admin-1 --\n")
  print(as.data.frame(LO |> mutate(p2 = is.finite(rho_admin2) & rho_admin2 > 0, p1 = is.finite(rho_admin1) & rho_admin1 > 0) |> group_by(set, target, variant) |>
    summarise(cells = dplyr::n(), rho_admin2_m = round(mean(rho_admin2, na.rm = TRUE), 3), pos_admin2 = sum(p2), rho_admin1_m = round(mean(rho_admin1, na.rm = TRUE), 3), pos_admin1 = sum(p1), .groups = "drop") |> arrange(set, target, desc(rho_admin2_m))), row.names = FALSE)
  if (!is.null(V2)) { v <- V2 |> filter(estimand == "country", arm == "domain_index") |> select(country, outcome, target, rho_district_v2 = spearman)
    d <- LO |> inner_join(v, by = c("heldout" = "country", "outcome", "target"))
    print(as.data.frame(d |> mutate(b = is.finite(rho_admin2) & rho_admin2 > rho_district_v2) |> group_by(set, target, variant) |> summarise(vs_district_full_index = round(mean(rho_district_v2, na.rm = TRUE), 3), cluster_better = sum(b), n = dplyr::n(), .groups = "drop")), row.names = FALSE) }
  cat("reference, district-fitted climate+soil index (DA-01/DA-03): Admin-2 0.368 level / 0.268 prev; Admin-1 0.452 level / 0.378 prev\n") }
cat("\nDONE\n")
