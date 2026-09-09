# =============================================================================
# scripts/cluster_level/03_cluster_benchmarks.R   [CL-03]
#
# THE PROTOCOL-V2 BENCHMARK, FITTED AT THE CLUSTER
#
# Same three estimands, same arms, same conventions as scripts/protocol_v2/02,
# but the unit of fitting is the survey cluster (323 across four countries)
# with covariates extracted at the cluster (CL-02). Folds are still cut by
# DISTRICT (in-fill, 5-fold, replicated) and by REGION (exhaustive), never by
# cluster, so the estimands are unchanged; transport is leave-one-country-out.
# Every fit is scored twice:
#   rho_cluster   Spearman across clusters in the held-out fold
#   rho_agg       Spearman across DISTRICTS after aggregating cluster
#                 predictions to the district (n_eff-weighted), against the
#                 district's own aggregated outcome -- the number directly
#                 comparable with the district-level benchmark
# The prevalence target is modelled on the empirical logit (continuity
# correction 0.5, because clusters of ten respondents are often 0 or 1) and
# scored on the natural scale. Two predictor sets: transportable columns only
# (static, slow and dynamic climatology / survey-year layers), and the same
# plus the fieldwork-window adjusters, which is an in-country sensitivity.
#
#   Rscript scripts/cluster_level/03_cluster_benchmarks.R      CL_REPS=10
# -> results/tables/cluster_level/benchmarks_cluster_raw.csv
# -> results/tables/cluster_level/benchmarks_cluster_cells.csv
# -> results/tables/cluster_level/benchmarks_cluster_loco.csv
# -> results/tables/cluster_level/cluster_vs_district.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
source("R/protocol_v2_mbg.R")   # MB-01: registers the geostatistical arms mbg / mbg_gp into ARMS_V2
OUTDIR <- Sys.getenv("CL_OUT_DIR", "results/tables/cluster_level"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE); CDIR <- "data/covariates/cluster"
PT <- Sys.getenv("CL_PRED_TAG", "")   # which predictor table: "" = 2/5 km buffers, "_r10" = 10 km rural
REPS <- as.integer(Sys.getenv("CL_REPS", "10")); MIN_TRAIN <- 20L; set.seed(20260904L)
TC <- read.csv("results/tables/cluster_level/targets_cluster.csv", stringsAsFactors = FALSE)
P  <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
MD <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, "_metadata.csv")), stringsAsFactors = FALSE)
V2 <- tryCatch(read.csv("results/tables/protocol_v2/benchmarks_v2_cells.csv", stringsAsFactors = FALSE), error = function(e) NULL)
A1T <- tryCatch(read.csv("results/tables/protocol_v2/admin1_transport.csv", stringsAsFactors = FALSE), error = function(e) NULL)
TC$cluster <- as.character(TC$cluster); P$cluster <- as.character(P$cluster)
MD <- MD[MD$column %in% names(P), ]
CS <- c("Climate and weather", "Soil characteristics")
SETS <- list(transportable = MD$column[MD$role != "fieldwork"], with_fieldwork = MD$column,
             climate_soil = MD$column[MD$role != "fieldwork" & MD$domain %in% CS])
# CL_SETS / CL_ARMS restrict a run (comma-separated); CL_TAG suffixes its output
# files so a partial run can sit beside the main one; CL_SUMMARY_ONLY=1 skips
# the fitting and re-prints the comparison from whatever CSVs exist.
.env_list <- function(k, default) { v <- Sys.getenv(k, ""); if (nzchar(v)) trimws(strsplit(v, ",")[[1]]) else default }
SETS <- SETS[.env_list("CL_SETS", c("transportable", "with_fieldwork"))]
TAG <- Sys.getenv("CL_TAG", ""); SUMMARY_ONLY <- Sys.getenv("CL_SUMMARY_ONLY", "0") == "1"
domain_of <- stats::setNames(MD$domain, MD$column)
pref <- make.names(substr(unique(MD$domain), 1, 12)); if (any(duplicated(pref))) stop("domain prefix collision: ", paste(unique(MD$domain)[duplicated(pref)], collapse = " | "))
COUNTRIES <- .env_list("CL_COUNTRIES", c("Gambia", "Ghana", "Malawi", "SierraLeone"))   # CL_COUNTRIES shards the in-fill / region estimands (MB-01); transport needs all four
ARMS <- .env_list("CL_ARMS", c("null_train_mean", "region_mean_jk", "spatial", "domain_index", "domain_enet", "spatial_plus_domain"))
LARMS <- intersect(c("null_train_mean", "domain_index", "domain_enet"), ARMS)
emp_logit <- function(y, n) log((y * n + 0.5) / ((1 - y) * n + 0.5))
group_folds <- function(groups, k, rep_id) { set.seed(20260951L + rep_id); g <- unique(groups); f <- sample(rep(seq_len(min(k, length(g))), length.out = length(g))); f[match(groups, g)] }
agg <- function(x, w, unit) { s <- tapply(x * w, unit, sum) / tapply(w, unit, sum); s[sort(names(s))] }

build <- function(cn, on, cols) {
  t <- TC[TC$country == cn & TC$outcome == on & is.finite(TC$lat), ]
  m <- inner_join(t, P[P$country == cn, setdiff(names(P), c("lat", "lon", "month_med"))], by = c("country", "cluster"))
  if (nrow(m) < 15) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, intersect(cols, names(m)), drop = FALSE])); if (ncol(Xr) < 10) return(NULL)
  list(country = cn, outcome = on, n = nrow(m), cluster = m$cluster, Admin1 = m$Admin1, Admin2 = paste(m$Admin1, m$Admin2, sep = "|"), lon = m$lon, lat = m$lat,
       y_prev = m$y_prev, y_level = m$y_level, n_raw = m$n_raw, w_prev = pmax(1, m$n_eff), w_level = pmax(1, m$n_eff_cont), X = Xr)
}
targets_of <- function(cl, target) {
  if (target == "prev") list(ok = is.finite(cl$y_prev) & is.finite(cl$w_prev) & cl$n_raw > 0, ymod = emp_logit(cl$y_prev, cl$n_raw), ynat = cl$y_prev, w = cl$w_prev, back = .v2_expit, scale = "prev")
  else list(ok = is.finite(cl$y_level) & is.finite(cl$w_level), ymod = cl$y_level, ynat = cl$y_level, w = cl$w_level, back = identity, scale = "level") }
score_both <- function(ynat, pred_nat, w, unit, scale) {
  s1 <- score_v2(ynat, pred_nat, w, scale = scale)
  ok <- is.finite(pred_nat); o <- agg(ynat[ok], w[ok], unit[ok]); p <- agg(pred_nat[ok], w[ok], unit[ok]); ww <- tapply(w[ok], unit[ok], sum)[names(o)]
  s2 <- score_v2(o, p, ww, scale = scale)
  data.frame(n_clusters = s1$n, rho_cluster = s1$spearman, wmae_cluster = s1$wmae, n_units = s2$n, rho_agg = s2$spearman, wmae_agg = s2$wmae)
}

# ── estimands A (in-fill) and B (region), per country ─────────────────────────
rows <- list(); cells <- TC |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
if (SUMMARY_ONLY) cells <- cells[0, ]
for (i in seq_len(nrow(cells))) for (set_name in names(SETS)) {
  cl <- tryCatch(build(cells$country[i], cells$outcome[i], SETS[[set_name]]), error = function(e) { cat("  build error", conditionMessage(e), "\n"); NULL }); if (is.null(cl)) next
  for (target in c("prev", "level")) { tv <- targets_of(cl, target); idx <- which(tv$ok); if (length(idx) < 15) next
    Y <- tv$ymod[idx]; X <- cl$X[idx, , drop = FALSE]; unit <- cl$Admin2[idx]; reg <- cl$Admin1[idx]
    aux <- list(lon = cl$lon[idx], lat = cl$lat[idx], Admin1 = reg, y_nat = Y,
                target = target, n_raw = cl$n_raw[idx], y_prev = cl$y_prev[idx])   # MB-01: counts for the geostatistical arm's binomial likelihood
    for (est in c("infill", "region")) { if (est == "region" && (dplyr::n_distinct(reg) < 2 || set_name == "with_fieldwork")) next
      set_arms <- if (set_name == "climate_soil") intersect(ARMS, c("null_train_mean", "region_mean_jk", "domain_index")) else ARMS
      for (r in seq_len(if (est == "infill") REPS else 1L)) {
        folds <- if (est == "infill") group_folds(unit, 5, r) else as.integer(factor(reg))
        arms <- if (est == "infill") set_arms else setdiff(set_arms, "region_mean_jk")
        preds <- lapply(arms, function(a) rep(NA_real_, length(Y))); names(preds) <- arms
        for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
          D <- domain_representation_v2(X, domain_of, sign_rows = tr)
          for (a in arms) { p <- tryCatch(ARMS_V2[[a]](tr, te, Y, X, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) preds[[a]][te] <- p } }
        for (a in arms) { s <- score_both(tv$ynat[idx], tv$back(preds[[a]]), tv$w[idx], unit, tv$scale)
          rows[[length(rows) + 1L]] <- cbind(data.frame(country = cl$country, outcome = cl$outcome, set = set_name, target = target, estimand = est, rep = r, arm = a, n_pred = ncol(X), stringsAsFactors = FALSE), s) }
      } }
  }
  cat("done", cells$country[i], cells$outcome[i], set_name, "\n")
}
if (!SUMMARY_ONLY) { R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, paste0("benchmarks_cluster_raw", TAG, ".csv")), row.names = FALSE)
  CELL <- R |> group_by(country, outcome, set, target, estimand, arm) |> summarise(reps = dplyr::n(), n_clusters = max(n_clusters), n_units = max(n_units),
    rho_cluster = median(rho_cluster, na.rm = TRUE), rho_agg = median(rho_agg, na.rm = TRUE), wmae_cluster = median(wmae_cluster, na.rm = TRUE), wmae_agg = median(wmae_agg, na.rm = TRUE), .groups = "drop")
  write.csv(CELL, file.path(OUTDIR, paste0("benchmarks_cluster_cells", TAG, ".csv")), row.names = FALSE) }

# ── estimand C (transport), pooled ───────────────────────────────────────────
lrows <- list(); LSETS <- intersect(names(SETS), c("transportable", "climate_soil")); if (SUMMARY_ONLY || length(COUNTRIES) < 4) LSETS <- character(0)
for (lset in LSETS) for (on in unique(TC$outcome)) for (target in c("prev", "level")) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on, SETS[[lset]]), error = function(e) NULL); if (is.null(z)) next
    tv <- targets_of(z, target); if (sum(tv$ok) >= 15) { z$tv <- tv; cl[[cn]] <- z } }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 10) next
  parts <- lapply(cl, function(z) { idx <- which(z$tv$ok); list(Y = as.numeric(scale(z$tv$ymod[idx])), X = z$X[idx, common, drop = FALSE], ynat = z$tv$ynat[idx], w = z$tv$w[idx],
    a2 = z$Admin2[idx], a1 = z$Admin1[idx], ctry = rep(z$country, length(idx)), lon = z$lon[idx], lat = z$lat[idx]) })
  g <- function(k) unlist(lapply(parts, `[[`, k))
  Y <- g("Y"); Xm <- do.call(rbind, lapply(parts, `[[`, "X")); ynat <- g("ynat"); wv <- g("w"); a2 <- g("a2"); a1 <- g("a1"); ctry <- g("ctry")
  aux <- list(lon = g("lon"), lat = g("lat"), Admin1 = paste(ctry, a1), y_nat = Y); scale <- if (target == "prev") "prev" else "level"
  for (h in unique(ctry)) { te <- which(ctry == h); tr <- which(ctry != h); if (length(tr) < MIN_TRAIN) next
    D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    for (a in if (lset == "climate_soil") intersect(LARMS, c("null_train_mean", "domain_index")) else LARMS) {
      p <- tryCatch(ARMS_V2[[a]](tr, te, Y, Xm, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) != length(te)) next
      s <- score_both(ynat[te], p, wv[te], a2[te], scale)
      o1 <- agg(ynat[te], wv[te], a1[te]); p1 <- agg(p, wv[te], a1[te]); s1 <- score_v2(o1, p1, tapply(wv[te], a1[te], sum)[names(o1)], scale = scale)
      lrows[[length(lrows) + 1L]] <- data.frame(set = lset, outcome = on, target = target, heldout = h, arm = a, n_clusters = s$n_clusters, rho_cluster = s$rho_cluster,
        n_admin2 = s$n_units, rho_admin2 = s$rho_agg, n_admin1 = s1$n, rho_admin1 = s1$spearman, stringsAsFactors = FALSE) } }
  cat("loco done", lset, on, target, "\n")
}
if (!SUMMARY_ONLY) { LO <- bind_rows(lrows); write.csv(LO, file.path(OUTDIR, paste0("benchmarks_cluster_loco", TAG, ".csv")), row.names = FALSE) }
# gather every run's tables (main + tagged partial runs) for the report
rd_all <- function(stem) { fs <- list.files(OUTDIR, pattern = paste0("^", stem, "(_[A-Za-z0-9]+)?[.]csv$"), full.names = TRUE); fs <- fs[file.size(fs) > 20]   # a shard without transport writes an empty LOCO file
  if (!length(fs)) return(NULL); bind_rows(lapply(fs, function(f) { d <- read.csv(f, stringsAsFactors = FALSE); if (!"set" %in% names(d)) d$set <- "transportable"; d })) |> distinct() }
CELL <- rd_all("benchmarks_cluster_cells"); LO <- rd_all("benchmarks_cluster_loco")

# ── side by side with the district-level benchmark ───────────────────────────
cat("\n===== CL-03: cluster-level benchmark =====\n")
# (counts are computed from the cell-level values BEFORE any mean is taken:
#  the first version redefined the column inside summarise and counted a scalar)
cat("\n-- in-fill and region: mean Spearman over cells, cluster-fitted (at the cluster / aggregated to districts) vs district-fitted (protocol v2) --\n")
v2c <- if (!is.null(V2)) V2 |> filter(estimand %in% c("infill", "region")) |> group_by(country, outcome, target, estimand, arm) |> summarise(rho_district_v2 = mean(spearman, na.rm = TRUE), .groups = "drop") else NULL
CMP <- if (!is.null(v2c)) CELL |> left_join(v2c, by = c("country", "outcome", "target", "estimand", "arm")) else CELL |> mutate(rho_district_v2 = NA_real_)
write.csv(CMP, file.path(OUTDIR, "cluster_vs_district.csv"), row.names = FALSE)
S1 <- CMP |> mutate(beats = is.finite(rho_agg) & is.finite(rho_district_v2) & rho_agg > rho_district_v2, paired = is.finite(rho_agg) & is.finite(rho_district_v2)) |>
  group_by(set, target, estimand, arm) |> summarise(cells = dplyr::n(), rho_cluster_m = round(mean(rho_cluster, na.rm = TRUE), 3), rho_agg_m = round(mean(rho_agg, na.rm = TRUE), 3),
    rho_district_v2_m = round(mean(rho_district_v2, na.rm = TRUE), 3), agg_beats_v2 = sum(beats), paired = sum(paired), .groups = "drop") |> arrange(estimand, target, set, desc(rho_agg_m))
print(as.data.frame(S1), row.names = FALSE)
cat("\n-- fieldwork adjusters (in-fill sensitivity), aggregated-to-district Spearman: with vs without --\n")
W <- CELL |> filter(estimand == "infill", set %in% c("transportable", "with_fieldwork"), arm %in% c("domain_index", "spatial_plus_domain")) |> select(country, outcome, target, arm, set, rho_agg) |> pivot_wider(names_from = set, values_from = rho_agg, values_fn = mean)   # tagged partial runs can duplicate a key (MB-01)
if (all(c("transportable", "with_fieldwork") %in% names(W))) print(as.data.frame(W |> mutate(better = is.finite(with_fieldwork) & is.finite(transportable) & with_fieldwork > transportable) |> group_by(target, arm) |>
  summarise(cells = dplyr::n(), transportable_m = round(mean(transportable, na.rm = TRUE), 3), with_fieldwork_m = round(mean(with_fieldwork, na.rm = TRUE), 3), better = sum(better), .groups = "drop")), row.names = FALSE)
if (!is.null(LO) && nrow(LO)) { cat("\n-- transport (leave-one-country-out), cluster-fitted, scored at the cluster / Admin-2 / Admin-1 --\n")
  print(as.data.frame(LO |> mutate(p2 = is.finite(rho_admin2) & rho_admin2 > 0, p1 = is.finite(rho_admin1) & rho_admin1 > 0) |> group_by(set, target, arm) |>
    summarise(cells = dplyr::n(), rho_cluster_m = round(mean(rho_cluster, na.rm = TRUE), 3), rho_admin2_m = round(mean(rho_admin2, na.rm = TRUE), 3), pos_admin2 = sum(p2),
              rho_admin1_m = round(mean(rho_admin1, na.rm = TRUE), 3), pos_admin1 = sum(p1), .groups = "drop") |> arrange(target, set, desc(rho_admin2_m))), row.names = FALSE)
  if (!is.null(V2)) { v <- V2 |> filter(estimand == "country", arm == "domain_index") |> select(country, outcome, target, rho_district_v2 = spearman)
    for (ls in unique(LO$set)) { d <- LO |> filter(set == ls, arm == "domain_index") |> inner_join(v, by = c("heldout" = "country", "outcome", "target"))
      for (tg in unique(d$target)) { dd <- d[d$target == tg, ]
        cat(sprintf("  [%s, %s] index transport at Admin-2: cluster-fitted %.3f vs district-fitted (full vocabulary) %.3f | cluster better in %d of %d\n", ls, tg, mean(dd$rho_admin2, na.rm = TRUE), mean(dd$rho_district_v2, na.rm = TRUE), sum(dd$rho_admin2 > dd$rho_district_v2, na.rm = TRUE), sum(is.finite(dd$rho_admin2 - dd$rho_district_v2)))) } } }
  if (!is.null(A1T) && all(c("country", "outcome", "arm", "spearman") %in% names(A1T))) { v1 <- A1T |> filter(arm == "domain_index") |> group_by(country, outcome) |> summarise(rho_a1_v2 = mean(spearman, na.rm = TRUE), .groups = "drop")
    for (ls in unique(LO$set)) { d <- LO |> filter(set == ls, arm == "domain_index", target == "prev") |> inner_join(v1, by = c("heldout" = "country", "outcome"))
      if (nrow(d)) cat(sprintf("  [%s, prev] index transport at Admin-1: cluster-fitted %.3f vs district-fitted %.3f | cluster better in %d of %d\n", ls, mean(d$rho_admin1, na.rm = TRUE), mean(d$rho_a1_v2, na.rm = TRUE), sum(d$rho_admin1 > d$rho_a1_v2, na.rm = TRUE), sum(is.finite(d$rho_admin1 - d$rho_a1_v2)))) } } }
cat("\nDONE\n")
