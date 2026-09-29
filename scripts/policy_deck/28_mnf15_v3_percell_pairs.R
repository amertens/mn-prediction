# =============================================================================
# scripts/policy_deck/28_mnf15_v3_percell_pairs.R
#
# "How often does the model put two districts in the survey's order?", per
# nutrient-country combination, with a 95% interval, for the v3 MNF15 talk
# (docs/slides/MNF15-talk-2026-09-v3.qmd). The pairwise metric is the share of
# district pairs whose order agrees with the survey (Kendall-type concordance;
# 50% = a coin toss), which is easier to explain than a Spearman correlation.
#
# Reproduces the committed benchmark before computing anything new:
#   in-country  = estimand A of scripts/protocol_v2/02_run_benchmarks_v2.R
#                 (headline tiers open + survey_public, domain components built
#                 once per cell, 5-fold x 10 draws, make_folds_v2 rep_id 1..10)
#   whole country held out = estimand C of 02b_merge_and_loco.R
# and stops unless every cell's mean Spearman matches benchmarks_v2_cells.csv
# to 0.005. Target: biomarker level (average status), as elsewhere in the talk.
#
# Interval: percentile bootstrap over districts (1,000 resamples). In-country,
# the statistic is computed on each district's held-out prediction averaged over
# the ten fold draws (point estimate and interval alike); the Spearman column is
# the benchmark's own per-draw mean, kept as the reproduction check.
#
#   Rscript scripts/policy_deck/28_mnf15_v3_percell_pairs.R
# -> results/tables/policy_deck/v3_percell_pairs.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R")
set.seed(20260927L)

P2 <- "results/tables/protocol_v2"; OUTT <- "results/tables/policy_deck"
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"))
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
TGT <- "level"; B <- 1000L

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2])
}))

# share of district pairs ordered as the survey orders them (ties in either are skipped)
concord <- function(o, p) {
  so <- sign(outer(o, o, "-")); sp <- sign(outer(p, p, "-"))
  m <- so * sp; ut <- upper.tri(m)
  agree <- sum(m[ut] > 0); tot <- sum(m[ut] != 0)
  if (tot == 0) NA_real_ else agree / tot
}

# as 02_run_benchmarks_v2.R build_cell (in-country) ------------------------------
cell_in <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(m = m, X = Xr, D = domain_representation_v2(Xr, domain_of), y = m$y_level,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = m$y_level, target = TGT))
}

rows <- list(); pts <- list()
for (cn in COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) {
  cc <- tryCatch(cell_in(cn, on), error = function(e) NULL); if (is.null(cc)) next
  n <- nrow(cc$m)
  P <- list(domain_index = matrix(NA_real_, n, 10), region_mean_jk = matrix(NA_real_, n, 10))
  sp <- list(domain_index = numeric(10), region_mean_jk = numeric(10))
  for (r in 1:10) {
    folds <- make_folds_v2("kfold_district", n, k = 5, rep_id = r)
    for (a in names(P)) {
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12) next
        P[[a]][te, r] <- ARMS_V2[[a]](tr, te, cc$y, cc$X, cc$D, cc$aux)
      }
      sp[[a]][r] <- score_v2(cc$y, P[[a]][, r], NULL, scale = "level")$spearman
    }
  }
  ref0 <- CELL$spearman[CELL$country == cn & CELL$outcome == on & CELL$target == TGT & CELL$estimand == "infill" & CELL$arm == "domain_index"]
  if (length(ref0) != 1 || !is.finite(ref0)) { cat("skipped (not scored in the benchmark):", cn, on, "\n"); next }
  for (a in names(P)) {
    ref <- CELL$spearman[CELL$country == cn & CELL$outcome == on & CELL$target == TGT & CELL$estimand == "infill" & CELL$arm == a]
    if (length(ref) != 1 || !is.finite(ref) || abs(mean(sp[[a]]) - ref) > 0.005)
      stop(sprintf("%s %s %s in-country does not reproduce: %.4f vs %s", cn, on, a, mean(sp[[a]]), paste(ref, collapse = ",")))
    pbar <- rowMeans(P[[a]])   # the prediction averaged over the ten fold draws, one per district
    stat <- function(idx) concord(cc$y[idx], pbar[idx])
    bt <- replicate(B, stat(sample.int(n, n, replace = TRUE)))
    rows[[length(rows) + 1]] <- data.frame(country = cn, outcome = on, estimand = "infill", arm = a, n = n,
      spearman = mean(sp[[a]]), pairs = stat(seq_len(n)), lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE))
  }
  cat(sprintf("in-country %-11s %-12s n=%2d  index %.3f (pairs %.0f%%)  regional %.3f\n", cn, on, n,
              mean(sp$domain_index), 100 * tail(rows, 2)[[1]]$pairs, mean(sp$region_mean_jk)))
}

# as 02b_merge_and_loco.R estimand C (whole country held out) --------------------
cell_lc <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), y = m$y_level, X = Xr, Admin1 = m$Admin1)
}
for (on in unique(TG$outcome)) {
  cl <- Filter(Negate(is.null), setNames(lapply(COUNTRIES, function(cn) tryCatch(cell_lc(cn, on), error = function(e) NULL)), COUNTRIES))
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y))))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y))
  aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  for (cn in names(cl)) {
    te <- which(ctry == cn); tr <- which(ctry != cn)
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    p <- ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, aux)
    s <- score_v2(ynat[te], p, NULL, scale = "level")$spearman
    ref <- CELL$spearman[CELL$country == cn & CELL$outcome == on & CELL$target == TGT & CELL$estimand == "country" & CELL$arm == "domain_index"]
    if (length(ref) != 1 || !is.finite(ref)) { cat("skipped held out (not scored):", cn, on, "\n"); next }
    if (abs(s - ref) > 0.005) stop(sprintf("%s %s held out does not reproduce: %.4f vs %s", cn, on, s, paste(ref, collapse = ",")))
    o <- ynat[te]; n <- length(te)
    bt <- replicate(B, { i <- sample.int(n, n, replace = TRUE); concord(o[i], p[i]) })
    rows[[length(rows) + 1]] <- data.frame(country = cn, outcome = on, estimand = "country", arm = "domain_index", n = n,
      spearman = s, pairs = concord(o, p), lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE))
  }
  cat("held out:", on, "\n")
}
R <- bind_rows(rows); rownames(R) <- NULL
write.csv(R, file.path(OUTT, "v3_percell_pairs.csv"), row.names = FALSE)
cat(sprintf("wrote v3_percell_pairs.csv: %d rows; all cells reproduce the committed benchmark\n", nrow(R)))
