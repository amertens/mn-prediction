# =============================================================================
# scripts/protocol_v2/23_domain_ablation_loco.R   [DA-01]
#
# WHICH DOMAINS CARRY TRANSPORT?
#
# Transport of district rankings to an unseen country (leave-one-country-out)
# is the project's headline capability. It is built on 20 domains of
# indicators reduced to principal components. Nobody has asked which domains
# the transported ranking actually depends on. Two ablations, on the identical
# LOCO folds, with the zero-tuning domain index (fast, and the arm the NCE
# cites):
#
#   drop-one   remove one domain's PCs; transport delta = full - without
#   only-one   keep one domain's PCs; how far does that domain get alone?
#
# A domain whose removal hurts and which alone transports is load-bearing; a
# domain that hurts when removed but is useless alone is complementary; a
# domain that neither hurts nor helps is dead weight for transport and is a
# candidate to stop collecting. This feeds the RA annotation priorities and
# the "what to collect" slide.
#
# PROTOCOL as in 02b_merge_and_loco.R: outcome z-scored within country,
# predictors rank-normalised within country and pooled on common columns,
# domain PCs oriented from the training countries only.
#
#   Rscript scripts/protocol_v2/23_domain_ablation_loco.R
# -> results/tables/protocol_v2/domain_ablation_loco.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2],
             stringsAsFactors = FALSE) }))

build_cell <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; if (!all(c(ycol, wcol) %in% names(t))) return(NULL)
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]; if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  yn <- m[[ycol]]
  list(country = cn, n = nrow(m), y_nat = yn, y_mod = if (target == "prev") .v2_logit(yn) else yn,
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}

# map a domain-PC column back to its domain: columns are named
# make.names(substr(domain, 1, 12))_PC<k> by build_domain_pcs_v2()
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) { pre <- sub("_PC[0-9]+$", "", cols); unname(prefix_of[pre]) }

loco_rho <- function(Dm_sel, Y, ctry, ynat, wv, aux, target) {
  folds <- as.integer(factor(ctry)); pred <- rep(NA_real_, length(Y))
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 20) next
    p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm_sel, aux), error = function(e) rep(NA_real_, length(te)))
    if (length(p) == length(te)) pred[te] <- p
  }
  vapply(unique(ctry), function(cn) { k <- which(ctry == cn)
    score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")$spearman }, numeric(1))
}

rows <- list()
for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
              Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  # orientation must be learned from training countries per fold; ablation
  # subsets columns of the SAME per-fold D, so build D per fold inside loco_rho
  # via a wrapper that recomputes it. To keep cost manageable, D is built once
  # per held-out country and the column subset applied afterwards.
  folds <- as.integer(factor(ctry))
  Dlist <- lapply(unique(folds), function(f) domain_representation_v2(Xm, domain_of, sign_rows = which(folds != f)))
  names(Dlist) <- unique(folds)
  # The number of PCs per domain is chosen to 80% variance ON THE TRAINING
  # ROWS, so the column set of D differs by held-out country. Subsetting by a
  # name list taken from one fold therefore fails on another ("subscript out
  # of bounds" in the first run). Intersect per fold, and never index by NA.
  run <- function(keep_cols) {
    pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 20) next
      Df <- Dlist[[as.character(f)]]
      kc <- intersect(keep_cols[!is.na(keep_cols)], colnames(Df))
      if (length(kc) < 1) next
      Dm <- Df[, kc, drop = FALSE]
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm, aux), error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) pred[te] <- p
    }
    vapply(unique(ctry), function(cn) { k <- which(ctry == cn)
      score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")$spearman }, numeric(1))
  }
  allcols <- unique(unlist(lapply(Dlist, colnames))); dom_of_col <- col_domain(allcols)
  if (any(is.na(dom_of_col)))
    cat("  unmapped D columns (excluded from ablation):", paste(allcols[is.na(dom_of_col)], collapse = ", "), "
")
  full <- run(allcols)
  for (cn in names(full)) rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn,
    domain = "ALL", variant = "full", spearman = full[[cn]], stringsAsFactors = FALSE)
  for (dm in unique(stats::na.omit(dom_of_col))) {
    drop <- run(allcols[which(is.na(dom_of_col) | dom_of_col != dm)])
    only <- run(allcols[which(!is.na(dom_of_col) & dom_of_col == dm)])
    for (cn in names(full)) {
      rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn, domain = dm,
        variant = "drop", spearman = drop[[cn]], stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn, domain = dm,
        variant = "only", spearman = only[[cn]], stringsAsFactors = FALSE)
    }
  }
  cat("ablation done", target, on, "\n")
}
R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "domain_ablation_loco.csv"), row.names = FALSE)

full <- R[R$variant == "full", ] |> group_by(target) |> summarise(full = mean(spearman, na.rm = TRUE), .groups = "drop")
A <- R[R$variant != "full", ] |> group_by(target, domain, variant) |>
  summarise(rho = mean(spearman, na.rm = TRUE), cells = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = variant, values_from = rho) |> left_join(full, by = "target") |>
  mutate(delta_drop = round(full - drop, 3), only = round(only, 3), full = round(full, 3), drop = round(drop, 3)) |>
  arrange(target, desc(delta_drop))
write.csv(A, file.path(OUTDIR, "domain_ablation_loco_summary.csv"), row.names = FALSE)
for (tg in unique(A$target)) {
  cat("\n===== DA-01 · transport (LOCO Spearman, mean over cells) · target =", tg, "=====\n")
  cat("full model:", A$full[A$target == tg][1], "\n")
  cat("delta_drop = loss when the domain is REMOVED (positive = load-bearing); only = the domain ALONE\n")
  print(as.data.frame(A[A$target == tg, c("domain", "cells", "drop", "delta_drop", "only")]), row.names = FALSE)
}
cat("\nDONE\n")
