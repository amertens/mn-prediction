# =============================================================================
# scripts/protocol_v2/43_source_ablation_loco.R   [DA-04]
#
# WHICH DATA SOURCES CARRY TRANSPORT? (the old deck's slide 17, redone)
#
# The January 2026 deck asked "which data sources contribute most to model
# performance" and answered with a change-in-MSE importance from the
# prescreened SuperLearner, per source (DHS, GEE, WFP, MAP, IHME, LSMS,
# FluNet). The corrected protocol replaced that with a domain ablation
# (DA-01, script 23). Sources and domains are not the same partition: the
# climate domain mixes GEE with Koppen/AEZ, the agriculture domain mixes
# MapSPAM with AlphaEarth, and DHS spans nine domains. This script is the
# source-level counterpart of DA-01: drop every column from one source (or
# keep only that source), rebuild the domain PCs on what is left, and score
# the zero-tuning index under leave-one-country-out. Same cells, same folds,
# same estimator as DA-01, so the two tables are directly comparable.
#
#   Rscript scripts/protocol_v2/43_source_ablation_loco.R
# -> results/tables/protocol_v2/source_ablation_loco.csv, _summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; MIN_TRAIN <- 20L
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
# collapse the long source labels to the names a slide can carry
src_short <- function(s) { s <- as.character(s)
  dplyr::case_when(grepl("^DHS", s) ~ "DHS", grepl("AlphaEarth", s) ~ "AlphaEarth", grepl("SoilGrids", s) ~ "SoilGrids/iSDA", grepl("^GEE", s) ~ "GEE",
    grepl("IHME", s) ~ "IHME", grepl("Malaria Atlas", s) ~ "Malaria Atlas", grepl("MapSPAM", s) ~ "MapSPAM", grepl("Koppen", s) ~ "Koppen/AEZ",
    grepl("WFP", s) ~ "WFP prices", grepl("FAOSTAT", s) ~ "FAOSTAT", grepl("WorldPop|GHS", s) ~ "WorldPop/GHS", grepl("RWI", s) ~ "RWI/GPW",
    grepl("HFID", s) ~ "HFID", grepl("GFDx", s) ~ "GFDx", TRUE ~ s) }
source_of <- stats::setNames(src_short(MD$source), MD$column)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
build_cell <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]; if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")); if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  yn <- m[[ycol]]; list(country = cn, n = nrow(m), y_nat = yn, y_mod = if (target == "prev") .v2_logit(yn) else yn, X = Xr, w = m[[wcol]], Admin1 = m$Admin1)
}
rows <- list()
for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y); folds <- as.integer(factor(ctry))
  run <- function(cols) { if (length(cols) < 2) return(stats::setNames(rep(NA_real_, length(unique(ctry))), unique(ctry)))
    pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < MIN_TRAIN) next
      Dm <- domain_representation_v2(Xm[, cols, drop = FALSE], domain_of, sign_rows = tr); if (!ncol(Dm)) next
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) pred[te] <- p }
    vapply(unique(ctry), function(cn) { k <- which(ctry == cn); score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")$spearman }, numeric(1)) }
  src <- source_of[common]; full <- run(common)
  for (cn in names(full)) rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn, source = "ALL", variant = "full", n_cols = length(common), spearman = full[[cn]], stringsAsFactors = FALSE)
  for (s in unique(src)) { drop <- run(common[src != s]); only <- run(common[src == s])
    for (cn in names(full)) { rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn, source = s, variant = "drop", n_cols = sum(src != s), spearman = drop[[cn]], stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = cn, source = s, variant = "only", n_cols = sum(src == s), spearman = only[[cn]], stringsAsFactors = FALSE) } }
  cat("source ablation done", target, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "source_ablation_loco.csv"), row.names = FALSE)
SUMM <- R |> group_by(target, source, variant) |> summarise(cells = dplyr::n(), n_cols = max(n_cols), spearman = mean(spearman, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = variant, values_from = c(spearman, n_cols)) |> filter(source != "ALL")
FULL <- R |> filter(variant == "full") |> group_by(target) |> summarise(full = mean(spearman, na.rm = TRUE), .groups = "drop")
SUMM <- left_join(SUMM, FULL, by = "target") |> mutate(delta_drop = full - spearman_drop, n_cols = n_cols_only) |>
  select(target, source, cells, n_cols, full, drop = spearman_drop, only = spearman_only, delta_drop) |> arrange(target, desc(delta_drop))
write.csv(SUMM, file.path(OUTDIR, "source_ablation_loco_summary.csv"), row.names = FALSE)
cat("\n===== DA-04: transport (LOCO Spearman, zero-tuning index) by DATA SOURCE =====\n")
for (tg in unique(SUMM$target)) { cat(sprintf("\n-- target = %s | full model %.3f --\n", tg, FULL$full[FULL$target == tg]))
  print(as.data.frame(SUMM[SUMM$target == tg, ] |> mutate(across(c(full, drop, only, delta_drop), ~ round(.x, 3)))), row.names = FALSE) }
cat("\nDONE\n")
