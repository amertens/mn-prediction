# =============================================================================
# scripts/protocol_v2/56_weight_sources.R   [WS-01]
#
# WHERE SHOULD THE DOMAIN INDEX'S WEIGHTS COME FROM?
#
# Scores the production index against twelve alternative weightings of the SAME
# domain axes (R/protocol_v2_weights.R) under the three estimands of protocol
# v2, on the same cells, folds and targets as scripts 02 / 02b:
#   A  in-fill      5-fold by district, V2_REPS replicated draws (default 5)
#   B  region       exhaustive leave-one-region-out
#   C  country      leave-one-country-out on the pooled, within-country
#                   rank-normalised matrix; PCs oriented from the training
#                   countries; outcome standardised within country; scored by
#                   rank correlation inside the held-out country
# aux$rep_block carries the replication / inner-CV block: Admin-1 regions
# in-country, the country under C.
#
#   Rscript -e "source('scripts/protocol_v2/56_weight_sources.R')"   (V2_REPS, V2_ESTIMANDS, V2_ARMS, V2_CELLS, V2_OUT_TAG, V2_MERGE)
#   (source() parses the whole file first; `Rscript file.R` reads it incrementally and an edit during the run breaks it)
# -> results/tables/protocol_v2/weight_sources_raw.csv      every draw
#    results/tables/protocol_v2/weight_sources_cells.csv    cell medians over draws
#    results/tables/protocol_v2/weight_sources_summary.csv  estimand x target x arm
#    results/tables/protocol_v2/weight_sources_paired.csv   each arm minus the index, paired on cells
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R")

OUTDIR <- "results/tables/protocol_v2"
REPS   <- as.integer(Sys.getenv("V2_REPS", "5"))
ESTS   <- trimws(strsplit(Sys.getenv("V2_ESTIMANDS", "infill,region,country"), ",")[[1]])
ARMS   <- c("domain_index", "domain_enet", WS_ARMS, WS_SPARSE_ARMS)
if (nzchar(Sys.getenv("V2_ARMS", ""))) ARMS <- union("domain_index", trimws(strsplit(Sys.getenv("V2_ARMS"), ",")[[1]]))   # a partial pass, always paired with the index
MERGE_ONLY <- Sys.getenv("V2_MERGE", "0") == "1"   # skip fitting; rebuild the summaries from every weight_sources_raw*.csv on disk
CELLS_ONLY <- trimws(strsplit(Sys.getenv("V2_CELLS", ""), ",")[[1]])   # smoke: "Ghana:child_iron,Malawi:women_iron"
OUT_TAG <- Sys.getenv("V2_OUT_TAG", "")
set.seed(20260961L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

build_cell <- function(cn, on, target, with_D = TRUE) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y_nat <- m[[ycol]]
  list(country = cn, outcome = on, target = target, n = nrow(m), y_nat = y_nat,
       y_mod = if (target == "prev") .v2_logit(y_nat) else y_nat,
       X = Xr, D = if (with_D) domain_representation_v2(Xr, domain_of) else NULL,
       w = m[[wcol]], Admin1 = m$Admin1, lon = m$lon, lat = m$lat,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, rep_block = m$Admin1, y_nat = y_nat))
}

score_draw <- function(cell, estimand, folds, rep_id) {
  out <- list(); ws_memo_clear()
  for (a in ARMS) {
    fn <- ARMS_V2[[a]]; pred_mod <- rep(NA_real_, cell$n)
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12 || !length(te)) next
      p <- tryCatch(fn(tr, te, cell$y_mod, cell$X, cell$D, cell$aux), error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) pred_mod[te] <- p
    }
    pred <- if (cell$target == "prev") .v2_expit(pred_mod) else pred_mod
    s <- score_v2(cell$y_nat, pred, cell$w, scale = if (cell$target == "prev") "prev" else "level")
    out[[a]] <- cbind(data.frame(country = cell$country, outcome = cell$outcome, target = cell$target, estimand = estimand,
                                 arm = a, rep = rep_id, n_areas = cell$n), s)
  }
  bind_rows(out)
}

rows <- list()
if (MERGE_ONLY) ESTS <- character(0)
cells_index <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
if (length(CELLS_ONLY) && nzchar(CELLS_ONLY[1])) cells_index <- cells_index |> filter(paste(country, outcome, sep = ":") %in% CELLS_ONLY)
LOCO_OUTCOMES <- if (length(CELLS_ONLY) && nzchar(CELLS_ONLY[1])) unique(sub("^.*:", "", CELLS_ONLY)) else unique(TG$outcome)
if (any(c("infill", "region") %in% ESTS)) {
  for (i in seq_len(nrow(cells_index))) {
    cn <- cells_index$country[i]; on <- cells_index$outcome[i]
    for (target in c("prev", "level")) {
      cell <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
      if (is.null(cell)) next
      if ("infill" %in% ESTS) for (r in seq_len(REPS)) {
        folds <- make_folds_v2("kfold_district", cell$n, k = 5, rep_id = r)
        rows[[paste(cn, on, target, "infill", r)]] <- score_draw(cell, "infill", folds, r)
      }
      if ("region" %in% ESTS && dplyr::n_distinct(cell$Admin1) >= 3) {
        folds <- make_folds_v2("loro", cell$n, blocks = cell$Admin1)
        rows[[paste(cn, on, target, "region", 1)]] <- score_draw(cell, "region", folds, 1L)
      }
      cat("done", cn, on, target, format(Sys.time(), "%H:%M:%S"), "\n"); flush.console()
    }
  }
}

# checkpoint: the in-country draws are written before the cross-country step so a
# failure there cannot lose them (the 2026-09-09 WS-01 run did exactly that)
if (length(rows)) write.csv(bind_rows(rows), file.path(OUTDIR, paste0("weight_sources_raw", OUT_TAG, "_incountry_checkpoint.csv")), row.names = FALSE)
# ── estimand C ────────────────────────────────────────────────────────────────
if ("country" %in% ESTS) {
  for (target in c("prev", "level")) for (on in LOCO_OUTCOMES) {
    cl <- list()
    for (cn in COUNTRIES) { cc <- tryCatch(build_cell(cn, on, target, with_D = FALSE), error = function(e) NULL); if (!is.null(cc)) cl[[cn]] <- cc }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next
    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
    aux  <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
                 Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), rep_block = ctry, y_nat = Y)
    folds <- as.integer(factor(ctry))
    for (a in ARMS) {
      fn <- ARMS_V2[[a]]; pred <- rep(NA_real_, length(Y)); ws_memo_clear()
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 20) next
        Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)        # orientation from the training countries only
        p <- tryCatch(fn(tr, te, Y, Xm, Dm, aux), error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p
      }
      for (cn in unique(ctry)) {
        k <- which(ctry == cn)
        s <- score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")
        s$mae <- NA_real_; s$wmae <- NA_real_; s$bias <- NA_real_; s$rmse_sd <- NA_real_   # transport is a ranking claim
        rows[[paste(target, on, a, cn, "country")]] <- cbind(data.frame(country = cn, outcome = on, target = target, estimand = "country",
                                                                        arm = a, rep = 1L, n_areas = length(k)), s)
      }
    }
    cat("loco done", target, on, format(Sys.time(), "%H:%M:%S"), "\n"); flush.console()
  }
}

if (!MERGE_ONLY) { RAW <- bind_rows(rows); write.csv(RAW, file.path(OUTDIR, paste0("weight_sources_raw", OUT_TAG, ".csv")), row.names = FALSE) }
if (MERGE_ONLY) {
  fs <- list.files(OUTDIR, pattern = "^weight_sources_raw.*[.]csv$", full.names = TRUE); fs <- fs[!grepl("smoke|checkpoint", fs)]
  cat("merging", length(fs), "raw files:", paste(basename(fs), collapse = " "), "
")
  RAW <- bind_rows(lapply(fs, read.csv, stringsAsFactors = FALSE)) |>
    distinct(country, outcome, target, estimand, arm, rep, .keep_all = TRUE)   # domain_index is in every pass on identical folds
  ARMS <- c(intersect(ARMS, unique(RAW$arm)), setdiff(unique(RAW$arm), ARMS)); OUT_TAG <- ""   # only arms that were actually run, canonical order first
}
CELLS <- RAW |> group_by(country, outcome, target, estimand, arm) |>
  summarise(reps = n(), n_areas = max(n_areas), spearman = median(spearman, na.rm = TRUE), wmae = median(wmae, na.rm = TRUE),
            topk = median(topk, na.rm = TRUE), na_draws = sum(!is.finite(spearman)), .groups = "drop")
CELLS$spearman[!is.finite(CELLS$spearman)] <- NA_real_
write.csv(CELLS, file.path(OUTDIR, paste0("weight_sources_cells", OUT_TAG, ".csv")), row.names = FALSE)
SUM <- CELLS |> group_by(estimand, target, arm) |>
  summarise(cells = n(), scored = sum(is.finite(spearman)), mean_spearman = mean(spearman, na.rm = TRUE),
            median_spearman = median(spearman, na.rm = TRUE), cells_positive = sum(spearman > 0, na.rm = TRUE),
            mean_wmae = mean(wmae, na.rm = TRUE), .groups = "drop") |>
  mutate(arm = factor(arm, levels = ARMS)) |> arrange(estimand, target, arm) |> mutate(arm = as.character(arm))
write.csv(SUM, file.path(OUTDIR, paste0("weight_sources_summary", OUT_TAG, ".csv")), row.names = FALSE)
W <- CELLS |> select(country, outcome, target, estimand, arm, spearman) |> pivot_wider(names_from = arm, values_from = spearman)
PAIR <- bind_rows(lapply(setdiff(ARMS, "domain_index"), function(a) {
  W |> group_by(estimand, target) |>
    summarise(arm = a, paired = sum(is.finite(.data[[a]] - domain_index)),
              mean_diff = mean(.data[[a]] - domain_index, na.rm = TRUE), median_diff = median(.data[[a]] - domain_index, na.rm = TRUE),
              cells_better = sum(.data[[a]] > domain_index, na.rm = TRUE), .groups = "drop") }))
write.csv(PAIR, file.path(OUTDIR, paste0("weight_sources_paired", OUT_TAG, ".csv")), row.names = FALSE)

cat("\n== mean Spearman over cells (cells positive / scored) ==\n")
for (est in unique(SUM$estimand)) {
  d <- SUM |> filter(estimand == est) |> mutate(v = sprintf("%.3f (%d/%d)", mean_spearman, cells_positive, scored)) |>
    select(target, arm, v) |> pivot_wider(names_from = target, values_from = v)
  cat("--", est, "--\n"); print(as.data.frame(d), row.names = FALSE)
}
cat("\n== each arm minus the production index, paired on cells (mean diff; cells better / paired) ==\n")
d <- PAIR |> mutate(v = sprintf("%+.3f (%d/%d)", mean_diff, cells_better, paired)) |> select(estimand, target, arm, v) |>
  unite("k", estimand, target) |> pivot_wider(names_from = k, values_from = v)
print(as.data.frame(d), row.names = FALSE)
