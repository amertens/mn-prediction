# =============================================================================
# scripts/protocol_v2/02b_merge_and_loco.R
#
# Merge the per-country shards written by 02_run_benchmarks_v2.R with
# V2_COUNTRY set, then run estimand C (transport), which needs every country
# at once and so cannot be sharded. Produces the same three output files the
# unsharded run would have produced.
#
#   Rscript scripts/protocol_v2/02b_merge_and_loco.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")],
                  by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y_nat <- m[[ycol]]
  # NOTE: domain scores are deliberately NOT built here. Building them per
  # country lets each country learn its own PC1 sign orientation, so a domain
  # score means the opposite thing in two countries and transport is destroyed.
  # For estimand C they are built once on the POOLED rank-normalised matrix,
  # with the orientation learned from the TRAINING countries only.
  list(country = cn, n = nrow(m), y_nat = y_nat,
       y_mod = if (target == "prev") .v2_logit(y_nat) else y_nat,
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1,
       lon = m$lon, lat = m$lat)
}

shards <- list.files(OUTDIR, pattern = "^benchmarks_v2_raw_.*\\.csv$",
                     full.names = TRUE)
if (!length(shards)) stop("no shards found - run 02 with V2_COUNTRY set")
cat("merging", length(shards), "shards\n")
RAW <- bind_rows(lapply(shards, read.csv, stringsAsFactors = FALSE))

# ── estimand C ──────────────────────────────────────────────────────────────
loco_rows <- list()
for (target in c("prev", "level")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      cc <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
      if (!is.null(cc)) cl[[cn]] <- cc
    }
    if (length(cl) < 3) next
    # pool the within-country rank-normalised matrices on their common columns
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next
    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat))
    wv   <- unlist(lapply(cl, function(z) z$w))
    aux  <- list(lon = unlist(lapply(cl, function(z) z$lon)),
                 lat = unlist(lapply(cl, function(z) z$lat)),
                 Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))),
                 y_nat = Y)
    folds <- as.integer(factor(ctry))
    for (a in c("null_train_mean", "domain_index", "domain_enet")) {
      fn <- ARMS_V2[[a]]
      pred <- rep(NA_real_, length(Y))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 20) next
        # orientation learned from TRAINING countries only, then applied to all
        Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        p <- tryCatch(fn(tr, te, Y, Xm, Dm, aux),
                      error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p
      }
      for (cn in unique(ctry)) {
        k <- which(ctry == cn)
        s <- score_v2(ynat[k], pred[k], wv[k],
                      scale = if (target == "prev") "prev" else "level")
        # transport is a RANKING claim: level metrics are not claimed
        s$mae <- NA_real_; s$wmae <- NA_real_; s$bias <- NA_real_
        s$rmse_sd <- NA_real_
        loco_rows[[paste(target, on, a, cn)]] <- cbind(
          data.frame(country = cn, outcome = on, target = target,
                     estimand = "country", arm = a, rep = 1L,
                     n_areas = length(k)), s)
      }
    }
    cat("loco done", target, on, "\n")
  }
}

RAW <- bind_rows(RAW, bind_rows(loco_rows))
write.csv(RAW, file.path(OUTDIR, "benchmarks_v2_raw.csv"), row.names = FALSE)

CELLS <- RAW |> group_by(country, outcome, target, estimand, arm) |>
  summarise(reps = n(), n_areas = max(n_areas),
            spearman = mean(spearman, na.rm = TRUE),
            spearman_sd = sd(spearman, na.rm = TRUE),
            spearman_lo = if (n() > 1) quantile(spearman, 0.1, na.rm = TRUE) else NA_real_,
            spearman_hi = if (n() > 1) quantile(spearman, 0.9, na.rm = TRUE) else NA_real_,
            pearson = mean(pearson, na.rm = TRUE),
            mae = mean(mae, na.rm = TRUE), wmae = mean(wmae, na.rm = TRUE),
            bias = mean(bias, na.rm = TRUE), topk = mean(topk, na.rm = TRUE),
            .groups = "drop")
write.csv(CELLS, file.path(OUTDIR, "benchmarks_v2_cells.csv"), row.names = FALSE)

SUMM <- CELLS |> group_by(estimand, target, arm) |>
  summarise(cells = n(),
            mean_spearman = round(mean(spearman, na.rm = TRUE), 3),
            median_spearman = round(median(spearman, na.rm = TRUE), 3),
            cells_positive = sum(spearman > 0, na.rm = TRUE),
            mean_wmae = round(mean(wmae, na.rm = TRUE), 2),
            mean_topk = round(mean(topk, na.rm = TRUE), 3), .groups = "drop") |>
  arrange(estimand, target, desc(mean_spearman))
write.csv(SUMM, file.path(OUTDIR, "benchmarks_v2_summary.csv"), row.names = FALSE)

cat("\n================ CORRECTED LEADERBOARD ================\n")
for (es in c("infill", "region", "country")) {
  if (!es %in% SUMM$estimand) next
  cat("\n---", ESTIMANDS_V2[[es]]$label, "---\n")
  print(as.data.frame(SUMM[SUMM$estimand == es, -1]), row.names = FALSE)
}
cat("\nNOTE: the null arm's correlation is mechanically negative (a leave-fold-out\n",
    "training mean is nearly constant across folds), so read its MAE, not its r.\n")
cat("\nDONE\n")
