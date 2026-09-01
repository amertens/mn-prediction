# =============================================================================
# scripts/protocol_v2/02_run_benchmarks_v2.R
#
# FIX 1, 3, 4, 5: the corrected leaderboard.
#
#   FIX 1  every scheme with a random component is REPLICATED (R draws) and
#          summarised over draws; the region scheme is exhaustive and therefore
#          has no draw luck at all. Scoring is precision-weighted by n_eff.
#   FIX 3  predictors are rank-normalised WITHIN country before anything else,
#          which is what makes cross-country pooling meaningful.
#   FIX 4  the primary predictor representation is 18 domain scores, not 373
#          raw columns; a raw-column arm is retained to measure the difference.
#   FIX 5  three estimands are scored separately, each against baselines that
#          saw the same information. In particular the covariate-free
#          comparator for in-fill is the JACKKNIFED regional mean, never the
#          withdrawn arm that reads the held-out district's own respondents.
#
# WHAT IS DELIBERATELY NOT CLAIMED
# --------------------------------
# Estimand C (transport to an unsurveyed country) is scored on RANK metrics
# only. Biomarker levels carry large cross-survey offsets - the project's own
# established finding - so a transported level is not a quantity this design
# can validate. Outcomes are within-country standardised before pooling for
# exactly that reason, which is the outcome-side twin of fix 3.
#
#   Rscript scripts/protocol_v2/02_run_benchmarks_v2.R
#   V2_REPS=20  override the number of in-fill replications (default 10)
#   PROFILE=smoke  Ghana only, 3 reps
# -> results/tables/protocol_v2/benchmarks_v2_raw.csv    one row per draw
# -> results/tables/protocol_v2/benchmarks_v2_cells.csv  summarised over draws
# -> results/tables/protocol_v2/benchmarks_v2_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR  <- "results/tables/protocol_v2"
PROFILE <- Sys.getenv("PROFILE", "full")
REPS    <- as.integer(Sys.getenv("V2_REPS", if (PROFILE == "smoke") "3" else "10"))
SEED    <- 20260961L
set.seed(SEED)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")

COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")
if (PROFILE == "smoke") COUNTRIES <- COUNTRIES["ghana"]

# V2_COUNTRY runs estimands A and B for one country only, writing a shard, so
# the four countries can run as parallel processes. Estimand C needs every
# country at once and is therefore skipped in shard mode; 02b merges the shards
# and runs it. Shards are deterministic: the fold seed depends on rep_id only.
.v2_args <- commandArgs(trailingOnly = TRUE)
SHARD <- if (length(.v2_args)) .v2_args[1] else Sys.getenv("V2_COUNTRY", "")
if (nzchar(SHARD)) COUNTRIES <- COUNTRIES[SHARD]
SUF <- if (nzchar(SHARD)) paste0("_", SHARD) else ""

# ── centroids on the pair key ───────────────────────────────────────────────
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

#' Assemble one cell: outcome, predictors (fix 3), domain scores (fix 4)
#' Predictors are prepared on the country's SURVEYED rows only, which is the
#' analysis set; rank-normalisation and imputation use no outcome information.
build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  ncol_eff <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[ncol_eff]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))          # PAIR KEY
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")],
                  by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))   # FIX 3
  if (ncol(Xr) < 20) return(NULL)
  D  <- domain_representation_v2(Xr, domain_of)                   # FIX 4
  y_nat <- m[[ycol]]
  y_mod <- if (target == "prev") .v2_logit(y_nat) else y_nat
  list(country = cn, outcome = on, target = target,
       y_nat = y_nat, y_mod = y_mod, X = Xr, D = D,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1,
                  y_nat = y_nat),
       w = m[[ncol_eff]], Admin1 = m$Admin1, n = nrow(m))
}

#' Run every arm over one fold assignment and score the pooled out-of-fold set
run_one_draw <- function(cell, estimand, folds, rep_id) {
  arms <- arms_for_estimand_v2(estimand)
  out <- list()
  for (a in arms) {
    fn <- ARMS_V2[[a]]
    pred_mod <- rep(NA_real_, cell$n)
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12 || !length(te)) next
      p <- tryCatch(fn(tr, te, cell$y_mod, cell$X, cell$D, cell$aux),
                    error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) pred_mod[te] <- p
    }
    pred <- if (cell$target == "prev") .v2_expit(pred_mod) else pred_mod
    # The null is defined on the NATURAL scale (the arithmetic training mean),
    # so that back-transforming a logit-scale mean cannot handicap it.
    if (a == "null_train_mean") {
      pred <- rep(NA_real_, cell$n)
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12 || !length(te)) next
        pred[te] <- mean(cell$y_nat[tr])
      }
    }
    s <- score_v2(cell$y_nat, pred, cell$w,
                  scale = if (cell$target == "prev") "prev" else "level")
    out[[a]] <- cbind(data.frame(country = cell$country, outcome = cell$outcome,
                                 target = cell$target, estimand = estimand,
                                 arm = a, rep = rep_id, n_areas = cell$n),
                      s)
  }
  bind_rows(out)
}

# ── estimands A and B, within country ───────────────────────────────────────
rows <- list()
cells_index <- TG |> distinct(country, outcome) |>
  filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells_index))) {
  cn <- cells_index$country[i]; on <- cells_index$outcome[i]
  for (target in c("prev", "level")) {
    cell <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
    if (is.null(cell)) next
    # A. in-fill: replicated K-fold over districts  (FIX 1)
    for (r in seq_len(REPS)) {
      folds <- make_folds_v2("kfold_district", cell$n, k = 5, rep_id = r)
      rows[[paste(cn, on, target, "infill", r)]] <-
        run_one_draw(cell, "infill", folds, r)
    }
    # B. region extrapolation: exhaustive LORO, deterministic
    if (dplyr::n_distinct(cell$Admin1) >= 3) {
      folds <- make_folds_v2("loro", cell$n, blocks = cell$Admin1)
      rows[[paste(cn, on, target, "region", 1)]] <-
        run_one_draw(cell, "region", folds, 1L)
    }
    cat("done", cn, on, target, "\n")
  }
}

# ── estimand C, transport across countries ──────────────────────────────────
# Outcomes are within-country standardised before pooling (the outcome-side
# twin of fix 3), so only rank metrics are meaningful and only those are read.
loco_rows <- list()
if (PROFILE != "smoke" && !nzchar(SHARD)) {
  for (target in c("prev", "level")) {
    for (on in unique(TG$outcome)) {
      cl <- list()
      for (cn in COUNTRIES) {
        cc <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
        if (!is.null(cc)) cl[[cn]] <- cc
      }
      if (length(cl) < 3) next
      common <- Reduce(intersect, lapply(cl, function(z) colnames(z$D)))
      if (length(common) < 5) next
      Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
      Dm <- do.call(rbind, lapply(cl, function(z) z$D[, common, drop = FALSE]))
      ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
      ynat <- unlist(lapply(cl, function(z) z$y_nat))
      wv <- unlist(lapply(cl, function(z) z$w))
      aux <- list(lon = unlist(lapply(cl, function(z) z$aux$lon)),
                  lat = unlist(lapply(cl, function(z) z$aux$lat)),
                  Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))),
                  y_nat = Y)
      Xdummy <- Dm
      folds <- as.integer(factor(ctry))
      for (a in c("null_train_mean", "domain_index", "domain_enet")) {
        fn <- ARMS_V2[[a]]
        pred <- rep(NA_real_, length(Y))
        for (f in unique(folds)) {
          te <- which(folds == f); tr <- which(folds != f)
          if (length(tr) < 20) next
          p <- tryCatch(fn(tr, te, Y, Xdummy, Dm, aux),
                        error = function(e) rep(NA_real_, length(te)))
          if (length(p) == length(te)) pred[te] <- p
        }
        # score WITHIN each held-out country: transport is a ranking claim
        for (cn in unique(ctry)) {
          k <- which(ctry == cn)
          s <- score_v2(ynat[k], pred[k], wv[k],
                        scale = if (target == "prev") "prev" else "level")
          s$mae <- NA_real_; s$wmae <- NA_real_; s$bias <- NA_real_
          s$rmse_sd <- NA_real_          # level is not transportable: not claimed
          loco_rows[[paste(target, on, a, cn)]] <- cbind(
            data.frame(country = cn, outcome = on, target = target,
                       estimand = "country", arm = a, rep = 1L,
                       n_areas = length(k)), s)
        }
      }
      cat("loco done", target, on, "\n")
    }
  }
}

RAW <- bind_rows(c(rows, loco_rows))
write.csv(RAW, file.path(OUTDIR, paste0("benchmarks_v2_raw", SUF, ".csv")),
          row.names = FALSE)
if (nzchar(SHARD)) { cat("shard", SHARD, "written\n"); quit(save = "no") }

# ── summarise over draws (FIX 1: never report a single draw) ────────────────
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
for (es in unique(SUMM$estimand)) {
  cat("\n---", ESTIMANDS_V2[[es]]$label, "---\n")
  print(as.data.frame(SUMM[SUMM$estimand == es, -1]), row.names = FALSE)
}
cat("\nDONE\n")
