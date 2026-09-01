# =============================================================================
# scripts/protocol_v2/04_fix3_normalisation.R
#
# FIX 3 in isolation: does within-country rank-normalisation before pooling
# change cross-country transport?
#
# The audit measured that 22 percent of the pooled cross-country covariate
# vocabulary carries >10x between-country mean offsets, plus pure country
# constants and Admin-1-broadcast columns, and that this enters the pooled
# LOCO fit uncentered. That is a scale defect: columns sharing a name but not
# a scale. Every arm in 02_run_benchmarks_v2.R already applies the fix, so the
# fix cannot be isolated there. This script runs leave-one-country-out
# transport three ways on IDENTICAL cells, folds and learner, varying only how
# predictors are put on a common scale:
#
#   raw_pooled     values as they come, pooled across countries (the defect)
#   z_pooled       z-scored over the POOLED rows (a global rescale; does not
#                  remove a between-country offset, it preserves it)
#   rank_within    rank-normalised WITHIN each country (the fix)
#
# Outcomes are within-country standardised in all three arms, so the contrast
# is purely about the predictor side. Transport is scored on rank correlation
# within each held-out country, for the reason given in 02: biomarker levels
# carry cross-survey offsets and a transported level is not a quantity this
# design can validate.
#
#   Rscript scripts/protocol_v2/04_fix3_normalisation.R
# -> results/tables/protocol_v2/fix3_normalisation.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
SEED <- 20260971L; set.seed(SEED)

#' Assemble a cell without any predictor scaling: raw values are kept
raw_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  X <- as.matrix(m[, PREDS, drop = FALSE])
  y <- m[[ycol]]
  list(country = cn, X = X, y_mod = if (target == "prev") .v2_logit(y) else y,
       y_nat = y, w = m[[wcol]], n = nrow(m))
}

#' Column-wise scaling schemes, all outcome-independent
scale_raw  <- function(X, ctry) { X[!is.finite(X)] <- NA; X }
scale_zall <- function(X, ctry) {
  mu <- colMeans(X, na.rm = TRUE)
  sdv <- apply(X, 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
}
scale_rank_within <- function(X, ctry) {
  out <- X
  for (cn in unique(ctry)) {
    k <- which(ctry == cn)
    out[k, ] <- apply(X[k, , drop = FALSE], 2, rank_normalize_v2)
  }
  out
}
SCHEMES <- list(raw_pooled = scale_raw, z_pooled = scale_zall,
                rank_within = scale_rank_within)

rows <- list()
for (target in c("prev", "level")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      cc <- tryCatch(raw_cell(cn, on, target), error = function(e) NULL)
      if (!is.null(cc)) cl[[cn]] <- cc
    }
    if (length(cl) < 3) next
    Xall  <- do.call(rbind, lapply(cl, function(z) z$X))
    ctry  <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat  <- unlist(lapply(cl, function(z) z$y_nat))
    wv    <- unlist(lapply(cl, function(z) z$w))
    # outcome standardised within country in EVERY arm
    Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    folds <- as.integer(factor(ctry))

    for (sch in names(SCHEMES)) {
      Xs <- SCHEMES[[sch]](Xall, ctry)
      # identical downstream handling for every scheme: drop unusable columns,
      # impute the rest to the column median, then build domain scores
      cov_j <- colMeans(is.finite(Xs))
      keep <- cov_j >= 0.70 & apply(Xs, 2, function(z)
        sum(is.finite(z)) > 2 && stats::sd(z[is.finite(z)]) > 0)
      keep[is.na(keep)] <- FALSE
      Xs <- Xs[, keep, drop = FALSE]
      med <- apply(Xs, 2, stats::median, na.rm = TRUE)
      for (j in seq_len(ncol(Xs))) Xs[!is.finite(Xs[, j]), j] <- med[j]
      colnames(Xs) <- PREDS[keep]
      if (ncol(Xs) < 20) next
      D <- build_domain_scores_v2(Xs, domain_of)
      # FAIRNESS: let glmnet standardise internally, so the raw_pooled arm is
      # judged on its information rather than punished for its units. Without
      # this a single penalty is applied to columns measured in thousands and
      # columns in [0,1] alike, and raw_pooled loses for a reason that has
      # nothing to do with cross-country comparability.
      aux <- list(lon = rep(NA_real_, length(Y)), lat = rep(NA_real_, length(Y)),
                  Admin1 = ctry, y_nat = Y, enet_standardize = TRUE)
      for (a in c("domain_index", "domain_enet")) {
        fn <- ARMS_V2[[a]]
        pred <- rep(NA_real_, length(Y))
        for (f in unique(folds)) {
          te <- which(folds == f); tr <- which(folds != f)
          if (length(tr) < 20) next
          p <- tryCatch(fn(tr, te, Y, Xs, D, aux),
                        error = function(e) rep(NA_real_, length(te)))
          if (length(p) == length(te)) pred[te] <- p
        }
        for (cn in unique(ctry)) {
          k <- which(ctry == cn)
          s <- score_v2(ynat[k], pred[k], wv[k],
                        scale = if (target == "prev") "prev" else "level")
          rows[[paste(target, on, sch, a, cn)]] <- data.frame(
            target = target, outcome = on, scheme = sch, arm = a,
            country = cn, n_areas = length(k),
            spearman = s$spearman, topk = s$topk)
        }
      }
    }
    cat("done", target, on, "\n")
  }
}

R3 <- bind_rows(rows)
write.csv(R3, file.path(OUTDIR, "fix3_normalisation.csv"), row.names = FALSE)

cat("\n=========== FIX 3: predictor scaling before pooling (LOCO transport) ===========\n")
summ <- R3 |> group_by(target, arm, scheme) |>
  summarise(cells = n(),
            mean_spearman = round(mean(spearman, na.rm = TRUE), 3),
            median_spearman = round(median(spearman, na.rm = TRUE), 3),
            positive = sum(spearman > 0, na.rm = TRUE),
            mean_topk = round(mean(topk, na.rm = TRUE), 3), .groups = "drop") |>
  arrange(target, arm, desc(mean_spearman))
print(as.data.frame(summ), row.names = FALSE)

cat("\n--- paired: rank_within minus raw_pooled, same cell/arm ---\n")
pw <- R3 |> filter(scheme %in% c("rank_within", "raw_pooled")) |>
  tidyr::pivot_wider(id_cols = c(target, outcome, arm, country),
                     names_from = scheme, values_from = spearman) |>
  filter(is.finite(rank_within), is.finite(raw_pooled))
cat("cells:", nrow(pw),
    "| mean gain:", round(mean(pw$rank_within - pw$raw_pooled), 3),
    "| rank_within better in:", sum(pw$rank_within > pw$raw_pooled), "of", nrow(pw), "\n")
cat("\nDONE\n")
