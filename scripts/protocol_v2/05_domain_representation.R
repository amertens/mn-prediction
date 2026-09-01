# =============================================================================
# scripts/protocol_v2/05_domain_representation.R
#
# Two questions, one run.
#
# Q1  Does a richer within-domain representation beat the single sign-aligned
#     mean the domain index currently uses? Diagnostics show PC1 explains only
#     0.13 to 0.17 of the variance in the two largest domains (infant/child
#     morbidity, 37 members; agriculture, 93 members) while reaching 0.84 in
#     Ruralness - so a single number is well matched to some domains and badly
#     matched to the biggest ones.
#
# Q2  Does leave-one-country-out transport work under the corrected protocol,
#     and which representation makes it work best?
#
# REPRESENTATIONS COMPARED (all built on TRAINING rows only, then applied)
#   mean1   sign-aligned mean of members             (the current domain score)
#   pc1     first principal component score per domain
#   pc12    first two PCs per domain
#   pc123   first three PCs per domain
#   pcvar   as many PCs per domain as reach 80% of that domain's variance
#   sup1    supervised: training-correlation-weighted mean of members
#   raw     no domain structure at all: the full pooled predictor matrix
#
# Every representation is learned inside the fold. PCA rotations come from the
# training countries and are applied to the held-out country, so no held-out
# information reaches the basis. sup1 uses the training outcome only.
#
# Predictors are rank-normalised WITHIN country before pooling, and outcomes are
# standardised within country, so this is a RANKING comparison - the same scope
# limit as estimand C in 02b.
#
#   Rscript scripts/protocol_v2/05_domain_representation.R
# -> results/tables/protocol_v2/domain_representation.csv
# -> results/tables/protocol_v2/domain_representation_summary.csv
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
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
set.seed(20260981L)

cell_X <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y <- m[[ycol]]
  list(country = cn, n = nrow(m), X = Xr, w = m[[wcol]], y_nat = y,
       y_mod = if (target == "prev") .v2_logit(y) else y)
}

#' Build a domain representation from training rows, apply to all rows
make_rep <- function(X, tr, ytr, kind) {
  cols <- colnames(X)
  domains <- sort(unique(stats::na.omit(domain_of[cols])))
  if (kind == "raw") return(X)
  if (kind == "mean1") return(build_domain_scores_v2(X, domain_of, sign_rows = tr))
  blocks <- list()
  for (dm in domains) {
    cc <- cols[which(domain_of[cols] == dm)]
    if (length(cc) < 2) next
    M <- X[, cc, drop = FALSE]
    Mtr <- M[tr, , drop = FALSE]
    if (kind == "sup1") {
      r <- suppressWarnings(apply(Mtr, 2, function(z)
        if (stats::sd(z) == 0) 0 else stats::cor(z, ytr, method = "spearman")))
      r[!is.finite(r)] <- 0
      if (all(r == 0)) next
      b <- matrix(M %*% r / sum(abs(r)), ncol = 1)
      colnames(b) <- paste0(substr(dm, 1, 12), "_sup")
      blocks[[dm]] <- b
      next
    }
    pc <- tryCatch(stats::prcomp(Mtr, center = TRUE, scale. = FALSE),
                   error = function(e) NULL)
    if (is.null(pc)) next
    ve <- pc$sdev^2 / sum(pc$sdev^2)
    npc <- switch(kind,
                  pc1 = 1L, pc12 = 2L, pc123 = 3L,
                  pcvar = max(1L, which(cumsum(ve) >= 0.80)[1]))
    npc <- min(npc, ncol(pc$rotation), length(cc))
    if (!is.finite(npc) || npc < 1) next
    sc <- scale(M, center = pc$center, scale = FALSE) %*%
      pc$rotation[, seq_len(npc), drop = FALSE]
    colnames(sc) <- paste0(substr(dm, 1, 12), "_PC", seq_len(npc))
    blocks[[dm]] <- sc
  }
  if (!length(blocks)) return(NULL)
  out <- do.call(cbind, blocks)
  out[!is.finite(out)] <- 0
  out
}

KINDS <- c("mean1", "pc1", "pc12", "pc123", "pcvar", "sup1", "raw")
rows <- list()
for (target in c("prev", "level")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      cc <- tryCatch(cell_X(cn, on, target), error = function(e) NULL)
      if (!is.null(cc)) cl[[cn]] <- cc
    }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat))
    wv   <- unlist(lapply(cl, function(z) z$w))
    folds <- as.integer(factor(ctry))
    for (kind in KINDS) {
      for (learner in c("enet", "index")) {
        if (kind == "raw" && learner == "index") next   # index needs blocks
        pred <- rep(NA_real_, length(Y))
        npred <- NA_integer_
        for (f in unique(folds)) {
          te <- which(folds == f); tr <- which(folds != f)
          if (length(tr) < 20) next
          Z <- make_rep(Xm, tr, Y[tr], kind)
          if (is.null(Z) || ncol(Z) < 2) next
          npred <- ncol(Z)
          aux <- list(lon = rep(NA_real_, length(Y)), lat = rep(NA_real_, length(Y)),
                      Admin1 = ctry, y_nat = Y, enet_standardize = TRUE)
          p <- if (learner == "enet")
            .v2_enet(Z[tr, , drop = FALSE], Y[tr], Z[te, , drop = FALSE],
                     standardize = TRUE)
          else arm_domain_index_v2(tr, te, Y, Z, Z, aux)
          if (length(p) == length(te)) pred[te] <- p
        }
        for (cn in unique(ctry)) {
          k <- which(ctry == cn)
          s <- score_v2(ynat[k], pred[k], wv[k],
                        scale = if (target == "prev") "prev" else "level")
          rows[[paste(target, on, kind, learner, cn)]] <- data.frame(
            target = target, outcome = on, representation = kind,
            learner = learner, country = cn, n_areas = length(k),
            n_predictors = npred, spearman = s$spearman, topk = s$topk)
        }
      }
    }
    cat("done", target, on, "\n")
  }
}

R5 <- bind_rows(rows)
write.csv(R5, file.path(OUTDIR, "domain_representation.csv"), row.names = FALSE)

SUMM <- R5 |> group_by(target, learner, representation) |>
  summarise(cells = n(), n_pred = round(mean(n_predictors, na.rm = TRUE)),
            mean_spearman = round(mean(spearman, na.rm = TRUE), 3),
            median_spearman = round(median(spearman, na.rm = TRUE), 3),
            positive = sum(spearman > 0, na.rm = TRUE),
            mean_topk = round(mean(topk, na.rm = TRUE), 3), .groups = "drop") |>
  arrange(target, learner, desc(mean_spearman))
write.csv(SUMM, file.path(OUTDIR, "domain_representation_summary.csv"),
          row.names = FALSE)

cat("\n====== LOCO TRANSPORT BY DOMAIN REPRESENTATION ======\n")
print(as.data.frame(SUMM), row.names = FALSE)
cat("\n--- best overall per target ---\n")
print(as.data.frame(SUMM |> group_by(target) |> slice_max(mean_spearman, n = 3)),
      row.names = FALSE)
cat("\nDONE\n")
