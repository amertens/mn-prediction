# =============================================================================
# R/protocol_v2_weights.R  [WS-01]
#
# WHERE SHOULD THE DOMAIN INDEX'S WEIGHTS COME FROM?
#
# The production index (arm_domain_index_v2) weights each domain axis by its
# pooled training-fold Fisher-z Spearman correlation with the outcome, scaled
# by sqrt(n - 3). That is the infinite-penalty limit of a ridge regression on
# the axes: beta(lambda) = (D'D + lambda I)^-1 D'y is proportional to D'y as
# lambda -> infinity, so "bivariate" means "maximal shrinkage of the joint fit
# toward the marginal one". This file supplies the other points of that family
# and three different SOURCES of evidence for the weights, all with the
# standard (tr, te, y, X, D, aux) arm interface so that scripts 02 / 02b /
# 56 can score them under identical folds:
#
#   index_soft1, index_soft2   soft-thresholded weights, sign(z)(|z| - c)_+,
#                              c = 1, 2: the sparse index. Falls back to the
#                              training mean (a constant, scored NA on rank)
#                              when no axis passes the threshold.
#   index_meta                 replication-weighted: the Fisher z is computed
#                              WITHIN each replication block of the training
#                              rows (aux$rep_block: the training countries in
#                              leave-one-country-out; Admin-1 regions with at
#                              least MIN_BLOCK training districts in-country)
#                              and combined by DerSimonian-Laird random-effects
#                              meta-analysis; the weight is the meta z, so an
#                              axis whose sign flips between blocks is pulled
#                              toward zero. With fewer than 2 usable blocks it
#                              is the production index.
#   index_meta_agree           index_meta with the weight set to zero unless
#                              every block agrees in sign.
#   index_decor                between-domain decorrelation: a global PCA of
#                              the axes on the training rows (components to
#                              95% of variance, at most n_tr - 2), then the
#                              production weighting on the orthogonal scores.
#                              Shared urban-rural variance is counted once.
#   ridge_min                  glmnet alpha = 0, lambda.min by inner random
#                              K-fold CV on squared error (the HP-01 ridge).
#   enet_1se, lasso_min, lasso_1se
#                              glmnet alpha = 0.5 / 1 at lambda.min or
#                              lambda.1se: the lasso "drops more variables",
#                              lambda.1se is the conventional reduced-tuning
#                              choice for small n.
#   ridge_nested, enet_nested, lasso_nested
#                              the penalty chosen by NESTED leave-one-block-out
#                              over the training rows (aux$rep_block), by the
#                              pooled inner out-of-block SPEARMAN rather than
#                              squared error, on a 30-value path; ties go to
#                              the larger penalty. This is "less CV, and CV on
#                              the decision metric".
#   ridge_nested_1se, enet_nested_1se, lasso_nested_1se
#                              the two fixes combined: nested block CV on rank,
#                              then the LARGEST penalty within one jackknife
#                              (over blocks) standard error of the best.
#
# Registered into ARMS_V2 when this file is sourced after R/protocol_v2.R
# (tar_source sorts by name: "protocol_v2_weights" sorts after "protocol_v2"
# and "protocol_v2_mbg"). aux$rep_block is optional: without it the meta arms
# reduce to the production index and the nested arms use aux$Admin1.
# =============================================================================

WS_MIN_BLOCK <- 5L      # smallest replication / inner block, in training districts

.ws_fisher_z <- function(x, y) {
  if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) return(0)
  r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  if (!is.finite(r)) return(0)
  r <- max(min(r, 0.999), -0.999)
  0.5 * log((1 + r) / (1 - r))
}

#' The production weights: Fisher z times sqrt(n - 3), pooled over the training rows
.ws_z_pooled <- function(tr, y, D) {
  ytr <- y[tr]; n <- length(tr)
  z <- apply(D[tr, , drop = FALSE], 2, function(x) .ws_fisher_z(x, ytr) * sqrt(max(n - 3, 1)))
  z[!is.finite(z)] <- 0
  z
}

#' Weighted sum of the axes, recentred and rescaled to the training outcome (step 5 of the index)
.ws_index_from_w <- function(tr, te, y, D, w) {
  ytr <- y[tr]
  if (all(w == 0)) return(rep(mean(ytr), length(te)))
  idx_tr <- as.numeric(D[tr, , drop = FALSE] %*% w)
  idx_te <- as.numeric(D[te, , drop = FALSE] %*% w)
  s <- stats::sd(idx_tr)
  if (!is.finite(s) || s == 0) return(rep(mean(ytr), length(te)))
  ((idx_te - mean(idx_tr)) / s) * stats::sd(ytr) + mean(ytr)
}

.ws_blocks <- function(tr, aux) {
  b <- aux$rep_block
  if (is.null(b)) b <- aux$Admin1
  if (is.null(b)) return(NULL)
  as.character(b)
}

# ── soft-thresholded weights ──────────────────────────────────────────────────
.ws_index_soft <- function(tr, te, y, D, cutoff) {
  z <- .ws_z_pooled(tr, y, D)
  w <- sign(z) * pmax(abs(z) - cutoff, 0)
  .ws_index_from_w(tr, te, y, D, w)
}
arm_index_soft1_v2 <- function(tr, te, y, X, D, aux) .ws_index_soft(tr, te, y, D, 1)
arm_index_soft2_v2 <- function(tr, te, y, X, D, aux) .ws_index_soft(tr, te, y, D, 2)

# ── replication-weighted (meta-analytic) weights ──────────────────────────────
#' DerSimonian-Laird random-effects meta z of per-block Fisher correlations
#' @return list(z = meta z per axis, agree = TRUE where every block shares a sign)
.ws_meta_z <- function(tr, y, D, blocks) {
  bt <- blocks[tr]
  tab <- table(bt)
  use <- names(tab)[tab >= WS_MIN_BLOCK]
  if (length(use) < 2) return(NULL)
  K <- ncol(D)
  Z <- matrix(0, nrow = length(use), ncol = K, dimnames = list(use, colnames(D)))
  V <- numeric(length(use))
  for (i in seq_along(use)) {
    rows <- tr[bt == use[i]]
    V[i] <- 1 / max(length(rows) - 3, 1)
    Z[i, ] <- apply(D[rows, , drop = FALSE], 2, function(x) .ws_fisher_z(x, y[rows]))
  }
  wf <- 1 / V
  zfe <- colSums(Z * wf) / sum(wf)
  Q <- colSums(wf * sweep(Z, 2, zfe, "-")^2)
  cval <- sum(wf) - sum(wf^2) / sum(wf)
  tau2 <- pmax(0, (Q - (length(use) - 1)) / cval)
  meta <- numeric(K); agree <- logical(K)
  for (k in seq_len(K)) {
    wr <- 1 / (V + tau2[k])
    mu <- sum(wr * Z[, k]) / sum(wr)
    se <- 1 / sqrt(sum(wr))
    meta[k] <- mu / se
    s <- sign(Z[, k]); s <- s[s != 0]
    agree[k] <- length(s) > 0 && (all(s > 0) || all(s < 0))
  }
  meta[!is.finite(meta)] <- 0
  list(z = meta, agree = agree)
}

.ws_index_meta <- function(tr, te, y, D, aux, gate) {
  blocks <- .ws_blocks(tr, aux)
  m <- if (is.null(blocks)) NULL else .ws_meta_z(tr, y, D, blocks)
  if (is.null(m)) return(.ws_index_from_w(tr, te, y, D, .ws_z_pooled(tr, y, D)))
  w <- m$z
  if (gate) w[!m$agree] <- 0
  .ws_index_from_w(tr, te, y, D, w)
}
arm_index_meta_v2       <- function(tr, te, y, X, D, aux) .ws_index_meta(tr, te, y, D, aux, gate = FALSE)
arm_index_meta_agree_v2 <- function(tr, te, y, X, D, aux) .ws_index_meta(tr, te, y, D, aux, gate = TRUE)

# ── between-domain decorrelation ──────────────────────────────────────────────
arm_index_decor_v2 <- function(tr, te, y, X, D, aux, var_target = 0.95) {
  Dtr <- D[tr, , drop = FALSE]
  keep <- apply(Dtr, 2, function(x) is.finite(stats::sd(x)) && stats::sd(x) > 0)
  if (sum(keep) < 2) return(rep(mean(y[tr]), length(te)))
  pc <- tryCatch(stats::prcomp(Dtr[, keep, drop = FALSE], center = TRUE, scale. = FALSE), error = function(e) NULL)
  if (is.null(pc)) return(rep(mean(y[tr]), length(te)))
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  npc <- max(1L, which(cumsum(ve) >= var_target)[1])
  npc <- min(npc, ncol(pc$rotation), max(1L, length(tr) - 2L))
  G <- scale(D[, keep, drop = FALSE], center = pc$center, scale = FALSE) %*% pc$rotation[, seq_len(npc), drop = FALSE]
  .ws_index_from_w(tr, te, y, G, .ws_z_pooled(tr, y, G))
}

# ── glmnet family ─────────────────────────────────────────────────────────────
.ws_memo <- new.env(parent = emptyenv())
ws_memo_clear <- function() rm(list = ls(.ws_memo), envir = .ws_memo)

#' One cv.glmnet per (alpha, training fold), shared by the lambda.min and lambda.1se arms
.ws_cv_fit <- function(tr, y, D, alpha) {
  key <- paste(alpha, ncol(D), length(tr), paste(tr, collapse = ","), signif(sum(y[tr]), 12), sep = "|")
  if (!is.null(.ws_memo[[key]])) return(.ws_memo[[key]])
  if (!requireNamespace("glmnet", quietly = TRUE) || ncol(D) < 2 || length(tr) < 12) return(NULL)
  nf <- max(3, min(5, floor(length(tr) / 5)))
  fit <- tryCatch(glmnet::cv.glmnet(D[tr, , drop = FALSE], y[tr], alpha = alpha, nfolds = nf,
                                    nlambda = 40, standardize = FALSE), error = function(e) NULL)
  .ws_memo[[key]] <- fit
  fit
}

.ws_glmnet_cv <- function(tr, te, y, D, alpha, s) {
  fit <- .ws_cv_fit(tr, y, D, alpha)
  if (is.null(fit)) return(rep(mean(y[tr]), length(te)))
  p <- tryCatch(as.numeric(stats::predict(fit, D[te, , drop = FALSE], s = s)), error = function(e) NULL)
  if (is.null(p) || !all(is.finite(p))) return(rep(mean(y[tr]), length(te)))
  p
}
arm_ridge_min_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_cv(tr, te, y, D, 0,   "lambda.min")
arm_enet_1se_v2  <- function(tr, te, y, X, D, aux) .ws_glmnet_cv(tr, te, y, D, 0.5, "lambda.1se")
arm_lasso_min_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_cv(tr, te, y, D, 1,   "lambda.min")
arm_lasso_1se_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_cv(tr, te, y, D, 1,   "lambda.1se")

#' Penalty by nested leave-one-block-out on the pooled inner out-of-block Spearman
.ws_glmnet_nested <- function(tr, te, y, D, aux, alpha, nlambda = 30, rule = c("max", "1se")) {
  rule <- match.arg(rule)
  if (!requireNamespace("glmnet", quietly = TRUE) || ncol(D) < 2 || length(tr) < 12) return(rep(mean(y[tr]), length(te)))
  blocks <- .ws_blocks(tr, aux)
  bt <- if (is.null(blocks)) NULL else blocks[tr]
  if (is.null(bt) || length(unique(bt)) < 3) {                  # fall back to 3 random district folds
    set.seed(20260956L + length(tr)); bt <- as.character(sample(rep(1:3, length.out = length(tr))))
  }
  full <- tryCatch(glmnet::glmnet(D[tr, , drop = FALSE], y[tr], alpha = alpha, nlambda = nlambda, standardize = FALSE),
                   error = function(e) NULL)
  if (is.null(full)) return(rep(mean(y[tr]), length(te)))
  lam <- full$lambda
  P <- matrix(NA_real_, nrow = length(tr), ncol = length(lam))
  for (b in unique(bt)) {
    inn <- which(bt != b); out <- which(bt == b)
    if (length(inn) < 8 || !length(out)) next
    f <- tryCatch(glmnet::glmnet(D[tr[inn], , drop = FALSE], y[tr[inn]], alpha = alpha, lambda = lam, standardize = FALSE),
                  error = function(e) NULL)
    if (is.null(f)) next
    p <- tryCatch(stats::predict(f, D[tr[out], , drop = FALSE], s = lam), error = function(e) NULL)
    if (is.null(p) || ncol(p) != length(lam)) next
    P[out, ] <- p
  }
  ok <- rowSums(is.finite(P)) == length(lam)
  if (sum(ok) < 8) return(rep(mean(y[tr]), length(te)))
  pooled_sp <- function(rows) apply(P[rows, , drop = FALSE], 2, function(p) if (stats::sd(p) == 0) -Inf else suppressWarnings(stats::cor(p, y[tr][rows], method = "spearman")))
  sc <- pooled_sp(which(ok)); sc[!is.finite(sc)] <- -Inf
  if (all(sc == -Inf)) return(rep(mean(y[tr]), length(te)))
  if (rule == "1se") {
    # jackknife over blocks: the SE of the pooled inner Spearman at each penalty,
    # then the LARGEST penalty whose score is within one SE of the best
    bl <- unique(bt[ok]); J <- sapply(bl, function(b) pooled_sp(which(ok & bt != b)))
    if (is.matrix(J) && ncol(J) >= 2) {
      J[!is.finite(J)] <- NA
      se <- sqrt((ncol(J) - 1) / ncol(J) * rowSums((J - rowMeans(J, na.rm = TRUE))^2, na.rm = TRUE))
      thr <- max(sc) - se[which.max(sc)]
      best <- lam[which(sc >= thr)[1]]                          # lambda is decreasing: first index = largest penalty
    } else best <- lam[which.max(sc)]
  } else best <- lam[which.max(sc)]                             # first max = largest penalty among ties
  p <- tryCatch(as.numeric(stats::predict(full, D[te, , drop = FALSE], s = best)), error = function(e) NULL)
  if (is.null(p) || !all(is.finite(p))) return(rep(mean(y[tr]), length(te)))
  p
}
arm_ridge_nested_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 0)
arm_enet_nested_v2  <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 0.5)
arm_lasso_nested_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 1)
arm_ridge_nested_1se_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 0,   rule = "1se")
arm_enet_nested_1se_v2  <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 0.5, rule = "1se")
arm_lasso_nested_1se_v2 <- function(tr, te, y, X, D, aux) .ws_glmnet_nested(tr, te, y, D, aux, 1,   rule = "1se")

WS_ARMS <- c("index_soft1", "index_soft2", "index_meta", "index_meta_agree", "index_decor",
             "ridge_min", "ridge_nested", "enet_1se", "enet_nested", "lasso_min", "lasso_1se", "lasso_nested",
             "ridge_nested_1se", "enet_nested_1se", "lasso_nested_1se")

if (exists("ARMS_V2")) {
  ARMS_V2$index_soft1      <- arm_index_soft1_v2
  ARMS_V2$index_soft2      <- arm_index_soft2_v2
  ARMS_V2$index_meta       <- arm_index_meta_v2
  ARMS_V2$index_meta_agree <- arm_index_meta_agree_v2
  ARMS_V2$index_decor      <- arm_index_decor_v2
  ARMS_V2$ridge_min        <- arm_ridge_min_v2
  ARMS_V2$ridge_nested     <- arm_ridge_nested_v2
  ARMS_V2$enet_1se         <- arm_enet_1se_v2
  ARMS_V2$enet_nested      <- arm_enet_nested_v2
  ARMS_V2$lasso_min        <- arm_lasso_min_v2
  ARMS_V2$lasso_1se        <- arm_lasso_1se_v2
  ARMS_V2$lasso_nested     <- arm_lasso_nested_v2
  ARMS_V2$ridge_nested_1se <- arm_ridge_nested_1se_v2
  ARMS_V2$enet_nested_1se  <- arm_enet_nested_1se_v2
  ARMS_V2$lasso_nested_1se <- arm_lasso_nested_1se_v2
}
