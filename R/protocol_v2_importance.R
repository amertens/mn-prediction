# =============================================================================
# R/protocol_v2_importance.R  [WS-02]
#
# WHAT DOES THE INDEX WEIGHT, IN TERMS OF THE ORIGINAL PREDICTORS?
#
# The zero-tuning index is linear in the rank-normalised predictors:
#     s_i = sum_k w_k d_ik,   d_ik = sum_j R*_jk (x_ij - c_j)
#  => s_i = sum_j beta_j x_ij + const,   beta_j = sum_k w_k R*_jk
# where R* is each domain's oriented PCA rotation (kept by build_domain_pcs_v2
# as attr(D, "basis")) and w the Fisher-z axis weights. Every x_j is on the
# same rank-normal scale, so beta_j (times the column's training SD, which is
# 1 minus imputation) is a comparable importance, and its sign is the
# direction. Two decompositions follow exactly:
#   share_j = beta_j cov(x_j, s) / var(s)     sums to 1 over predictors
#   share_g = cov(s_g, s) / var(s)            sums to 1 over domains
# (s_g the part of the index carried by domain g's axes).
#
# The sparse few-predictor arms use the same projection INSIDE the training
# fold: rank predictors by |beta_j| sd_j, keep the top m, and form an
# equal-weight sign-aligned composite (Dawes' improper linear model), or a
# beta-weighted one, rescaled to the training outcome like the index. Nothing
# is chosen on held-out rows. Registered into ARMS_V2 when sourced after
# R/protocol_v2.R and R/protocol_v2_weights.R:
#   sparse5, sparse10, sparse20   equal weights, m = 5 / 10 / 20
#   sparse10w                     beta weights, m = 10
# =============================================================================

#' Project axis weights back onto predictors. Returns a named numeric vector over colnames(X).
index_backproject_v2 <- function(w, basis, cols) {
  beta <- stats::setNames(numeric(length(cols)), cols)
  for (g in names(basis)) {
    b <- basis[[g]]
    k <- intersect(colnames(b$rot), names(w))
    if (!length(k)) next
    contrib <- as.numeric(b$rot[, k, drop = FALSE] %*% w[k])
    j <- match(b$cols, cols)
    ok <- is.finite(j)
    beta[j[ok]] <- beta[j[ok]] + contrib[ok]
  }
  beta
}

#' Importance table for one fit: the index on the training rows, projected onto predictors
#' @return list(columns = data.frame(column, beta, sd_tr, beta_std, share), domains = data.frame(domain, share, n_axes))
index_importance_v2 <- function(tr, y, X, D) {
  basis <- attr(D, "basis")
  if (is.null(basis)) stop("D carries no basis attribute; build it with build_domain_pcs_v2()")
  w <- .ws_z_pooled(tr, y, D)
  beta <- index_backproject_v2(w, basis, colnames(X))
  s <- as.numeric(X[tr, , drop = FALSE] %*% beta)
  vs <- stats::var(s)
  sd_tr <- apply(X[tr, , drop = FALSE], 2, stats::sd)
  share <- if (is.finite(vs) && vs > 0) beta * apply(X[tr, , drop = FALSE], 2, function(x) stats::cov(x, s)) / vs else beta * 0
  cols <- data.frame(column = names(beta), beta = as.numeric(beta), sd_tr = as.numeric(sd_tr),
                     beta_std = as.numeric(beta * sd_tr), share = as.numeric(share), stringsAsFactors = FALSE)
  dom <- do.call(rbind, lapply(names(basis), function(g) {
    k <- intersect(colnames(basis[[g]]$rot), names(w))
    sg <- if (length(k)) as.numeric(D[tr, k, drop = FALSE] %*% w[k]) else rep(0, length(tr))
    data.frame(domain = g, share = if (is.finite(vs) && vs > 0) stats::cov(sg, s) / vs else 0,
               n_axes = length(k), stringsAsFactors = FALSE)
  }))
  list(columns = cols, domains = dom)
}

# ── sparse few-predictor arms ─────────────────────────────────────────────────
.ws_sparse <- function(tr, te, y, X, D, m, weighted = FALSE) {
  basis <- attr(D, "basis")
  ytr <- y[tr]
  if (is.null(basis)) {                                   # fallback: marginal Fisher z of the raw columns
    beta <- apply(X[tr, , drop = FALSE], 2, function(x) .ws_fisher_z(x, ytr))
  } else {
    beta <- index_backproject_v2(.ws_z_pooled(tr, y, D), basis, colnames(X))
  }
  score <- abs(beta) * apply(X[tr, , drop = FALSE], 2, stats::sd)
  score[!is.finite(score)] <- 0
  if (all(score == 0)) return(rep(mean(ytr), length(te)))
  top <- order(score, decreasing = TRUE)[seq_len(min(m, sum(score > 0)))]
  wv <- stats::setNames(numeric(ncol(X)), colnames(X))
  wv[top] <- if (weighted) beta[top] else sign(beta[top])
  .ws_index_from_w(tr, te, y, X, wv)
}
arm_sparse5_v2   <- function(tr, te, y, X, D, aux) .ws_sparse(tr, te, y, X, D, 5)
arm_sparse10_v2  <- function(tr, te, y, X, D, aux) .ws_sparse(tr, te, y, X, D, 10)
arm_sparse20_v2  <- function(tr, te, y, X, D, aux) .ws_sparse(tr, te, y, X, D, 20)
arm_sparse10w_v2 <- function(tr, te, y, X, D, aux) .ws_sparse(tr, te, y, X, D, 10, weighted = TRUE)

WS_SPARSE_ARMS <- c("sparse5", "sparse10", "sparse20", "sparse10w")

if (exists("ARMS_V2")) {
  ARMS_V2$sparse5   <- arm_sparse5_v2
  ARMS_V2$sparse10  <- arm_sparse10_v2
  ARMS_V2$sparse20  <- arm_sparse20_v2
  ARMS_V2$sparse10w <- arm_sparse10w_v2
}
