# =============================================================================
# R/protocol_v2_hapc.R   [HP-02 / P8]
#
# PRINCIPAL-COMPONENT HIGHLY ADAPTIVE LASSO / RIDGE AS STANDALONE ARMS
#
# The hapc kernel learners (Wang, Schuler, van der Laan, Garcia Meixide; wrapped
# in R/sl_hapc.R) with the (tr, te, y, X, D, aux) interface of every protocol-v2
# arm, fitted on the RAW rank-normalised columns X — not the domain PCs D: HP-02
# showed the kernel's own reduction beats a second one in 11 of 12 comparisons,
# and SL-06 that PCHAL on the raw columns ties the index across borders on the
# level target. `hapc_lasso` (norm 1, PCHAL) is the pre-registered transport
# candidate P8 of docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md;
# `hapc_ridge` (norm 2, PCHAR) is kept beside it for the record. Degree 1,
# lambda by hapc's own inner 5-fold CV inside the training rows, survey weights
# not used (the package has no weights argument). Requires the reticulate venv
# of R/sl_hapc.R; if it is missing the arm returns the training mean and says
# so once. Registered into ARMS_V2 when sourced after R/protocol_v2.R.
#
#   Score (district rung, leave-one-country-out, 22 cells):
#   V2_ARMS=hapc_lasso,hapc_ridge V2_ESTIMANDS=country V2_OUT_TAG=_hapc \
#     Rscript -e "source('scripts/protocol_v2/56_weight_sources.R')"
# =============================================================================

.hapc_arm_env <- new.env(parent = emptyenv())
.hapc_arm_v2 <- function(tr, te, y, X, D, aux, norm) {
  fallback <- rep(mean(y[tr]), length(te))
  if (!exists("SL.hapc", mode = "function")) {
    if (!isTRUE(.hapc_arm_env$warned)) { cat("    [hapc arm] R/sl_hapc.R not sourced; returning the training mean\n"); .hapc_arm_env$warned <- TRUE }
    return(fallback)
  }
  if (is.null(X) || ncol(X) < 2 || length(tr) < 12) return(fallback)
  Xtr <- as.data.frame(X[tr, , drop = FALSE]); Xte <- as.data.frame(X[te, , drop = FALSE])
  r <- tryCatch(SL.hapc(y[tr], Xtr, Xte, stats::gaussian(), rep(1, length(tr)), seq_along(tr), norm = norm),
                error = function(e) { if (!isTRUE(.hapc_arm_env$warned_err)) { cat("    [hapc arm] fit failed:", conditionMessage(e), "\n"); .hapc_arm_env$warned_err <- TRUE }; NULL })
  if (is.null(r) || length(r$pred) != length(te) || !all(is.finite(r$pred))) return(fallback)
  as.numeric(r$pred)
}
arm_hapc_lasso_v2 <- function(tr, te, y, X, D, aux) .hapc_arm_v2(tr, te, y, X, D, aux, norm = "1")   # PCHAL
arm_hapc_ridge_v2 <- function(tr, te, y, X, D, aux) .hapc_arm_v2(tr, te, y, X, D, aux, norm = "2")   # PCHAR
HAPC_ARMS <- c("hapc_lasso", "hapc_ridge")

if (exists("ARMS_V2") && is.list(ARMS_V2)) {
  ARMS_V2$hapc_lasso <- arm_hapc_lasso_v2
  ARMS_V2$hapc_ridge <- arm_hapc_ridge_v2
}

# ── HP-03 sandbox arms (2026-09-19) ──────────────────────────────────────────
# Two questions from the P8 discussion:
#   (1) does a HAL with the DOMAIN STRUCTURE built in do better than PCHAL on
#       the raw columns? `hal_group`: an explicit zero-order HAL basis on the
#       domain PCs (indicator I(x >= knot) at the training deciles of every
#       axis, ~9 knots x ~100 axes) with a GROUP lasso over domains (gglasso,
#       lambda.min by 5-fold CV); `hal_lasso`: the same basis with a plain
#       lasso, to isolate the grouping. Degree 1 only (no interactions), which
#       is what hapc was run at.
#   (2) does the kernel learner gain from COORDINATES? `hapc_lasso_xy`: PCHAL
#       on the raw columns plus lon / lat; `hapc_xy_only`: PCHAL on lon / lat
#       alone (HAL as a spatial smoother, the kernel analogue of the GAM).
#       Coordinates only mean something in-country; under transport they
#       extrapolate, and the arm is scored there only to show it.
.hal_basis_v2 <- function(Dtr, Dall, n_knots = 9L) {
  cols <- colnames(Dtr); B <- list(); grp <- character(0)
  for (j in cols) {
    q <- unique(stats::quantile(Dtr[, j], probs = seq(0.1, 0.9, length.out = n_knots), names = FALSE, type = 7))
    q <- q[is.finite(q)]; if (!length(q)) next
    M <- sapply(q, function(k) as.numeric(Dall[, j] >= k)); if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    keep <- apply(M[seq_len(nrow(Dtr)), , drop = FALSE], 2, function(z) length(unique(z)) > 1)
    if (!any(keep)) next
    M <- M[, keep, drop = FALSE]; colnames(M) <- paste0(j, "_k", seq_len(ncol(M)))
    B[[j]] <- M; grp <- c(grp, rep(sub("_PC[0-9]+$", "", j), ncol(M)))
  }
  if (!length(B)) return(NULL)
  list(X = do.call(cbind, B), group = grp)
}
.hal_domain_arm_v2 <- function(tr, te, y, D, grouped) {
  fallback <- rep(mean(y[tr]), length(te))
  if (is.null(D) || ncol(D) < 2 || length(tr) < 12) return(fallback)
  b <- .hal_basis_v2(D[tr, , drop = FALSE], D)
  if (is.null(b) || ncol(b$X) < 2) return(fallback)
  Xtr <- b$X[tr, , drop = FALSE]; Xte <- b$X[te, , drop = FALSE]
  nf <- max(3, min(5, floor(length(tr) / 5)))
  if (grouped) {
    if (!requireNamespace("gglasso", quietly = TRUE)) return(fallback)
    g <- as.integer(factor(b$group)); o <- order(g)
    fit <- tryCatch(gglasso::cv.gglasso(Xtr[, o, drop = FALSE], y[tr], group = g[o], loss = "ls", pred.loss = "L2", nfolds = nf),
                    error = function(e) NULL)
    if (is.null(fit)) return(fallback)
    p <- tryCatch(as.numeric(stats::predict(fit$gglasso.fit, newx = Xte[, o, drop = FALSE], s = fit$lambda.min, type = "link")), error = function(e) NULL)
  } else {
    fit <- tryCatch(glmnet::cv.glmnet(Xtr, y[tr], alpha = 1, nfolds = nf, nlambda = 40, standardize = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(fallback)
    p <- tryCatch(as.numeric(stats::predict(fit, Xte, s = "lambda.min")), error = function(e) NULL)
  }
  if (is.null(p) || !all(is.finite(p))) fallback else p
}
arm_hal_group_v2 <- function(tr, te, y, X, D, aux) .hal_domain_arm_v2(tr, te, y, D, grouped = TRUE)
arm_hal_lasso_v2 <- function(tr, te, y, X, D, aux) .hal_domain_arm_v2(tr, te, y, D, grouped = FALSE)
.xy_v2 <- function(aux, n) if (!is.null(aux$lon) && !is.null(aux$lat) && all(is.finite(aux$lon)) && all(is.finite(aux$lat))) cbind(lon = aux$lon, lat = aux$lat) else NULL
arm_hapc_lasso_xy_v2 <- function(tr, te, y, X, D, aux) {
  xy <- .xy_v2(aux, nrow(X)); if (is.null(xy)) return(rep(mean(y[tr]), length(te)))
  .hapc_arm_v2(tr, te, y, cbind(X, xy), D, aux, norm = "1")
}
arm_hapc_xy_only_v2 <- function(tr, te, y, X, D, aux) {
  xy <- .xy_v2(aux, nrow(X)); if (is.null(xy)) return(rep(mean(y[tr]), length(te)))
  .hapc_arm_v2(tr, te, y, xy, D, aux, norm = "1")
}
HP03_ARMS <- c("hal_group", "hal_lasso", "hapc_lasso_xy", "hapc_xy_only")
if (exists("ARMS_V2") && is.list(ARMS_V2)) {
  ARMS_V2$hal_group     <- arm_hal_group_v2
  ARMS_V2$hal_lasso     <- arm_hal_lasso_v2
  ARMS_V2$hapc_lasso_xy <- arm_hapc_lasso_xy_v2
  ARMS_V2$hapc_xy_only  <- arm_hapc_xy_only_v2
}
