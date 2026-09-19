# =============================================================================
# R/protocol_v2_sl.R   [SL-06]
#
# ONE CROSS-VALIDATED SUPERLEARNER OVER EVERY CANDIDATE, AS PROTOCOL-V2 ARMS
#
# Every candidate estimator the project has tried is scored inside one classic
# SuperLearner (fit_area_superlearner, R/area_superlearner.R) on the protocol's
# own cells, folds and targets, so the comparison between them is honest and
# paired: the same survey weights, the same Admin-1-blocked inner folds, the
# same rank-normalised predictor matrix with the domain map of the current
# harmonised vocabulary. The library:
#
#   mean                       constant
#   lasso / enet / ridge       penalised regression on the raw columns
#   ranger / xgb               forest and boosted trees on the raw columns
#   domain_index               zero-tuning index on the per-domain PCs
#   domain_pc_ridge / _enet    penalised regression on the per-domain PCs
#   domain_pc1_ols             OLS on the first PC of every domain
#   hapc / hapc_lasso          principal-component Highly Adaptive Ridge /
#                              Lasso on the raw columns (R/sl_hapc.R; off with
#                              V2_SL_HAPC=0)
#   spatial_gam                GAM on district centroids (in-country only)
#   spatial_plus_domain        GAM + elastic net on the domain PCs of the
#                              residual (in-country only)
#   mbg                        SPDE field + domain PC1s, INLA (in-country
#                              only; opt in with V2_SL_MBG=1, seconds per fit)
#
# Each learner rebuilds its own representation from the training rows it is
# handed, so nothing leaks across the outer fold. One SuperLearner fit per
# (cell, training fold) is memoised and serves every arm below:
#
#   sl_disc        the learner with the lowest cross-validated squared error
#   sl_nnls        the NNLS ensemble (squared-error loss)
#   sl_rank_disc   the learner with the highest cross-validated Spearman
#   sl_rank_nnls   NNLS on ranks (method.asl_rank, script 20)
#   sl_lrn_<name>  every library member on its own, refit on the full
#                  training fold -- the honest per-candidate comparison
#
# Registered into ARMS_V2 when sourced after R/protocol_v2.R (tar_source
# sorts by name). The arm reads from aux: domain_of (column -> domain),
# w (survey weights), rep_block (inner-fold blocks), lon / lat, estimand
# ("country" drops the spatial learners). Run through
# scripts/protocol_v2/56_weight_sources.R with V2_ARMS=<arms>.
# =============================================================================

.SL_BASE_LEARNERS_V2    <- c("mean", "lasso", "enet", "ridge", "ranger", "xgb")
.SL_DOMAIN_LEARNERS_V2  <- c("domain_index", "domain_pc_ridge", "domain_pc_enet", "domain_pc1_ols")
.SL_HAPC_LEARNERS_V2    <- c("hapc", "hapc_lasso")
.SL_SPATIAL_LEARNERS_V2 <- c("spatial_gam", "spatial_plus_domain", "mbg")

sl_library_v2 <- function(estimand = "infill", has_coords = TRUE) {
  lib <- c(.SL_BASE_LEARNERS_V2, .SL_DOMAIN_LEARNERS_V2)
  if (!identical(Sys.getenv("V2_SL_HAPC", "1"), "0") && exists("SL.hapc")) lib <- c(lib, .SL_HAPC_LEARNERS_V2)
  if (has_coords && !identical(estimand, "country")) {
    lib <- c(lib, "spatial_gam", "spatial_plus_domain")
    if (identical(Sys.getenv("V2_SL_MBG", "0"), "1") && exists(".mbg_ok") && .mbg_ok()) lib <- c(lib, "mbg")
  }
  lib
}
#' coordinates reach the spatial learners only
sl_screens_v2 <- function(lib) {
  sc <- stats::setNames(rep("screen.asl_no_coords", length(lib)), lib)
  sc["spatial_gam"] <- "screen.asl_coords"
  sc[intersect(lib, c("spatial_plus_domain", "mbg"))] <- "All"
  sc[lib]
}
SL_LEARNERS_V2 <- c(.SL_BASE_LEARNERS_V2, .SL_DOMAIN_LEARNERS_V2, .SL_HAPC_LEARNERS_V2, .SL_SPATIAL_LEARNERS_V2)

# ── memoised fit ─────────────────────────────────────────────────────────────
.sl_memo <- new.env(parent = emptyenv())
sl_memo_clear <- function() rm(list = ls(.sl_memo), envir = .sl_memo)
.sl_selection <- new.env(parent = emptyenv()); .sl_selection$rows <- list()
#' One row per SuperLearner fit: which learner each meta-learner chose and the
#' NNLS weights, for the selection-frequency table of the log.
sl_domain_selection_table <- function() {
  if (!length(.sl_selection$rows)) return(data.frame())
  dplyr::bind_rows(.sl_selection$rows)
}
sl_selection_clear <- function() .sl_selection$rows <- list()

.sl_fit_v2 <- function(tr, te, y, X, aux) {
  estimand <- if (is.null(aux$estimand)) "infill" else aux$estimand
  key <- paste(estimand, ncol(X), length(tr), paste(tr, collapse = ","), paste(te, collapse = ","),
               signif(sum(y[tr]), 12), sep = "|")
  hit <- .sl_memo[[key]]
  if (!is.null(hit)) return(hit)
  orig <- colnames(X)
  Xdf <- as.data.frame(X); names(Xdf) <- make.names(orig, unique = TRUE)
  dom <- if (is.null(aux$domain_of)) NULL else stats::setNames(unname(aux$domain_of[orig]), names(Xdf))
  has_coords <- !is.null(aux$lon) && !is.null(aux$lat) && all(is.finite(aux$lon[tr])) && all(is.finite(aux$lat[tr]))
  if (has_coords) { Xdf$lon <- aux$lon; Xdf$lat <- aux$lat }
  lib <- sl_library_v2(estimand, has_coords)
  fit <- fit_area_superlearner(
    Y = y[tr], X = Xdf[tr, , drop = FALSE], newX = Xdf[te, , drop = FALSE],
    weights = if (is.null(aux$w)) NULL else aux$w[tr],
    block = if (is.null(aux$rep_block)) NULL else aux$rep_block[tr],
    library = lib, screens = sl_screens_v2(lib), V = 5L, discrete = TRUE,
    meta = "mse", domain_of = dom)
  if (is.null(fit)) { .sl_memo[[key]] <- list(); return(list()) }
  lpn <- fit$library_pred_new
  ok <- apply(lpn, 2, function(z) all(is.finite(z)))
  coef <- fit$coef; coef[!is.finite(coef) | !ok] <- 0
  nnls <- if (sum(coef) > 0) as.numeric(lpn[, coef > 0, drop = FALSE] %*% coef[coef > 0]) else rep(mean(y[tr]), length(te))
  rk <- .asl_rank_coef(fit$Z, y[tr])
  rc <- rk$coef_nnls; rc[!ok] <- 0
  rank_nnls <- if (sum(rc) > 0) as.numeric(lpn[, rc > 0, drop = FALSE] %*% rc[rc > 0]) else rep(mean(y[tr]), length(te))
  rank_pick <- colnames(lpn)[which.max(ifelse(ok, rk$rho, -Inf))]
  preds <- list(sl_disc = fit$pred_new, sl_nnls = nnls,
                sl_rank_disc = lpn[, rank_pick], sl_rank_nnls = rank_nnls)
  for (nm in colnames(lpn)) preds[[paste0("sl_lrn_", nm)]] <- lpn[, nm]
  # selection record
  row <- data.frame(estimand = estimand, n_train = length(tr), n_test = length(te),
                    country = if (is.null(aux$country)) NA_character_ else aux$country,
                    outcome = if (is.null(aux$outcome)) NA_character_ else aux$outcome,
                    target = if (is.null(aux$target)) NA_character_ else aux$target,
                    rep = if (is.null(aux$rep)) NA_integer_ else aux$rep,
                    mse_pick = fit$pick, rank_pick = rank_pick, library = paste(lib, collapse = ","),
                    stringsAsFactors = FALSE)
  for (nm in colnames(lpn)) { row[[paste0("w_", nm)]] <- unname(coef[nm]); row[[paste0("wr_", nm)]] <- unname(rc[nm]); row[[paste0("rho_", nm)]] <- rk$rho[match(nm, colnames(lpn))] }
  .sl_selection$rows[[length(.sl_selection$rows) + 1L]] <- row
  .sl_memo[[key]] <- preds
  preds
}

.sl_pred_v2 <- function(tr, te, y, X, aux, what) {
  p <- .sl_fit_v2(tr, te, y, X, aux)[[what]]
  if (is.null(p) || length(p) != length(te)) rep(NA_real_, length(te)) else as.numeric(p)
}
arm_sl_disc_v2      <- function(tr, te, y, X, D, aux) .sl_pred_v2(tr, te, y, X, aux, "sl_disc")
arm_sl_nnls_v2      <- function(tr, te, y, X, D, aux) .sl_pred_v2(tr, te, y, X, aux, "sl_nnls")
arm_sl_rank_disc_v2 <- function(tr, te, y, X, D, aux) .sl_pred_v2(tr, te, y, X, aux, "sl_rank_disc")
arm_sl_rank_nnls_v2 <- function(tr, te, y, X, D, aux) .sl_pred_v2(tr, te, y, X, aux, "sl_rank_nnls")

SL_ENSEMBLE_ARMS <- c("sl_disc", "sl_nnls", "sl_rank_disc", "sl_rank_nnls")
SL_LEARNER_ARMS  <- paste0("sl_lrn_", SL_LEARNERS_V2)
SL_ARMS <- c(SL_ENSEMBLE_ARMS, SL_LEARNER_ARMS)

if (exists("ARMS_V2") && is.list(ARMS_V2)) {
  ARMS_V2$sl_disc      <- arm_sl_disc_v2
  ARMS_V2$sl_nnls      <- arm_sl_nnls_v2
  ARMS_V2$sl_rank_disc <- arm_sl_rank_disc_v2
  ARMS_V2$sl_rank_nnls <- arm_sl_rank_nnls_v2
  for (.nm in SL_LEARNERS_V2) {
    ARMS_V2[[paste0("sl_lrn_", .nm)]] <- local({
      what <- paste0("sl_lrn_", .nm)
      function(tr, te, y, X, D, aux) .sl_pred_v2(tr, te, y, X, aux, what)
    })
  }
  rm(.nm)
}
