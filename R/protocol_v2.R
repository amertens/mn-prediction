# =============================================================================
# R/protocol_v2.R
#
# Corrected evaluation protocol ("v2"). Reusable functions for the five fixes
# adopted after the 2026-09-01 methods audit and signal-probe battery
# (docs/findings/TWO_READINGS_2026-09c.md):
#
#   1 replicated folds + effective sample size    make_folds_v2(), effective_n_v2()
#   2 continuous biomarker as a first-class target build_targets_v2()
#   3 within-country rank-normalisation before any pooling  rank_normalize_v2()
#   4 domain scores instead of raw columns        build_domain_scores_v2()
#   5 three estimands, information-matched nulls   ESTIMANDS_V2, arm_*()
#
# NOTHING HERE IS IN THE TARGETS DAG. These functions are auto-sourced by
# tar_source("R/") but no target depends on them, so adding this file does not
# invalidate any cached target. The drivers live in scripts/protocol_v2/.
#
# WHY A LAYER RATHER THAN AN EDIT TO R/admin2_analysis.R
# ------------------------------------------------------
# The upstream fix for n_svy semantics is a two-line addition to
# compute_svy_admin2() (carry n_eff alongside n_svy). That edit invalidates the
# whole DAG and costs a multi-hour rebuild, so it is deliberately NOT made here;
# this layer recomputes the corrected quantities from the cached store instead.
# When a rebuild is next acceptable, move effective_n_v2() upstream and delete
# the recomputation.
# =============================================================================

# ── small helpers ────────────────────────────────────────────────────────────

.v2_num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

.v2_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

.v2_wvar <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  m <- .v2_wmean(x, w)
  v <- sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
  v * sum(ok) / max(1, sum(ok) - 1)
}

#' Logit with an explicit, reported clamp
#' @param eps clamp; districts pinned at the clamp are counted by the caller
.v2_logit <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.v2_expit <- function(z) stats::plogis(z)

# ── FIX 1. Effective sample size ─────────────────────────────────────────────

#' Kish effective sample size from survey weights
#'
#' `sum(w)^2 / sum(w^2)`. This is the weighting component only; it does NOT
#' account for clustering. Combine with a design effect via `effective_n_v2()`.
kish_n_v2 <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(0)
  sum(w)^2 / sum(w^2)
}

#' Total design effect for one country x outcome, estimated where it is estimable
#'
#' The audit established that svy_prev_se is degenerate in the majority of
#' districts (Ghana 64 of 75, Malawi 74 of 87, Gambia 17 of 30) because those
#' districts contain a single PSU, so a district-level design effect cannot be
#' estimated. This estimates ONE deff per country x outcome at the national
#' level, where every country has 60-103 PSUs, and applies it to districts.
#'
#' deff = Var_design(estimate) / Var_SRS(estimate), computed with a
#' cluster-robust (linearisation) variance over PSUs.
#'
#' @param y numeric outcome (0/1 for a proportion, continuous otherwise)
#' @param w survey weights
#' @param psu PSU / cluster identifier
#' @return list(deff, n_raw, n_kish, method); deff is floored at 1 only when the
#'   estimate is degenerate, never silently
deff_national_v2 <- function(y, w, psu) {
  ok <- is.finite(y) & is.finite(w) & w > 0 & !is.na(psu)
  y <- y[ok]; w <- w[ok]; psu <- as.character(psu)[ok]
  n <- length(y)
  out <- list(deff = NA_real_, n_raw = n, n_kish = kish_n_v2(w),
              n_psu = length(unique(psu)), method = "unestimable")
  if (n < 10 || out$n_psu < 4) return(out)

  # Cluster-robust variance of the weighted mean (Hansen-Hurwitz / linearised).
  m   <- .v2_wmean(y, w)
  W   <- sum(w)
  # per-PSU totals of the residual-weighted contributions
  z   <- w * (y - m)
  agg <- tapply(z, psu, sum)
  k   <- length(agg)
  v_design <- (k / (k - 1)) * sum((agg - mean(agg))^2) / W^2
  v_srs    <- .v2_wvar(y, w) / n
  if (!is.finite(v_design) || !is.finite(v_srs) || v_srs <= 0) return(out)

  deff <- v_design / v_srs
  # A deff below 1 is possible under stratification but is usually noise at
  # these sizes; keep the raw value and let the caller see it.
  out$deff <- deff
  out$method <- "cluster_robust_national"
  out
}

#' Effective sample size for a district
#'
#' n_eff = n_raw / deff. Uses the country x outcome deff from
#' `deff_national_v2()`. Returns n_raw when deff is unestimable, flagged.
effective_n_v2 <- function(n_raw, deff, fallback_deff = 1.5) {
  d <- ifelse(is.finite(deff) & deff > 0, deff, fallback_deff)
  pmax(1, n_raw / d)
}

#' Intra-cluster correlation implied by a national design effect
#'
#' AU-01 finding 6 (2026-09-07): dividing every district by the national deff
#' weights a three-PSU district like a single-PSU one. The national total deff
#' is split Kish-fashion into a weighting part, deff_w = n_raw / n_kish, and a
#' clustering part, deff_c = deff / deff_w = 1 + (b - 1) rho with b = n_raw /
#' n_psu the mean PSU take; rho is what transfers to a district, whose own
#' PSU count and Kish n then give its design effect (effective_n_district_v2).
#' @return rho clipped to [0, 0.95], or NA when the national deff is unestimable
icc_from_deff_v2 <- function(deff, n_raw, n_psu, n_kish) {
  if (!is.finite(deff) || deff <= 0 || !is.finite(n_raw) || !is.finite(n_psu) || n_psu < 2 ||
      !is.finite(n_kish) || n_kish <= 0) return(NA_real_)
  b <- n_raw / n_psu
  if (b <= 1) return(NA_real_)
  deff_w <- n_raw / n_kish
  deff_c <- deff / deff_w
  min(max((deff_c - 1) / (b - 1), 0), 0.95)
}

#' Effective sample size of a district from its own PSU count
#'
#' n_eff_d = n_kish_d / (1 + (n_d / k_d - 1) rho): the district's Kish n for the
#' weighting part and its own mean PSU take for the clustering part. Falls back
#' to effective_n_v2() (national deff) where rho is NA.
effective_n_district_v2 <- function(n_raw, n_kish, n_psu, rho, deff_national = NA_real_, fallback_deff = 1.5) {
  b <- ifelse(is.finite(n_psu) & n_psu > 0, n_raw / pmax(n_psu, 1), n_raw)
  deff_c <- 1 + pmax(b - 1, 0) * rho
  out <- if (is.finite(rho)) pmax(1, n_kish / deff_c) else rep(NA_real_, length(n_raw))
  fb <- effective_n_v2(n_raw, deff_national, fallback_deff)
  ifelse(is.finite(out), out, fb)
}

# ── FIX 3. Within-country rank-normalisation ────────────────────────────────

#' Rank-normal transform (van der Waerden), outcome-independent
#'
#' Applied WITHIN country before any pooling. This is fix 3: the audit found
#' that 22 percent of the pooled cross-country covariate set carries >10x
#' between-country mean offsets and enters an uncentered pooled fit, which is
#' a scale defect rather than a signal problem. Ranking within country removes
#' the offset while preserving within-country ordering. It uses no outcome
#' information, so it is safe to apply before folds are drawn.
rank_normalize_v2 <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && stats::sd(x[ok]) > 0) {
    out[ok] <- stats::qnorm((rank(x[ok]) - 0.5) / sum(ok))
  }
  out
}

#' Rank-normalise a predictor matrix within country, then median-impute
#'
#' Median imputation on the rank-normal scale is imputation to 0. Columns below
#' `min_cov` coverage are dropped rather than imputed. Both steps are
#' outcome-independent (fix 4's companion to complete-case removal, which the
#' audit found silently deletes Ghana's entire 153-column DHS block).
prep_predictors_v2 <- function(X, min_cov = 0.70) {
  X <- as.matrix(X)
  Xr <- apply(X, 2, rank_normalize_v2)
  if (is.null(dim(Xr))) Xr <- matrix(Xr, nrow = nrow(X))
  cov_j <- colMeans(is.finite(Xr))
  keep <- cov_j >= min_cov &
    apply(Xr, 2, function(z) sum(is.finite(z)) > 2 &&
            stats::sd(z[is.finite(z)]) > 0)
  keep[is.na(keep)] <- FALSE
  Xr <- Xr[, keep, drop = FALSE]
  Xr[!is.finite(Xr)] <- 0
  colnames(Xr) <- colnames(X)[keep]
  Xr
}

# ── FIX 4. Domain scores ────────────────────────────────────────────────────

#' Sign-aligned domain scores from rank-normalised predictors
#'
#' Fix 4: the audit measured the effective dimensionality of the 373-predictor
#' matrix at roughly 30-55 (participation ratio ~9-10), and the project's own
#' ablation shows the top-20 screen makes most predictors unreachable. Domain
#' scores collapse each of the 18 conceptual domains to one number, orienting
#' members by the sign of the first principal component so that members do not
#' cancel.
#'
#' Signs are learned from `sign_rows` only (pass training-country rows for a
#' leave-one-country-out fit) and never from the outcome.
#'
#' @param Xr rank-normalised predictor matrix (rows = areas)
#' @param domain_of named character vector mapping column -> domain
#' @param sign_rows logical/integer index of rows used to learn PC1 signs
build_domain_scores_v2 <- function(Xr, domain_of, sign_rows = NULL) {
  cols <- colnames(Xr)
  domains <- sort(unique(stats::na.omit(domain_of[cols])))
  if (is.null(sign_rows)) sign_rows <- seq_len(nrow(Xr))
  sgn <- stats::setNames(rep(1, length(cols)), cols)
  for (dm in domains) {
    cc <- cols[which(domain_of[cols] == dm)]
    if (length(cc) < 2) next
    M <- Xr[sign_rows, cc, drop = FALSE]
    M[!is.finite(M)] <- 0
    if (nrow(M) < 3 || all(apply(M, 2, stats::sd) == 0)) next
    pc <- tryCatch(stats::prcomp(M, center = TRUE), error = function(e) NULL)
    if (is.null(pc)) next
    ld <- pc$rotation[, 1]
    if (mean(ld > 0) < 0.5) ld <- -ld   # majority-positive orientation
    sgn[cc] <- ifelse(ld >= 0, 1, -1)
  }
  if (!length(domains)) return(matrix(numeric(0), nrow = nrow(Xr), ncol = 0))
  # NOTE: this single-score representation is retained for comparability with
  # earlier runs, but it is NOT the default any more. Measured over 22
  # leave-one-country-out cells it is the WORST of seven representations
  # (index 0.151 on the level target); build_domain_pcs_v2() reaches 0.255.
  # See docs/findings/PROTOCOL_V2.md, "domain representation".
  S <- sweep(Xr, 2, sgn[cols], "*")
  D <- matrix(0, nrow = nrow(Xr), ncol = length(domains),
              dimnames = list(NULL, domains))
  for (dm in domains) {
    cc <- which(domain_of[cols] == dm)
    if (!length(cc)) next
    D[, dm] <- rowMeans(S[, cc, drop = FALSE], na.rm = TRUE)
  }
  D[!is.finite(D)] <- 0
  D
}

#' Principal-component representation of each domain — THE DEFAULT
#'
#' One score per domain throws away most of the largest domains: PC1 explains
#' 0.13 of the variance in infant/child morbidity (37 members) and 0.17 in
#' agriculture (93), against 0.84 in Ruralness. Keeping as many components as
#' reach `var_target` of each domain's variance adapts the representation to
#' the domain instead of assuming every domain is one factor.
#'
#' Measured over 22 leave-one-country-out cells, mean Spearman on the level
#' target: sign-aligned mean 0.151, PC1 0.218, PC1-2 0.244, PCs-to-80% 0.255.
#' Supervised weighting adds nothing (0.196), so the gain comes from
#' representing more of each domain's variance, not from aiming it at the
#' outcome — which also keeps this step outcome-independent and leak-free.
#'
#' Rotations are learned from `sign_rows` (pass the training rows) and applied
#' to every row, so no held-out information reaches the basis.
#'
#' @param var_target cumulative variance share to retain per domain
#' @param max_pc hard cap per domain, so one huge domain cannot dominate
build_domain_pcs_v2 <- function(Xr, domain_of, sign_rows = NULL,
                                var_target = 0.80, max_pc = 12L) {
  cols <- colnames(Xr)
  domains <- sort(unique(stats::na.omit(domain_of[cols])))
  if (!length(domains)) return(matrix(numeric(0), nrow = nrow(Xr), ncol = 0))
  # PC columns are named make.names(substr(domain, 1, 12))_PC<k>, and scripts 23 / 25 / 28 / 30 / 35 / 36 / 39 / 44 map
  # them back to the domain by that prefix. Two labels sharing a prefix would be merged silently (found 2026-09-16:
  # "Food prices (RTFP)" against "Food prices and supply"), so the collision is an error here.
  pref <- make.names(substr(domains, 1, 12))
  if (anyDuplicated(pref))
    stop("domain labels collide on their first 12 characters: ",
         paste(domains[pref %in% pref[duplicated(pref)]], collapse = " | "), call. = FALSE)
  if (is.null(sign_rows)) sign_rows <- seq_len(nrow(Xr))
  blocks <- list(); basis <- list()
  for (dm in domains) {
    cc <- cols[which(domain_of[cols] == dm)]
    if (!length(cc)) next
    M <- Xr[, cc, drop = FALSE]
    if (length(cc) == 1) {
      b <- M
      colnames(b) <- paste0(make.names(substr(dm, 1, 12)), "_PC1")
      blocks[[dm]] <- b
      basis[[dm]] <- list(cols = cc, center = 0, rot = matrix(1, 1, 1, dimnames = list(cc, colnames(b))))
      next
    }
    Mtr <- M[sign_rows, , drop = FALSE]
    if (nrow(Mtr) < 3) next
    pc <- tryCatch(stats::prcomp(Mtr, center = TRUE, scale. = FALSE),
                   error = function(e) NULL)
    if (is.null(pc)) next
    ve <- pc$sdev^2 / sum(pc$sdev^2)
    npc <- max(1L, which(cumsum(ve) >= var_target)[1])
    if (!is.finite(npc)) npc <- 1L
    npc <- min(npc, max_pc, ncol(pc$rotation), length(cc), nrow(Mtr) - 1L)
    if (npc < 1) next
    sc <- scale(M, center = pc$center, scale = FALSE) %*%
      pc$rotation[, seq_len(npc), drop = FALSE]
    # orient each component so the majority of its loadings are positive, for
    # cross-country comparability of the sign (the per-country orientation
    # hazard documented in PROTOCOL_V2.md)
    flip <- ifelse(colMeans(pc$rotation[, seq_len(npc), drop = FALSE] > 0) < 0.5,
                   -1, 1)
    sc <- sweep(sc, 2, flip, "*")
    colnames(sc) <- paste0(make.names(substr(dm, 1, 12)), "_PC", seq_len(npc))
    blocks[[dm]] <- sc
    # the linear map from centred predictors to oriented axes, kept so that an
    # axis-weighted index can be projected back onto its columns (WS-02)
    rot <- sweep(pc$rotation[, seq_len(npc), drop = FALSE], 2, flip, "*")
    dimnames(rot) <- list(cc, colnames(sc))
    basis[[dm]] <- list(cols = cc, center = pc$center, rot = rot)
  }
  if (!length(blocks)) return(matrix(numeric(0), nrow = nrow(Xr), ncol = 0))
  out <- do.call(cbind, blocks)
  out[!is.finite(out)] <- 0
  attr(out, "basis") <- basis
  out
}

#' The project's default domain representation.
#'
#' Switch back to the single-score version with
#' `Sys.setenv(V2_DOMAIN_REP = "mean1")` when reproducing pre-2026-09-01 runs.
domain_representation_v2 <- function(Xr, domain_of, sign_rows = NULL) {
  if (identical(Sys.getenv("V2_DOMAIN_REP", "pcvar"), "mean1"))
    build_domain_scores_v2(Xr, domain_of, sign_rows = sign_rows)
  else
    build_domain_pcs_v2(Xr, domain_of, sign_rows = sign_rows)
}

# ── FIX 1 & 5. Folds for the three estimands ────────────────────────────────

#' The three estimands, each with the baseline that saw the same information
#'
#' Fix 5. The audit found the project's leaderboard scored transport-style
#' models against an in-fill-style baseline (the flat regional mean, which reads
#' the held-out district's own survey responses and which the project's own
#' jackknife control had already withdrawn).
ESTIMANDS_V2 <- list(
  infill = list(
    label = "A. in-fill: unsurveyed district inside a surveyed region",
    scheme = "kfold_district",
    baselines = c("null_train_mean", "region_mean_jk", "spatial"),
    replicated = TRUE
  ),
  region = list(
    label = "B. extrapolation: whole unsurveyed region",
    scheme = "loro",
    baselines = c("null_train_mean", "spatial"),
    replicated = FALSE
  ),
  country = list(
    label = "C. transport: unsurveyed country",
    scheme = "loco",
    baselines = c("null_train_mean"),
    replicated = FALSE
  )
)

#' Build fold assignments for one cell
#'
#' Fix 1: `kfold_district` is REPLICATED. The audit traced the project's
#' headline within-country negative (median r 0.058) to a single unreplicated
#' 3-fold region draw sitting in the bottom decile of its own protocol's
#' distribution; re-randomising the identical protocol gave 0.217. Any scheme
#' with a random component must therefore be run R times and summarised over
#' draws, never reported from one draw.
#'
#' `loro` is exhaustive and deterministic, which removes draw luck entirely and
#' is the preferred within-country scheme for that reason.
make_folds_v2 <- function(scheme, n, blocks = NULL, k = 5, rep_id = 1) {
  if (scheme == "kfold_district") {
    set.seed(20260951L + rep_id)
    return(sample(rep(seq_len(min(k, n)), length.out = n)))
  }
  if (scheme %in% c("loro", "loco")) {
    if (is.null(blocks)) stop("blocks required for ", scheme)
    return(as.integer(factor(blocks)))
  }
  stop("unknown scheme: ", scheme)
}

# ── arms ────────────────────────────────────────────────────────────────────
# Every arm has the signature (tr, te, y, X, D, aux) and returns predictions on
# the modelling scale for the test rows. `aux` carries lon/lat, Admin1 and the
# n_eff weights so that arms which need them can see them, and so that an arm
# which peeks can be identified by reading one function.

arm_null_train_mean_v2 <- function(tr, te, y, X, D, aux) {
  rep(mean(y[tr]), length(te))
}

#' Region mean, jackknifed: the honest covariate-free comparator for in-fill
#'
#' Assigns each held-out district the mean of the OTHER surveyed districts in
#' its own region. This is the arm the project's rule 3.6 baseline should have
#' been: the withdrawn version used the held-out district's own respondents.
#' Falls back to the training mean for a region with no other training district.
arm_region_mean_jk_v2 <- function(tr, te, y, X, D, aux) {
  a1 <- aux$Admin1
  gm <- mean(y[tr])
  vapply(te, function(i) {
    same <- tr[a1[tr] == a1[i]]
    if (!length(same)) gm else mean(y[same])
  }, 0)
}

#' Fit the spatial smoother once on a training fold; predict any rows
#'
#' Returned as a closure so that an arm needing both training and test
#' predictions (spatial_plus_domain) fits one GAM rather than two.
.v2_spatial_fit <- function(tr, y, aux) {
  fallback <- mean(y[tr])
  if (!requireNamespace("mgcv", quietly = TRUE))
    return(function(idx) rep(fallback, length(idx)))
  dtr <- data.frame(y = y[tr], lon = aux$lon[tr], lat = aux$lat[tr])
  dtr <- dtr[stats::complete.cases(dtr), ]
  if (nrow(dtr) < 12 || dplyr::n_distinct(dtr$lon) < 5)
    return(function(idx) rep(fallback, length(idx)))
  kk <- max(3, min(25, floor(nrow(dtr) / 3)))
  fit <- tryCatch(mgcv::gam(y ~ s(lon, lat, k = kk), data = dtr),
                  error = function(e) NULL)
  if (is.null(fit)) return(function(idx) rep(fallback, length(idx)))
  function(idx) {
    p <- tryCatch(as.numeric(stats::predict(fit, data.frame(
      lon = aux$lon[idx], lat = aux$lat[idx]))), error = function(e) NULL)
    if (is.null(p)) return(rep(fallback, length(idx)))
    p[!is.finite(p)] <- fallback
    p
  }
}

arm_spatial_v2 <- function(tr, te, y, X, D, aux) {
  .v2_spatial_fit(tr, y, aux)(te)
}

#' @param standardize glmnet's internal column standardisation. FALSE is right
#'   when the caller has already put every column on a common scale (all arms
#'   in 02/02b do). It must be TRUE when comparing scaling schemes, or the
#'   unscaled arm is penalised for its units rather than judged on its
#'   information and the comparison is not fair.
.v2_enet <- function(Xtr, ytr, Xte, alpha = 0.5, standardize = FALSE) {
  if (!requireNamespace("glmnet", quietly = TRUE) || ncol(Xtr) < 2)
    return(rep(mean(ytr), nrow(Xte)))
  if (nrow(Xtr) < 12) return(rep(mean(ytr), nrow(Xte)))
  nf <- max(3, min(5, floor(nrow(Xtr) / 5)))
  # nlambda 40 rather than the default 100: at n <= 87 the extra path
  # resolution changes lambda.min negligibly and costs ~2.5x the runtime.
  fit <- tryCatch(glmnet::cv.glmnet(Xtr, ytr, alpha = alpha, nfolds = nf,
                                    nlambda = 40, standardize = standardize),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(mean(ytr), nrow(Xte)))
  p <- tryCatch(as.numeric(stats::predict(fit, Xte, s = "lambda.min")),
                error = function(e) NULL)
  if (is.null(p) || !all(is.finite(p))) return(rep(mean(ytr), nrow(Xte)))
  p
}

arm_domain_enet_v2 <- function(tr, te, y, X, D, aux) {
  .v2_enet(D[tr, , drop = FALSE], y[tr], D[te, , drop = FALSE],
           standardize = isTRUE(aux$enet_standardize))
}

arm_raw_enet_v2 <- function(tr, te, y, X, D, aux) {
  .v2_enet(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
}

#' Zero-tuning sign-aligned domain index (the probe P3a estimator)
#'
#' Weight each domain by its training-fold Fisher-z Spearman association with
#' the outcome, then take the weighted sum. No hyperparameters, nothing to
#' overfit, and the arm that transported at Spearman 0.309 across countries.
#'
#' The score is on an arbitrary scale. It is standardised on the training rows
#' and mapped to the outcome scale as mean(y_tr) + rho * sd(y_tr) * z. For the
#' RANKING product (`domain_index`) rho = 1: the map is a fixed affine rescale
#' that keeps correlation meaningful without any test information. For the
#' LEVEL product (`domain_index_cal`, IS-01, 2026-09-18) rho is the score's
#' out-of-sample correlation with the outcome, estimated honestly by nested
#' 5-fold predictions inside the training rows: a score that explains rho^2 of
#' the outcome's variance must be given rho times its spread, otherwise the
#' levels are over-dispersed by 1/rho (in-fill prevalence MAE 14.5 pp against
#' 12.4 for the training mean; calibrated 11.5, RR-12 cells). The two arms have
#' IDENTICAL rankings inside any one fit; they are kept apart because the
#' protocol pools out-of-fold predictions across folds before ranking them,
#' and a shrunk spread lets fold-to-fold differences in the training mean into
#' the pooled order (Spearman 0.39 -> 0.34 on the same rankings,
#' index_shrinkage_variants.csv). Read Spearman / top-k from `domain_index`
#' and MAE / bias from `domain_index_cal`. V2_INDEX_SHRINK overrides the
#' default of `domain_index` ("none" | "nested" | "insample"); in-sample rho
#' under-shrinks (MAE 0.189 / 11.7 vs nested 0.187 / 11.5) and is kept only
#' for comparison. rho is floored at 0.001 so a negative estimate collapses
#' the levels to the mean without destroying the ranking.
.index_weights_v2 <- function(Dtr, ytr) {
  n <- length(ytr)
  z <- apply(Dtr, 2, function(x) {
    if (stats::sd(x) == 0) return(0)
    r <- suppressWarnings(stats::cor(x, ytr, method = "spearman"))
    if (!is.finite(r)) return(0)
    r <- max(min(r, 0.999), -0.999)
    0.5 * log((1 + r) / (1 - r)) * sqrt(max(n - 3, 1))
  })
  z[!is.finite(z)] <- 0
  z
}
.index_rho_v2 <- function(tr, y, D, how = c("nested", "insample", "none"), k = 5L) {
  how <- match.arg(how)
  if (how == "none") return(1)
  ytr <- y[tr]
  if (how == "insample") {
    s <- as.numeric(D[tr, , drop = FALSE] %*% .index_weights_v2(D[tr, , drop = FALSE], ytr))
    r <- suppressWarnings(stats::cor(s, ytr))
  } else {
    # deterministic inner folds (systematic 1..k over the training rows): no
    # call to the RNG, so the caller's random stream is untouched
    f <- rep_len(seq_len(min(k, length(tr))), length(tr))
    oof <- rep(NA_real_, length(tr))
    for (j in unique(f)) {
      itr <- tr[f != j]; ite <- tr[f == j]
      if (length(itr) < 5) next
      w <- .index_weights_v2(D[itr, , drop = FALSE], y[itr])
      s_itr <- as.numeric(D[itr, , drop = FALSE] %*% w)
      s_ite <- as.numeric(D[ite, , drop = FALSE] %*% w)
      if (stats::sd(s_itr) == 0) { oof[f == j] <- mean(y[itr]); next }
      oof[f == j] <- ((s_ite - mean(s_itr)) / stats::sd(s_itr)) * stats::sd(y[itr]) + mean(y[itr])
    }
    r <- suppressWarnings(stats::cor(oof, ytr, use = "complete.obs"))
  }
  if (!is.finite(r)) r <- 0
  max(min(r, 1), 0.001)
}
arm_domain_index_v2 <- function(tr, te, y, X, D, aux) {
  ytr <- y[tr]
  z <- .index_weights_v2(D[tr, , drop = FALSE], ytr)
  idx_tr <- as.numeric(D[tr, , drop = FALSE] %*% z)
  idx_te <- as.numeric(D[te, , drop = FALSE] %*% z)
  if (stats::sd(idx_tr) == 0) return(rep(mean(ytr), length(te)))
  rho <- .index_rho_v2(tr, y, D, how = Sys.getenv("V2_INDEX_SHRINK", "none"))
  ((idx_te - mean(idx_tr)) / stats::sd(idx_tr)) * rho * stats::sd(ytr) + mean(ytr)
}
#' The calibrated-level index: the same ranking, spread shrunk by the nested
#' out-of-sample correlation (IS-01). The arm to read levels, MAE and bias from.
#'
#' On the prevalence target the shrinkage happens on the logit scale, and the
#' back-transformed logit mean sits below the arithmetic mean (-2.6 pp in
#' RR-13; the spatial arm shows the same). When `aux$target == "prev"` and
#' `aux$y_nat` is present the logit predictions are therefore shifted by the
#' constant that makes the mean back-transformed TRAINING prediction equal
#' the training rows' mean prevalence (the null arm's own definition), which
#' is exactly the dashboard's national anchor applied to the training mean.
#' A shift is a monotone map, so the ranking is untouched.
arm_domain_index_cal_v2 <- function(tr, te, y, X, D, aux) {
  old <- Sys.getenv("V2_INDEX_SHRINK", unset = NA)
  Sys.setenv(V2_INDEX_SHRINK = "nested"); on.exit(if (is.na(old)) Sys.unsetenv("V2_INDEX_SHRINK") else Sys.setenv(V2_INDEX_SHRINK = old))
  p_te <- arm_domain_index_v2(tr, te, y, X, D, aux)
  if (identical(aux$target, "prev") && !is.null(aux$y_nat) && all(is.finite(aux$y_nat[tr]))) {
    p_tr <- arm_domain_index_v2(tr, tr, y, X, D, aux)
    target_mean <- mean(aux$y_nat[tr])
    f <- function(c) mean(.v2_expit(p_tr + c)) - target_mean
    shift <- tryCatch(stats::uniroot(f, c(-8, 8))$root, error = function(e) 0)
    p_te <- p_te + shift
  }
  p_te
}

arm_spatial_plus_domain_v2 <- function(tr, te, y, X, D, aux) {
  predict_sp <- .v2_spatial_fit(tr, y, aux)     # one GAM, used for both sets
  r_tr <- y[tr] - predict_sp(tr)
  add <- .v2_enet(D[tr, , drop = FALSE], r_tr, D[te, , drop = FALSE])
  predict_sp(te) + add
}

ARMS_V2 <- list(
  null_train_mean     = arm_null_train_mean_v2,
  region_mean_jk      = arm_region_mean_jk_v2,
  spatial             = arm_spatial_v2,
  domain_index        = arm_domain_index_v2,
  domain_index_cal    = arm_domain_index_cal_v2,
  domain_enet         = arm_domain_enet_v2,
  raw_enet            = arm_raw_enet_v2,
  spatial_plus_domain = arm_spatial_plus_domain_v2
)

#' Which arms may be run under which estimand, and what each one saw
#'
#' `region_mean_jk` needs other surveyed districts inside the held-out unit's
#' own region, which exist only under in-fill. Running it under `region` or
#' `country` would be exactly the information asymmetry fix 5 exists to remove.
arms_for_estimand_v2 <- function(estimand) {
  base <- c("null_train_mean", "spatial", "domain_index", "domain_enet",
            "raw_enet", "spatial_plus_domain")
  # the calibrated-level index only where levels are scored (transport blanks MAE)
  if (estimand != "country") base <- c(base, "domain_index_cal")
  if (estimand == "infill") c(base, "region_mean_jk") else base
}

# ── scoring ─────────────────────────────────────────────────────────────────

#' Precision-weighted metrics
#'
#' Fix 1: the audit found the primary benchmark scores districts unweighted,
#' so a district resting on a single respondent counts as much as one resting
#' on 87, and no minimum-n rule exists anywhere. Weighted metrics use n_eff.
score_v2 <- function(obs, pred, w = NULL, scale = c("prev", "level")) {
  scale <- match.arg(scale)
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  w <- if (is.null(w)) rep(1, length(obs)) else w[ok]
  good <- w[is.finite(w) & w > 0]
  w[!is.finite(w) | w <= 0] <- if (length(good)) min(good) else 1
  n <- length(obs)
  if (n < 4 || stats::sd(obs) == 0) {
    return(data.frame(n = n, spearman = NA_real_, pearson = NA_real_,
                      mae = NA_real_, wmae = NA_real_, bias = NA_real_,
                      rmse_sd = NA_real_, topk = NA_real_))
  }
  sp <- if (stats::sd(pred) == 0) NA_real_ else
    suppressWarnings(stats::cor(obs, pred, method = "spearman"))
  pe <- if (stats::sd(pred) == 0) NA_real_ else
    suppressWarnings(stats::cor(obs, pred))
  mult <- if (scale == "prev") 100 else 1
  ae <- abs(obs - pred)
  # top-quartile capture: share of the truly worst quarter also ranked worst
  q <- max(2, ceiling(n / 4))
  worst_true <- order(obs, decreasing = TRUE)[seq_len(q)]
  worst_pred <- order(pred, decreasing = TRUE)[seq_len(q)]
  data.frame(
    n = n,
    spearman = sp, pearson = pe,
    mae  = mean(ae) * mult,
    wmae = stats::weighted.mean(ae, w) * mult,
    bias = mean(pred - obs) * mult,
    rmse_sd = sqrt(mean((obs - pred)^2)) / stats::sd(obs),
    topk = length(intersect(worst_true, worst_pred)) / q
  )
}

# ── modelled-surface handling ───────────────────────────────────────────────

#' Optionally drop predictors that are themselves modelled surfaces
#'
#' These columns are INCLUDED BY DEFAULT. An earlier version of this function
#' excluded them on the suspicion that they were trained on the same surveys
#' this project scores against. That suspicion was checked against source
#' documentation and does not hold:
#'
#'   GFDx anaemia  = WHO, "The global prevalence of anaemia in 2011" (2015).
#'                   It PRE-DATES every survey here (Gambia 2021, Ghana 2017,
#'                   Malawi 2015-16, Sierra Leone 2013), so it cannot have been
#'                   fitted to them.
#'   GFDx zinc     = Wessells & Brown (2012), estimated from FAO food-balance
#'                   -sheet zinc availability and stunting. It uses no
#'                   biomarker survey at all, which makes it a food-supply
#'                   proxy - the mechanistically desirable kind of predictor.
#'   IHME stunting/wasting/underweight = anthropometry. A different construct
#'                   from every micronutrient biomarker modelled here.
#'   IHME anaemia  = the LBD "global anemia prev geospatial estimates
#'                   2000-2019" surface for women 15-49, defined on
#'                   HAEMOGLOBIN. This project's outcomes are ferritin-based
#'                   iron deficiency, RBP-based vitamin A deficiency, folate,
#'                   B12 and zinc. Haemoglobin is a different biomarker and
#'                   anaemia a different condition - the project's own audit
#'                   makes exactly this point when criticising the WHO anaemia
#'                   bands applied to ferritin-based iron deficiency.
#'
#' The residual caveat is narrow and worth stating rather than acting on: IHME
#' anaemia is fitted partly to DHS haemoglobin from survey rounds that in some
#' countries ran alongside the micronutrient surveys used here, and haemoglobin
#' correlates with iron status. So for the IRON outcomes specifically there is a
#' weak shared-information channel; for vitamin A, folate, B12 and zinc there is
#' effectively none. That is a sensitivity analysis, not grounds for exclusion -
#' and a modelled haemoglobin surface is in fact an attractive predictor, since
#' the project has measured a +0.174 gain from FIELD haemoglobin, which needs a
#' blood draw, whereas this surface needs no survey at all.
#'
#' Set V2_DROP_MODELLED=1 to exclude them for a sensitivity run.
#' Exclusion rules (metadata/covariates/exclusions.csv), found from the project root
.v2_exclusions <- function() {
  d <- getwd()
  for (i in 1:6) {
    f <- file.path(d, "metadata", "covariates", "exclusions.csv")
    if (file.exists(f)) { ex <- read.csv(f, stringsAsFactors = FALSE); if (!"policy" %in% names(ex)) ex$policy <- "data_defect"; return(ex) }
    d <- dirname(d)
  }
  NULL
}

#' Find a project file upward from the working directory (scripts setwd() to
#' the project root, tests run from tests/testthat).
.v2_project_file <- function(...) {
  d <- getwd()
  for (i in 1:6) {
    f <- file.path(d, ...)
    if (file.exists(f)) return(f)
    d <- dirname(d)
  }
  NULL
}

# ── Predictor tiers (TP-01, 2026-09-15) ─────────────────────────────────────
#
# Every column of the shared set belongs to one of three tiers, assigned from
# its `source` by the first matching row of metadata/covariates/predictor_tiers.csv:
#   open           gridded / modelled / administrative products obtainable for
#                  any country without fielding a household survey
#   survey_public  derived from public household surveys other than DHS
#                  (HCES microdata, MICS regional estimates via WHO HEAT)
#   survey_dhs     derived from DHS microdata
# V2_PREDICTOR_TIERS (comma-separated) selects the tiers a run may use; the
# default is every tier. The transport driver (02b) defaults to
# "open,survey_public" - the pre-registered no-DHS arm, which DA-04 showed
# transports better - so the headline transport figure never depends on DHS
# unless a run asks for it.
V2_TIERS_ALL <- c("open", "survey_public", "survey_dhs")

tier_rules_v2 <- function() {
  f <- .v2_project_file("metadata", "covariates", "predictor_tiers.csv")
  if (is.null(f)) stop("metadata/covariates/predictor_tiers.csv not found above ", getwd())
  r <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(c("source_regex", "tier") %in% names(r)), all(r$tier %in% V2_TIERS_ALL))
  r
}

#' Tier of every row of a metadata table, from its `source`
assign_tier_v2 <- function(meta) {
  rules <- tier_rules_v2()
  src <- as.character(meta$source); src[is.na(src)] <- ""
  tier <- rep(NA_character_, length(src))
  for (i in seq_len(nrow(rules))) {
    hit <- is.na(tier) & grepl(rules$source_regex[i], src, perl = TRUE)
    tier[hit] <- rules$tier[i]
  }
  if (anyNA(tier)) stop("predictor_tiers.csv has no rule for source(s): ", paste(unique(src[is.na(tier)]), collapse = "; "))
  tier
}

#' The tiers this run may use (V2_PREDICTOR_TIERS, default all)
predictor_tiers_v2 <- function() {
  s <- Sys.getenv("V2_PREDICTOR_TIERS", "")
  if (!nzchar(s)) return(V2_TIERS_ALL)
  t <- trimws(strsplit(s, ",")[[1]])
  bad <- setdiff(t, V2_TIERS_ALL)
  if (length(bad)) stop("V2_PREDICTOR_TIERS: unknown tier(s) ", paste(bad, collapse = ", "), "; valid: ", paste(V2_TIERS_ALL, collapse = ", "))
  t
}

.v2_as_logical <- function(x) { x <- toupper(trimws(as.character(x))); out <- x %in% c("TRUE", "T", "1"); out[is.na(x) | x == "" | x == "NA"] <- NA; out }

#' The fit-time predictor policy. Applied by every script to the column list
#' it takes from the shared set, so a column that reaches the file by any route
#' is still governed here. In order:
#'   1. leakage rules in exclusions.csv (policy == "leakage")           always
#'   2. national constants (metadata subnational == FALSE)               unless V2_KEEP_NATIONAL=1
#'      A value broadcast to every district of a country cannot rank its
#'      districts; in the four-country transport it is a three-value country
#'      effect. (The per-country zero-variance filter in prep_predictors_v2()
#'      already removed them silently; this makes the policy explicit.)
#'   3. predictor tiers (V2_PREDICTOR_TIERS)                             default: all tiers
#'   4. modelled surfaces                                                only with V2_DROP_MODELLED=1
drop_near_outcome_v2 <- function(preds, meta) {
  n0 <- length(preds)
  # Leakage policy (LK-01, 2026-09-07): the rows of exclusions.csv flagged policy == "leakage"
  # are enforced at fit time as well as in the builder, so a column that reaches the shared
  # set by any route is still kept out of the design matrix. Always on.
  ex <- .v2_exclusions()
  if (!is.null(ex) && any(ex$policy == "leakage") && length(preds)) {
    rx <- ex$canonical_regex[ex$policy == "leakage"]
    bad <- preds[Reduce(`|`, lapply(rx, function(r) grepl(r, preds, perl = TRUE)), init = rep(FALSE, length(preds)))]
    if (length(bad)) {
      message(sprintf("[protocol v2] leakage policy: excluded %d predictor(s): %s", length(bad), paste(bad, collapse = ", ")))
      preds <- setdiff(preds, bad)
    }
  }
  # National-constant policy (TP-01)
  n_nat <- 0L
  if (!identical(Sys.getenv("V2_KEEP_NATIONAL", "0"), "1") && "subnational" %in% names(meta)) {
    sub <- .v2_as_logical(meta$subnational)
    nat <- meta$column[!is.na(sub) & !sub]
    bad <- intersect(preds, nat); n_nat <- length(bad)
    preds <- setdiff(preds, bad)
  }
  # Tier policy (TP-01)
  tiers <- predictor_tiers_v2(); n_tier <- 0L
  if (!setequal(tiers, V2_TIERS_ALL) && length(preds)) {
    tier <- if ("tier" %in% names(meta)) as.character(meta$tier) else assign_tier_v2(meta)
    tier_of <- stats::setNames(tier, meta$column)
    bad <- preds[!(tier_of[preds] %in% tiers)]; n_tier <- length(bad)
    preds <- setdiff(preds, bad)
  }
  message(sprintf("[protocol v2] predictor policy: %d of %d kept | tiers %s | %d national constant(s) %s | %d outside the tiers",
                  length(preds), n0, paste(tiers, collapse = ","), n_nat,
                  if (identical(Sys.getenv("V2_KEEP_NATIONAL", "0"), "1")) "kept (V2_KEEP_NATIONAL=1)" else "dropped", n_tier))
  if (!identical(Sys.getenv("V2_DROP_MODELLED", "0"), "1")) return(preds)
  bad <- if ("modelled_surface" %in% names(meta)) meta$column[.v2_as_logical(meta$modelled_surface) %in% TRUE]
         else if ("domain" %in% names(meta)) meta$column[grepl("MODELLED SURFACE", meta$domain, fixed = TRUE)] else character(0)
  keep <- setdiff(preds, bad)
  n <- length(preds) - length(keep)
  if (n) message(sprintf("[protocol v2] sensitivity: excluded %d modelled-surface predictors", n))
  keep
}

# ── JK-01: checked Admin-2 joins travel with the protocol ────────────────────
# Scripts that source only this file get join_admin2_v2(), admin2_population_v2()
# and admin2_spine() from R/admin2_keys.R.
if (!exists("join_admin2_v2", mode = "function")) {
  .a2 <- .v2_project_file("R", "admin2_keys.R")
  if (!is.null(.a2)) source(.a2)
}
