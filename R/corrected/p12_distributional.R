# =============================================================================
# R/corrected/p12_distributional.R  —  P12: the distributional estimator, for
# every outcome with a continuous biomarker
#
# THE IDEA
# --------
# Deficiency is a threshold on a continuous biomarker. The production path
# dichotomises first and then models the binary. The distributional estimator
# models the continuous biomarker and integrates the fitted distribution past
# the cut-point instead, which keeps the information in how far a value sits
# from the threshold.
#
# The prototype (sensitivity/10_distributional_vs_binary_prototype.R, one
# country, three cells, documented in docs/exploratory_distributional_prototype.md)
# found the gain concentrates where the cut-point sits off the median: for
# Ghana women's iron at 12.3 percent prevalence it roughly tripled the area-level
# correlation, and for mid-prevalence cells it was a wash. That is the expected
# shape, since the efficiency of continuous over binary grows as the threshold
# moves into a tail.
#
# WHAT THIS FILE ADDS
# -------------------
# The prototype ran one country and three cells with no calibration or
# discrimination metrics. This runs every country and every outcome that has a
# continuous biomarker, and scores both arms on the four families the corrected
# methods layer reports: discrimination, calibration, Brier skill, and Admin-2
# correlation read against the reliability ceiling.
#
# Both arms produce an individual-level predicted probability, which is what
# makes AUC and Brier comparable between them:
#
#   binary_classifier      ridge logistic on the individual binary outcome
#   distributional_default ridge gaussian on the individual log biomarker, then
#                          p_i = P(mu_i + e < log(threshold)) taken over the
#                          empirical TRAINING residual distribution e
#
# The residual distribution is empirical rather than Gaussian on purpose: the
# heteroskedasticity check in sensitivity/11_distributional_heterosked_transport.R
# is the reason not to assume a constant-variance normal.
#
# HARMONIZED INPUTS
# -----------------
# Where WS3 provides a harmonized biomarker it is used, so this workstream is
# not silently re-measuring adjustment heterogeneity:
#   vitamin A  BRINDA CRP+AGP adjusted RBP
#   iron       BRINDA CRP+AGP adjusted ferritin
# Folate, B12 and zinc have no inflammation adjustment defined in this pipeline,
# so their configured continuous column is used as-is and the row records that.
#
# CROSS-VALIDATION
# ----------------
# Leave-one-Admin-2-out within country, matching the prototype. Every fold fits
# its own imputation, standardization and glmnet on the training areas only, and
# the held-out area contributes to no training quantity.
# =============================================================================

#' Covariate set. Deliberately the prototype's fixed a-priori eight, so the
#' comparison against docs/exploratory_distributional_prototype.md is
#' like-for-like and the two arms differ only in how the outcome is modelled.
#' Widening it would confound "distributional beats binary" with "more
#' covariates beat fewer".
DIST_COVARIATE_PATTERNS <- c(
  "gee_soilzinc_mean_0_20", "gee_soiliron_mean_0_20",
  "gee_soilnitrogen_mean_0_20", "gee_soilphosphorus_mean_0_20",
  "gee_accessibility_2019", "gee_elevation_2000",
  "gee_popdensity_2010", "gee_ghsbuilts_2015_built_surface"
)

DIST_SCHEMES <- c("binary_classifier", "distributional_default")

#' Areas with fewer than this many biomarker reads are dropped from the
#' EVALUATION only, never from training. Matches the prototype's MIN_NSVY.
DIST_MIN_NSVY <- 8L

#' In-fold median imputation and standardization, training statistics only.
.dist_prep <- function(Xtr, Xte) {
  med <- apply(Xtr, 2, stats::median, na.rm = TRUE); med[!is.finite(med)] <- 0
  ctr <- colMeans(Xtr, na.rm = TRUE);                ctr[!is.finite(ctr)] <- 0
  sdv <- apply(Xtr, 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  imp <- function(M) {
    for (j in seq_len(ncol(M))) {
      x <- M[, j]; x[!is.finite(x)] <- med[j]; M[, j] <- (x - ctr[j]) / sdv[j]
    }
    M
  }
  list(tr = imp(Xtr), te = imp(Xte))
}

.dist_ridge <- function(x, y, family, w = NULL, seed = 12345L) {
  if (is.null(w)) w <- rep(1, length(y))
  if (!is.null(seed)) set.seed(seed)
  tryCatch(glmnet::cv.glmnet(x, y, family = family, alpha = 0, weights = w, nfolds = 5),
           error = function(e) NULL)
}

#' Resolve the harmonized continuous biomarker for one country-outcome.
#'
#' @return list(value, threshold, marker, harmonized) or NULL when the outcome
#'   has no usable continuous biomarker. Never falls back silently: the
#'   `harmonized` flag records which path was taken.
dist_continuous_biomarker <- function(merged, cc, oc) {
  pop <- if (grepl("^child", oc$tag)) "child" else "women"

  if (grepl("vitA", oc$tag, ignore.case = TRUE)) {
    spec <- brinda_rbp_cols(cc$country)[[pop]]
    if (!is.null(spec) && all(unlist(spec[c("rbp", "crp")]) %in% colnames(merged)))
      return(list(value = brinda_adjust_rbp(merged[[spec$rbp]], merged[[spec$crp]],
                                            if (is.null(spec$agp)) NULL else merged[[spec$agp]]),
                  threshold = 0.70, marker = "BRINDA-adjusted RBP (umol/L)",
                  harmonized = TRUE))
  }
  if (grepl("iron", oc$tag, ignore.case = TRUE)) {
    spec <- brinda_ferritin_cols(cc$country)[[pop]]
    if (!is.null(spec) && all(unlist(spec[c("fer", "crp")]) %in% colnames(merged)))
      return(list(value = brinda_adjust_ferritin(merged[[spec$fer]], merged[[spec$crp]],
                                                 if (is.null(spec$agp)) NULL else merged[[spec$agp]]),
                  threshold = if (pop == "child") 12 else 15,
                  marker = "BRINDA-adjusted ferritin (ug/L)", harmonized = TRUE))
  }

  # Folate, B12, zinc: no inflammation adjustment is defined in this pipeline,
  # so the configured continuous column is used as-is.
  cont <- oc$continuous
  if (is.null(cont) || !cont %in% colnames(merged) || is.null(oc$cutoff)) return(NULL)
  thr <- oc$cutoff
  if (identical(oc$cutoff_scale, "log")) thr <- exp(thr)
  list(value = suppressWarnings(as.numeric(merged[[cont]])), threshold = thr,
       marker = paste0("configured column ", cont, ", no inflammation adjustment"),
       harmonized = FALSE)
}

#' Metrics shared by both arms.
#'
#' Individual level: AUC and Brier, with Brier skill taken against the TRAINING
#' prevalence, which is the honest null a transported or held-out prediction has
#' to beat. Area level: Pearson and Spearman against the design-based Admin-2
#' prevalence, MAE, and the calibration slope of truth on prediction.
.dist_metrics <- function(p_ind, y_ind, train_prev, area_pred, area_truth) {
  ok <- is.finite(p_ind) & is.finite(y_ind)
  auc <- if (sum(ok) > 10 && length(unique(y_ind[ok])) > 1)
    tryCatch(as.numeric(pROC::auc(pROC::roc(y_ind[ok], p_ind[ok], quiet = TRUE))),
             error = function(e) NA_real_) else NA_real_
  brier      <- mean((y_ind[ok] - p_ind[ok])^2)
  brier_null <- mean((y_ind[ok] - train_prev)^2)
  ok2 <- is.finite(area_pred) & is.finite(area_truth)
  cal <- if (sum(ok2) > 3)
    tryCatch(unname(stats::coef(stats::lm(area_truth[ok2] ~ area_pred[ok2]))[2]),
             error = function(e) NA_real_) else NA_real_
  data.frame(
    n_ind        = sum(ok),
    n_areas_eval = sum(ok2),
    auc          = auc,
    brier        = brier,
    brier_null_train = brier_null,
    brier_skill  = if (is.finite(brier_null) && brier_null > 0) 1 - brier / brier_null else NA_real_,
    pearson_r    = if (sum(ok2) > 3) suppressWarnings(stats::cor(area_pred[ok2], area_truth[ok2])) else NA_real_,
    spearman_r   = if (sum(ok2) > 3) suppressWarnings(stats::cor(area_pred[ok2], area_truth[ok2], method = "spearman")) else NA_real_,
    mae_pp       = if (sum(ok2) > 3) 100 * mean(abs(area_pred[ok2] - area_truth[ok2])) else NA_real_,
    calib_slope  = cal,
    stringsAsFactors = FALSE
  )
}

#' Run both arms for one country-outcome under leave-one-Admin-2-out CV.
#'
#' @param merged per-country merged data.frame
#' @param area_cov data.frame with Admin2 plus the covariate columns
#' @param svy compute_svy_admin2() output for this cell
#' @param cc,oc country and outcome configs
#' @param seed integer
#' @return data.frame with one row per scheme, or NULL
run_distributional_cell <- function(merged, area_cov, svy, cc, oc, seed = 12345L) {
  bio <- dist_continuous_biomarker(merged, cc, oc)
  if (is.null(bio)) return(NULL)

  d <- merged
  pop_col <- cc$child_flag
  val <- bio$value
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    keep <- !is.na(d[[pop_col]]) & d[[pop_col]] == oc$child_flag_val
    d <- d[keep, , drop = FALSE]; val <- val[keep]
  }

  w <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  ind <- data.frame(Admin2 = as.character(d[[cc$admin2_col]]),
                    cont = suppressWarnings(as.numeric(val)), wt = w,
                    stringsAsFactors = FALSE)
  ind <- ind[is.finite(ind$cont) & ind$cont > 0 &
             !is.na(ind$Admin2) & nzchar(ind$Admin2), , drop = FALSE]
  if (nrow(ind) < 100) return(NULL)
  ind$def  <- as.integer(ind$cont < bio$threshold)
  ind$logc <- log(ind$cont)
  if (length(unique(ind$def)) < 2) return(NULL)

  cov_names <- intersect(DIST_COVARIATE_PATTERNS, names(area_cov))
  if (length(cov_names) < 3) return(NULL)
  acov <- area_cov[, c("Admin2", cov_names), drop = FALSE]
  acov$Admin2 <- as.character(acov$Admin2)
  for (k in cov_names) acov[[k]] <- suppressWarnings(as.numeric(acov[[k]]))

  sv <- as.data.frame(svy); sv$Admin2 <- as.character(sv$Admin2)
  areas <- sort(unique(ind$Admin2))
  areas <- areas[areas %in% acov$Admin2 & areas %in% sv$Admin2]
  if (length(areas) < 8) return(NULL)

  A   <- acov[match(areas, acov$Admin2), , drop = FALSE]
  sv  <- sv[match(areas, sv$Admin2), , drop = FALSE]
  ind <- ind[ind$Admin2 %in% areas, , drop = FALSE]
  X   <- as.matrix(A[, cov_names, drop = FALSE])

  # Individual-level out-of-fold predictions, one per scheme.
  p_bin <- rep(NA_real_, nrow(ind))
  p_dis <- rep(NA_real_, nrow(ind))
  train_prev <- rep(NA_real_, nrow(ind))
  thr_log <- log(bio$threshold)

  for (i in seq_along(areas)) {
    te_rows <- which(ind$Admin2 == areas[i])
    tr_rows <- which(ind$Admin2 != areas[i])
    if (!length(te_rows) || length(tr_rows) < 50) next
    itr <- ind[tr_rows, , drop = FALSE]; ite <- ind[te_rows, , drop = FALSE]
    if (length(unique(itr$def)) < 2) next

    P <- .dist_prep(X[match(itr$Admin2, areas), , drop = FALSE],
                    X[match(ite$Admin2, areas), , drop = FALSE])

    fb <- .dist_ridge(P$tr, itr$def, "binomial", itr$wt, seed)
    if (!is.null(fb))
      p_bin[te_rows] <- as.numeric(stats::predict(fb, P$te, s = "lambda.min",
                                                  type = "response"))

    fd <- .dist_ridge(P$tr, itr$logc, "gaussian", itr$wt, seed)
    if (!is.null(fd)) {
      mu_tr <- as.numeric(stats::predict(fd, P$tr, s = "lambda.min"))
      mu_te <- as.numeric(stats::predict(fd, P$te, s = "lambda.min"))
      resid <- itr$logc - mu_tr
      p_dis[te_rows] <- vapply(mu_te, function(m)
        stats::weighted.mean((m + resid) < thr_log, itr$wt), numeric(1))
    }
    train_prev[te_rows] <- stats::weighted.mean(itr$def, itr$wt)
  }

  # Area aggregation happens AFTER the individual fit, which is the ordering the
  # prototype identified as the dominant lever.
  .area <- function(p) vapply(areas, function(a) {
    k <- ind$Admin2 == a & is.finite(p)
    if (!any(k)) NA_real_ else stats::weighted.mean(p[k], ind$wt[k])
  }, numeric(1))

  truth <- sv$svy_prev
  n_svy <- sv$n_svy
  evalable <- is.finite(truth) & is.finite(n_svy) & n_svy >= DIST_MIN_NSVY

  out <- list()
  for (scheme in DIST_SCHEMES) {
    p <- if (scheme == "binary_classifier") p_bin else p_dis
    ap <- .area(p)
    m <- .dist_metrics(p, ind$def, stats::weighted.mean(ind$def, ind$wt),
                       ap[evalable], truth[evalable])
    m$scheme            <- scheme
    m$country           <- cc$country
    m$outcome           <- oc$tag
    m$marker            <- bio$marker
    m$harmonized_marker <- bio$harmonized
    m$threshold         <- bio$threshold
    m$national_prev     <- stats::weighted.mean(ind$def, ind$wt)
    m$n_areas_total     <- length(areas)
    m$seed              <- seed
    # r_share is the correlation read against what sampling noise permits.
    m <- add_reliability_columns(m, sv[evalable, , drop = FALSE])
    out[[scheme]] <- m
  }
  dplyr::bind_rows(out)
}
