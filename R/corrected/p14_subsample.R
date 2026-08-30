# =============================================================================
# R/corrected/p14_subsample.R  —  P14: the survey-subsample cost-of-accuracy
# simulation, rebuilt
#
# WHY THIS IS A REBUILD
# ---------------------
# archive/docs/PROJECT_STATUS_2026-08_UPDATE.md section 5 reports this project's
# flagship result: simulate a reduced survey budget, and ask whether a
# covariate model can recover the subnational accuracy the cut clusters would
# have provided. Its stated conclusions are that for any district retaining even
# a few clusters the direct estimate beats the model by a wide margin, and that
# the model's only reliable value is for districts with ZERO retained clusters,
# where it beats the national mean by about 0.7 percentage points.
#
# No code or output for that simulation exists anywhere in this repository. The
# Stage 0 audit searched for *subsample*, *cost_of_accuracy*, *retention* and
# *sentinel* and found nothing, and results/ contains no retention table. The
# result lives only as prose. Under this branch's anti-fabrication rule its
# figures could not be cited, so the simulation is rebuilt here from the
# specification in that section before anything is added to it.
#
# SPECIFICATION, as stated in that section
# ----------------------------------------
#   * retention fractions 50 to 90 percent of clusters
#   * 25 replicates per fraction per country-outcome
#   * compare a covariate-based area model against "use the local direct
#     estimate where you have one, else the national average"
#   * benchmark both against each country's own FULL survey
#   * genuine k-fold cross-validation, so a district's model prediction is
#     never trained on that district's own subsampled estimate
#
# The last point is the correction the note describes: an initial pass had
# in-sample leakage and looked artificially good.
#
# WHAT IS ADDED
# -------------
#   spatial_cv   the WS5 cluster-level smoother supplying predictions for
#                zero-cluster districts, so skipped districts borrow from
#                surveyed neighbours rather than from covariates
#
# The reported quantity is RMSE against the full-survey district prevalence,
# split by whether a district retained any clusters. That split is the whole
# point: the two strata answer different questions and averaging them hides the
# result.
# =============================================================================

SUBSAMPLE_RETENTION <- c(0.5, 0.6, 0.7, 0.8, 0.9)
SUBSAMPLE_REPLICATES <- 25L
SUBSAMPLE_STRATEGIES <- c("direct_or_national", "model_cv", "spatial_cv")

#' K-fold assignment over districts, so a district's prediction is never trained
#' on its own subsampled estimate.
.ss_folds <- function(n, k = 5L) {
  k <- max(2L, min(k, n))
  sample(rep_len(seq_len(k), n))
}

#' Cross-validated covariate model over districts with a retained estimate.
#'
#' Ridge on the logit of the district's subsampled prevalence. Fitted only on
#' districts that HAVE a retained estimate, then used to predict every district,
#' including the zero-cluster ones. Predictions for training districts are
#' out-of-fold; predictions for zero-cluster districts come from the full fit,
#' since they contribute nothing to training and no leakage is possible.
#'
#' @param X covariate matrix, all districts, rows aligned to `have`
#' @param y subsampled prevalence, NA where a district retained no clusters
#' @param w weights (retained n)
#' @return numeric prediction for every district
.ss_model_cv <- function(X, y, w, k = 5L) {
  have <- which(is.finite(y))
  if (length(have) < 10) return(rep(NA_real_, length(y)))
  pred <- rep(NA_real_, length(y))
  yl <- stats::qlogis(pmin(pmax(y[have], 0.01), 0.99))
  Xh <- X[have, , drop = FALSE]
  wh <- w[have]; wh[!is.finite(wh) | wh <= 0] <- 1

  folds <- .ss_folds(length(have), k)
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    if (sum(tr) < 5) next
    m <- tryCatch(glmnet::cv.glmnet(Xh[tr, , drop = FALSE], yl[tr], alpha = 0,
                                    weights = wh[tr], nfolds = 5),
                  error = function(e) NULL)
    if (is.null(m)) next
    pred[have[te]] <- stats::plogis(
      as.numeric(stats::predict(m, Xh[te, , drop = FALSE], s = "lambda.min")))
  }
  # Zero-cluster districts: fit on everything that has an estimate. They are not
  # in the training set, so this is out-of-sample for them by construction.
  miss <- which(!is.finite(y))
  if (length(miss)) {
    m <- tryCatch(glmnet::cv.glmnet(Xh, yl, alpha = 0, weights = wh, nfolds = 5),
                  error = function(e) NULL)
    if (!is.null(m))
      pred[miss] <- stats::plogis(
        as.numeric(stats::predict(m, X[miss, , drop = FALSE], s = "lambda.min")))
  }
  pred
}

#' Cross-validated cluster-level spatial smoother, same contract as .ss_model_cv.
#'
#' Fitted on the RETAINED clusters, so a skipped district's prediction comes
#' from its surveyed neighbours. District-level folds again, so a retained
#' district's prediction is not trained on its own clusters.
#'
#' @param cl retained clusters (lon, lat, p_w, n, Admin2)
#' @param districts district names in output order
#' @param centroid data.frame(Admin2, lon, lat) to predict at
#' @return numeric prediction for every district
.ss_spatial_cv <- function(cl, districts, centroid, k = 5L) {
  pred <- rep(NA_real_, length(districts))
  if (nrow(cl) < 20) return(pred)
  d_have <- unique(cl$Admin2)
  folds <- .ss_folds(length(d_have), k)
  names(folds) <- d_have

  fit_predict <- function(train_cl, newd) {
    if (nrow(train_cl) < 20 || !nrow(newd)) return(rep(NA_real_, nrow(newd)))
    kk <- max(5L, min(30L, floor(nrow(train_cl) / 3)))
    f <- stats::as.formula(sprintf("p_w ~ s(lon, lat, bs = 'tp', k = %d)", kk))
    m <- tryCatch(mgcv::gam(f, data = train_cl, family = stats::binomial(),
                            weights = train_cl$n, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) return(rep(NA_real_, nrow(newd)))
    tryCatch(as.numeric(mgcv::predict.gam(m, newdata = newd, type = "response")),
             error = function(e) rep(NA_real_, nrow(newd)))
  }

  for (f in unique(folds)) {
    out_d <- names(folds)[folds == f]
    tr <- cl[!(cl$Admin2 %in% out_d), , drop = FALSE]
    idx <- match(out_d, districts)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) next
    pred[idx] <- fit_predict(tr, centroid[idx, , drop = FALSE])
  }
  # Zero-cluster districts: fit on every retained cluster.
  miss <- which(!(districts %in% d_have))
  if (length(miss)) pred[miss] <- fit_predict(cl, centroid[miss, , drop = FALSE])
  pred
}

#' One replicate of the retention simulation.
#'
#' @param cl_all all clusters for this country-outcome
#' @param truth full-survey district prevalence, aligned to `districts`
#' @param X district covariate matrix
#' @param centroid district centroids for the spatial arm
#' @param districts district names
#' @param frac retention fraction
#' @return data.frame, one row per strategy per coverage stratum
.ss_one_replicate <- function(cl_all, truth, X, centroid, districts, frac) {
  keep <- sample(seq_len(nrow(cl_all)), size = max(2L, round(frac * nrow(cl_all))))
  cl <- cl_all[keep, , drop = FALSE]

  # Subsampled district estimate from retained clusters only.
  agg <- stats::aggregate(cbind(num = cl$p_w * cl$n, den = cl$n),
                          by = list(Admin2 = cl$Admin2), FUN = sum)
  y <- rep(NA_real_, length(districts))
  w <- rep(0, length(districts))
  i <- match(agg$Admin2, districts)
  ok <- !is.na(i)
  y[i[ok]] <- agg$num[ok] / agg$den[ok]
  w[i[ok]] <- agg$den[ok]

  national <- sum(cl$p_w * cl$n) / sum(cl$n)
  has_clusters <- is.finite(y)

  preds <- list(
    # The baseline the note describes: local direct estimate where one exists,
    # national average elsewhere.
    direct_or_national = ifelse(has_clusters, y, national),
    model_cv   = .ss_model_cv(X, y, w),
    spatial_cv = .ss_spatial_cv(cl, districts, centroid)
  )

  rows <- list()
  for (s in names(preds)) {
    p <- preds[[s]]
    for (stratum in c("has_clusters", "zero_clusters")) {
      sel <- if (stratum == "has_clusters") has_clusters else !has_clusters
      sel <- sel & is.finite(p) & is.finite(truth)
      if (sum(sel) < 1) next
      rows[[length(rows) + 1L]] <- data.frame(
        strategy = s, stratum = stratum, retention = frac,
        n_districts = sum(sel),
        rmse_pp = 100 * sqrt(mean((p[sel] - truth[sel])^2)),
        mae_pp  = 100 * mean(abs(p[sel] - truth[sel])),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

#' Run the retention simulation for one country-outcome.
#'
#' @param cl_all cluster table from build_cluster_dataset(), with covariates joined
#' @param svy full-survey Admin-2 table
#' @param covs covariate column names
#' @param fracs retention fractions
#' @param reps replicates per fraction
#' @param seed integer
#' @return data.frame of per-replicate rows
run_subsample_simulation <- function(cl_all, svy, covs,
                                     fracs = SUBSAMPLE_RETENTION,
                                     reps = SUBSAMPLE_REPLICATES,
                                     seed = 12345L) {
  sv <- as.data.frame(svy); sv$Admin2 <- as.character(sv$Admin2)
  districts <- sort(unique(cl_all$Admin2))
  districts <- districts[districts %in% sv$Admin2]
  if (length(districts) < 10) return(NULL)
  cl_all <- cl_all[cl_all$Admin2 %in% districts, , drop = FALSE]

  truth <- sv$svy_prev[match(districts, sv$Admin2)]
  n_svy <- sv$n_svy[match(districts, sv$Admin2)]
  # Districts whose own full-survey estimate is itself too thin to be a target
  # are dropped from the evaluation, matching every other workstream here.
  usable <- is.finite(truth) & is.finite(n_svy) & n_svy >= DIST_MIN_NSVY
  truth[!usable] <- NA_real_

  cov_ok <- intersect(covs, names(cl_all))
  if (length(cov_ok) < 3) return(NULL)
  dcov <- cl_all[!duplicated(cl_all$Admin2), c("Admin2", cov_ok), drop = FALSE]
  X <- as.matrix(dcov[match(districts, dcov$Admin2), cov_ok, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE); if (!is.finite(m)) m <- 0
    X[!is.finite(X[, j]), j] <- m
    s <- stats::sd(X[, j]); if (!is.finite(s) || s == 0) s <- 1
    X[, j] <- (X[, j] - mean(X[, j])) / s
  }

  centroid <- stats::aggregate(cbind(lon = cl_all$lon, lat = cl_all$lat),
                               by = list(Admin2 = cl_all$Admin2), FUN = mean)
  centroid <- centroid[match(districts, centroid$Admin2), , drop = FALSE]

  set.seed(seed)
  out <- list()
  for (f in fracs) for (r in seq_len(reps)) {
    one <- .ss_one_replicate(cl_all, truth, X, centroid, districts, f)
    if (is.null(one)) next
    one$replicate <- r
    out[[length(out) + 1L]] <- one
  }
  if (!length(out)) return(NULL)
  res <- dplyr::bind_rows(out)
  res$n_clusters_total <- nrow(cl_all)
  res$n_districts_total <- length(districts)
  res$seed <- seed
  res
}

# ── Calibration learning curve ───────────────────────────────────────────────

#' How many sentinel districts does a transported map need?
#'
#' Simulates deploying a model trained on OTHER countries into a target country
#' that can afford k districts of its own. Calibration is a location shift on the
#' logit scale fitted on those k districts only, which is the same monotone
#' correction WS7 uses, and error is measured on the districts NOT used for
#' calibration.
#'
#' @param pred transported prediction per district
#' @param truth full-survey prevalence per district
#' @param ks sentinel counts to trace
#' @param reps replicates per k
#' @param seed integer
#' @return data.frame, one row per k per replicate
calibration_learning_curve <- function(pred, truth, ks = c(0, 3, 5, 10, 15, 20),
                                       reps = 25L, seed = 12345L) {
  ok <- is.finite(pred) & is.finite(truth) & pred > 0 & pred < 1
  pred <- pred[ok]; truth <- truth[ok]
  n <- length(pred)
  if (n < 12) return(NULL)
  lp <- stats::qlogis(pmin(pmax(pred, 1e-4), 1 - 1e-4))
  lt <- stats::qlogis(pmin(pmax(truth, 1e-4), 1 - 1e-4))

  set.seed(seed)
  rows <- list()
  for (k in ks) {
    if (k >= n - 4) next
    for (r in seq_len(reps)) {
      sent <- if (k == 0) integer(0) else sample.int(n, k)
      held <- setdiff(seq_len(n), sent)
      # k = 0 is the uncalibrated transported map, the reference point.
      shift <- if (k == 0) 0 else mean(lt[sent] - lp[sent])
      p <- stats::plogis(lp[held] + shift)
      rows[[length(rows) + 1L]] <- data.frame(
        k_sentinel = k, replicate = r, n_held = length(held),
        rmse_pp = 100 * sqrt(mean((p - truth[held])^2)),
        mae_pp  = 100 * mean(abs(p - truth[held])),
        shift_logit = shift, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}
