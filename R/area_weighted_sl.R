# =============================================================================
# R/area_weighted_sl.R
#
# Precision-weighted area-level ensemble.
#
# WHY IT EXISTS
# -------------
# Districts contribute wildly unequal amounts of information -- 6 biomarker
# reads in one, 124 in another -- and the unweighted area-level SL treats them
# as equally informative. Weighting each area by the inverse of its sampling
# variance is standard small-area practice (it is what makes Fay-Herriot work)
# and in the ablation it was the single largest estimator gain measured:
#
#   config                          median r   r_share   MAE pp   beats base (r)
#   baseline (current 5-learner)      0.098      0.30     10.44         --
#   baseline + precision weights      0.224      0.66      9.78        17/24
#
# (scripts/covariates/09_estimator_ablation.R, leave-one-region-out over 24
# country x outcome cells; results/tables/estimator_ablation.csv.)
#
# WHY IT IS A SEPARATE FUNCTION RATHER THAN A FLAG ON fit_area_sl()
# -----------------------------------------------------------------
# mlr3superlearner's entry point takes no weights argument at all --
# (data, target, library, outcome_type, folds, discrete, newdata, group, info)
# -- so case weights cannot be passed through it, and forcing them would mean
# patching a third-party package inside a pipeline that has to stay
# reproducible. This implements the same five learner families with the same
# hyperparameters, stacked by non-negative least squares on inner out-of-fold
# predictions, which is what the SuperLearner ensemble step does.
#
# It is ADDITIVE: fit_area_sl() is untouched, so the existing analysis of record
# is unaffected and the two appear side by side in the comparison table. Which
# one becomes primary is then an evidence question answered by that table, not a
# decision buried in a refactor.
#
# The same ablation found that TRIMMING the library is harmful (enet alone:
# median r -0.088; linear-only: 0.028), so the five families are kept.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

#' Inverse-sampling-variance weights for area prevalences, on the logit scale.
#'
#' The model is fitted on logit(prevalence), so the weight has to be the
#' precision of the LOGIT, not of the proportion. Delta method:
#'   Var(logit(p)) ~= Var(p) / (p(1-p))^2
#' Uses the Agresti-Coull shrunk sampling variance so a district observed at 0%
#' or 100% does not receive infinite weight. Normalised to mean 1 so the
#' effective sample size, and hence any learner's internal regularisation, is
#' comparable to the unweighted fit.
area_precision_weights <- function(prev, n_svy, deff = 1.5) {
  p  <- pmin(pmax(as.numeric(prev), 0.005), 0.995)
  vp <- area_sampling_var_shrunk(as.numeric(prev), as.numeric(n_svy), deff)
  vl <- vp / (p * (1 - p))^2
  w  <- 1 / pmax(vl, .Machine$double.eps)
  w[!is.finite(w) | w <= 0] <- stats::median(w[is.finite(w) & w > 0])
  if (!any(is.finite(w))) return(rep(1, length(prev)))
  w / mean(w, na.rm = TRUE)
}

.awsl_glmnet <- function(X, y, w, alpha) {
  f <- tryCatch(glmnet::cv.glmnet(X, y, alpha = alpha, weights = w,
                                  nfolds = max(3L, min(8L, floor(length(y) / 4)))),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, newx = Z, s = "lambda.min")))
}
.awsl_ranger <- function(X, y, w) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  f <- tryCatch(ranger::ranger(.y ~ ., data = data.frame(.y = y, X), num.trees = 500,
                               min.node.size = 5, case.weights = w, num.threads = 1),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, data = as.data.frame(Z))$predictions))
}
.awsl_xgb <- function(X, y, w) {
  if (!requireNamespace("xgboost", quietly = TRUE)) return(NULL)
  f <- tryCatch(xgboost::xgboost(data = X, label = y, weight = w, max_depth = 2,
                                 eta = 0.05, nrounds = 200, min_child_weight = 3,
                                 subsample = 0.8, verbose = 0, nthread = 1),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, Z)))
}
.AWSL_LEARNERS <- list(
  lasso  = function(X, y, w) .awsl_glmnet(X, y, w, 1),
  ridge  = function(X, y, w) .awsl_glmnet(X, y, w, 0),
  enet   = function(X, y, w) .awsl_glmnet(X, y, w, 0.5),
  ranger = .awsl_ranger,
  xgb    = .awsl_xgb)

#' NNLS-stacked ensemble; weights learned from inner out-of-fold predictions.
.awsl_stack <- function(X, y, w, seed = 20260829L) {
  set.seed(seed)
  ls <- .AWSL_LEARNERS; n <- length(y)
  k <- max(3L, min(5L, floor(n / 5)))
  fold <- sample(rep_len(seq_len(k), n))
  P <- matrix(NA_real_, n, length(ls), dimnames = list(NULL, names(ls)))
  for (f in seq_len(k)) {
    i <- which(fold == f)
    for (j in seq_along(ls)) {
      m <- ls[[j]](X[-i, , drop = FALSE], y[-i], w[-i])
      if (!is.null(m)) P[i, j] <- m$predict(X[i, , drop = FALSE])
    }
  }
  ok <- stats::complete.cases(P)
  wts <- rep(1 / length(ls), length(ls))
  if (sum(ok) > length(ls) + 2) {
    nn <- if (requireNamespace("nnls", quietly = TRUE))
      tryCatch(nnls::nnls(P[ok, , drop = FALSE], y[ok]), error = function(e) NULL) else NULL
    if (!is.null(nn) && sum(nn$x) > 0) wts <- nn$x / sum(nn$x) else {
      cf <- tryCatch(stats::coef(stats::lm.fit(P[ok, , drop = FALSE], y[ok])),
                     error = function(e) NULL)
      if (!is.null(cf)) { cf[!is.finite(cf) | cf < 0] <- 0
                          if (sum(cf) > 0) wts <- cf / sum(cf) }
    }
  }
  models <- lapply(ls, function(L) L(X, y, w))
  keep <- !vapply(models, is.null, logical(1))
  if (!any(keep)) return(NULL)
  list(models = models[keep], wts = wts[keep] / max(sum(wts[keep]), 1e-12),
       learner_weights = stats::setNames(wts, names(ls)))
}
.awsl_predict <- function(s, Z) {
  if (is.null(s)) return(rep(NA_real_, nrow(Z)))
  P <- vapply(s$models, function(m) m$predict(Z), numeric(nrow(Z)))
  if (is.null(dim(P))) P <- matrix(P, nrow = nrow(Z))
  as.numeric(P %*% s$wts)
}
.awsl_matrix <- function(df, vars) {
  X <- as.matrix(df[, vars, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  X
}

#' Canonical covariates named by the literature set.
#' @return character vector; empty when the metadata file is absent
literature_covariate_cols <- function(df) {
  f <- here::here("metadata", "covariates", "literature_covariates.csv")
  if (!file.exists(f)) return(character(0))
  L <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  cols <- area_covariate_cols(df)
  unique(unlist(lapply(L$canonical_regex, function(r) grep(r, cols, value = TRUE))))
}

#' Keep the k covariates most correlated with the area outcome.
.awsl_screen <- function(X, y, k) {
  if (k <= 0 || k >= ncol(X)) return(colnames(X))
  r <- abs(suppressWarnings(stats::cor(X, y)))
  r[!is.finite(r)] <- 0
  colnames(X)[order(r, decreasing = TRUE)][seq_len(min(k, ncol(X)))]
}

#' Precision-weighted area-level ensemble with honest out-of-fold predictions.
#'
#' @param svy_admin2 district survey prevalences (Admin2, svy_prev, n_svy)
#' @param gee_admin2 admin-2 covariates (Admin1/Admin2 + numeric columns)
#' @param weighted apply precision weights (FALSE reproduces the unweighted
#'   ensemble, which is what makes the weighting effect measurable in-pipeline)
#' @param fold_by "region" uses leave-one-Admin1-out -- no district is predicted
#'   by a model that has seen its neighbours; falls back to 10-fold if Admin1 is
#'   unavailable or too coarse
#' @return list(area_preds, metrics, learner_weights, n_areas, n_vars)
#' @param covariate_set "all" (every canonical covariate) or "literature" (the
#'   hypothesis-driven set in metadata/covariates/literature_covariates.csv)
#' @param screen screening protocol:
#'   "none"      no screening
#'   "in_fold"   re-screen inside every training fold -- the HONEST protocol
#'   "all_data"  screen ONCE on the full data, then cross-validate the screened
#'               set. Deliberately OPTIMISTIC: the held-out areas helped choose
#'               the covariates, so the metric is an UPPER BOUND, not a
#'               performance estimate. It exists to QUANTIFY that optimism --
#'               the manuscript already documents one instance of exactly this
#'               bias (a recipe selected on the LOCO correlations it was then
#'               scored against), and the in_fold-vs-all_data gap measures how
#'               large the effect is here. Never quote it as accuracy.
#' @param k_screen covariates kept when screening
fit_area_weighted_sl <- function(svy_admin2, gee_admin2, weighted = TRUE,
                                 fold_by = "region", deff = 1.5,
                                 covariate_set = c("all", "literature"),
                                 screen = c("none", "in_fold", "all_data"),
                                 k_screen = 20L) {
  covariate_set <- match.arg(covariate_set); screen <- match.arg(screen)
  # Join on the PAIR key where both sides carry Admin1. `keys` was computed
  # here but never used, and the join ran on Admin2 alone -- which fans out
  # against a covariate table whose Admin-2 names are not unique (Malawi: 243
  # rows, 239 distinct names), turning 87 surveyed districts into 90 training
  # rows with three duplicate districts sharing one survey value.
  keys <- admin2_join_by(svy_admin2, gee_admin2)
  tr <- dplyr::inner_join(
    svy_admin2[, c(keys, "svy_prev", intersect("n_svy", names(svy_admin2)))],
    gee_admin2, by = keys)
  tr <- tr[is.finite(tr$svy_prev), , drop = FALSE]
  if (nrow(tr) < 10) {
    cat("  [area_wsl] too few areas (", nrow(tr), ") - skipping\n"); return(NULL)
  }
  if (!"n_svy" %in% names(tr)) tr$n_svy <- 30

  vars <- setdiff(names(tr), c(keys, "svy_prev", "n_svy", "country"))
  vars <- vars[vapply(vars, function(v) {
    x <- suppressWarnings(as.numeric(tr[[v]]))
    sum(is.finite(x)) > 2 && stats::sd(x, na.rm = TRUE) > 0 }, logical(1))]
  if (identical(covariate_set, "literature")) {
    lit <- literature_covariate_cols(tr)
    if (length(lit) >= 5) vars <- intersect(vars, lit) else
      cat("  [area_wsl] literature set fell back to all
")
  }
  if (length(vars) < 2) { cat("  [area_wsl] too few usable covariates\n"); return(NULL) }

  X <- .awsl_matrix(tr, vars)
  y <- stats::qlogis(pmin(pmax(tr$svy_prev, 0.005), 0.995))
  w <- if (weighted) area_precision_weights(tr$svy_prev, tr$n_svy, deff) else rep(1, nrow(tr))

  folds <- if (fold_by == "region" && "Admin1" %in% names(tr) &&
               dplyr::n_distinct(tr$Admin1) >= 3) as.character(tr$Admin1) else {
    set.seed(20260829L); as.character(sample(rep_len(seq_len(min(10L, nrow(tr))), nrow(tr))))
  }
  # OPTIMISTIC: choose covariates once using EVERY area, including the ones
  # that will be held out, then cross-validate that fixed set.
  if (identical(screen, "all_data")) {
    sel <- .awsl_screen(X, y, k_screen)
    X <- X[, sel, drop = FALSE]; vars <- sel
  }

  oof <- rep(NA_real_, nrow(tr))
  for (f in unique(folds)) {
    i <- which(folds == f)
    if (length(i) >= nrow(tr)) next
    Xtr <- X[-i, , drop = FALSE]; Xte <- X[i, , drop = FALSE]
    if (identical(screen, "in_fold")) {
      sel <- .awsl_screen(Xtr, y[-i], k_screen)
      Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
    }
    s <- .awsl_stack(Xtr, y[-i], w[-i])
    oof[i] <- .awsl_predict(s, Xte)
  }
  full_X <- if (identical(screen, "in_fold"))
    X[, .awsl_screen(X, y, k_screen), drop = FALSE] else X
  full <- .awsl_stack(full_X, y, w)

  preds <- data.frame(tr[, keys, drop = FALSE],
                      area_pred = stats::plogis(.awsl_predict(full, full_X)),
                      area_pred_oof = stats::plogis(oof),
                      stringsAsFactors = FALSE)
  ok <- is.finite(preds$area_pred_oof)
  metrics <- data.frame(
    n_areas = sum(ok), n_vars = length(vars), weighted = weighted,
    covariate_set = covariate_set, screen = screen,
    optimistic = identical(screen, "all_data"),
    fold_by = if (fold_by == "region" && dplyr::n_distinct(folds) >= 3) "region" else "random",
    pearson_r = if (sum(ok) > 3) suppressWarnings(stats::cor(tr$svy_prev[ok], preds$area_pred_oof[ok])) else NA_real_,
    mae_pp = if (any(ok)) 100 * mean(abs(tr$svy_prev[ok] - preds$area_pred_oof[ok])) else NA_real_,
    stringsAsFactors = FALSE)
  cat(sprintf("  [area_wsl] n=%d, p=%d, %s, set=%s, screen=%s%s | OOF (%s) r=%.3f MAE=%.2f pp
",
              nrow(tr), length(vars), if (weighted) "precision-weighted" else "unweighted",
              covariate_set, screen,
              if (identical(screen, "all_data")) " [OPTIMISTIC]" else "",
              metrics$fold_by, metrics$pearson_r, metrics$mae_pp))
  list(area_preds = preds, metrics = metrics,
       learner_weights = full$learner_weights, n_areas = nrow(tr), n_vars = length(vars))
}
