# =============================================================================
# R/fit_superlearner.R
#
# SuperLearner ensemble modeling with cluster-blocked cross-validation.
#
# Key functions:
#   make_folds_cluster    – cluster-blocked CV fold assignments
#   SL.mgcv_gam           – custom mgcv GAM wrapper (CRAN-only)
#   SL.gam_spatial        – GAM with 2-D thin-plate spline on coordinates
#   fit_sl_continuous     – SuperLearner for continuous biomarker Y
#   fit_sl_binary         – SuperLearner for binary deficiency indicator D
#   get_cv_preds          – extract cross-validated predictions
#   summarise_sl_weights  – extract and display metalearner weights
# =============================================================================

# ── Cluster-blocked cross-validation ─────────────────────────────────────────

#' Create cluster-blocked CV fold assignments for SuperLearner.
#'
#' All individuals in a cluster are assigned to the same validation fold,
#' preventing data leakage across the cluster boundary.
#'
#' @param cluster_id Integer or character vector of cluster memberships.
#' @param V          Number of folds (default 5).
#' @param seed       RNG seed for fold assignment.
#' @return A list of length V, each element a numeric vector of row indices
#'         comprising the VALIDATION set for that fold.  Suitable for passing
#'         to SuperLearner via \code{cvControl = list(validRows = ...)}.
make_folds_cluster <- function(cluster_id, V = 5, seed = 42) {
  set.seed(seed)
  unique_clusters <- unique(cluster_id)
  n_cl            <- length(unique_clusters)

  # Randomly assign clusters to folds (balanced as possible)
  fold_labels <- sample(rep_len(seq_len(V), n_cl))
  names(fold_labels) <- as.character(unique_clusters)

  # Translate cluster fold to individual fold
  indiv_fold <- fold_labels[as.character(cluster_id)]

  # Return validation row-index lists
  lapply(seq_len(V), function(v) which(indiv_fold == v))
}


# ── Custom SuperLearner wrappers ──────────────────────────────────────────────

#' SuperLearner wrapper: mgcv GAM with penalised smooth terms.
#'
#' Fits a GAM using mgcv::gam() with s(., k=4) for each numeric predictor
#' (up to max_smooth terms to prevent over-parameterisation).
#' Falls back to a GLM if fewer than 2 predictors are available.
#'
#' @inheritParams SuperLearner::SL.glmnet
#' @param max_smooth Maximum number of smooth terms (default 8).
SL.mgcv_gam <- function(Y, X, newX, family, obsWeights,
                          max_smooth = 8, ...) {
  require(mgcv, quietly = TRUE)

  num_cols  <- names(X)[sapply(X, function(x) is.numeric(x) && length(unique(x)) > 5)]
  bin_cols  <- setdiff(names(X), num_cols)

  smooth_terms  <- paste0("s(", head(num_cols, max_smooth), ", k = 4)",
                           collapse = " + ")
  linear_terms  <- if (length(bin_cols) > 0) paste(bin_cols, collapse = " + ")
                   else ""
  all_terms     <- paste(c(smooth_terms, linear_terms)[nchar(c(smooth_terms, linear_terms)) > 0],
                          collapse = " + ")
  if (nchar(all_terms) == 0) all_terms <- "1"

  formula_str <- paste("Y ~", all_terms)
  fit_data    <- cbind(Y = Y, as.data.frame(X))

  fit <- tryCatch(
    mgcv::gam(as.formula(formula_str), data = fit_data,
               family = family, weights = obsWeights, method = "REML"),
    error = function(e) {
      # Fallback: linear model
      glm(as.formula(paste("Y ~", paste(names(X), collapse = " + "))),
           data = fit_data, family = family, weights = obsWeights)
    }
  )

  pred <- as.numeric(predict(fit, newdata = as.data.frame(newX),
                              type = "response"))
  pred <- pmax(0, pmin(if (identical(family, binomial())) 1 else Inf, pred))

  list(pred = pred, fit = fit)
}

#' Prediction method for SL.mgcv_gam.
predict.SL.mgcv_gam <- function(object, newdata, ...) {
  as.numeric(predict(object$fit, newdata = as.data.frame(newdata),
                     type = "response"))
}


#' SuperLearner wrapper: spatially-aware GAM.
#'
#' Adds a 2-D thin-plate spline on (x_coord, y_coord) to capture residual
#' spatial structure that other learners may miss.
#'
#' @inheritParams SuperLearner::SL.glmnet
SL.gam_spatial <- function(Y, X, newX, family, obsWeights, ...) {
  require(mgcv, quietly = TRUE)

  has_coords <- all(c("x_coord", "y_coord") %in% names(X))
  other_cols <- setdiff(names(X), c("x_coord", "y_coord"))

  if (has_coords && length(other_cols) > 0) {
    formula_str <- paste("Y ~ s(x_coord, y_coord, bs = 'tp', k = 15) +",
                          paste(other_cols, collapse = " + "))
  } else if (has_coords) {
    formula_str <- "Y ~ s(x_coord, y_coord, bs = 'tp', k = 15)"
  } else if (length(other_cols) > 0) {
    formula_str <- paste("Y ~", paste(other_cols, collapse = " + "))
  } else {
    formula_str <- "Y ~ 1"
  }

  fit_data <- cbind(Y = Y, as.data.frame(X))
  fit <- tryCatch(
    mgcv::gam(as.formula(formula_str), data = fit_data,
               family = family, weights = obsWeights, method = "REML"),
    error = function(e) {
      lm(Y ~ ., data = fit_data)
    }
  )

  pred <- as.numeric(predict(fit, newdata = as.data.frame(newX),
                              type = "response"))
  pred <- pmax(0, pmin(if (identical(family, binomial())) 1 else Inf, pred))

  list(pred = pred, fit = fit)
}

#' Prediction method for SL.gam_spatial.
predict.SL.gam_spatial <- function(object, newdata, ...) {
  as.numeric(predict(object$fit, newdata = as.data.frame(newdata),
                     type = "response"))
}


# ── Fast ranger wrapper (fewer trees for tutorial speed) ─────────────────────

#' SuperLearner wrapper: ranger random forest (speed-optimised for tutorial).
SL.ranger_fast <- function(Y, X, newX, family, obsWeights,
                             num.trees = 200, min.node.size = 5, ...) {
  require(ranger, quietly = TRUE)
  df_fit <- cbind(Y = Y, as.data.frame(X))
  fit    <- ranger::ranger(Y ~ ., data = df_fit,
                             num.trees     = num.trees,
                             min.node.size = min.node.size,
                             case.weights  = obsWeights,
                             probability   = identical(family, binomial()),
                             seed          = 1)
  if (identical(family, binomial())) {
    pred <- predict(fit, data = as.data.frame(newX))$predictions[, 2]
  } else {
    pred <- predict(fit, data = as.data.frame(newX))$predictions
  }
  list(pred = as.numeric(pred), fit = fit)
}

#' Prediction method for SL.ranger_fast.
predict.SL.ranger_fast <- function(object, newdata, family, ...) {
  if (inherits(object$fit, "ranger") && object$fit$treetype == "Probability estimation") {
    predict(object$fit, data = as.data.frame(newdata))$predictions[, 2]
  } else {
    predict(object$fit, data = as.data.frame(newdata))$predictions
  }
}


# ── Core SL fitting functions ─────────────────────────────────────────────────

#' Fit a SuperLearner ensemble for a CONTINUOUS biomarker outcome.
#'
#' Uses squared-error loss and a cluster-blocked CV scheme defined by
#' \code{foldid} (output of \code{make_folds_cluster}).
#'
#' @param df         Data frame containing outcome and predictors.
#' @param outcome    Name of the continuous outcome column (e.g., "Y").
#' @param X_cols     Character vector of predictor column names.
#' @param foldid     Validation row-index list from \code{make_folds_cluster}.
#' @param sl_library Character vector of SuperLearner wrapper names.
#' @param ...        Additional arguments passed to \code{SuperLearner::SuperLearner}.
#' @return A fitted SuperLearner object.
fit_sl_continuous <- function(df,
                               outcome    = "Y",
                               X_cols     = NULL,
                               foldid     = NULL,
                               sl_library = c("SL.mean", "SL.glmnet",
                                              "SL.ranger_fast",
                                              "SL.mgcv_gam",
                                              "SL.gam_spatial"),
                               ...) {
  if (is.null(X_cols))
    X_cols <- get_predictor_cols(df)

  # Drop all-NA columns
  X_cols <- X_cols[X_cols %in% colnames(df)]
  X_cols <- X_cols[!sapply(df[, X_cols, drop = FALSE],
                            function(x) all(is.na(x)))]

  Y <- df[[outcome]]
  X <- df[, X_cols, drop = FALSE]

  cv_control <- if (!is.null(foldid)) list(validRows = foldid) else list(V = 5L)

  SuperLearner::SuperLearner(
    Y          = Y,
    X          = X,
    family     = gaussian(),
    SL.library = sl_library,
    cvControl  = cv_control,
    ...
  )
}


#' Fit a SuperLearner ensemble for a BINARY deficiency outcome.
#'
#' Uses binomial negative log-likelihood (logistic) loss.  The metalearner
#' is constrained to the probability simplex (non-negative weights summing
#' to 1) via the default NNLS metalearner.
#'
#' @param df         Data frame.
#' @param outcome    Name of the binary outcome column (e.g., "D").
#' @param X_cols     Predictor column names.
#' @param foldid     Validation row-index list from \code{make_folds_cluster}.
#' @param sl_library SuperLearner wrapper names.
#' @param ...        Additional arguments to SuperLearner.
#' @return A fitted SuperLearner object.
fit_sl_binary <- function(df,
                           outcome    = "D",
                           X_cols     = NULL,
                           foldid     = NULL,
                           sl_library = c("SL.mean", "SL.glmnet",
                                          "SL.ranger_fast",
                                          "SL.mgcv_gam",
                                          "SL.gam_spatial"),
                           ...) {
  if (is.null(X_cols))
    X_cols <- get_predictor_cols(df)

  X_cols <- X_cols[X_cols %in% colnames(df)]
  X_cols <- X_cols[!sapply(df[, X_cols, drop = FALSE],
                            function(x) all(is.na(x)))]

  Y <- as.numeric(df[[outcome]])
  X <- df[, X_cols, drop = FALSE]

  cv_control <- if (!is.null(foldid)) list(validRows = foldid) else list(V = 5L)

  SuperLearner::SuperLearner(
    Y          = Y,
    X          = X,
    family     = binomial(),
    SL.library = sl_library,
    cvControl  = cv_control,
    ...
  )
}


# ── Cross-validated prediction extraction ─────────────────────────────────────

#' Extract cross-validated (out-of-fold) predictions from a fitted SuperLearner.
#'
#' @param sl_fit A fitted SuperLearner object.
#' @return Numeric vector of cross-validated predictions (length = n).
get_cv_preds <- function(sl_fit) {
  as.numeric(sl_fit$SL.predict)
}

#' Extract individual base-learner CV predictions from a SuperLearner.
#' @param sl_fit A fitted SuperLearner object.
#' @return Data frame: one column per base learner, one row per individual.
get_cv_preds_all <- function(sl_fit) {
  as.data.frame(sl_fit$library.predict)
}


# ── Metalearner weights ───────────────────────────────────────────────────────

#' Extract and display metalearner (ensemble) coefficient weights.
#'
#' @param sl_fit A fitted SuperLearner object.
#' @return Data frame: learner name, weight, cv_risk.
summarise_sl_weights <- function(sl_fit) {
  risk_tbl <- sl_fit$cvRisk
  weights   <- coef(sl_fit)
  learners  <- names(weights)

  data.frame(
    learner  = learners,
    weight   = round(weights, 4),
    cv_risk  = round(risk_tbl[learners], 6),
    row.names = NULL
  ) %>%
    dplyr::arrange(dplyr::desc(weight))
}


# ── Nested CV for threshold tuning ───────────────────────────────────────────

#' Inner-loop threshold tuner: find T* that maximises a criterion.
#'
#' @param obs    Observed binary outcomes (training fold).
#' @param pred   Predicted probabilities (training fold).
#' @param criterion "f1", "youden", or "cost".
#' @param cfp    Cost of false positive (used when criterion = "cost").
#' @param cfn    Cost of false negative (used when criterion = "cost").
#' @return Scalar optimal threshold.
tune_threshold <- function(obs, pred,
                            criterion = c("f1", "youden", "cost"),
                            cfp = 1, cfn = 1) {
  criterion <- match.arg(criterion)
  thresholds <- seq(0.05, 0.95, by = 0.01)

  scores <- sapply(thresholds, function(t) {
    pred_class <- as.integer(pred >= t)
    tp <- sum(pred_class == 1 & obs == 1)
    fp <- sum(pred_class == 1 & obs == 0)
    tn <- sum(pred_class == 0 & obs == 0)
    fn <- sum(pred_class == 0 & obs == 1)
    sens <- tp / (tp + fn + 1e-9)
    spec <- tn / (tn + fp + 1e-9)
    prec <- tp / (tp + fp + 1e-9)

    if (criterion == "f1") {
      2 * prec * sens / (prec + sens + 1e-9)
    } else if (criterion == "youden") {
      sens + spec - 1
    } else {  # cost-minimising
      -(cfp * fp + cfn * fn) / length(obs)
    }
  })

  thresholds[which.max(scores)]
}

#' Compute sensitivity, specificity, PPV, NPV at a given threshold.
#'
#' @param obs   Binary observed outcomes.
#' @param pred  Predicted probabilities.
#' @param T     Classification threshold.
#' @return Named numeric vector.
classification_metrics_at_T <- function(obs, pred, T = 0.5) {
  pred_class <- as.integer(pred >= T)
  tp <- sum(pred_class == 1 & obs == 1)
  fp <- sum(pred_class == 1 & obs == 0)
  tn <- sum(pred_class == 0 & obs == 0)
  fn <- sum(pred_class == 0 & obs == 1)
  c(
    sensitivity = tp / (tp + fn + 1e-9),
    specificity = tn / (tn + fp + 1e-9),
    ppv         = tp / (tp + fp + 1e-9),
    npv         = tn / (tn + fn + 1e-9),
    f1          = 2 * tp / (2 * tp + fp + fn + 1e-9)
  )
}
