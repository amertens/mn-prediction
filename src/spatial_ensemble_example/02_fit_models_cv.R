# 02_fit_models_cv.R
# Spatial block CV + comparison of alternative spatial-ensemble strategies.
#
# Methods:
# (A) Non-spatial stacking (Super Learner-style convex combination of base learners)
# (B) Stacking with coordinates available to the base learners
# (C) Two-stage: non-spatial stacking + spatial meta-learner (GAM "gp" smooth)
# (D) Single spatial model with linear covariate mean + spatial random field (GAM "gp")
#
# NOTES ON BEST PRACTICE:
# - Predictive accuracy metrics are computed using cross-validated out-of-fold predictions.
# - Mapping models (fit on the full data for wall-to-wall prediction) should not be used to
#   compute CV metrics unless they are refit within each CV split. This script adheres to that.
#
# Outputs:
# - outputs/cv_predictions.csv
# - outputs/cv_metrics.csv

library(here)
library(tidyverse)
stopifnot(requireNamespace("data.table", quietly = TRUE))
stopifnot(requireNamespace("sf", quietly = TRUE))
stopifnot(requireNamespace("blockCV", quietly = TRUE))
stopifnot(requireNamespace("mgcv", quietly = TRUE))
stopifnot(requireNamespace("glmnet", quietly = TRUE))
stopifnot(requireNamespace("ranger", quietly = TRUE))
stopifnot(requireNamespace("xgboost", quietly = TRUE))

dt_sim <- data.table::fread("data/sim/cluster_data.csv")
head(dt_sim)


dt <- readRDS(here("data", "sim", "Ghana_gee_wfp_prescreened_data.rds"))
dt <- dt %>% rename(cluster_id=gw_cn,
                    x1=wfp_rice_local,
                    x2=gee_built_surface,
                    x3=gee_smod_code,
                    rs_cov1=gee_trmm_mean_50km,
                    rs_cov2=gee_ndvi_mean_50km,
                    rs_cov3=gee_temp_mean_50km)

head(dt)
table(is.na(dt))

dt$N2 <- dt$N
dt$N <- round(dt$Y / dt$prev)
dt$N[is.na(dt$N)] <- dt$N2[is.na(dt$N)]



# ---- Spatial block cross-validation folds ----
# Inputs required:
# - coordinates (lon, lat)
# - a block size (in the same units as coordinates)
# - desired number of folds
# ------------------------------------------------------------------------------
# Spatial block cross-validation for point-level survey data (sf-based)
#
# PURPOSE
#   Create spatially structured CV folds for cluster-level survey data
#   without rasterization or implicit CRS assumptions.
#
# KEY DESIGN PRINCIPLES
#   - Operates on point geometries only (no rasters created)
#   - Forces explicit CRS definition and projection
#   - Block size is defined in meters
#   - Returns a simple fold-id vector for drop-in compatibility
#
# DEPENDENCIES
#   sf
#   blockCV
#
# INPUTS
#   dt          : data.frame with longitude/latitude columns
#   lon         : longitude column name (character)
#   lat         : latitude column name (character)
#   k           : number of folds
#   block_size  : spatial block size in meters
#   crs_ll      : CRS of input coordinates (default WGS84)
#   crs_proj    : projected CRS for distance-based blocking (default Web Mercator)
#   seed        : RNG seed for reproducibility
#
# OUTPUT
#   Integer vector of length nrow(dt) with fold IDs in {1, ..., k}
# ------------------------------------------------------------------------------
make_spatial_folds_sf <- function(
    dt,
    lon = "lon",
    lat = "lat",
    k = 5,
    block_size = 200000,
    crs_ll = 4326,
    crs_proj = 3857,
    seed = 123
) {

  stopifnot(
    is.data.frame(dt),
    lon %in% names(dt),
    lat %in% names(dt)
  )

  requireNamespace("sf", quietly = TRUE)
  requireNamespace("blockCV", quietly = TRUE)

  # Convert to sf with explicit CRS
  sf_pts <- sf::st_as_sf(
    dt,
    coords = c(lon, lat),
    crs = crs_ll,
    remove = FALSE
  )

  # Project to planar CRS (meters)
  sf_pts <- sf::st_transform(sf_pts, crs_proj)

  # Compute spatial extent
  bbox <- sf::st_bbox(sf_pts)
  domain_size <- max(
    bbox$xmax - bbox$xmin,
    bbox$ymax - bbox$ymin
  )

  # Sanity check on block size
  if (block_size < 0.02 * domain_size) {
    stop(
      paste0(
        "block_size is too small relative to the spatial extent of the data.\n",
        "Provided block_size = ", block_size, " meters.\n",
        "Spatial extent ≈ ", round(domain_size), " meters.\n",
        "Use block sizes on the order of tens to hundreds of kilometers ",
        "(e.g., 50000–300000 meters)."
      ),
      call. = FALSE
    )
  }

  set.seed(seed)

  cv_obj <- tryCatch(
    blockCV::cv_spatial(
      x = sf_pts,
      k = k,
      size = block_size,
      selection = "random",
      iteration = 100,
      progress = FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "Spatial blocking failed.\n",
          "This is usually caused by an incompatible block_size.\n",
          "Original error: ", e$message
        ),
        call. = FALSE
      )
    }
  )

  folds <- cv_obj$folds_ids

  if (length(folds) != nrow(dt)) {
    stop(
      "Fold assignment length does not match number of observations.",
      call. = FALSE
    )
  }

  return(folds)
}


folds <- make_spatial_folds_sf(
  dt,
  k = 5,
  block_size = 100000
)
length(folds )
folds


# ---- Helpers ----
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))
clip01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

# Weighted Brier score (probability scale)
brier_score <- function(y, n, p_hat) {
  w <- n
  mean(w * (p_hat - (y / n))^2) / mean(w)
}

# Weighted CV R^2 on prevalence scale (cluster proportions), using weights n
cv_r2 <- function(y, n, p_hat) {
  w <- n
  ybar <- sum(w * (y / n)) / sum(w)
  sse <- sum(w * ((y / n) - p_hat)^2)
  sst <- sum(w * ((y / n) - ybar)^2)
  1 - sse / sst
}

# Calibration diagnostics:
# - intercept and slope from a binomial GLM: cbind(Y, N-Y) ~ offset(0) + logit(p_hat)
#   In practice, a flexible calibration curve can be used; here we also provide a simple ECE.
calibration_diagnostics <- function(y, n, p_hat, n_bins = 10) {
  p_hat <- clip01(p_hat)
  z <- logit(p_hat)

  # Fit: logit(E[Y/N]) = a + b * logit(p_hat)
  fit <- suppressWarnings(stats::glm(
    cbind(y, n - y) ~ z,
    family = stats::binomial(),
    weights = n
  ))
  a <- stats::coef(fit)[1]
  b <- stats::coef(fit)[2]

  # Expected calibration error (ECE) with weighted bins
  bins <- cut(p_hat, breaks = quantile(p_hat, probs = seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  dtb <- data.table::data.table(bin = bins, y = y, n = n, p = p_hat)
  agg <- dtb[, .(
    p_mean = sum(n * p) / sum(n),
    y_mean = sum(y) / sum(n),
    w = sum(n)
  ), by = bin]
  ece <- sum(agg$w * abs(agg$y_mean - agg$p_mean)) / sum(agg$w)

  list(intercept = unname(a), slope = unname(b), ece = unname(ece))
}

# Convex stacking weights: minimize weighted Brier score under w>=0 and sum(w)=1
# Uses quadratic programming via base 'quadprog' if available; otherwise falls back to NNLS + renormalization.
stack_weights_convex <- function(P, y, n) {
  # P: matrix (n_obs x n_learners) of predictions on probability scale
  # y, n: observed counts
  stopifnot(is.matrix(P))
  p <- ncol(P)
  wts <- as.numeric(n)

  # Objective: minimize sum_i w_i ( (P_i %*% w) - y_i/n_i )^2
  # Equivalent QP: (1/2) w' D w - d' w, with constraints w>=0, sum(w)=1
  yprop <- y / n

  Dmat <- t(P) %*% (P * wts)  # p x p
  dvec <- t(P) %*% (wts * yprop)

  # Add small ridge for numerical stability
  Dmat <- Dmat + diag(1e-8, p)

  if (requireNamespace("quadprog", quietly = TRUE)) {
    # Constraints: A' w >= b
    # w >= 0  => I' w >= 0
    # sum(w)=1 => enforce with two inequalities: sum(w) >= 1 and -sum(w) >= -1
    Amat <- cbind(diag(p), rep(1, p), -rep(1, p))
    bvec <- c(rep(0, p), 1, -1)

    sol <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0)
    w <- sol$solution
    w[w < 0] <- 0
    w <- w / sum(w)
    return(w)
  }

  # Fallback: non-negative least squares then renormalize
  if (requireNamespace("nnls", quietly = TRUE)) {
    fit <- nnls::nnls(A = P * sqrt(wts), b = yprop * sqrt(wts))
    w <- pmax(coef(fit), 0)
    w <- w / sum(w)
    return(w)
  }

  # Last resort: equal weights
  rep(1 / p, p)
}

# ---- Base learners ----
# Each learner returns predicted probabilities.
fit_predict_glm <- function(train, test, use_coords = FALSE) {
  f <- if (use_coords) {
    cbind(Y, N - Y) ~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3 + lon + lat
  } else {
    cbind(Y, N - Y) ~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3
  }
  m <- suppressWarnings(stats::glm(f, family = stats::binomial(), data = train))
  as.numeric(stats::predict(m, newdata = test, type = "response"))
}

fit_predict_gam <- function(train, test, use_coords = FALSE) {
  # Spatial info as smooths only when use_coords=TRUE (non-spatial case uses smooths on covariates)
  if (use_coords) {
    f <- cbind(Y, N - Y) ~ s(x1, k = 6) + x2 + s(x3, k = 6) +
      s(rs_cov1, k = 6) + s(rs_cov2, k = 6) + s(rs_cov3, k = 6) +
      s(lon, lat, bs = "tp", k = 60)
  } else {
    f <- cbind(Y, N - Y) ~ s(x1, k = 6) + x2 + s(x3, k = 6) +
      s(rs_cov1, k = 6) + s(rs_cov2, k = 6) + s(rs_cov3, k = 6)
  }
  m <- mgcv::gam(f, family = stats::binomial(), data = train, method = "REML")
  as.numeric(stats::predict(m, newdata = test, type = "response"))
}

fit_predict_glmnet <- function(train, test, use_coords = FALSE) {
  X_train <- model.matrix(~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3 + (if (use_coords) lon + lat else 0),
                          data = train)[, -1, drop = FALSE]
  X_test <- model.matrix(~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3 + (if (use_coords) lon + lat else 0),
                         data = test)[, -1, drop = FALSE]

  y <- train$Y / train$N
  w <- train$N

  m <- glmnet::cv.glmnet(
    x = X_train, y = y,
    family = "binomial",
    weights = w,
    alpha = 1,
    nfolds = 5,
    type.measure = "deviance"
  )

  p <- as.numeric(stats::predict(m, newx = X_test, s = "lambda.min", type = "response"))
  clip01(p)
}

fit_predict_ranger <- function(train, test, use_coords = FALSE) {
  feats <- c("x1", "x2", "x3", "rs_cov1", "rs_cov2", "rs_cov3", if (use_coords) c("lon", "lat") else NULL)
  # Train on proportions with case weights N (approximation)
  df_tr <- data.frame(train[, c("Y", "N", feats), with = FALSE])
  df_tr$prev <- df_tr$Y / df_tr$N

  rf <- ranger::ranger(
    dependent.variable.name = "prev",
    data = df_tr[, c("prev", feats), drop = FALSE],
    num.trees = 500,
    min.node.size = 5,
    respect.unordered.factors = TRUE,
    case.weights = df_tr$N,
    seed = 123
  )
  df_te <- data.frame(test[, feats, with = FALSE])
  p <- as.numeric(predict(rf, data = df_te)$predictions)
  clip01(p)
}

fit_predict_xgb <- function(train, test, use_coords = FALSE) {
  feats <- c("x1", "x2", "x3", "rs_cov1", "rs_cov2", "rs_cov3", if (use_coords) c("lon", "lat") else NULL)

  X_tr <- as.matrix(train[, feats, with = FALSE])
  X_te <- as.matrix(test[, feats, with = FALSE])
  y_tr <- train$Y / train$N
  w_tr <- train$N

  dtr <- xgboost::xgb.DMatrix(data = X_tr, label = y_tr, weight = w_tr)
  dte <- xgboost::xgb.DMatrix(data = X_te)

  params <- list(
    objective = "reg:squarederror",
    eta = 0.05,
    max_depth = 4,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1
  )

  m <- xgboost::xgb.train(
    params = params,
    data = dtr,
    nrounds = 400,
    verbose = 0
  )

  p <- as.numeric(predict(m, dte))
  clip01(p)
}

# ---- Method implementations (CV) ----
get_oof_predictions <- function(dt, fold_id, use_coords = FALSE) {

  stopifnot(length(fold_id) == nrow(dt))

  preds <- data.table::data.table(cluster_id = dt$cluster_id)
  preds[, `:=`(
    p_glm = NA_real_,
    p_gam = NA_real_,
    p_rf  = NA_real_,
    p_xgb = NA_real_
  )]

  fold_levels <- sort(unique(fold_id))

  for (k in fold_levels) {

    test_idx  <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    dt <- as.data.frame(dt)
    tr <- dt[dt$cluster_id %in% train_idx,]
    te <- dt[dt$cluster_id %in% test_idx,]
    tr <- data.table::data.table(tr)
    te <- data.table::data.table(te)

    preds[test_idx, p_glm := fit_predict_glm(tr, te, use_coords = use_coords)]
    preds[test_idx, p_gam := fit_predict_gam(tr, te, use_coords = use_coords)]
    preds[test_idx, p_rf  := fit_predict_ranger(tr, te, use_coords = use_coords)]
    preds[test_idx, p_xgb := fit_predict_xgb(tr, te, use_coords = use_coords)]
  }

  return(preds)
}

# (A) Non-spatial stacking: base learners ignore coordinates; stacking weights estimated within each fold
fit_method_A <- function(dt, fold_id) {

  oof <- get_oof_predictions(dt, fold_id, use_coords = FALSE)
  Pcols <- c("p_glm", "p_gam", "p_rf", "p_xgb")

  out <- data.table::data.table(
    cluster_id = dt$cluster_id,
    p_hat = NA_real_,
    method = "A_nonspatial_stack"
  )

  fold_levels <- sort(unique(fold_id))

  for (k in fold_levels) {

    test_idx  <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    P_tr <- as.matrix(oof[train_idx, ..Pcols])
    P_te <- as.matrix(oof[test_idx,  ..Pcols])

    w <- stack_weights_convex(
      P_tr,
      y = dt$Y[train_idx],
      n = dt$N[train_idx]
    )

    out[test_idx, p_hat := as.numeric(P_te %*% w)]
  }

  list(pred = out, base_oof = oof)
}

# (B) Coordinates in base learners; stacking weights estimated within each fold
fit_method_B <- function(dt, fold_id) {

  oof <- get_oof_predictions(dt, fold_id, use_coords = TRUE)
  Pcols <- c("p_glm", "p_gam", "p_rf", "p_xgb")

  out <- data.table::data.table(
    cluster_id = dt$cluster_id,
    p_hat = NA_real_,
    method = "B_stack_with_coords"
  )

  fold_levels <- sort(unique(fold_id))

  for (k in fold_levels) {

    test_idx  <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    P_tr <- as.matrix(oof[train_idx, ..Pcols])
    P_te <- as.matrix(oof[test_idx,  ..Pcols])

    w <- stack_weights_convex(
      P_tr,
      y = dt$Y[train_idx],
      n = dt$N[train_idx]
    )

    out[test_idx, p_hat := as.numeric(P_te %*% w)]
  }

  list(pred = out, base_oof = oof)
}

fit_method_C <- function(dt, fold_id) {

  requireNamespace("mgcv", quietly = TRUE)
  requireNamespace("data.table", quietly = TRUE)

  # First-stage OOF predictions (non-spatial base learners)
  oof <- get_oof_predictions(dt, fold_id, use_coords = FALSE)

  Pcols <- c("p_glm", "p_gam", "p_rf", "p_xgb")

  dt2 <- cbind(dt, oof)

  out <- data.table::data.table(
    cluster_id = dt$cluster_id,
    p_hat = NA_real_,
    method = "C_stack_plus_spatial_meta"
  )

  fold_levels <- sort(unique(fold_id))

  for (k in fold_levels) {

    train_idx <- which(fold_id != k)
    test_idx  <- which(fold_id == k)

    dt2 <- as.data.frame(dt2)
    tr <- dt2[dt$cluster_id %in% train_idx,]
    te <- dt2[dt$cluster_id %in% test_idx,]
    tr <- data.table::data.table(tr)
    te <- data.table::data.table(te)

    # Number of unique spatial locations in training data
    n_loc <- length(unique(paste(tr$lon, tr$lat)))

    # Adaptive basis size for spatial smooth
    k_gp <- max(5, min(50, n_loc - 1))

    fit <- mgcv::gam(
      cbind(Y, N - Y) ~
        p_glm + p_gam + p_rf + p_xgb +
        s(lon, lat, bs = "gp", k = k_gp),
      family = stats::binomial(),
      data = tr,
      method = "REML",
      select = TRUE
    )

    p <- as.numeric(stats::predict(fit, newdata = te, type = "response"))
    out[test_idx, p_hat := clip01(p)]
  }

  list(
    pred = out,
    base_oof = oof
  )
}

# (D) Single spatial model (no ensemble): linear mean + spatial random field
fit_method_D <- function(dt, fold_id) {

  out <- data.table::data.table(
    cluster_id = dt$cluster_id,
    p_hat = NA_real_,
    method = "D_spatial_gam_only"
  )

  fold_levels <- sort(unique(fold_id))

  for (k in fold_levels) {

    test_idx  <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    dt <- as.data.frame(dt)
    tr <- dt[dt$cluster_id %in% train_idx,]
    te <- dt[dt$cluster_id %in% test_idx,]
    tr <- data.table::data.table(tr)
    te <- data.table::data.table(te)


    m <- mgcv::gam(
      cbind(Y, N - Y) ~
        x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3 +
        s(lon, lat, bs = "gp", k = 5),
      family = stats::binomial(),
      data = tr,
      method = "REML"
    )

    p <- as.numeric(stats::predict(m, newdata = te, type = "response"))
    out[test_idx, p_hat := clip01(p)]
  }

  list(pred = out)
}

# ---- Fit all methods ----
resA <- fit_method_A(dt, folds)
resB <- fit_method_B(dt, folds)
resC <- fit_method_C(dt, folds)
resD <- fit_method_D(dt, folds)

pred_all <- data.table::rbindlist(list(resA$pred, resB$pred, resC$pred, resD$pred))
dt <- data.table::data.table(dt)

# Attach observed
pred_all <- merge(pred_all, dt[, .(cluster_id, Y, N, prev, lon, lat)], by = "cluster_id", all.x = TRUE)
pred_all


# ---- Compute metrics ----
metric_one <- function(dsub) {
  p_hat <- clip01(dsub$p_hat)
  list(
    brier = brier_score(dsub$Y, dsub$N, p_hat),
    r2 = cv_r2(dsub$Y, dsub$N, p_hat),
    cal_intercept = calibration_diagnostics(dsub$Y, dsub$N, p_hat)$intercept,
    cal_slope = calibration_diagnostics(dsub$Y, dsub$N, p_hat)$slope,
    ece = calibration_diagnostics(dsub$Y, dsub$N, p_hat)$ece
  )
}

metrics <- pred_all[, {
  m <- metric_one(.SD)
  .(brier = m$brier, cv_r2 = m$r2, cal_intercept = m$cal_intercept, cal_slope = m$cal_slope, ece = m$ece)
}, by = method]

data.table::fwrite(pred_all, here("outputs/cv_predictions.csv"))
data.table::fwrite(metrics, here("outputs/cv_metrics.csv"))

message("CV complete. See outputs/cv_metrics.csv and outputs/cv_predictions.csv")


# ---- Optional: quick calibration plot for one method ----
# (Uncomment to view)
# library(ggplot2)
# plot_calibration <- function(dsub, n_bins = 10) {
#   dsub[, p_hat := clip01(p_hat)]
#   bins <- cut(p_hat, breaks = quantile(p_hat, probs = seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
#   agg <- dsub[, .(p_mean = sum(N * p_hat) / sum(N), y_mean = sum(Y) / sum(N), w = sum(N)), by = bins]
#   ggplot(agg, aes(x = p_mean, y = y_mean, size = w)) +
#     geom_point(alpha = 0.8) +
#     geom_abline(intercept = 0, slope = 1, linetype = 2) +
#     coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(x = "Mean predicted prevalence", y = "Mean observed prevalence", title = unique(dsub$method))
# }
# print(plot_calibration(pred_all[method == "C_two_stage_spatial_meta"]))
