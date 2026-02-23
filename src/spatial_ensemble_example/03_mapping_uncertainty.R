# 03_mapping_uncertainty.R
# Mapping + uncertainty demonstration using Method (C): two-stage spatial meta-learner.
#
# Key principle:
# - CV metrics were computed in 02_fit_models_cv.R using split-specific refits.
# - For mapping, we refit the selected approach on ALL available data to maximize information
#   and then predict wall-to-wall. Uncertainty shown here pertains to the spatial meta-learner,
#   using mgcv's approximate Bayesian covariance of coefficients.
#
# Outputs:
# - outputs/map_mean.tif, outputs/map_sd.tif
# - outputs/map_mean.png, outputs/map_sd.png

library(here)
stopifnot(requireNamespace("data.table", quietly = TRUE))
stopifnot(requireNamespace("terra", quietly = TRUE))
stopifnot(requireNamespace("mgcv", quietly = TRUE))
stopifnot(requireNamespace("glmnet", quietly = TRUE))
stopifnot(requireNamespace("ranger", quietly = TRUE))
stopifnot(requireNamespace("xgboost", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))

dt <- data.table::fread(here("data/sim/cluster_data.csv"))
ras <- terra::rast(here("data/sim/domain_covariates.tif"))

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))
clip01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

# ---- Stage-1: refit ensemble on full data (non-spatial learners) ----
# For mapping, you typically refit each base learner on the full dataset.
# Inputs required:
# - cluster-level covariates for training
# - (if using raster covariates) extraction of raster values at clusters (done in simulation)
# - outcome (Y,N)

fit_base_full <- function(dt) {
  # Return a list of fitted models and a helper to predict on newdata.
  # Base learners match those in 02_fit_models_cv.R (non-spatial: no lon/lat).
  m_glm <- suppressWarnings(stats::glm(
    cbind(Y, N - Y) ~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3,
    family = stats::binomial(),
    data = dt
  ))

  X <- model.matrix(~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3, data = dt)[, -1, drop = FALSE]
  y <- dt$Y / dt$N
  w <- dt$N
  # m_glmnet <- glmnet::cv.glmnet(
  #   x = X, y = y, family = "binomial", weights = w, alpha = 1, nfolds = 5, type.measure = "deviance"
  # )

  feats <- c("x1", "x2", "x3", "rs_cov1", "rs_cov2", "rs_cov3")
  df_tr <- data.frame(dt[, c("Y", "N", feats), with = FALSE])
  df_tr$prev <- df_tr$Y / df_tr$N
  m_rf <- ranger::ranger(
    dependent.variable.name = "prev",
    data = df_tr[, c("prev", feats), drop = FALSE],
    num.trees = 500,
    min.node.size = 5,
    case.weights = df_tr$N,
    seed = 123
  )

  X_tr <- as.matrix(dt[, feats, with = FALSE])
  dtr <- xgboost::xgb.DMatrix(data = X_tr, label = y, weight = w)
  params <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 4, subsample = 0.8, colsample_bytree = 0.8)
  m_xgb <- xgboost::xgb.train(params = params, data = dtr, nrounds = 400, verbose = 0)

  m_gam <- mgcv::gam(
    cbind(Y, N - Y) ~ s(x1, k = 6) + x2 + s(x3, k = 6) +
      s(rs_cov1, k = 6) + s(rs_cov2, k = 6) + s(rs_cov3, k = 6),
    family = stats::binomial(),
    data = dt,
    method = "REML"
  )

  predict_all <- function(newdata) {
    p_glm <- as.numeric(stats::predict(m_glm, newdata = newdata, type = "response"))

    # Xn <- model.matrix(~ x1 + x2 + x3 + rs_cov1 + rs_cov2 + rs_cov3, data = newdata)[, -1, drop = FALSE]
    # p_glmnet <- as.numeric(stats::predict(m_glmnet, newx = Xn, s = "lambda.min", type = "response"))

    df_te <- data.frame(newdata[, feats, with = FALSE])
    p_rf <- as.numeric(predict(m_rf, data = df_te)$predictions)

    dte <- xgboost::xgb.DMatrix(data = as.matrix(newdata[, feats, with = FALSE]))
    p_xgb <- as.numeric(predict(m_xgb, dte))

    p_gam <- as.numeric(stats::predict(m_gam, newdata = newdata, type = "response"))

    cbind(p_glm = clip01(p_glm), #p_glmnet = clip01(p_glmnet),
          p_gam = clip01(p_gam), p_rf = clip01(p_rf), p_xgb = clip01(p_xgb))
  }

  list(
    models = list(glm = m_glm, #glmnet = m_glmnet,
                  gam = m_gam, rf = m_rf, xgb = m_xgb),
    predict = predict_all
  )
}

head(dt)
base_full <- fit_base_full(dt)

# Estimate stacking weights on the full dataset using in-sample fit for illustration.
# In real applications, weights are usually learned via CV; for mapping, you can:
# - fix weights learned in CV, OR
# - re-estimate weights on full data using CV predictions.
# Here we reuse the CV-trained weights approach indirectly by fitting a convex combination
# against observed prevalence on the full data as a simple demonstration.
stack_weights_convex <- function(P, y, n) {
  wts <- as.numeric(n)
  yprop <- y / n
  p <- ncol(P)
  Dmat <- t(P) %*% (P * wts)
  dvec <- t(P) %*% (wts * yprop)
  Dmat <- Dmat + diag(1e-8, p)

  if (requireNamespace("quadprog", quietly = TRUE)) {
    Amat <- cbind(diag(p), rep(1, p), -rep(1, p))
    bvec <- c(rep(0, p), 1, -1)
    sol <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0)
    w <- sol$solution
    w[w < 0] <- 0
    w <- w / sum(w)
    return(w)
  }
  rep(1 / p, p)
}

P_dt <- base_full$predict(dt)
w_full <- stack_weights_convex(P_dt, y = dt$Y, n = dt$N)

dt[, p_stage1 := as.numeric(P_dt %*% w_full)]
dt[, offset_eta := logit(clip01(p_stage1))]

# ---- Stage-2 spatial meta-learner: binomial with offset + spatial GP smooth ----
# Inputs required:
# - coordinates (lon, lat)
# - stage-1 predictions (as offset on logit scale)
# - outcome (Y,N)
m2_full <- mgcv::gam(
  cbind(Y, N - Y) ~ offset(offset_eta) + s(lon, lat, bs = "gp", k = 120),
  family = stats::binomial(),
  data = dt,
  method = "REML"
)

# ---- Prediction grid ----
# Here we use the raster grid as the prediction grid (pixel centers).
# Inputs required:
# - wall-to-wall raster covariates aligned on a grid
# - if additional covariates exist (e.g., administrative effects), they must be available on the grid

xy <- terra::crds(ras, df = TRUE)
grid <- data.table::data.table(lon = xy[,1], lat = xy[,2])

# Extract raster covariates for the grid (they already exist as ras layers)
vals <- terra::as.data.frame(ras, xy = FALSE, na.rm = FALSE)
grid[, `:=`(
  rs_cov1 = vals$rs_cov1,
  rs_cov2 = vals$rs_cov2,
  rs_cov3 = vals$rs_cov3
)]

# For the non-spatial cluster covariates (x1/x2/x3), a real workflow would use
# covariates available wall-to-wall (or modeled surfaces). Here we create simple proxies:
grid[, x1 := scale(lon)[,1] * 0.3 + rnorm(.N, 0, 0.1)]
grid[, x2 := as.integer(rs_cov2 > 0.4)]
grid[, x3 := scale(lat)[,1] * -0.2 + rnorm(.N, 0, 0.1)]

# Stage-1 predictions on the grid
P_grid <- base_full$predict(grid)
grid[, p_stage1 := as.numeric(P_grid %*% w_full)]
grid[, offset_eta := logit(clip01(p_stage1))]

# ---- Uncertainty: approximate posterior draws from mgcv ----
# We obtain uncertainty from:
# - the fitted spatial meta-learner coefficients and covariance matrix
# - simulating linear predictors, then transforming to probability scale
#
# This yields an approximate posterior predictive SD for prevalence, conditional on stage-1 predictions.

n_sims <- 400
Xp <- mgcv::predict.gam(m2_full, newdata = grid, type = "lpmatrix")
beta_hat <- coef(m2_full)
Vb <- vcov(m2_full)

# Draw coefficients
set.seed(202)
B <- MASS::mvrnorm(n = n_sims, mu = beta_hat, Sigma = Vb)

# Linear predictor draws (includes offset through lpmatrix)
eta_draws <- Xp %*% t(B)
p_draws <- inv_logit(eta_draws)

p_mean <- rowMeans(p_draws)
p_sd <- apply(p_draws, 1, sd)

grid[, `:=`(p_mean = p_mean, p_sd = p_sd)]

# ---- Write rasters and quick plots ----
r_mean <- ras[[1]]; terra::values(r_mean) <- grid$p_mean; names(r_mean) <- "p_mean"
r_sd <- ras[[1]]; terra::values(r_sd) <- grid$p_sd; names(r_sd) <- "p_sd"

terra::writeRaster(r_mean, "outputs/map_mean.tif", overwrite = TRUE)
terra::writeRaster(r_sd, "outputs/map_sd.tif", overwrite = TRUE)

# PNG outputs (simple)
png(here("outputs/map_mean.png"), width = 900, height = 700)
terra::plot(r_mean, main = "Posterior-approximate mean prevalence (Method C)")
dev.off()

png(here("outputs/map_sd.png"), width = 900, height = 700)
terra::plot(r_sd, main = "Posterior-approximate SD of prevalence (Method C)")
dev.off()

message("Mapping complete. See outputs/map_mean.* and outputs/map_sd.*")

# ---- Documentation of required inputs (SPDE-INLA placeholder) ----
# If you implement SPDE-INLA for uncertainty/mapping instead of mgcv:
# Required inputs after fitting stage-1 predictions include:
# - A mesh over the domain (coordinates + tuning parameters controlling triangle sizes)
# - A spatial index linking observation locations to the mesh
# - A likelihood for binomial counts (Y, N) with offset=logit(p_stage1)
# - Priors for the spatial field (range, marginal variance) and (optionally) iid nugget
# - Prediction locations (grid) with projection matrices from mesh->locations
# SPDE-INLA would then return posterior marginals at grid points, enabling mean/SD/interval maps.
