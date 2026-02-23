# -------------------- Packages --------------------
rm(list=ls())
library(sf)
library(dplyr)
library(SuperLearner)
library(ranger)      # SL.ranger
library(gam)         # SL.gam
library(xgboost)     # SL.xgboost
library(blockCV)
library(INLA)
library(rnaturalearth)
library(ggplot2)
library(here)
library(future)
library(future.apply)

options(future.globals.maxSize = 5 * 1024^3)

prev_gps <- readRDS(here("data","IPD","Ghana","GMS_cluster_prevalence.rds"))
prev <- readRDS(here("data","IPD","Ghana","GMS_cluster_prevalence_points.rds"))
head(prev)

#"C:\Users\andre\Downloads\TRMM_Ghana_2017.tif"

prev %>% group_by(outcome) %>% summarise(mean(prev), sd(prev))

# temp demo covars (replace with your cluster-level summaries when ready)
prev_gps <- prev_gps %>% select(gw_cnum, outcome, prev, lat, lon) #%>% mutate(sex = 1, age = 55)
prev <- prev %>% select(gw_cnum, outcome, prev, geometry) #%>% mutate(sex = 1, age = 55)
df <- prev %>% filter(outcome == "Child iron deficiency (adj)")

# -------------------- 1) Prepare data --------------------
dat_sf <- df %>%
  st_as_sf() %>%
  st_make_valid() %>%
  st_transform(4326)

stopifnot(inherits(dat_sf$geometry, "sfc_POINT"))
stopifnot(all(dat_sf$prev >= 0 & dat_sf$prev <= 1))
dat_sf <- dat_sf %>% distinct(gw_cnum, .keep_all = TRUE)

coords_ll <- st_coordinates(dat_sf)
dat_sf <- dat_sf %>% mutate(lon = coords_ll[,1], lat = coords_ll[,2])

pred_cols <- c("lon", "lat")
X <- dat_sf %>% st_drop_geometry() %>% select(any_of(pred_cols)) %>% as.data.frame()
Y <- dat_sf$prev

eps   <- 1e-4
Y_eta <- qlogis(pmin(pmax(Y, eps), 1 - eps))

# -------------------- 2) Spatial CV folds (blocking) --------------------
dat_utm      <- st_transform(dat_sf, 32630)  # UTM 30N
dat_utm$resp <- dat_sf$prev                  # harmless (we won't pass 'column')

set.seed(2025)
cv <- blockCV::cv_spatial(
  x         = dat_utm,
  size      = 50000,   # ~50 km
  k         = 5,
  selection = "random",
  iteration = 200,
  hexagon   = FALSE
)

# robustly extract fold IDs across blockCV versions
fold_id <- if (!is.null(cv$foldID)) cv$foldID else cv$folds_ids
stopifnot(!is.null(fold_id))
validRows <- lapply(seq_len(max(fold_id)), function(i) which(fold_id == i))

# -------------------- 3) Super Learner on logit(prevalence) --------------------
#SL.lib <- c("SL.glm", "SL.gam", "SL.ranger", "SL.xgboost")
#SL.lib <- c("SL.glm", "SL.glmnet")
SL.lib <- c("SL.glm", "SL.ranger")

set.seed(2025)
fit_sl <- SuperLearner(
  Y = Y_eta,
  X = X,
  family     = gaussian(),     # logit(y) modeled as Gaussian
  SL.library = SL.lib,
  method     = "method.NNLS",
  cvControl  = list(V = length(validRows), validRows = validRows, shuffle = FALSE)
)

oof_eta   <- fit_sl$SL.predict
resid_eta <- as.numeric(Y_eta - oof_eta)

# -------------------- 4) Spatial residual GP via INLA/SPDE --------------------
coords_utm <- st_coordinates(dat_utm)
mesh <- inla.mesh.2d(
  loc      = coords_utm,
  max.edge = c(20e3, 60e3),
  cutoff   = 10e3,
  offset   = c(50e3, 100e3)
)
spde <- inla.spde2.pcmatern(
  mesh,
  alpha       = 2,
  prior.range = c(100e3, 0.5),   # P(range < 100 km) = 0.5
  prior.sigma = c(1.0, 0.01)     # P(sigma > 1) = 0.01
)

A_obs   <- inla.spde.make.A(mesh, loc = coords_utm)
s.index <- inla.spde.make.index("s", n.spde = spde$n.spde)

n_obs <- length(resid_eta)
X_obs <- data.frame(intercept = rep(1, n_obs))

stk_obs <- inla.stack(
  data    = list(y = resid_eta),
  A       = list(A_obs, 1),
  effects = list(s = s.index, X_obs),
  tag     = "obs"
)

# -------------------- 5) Build a 5-km prediction grid over Ghana --------------------
ghana_ll  <- rnaturalearth::ne_countries(scale = "medium", country = "Ghana", returnclass = "sf") %>% st_transform(4326)
ghana_utm <- st_transform(ghana_ll, 32630)

grid_centers <- st_make_grid(ghana_utm, cellsize = 5000, what = "centers") %>% st_sf()
grid_sf      <- st_intersection(grid_centers, ghana_utm)  # keep inside Ghana

# predictors for the grid (lon/lat + any area-level covariates you later add)
grid_ll        <- st_transform(grid_sf, 4326)
grid_coords_ll <- st_coordinates(grid_ll)
X_grid <- data.frame(
  lon = grid_coords_ll[,1],
  lat = grid_coords_ll[,2]
)

# ensemble mean on grid (logit scale)
eta_sl_grid <- as.numeric(predict(fit_sl, newdata = X_grid, onlySL = TRUE)$pred)

# spatial residuals on grid
A_pred  <- inla.spde.make.A(mesh, loc = st_coordinates(grid_sf))
n_pred  <- nrow(grid_sf)
X_pred  <- data.frame(intercept = rep(1, n_pred))

stk_pred <- inla.stack(
  data    = list(y = NA),
  A       = list(A_pred, 1),
  effects = list(s = s.index, X_pred),
  tag     = "pred"
)

stk_all <- inla.stack(stk_obs, stk_pred)

fit_gp <- inla(
  y ~ 1 + f(s, model = spde),
  data = inla.stack.data(stk_all),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_all), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

idx_pred        <- inla.stack.index(stk_all, tag = "pred")$data
eta_resid_grid  <- as.numeric(fit_gp$summary.fitted.values[idx_pred, "mean"])

eta_grid <- eta_sl_grid + eta_resid_grid
p_grid   <- plogis(eta_grid)
pred_sf  <- st_set_geometry(grid_ll, "geometry") %>% mutate(pred_prev = p_grid)

# -------------------- 6) Quick map --------------------
ggplot() +
  geom_sf(data = ghana_ll, fill = NA, linewidth = 0.3) +
  geom_sf(data = pred_sf, aes(color = pred_prev), size = 3) +
  scale_color_viridis_c(name = "Predicted prevalence", limits = c(0, 1)) +
  coord_sf(expand = FALSE) +
  theme_minimal()

# -------------------- 7) Sanity checks --------------------
stopifnot(n_obs == nrow(A_obs),
          n_obs == nrow(X_obs),
          n_pred == nrow(A_pred),
          n_pred == nrow(X_pred),
          length(resid_eta) == n_obs)



# --- assumes you already ran your working script up to building:
# dat_sf (clusters; WGS84 with lon/lat and prev, sex, age)
# dat_utm (same in EPSG:32630)
# mesh, spde, ghana_ll, ghana_utm, grid_sf, grid_ll, X_grid
# SL library definition (SL.lib)
# -------------------------------------------------------------
# -------------------- 0) Setup --------------------

# Toggle parallel vs sequential in one place while debugging
USE_PARALLEL <- TRUE
if (USE_PARALLEL) {
  plan(multisession, workers = max(1, parallel::detectCores()/2))
} else {
  plan(sequential)
}

set.seed(4242)
B <- 10  # start small when debugging

pkgs_needed <- c("SuperLearner","ranger","xgboost","gam","blockCV","sf","INLA")

# Make sure these globals exist before we spawn workers
stopifnot(
  exists("dat_sf"), exists("mesh"), exists("spde"),
  exists("grid_sf"), exists("X_grid"), exists("SL.lib")
)

# -------------------- 1) One replicate function --------------------
boot_fit_once <- function(dat_sf_boot) {
  # ----- Ensure a clean worker environment -----
  invisible(lapply(pkgs_needed, require, character.only = TRUE))
  grDevices::pdf(NULL)                               # avoid device warnings
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  INLA::inla.setOption(num.threads = "1:1")          # INLA single-thread per worker
  dummy <- SuperLearner::All                         # ensure SL screener 'All' is visible

  # ----- 1) Spatial CV folds for this bootstrap sample -----
  dat_utm_b <- sf::st_transform(dat_sf_boot, 32630)
  cv_b <- blockCV::cv_spatial(
    x         = dat_utm_b,
    size      = 50000,
    k         = 5,
    selection = "random",
    iteration = 10, #less for speed
    hexagon   = FALSE,
    verbose   = FALSE
  )
  fold_id_b <- if (!is.null(cv_b$foldID)) cv_b$foldID else cv_b$folds_ids
  if (is.null(fold_id_b)) stop("blockCV returned no fold IDs.")
  validRows_b <- lapply(seq_len(max(fold_id_b)), function(i) which(fold_id_b == i))

  # ----- 2) SuperLearner on logit(prevalence) -----
  X_b <- dat_sf_boot |>
    sf::st_drop_geometry() |>
    dplyr::select(lon, lat) |>
    as.data.frame()

  eps <- 1e-4
  Y_b_eta <- qlogis(pmin(pmax(dat_sf_boot$prev, eps), 1 - eps))

  fit_sl_b <- SuperLearner::SuperLearner(
    Y = Y_b_eta, X = X_b,
    family     = gaussian(),
    SL.library = SL.lib,
    method     = "method.NNLS",
    cvControl  = list(V = length(validRows_b), validRows = validRows_b, shuffle = FALSE)
  )

  oof_eta_b   <- fit_sl_b$SL.predict
  resid_eta_b <- as.numeric(Y_b_eta - oof_eta_b)

  # ----- 3) Residual GP via INLA/SPDE -----
  coords_utm_b <- sf::st_coordinates(dat_utm_b)
  A_obs_b      <- INLA::inla.spde.make.A(mesh, loc = coords_utm_b)
  s.index      <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)

  # stacks with explicit intercept column named "intercept"
  X_obs_b <- data.frame(intercept = rep(1, length(resid_eta_b)))
  stk_obs_b <- INLA::inla.stack(
    data    = list(y = resid_eta_b),
    A       = list(A_obs_b, 1),
    effects = list(s = s.index, X_obs_b),
    tag     = "obs"
  )

  A_pred  <- INLA::inla.spde.make.A(mesh, loc = sf::st_coordinates(grid_sf))
  X_pred  <- data.frame(intercept = rep(1, nrow(grid_sf)))
  stk_pred_b <- INLA::inla.stack(
    data    = list(y = NA),
    A       = list(A_pred, 1),
    effects = list(s = s.index, X_pred),
    tag     = "pred"
  )

  stk_all_b <- INLA::inla.stack(stk_obs_b, stk_pred_b)

  fit_gp_b <- INLA::inla(
    y ~ 0 + intercept + f(s, model = spde),   # << remove default intercept, use our column
    data = INLA::inla.stack.data(stk_all_b),
    family = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stk_all_b), compute = TRUE),
    control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE)
  )

  idx_pred_b       <- INLA::inla.stack.index(stk_all_b, tag = "pred")$data
  eta_resid_grid_b <- as.numeric(fit_gp_b$summary.fitted.values[idx_pred_b, "mean"])

  # ----- 4) Combine SL mean + GP residuals on grid -----
  eta_sl_grid_b <- as.numeric(predict(fit_sl_b, newdata = X_grid, onlySL = TRUE)$pred)
  p_grid_b      <- plogis(eta_sl_grid_b + eta_resid_grid_b)

  return(p_grid_b)
}

# -------------------- 2) Bootstrap loop with robust error handling --------------------
library(future.apply)

USE_PARALLEL <- TRUE
if (USE_PARALLEL) {
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
} else {
  plan(sequential)
}

set.seed(4242)
B <- 100

boot_draws <- future_lapply(
  seq_len(B),
  function(b) {
    # ---- worker hygiene ----
    pkgs_needed <- c("SuperLearner","ranger","xgboost","gam","blockCV","sf","INLA")
    invisible(lapply(pkgs_needed, require, character.only = TRUE))
    grDevices::pdf(NULL); on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    INLA::inla.setOption(num.threads = "1:1")
    dummy <- SuperLearner::All

    # ---- draw bootstrap sample (compute N *inside*) ----
    N_b  <- nrow(dat_sf)                 # size of the original cluster set
    idx  <- sample.int(N_b, size = N_b, replace = TRUE)
    dat_b <- dat_sf[idx, ]               # bootstrap sample (sf)

    # ---- spatial folds for this replicate ----
    dat_utm_b <- sf::st_transform(dat_b, 32630)
    cv_b <- blockCV::cv_spatial(
      x = dat_utm_b, size = 50000, k = 5,
      selection = "random",
      iteration = 10,
      hexagon = FALSE
    )
    fold_id_b <- if (!is.null(cv_b$foldID)) cv_b$foldID else cv_b$folds_ids
    validRows_b <- lapply(seq_len(max(fold_id_b)), function(i) which(fold_id_b == i))

    # ---- SuperLearner on logit(prev) ----
    X_b <- dat_b |>
      sf::st_drop_geometry() |>
      dplyr::select(lon, lat) |>
      as.data.frame()

    eps <- 1e-4
    Y_eta_b <- qlogis(pmin(pmax(dat_b$prev, eps), 1 - eps))

    fit_sl_b <- SuperLearner::SuperLearner(
      Y = Y_eta_b, X = X_b,
      family     = gaussian(),
      SL.library = SL.lib,
      method     = "method.NNLS",
      cvControl  = list(V = length(validRows_b), validRows = validRows_b, shuffle = FALSE)
    )

    oof_eta_b   <- fit_sl_b$SL.predict
    resid_eta_b <- as.numeric(Y_eta_b - oof_eta_b)

    # ---- INLA residual GP (with explicit intercept, no duplicate) ----
    coords_utm_b <- sf::st_coordinates(dat_utm_b)
    A_obs_b      <- INLA::inla.spde.make.A(mesh, loc = coords_utm_b)
    s.index      <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)

    X_obs_b <- data.frame(intercept = rep(1, length(resid_eta_b)))
    stk_obs_b <- INLA::inla.stack(
      data    = list(y = resid_eta_b),
      A       = list(A_obs_b, 1),
      effects = list(s = s.index, X_obs_b),
      tag     = "obs"
    )

    A_pred  <- INLA::inla.spde.make.A(mesh, loc = sf::st_coordinates(grid_sf))
    X_pred  <- data.frame(intercept = rep(1, nrow(grid_sf)))
    stk_pred_b <- INLA::inla.stack(
      data    = list(y = NA),
      A       = list(A_pred, 1),
      effects = list(s = s.index, X_pred),
      tag     = "pred"
    )

    stk_all_b <- INLA::inla.stack(stk_obs_b, stk_pred_b)

    fit_gp_b <- INLA::inla(
      y ~ 0 + intercept + f(s, model = spde),
      data = INLA::inla.stack.data(stk_all_b),
      family = "gaussian",
      control.predictor = list(A = INLA::inla.stack.A(stk_all_b), compute = TRUE),
      control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE)
    )

    idx_pred <- INLA::inla.stack.index(stk_all_b, tag = "pred")$data
    eta_resid_grid_b <- as.numeric(fit_gp_b$summary.fitted.values[idx_pred, "mean"])

    # ---- grid predictions: SL mean + GP residuals ----
    eta_sl_grid_b <- as.numeric(predict(fit_sl_b, newdata = X_grid, onlySL = TRUE)$pred)
    p_grid_b      <- plogis(eta_sl_grid_b + eta_resid_grid_b)

    p_grid_b
  },
  future.seed     = TRUE,
  future.packages = pkgs_needed,
  future.globals  = c("dat_sf","SL.lib","mesh","spde","grid_sf","X_grid")
)

# Summaries (unchanged)
pred_mat <- do.call(rbind, boot_draws)

saveRDS(pred_mat, here("data","IPD","Ghana","GMS_cluster_prevalence_bootdraws.rds"))
