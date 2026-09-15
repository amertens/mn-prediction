# 01_simulate_data.R
# Simulate cluster-level survey prevalence data + wall-to-wall raster covariates.
#
# Key design features:
# - Cluster locations (lon/lat) over a rectangular domain
# - Binomial outcomes: Y successes out of N sampled per cluster
# - Non-spatial covariates (measured at clusters)
# - Raster-like covariates available wall-to-wall (terra SpatRaster)
# - A latent spatially correlated surface (Gaussian process) that induces residual
#   spatial autocorrelation after accounting for covariates.
#
# Outputs:
# - data/sim/cluster_data.csv
# - data/sim/domain_covariates.tif  (multi-layer raster)
# - data/sim/true_surface.tif       (latent GP surface on grid; for diagnostics only)

library(here)
stopifnot(requireNamespace("terra", quietly = TRUE))
stopifnot(requireNamespace("data.table", quietly = TRUE))

# ---- Domain ----
# Approximate subnational to regional scale (≈ 500–700 km)
domain <- list(
  xmin = 0,  xmax = 6,   # degrees longitude
  ymin = 0,  ymax = 6    # degrees latitude
)

# ---- Utility: pairwise distance ----
pairwise_dist <- function(x, y) {
  # x, y are numeric vectors of length n (coordinates)
  # returns an n x n Euclidean distance matrix
  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  sqrt(dx^2 + dy^2)
}

# ---- Simulate latent spatial GP at cluster locations ----
simulate_gp <- function(lon, lat, sigma = 0.8, range = 20, nugget = 1e-6) {
  # Exponential covariance: Cov(h) = sigma^2 * exp(-h / range)
  D <- pairwise_dist(lon, lat)
  K <- sigma^2 * exp(-D / range)
  diag(K) <- diag(K) + nugget
  # Cholesky for simulation
  L <- chol(K)
  as.numeric(t(L) %*% rnorm(length(lon)))
}

# ---- Create wall-to-wall raster covariates ----
make_rasters <- function(domain, res = 1) {
  r <- terra::rast(
    xmin = domain$xmin, xmax = domain$xmax,
    ymin = domain$ymin, ymax = domain$ymax,
    resolution = res,
    crs = NA
  )

  # Create synthetic spatial covariates with different structures
  xy <- terra::crds(r, df = TRUE)
  lon <- xy[,1]; lat <- xy[,2]

  # Covariate 1: smooth gradient + sinusoid
  z1 <- scale(lon) + 0.6 * sin(lat / 8)

  # Covariate 2: radial bump
  z2 <- exp(-((lon - 70)^2 + (lat - 30)^2) / (2 * 18^2))

  # Covariate 3: patchy field (smoothed noise)
  set.seed(123)  # stable raster realization across runs
  z3 <- rnorm(nrow(xy))
  # Smooth by averaging in a moving window
  r3 <- r
  terra::values(r3) <- z3
  r3s <- terra::focal(r3, w = matrix(1, 7, 7), fun = mean, na.rm = TRUE, fillvalue = NA)
  z3s <- as.numeric(terra::values(r3s))
  z3s <- (z3s - mean(z3s, na.rm = TRUE)) / sd(z3s, na.rm = TRUE)

  # Assign
  r1 <- r; terra::values(r1) <- as.numeric(z1)
  r2 <- r; terra::values(r2) <- as.numeric(z2)
  r3 <- r; terra::values(r3) <- z3s

  names(r1) <- "rs_cov1"
  names(r2) <- "rs_cov2"
  names(r3) <- "rs_cov3"

  c(r1, r2, r3)
}

# ---- Simulate cluster survey data ----
simulate_cluster_data <- function(
    n_clusters = 300,
    domain = domain,
    min_n = 20, max_n = 80,
    gp_sigma = 0.9, gp_range = 18,
    beta = c(intercept = -1.2, x1 = 0.9, x2 = -0.6, x3 = 0.7, z1 = 0.8, z2 = -1.0, z3 = 0.5),
    seed = 999
) {
  set.seed(seed)

  dt <- data.table::data.table(
    cluster_id = 1:n_clusters,
    lon = runif(n_clusters, domain$xmin, domain$xmax),
    lat = runif(n_clusters, domain$ymin, domain$ymax)
  )

  # Non-spatial covariates (measured at clusters)
  dt[, x1 := rnorm(.N)]
  dt[, x2 := rbinom(.N, 1, 0.45)]
  dt[, x3 := runif(.N, -1, 1)]

  # Raster covariates extracted at cluster locations
  ras <- make_rasters(domain, res = 1)
  pts <- terra::vect(dt[, .(lon, lat)], geom = c("lon", "lat"), crs = NA)
  rc <- terra::extract(ras, pts)
  dt[, c("rs_cov1", "rs_cov2", "rs_cov3") := rc[, -1]]

  # Cluster sample sizes
  dt[, N := sample(min_n:max_n, .N, replace = TRUE)]

  # Latent spatial surface (unobserved)
  dt[, gp := simulate_gp(lon, lat, sigma = gp_sigma, range = gp_range)]

  # Linear predictor for prevalence
  eta <- beta["intercept"] +
    beta["x1"] * dt$x1 +
    beta["x2"] * dt$x2 +
    beta["x3"] * dt$x3 +
    beta["z1"] * dt$rs_cov1 +
    beta["z2"] * dt$rs_cov2 +
    beta["z3"] * dt$rs_cov3 +
    dt$gp

  p <- plogis(eta)
  dt[, p_true := p]

  # Observed binomial outcomes
  dt[, Y := rbinom(.N, size = N, prob = p_true)]
  dt[, prev := Y / N]

  dt[]
}

# ---- Run simulation and save ----
rasters <- make_rasters(domain, res = 1)
terra::writeRaster(rasters, filename = here("data/sim/domain_covariates.tif"), overwrite = TRUE)

dt <- simulate_cluster_data(  n_clusters = 300,
                              domain = domain,
                              min_n = 20, max_n = 80,
                              gp_sigma = 0.9, gp_range = 18,
                              beta = c(intercept = -1.2, x1 = 0.9, x2 = -0.6, x3 = 0.7, z1 = 0.8, z2 = -1.0, z3 = 0.5),
                              seed = 999)

data.table::fwrite(dt, here("data/sim/cluster_data.csv"))

# Save latent GP surface on the same grid for diagnostics (not used for fitting)
# For mapping exercises, the latent field is unknown in real problems; here it is for validation only.
grid_xy <- terra::crds(rasters, df = TRUE)
gp_grid <- simulate_gp(grid_xy[,1], grid_xy[,2], sigma = 0.9, range = 18)
gp_r <- rasters[[1]]
terra::values(gp_r) <- gp_grid
names(gp_r) <- "gp_true"
terra::writeRaster(gp_r, filename = here("data/sim/true_surface.tif"), overwrite = TRUE)

message("Simulation complete: data/sim/cluster_data.csv and raster covariates written.")
