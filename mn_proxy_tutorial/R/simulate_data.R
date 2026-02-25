# =============================================================================
# R/simulate_data.R
#
# Simulation of a realistic micronutrient deficiency survey.
#
# DGP for continuous biomarker Y:
#
#   Y_i = mu_{a(i)}                 (Admin1 random effect)
#         + S(s_{c(i)})             (spatial latent field at cluster location)
#         + u_{c(i)}                (non-spatial cluster random effect)
#         + f(X_i)                  (nonlinear individual/cluster covariates)
#         + epsilon_i               (measurement noise)
#
#   D_i = I(Y_i < deficiency_cutoff)
#
# Spatial field S(s) is sampled from an isotropic MVN with exponential
# covariance K(d) = sigma_s^2 * exp(-d / rho) at cluster centroids.
# =============================================================================

#' Simulate a micronutrient survey dataset.
#'
#' @param n_admin1          Number of Admin1 regions (default 12).
#' @param clusters_per_admin1 Clusters per Admin1 (default 8).
#' @param n_per_cluster     Individuals per cluster (default 15).
#' @param seed              RNG seed for reproducibility (default 42).
#' @param deficiency_cutoff Threshold on Y below which an individual is
#'                          classified as deficient (default -0.5).
#' @param sigma_admin1      SD of Admin1 random effects (default 0.4).
#' @param sigma_cluster     SD of cluster random effects (default 0.3).
#' @param sigma_spatial     SD of spatial latent field (default 0.5).
#' @param spatial_range     Range parameter rho for exponential covariance;
#'                          units match the [0,1]^2 coordinate space (default 0.25).
#' @param sigma_epsilon     SD of individual-level noise (default 0.6).
#' @return A data frame with one row per individual.
simulate_mn_survey <- function(n_admin1           = 12,
                                clusters_per_admin1 = 8,
                                n_per_cluster       = 15,
                                seed                = 42,
                                deficiency_cutoff   = -0.5,
                                sigma_admin1        = 0.4,
                                sigma_cluster       = 0.3,
                                sigma_spatial       = 0.5,
                                spatial_range       = 0.25,
                                sigma_epsilon       = 0.6) {

  set.seed(seed)
  n_clusters <- n_admin1 * clusters_per_admin1
  N          <- n_clusters * n_per_cluster

  # ── Admin1 layout ───────────────────────────────────────────────────────────
  admin1_id   <- rep(seq_len(n_admin1), each = clusters_per_admin1)
  admin1_name <- paste0("Region_", formatC(admin1_id, width = 2, flag = "0"))
  cluster_id  <- seq_len(n_clusters)

  # Admin1 random effects (centred near 0 so mean Y is interpretable)
  mu_admin1 <- rnorm(n_admin1, mean = 0, sd = sigma_admin1)

  # ── Cluster-level spatial coordinates ──────────────────────────────────────
  # Place cluster centroids on a quasi-regular jittered grid within [0,1]^2,
  # grouped by Admin1 "territory" to give geographic coherence.
  admin1_centres <- data.frame(
    cx = seq(0.1, 0.9, length.out = n_admin1),
    cy = rep(c(0.25, 0.5, 0.75), length.out = n_admin1)
  )
  x_coord <- numeric(n_clusters)
  y_coord <- numeric(n_clusters)
  for (a in seq_len(n_admin1)) {
    idx <- which(admin1_id == a)
    x_coord[idx] <- admin1_centres$cx[a] + runif(length(idx), -0.07, 0.07)
    y_coord[idx] <- admin1_centres$cy[a] + runif(length(idx), -0.07, 0.07)
  }
  x_coord <- pmax(0.01, pmin(0.99, x_coord))
  y_coord <- pmax(0.01, pmin(0.99, y_coord))

  # ── Spatial latent field at cluster locations ───────────────────────────────
  coords_mat <- cbind(x_coord, y_coord)
  dists      <- as.matrix(dist(coords_mat))
  Sigma_s    <- exp_cov(dists, sigma2 = sigma_spatial^2, rho = spatial_range)
  # Add small nugget for numerical stability
  diag(Sigma_s) <- diag(Sigma_s) + 1e-6
  spatial_field <- as.numeric(MASS::mvrnorm(1, mu = rep(0, n_clusters),
                                             Sigma = Sigma_s))

  # ── Cluster-level random effects ────────────────────────────────────────────
  cluster_re <- rnorm(n_clusters, 0, sigma_cluster)

  # ── Cluster-level covariates (proxy "raster" variables) ────────────────────
  # WASH domain
  wash_water      <- plogis(rnorm(n_clusters, mean = 0.5 - x_coord, sd = 0.5))
  wash_sanitation <- plogis(rnorm(n_clusters, mean = 0.4 + y_coord * 0.6, sd = 0.4))

  # Climate domain (smooth gradient + noise)
  clim_rainfall <- 400 + 600 * x_coord + rnorm(n_clusters, 0, 80)
  clim_temp     <-  22 + 8   * y_coord  + rnorm(n_clusters, 0, 1.5)

  # Disease ecology (Admin1-level HIV + cluster-level malaria)
  hiv_by_admin1    <- rbeta(n_admin1, 2, 10)
  dis_hiv          <- hiv_by_admin1[admin1_id]
  dis_malaria      <- plogis(rnorm(n_clusters, mean = -0.5 + 2 * y_coord, sd = 0.8))

  # SES domain (cluster-level wealth index — used as cluster mean)
  ses_wealth_clust <- rnorm(n_clusters, mean = 0.4 * wash_water - 0.3, sd = 0.3)

  # ── Expand cluster data to individual level ─────────────────────────────────
  cluster_idx <- rep(seq_len(n_clusters), each = n_per_cluster)
  indiv_id    <- seq_len(N)

  # Individual-level covariates
  indiv_age     <- runif(N, 15, 49)        # women of reproductive age
  indiv_sex     <- rbinom(N, 1, 0.5)       # 0 = female (majority of target pop)
  indiv_suppl   <- rbinom(N, 1,
                           plogis(-1 + 0.5 * wash_water[cluster_idx]))

  # Individual SES (cluster mean + individual noise)
  ses_wealth    <- ses_wealth_clust[cluster_idx] + rnorm(N, 0, 0.3)
  ses_educ      <- pmax(0, round(8 + 2 * ses_wealth + rnorm(N, 0, 2)))

  # ── Nonlinear covariate effects on Y ───────────────────────────────────────
  # f(X) captures realistic associations:
  #   - age: U-shaped deficiency (young and old more deficient)
  #   - wealth: protective (higher wealth → higher biomarker)
  #   - rainfall: moderate rainfall associated with better nutrition
  #   - malaria: strong negative effect (inflammation depletes ferritin/RBP)
  #   - supplementation: protective

  age_z      <- (indiv_age - 30) / 10
  f_X <- (
    -0.5 * dis_malaria[cluster_idx]                  +   # main driver
    0.4  * ses_wealth                                 +   # SES
    0.3  * wash_water[cluster_idx]                    +   # WASH
    -0.2 * age_z^2                                    +   # U-shape age
    0.3  * indiv_suppl                                +   # supplementation
    0.15 * (clim_rainfall[cluster_idx] - 700) / 300   +   # moderate rainfall
    -0.2 * clim_temp[cluster_idx] / 30                    # temperature
  )

  # ── Compose Y ──────────────────────────────────────────────────────────────
  Y <- (mu_admin1[admin1_id[cluster_idx]]
      + spatial_field[cluster_idx]
      + cluster_re[cluster_idx]
      + f_X
      + rnorm(N, 0, sigma_epsilon))

  # ── Deficiency indicator ────────────────────────────────────────────────────
  D <- as.integer(Y < deficiency_cutoff)

  # ── Survey weights (unequal inclusion probability) ─────────────────────────
  # Clusters in rural areas (low wash_water) have lower sampling probability.
  # Weights are the inverse of the normalised inclusion probability.
  incl_prob_clust <- plogis(0.5 + wash_water)   # 0.5–0.7 range
  weight_raw      <- 1 / incl_prob_clust[cluster_idx]
  svy_weight      <- weight_raw / mean(weight_raw)  # normalised to mean 1

  # ── Assemble output data frame ──────────────────────────────────────────────
  df <- data.frame(
    # Identifiers
    indiv_id   = indiv_id,
    cluster_id = cluster_idx,
    admin1_id  = admin1_id[cluster_idx],
    admin1_name = admin1_name[cluster_idx],
    # Coordinates (cluster centroid, repeated for individuals)
    x_coord    = x_coord[cluster_idx],
    y_coord    = y_coord[cluster_idx],
    # WASH covariates (cluster-level)
    wash_water      = wash_water[cluster_idx],
    wash_sanitation = wash_sanitation[cluster_idx],
    # SES
    ses_wealth  = ses_wealth,
    ses_educ    = ses_educ,
    # Climate (cluster-level)
    clim_rainfall = clim_rainfall[cluster_idx],
    clim_temp     = clim_temp[cluster_idx],
    # Disease ecology
    dis_malaria = dis_malaria[cluster_idx],
    dis_hiv     = dis_hiv[cluster_idx],
    # Individual-level
    indiv_age   = indiv_age,
    indiv_sex   = indiv_sex,
    indiv_suppl = indiv_suppl,
    # Outcomes
    Y           = Y,
    D           = D,
    svy_weight  = svy_weight
  )

  # Store simulation metadata as attributes for reference
  attr(df, "deficiency_cutoff") <- deficiency_cutoff
  attr(df, "true_national_prev") <- mean(D)
  df
}


# ── True Admin1 prevalence ────────────────────────────────────────────────────

#' Compute the true (population) Admin1 prevalence from a simulated dataset.
#'
#' @param df          Simulated data frame from simulate_mn_survey.
#' @param admin1_col  Name of Admin1 column.
#' @return Named numeric vector: Admin1 name → true prevalence.
true_admin1_prevalence <- function(df, admin1_col = "admin1_name") {
  tapply(df$D, df[[admin1_col]], mean)
}


# ── Covariate grid simulation ─────────────────────────────────────────────────

#' Simulate a wall-to-wall covariate grid over [0,1]^2.
#'
#' Creates a data frame of grid points with the same covariate names as the
#' survey data, drawn from the same generative model (without individual-level
#' covariates or outcomes).
#'
#' @param resolution  Grid spacing in each dimension (default 0.05 → 400 cells).
#' @param survey_df   Original survey data frame (used to infer covariate means).
#' @param seed        RNG seed.
#' @return Data frame with one row per grid cell.
simulate_covariate_grid <- function(resolution = 0.05,
                                     survey_df   = NULL,
                                     seed        = 99) {
  set.seed(seed)
  grid_x <- seq(resolution / 2, 1 - resolution / 2, by = resolution)
  grid_y <- seq(resolution / 2, 1 - resolution / 2, by = resolution)
  grid   <- expand.grid(x_coord = grid_x, y_coord = grid_y)
  n_grid <- nrow(grid)

  # Covariates: same functional form as DGP
  wash_water      <- plogis(rnorm(n_grid, mean = 0.5 - grid$x_coord, sd = 0.5))
  wash_sanitation <- plogis(rnorm(n_grid, mean = 0.4 + grid$y_coord * 0.6, sd = 0.4))
  clim_rainfall   <- 400 + 600 * grid$x_coord + rnorm(n_grid, 0, 80)
  clim_temp       <-  22 + 8   * grid$y_coord  + rnorm(n_grid, 0, 1.5)
  dis_malaria     <- plogis(rnorm(n_grid, mean = -0.5 + 2 * grid$y_coord, sd = 0.8))
  dis_hiv         <- rbeta(n_grid, 2, 10)
  ses_wealth      <- rnorm(n_grid, mean = 0.4 * wash_water - 0.3, sd = 0.4)
  ses_educ        <- pmax(0, round(8 + 2 * ses_wealth + rnorm(n_grid, 0, 2)))

  # Use population-mean individual-level values for grid
  data.frame(
    x_coord         = grid$x_coord,
    y_coord         = grid$y_coord,
    wash_water      = wash_water,
    wash_sanitation = wash_sanitation,
    clim_rainfall   = clim_rainfall,
    clim_temp       = clim_temp,
    dis_malaria     = dis_malaria,
    dis_hiv         = dis_hiv,
    ses_wealth      = ses_wealth,
    ses_educ        = ses_educ,
    # Use mean individual values for age/sex/supplement
    indiv_age       = 30,
    indiv_sex       = 0.5,
    indiv_suppl     = 0.3
  )
}
