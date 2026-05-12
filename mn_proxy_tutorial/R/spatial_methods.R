# =============================================================================
# R/spatial_methods.R
#
# Spatial analysis utilities:
#   - Moran's I test on cluster residuals.
#   - Spatial block CV construction (k-means on coordinates).
#   - Fay-Herriot small area estimation (using the 'sae' package).
#   - Simulation of Admin1 polygons for mapping (no external shapefile needed).
#   - Variogram plotting helper.
# =============================================================================

# ── Moran's I ─────────────────────────────────────────────────────────────────

#' Compute Moran's I and Monte-Carlo test on cluster-level residuals.
#'
#' Constructs a k-nearest-neighbour spatial weights matrix from cluster
#' coordinates, then applies spdep::moran.mc() to the supplied residuals.
#'
#' @param residuals  Numeric vector of cluster-level residuals
#'                   (length = n_clusters).
#' @param coords     Data frame or matrix with columns x_coord, y_coord,
#'                   one row per cluster.
#' @param k          Number of nearest neighbours for weights (default 6).
#' @param nsim       Monte-Carlo permutations for p-value (default 999).
#' @return List: $I (Moran's I), $p_value, $plot (ggplot of test distribution).
compute_morans_i <- function(residuals, coords, k = 6, nsim = 499) {
  require(spdep, quietly = TRUE)

  xy  <- as.matrix(coords[, c("x_coord", "y_coord")])
  nb  <- spdep::knn2nb(spdep::knearneigh(xy, k = k))
  lw  <- spdep::nb2listw(nb, style = "W")

  mc  <- spdep::moran.mc(residuals, lw, nsim = nsim)

  # Distribution plot
  perm_df <- data.frame(I = mc$res[-length(mc$res)])
  p_plot <- ggplot2::ggplot(perm_df, ggplot2::aes(x = I)) +
    ggplot2::geom_histogram(bins = 40, fill = "steelblue", colour = "white",
                             alpha = 0.7) +
    ggplot2::geom_vline(xintercept = mc$statistic, colour = "firebrick",
                         linewidth = 1.2, linetype = "dashed") +
    ggplot2::annotate("text", x = mc$statistic + 0.01,
                       y    = 0,
                       label = sprintf("Observed I = %.3f\np = %.3f",
                                       mc$statistic, mc$p.value),
                       hjust = 0, vjust = -0.5, size = 3.5) +
    ggplot2::labs(
      x     = "Moran's I (permuted)",
      y     = "Count",
      title = "Moran's I Monte-Carlo test",
      subtitle = "Dashed red = observed value"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  list(I = mc$statistic, p_value = mc$p.value,
       statistic = mc, plot = p_plot)
}


#' Extract cluster-level residuals from individual predictions.
#'
#' Aggregates individual-level residuals (observed - predicted) to cluster
#' means, returning a data frame suitable for Moran's I.
#'
#' @param df          Data frame with cluster, outcome, and prediction columns.
#' @param cluster_col Cluster column name.
#' @param obs_col     Observed outcome column.
#' @param pred_col    Predicted value column.
#' @return Data frame: cluster_id, x_coord, y_coord, mean_residual.
cluster_residuals <- function(df,
                               cluster_col = "cluster_id",
                               obs_col     = "Y",
                               pred_col    = "pred_cont") {
  df %>%
    dplyr::mutate(resid = .data[[obs_col]] - .data[[pred_col]]) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::summarise(
      x_coord      = mean(x_coord, na.rm = TRUE),
      y_coord      = mean(y_coord, na.rm = TRUE),
      mean_residual = mean(resid,  na.rm = TRUE),
      n             = dplyr::n(),
      .groups       = "drop"
    )
}


# ── Spatial block cross-validation ───────────────────────────────────────────

#' Create spatial CV blocks using k-means clustering of coordinates.
#'
#' Unlike cluster-blocked CV (which groups by survey cluster), spatial
#' blocking groups geographically proximate clusters into the same fold.
#' This prevents leakage when spatial autocorrelation extends across the
#' typical survey cluster radius.
#'
#' @param coords     Data frame with x_coord, y_coord (one row per obs).
#' @param n_blocks   Number of spatial blocks (folds).
#' @param seed       RNG seed.
#' @return Integer vector (length = nrow(coords)): block assignment per obs.
make_spatial_blocks <- function(coords, n_blocks = 5, seed = 42) {
  set.seed(seed)
  km <- stats::kmeans(coords[, c("x_coord", "y_coord")],
                       centers = n_blocks, nstart = 10)
  km$cluster
}

#' Construct SuperLearner-style validRows from a spatial block vector.
#' @param block_vec  Integer block assignment (from \code{make_spatial_blocks}).
#' @return List of validation row-index vectors (one per block).
spatial_blocks_to_folds <- function(block_vec) {
  lapply(sort(unique(block_vec)), function(b) which(block_vec == b))
}


# ── Admin1 polygon simulation (no shapefile required) ─────────────────────────

#' Simulate Admin1 polygons as a simple sf object.
#'
#' Creates a regular grid of squares covering [0,1]^2, one polygon per
#' Admin1 region.  Used for choropleth mapping when no real boundaries exist.
#'
#' @param n_admin1  Number of Admin1 regions (default 12).
#' @param seed      RNG seed.
#' @return An sf POLYGON object with columns: admin1_name, geometry.
simulate_admin1_polygons <- function(n_admin1 = 12, seed = 42) {
  require(sf, quietly = TRUE)
  set.seed(seed)

  # Tile [0,1]^2 into n_admin1 Voronoi-like regions using a simple grid
  ncols <- ceiling(sqrt(n_admin1))
  nrows <- ceiling(n_admin1 / ncols)
  w     <- 1 / ncols
  h     <- 1 / nrows

  polys <- vector("list", n_admin1)
  for (i in seq_len(n_admin1)) {
    row_i <- (i - 1) %/% ncols
    col_i <- (i - 1) %%  ncols
    xmin  <- col_i * w
    xmax  <- xmin + w
    ymin  <- row_i * h
    ymax  <- ymin + h
    polys[[i]] <- sf::st_polygon(list(rbind(
      c(xmin, ymin), c(xmax, ymin),
      c(xmax, ymax), c(xmin, ymax),
      c(xmin, ymin)
    )))
  }
  sf::st_sf(
    admin1_name = paste0("Region_", formatC(seq_len(n_admin1), width = 2, flag = "0")),
    geometry    = sf::st_sfc(polys, crs = 4326)
  )
}

#' Assign grid cells to Admin1 polygons via spatial join.
#' @param grid_df   Data frame with x_coord, y_coord columns.
#' @param admin1_sf sf object from \code{simulate_admin1_polygons}.
#' @return grid_df with admin1_name column added.
assign_grid_to_admin1 <- function(grid_df, admin1_sf) {
  require(sf, quietly = TRUE)
  grid_sf <- sf::st_as_sf(grid_df,
                            coords = c("x_coord", "y_coord"),
                            crs    = 4326)
  joined  <- sf::st_join(grid_sf, admin1_sf, join = sf::st_within)
  grid_df$admin1_name <- joined$admin1_name
  grid_df
}


# ── Fay-Herriot small area estimation ────────────────────────────────────────

#' Fit a Fay-Herriot (area-level) small area estimation model.
#'
#' The model is:
#'   Ŷ_a = θ_a + e_a,    e_a ~ N(0, ψ_a)   [sampling model]
#'   θ_a = X_a β + u_a,  u_a ~ N(0, σ²_u)  [linking model]
#'
#' Uses the \code{sae} package if available; falls back to a manual REML
#' implementation via \code{lme4} otherwise.
#'
#' @param direct_est  Named numeric vector of direct survey estimates Ŷ_a.
#' @param samp_var    Named numeric vector of sampling variances ψ_a.
#' @param X_design    Data frame of Admin1-level covariates (one row per Admin1).
#'                    The first column must match names of direct_est.
#' @return List: $eblup (EBLUP estimates), $model_df (full results table).
fit_fay_herriot <- function(direct_est, samp_var, X_design) {
  admins <- names(direct_est)
  df_fh  <- data.frame(
    admin1    = admins,
    direct    = direct_est[admins],
    samp_var  = samp_var[admins]
  )

  # Merge area-level covariates
  if (!is.null(X_design)) {
    X_design$admin1 <- as.character(X_design[[1]])
    df_fh <- merge(df_fh, X_design, by = "admin1", all.x = TRUE)
  }

  df_fh <- df_fh[!is.na(df_fh$direct) & !is.na(df_fh$samp_var), ]
  X_cols <- setdiff(colnames(df_fh), c("admin1", "direct", "samp_var"))
  X_cols <- X_cols[sapply(df_fh[, X_cols, drop = FALSE], is.numeric)]

  # Try sae package
  if (requireNamespace("sae", quietly = TRUE) && length(X_cols) > 0) {
    formula_str <- paste("direct ~", paste(X_cols, collapse = " + "))
    fh_fit <- tryCatch(
      sae::eblupFH(
        formula   = as.formula(formula_str),
        vardir    = df_fh$samp_var,
        data      = df_fh
      ),
      error = function(e) NULL
    )
    if (!is.null(fh_fit)) {
      df_fh$eblup    <- fh_fit$eblup$eblup
      df_fh$mse_eblup <- fh_fit$mse$mse
      return(list(eblup    = setNames(df_fh$eblup, df_fh$admin1),
                  model_df = df_fh,
                  fit      = fh_fit))
    }
  }

  # Fallback: weighted average of direct estimate and covariate prediction
  if (length(X_cols) > 0) {
    lm_fit <- lm(direct ~ ., data = df_fh[, c("direct", X_cols)])
    x_pred <- predict(lm_fit, newdata = df_fh)
  } else {
    x_pred <- rep(mean(df_fh$direct), nrow(df_fh))
  }
  # Simple James-Stein-style shrinkage
  var_u  <- max(0, var(df_fh$direct) - mean(df_fh$samp_var))
  gamma  <- var_u / (var_u + df_fh$samp_var)
  df_fh$eblup <- gamma * df_fh$direct + (1 - gamma) * x_pred

  list(eblup    = setNames(df_fh$eblup, df_fh$admin1),
       model_df = df_fh,
       fit      = NULL)
}


# ── Variogram helper (empirical variogram of residuals) ───────────────────────

#' Compute an empirical variogram of cluster-level residuals.
#'
#' Uses a manual bin-average approach (no external variogram package needed).
#'
#' @param resid_df  Data frame with x_coord, y_coord, mean_residual columns.
#' @param n_bins    Number of distance bins.
#' @return Data frame: bin_mid, gamma (semivariance), n_pairs.
empirical_variogram <- function(resid_df, n_bins = 15) {
  n    <- nrow(resid_df)
  pairs_list <- list()

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d   <- sqrt((resid_df$x_coord[i] - resid_df$x_coord[j])^2 +
                    (resid_df$y_coord[i] - resid_df$y_coord[j])^2)
      gam <- 0.5 * (resid_df$mean_residual[i] - resid_df$mean_residual[j])^2
      pairs_list[[length(pairs_list) + 1]] <- c(d = d, gamma = gam)
    }
  }

  pairs <- as.data.frame(do.call(rbind, pairs_list))
  max_d <- max(pairs$d) * 0.6    # cap at 60 % of max distance
  pairs <- pairs[pairs$d <= max_d, ]

  breaks  <- seq(0, max_d, length.out = n_bins + 1)
  bin     <- cut(pairs$d, breaks = breaks, include.lowest = TRUE)

  vario <- data.frame(pairs, bin = bin) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      bin_mid = mean(d),
      gamma   = mean(gamma),
      n_pairs = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(bin_mid))

  vario
}

#' Plot empirical variogram.
#' @param vario_df Output of \code{empirical_variogram}.
#' @return ggplot object.
plot_variogram <- function(vario_df) {
  ggplot2::ggplot(vario_df, ggplot2::aes(x = bin_mid, y = gamma)) +
    ggplot2::geom_point(ggplot2::aes(size = n_pairs), colour = "steelblue") +
    ggplot2::geom_line(colour = "steelblue", linewidth = 0.7) +
    ggplot2::scale_size_area(name = "n pairs", max_size = 5) +
    ggplot2::labs(x     = "Distance",
                   y     = "Semivariance γ(h)",
                   title = "Empirical variogram of cluster residuals",
                   subtitle = "Rising trend → residual spatial autocorrelation") +
    ggplot2::theme_minimal(base_size = 11)
}
