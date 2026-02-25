# =============================================================================
# R/uncertainty.R
#
# Functions for uncertainty quantification of Admin1 prevalence estimates:
#   1. Nonparametric cluster bootstrap (resamples clusters, refits SL).
#   2. Fast bootstrap using stored CV predictions (no refit).
#   3. Delta-method variance for plug-in prevalence.
#   4. GAM-based Monte Carlo uncertainty propagation.
#   5. Cross-fold variability as a heuristic.
# =============================================================================

# ── 1. Full cluster bootstrap (resamples clusters + refits model) ─────────────

#' Nonparametric cluster bootstrap for Admin1 prevalence CIs.
#'
#' For each replicate:
#'   1. Resample clusters WITH replacement (preserving cluster structure).
#'   2. Refit the chosen model on the resampled data using \code{fit_fun}.
#'   3. Predict on the ORIGINAL full dataset using \code{pred_fun}.
#'   4. Aggregate to Admin1 prevalence.
#' Returns B Admin1 prevalence vectors from which percentile CIs are computed.
#'
#' NOTE: For large SL libraries this can be slow (minutes).  Use
#' \code{bootstrap_fast_admin1} for a faster alternative that reuses CV preds.
#'
#' @param df          Data frame (individual level).
#' @param cluster_col Column name for cluster identifier.
#' @param admin1_col  Column name for Admin1 identifier.
#' @param B           Number of bootstrap replicates.
#' @param fit_fun     Function(df_boot) → model object.
#' @param pred_fun    Function(model, df_orig) → numeric vector of predictions.
#' @param seed        RNG seed.
#' @param verbose     Print progress?
#' @return List: $boot_mat (B × n_admin1 matrix), $ci (data frame with CIs).
bootstrap_admin1 <- function(df,
                               cluster_col = "cluster_id",
                               admin1_col  = "admin1_name",
                               B           = 100,
                               fit_fun,
                               pred_fun,
                               seed        = 42,
                               verbose     = TRUE) {
  set.seed(seed)
  admin1_levels <- sort(unique(df[[admin1_col]]))
  clusters      <- unique(df[[cluster_col]])
  boot_mat      <- matrix(NA_real_, nrow = B, ncol = length(admin1_levels),
                           dimnames = list(NULL, admin1_levels))

  for (b in seq_len(B)) {
    if (verbose && b %% 10 == 0)
      cat(sprintf("  Bootstrap replicate %d / %d\n", b, B))

    # Resample clusters with replacement
    boot_clusters <- sample(clusters, size = length(clusters), replace = TRUE)
    df_b <- do.call(rbind, lapply(boot_clusters, function(cl)
      df[df[[cluster_col]] == cl, , drop = FALSE]))

    # Refit and predict on original data
    model_b <- tryCatch(fit_fun(df_b), error = function(e) NULL)
    if (is.null(model_b)) next

    preds_b <- tryCatch(pred_fun(model_b, df), error = function(e) NULL)
    if (is.null(preds_b) || length(preds_b) != nrow(df)) next

    # Admin1 aggregation
    df_tmp  <- df; df_tmp$pred_b <- preds_b
    a1_prev <- aggregate_to_admin1(df_tmp, admin1_col, "pred_b")
    idx     <- match(a1_prev[[admin1_col]], admin1_levels)
    boot_mat[b, idx] <- a1_prev$pred_prevalence
  }

  # Compute percentile CIs
  ci <- data.frame(
    admin1    = admin1_levels,
    boot_mean = apply(boot_mat, 2, mean,            na.rm = TRUE),
    ci_lo     = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
    ci_hi     = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE),
    n_valid   = apply(boot_mat, 2, function(x) sum(!is.na(x)))
  )
  ci$ci_width <- ci$ci_hi - ci$ci_lo

  list(boot_mat = boot_mat, ci = ci)
}


# ── 2. Fast bootstrap using stored CV predictions ─────────────────────────────

#' Fast cluster bootstrap using STORED cross-validated predictions.
#'
#' Instead of refitting, this resamples individuals (cluster-blocked)
#' from their stored CV predictions.  This is a VARIANCE ESTIMATE ONLY
#' (it does not account for model refitting variance) and should be
#' interpreted as a lower bound on true bootstrap uncertainty.
#'
#' Use as an approximation when full refitting is too slow.
#'
#' @param df          Data frame with cv_pred_col and admin1/cluster columns.
#' @param cluster_col Cluster column name.
#' @param admin1_col  Admin1 column name.
#' @param pred_col    Column of cross-validated predicted probabilities.
#' @param B           Number of bootstrap replicates.
#' @param seed        RNG seed.
#' @return List: $boot_mat (B × n_admin1), $ci (data frame with CIs).
bootstrap_fast_admin1 <- function(df,
                                   cluster_col = "cluster_id",
                                   admin1_col  = "admin1_name",
                                   pred_col    = "pred_bin",
                                   B           = 200,
                                   seed        = 42) {
  set.seed(seed)
  admin1_levels <- sort(unique(df[[admin1_col]]))
  clusters      <- unique(df[[cluster_col]])
  boot_mat      <- matrix(NA_real_, nrow = B, ncol = length(admin1_levels),
                           dimnames = list(NULL, admin1_levels))

  for (b in seq_len(B)) {
    boot_clusters <- sample(clusters, size = length(clusters), replace = TRUE)
    df_b <- do.call(rbind, lapply(boot_clusters, function(cl)
      df[df[[cluster_col]] == cl, , drop = FALSE]))

    a1 <- df_b %>%
      dplyr::group_by(.data[[admin1_col]]) %>%
      dplyr::summarise(prev = mean(.data[[pred_col]], na.rm = TRUE),
                       .groups = "drop")
    idx <- match(a1[[admin1_col]], admin1_levels)
    boot_mat[b, idx] <- a1$prev
  }

  ci <- data.frame(
    admin1    = admin1_levels,
    boot_mean = apply(boot_mat, 2, mean,            na.rm = TRUE),
    ci_lo     = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
    ci_hi     = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
  ci$ci_width <- ci$ci_hi - ci$ci_lo

  list(boot_mat = boot_mat, ci = ci)
}


# ── 3. Delta-method variance for plug-in prevalence ──────────────────────────

#' Delta-method variance for plug-in Admin1 prevalence.
#'
#' For an Admin1 region a with N_a observations, the plug-in estimator is:
#'   p̂_a = (1/N_a) Σ p̂_i
#' Under weak dependence assumptions, Var(p̂_a) ≈ (1/N_a²) Σ Var(p̂_i).
#' Here we use the design-based (cluster-mean) variance approximation.
#'
#' @param df          Data frame with pred_col, admin1_col, cluster_col.
#' @param admin1_col  Admin1 column name.
#' @param cluster_col Cluster column name.
#' @param pred_col    Predicted probability column.
#' @param alpha       Significance level for CI (default 0.05 → 95% CI).
#' @return Data frame: admin1, mean_pred, var_delta, se_delta, ci_lo, ci_hi.
delta_method_admin1 <- function(df,
                                 admin1_col  = "admin1_name",
                                 cluster_col = "cluster_id",
                                 pred_col    = "pred_bin",
                                 alpha       = 0.05) {
  z <- qnorm(1 - alpha / 2)

  df %>%
    dplyr::group_by(.data[[admin1_col]], .data[[cluster_col]]) %>%
    dplyr::summarise(clust_mean = mean(.data[[pred_col]], na.rm = TRUE),
                     n_clust    = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(
      mean_pred   = mean(clust_mean),
      n_clusters  = dplyr::n(),
      var_between = var(clust_mean) / dplyr::n(),  # cluster-mean variance / n
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      se_delta = sqrt(var_between),
      ci_lo    = pmax(0, mean_pred - z * se_delta),
      ci_hi    = pmin(1, mean_pred + z * se_delta)
    )
}


# ── 4. GAM-based Monte Carlo uncertainty ──────────────────────────────────────

#' Monte Carlo uncertainty from a fitted GAM (mgcv).
#'
#' Draws n_sim realisations from the posterior of GAM coefficients
#' (Gaussian approximation from mgcv), then aggregates to Admin1 CIs.
#'
#' @param gam_fit     A fitted mgcv::gam object.
#' @param newdata     Data frame for prediction (with Admin1 column).
#' @param admin1_col  Admin1 column in newdata.
#' @param n_sim       Number of Monte Carlo draws (default 500).
#' @param alpha       CI significance level.
#' @return Data frame: admin1, gam_mean, ci_lo, ci_hi.
gam_uncertainty_admin1 <- function(gam_fit, newdata, admin1_col = "admin1_name",
                                    n_sim = 500, alpha = 0.05) {
  require(mgcv, quietly = TRUE)

  # Simulate from posterior of linear predictor
  Xp    <- predict(gam_fit, newdata = newdata, type = "lpmatrix")
  beta  <- coef(gam_fit)
  Vb    <- vcov(gam_fit)

  # Sample from N(beta, Vb)
  beta_sim <- MASS::mvrnorm(n_sim, beta, Vb)
  lp_sim   <- Xp %*% t(beta_sim)               # n_pred × n_sim

  # Apply link inverse
  inv_link <- gam_fit$family$linkinv
  p_sim    <- apply(lp_sim, 2, inv_link)        # n_pred × n_sim

  # Aggregate to Admin1 for each draw
  admin1_vec  <- newdata[[admin1_col]]
  admin1_levs <- sort(unique(admin1_vec))
  prev_mat    <- matrix(NA_real_, nrow = n_sim, ncol = length(admin1_levs),
                         dimnames = list(NULL, admin1_levs))

  for (s in seq_len(n_sim)) {
    for (a in admin1_levs) {
      idx <- which(admin1_vec == a)
      if (length(idx) > 0)
        prev_mat[s, a] <- mean(p_sim[idx, s])
    }
  }

  data.frame(
    admin1   = admin1_levs,
    gam_mean = apply(prev_mat, 2, mean,            na.rm = TRUE),
    ci_lo    = apply(prev_mat, 2, quantile, alpha / 2,   na.rm = TRUE),
    ci_hi    = apply(prev_mat, 2, quantile, 1 - alpha / 2, na.rm = TRUE)
  )
}


# ── 5. Cross-fold variability heuristic ──────────────────────────────────────

#' Summarise fold-to-fold Admin1 prevalence variability.
#'
#' Computes predictions within each CV fold validation set and measures
#' the spread across folds.  This is a heuristic for model instability
#' (NOT a valid confidence interval).
#'
#' @param df          Data frame with fold_col, pred_col, admin1_col.
#' @param fold_col    Column of fold assignments (integer 1:V).
#' @param pred_col    Predicted probability column.
#' @param admin1_col  Admin1 column name.
#' @return Data frame: admin1, fold_sd (SD across folds), min, max.
cv_fold_variability <- function(df,
                                 fold_col   = "fold",
                                 pred_col   = "pred_bin",
                                 admin1_col = "admin1_name") {
  fold_prevs <- df %>%
    dplyr::group_by(.data[[admin1_col]], .data[[fold_col]]) %>%
    dplyr::summarise(prev = mean(.data[[pred_col]], na.rm = TRUE),
                     .groups = "drop")

  fold_prevs %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(
      fold_mean = mean(prev),
      fold_sd   = sd(prev),
      fold_min  = min(prev),
      fold_max  = max(prev),
      .groups   = "drop"
    )
}


# ── Plotting ──────────────────────────────────────────────────────────────────

#' Plot Admin1 prevalence with uncertainty intervals.
#'
#' @param obs_prev   Named numeric vector: true/observed Admin1 prevalence.
#' @param ci_df      Data frame with columns: admin1, boot_mean/gam_mean, ci_lo, ci_hi.
#' @param mean_col   Column name of the point estimate in ci_df.
#' @param title      Plot title.
#' @return ggplot object.
plot_prevalence_ci <- function(obs_prev, ci_df,
                                mean_col = "boot_mean",
                                title    = "Admin1 prevalence with 95% CIs") {
  ci_df$true_prev <- obs_prev[ci_df$admin1]
  ci_df <- ci_df[order(ci_df[[mean_col]]), ]
  ci_df$admin1 <- factor(ci_df$admin1, levels = ci_df$admin1)

  ggplot2::ggplot(ci_df,
                   ggplot2::aes(y = admin1, x = .data[[mean_col]] * 100)) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_lo * 100, xmax = ci_hi * 100),
      height = 0.4, colour = "steelblue", linewidth = 0.7
    ) +
    ggplot2::geom_point(colour = "steelblue", size = 3) +
    ggplot2::geom_point(
      ggplot2::aes(x = true_prev * 100),
      shape = 4, colour = "firebrick", size = 3, stroke = 1.2
    ) +
    ggplot2::labs(
      x       = "Prevalence (%)",
      y       = NULL,
      title   = title,
      caption = "● = estimated (95% CI); ✕ = true simulated prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
