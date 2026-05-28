# =============================================================================
# R/feature_engineering_bucket_a.R
#
# Bucket A feature-engineering wins from existing data (no external
# dependencies, no downloads, no auth).  Four helpers:
#
#   1. add_spatial_smoothing()      — for each covariate, append the
#       neighbourhood mean within the same country (Queen contiguity).
#
#   2. add_within_domain_pca()      — group covariates by data-source
#       prefix, run PCA per group, take the top K components per group.
#       Replaces the within-group multicollinearity that bites Fay-
#       Herriot and OLS rank-deficiency.
#
#   3. compute_admin2_dhs_aggregates() / add_dhs_admin2_features()
#       — admin-2 means of cluster-level DHS variables.  These are
#       DHS-derived so are tagged with the `dhs_` prefix and the
#       AREA-level benchmark runner should exclude them from LOCO
#       transportability evaluation (held-out countries lack DHS data
#       by construction).  Useful for within-country admin-2 prediction.
#
#   4. compute_admin2_wash_and_gini() — admin-2 WASH composite
#       (mean of improved-water + improved-sanitation indicators) and
#       intra-admin-2 wealth Gini (income inequality proxy within each
#       district).  Also DHS-derived → `dhs_`-tagged.
#
# All four functions return data frames that can be merged into the
# pooled area-level dataset and consumed by the existing benchmark
# methods in R/benchmark_models.R.
#
# Distinct from R/feature_engineering.R (which does temporal anomalies
# and seasonality features from year-suffixed columns).
# =============================================================================


`%||%` <- function(a, b) if (is.null(a)) b else a


# ── 1. Spatial smoothing of admin-2 covariates ───────────────────────────────
# For each covariate in `vars`, add a `<var>_smooth` column whose value at
# admin-2 polygon i is the mean of that covariate over polygon i's
# Queen-contiguity neighbours within the same country.  Polygons with no
# neighbours (islands etc.) get their own value.
#
# Inputs:
#   pooled_data       — data.frame with columns `country`, `Admin2`, and the
#                       variables in `vars`.
#   vars              — character vector of covariate column names.
#   adjacency_list    — named list (country -> spdep::nb) as produced by
#                       build_adjacency_list() in R/benchmark_models.R.
#   svy_admin2_list   — named list of svy_admin2 data.frames (country ->
#                       data.frame with at least an `Admin2` column).
#                       Used to align nb indices with Admin2 names.
#   suffix            — suffix appended to smoothed-variable column names
#                       (default "_smooth").
#
# Returns the input data.frame with new `<var><suffix>` columns appended.
add_spatial_smoothing <- function(pooled_data, vars, adjacency_list,
                                    svy_admin2_list, suffix = "_smooth") {
  if (length(vars) == 0 || length(adjacency_list) == 0) return(pooled_data)
  smooth_mat <- matrix(NA_real_, nrow = nrow(pooled_data), ncol = length(vars))
  colnames(smooth_mat) <- paste0(vars, suffix)

  for (ctry in names(adjacency_list)) {
    nb <- adjacency_list[[ctry]]
    svy <- svy_admin2_list[[ctry]]
    if (is.null(svy) || is.null(nb) || length(nb) != nrow(svy)) next
    pdrows <- which(pooled_data$country == ctry)
    if (length(pdrows) == 0) next
    poly_idx <- match(pooled_data$Admin2[pdrows], svy$Admin2)
    for (j in seq_along(vars)) {
      v <- vars[j]
      vals <- pooled_data[[v]][pdrows]
      for (i in seq_along(pdrows)) {
        pi <- poly_idx[i]
        if (is.na(pi)) { smooth_mat[pdrows[i], j] <- vals[i]; next }
        nbi <- nb[[pi]]
        if (length(nbi) == 0 || identical(nbi, 0L)) {
          smooth_mat[pdrows[i], j] <- vals[i]
        } else {
          nb_names <- svy$Admin2[nbi]
          nb_rows  <- which(pooled_data$country == ctry &
                              pooled_data$Admin2 %in% nb_names)
          nb_vals <- pooled_data[[v]][nb_rows]
          smooth_mat[pdrows[i], j] <- mean(nb_vals, na.rm = TRUE)
        }
      }
    }
  }
  for (j in seq_along(vars)) {
    bad <- !is.finite(smooth_mat[, j])
    if (any(bad)) smooth_mat[bad, j] <- pooled_data[[vars[j]]][bad]
  }
  out <- cbind(pooled_data, as.data.frame(smooth_mat))
  attr(out, "smooth_vars") <- paste0(vars, suffix)
  out
}


# ── 2. Within-domain PCA ────────────────────────────────────────────────────
# Groups covariates by data-source prefix (via .assign_variable_groups() in
# R/benchmark_models.R), runs prcomp per group on the pooled training data,
# extracts the top `top_k` principal components, and appends them as new
# `pc_<group>_<i>` columns.
add_within_domain_pca <- function(pooled_data, vars,
                                    top_k = 3L, min_group_size = 3L,
                                    drop_originals = FALSE) {
  if (length(vars) < min_group_size) return(pooled_data)
  if (!exists(".assign_variable_groups")) {
    # Source benchmark_models.R for the grouping helper if not loaded.
    return(pooled_data)
  }
  groups <- .assign_variable_groups(vars)
  group_levels <- unique(groups)
  new_cols <- character()
  for (g in group_levels) {
    gv <- vars[groups == g]
    if (length(gv) < min_group_size) next
    X <- as.matrix(pooled_data[, gv, drop = FALSE])
    for (j in seq_len(ncol(X))) {
      bad <- !is.finite(X[, j])
      if (any(bad)) X[bad, j] <- median(X[, j], na.rm = TRUE)
    }
    pr <- tryCatch(stats::prcomp(X, center = TRUE, scale. = TRUE),
                    error = function(e) NULL)
    if (is.null(pr)) next
    k <- min(top_k, ncol(pr$x))
    cols <- pr$x[, seq_len(k), drop = FALSE]
    # Build a readable label from the variable prefix of the first
    # variable in the group.
    label_token <- sub("^(gee|ihme|mal_atlas|gdl|spam|chirps|wfp|map|wp|esa)_",
                        "", gv[1])
    tokens <- strsplit(label_token, "_")[[1]]
    label_pieces <- tokens[seq_len(min(2, length(tokens)))]
    label <- gsub("[^a-zA-Z0-9]+", "_", paste(label_pieces, collapse = "_"))
    if (label == "") label <- paste0("g", g)
    new_names <- paste0("pc_", label, "_", seq_len(k))
    # Guard against name collisions across iterations.
    while (any(new_names %in% colnames(pooled_data))) {
      new_names <- paste0(new_names, "x")
    }
    colnames(cols) <- new_names
    pooled_data <- cbind(pooled_data, as.data.frame(cols))
    new_cols <- c(new_cols, new_names)
    if (drop_originals) pooled_data[, gv] <- NULL
  }
  attr(pooled_data, "pca_vars") <- new_cols
  pooled_data
}


# ── 3. Cluster-level DHS aggregates at admin-2 ───────────────────────────────
# Default set of DHS variables to aggregate. Restricted to the clean
# binary WASH indicators that exist consistently across the 4 DHS surveys.
# (The DHS-derived wealth quintile is harder to identify reliably across
# countries — gw_hhwq turned out to be a binary indicator, gw_WQ* are
# country-specific aggregates rather than quintile dummies. Pass a custom
# `vars` argument with the correct wealth column if available.)
.default_dhs_vars <- c(
  "gw_water_12",     # piped into dwelling
  "gw_water_13",     # piped to yard / public tap
  "gw_toilet_11",    # flush toilet
  "gw_toilet_22"     # improved pit latrine
)
compute_admin2_dhs_aggregates <- function(merged_country_list,
                                            vars = .default_dhs_vars,
                                            id_col = "Admin2") {
  out_rows <- list()
  for (ctry in names(merged_country_list)) {
    m <- merged_country_list[[ctry]]
    if (is.null(m) || nrow(m) == 0 || !id_col %in% colnames(m)) next
    vars_present <- intersect(vars, colnames(m))
    if (length(vars_present) == 0) next
    # Numeric-coerce the wealth quintile (it's often haven_labelled).
    df <- m[, c(id_col, vars_present), drop = FALSE]
    for (v in vars_present) {
      df[[v]] <- suppressWarnings(as.numeric(unclass(df[[v]])))
    }
    agg <- df |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(vars_present),
                                      ~ mean(.x, na.rm = TRUE)),
                        n_dhs = dplyr::n(),
                        .groups = "drop") |>
      dplyr::rename_with(.fn = ~ paste0("dhs_", .),
                          .cols = -dplyr::all_of(c(id_col, "n_dhs")))
    agg$country <- ctry
    out_rows[[ctry]] <- agg
  }
  if (length(out_rows) == 0) return(data.frame())
  dplyr::bind_rows(out_rows)
}

add_dhs_admin2_features <- function(pooled_data,
                                      merged_country_list,
                                      vars = .default_dhs_vars,
                                      id_col = "Admin2") {
  dhs <- compute_admin2_dhs_aggregates(merged_country_list, vars, id_col)
  if (nrow(dhs) == 0) return(pooled_data)
  out <- dplyr::left_join(pooled_data, dhs, by = c("country", id_col))
  attr(out, "dhs_vars") <- grep("^dhs_", colnames(out), value = TRUE)
  out
}


# ── 4. Admin-2 WASH composite + intra-admin-2 wealth Gini ───────────────────
gini_coefficient <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x); if (n < 2) return(NA_real_)
  x_sorted <- sort(x)
  sum_x <- sum(x_sorted)
  if (sum_x == 0) return(0)
  i <- seq_len(n)
  (2 * sum(i * x_sorted) / (n * sum_x)) - (n + 1) / n
}

compute_admin2_wash_and_gini <- function(merged_country_list,
                                          id_col = "Admin2",
                                          wealth_var = NULL) {
  # wealth_var: optional name of a continuous-or-ordinal wealth column.
  # If NULL or absent in the data, wealth-Gini / wealth-mean columns are
  # omitted. (No reliable cross-country wealth-quintile column was
  # identified in the gw_* set; pass one explicitly when known.)
  out_rows <- list()
  for (ctry in names(merged_country_list)) {
    m <- merged_country_list[[ctry]]
    if (is.null(m) || nrow(m) == 0 || !id_col %in% colnames(m)) next
    water_vars  <- intersect(c("gw_water_12", "gw_water_13"),  colnames(m))
    toilet_vars <- intersect(c("gw_toilet_11", "gw_toilet_22"), colnames(m))
    wv          <- if (!is.null(wealth_var)) intersect(wealth_var, colnames(m))
                   else character()
    if (length(c(water_vars, toilet_vars, wv)) == 0) next
    df <- m[, c(id_col, water_vars, toilet_vars, wv), drop = FALSE]
    for (v in c(water_vars, toilet_vars)) {
      df[[v]] <- suppressWarnings(as.numeric(unclass(df[[v]])))
      df[[v]][is.na(df[[v]])] <- 0
    }
    if (length(wv) == 1) {
      df[[wv]] <- suppressWarnings(as.numeric(unclass(df[[wv]])))
    }
    has_water  <- length(water_vars)  > 0
    has_toilet <- length(toilet_vars) > 0
    df$.water  <- if (has_water)  rowSums(df[, water_vars,  drop = FALSE]) > 0 else FALSE
    df$.toilet <- if (has_toilet) rowSums(df[, toilet_vars, drop = FALSE]) > 0 else FALSE
    df$.wash   <- (as.numeric(df$.water) + as.numeric(df$.toilet)) /
                   max(1, has_water + has_toilet)
    agg <- df |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(
        dhs_wash_composite = mean(.data$.wash, na.rm = TRUE),
        dhs_water_pct      = mean(as.numeric(.data$.water),  na.rm = TRUE),
        dhs_sanit_pct      = mean(as.numeric(.data$.toilet), na.rm = TRUE),
        dhs_wealth_gini    = if (length(wv) == 1)
                                gini_coefficient(.data[[wv]])
                              else NA_real_,
        dhs_wealth_mean    = if (length(wv) == 1)
                                mean(.data[[wv]], na.rm = TRUE)
                              else NA_real_,
        n_wash             = dplyr::n(),
        .groups = "drop")
    if (length(wv) == 0) {
      agg$dhs_wealth_gini <- NULL; agg$dhs_wealth_mean <- NULL
    }
    agg$country <- ctry
    out_rows[[ctry]] <- agg
  }
  if (length(out_rows) == 0) return(data.frame())
  dplyr::bind_rows(out_rows)
}

add_wash_and_gini_features <- function(pooled_data,
                                         merged_country_list,
                                         id_col = "Admin2") {
  add <- compute_admin2_wash_and_gini(merged_country_list, id_col)
  if (nrow(add) == 0) return(pooled_data)
  # Avoid duplicate columns if user has already called add_dhs_admin2_features.
  drop_cols <- intersect(colnames(add), colnames(pooled_data))
  drop_cols <- setdiff(drop_cols, c("country", id_col))
  if (length(drop_cols) > 0) pooled_data[, drop_cols] <- NULL
  dplyr::left_join(pooled_data, add, by = c("country", id_col))
}


# ── Convenience: apply all four feature-engineering augmentations at once ────
# Returns the augmented data.frame plus attributes recording the new
# variable names added (so downstream methods can screen or exclude as
# needed):
#   attr(out, "gee_vars")            — augmented GEE-only set (safe for LOCO)
#   attr(out, "dhs_vars")            — DHS-derived (NOT safe for LOCO)
apply_bucket_a_features <- function(pooled_data, gee_vars,
                                      adjacency_list = NULL,
                                      svy_admin2_list = NULL,
                                      merged_country_list = NULL,
                                      spatial_smooth = TRUE,
                                      domain_pca = TRUE,
                                      dhs_aggregates = TRUE,
                                      wash_gini = TRUE) {
  augmented <- pooled_data
  new_gee_vars <- gee_vars
  new_dhs_vars <- character()

  if (isTRUE(spatial_smooth) && !is.null(adjacency_list) &&
      !is.null(svy_admin2_list) && length(gee_vars) > 0) {
    augmented <- add_spatial_smoothing(augmented, gee_vars,
                                         adjacency_list, svy_admin2_list)
    smooth_vars <- attr(augmented, "smooth_vars") %||% character()
    new_gee_vars <- c(new_gee_vars, smooth_vars)
    cat(sprintf("  [features] added %d smoothed covariates\n", length(smooth_vars)))
  }
  if (isTRUE(domain_pca) && length(gee_vars) >= 3) {
    augmented <- add_within_domain_pca(augmented, gee_vars)
    pca_vars <- attr(augmented, "pca_vars") %||% character()
    new_gee_vars <- c(new_gee_vars, pca_vars)
    cat(sprintf("  [features] added %d within-domain PCs\n", length(pca_vars)))
  }
  if (isTRUE(dhs_aggregates) && !is.null(merged_country_list)) {
    augmented <- add_dhs_admin2_features(augmented, merged_country_list)
    dhs_vars <- attr(augmented, "dhs_vars") %||% character()
    new_dhs_vars <- c(new_dhs_vars, dhs_vars)
    cat(sprintf("  [features] added %d DHS admin-2 aggregates (dhs_ tagged)\n",
                length(dhs_vars)))
  }
  if (isTRUE(wash_gini) && !is.null(merged_country_list)) {
    augmented <- add_wash_and_gini_features(augmented, merged_country_list)
    wash_cols <- grep("^dhs_wash|^dhs_wealth|^dhs_water_pct|^dhs_sanit_pct",
                       colnames(augmented), value = TRUE)
    new_dhs_vars <- c(new_dhs_vars, wash_cols)
    cat(sprintf("  [features] added %d WASH/wealth composites (dhs_ tagged)\n",
                length(wash_cols)))
  }

  attr(augmented, "gee_vars")           <- new_gee_vars
  attr(augmented, "dhs_vars")           <- unique(new_dhs_vars)
  attr(augmented, "augmented_gee_vars") <- setdiff(new_gee_vars, gee_vars)
  attr(augmented, "augmented_dhs_vars") <- unique(new_dhs_vars)
  augmented
}
