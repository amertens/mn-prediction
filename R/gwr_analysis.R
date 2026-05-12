# =============================================================================
# R/gwr_analysis.R
#
# Supplementary spatial analysis at Admin2:
#   - Moran's I on SL residuals (spatial autocorrelation diagnostic)
#   - Geographically Weighted Regression (GWR) of survey-weighted prevalence
#     on aggregated proxy predictors, to map locally-varying determinants.
#
# Purpose: interpretive supplement to the global SuperLearner pipeline;
# does NOT introduce gw_ survey variables (proxy-only constraint preserved).
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(spdep)
})


#' Load Admin2 polygons + centroids for a country
#'
#' Uses the same GADM cache as scripts/06_admin2_predictions_map.R.
#'
#' @param gadm_code 3-letter ISO code (e.g. "GMB", "GHA")
#' @return sf with columns Admin1, Admin2, geometry, plus lon/lat centroids
load_admin2_centroids <- function(gadm_code) {
  gadm_raw <- geodata::gadm(gadm_code, level = 2,
                            path = here::here("data", "gadm"))
  polys <- sf::st_as_sf(gadm_raw)
  polys$Admin1 <- polys$NAME_1
  polys$Admin2 <- polys$NAME_2

  # Centroids in WGS84 lon/lat — sufficient for kNN weights and adaptive GWR.
  ctr <- suppressWarnings(sf::st_centroid(polys))
  coords <- sf::st_coordinates(ctr)
  polys$lon <- coords[, 1]
  polys$lat <- coords[, 2]

  polys[, c("Admin1", "Admin2", "lon", "lat", "geometry")]
}


#' Moran's I on Admin2 residuals
#'
#' Builds a row-standardised k-nearest-neighbour weights matrix on Admin2
#' centroids and tests for spatial autocorrelation in (svy_prev - sl_prev).
#'
#' @param admin2_df data.frame with columns Admin1, Admin2, svy_prev, sl_prev
#' @param polys sf from load_admin2_centroids() (must include same Admin2s)
#' @param k_neighbors integer (default 5)
#' @return list with $morans_i (data.frame, 1 row), $listw (for reuse)
morans_i_residuals <- function(admin2_df, polys, k_neighbors = 5) {
  d <- admin2_df %>%
    dplyr::filter(!is.na(svy_prev), !is.na(sl_prev)) %>%
    dplyr::mutate(residual = svy_prev - sl_prev)

  joined <- polys %>%
    sf::st_drop_geometry() %>%
    dplyr::inner_join(d, by = c("Admin1", "Admin2"))

  if (nrow(joined) < (k_neighbors + 2)) {
    return(list(
      morans_i = data.frame(
        n = nrow(joined), morans_i = NA_real_, expectation = NA_real_,
        variance = NA_real_, z = NA_real_, p_value = NA_real_,
        note = "too few admin2 with both survey and SL predictions"
      ),
      listw = NULL
    ))
  }

  coords <- as.matrix(joined[, c("lon", "lat")])
  k_use  <- min(k_neighbors, nrow(coords) - 1)
  knn    <- spdep::knearneigh(coords, k = k_use)
  nb     <- spdep::knn2nb(knn)
  lw     <- spdep::nb2listw(nb, style = "W")

  mt <- spdep::moran.test(joined$residual, lw, randomisation = TRUE,
                          alternative = "two.sided")

  list(
    morans_i = data.frame(
      n           = nrow(joined),
      k_neighbors = k_use,
      morans_i    = unname(mt$estimate["Moran I statistic"]),
      expectation = unname(mt$estimate["Expectation"]),
      variance    = unname(mt$estimate["Variance"]),
      z           = unname(mt$statistic),
      p_value     = mt$p.value,
      note        = NA_character_
    ),
    listw = lw
  )
}


#' Aggregate individual-level proxy predictors to Admin2 means
#'
#' Strips gw_ columns to enforce proxy-only constraint, then computes
#' simple (n-weighted) means per Admin2 for numeric predictors.
#'
#' @param ind_df individual-level data.frame (must have Admin1, Admin2)
#' @param keep_vars character vector of predictor column names to aggregate
#' @return data.frame with Admin1, Admin2, n_ind, and one column per keep_var
aggregate_predictors_admin2 <- function(ind_df, keep_vars) {
  drop_gw <- grep("^gw_", keep_vars, value = TRUE)
  if (length(drop_gw) > 0) {
    warning(sprintf("[aggregate_predictors_admin2] Dropping %d gw_ vars (proxy-only): %s",
                    length(drop_gw), paste(drop_gw, collapse = ", ")))
    keep_vars <- setdiff(keep_vars, drop_gw)
  }

  keep_vars <- intersect(keep_vars, colnames(ind_df))
  if (length(keep_vars) == 0) {
    stop("No usable predictor columns found in individual data.")
  }

  ind_df %>%
    dplyr::filter(!is.na(Admin1), !is.na(Admin2)) %>%
    dplyr::group_by(Admin1, Admin2) %>%
    dplyr::summarise(
      n_ind = dplyr::n(),
      dplyr::across(dplyr::all_of(keep_vars),
                    ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}


#' Pick GWR predictors
#'
#' Selects up to `top_k` numeric predictors with VIF <= max_vif,
#' ranked by absolute Pearson correlation with the outcome.
#' (Used in the absence of pre-computed SHAP rankings.)
#'
#' @param df data.frame containing y_col and candidate columns
#' @param y_col outcome column name
#' @param candidates character vector of candidate predictors
#' @param top_k integer (default 6)
#' @param max_vif numeric (default 5)
#' @return character vector of selected predictor names
select_gwr_predictors <- function(df, y_col, candidates, top_k = 6, max_vif = 5) {
  candidates <- intersect(candidates, colnames(df))
  candidates <- candidates[vapply(candidates,
                                  function(v) is.numeric(df[[v]]) &&
                                    sum(!is.na(df[[v]])) >= 0.5 * nrow(df) &&
                                    stats::var(df[[v]], na.rm = TRUE) > 0,
                                  logical(1))]
  if (length(candidates) == 0) return(character(0))

  y <- df[[y_col]]
  cors <- vapply(candidates, function(v) {
    suppressWarnings(stats::cor(df[[v]], y, use = "pairwise.complete.obs"))
  }, numeric(1))
  ord <- candidates[order(abs(cors), decreasing = TRUE)]

  selected <- character(0)
  for (v in ord) {
    trial <- c(selected, v)
    if (length(trial) < 2) { selected <- trial; next }
    fmla <- stats::as.formula(paste(y_col, "~", paste(trial, collapse = " + ")))
    sub  <- df[, c(y_col, trial), drop = FALSE]
    sub  <- sub[stats::complete.cases(sub), , drop = FALSE]
    if (nrow(sub) <= length(trial) + 1) next
    fit  <- tryCatch(stats::lm(fmla, data = sub), error = function(e) NULL)
    if (is.null(fit)) next
    vifs <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(vifs) || any(vifs > max_vif, na.rm = TRUE)) next
    selected <- trial
    if (length(selected) >= top_k) break
  }
  selected
}


#' Fit Admin2 GWR
#'
#' Runs global OLS, selects adaptive bandwidth via CV, fits GWR with a
#' Gaussian kernel, then compares global vs local fit.
#'
#' @param df data.frame including Admin1, Admin2, lon, lat, y_col, predictors
#' @param y_col outcome column (e.g. "svy_prev")
#' @param predictors character vector of predictor column names
#' @return list: $global_fit, $gwr_fit, $bandwidth, $summary, $local_coefs
fit_admin2_gwr <- function(df, y_col, predictors) {
  stopifnot(all(c("lon", "lat", y_col) %in% colnames(df)))
  predictors <- intersect(predictors, colnames(df))
  if (length(predictors) == 0) stop("No usable predictors for GWR.")

  d <- df[stats::complete.cases(df[, c(y_col, predictors, "lon", "lat")]), , drop = FALSE]
  if (nrow(d) < length(predictors) + 5) {
    stop(sprintf("Too few admin2 (%d) for GWR with %d predictors.",
                 nrow(d), length(predictors)))
  }

  fmla <- stats::as.formula(paste(y_col, "~", paste(predictors, collapse = " + ")))

  global_fit <- stats::lm(fmla, data = d)

  sp_df <- as.data.frame(d)
  sp::coordinates(sp_df) <- ~ lon + lat

  bw <- GWmodel::bw.gwr(formula = fmla, data = sp_df,
                       approach = "CV", kernel = "gaussian",
                       adaptive = TRUE)

  gwr_fit <- GWmodel::gwr.basic(formula = fmla, data = sp_df,
                               bw = bw, kernel = "gaussian",
                               adaptive = TRUE)

  sdf <- as.data.frame(gwr_fit$SDF)
  beta_cols <- predictors
  tv_cols   <- paste0(predictors, "_TV")

  local_coefs <- data.frame(
    Admin1 = d$Admin1,
    Admin2 = d$Admin2,
    lon    = d$lon,
    lat    = d$lat,
    Intercept     = sdf[["Intercept"]],
    local_R2      = sdf[["Local_R2"]],
    residual      = sdf[["residual"]]
  )
  for (p in predictors) {
    local_coefs[[paste0("beta_", p)]] <- sdf[[p]]
    if (paste0(p, "_TV") %in% colnames(sdf)) {
      tv <- sdf[[paste0(p, "_TV")]]
      local_coefs[[paste0("t_", p)]]    <- tv
      local_coefs[[paste0("pval_", p)]] <- 2 * stats::pnorm(-abs(tv))
    }
  }
  for (p in predictors) {
    pcol <- paste0("pval_", p)
    if (pcol %in% colnames(local_coefs)) {
      local_coefs[[paste0("padj_", p)]] <-
        stats::p.adjust(local_coefs[[pcol]], method = "BH")
    }
  }

  aic_global <- stats::AIC(global_fit)
  aic_gwr    <- gwr_fit$GW.diagnostic$AIC
  r2_global  <- summary(global_fit)$r.squared
  r2_gwr     <- gwr_fit$GW.diagnostic$gw.R2

  ftest <- tryCatch(GWmodel::BFC02.gwr.test(gwr_fit), error = function(e) NULL)

  summary_df <- data.frame(
    n             = nrow(d),
    n_predictors  = length(predictors),
    bandwidth     = bw,
    kernel        = "gaussian",
    adaptive      = TRUE,
    aic_global    = aic_global,
    aic_gwr       = aic_gwr,
    aic_delta     = aic_global - aic_gwr,
    r2_global     = r2_global,
    r2_gwr        = r2_gwr,
    f_statistic   = if (!is.null(ftest)) unname(ftest$statistic) else NA_real_,
    f_p_value     = if (!is.null(ftest)) ftest$p.value else NA_real_,
    predictors    = paste(predictors, collapse = ";")
  )

  list(
    global_fit  = global_fit,
    gwr_fit     = gwr_fit,
    bandwidth   = bw,
    summary     = summary_df,
    local_coefs = local_coefs
  )
}


#' Run the full Admin2 GWR supplementary pipeline for one country/outcome
#'
#' @param admin2_pred_csv path to results/shared/admin2_predictions_*.csv
#' @param ind_data individual-level data.frame (with Admin1/Admin2) OR NULL
#'                 (if NULL, GWR is skipped and only Moran's I is returned)
#' @param gadm_code GADM 3-letter code
#' @param candidate_predictors character vector (proxy-only) for GWR
#' @param k_neighbors for Moran's I
#' @param top_k max GWR predictors
#' @return list with $morans, $gwr_summary, $gwr_local_coefs, $polys
run_admin2_gwr_pipeline <- function(admin2_pred_csv,
                                    ind_data = NULL,
                                    gadm_code,
                                    candidate_predictors = NULL,
                                    k_neighbors = 5,
                                    top_k = 6,
                                    skip_gwr = FALSE) {
  admin2_df <- utils::read.csv(admin2_pred_csv, stringsAsFactors = FALSE)
  polys <- load_admin2_centroids(gadm_code)

  morans <- morans_i_residuals(admin2_df, polys, k_neighbors = k_neighbors)$morans_i

  out <- list(
    morans            = morans,
    gwr_summary       = NULL,
    gwr_local_coefs   = NULL,
    polys             = polys,
    selected_predictors = character(0)
  )

  if (skip_gwr || is.null(ind_data) || is.null(candidate_predictors) ||
      length(candidate_predictors) == 0) {
    return(out)
  }

  # Drop coordinate-like columns from candidates: they would clash with the
  # polygon centroids and don't make sense as Admin2-aggregated predictors.
  candidate_predictors <- setdiff(
    candidate_predictors,
    c("lon", "lat", "longitude", "latitude", "LON", "LAT")
  )

  agg <- aggregate_predictors_admin2(ind_data, candidate_predictors)

  centroids <- sf::st_drop_geometry(polys[, c("Admin1", "Admin2", "lon", "lat")])

  d <- admin2_df %>%
    dplyr::filter(!is.na(svy_prev)) %>%
    dplyr::inner_join(agg, by = c("Admin1", "Admin2")) %>%
    dplyr::inner_join(centroids, by = c("Admin1", "Admin2"))

  if (!all(c("lon", "lat") %in% colnames(d))) {
    stop("Internal error: lon/lat missing after centroid join.")
  }

  selected <- select_gwr_predictors(
    d, y_col = "svy_prev",
    candidates = setdiff(colnames(agg), c("Admin1", "Admin2", "n_ind",
                                           "lon", "lat")),
    top_k = top_k
  )
  out$selected_predictors <- selected
  if (length(selected) < 2) {
    warning("Fewer than 2 usable predictors after VIF filtering — skipping GWR.")
    return(out)
  }

  fit <- fit_admin2_gwr(d, y_col = "svy_prev", predictors = selected)
  out$gwr_summary     <- fit$summary
  out$gwr_local_coefs <- fit$local_coefs
  out
}
