# =============================================================================
# R/oos_prediction.R
#
# Out-of-sample prediction for unsurveyed countries using GEE rasters.
# Takes an area-level model trained on surveyed countries and applies it
# to GADM Admin-2 polygons in a new country, extracting the same GEE
# covariates from locally stored rasters.
# =============================================================================

#' Normalize a raster filename into a canonical GEE variable name
#'
#' Strips country names, standardizes separators, and creates names that
#' match across countries. E.g.:
#'   "CCNL_2013_Gambia.tif"         -> "gee_ccnl_2013"
#'   "CCNL_2013_Cote_dIvoire.tif"   -> "gee_ccnl_2013"
#'   "NDVI_Ghana_2017.tif"          -> "gee_ndvi_2017"
#'   "FLDAS_2017_Annual_Mean_Gambia" -> "gee_fldas_2017_annual_mean"
#'   "FLDAS_2017_Monthly_Cote_dIvoire" -> "gee_fldas_2017_monthly"
#'
#' @param filename Basename of .tif file (without extension)
#' @return Canonical variable name starting with "gee_"
canonicalize_gee_varname <- function(filename) {
  # Remove file extension
  base <- tools::file_path_sans_ext(basename(filename))

  # Country names to strip (order matters — longer first)
  country_patterns <- c(
    "Cote_dIvoire", "Cote_d_Ivoire", "CotedIvoire",
    "Sierra_Leone", "Sierra Leone", "SierraLeone",
    "Gambia", "Ghana", "Malawi", "Nigeria", "Senegal",
    "Burkina_Faso", "BurkinaFaso", "Guinea", "Mali", "Niger",
    "Togo", "Benin", "Liberia"
  )

  cleaned <- base
  for (pat in country_patterns) {
    # Remove country name and any surrounding underscores/spaces
    cleaned <- gsub(paste0("_?", pat, "_?"), "_", cleaned, ignore.case = TRUE)
  }

  # Clean up: lowercase, replace non-alphanumeric with _, collapse multiple _
  cleaned <- tolower(cleaned)
  cleaned <- gsub("[^a-z0-9]+", "_", cleaned)
  cleaned <- gsub("_+", "_", cleaned)
  cleaned <- sub("^_|_$", "", cleaned)

  paste0("gee_", cleaned)
}


#' Extract GEE covariates for a new country's Admin-2 polygons
#'
#' @param gadm_code ISO3 code for geodata::gadm()
#' @param raster_dir Path to folder of .tif rasters
#' @param gadm_path Path for caching GADM downloads
#' @return data.frame with Admin1, Admin2, and gee_* columns
extract_gee_for_country <- function(gadm_code, raster_dir, gadm_path = here::here("data", "gadm")) {

  gadm_raw <- geodata::gadm(gadm_code, level = 2, path = gadm_path)
  gadm <- sf::st_as_sf(gadm_raw)
  gadm$Admin2 <- gadm$NAME_2
  gadm$Admin1 <- gadm$NAME_1

  tif_files <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  cat(sprintf("[oos] %s: %d Admin2 polygons, %d .tif files\n",
              gadm_code, nrow(gadm), length(tif_files)))

  gee_df <- data.frame(Admin1 = gadm$Admin1, Admin2 = gadm$Admin2)

  for (tif in tif_files) {
    varname <- canonicalize_gee_varname(tif)

    r <- tryCatch(terra::rast(tif), error = function(e) NULL)
    if (is.null(r)) next
    if (terra::nlyr(r) > 1) r <- r[[1]]

    vals <- tryCatch(
      exactextractr::exact_extract(r, gadm, fun = "mean"),
      error = function(e) rep(NA_real_, nrow(gadm))
    )

    # If duplicate canonical name (e.g., multiple years), keep first
    if (!varname %in% colnames(gee_df)) {
      gee_df[[varname]] <- vals
    }
  }

  gee_vars <- setdiff(colnames(gee_df), c("Admin1", "Admin2"))
  cat(sprintf("[oos] Extracted %d GEE variables for %s\n", length(gee_vars), gadm_code))

  attr(gee_df, "polygons") <- gadm
  gee_df
}


#' Predict micronutrient deficiency in an unsurveyed country
#'
#' Uses an area-level glmnet model (from fit_area_level_model) trained on
#' surveyed countries' Admin-2 prevalences, and applies it to a new country's
#' Admin-2 GEE covariates.
#'
#' @param area_model Output from fit_area_level_model() for a training country
#'   (must contain fit_area and gee_vars)
#' @param oos_gee Output from extract_gee_for_country() for the target country
#' @param outcome_tag Character label for the outcome
#' @return list with predictions data.frame and polygons
predict_oos_country <- function(area_model, oos_gee, outcome_tag) {

  if (is.null(area_model$fit_area) || is.null(area_model$gee_vars)) {
    warning("Area model has no fit or gee_vars — cannot predict OOS")
    return(NULL)
  }

  model_vars <- area_model$gee_vars
  fit        <- area_model$fit_area
  polygons   <- attr(oos_gee, "polygons")

  # Build predictor matrix, padding missing columns with 0
  available <- intersect(model_vars, colnames(oos_gee))
  missing   <- setdiff(model_vars, colnames(oos_gee))
  cat(sprintf("[oos_predict] %d/%d model variables found in target country (%d missing)\n",
              length(available), length(model_vars), length(missing)))

  if (length(missing) > 0) {
    cat(sprintf("  Missing vars: %s\n", paste(head(missing, 10), collapse = ", ")))
  }

  X_oos <- as.matrix(oos_gee[, available, drop = FALSE])
  # Pad missing columns with 0
  if (length(missing) > 0) {
    pad <- matrix(0, nrow = nrow(X_oos), ncol = length(missing))
    colnames(pad) <- missing
    X_oos <- cbind(X_oos, pad)
  }
  # Reorder to match model
  X_oos <- X_oos[, model_vars, drop = FALSE]

  # Impute any NAs with column medians
  col_medians <- apply(X_oos, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_oos))) {
    na_idx <- is.na(X_oos[, j])
    if (any(na_idx)) X_oos[na_idx, j] <- col_medians[j]
  }

  # Predict
  yhat <- pmin(pmax(
    as.numeric(predict(fit, X_oos, s = "lambda.min")),
    0), 1)

  preds <- data.frame(
    Admin1   = oos_gee$Admin1,
    Admin2   = oos_gee$Admin2,
    pred_prev = yhat,
    stringsAsFactors = FALSE
  )

  cat(sprintf("[oos_predict] %s: %d Admin2 predictions, range [%.1f%%, %.1f%%]\n",
              outcome_tag, nrow(preds),
              100 * min(yhat, na.rm = TRUE), 100 * max(yhat, na.rm = TRUE)))

  # Admin1 summaries
  admin1_preds <- preds %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n_admin2  = dplyr::n(),
      pred_prev = mean(pred_prev, na.rm = TRUE),
      .groups   = "drop"
    )
  cat(sprintf("[oos_predict] Admin1 summary: %d regions\n", nrow(admin1_preds)))

  list(
    admin2_preds = preds,
    admin1_preds = admin1_preds,
    polygons     = polygons,
    outcome      = outcome_tag,
    n_model_vars = length(model_vars),
    n_matched    = length(available),
    n_missing    = length(missing)
  )
}


#' Run out-of-sample predictions for all outcomes using pooled area models
#'
#' Trains area-level models on ALL surveyed countries' Admin-2 data combined,
#' then predicts for the target unsurveyed country.
#'
#' @param svy_admin2_list Named list of svy_admin2 data.frames (one per country)
#' @param country_configs Country configs for surveyed countries
#' @param oos_gadm_code GADM ISO3 code for target country
#' @param oos_raster_dir Path to GEE rasters for target country
#' @param oc Outcome config
#' @param params Pipeline parameters
#' @return Output from predict_oos_country()
predict_oos_pooled <- function(svy_admin2_list, country_configs,
                                oos_gadm_code, oos_raster_dir,
                                oc, params) {

  # 1. Pool all surveyed countries' Admin-2 data
  pooled_svy <- dplyr::bind_rows(lapply(names(svy_admin2_list), function(cn) {
    svy <- svy_admin2_list[[cn]]
    if (is.null(svy) || !is.data.frame(svy) || nrow(svy) == 0) return(NULL)
    svy$source_country <- cn
    svy
  }))

  if (nrow(pooled_svy) < 5) {
    warning("Too few pooled Admin2 observations for OOS prediction")
    return(NULL)
  }

  cat(sprintf("[oos_pooled] Pooled %d Admin2 areas from %d countries\n",
              nrow(pooled_svy), length(unique(pooled_svy$source_country))))

  # 2. Extract GEE for all surveyed countries' polygons
  all_gee <- list()
  for (cn in names(country_configs)) {
    cc <- country_configs[[cn]]
    rdir <- cc$raster_dir
    if (!dir.exists(rdir)) {
      alt <- gsub("_GEE_rasters$", "_GEE rasters", rdir)
      if (dir.exists(alt)) rdir <- alt
    }
    if (!dir.exists(rdir)) next

    gadm_raw <- tryCatch(
      geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
      error = function(e) NULL
    )
    if (is.null(gadm_raw)) next
    gadm <- sf::st_as_sf(gadm_raw)
    gadm$Admin2 <- gadm$NAME_2

    tif_files <- sort(list.files(rdir, pattern = "\\.tif$", full.names = TRUE))
    gee_df <- data.frame(Admin2 = gadm$Admin2)
    for (tif in tif_files) {
      varname <- canonicalize_gee_varname(tif)
      r <- tryCatch(terra::rast(tif), error = function(e) NULL)
      if (is.null(r)) next
      if (terra::nlyr(r) > 1) r <- r[[1]]
      vals <- tryCatch(exactextractr::exact_extract(r, gadm, fun = "mean"),
                        error = function(e) rep(NA_real_, nrow(gadm)))
      if (!varname %in% colnames(gee_df)) {
        gee_df[[varname]] <- vals
      }
    }
    all_gee[[cn]] <- gee_df
  }

  # Find common GEE variables across all countries
  gee_var_lists <- lapply(all_gee, function(df) setdiff(colnames(df), "Admin2"))
  common_gee <- Reduce(intersect, gee_var_lists)
  cat(sprintf("[oos_pooled] %d common GEE variables across %d countries\n",
              length(common_gee), length(all_gee)))

  # 3. Build training data: merge survey prevalences with GEE covariates
  train_rows <- list()
  for (cn in names(all_gee)) {
    gee <- all_gee[[cn]]
    svy <- svy_admin2_list[[cn]]
    if (is.null(svy) || nrow(svy) == 0) next
    merged <- dplyr::inner_join(
      svy %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      gee,
      by = "Admin2"
    ) %>% dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))
    if (nrow(merged) > 0) train_rows[[cn]] <- merged
  }
  train_df <- dplyr::bind_rows(train_rows)
  cat(sprintf("[oos_pooled] Training on %d Admin2 areas with survey data\n", nrow(train_df)))

  if (nrow(train_df) < 5) {
    warning("Too few training areas for pooled OOS model")
    return(NULL)
  }

  # 4. Fit pooled glmnet
  valid_vars <- common_gee[sapply(common_gee, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  X_train <- as.matrix(train_df[, valid_vars])
  Y_train <- train_df$svy_prev
  col_medians <- apply(X_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_train))) {
    na_idx <- is.na(X_train[, j])
    if (any(na_idx)) X_train[na_idx, j] <- col_medians[j]
  }

  set.seed(params$seed)
  fit <- tryCatch(
    glmnet::cv.glmnet(x = X_train, y = Y_train, alpha = 0.5,
                       nfolds = max(min(nrow(X_train), 10), 3),
                       type.measure = "mse"),
    error = function(e) { warning("Pooled glmnet failed: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  area_model <- list(fit_area = fit, gee_vars = valid_vars)

  # 5. Extract GEE for target country and predict
  oos_gee <- extract_gee_for_country(oos_gadm_code, oos_raster_dir)
  predict_oos_country(area_model, oos_gee, oc$tag)
}
