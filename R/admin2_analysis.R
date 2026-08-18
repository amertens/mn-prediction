# =============================================================================
# R/admin2_analysis.R
#
# Functions for Admin2-level analysis:
#   - Extract GEE raster zonal means for Admin-2 polygons
#   - Aggregate individual SL predictions to Admin2
#   - Survey-weighted Admin2 prevalence (srvyr)
#   - Area-level model: GEE rasters -> Admin2 prevalence (for unsurveyed areas)
#   - Error analysis
# =============================================================================

# Outcomes whose Admin-2 prevalence is derived from the adjusted continuous
# biomarker + a uniform WHO cutoff (apply_threshold) rather than the survey-
# provided binary column. This makes the cross-country transportability
# outcome COMPARABLE: the survey binaries otherwise mix iron deficiency (ID,
# ferritin only) with iron-deficiency anaemia (IDA, ferritin + anaemia) and
# different inflammation-adjustment methods across (and within) countries.
# Vitamin A is already cutoff-uniform (RBP < 0.70); including it harmonises
# Malawi women (whose survey VAD was unadjusted) onto the adjusted biomarker.
UNIFORM_TRANSPORT_TAGS <- c("child_iron", "women_iron", "child_vitA", "women_vitA")


#' Resolve a country's GEE raster directory (handles space/underscore variants)
#'
#' @param raster_dir Configured raster directory path
#' @return Existing directory path, or NULL if neither variant exists
.resolve_raster_dir <- function(raster_dir) {
  if (dir.exists(raster_dir)) return(raster_dir)
  alt <- gsub("_GEE_rasters$", "_GEE rasters", raster_dir)
  if (dir.exists(alt)) return(alt)
  NULL
}

#' Path of a country's legacy-parity GEE CSV.
#'
#' Keyed on the ISO3 `gadm_code`, NOT on `cc$country`: the config's display name
#' ("Sierra Leone") differs from the config key and the CLI argument
#' ("SierraLeone"), so a name-based path silently misses the file for any
#' multi-word country. Written by scripts/build_gee_legacy_parity.R.
legacy_parity_csv_path <- function(cc) {
  here::here("data", "GEE",
             paste0(cc$gadm_code, "_legacy_parity_admin2_gee.csv"))
}

#' Attach Earth-Engine-extracted admin-2 covariates under the LEGACY gee_* names
#'
#' Fallback for a country that has no local .tif exports in
#' data/<Country>_GEE_rasters/. scripts/build_gee_legacy_parity.R writes
#' data/GEE/<ISO3>_legacy_parity_admin2_gee.csv with exactly the column names
#' .append_gee_zonal_cols() would have produced from rasters, so such a country
#' joins the shared cross-country predictor vocabulary instead of silently
#' contributing nothing (and collapsing every pooled intersection to zero).
#'
#' Values are matched by Admin2 name and returned in `base`'s row order, so the
#' result lines up with the polygon order callers rely on.
#'
#' @param base data.frame with at least an Admin2 column, in polygon order
#' @param cc Country config (needs gadm_code, country)
#' @return `base` with gee_* columns appended, or `base` unchanged if no CSV
.append_legacy_parity_cols <- function(base, cc) {
  csv <- legacy_parity_csv_path(cc)
  if (!file.exists(csv)) return(base)

  parity <- tryCatch(utils::read.csv(csv, check.names = FALSE,
                                     stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(parity) || !"Admin2" %in% names(parity)) {
    warning(sprintf("[gee_parity] unreadable or Admin2-less CSV: %s", csv))
    return(base)
  }
  gee_cols <- grep("^gee_", names(parity), value = TRUE)
  if (!length(gee_cols)) return(base)

  parity$Admin2 <- trimws(as.character(parity$Admin2))
  parity <- parity[!duplicated(parity$Admin2), , drop = FALSE]
  idx <- match(trimws(as.character(base$Admin2)), parity$Admin2)
  for (v in gee_cols) base[[v]] <- parity[[v]][idx]

  cat(sprintf("[gee_parity] %s: %d legacy-named GEE columns from %s (%d/%d areas matched)\n",
              cc$country, length(gee_cols), basename(csv),
              sum(!is.na(idx)), nrow(base)))
  base
}

#' Append GEE raster zonal-mean columns to a base Admin-2 data frame
#'
#' Shared by extract_gee_admin2() (individual-level merge) and
#' extract_area_covariates() (area-level model) so the extraction logic — and
#' therefore the gee_* column naming — stays identical across both. `gee_admin2`
#' must have one row per feature in `all_polys`, in the same order; gee_* columns
#' are appended in place. Distinct-variable multi-band rasters (FLDAS, Soil, ESA,
#' 12-month climatologies) with up to `max_bands_individual` bands are kept
#' individually plus annual summaries. Big temporal stacks above that threshold
#' (e.g. Atmosphere = 628 daily bands, GPW = 77 demographic bins, WAPOR = 36
#' dekadal slices) are COLLAPSED across bands with terra FIRST and then extracted
#' once, mirroring extract_gee_cluster_buffers(): one mean/sd layer apiece instead
#' of one exact_extract per band. The old per-band path looped exact_extract once
#' PER BAND, which segfaults on huge stacks and spawns hundreds of bloat columns.
#'
#' @param gee_admin2 Base data frame (one row per polygon, in polygon order)
#' @param all_polys sf polygons passed to exactextractr::exact_extract
#' @param raster_dir Resolved directory containing the country's .tif rasters
#' @param max_bands_individual Bands at/below this are kept individually; above,
#'   the stack is collapsed across bands with terra before a single extract
#' @return `gee_admin2` with gee_* columns appended
.append_gee_zonal_cols <- function(gee_admin2, all_polys, raster_dir,
                                    max_bands_individual = 24L) {
  tif_files <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  cat(sprintf("[gee_zonal] %d .tif files, %d polygons\n",
              length(tif_files), nrow(all_polys)))

  for (tif in tif_files) {
    base <- tools::file_path_sans_ext(basename(tif))
    # Strip country name for cross-country variable-name consistency
    base_clean <- base
    for (cname in c("Gambia", "Ghana", "Sierra_Leone", "Sierra Leone",
                    "Malawi", "Cote_dIvoire", "Cote d'Ivoire")) {
      base_clean <- gsub(paste0("_?", gsub("[' ]", ".", cname), "_?"), "_", base_clean)
    }
    base_varname <- paste0("gee_", tolower(gsub("[^A-Za-z0-9]+", "_", base_clean)))
    base_varname <- gsub("_+", "_", base_varname)
    base_varname <- sub("^_|_$", "", base_varname)

    r <- tryCatch(terra::rast(tif), error = function(e) NULL)
    if (is.null(r)) next

    n_layers <- terra::nlyr(r)

    if (n_layers == 1) {
      vals <- tryCatch(
        exactextractr::exact_extract(r, all_polys, fun = "mean"),
        error = function(e) rep(NA_real_, nrow(all_polys))
      )
      gee_admin2[[base_varname]] <- vals
    } else if (n_layers <= max_bands_individual) {
      # Distinct-variable / small temporal stack: keep each band individually,
      # then add annual summaries. Extract every band once into a matrix so the
      # summaries reuse those values instead of re-extracting per band.
      layer_names <- tryCatch(names(r), error = function(e) NULL)
      all_layer_vals <- tryCatch(
        vapply(seq_len(n_layers),
               function(i) exactextractr::exact_extract(r[[i]], all_polys, fun = "mean"),
               numeric(nrow(all_polys))),
        error = function(e) NULL
      )

      if (!is.null(all_layer_vals) && is.matrix(all_layer_vals)) {
        for (lyr_idx in seq_len(n_layers)) {
          lyr_name <- if (!is.null(layer_names) && nchar(layer_names[lyr_idx]) > 0) {
            ln <- tolower(gsub("[^A-Za-z0-9]+", "_", layer_names[lyr_idx]))
            paste0(base_varname, "_", ln)
          } else {
            month_label <- if (n_layers == 12) tolower(month.abb[lyr_idx])
                           else sprintf("b%02d", lyr_idx)
            paste0(base_varname, "_", month_label)
          }
          lyr_name <- gsub("_+", "_", lyr_name)
          lyr_name <- sub("_$", "", lyr_name)
          gee_admin2[[lyr_name]] <- all_layer_vals[, lyr_idx]
        }

        gee_admin2[[paste0(base_varname, "_annual_mean")]] <- rowMeans(all_layer_vals, na.rm = TRUE)
        gee_admin2[[paste0(base_varname, "_annual_sd")]]   <- apply(all_layer_vals, 1, sd, na.rm = TRUE)
        gee_admin2[[paste0(base_varname, "_annual_min")]]  <- apply(all_layer_vals, 1, min, na.rm = TRUE)
        gee_admin2[[paste0(base_varname, "_annual_max")]]  <- apply(all_layer_vals, 1, max, na.rm = TRUE)
        gee_admin2[[paste0(base_varname, "_annual_range")]] <- gee_admin2[[paste0(base_varname, "_annual_max")]] -
                                                                gee_admin2[[paste0(base_varname, "_annual_min")]]
      }
    } else {
      # Big temporal stack (Atmosphere 628 bands, GPW 77, WAPOR 36): collapse
      # across bands with terra FIRST, then exact_extract the single collapsed
      # layer ONCE. Avoids the per-band exact_extract blow-up (which segfaults on
      # Atmosphere) and the hundreds of country-specific per-band bloat columns.
      cat(sprintf("    [collapse] %s (%d bands -> annual mean+sd)\n", base_varname, n_layers))
      r_mean <- tryCatch(terra::app(r, fun = mean, na.rm = TRUE), error = function(e) NULL)
      if (!is.null(r_mean)) {
        gee_admin2[[paste0(base_varname, "_annual_mean")]] <- tryCatch(
          exactextractr::exact_extract(r_mean, all_polys, fun = "mean"),
          error = function(e) rep(NA_real_, nrow(all_polys))
        )
      }
      r_sd <- tryCatch(terra::app(r, fun = sd, na.rm = TRUE), error = function(e) NULL)
      if (!is.null(r_sd)) {
        gee_admin2[[paste0(base_varname, "_annual_sd")]] <- tryCatch(
          exactextractr::exact_extract(r_sd, all_polys, fun = "mean"),
          error = function(e) rep(NA_real_, nrow(all_polys))
        )
      }
    }
  }

  gee_admin2
}

#' Extract GEE raster zonal means for all Admin-2 polygons in a country
#'
#' This is a shared extraction step used by both the individual-level model
#' (merged into individual records via Admin-2 join) and the area-level model.
#' Multi-band rasters with >20 temporal/demographic bands (GPW, WAPOR) are
#' reduced to summary statistics; distinct-variable bands (FLDAS, Soil) are
#' kept individually.
#'
#' @param cc Country config (needs gadm_code, raster_dir)
#' @return data.frame with Admin1, Admin2, and gee_* columns (one row per Admin-2)
extract_gee_admin2 <- function(cc) {

  gadm_raw <- tryCatch(
    geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
    error = function(e) {
      tryCatch(
        geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(gadm_raw)) {
    warning(sprintf("[extract_gee_admin2] GADM download failed for %s", cc$gadm_code))
    return(NULL)
  }
  all_polys <- sf::st_as_sf(gadm_raw)
  all_polys$Admin2 <- all_polys$NAME_2
  all_polys$Admin1 <- all_polys$NAME_1

  gee_admin2 <- data.frame(Admin1 = all_polys$Admin1, Admin2 = all_polys$Admin2)
  raster_dir <- .resolve_raster_dir(cc$raster_dir)

  if (is.null(raster_dir)) {
    # No local rasters — fall back to the Earth-Engine legacy-parity CSV, which
    # carries the same column names the raster path would have produced.
    cat(sprintf("[extract_gee_admin2] No raster dir for %s — trying legacy-parity CSV\n",
                cc$country))
    gee_admin2 <- .append_legacy_parity_cols(gee_admin2, cc)
    if (ncol(gee_admin2) == 2L)
      warning(sprintf(paste0("[extract_gee_admin2] %s has neither GEE rasters nor a ",
                             "legacy-parity CSV: it will contribute NO GEE predictors ",
                             "and will empty the pooled/LOCO GEE intersection. Run ",
                             "scripts/build_gee_legacy_parity.R %s."),
                      cc$country, cc$country))
    return(gee_admin2)
  }

  cat(sprintf("[extract_gee_admin2] %s: %d Admin-2 polygons\n",
              cc$country, nrow(all_polys)))

  gee_admin2 <- .append_gee_zonal_cols(gee_admin2, all_polys, raster_dir)

  n_gee <- ncol(gee_admin2) - 2  # minus Admin1, Admin2
  cat(sprintf("[extract_gee_admin2] %s: %d GEE variables extracted\n", cc$country, n_gee))

  # Deduplicate Admin2 names (e.g., "Lake Malawi" appears multiple times)
  dup_names <- gee_admin2$Admin2[duplicated(gee_admin2$Admin2)]
  if (length(dup_names) > 0) {
    cat(sprintf("[extract_gee_admin2] Deduplicating %d duplicate Admin2 names\n",
                length(unique(dup_names))))
    num_cols <- names(gee_admin2)[sapply(gee_admin2, is.numeric)]
    gee_admin2 <- gee_admin2 |>
      dplyr::group_by(Admin2) |>
      dplyr::summarise(
        Admin1 = dplyr::first(Admin1),
        dplyr::across(dplyr::all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) |>
      as.data.frame()
  }

  gee_admin2
}

#' Extract per-country area-level covariates (GADM polygons + GEE zonal means)
#'
#' Runs the expensive GADM download + raster zonal extraction ONCE per country
#' so all of that country's area_model_<outcome> targets can reuse it — the
#' extraction depends only on the country, not the outcome. Returns both the
#' gee_* covariate table (one row per Admin-2 polygon, in polygon order and NOT
#' deduplicated, so it lines up with the area model's join/prediction logic) and
#' the sf polygons (needed for the coverage map). Column naming is shared with
#' extract_gee_admin2() via .append_gee_zonal_cols(), so downstream predictor
#' screening and outputs are unchanged.
#'
#' @param cc Country config (needs gadm_code, raster_dir)
#' @return list(gee_admin2 = data.frame, polygons = sf); both NULL if GADM fails
extract_area_covariates <- function(cc) {
  gadm_raw <- tryCatch(
    geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
    error = function(e) {
      # Retry once — GADM server can be flaky
      tryCatch(
        geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(gadm_raw)) {
    warning(sprintf("[extract_area_covariates] GADM download failed for %s", cc$gadm_code))
    return(list(gee_admin2 = NULL, polygons = NULL))
  }
  all_polys <- sf::st_as_sf(gadm_raw)
  all_polys$Admin2 <- all_polys$NAME_2

  gee_admin2 <- data.frame(Admin2 = all_polys$Admin2)
  raster_dir <- .resolve_raster_dir(cc$raster_dir)

  if (is.null(raster_dir)) {
    # Same fallback as extract_gee_admin2(): Earth-Engine covariates under the
    # legacy column names, matched by Admin2 and returned in polygon order.
    gee_admin2 <- .append_legacy_parity_cols(gee_admin2, cc)
    if (ncol(gee_admin2) == 1L)
      warning(sprintf(paste0("[extract_area_covariates] no GEE rasters (%s) and no ",
                             "legacy-parity CSV for %s — the area model has no ",
                             "covariates. Run scripts/build_gee_legacy_parity.R %s."),
                      cc$raster_dir, cc$country, cc$country))
    return(list(gee_admin2 = gee_admin2, polygons = all_polys))
  }

  gee_admin2 <- .append_gee_zonal_cols(gee_admin2, all_polys, raster_dir)
  cat(sprintf("[extract_area_covariates] %s: %d GEE variables extracted\n",
              cc$country, ncol(gee_admin2) - 1))

  list(gee_admin2 = gee_admin2, polygons = all_polys)
}

#' Aggregate individual-level SL predictions to Admin2 prevalence
#'
#' @param sl_fit Output from fit_sl_models()
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with Admin2, n, sl_prev
aggregate_admin2_sl <- function(sl_fit, outcome_data, cc, oc) {
  if (is.null(sl_fit)) return(NULL)

  fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
  if (is.null(fit)) return(NULL)

  d   <- outcome_data$data
  res <- fit$res

  # Predicted prevalence: binary model gives P(deficient), continuous uses threshold
  if (!is.null(sl_fit$bin_fit)) {
    d$deficient <- res$yhat_full
  } else {
    d$deficient <- apply_threshold(res$yhat_full, oc$cutoff, oc$cutoff_dir)
  }

  admin2_col <- cc$admin2_col

  # Use survey weights for predicted prevalence aggregation
  wt_col <- cc$weight_col
  has_weights <- !is.null(wt_col) && wt_col %in% colnames(d) &&
                 any(!is.na(d[[wt_col]]) & d[[wt_col]] > 0)
  if (has_weights) {
    d$`.wt` <- as.numeric(d[[wt_col]])
    d$`.wt`[is.na(d$`.wt`) | d$`.wt` <= 0] <- 1
  } else {
    d$`.wt` <- 1
  }

  d %>%
    dplyr::filter(!is.na(.data[[admin2_col]])) %>%
    dplyr::group_by(.data[[admin2_col]]) %>%
    dplyr::summarise(
      n_sl    = dplyr::n(),
      sl_prev = stats::weighted.mean(deficient, `.wt`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(Admin2 = dplyr::all_of(admin2_col))
}


#' Compute survey-weighted Admin2 prevalence
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with Admin2, svy_prev, svy_prev_se, n_svy, etc.
compute_svy_admin2 <- function(outcome_data, cc, oc) {
  d <- outcome_data$data

  # ── Uniform, cross-country-comparable outcome ───────────────────────────
  # For the harmonized transportability outcomes, overwrite the survey binary
  # with one derived from the adjusted continuous biomarker + a uniform WHO
  # cutoff, so every country uses the SAME definition (e.g. ferritin < 12/15 =
  # ID, RBP < 0.70 = VAD). Falls back to the survey binary if the continuous
  # column is unavailable. Downstream survey weighting is unchanged.
  if (!is.null(oc$tag) && oc$tag %in% UNIFORM_TRANSPORT_TAGS && !is.null(oc$binary)) {
    derived <- NULL
    if (grepl("vitA", oc$tag, ignore.case = TRUE)) {
      # 2026-06-23 (DC-H2): uniform BRINDA inflammation adjustment of RBP, so every
      # country's VAD outcome uses ONE method (R/brinda_adjustment.R), validated in
      # docs/dc_h2_brinda_validation.md. Replaces the prior Thurnham/raw RBP mix.
      # brinda_vad_binary() is the single source of truth, shared with
      # apply_brinda_vita_binary(); it warns and returns NULL when a country's
      # biomarker columns are unavailable, in which case the configured binary
      # is kept.
      newbin <- brinda_vad_binary(d, cc, oc, label = "[svy_admin2]")
      if (!is.null(newbin)) derived <- as.numeric(newbin)
    } else if (!is.null(oc$continuous) && oc$continuous %in% colnames(d) && !is.null(oc$cutoff)) {
      # Non-VitA uniform outcomes (iron): threshold the configured continuous.
      derived <- apply_threshold(suppressWarnings(as.numeric(d[[oc$continuous]])),
                                 oc$cutoff, oc$cutoff_dir %||% "less")
      cat(sprintf("  [svy_admin2] %s — %s: uniform outcome %s %s %g => %d/%d (%.1f%%) deficient\n",
                  cc$country, oc$tag, oc$continuous, oc$cutoff_dir %||% "less", oc$cutoff,
                  sum(derived == 1, na.rm = TRUE), sum(!is.na(derived)), 100 * mean(derived, na.rm = TRUE)))
    }
    if (!is.null(derived)) { d[[oc$binary]] <- derived; outcome_data$data <- d }
  }

  svy_cols <- c(cc$psu_col, cc$weight_col, cc$admin2_col,
                cc$admin1_col, oc$binary)
  if (!is.null(cc$strata_col)) svy_cols <- c(svy_cols, cc$strata_col)
  svy_cols <- unique(svy_cols[svy_cols %in% colnames(d)])

  svy_df <- d %>%
    dplyr::select(dplyr::all_of(svy_cols)) %>%
    dplyr::filter(!is.na(.data[[oc$binary]]),
                  !is.na(.data[[cc$admin2_col]]),
                  !is.na(.data[[cc$weight_col]])) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(cc$psu_col), as.factor),
      !!oc$binary := as.numeric(.data[[oc$binary]]),
      !!cc$weight_col := as.numeric(.data[[cc$weight_col]])
    )

  if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(svy_df))
    svy_df <- svy_df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(cc$strata_col), as.factor))

  options(survey.lonely.psu = "adjust")

  if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(svy_df)) {
    svy_des <- srvyr::as_survey_design(
      svy_df,
      ids     = !!rlang::sym(cc$psu_col),
      strata  = !!rlang::sym(cc$strata_col),
      weights = !!rlang::sym(cc$weight_col),
      nest    = TRUE
    )
  } else {
    svy_des <- srvyr::as_survey_design(
      svy_df,
      ids     = !!rlang::sym(cc$psu_col),
      weights = !!rlang::sym(cc$weight_col)
    )
  }

  # Some Admin2 areas have too few clusters for survey_mean to compute SEs.
  # Compute per-group safely: first get simple counts, then try survey_mean.
  svy_admin2 <- tryCatch({
    svy_des %>%
      dplyr::group_by(.data[[cc$admin2_col]]) %>%
      dplyr::summarise(
        svy_prev = srvyr::survey_mean(!!rlang::sym(oc$binary),
                                       vartype = c("se", "ci"), na.rm = TRUE),
        n_svy    = dplyr::n(),
        .groups  = "drop"
      ) %>%
      dplyr::rename(Admin2 = dplyr::all_of(cc$admin2_col)) %>%
      dplyr::mutate(
        svy_cv       = dplyr::if_else(svy_prev > 0, svy_prev_se / svy_prev, NA_real_),
        svy_ci_width = svy_prev_upp - svy_prev_low
      )
  }, error = function(e) {
    # Fallback: compute per Admin2 area individually, skipping failures
    cat(sprintf("  [svy_admin2] Grouped survey_mean failed (%s), computing per-area...\n",
                conditionMessage(e)))
    admin2_vals <- unique(svy_df[[cc$admin2_col]])
    admin2_vals <- admin2_vals[!is.na(admin2_vals)]

    results <- lapply(admin2_vals, function(a2) {
      sub <- svy_df[svy_df[[cc$admin2_col]] == a2, ]
      n_clusters <- length(unique(sub[[cc$psu_col]]))
      n_obs <- nrow(sub)

      if (n_obs < 2) return(NULL)

      res <- tryCatch({
        if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(sub)) {
          des <- srvyr::as_survey_design(sub,
            ids = !!rlang::sym(cc$psu_col),
            strata = !!rlang::sym(cc$strata_col),
            weights = !!rlang::sym(cc$weight_col), nest = TRUE)
        } else {
          des <- srvyr::as_survey_design(sub,
            ids = !!rlang::sym(cc$psu_col),
            weights = !!rlang::sym(cc$weight_col))
        }
        sm <- des %>%
          dplyr::summarise(
            svy_prev = srvyr::survey_mean(!!rlang::sym(oc$binary),
                                           vartype = c("se", "ci"), na.rm = TRUE),
            .groups = "drop")
        data.frame(Admin2 = a2, svy_prev = sm$svy_prev,
                   svy_prev_se = sm$svy_prev_se,
                   svy_prev_low = sm$svy_prev_low,
                   svy_prev_upp = sm$svy_prev_upp,
                   n_svy = n_obs, stringsAsFactors = FALSE)
      }, error = function(e2) {
        # Last resort: unweighted mean, no SE
        prev <- mean(sub[[oc$binary]], na.rm = TRUE)
        data.frame(Admin2 = a2, svy_prev = prev,
                   svy_prev_se = NA_real_,
                   svy_prev_low = NA_real_,
                   svy_prev_upp = NA_real_,
                   n_svy = n_obs, stringsAsFactors = FALSE)
      })
      res
    })

    valid_results <- results[!sapply(results, is.null)]
    if (length(valid_results) == 0) {
      return(data.frame(Admin2 = character(), svy_prev = numeric(),
                        svy_prev_se = numeric(), svy_prev_low = numeric(),
                        svy_prev_upp = numeric(), n_svy = integer(),
                        svy_cv = numeric(), svy_ci_width = numeric(),
                        stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, valid_results)
    out$svy_cv <- ifelse(out$svy_prev > 0 & !is.na(out$svy_prev_se),
                         out$svy_prev_se / out$svy_prev, NA_real_)
    out$svy_ci_width <- out$svy_prev_upp - out$svy_prev_low
    out
  })

  svy_admin2
}


#' Fit area-level model (GEE rasters -> Admin2 prevalence) and predict everywhere
#'
#' This trains an elastic net on surveyed Admin2 units and predicts to ALL
#' polygons (including unsurveyed areas with no biomarker data).
#'
#' @param svy_admin2 Survey-weighted Admin2 prevalence (from compute_svy_admin2)
#' @param area_cov Pre-extracted per-country covariates from
#'   extract_area_covariates(): list(gee_admin2 = data.frame, polygons = sf).
#'   Shared across all of a country's outcomes so the expensive GADM download +
#'   raster zonal extraction runs only once per country.
#' @param cc Country config
#' @param oc Outcome config
#' @param params Pipeline parameters
#' @return list with area_preds, loo_summary, fit_area (glmnet object)
fit_area_level_model <- function(svy_admin2, area_cov, cc, oc, params) {
  if (is.null(svy_admin2) || !is.data.frame(svy_admin2) ||
      nrow(svy_admin2) < 3 || all(is.na(svy_admin2$svy_prev))) {
    warning("Too few survey Admin2 estimates for area-level model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = NULL))
  }

  # Per-country GADM polygons + GEE raster zonal means are pre-extracted ONCE
  # per country by the area_covariates_<country> target (extract_area_covariates)
  # and shared across all of that country's outcomes — the extraction depends
  # only on the country, not the outcome.
  if (is.null(area_cov) || is.null(area_cov$polygons)) {
    warning(sprintf("Area covariates unavailable for %s — skipping area model",
                    cc$gadm_code))
    return(list(area_preds = NULL, loo_summary = NULL, polygons = NULL))
  }
  all_polys  <- area_cov$polygons
  gee_admin2 <- area_cov$gee_admin2

  # Identify surveyed/unsurveyed
  surveyed_names <- svy_admin2$Admin2
  n_surveyed     <- length(surveyed_names)
  n_unsurveyed   <- sum(!all_polys$Admin2 %in% surveyed_names)

  cat(sprintf("[area_model] %d surveyed, %d unsurveyed Admin2 units\n",
              n_surveyed, n_unsurveyed))

  gee_vars <- setdiff(colnames(gee_admin2), c("Admin1", "Admin2"))
  if (length(gee_vars) == 0) {
    warning("No GEE covariates available for area-level model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }
  cat(sprintf("  %d GEE variables available\n", length(gee_vars)))

  # Training set: surveyed areas with GEE covariates

  train_df <- gee_admin2 %>%
    dplyr::inner_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    ) %>%
    dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))

  if (nrow(train_df) < 3) {
    warning("Too few surveyed Admin2 areas with valid prevalence for area model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }

  # Drop zero-variance or all-NA predictors
  valid_vars <- gee_vars[sapply(gee_vars, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  if (length(valid_vars) < 2) {
    warning("Too few valid GEE predictors for area-level model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }

  X_train <- as.matrix(train_df[, valid_vars])
  Y_train <- train_df$svy_prev
  col_medians <- apply(X_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_train))) {
    na_idx <- is.na(X_train[, j])
    if (any(na_idx)) X_train[na_idx, j] <- col_medians[j]
  }

  X_all <- as.matrix(gee_admin2[, valid_vars])
  for (j in seq_len(ncol(X_all))) {
    na_idx <- is.na(X_all[, j])
    if (any(na_idx)) X_all[na_idx, j] <- col_medians[j]
  }

  # ── HAL (Highly Adaptive Lasso) for area-level model ──────────────────
  # HAL is a nonparametric estimator that:
  #   - Converges at n^{-1/3} rate (faster than most nonparametric methods)
  #   - Has built-in L1 regularization via the lasso penalty
  #   - Discovers interactions automatically via indicator basis functions
  #   - Provides valid semiparametric inference
  #   - Works well in small-n settings (n=30-100 Admin2 areas)
  #
  # We use max_degree=2 (main effects + pairwise interactions) and
  # smoothness_orders=0 (piecewise constant, appropriate for geographic data).
  # The number of knots is kept small to control computation.

  # HAL parameters — conservative for small n
  hal_max_degree   <- if (nrow(X_train) < 20) 1L else 2L
  hal_num_knots    <- if (nrow(X_train) < 30) c(3, 1) else c(5, 2)
  if (hal_max_degree == 1L) hal_num_knots <- hal_num_knots[1]

  # Helper: fit HAL with fallback to glmnet
  fit_hal_safe <- function(X, Y, ...) {
    tryCatch({
      hal9001::fit_hal(
        X = X, Y = Y,
        max_degree       = hal_max_degree,
        num_knots        = hal_num_knots,
        smoothness_orders = 0L,
        family           = "gaussian",
        ...
      )
    }, error = function(e) {
      cat(sprintf("  HAL failed (%s), falling back to glmnet\n", e$message))
      NULL
    })
  }

  # Helper: predict from HAL or glmnet
  predict_hal_safe <- function(fit, newdata) {
    if (inherits(fit, "hal9001")) {
      as.numeric(predict(fit, new_data = newdata))
    } else {
      # glmnet fallback
      as.numeric(predict(fit, newdata, s = "lambda.min"))
    }
  }

  # ── Screen predictors for the area-level model ────────────────────────
  # HAL with max_degree=2 builds pairwise-interaction basis functions over
  # every predictor, so both runtime and memory explode with the number of
  # columns. With only n ≈ 30–100 surveyed Admin2 areas, a large covariate
  # set is also statistically unidentifiable. Cap the predictor set to the
  # top-k columns most correlated (marginally) with the outcome. Screening is
  # redone inside each LOO fold so the cross-validated metrics stay honest.
  AREA_MAX_PREDICTORS <- 20L
  screen_top_k <- function(Xm, Yv, k = AREA_MAX_PREDICTORS) {
    if (ncol(Xm) <= k) return(seq_len(ncol(Xm)))
    score <- vapply(seq_len(ncol(Xm)), function(j) {
      x <- Xm[, j]
      if (length(unique(x)) < 2) return(0)
      v <- abs(suppressWarnings(stats::cor(x, Yv)))
      if (is.na(v)) 0 else v
    }, numeric(1))
    order(score, decreasing = TRUE)[seq_len(k)]
  }

  keep_final <- screen_top_k(X_train, Y_train)
  X_train_s  <- X_train[, keep_final, drop = FALSE]
  X_all_s    <- X_all[,   keep_final, drop = FALSE]
  cat(sprintf("  Screened %d -> %d area-level predictors (top-k by |cor|)\n",
              length(valid_vars), length(keep_final)))

  # LOO-CV (screening repeated within each fold to avoid leakage)
  cat(sprintf("  Running LOO-CV (HAL, max_degree=%d)...\n", hal_max_degree))
  loo_preds <- rep(NA_real_, nrow(X_train))
  for (i in seq_len(nrow(X_train))) {
    n_remain <- nrow(X_train) - 1
    if (n_remain < 3) next
    keep_i <- screen_top_k(X_train[-i, , drop = FALSE], Y_train[-i])
    X_in   <- X_train[-i, keep_i, drop = FALSE]
    fit_loo <- fit_hal_safe(X_in, Y_train[-i])
    # Fallback to glmnet if HAL fails
    if (is.null(fit_loo)) {
      fit_loo <- tryCatch(
        glmnet::cv.glmnet(
          x = X_in, y = Y_train[-i],
          alpha = 0.5, nfolds = min(n_remain, 10), type.measure = "mse"
        ), error = function(e) NULL
      )
    }
    if (!is.null(fit_loo))
      loo_preds[i] <- predict_hal_safe(fit_loo, X_train[i, keep_i, drop = FALSE])
  }
  loo_preds <- pmin(pmax(loo_preds, 0), 1)

  loo_valid <- !is.na(loo_preds)
  loo_summary <- if (sum(loo_valid) >= 3) {
    data.frame(
      outcome      = oc$tag,
      country      = cc$country,
      n_train      = sum(loo_valid),
      n_predictors = length(keep_final),
      model_type   = "HAL",
      loo_mae_pp   = round(mean(abs(Y_train[loo_valid] - loo_preds[loo_valid])) * 100, 2),
      loo_rmse_pp  = round(sqrt(mean((Y_train[loo_valid] - loo_preds[loo_valid])^2)) * 100, 2),
      loo_r        = round(cor(Y_train[loo_valid], loo_preds[loo_valid]), 3)
    )
  } else NULL

  if (!is.null(loo_summary))
    cat(sprintf("  LOO-CV: MAE = %.1f pp, r = %.3f\n",
                loo_summary$loo_mae_pp, loo_summary$loo_r))

  # Final model on all training data
  set.seed(params$seed)
  fit_area <- fit_hal_safe(X_train_s, Y_train)
  # Fallback to glmnet if HAL fails on full data
  if (is.null(fit_area)) {
    fit_area <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train_s, y = Y_train, alpha = 0.5,
        nfolds = max(min(nrow(X_train_s), 10), 3), type.measure = "mse"
      ),
      error = function(e) {
        warning("Both HAL and glmnet failed for area model: ", e$message)
        NULL
      }
    )
  }

  if (is.null(fit_area)) {
    return(list(area_preds = NULL, loo_summary = loo_summary, polygons = all_polys))
  }

  yhat_all <- pmin(pmax(predict_hal_safe(fit_area, X_all_s), 0), 1)

  area_preds <- data.frame(
    Admin2         = gee_admin2$Admin2,
    area_pred_prev = yhat_all,
    has_survey     = gee_admin2$Admin2 %in% surveyed_names
  ) %>%
    dplyr::left_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    )

  # Bootstrap CIs
  B_AREA <- params$B_area
  cat(sprintf("  Bootstrapping area model (B = %d)...\n", B_AREA))
  set.seed(params$seed)
  boot_all <- matrix(NA_real_, nrow = nrow(X_all), ncol = B_AREA)

  for (b in seq_len(B_AREA)) {
    idx_b <- sample(nrow(X_train_s), replace = TRUE)
    fit_b <- fit_hal_safe(X_train_s[idx_b, , drop = FALSE], Y_train[idx_b])
    if (is.null(fit_b)) {
      # glmnet fallback for this bootstrap rep
      fit_b <- tryCatch(
        glmnet::cv.glmnet(
          x = X_train_s[idx_b, , drop = FALSE], y = Y_train[idx_b],
          alpha = 0.5, nfolds = min(length(unique(idx_b)), 10),
          type.measure = "mse"
        ), error = function(e) NULL
      )
    }
    if (!is.null(fit_b))
      boot_all[, b] <- pmin(pmax(predict_hal_safe(fit_b, X_all_s), 0), 1)
  }

  valid_boots <- sum(!is.na(boot_all[1, ]))
  if (valid_boots >= 3) {
    area_preds$boot_mean <- apply(boot_all, 1, mean, na.rm = TRUE)
    area_preds$ci_lo     <- apply(boot_all, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
    area_preds$ci_hi     <- apply(boot_all, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
    area_preds$ci_width  <- area_preds$ci_hi - area_preds$ci_lo
  }

  list(
    area_preds  = area_preds,
    loo_summary = loo_summary,
    polygons    = all_polys,
    fit_area    = fit_area,
    gee_vars    = valid_vars[keep_final]
  )
}


#' Compute error metrics between SL predictions and survey estimates
#'
#' @param sl_admin2 Admin2 SL predictions (from aggregate_admin2_sl)
#' @param svy_admin2 Survey-weighted Admin2 prevalence (from compute_svy_admin2)
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with error metrics, or NULL if inputs are missing
compute_admin2_error <- function(sl_admin2, svy_admin2, cc, oc) {
  if (is.null(sl_admin2) || is.null(svy_admin2)) return(NULL)

  merged <- sl_admin2 %>%
    dplyr::inner_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    ) %>%
    dplyr::filter(!is.na(sl_prev), !is.na(svy_prev))

  if (nrow(merged) < 2) return(NULL)

  data.frame(
    outcome      = oc$tag,
    country      = cc$country,
    n_admin2     = nrow(merged),
    mae_pp       = round(mean(abs(merged$sl_prev - merged$svy_prev)) * 100, 2),
    rmse_pp      = round(sqrt(mean((merged$sl_prev - merged$svy_prev)^2)) * 100, 2),
    pearson_r    = round(cor(merged$sl_prev, merged$svy_prev), 3),
    mean_bias_pp = round(mean(merged$sl_prev - merged$svy_prev) * 100, 2)
  )
}
