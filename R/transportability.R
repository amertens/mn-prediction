# =============================================================================
# R/transportability.R
#
# Cross-country transportability analysis via leave-one-country-out (LOCO) CV.
#
# Design: Pool data from all countries using only covariates available in every
# country (IHME, MAP, GEE). For each held-out country, train a SuperLearner on
# the remaining countries and evaluate on the held-out one. This tests whether
# proxy-based deficiency models generalize across populations.
# =============================================================================

# Prefixes of domains expected to be present in ALL countries.
# DHS, MICS, LSMS, WFP, FluNet vary by country and are excluded.
# FSEC (food security) is included because HFID covers all countries;
# CH-only columns will be dropped by the intersection step if not universal.
COMMON_DOMAIN_PREFIXES <- c("ihme_", "MAP_", "gee_", "fsec_")


#' Build a pooled multi-country dataset for one outcome
#'
#' Loads per-country merged data, filters to the correct population, harmonizes
#' outcome column names, finds the intersection of common proxy predictors, and
#' stacks into a single dataset.
#'
#' @param all_merged Named list of merged data.frames (one per country)
#' @param all_configs Output of get_country_configs()
#' @param outcome_tag Character, e.g. "child_vitA"
#' @return list with: data, Xvars_common, outcome_tag
build_pooled_dataset <- function(all_merged, all_configs, outcome_tag) {

  country_frames <- list()
  country_proxy_cols <- list()

  for (cn in names(all_merged)) {
    cc <- all_configs[[cn]]
    oc <- cc$outcomes[[outcome_tag]]
    if (is.null(oc)) {
      cat(sprintf("  [pool] Skipping %s — outcome %s not defined\n", cn, outcome_tag))
      next
    }

    d <- all_merged[[cn]]

    # Filter to population
    pop_col <- cc$child_flag
    if (!is.null(pop_col) && pop_col %in% colnames(d)) {
      d <- d[d[[pop_col]] == oc$child_flag_val, ]
    }

    # Require non-missing binary outcome
    bin_col <- oc$binary
    if (is.null(bin_col) || !(bin_col %in% colnames(d))) {
      cat(sprintf("  [pool] Skipping %s — binary column %s not found\n", cn, bin_col))
      next
    }
    d <- d[!is.na(d[[bin_col]]), ]

    if (nrow(d) == 0) {
      cat(sprintf("  [pool] Skipping %s — 0 rows after filtering\n", cn))
      next
    }

    # Identify proxy columns from common domains only
    all_cols <- colnames(d)
    proxy_cols <- character()
    for (pfx in COMMON_DOMAIN_PREFIXES) {
      proxy_cols <- c(proxy_cols, all_cols[startsWith(all_cols, pfx)])
    }
    country_proxy_cols[[cn]] <- proxy_cols

    # Harmonize: rename binary outcome to Y_binary
    # Handle 1/2 factor coding (1=deficient, 2=not) same as build_outcome_dataset
    yb <- d[[bin_col]]
    yb_vals <- sort(unique(as.numeric(yb[!is.na(yb)])))
    if (length(yb_vals) == 2 && all(yb_vals == c(1, 2))) {
      cat(sprintf("  [pool] Recoding %s from {1,2} to {1,0}\n", bin_col))
      d[[bin_col]] <- ifelse(as.numeric(d[[bin_col]]) == 1, 1L, 0L)
    }
    # Ensure values are in {0, 1}
    d$Y_binary <- as.integer(d[[bin_col]])
    bad <- !d$Y_binary %in% c(0L, 1L) & !is.na(d$Y_binary)
    if (any(bad)) {
      cat(sprintf("  [pool] WARNING: %d values of %s not in {0,1} — setting to NA\n",
                  sum(bad), bin_col))
      d$Y_binary[bad] <- NA_integer_
    }

    # Add country identifier and make cluster IDs globally unique
    d$country <- cn
    cluster_col <- cc$cluster_id
    if (cluster_col %in% colnames(d)) {
      d$pooled_cluster <- paste0(tolower(cn), "_", d[[cluster_col]])
    } else {
      d$pooled_cluster <- paste0(tolower(cn), "_", seq_len(nrow(d)))
    }

    # Carry Admin1 for aggregation (rename if needed)
    if (cc$admin1_col %in% colnames(d)) {
      d$Admin1 <- d[[cc$admin1_col]]
    }

    country_frames[[cn]] <- d
  }

  if (length(country_frames) < 2) {
    stop("[pool] Need at least 2 countries with valid data for transportability")
  }

  # Find intersection of proxy predictor columns across ALL countries
  Xvars_common <- Reduce(intersect, country_proxy_cols)
  cat(sprintf("[pool] %s: %d common proxy predictors across %d countries\n",
              outcome_tag, length(Xvars_common), length(country_frames)))

  if (length(Xvars_common) < 5) {
    warning(sprintf("Only %d common predictors — model may be unreliable",
                    length(Xvars_common)))
  }

  # Select common columns and stack
  keep_cols <- unique(c("Y_binary", "country", "pooled_cluster", "Admin1",
                        Xvars_common))

  stacked <- dplyr::bind_rows(lapply(country_frames, function(d) {
    d_sub <- d[, intersect(keep_cols, colnames(d)), drop = FALSE]
    # Ensure all columns are atomic (haven labels, list columns can cause
    # is.nan() errors in sl3). Convert factors/labelled to numeric.
    for (col in colnames(d_sub)) {
      v <- d_sub[[col]]
      if (is.list(v)) {
        d_sub[[col]] <- as.numeric(unlist(v))
      } else if (inherits(v, "haven_labelled")) {
        d_sub[[col]] <- as.numeric(v)
      } else if (is.factor(v)) {
        d_sub[[col]] <- as.numeric(as.character(v))
      }
    }
    d_sub
  }))

  cat(sprintf("[pool] Pooled dataset: %d rows x %d cols\n",
              nrow(stacked), ncol(stacked)))
  for (cn in unique(stacked$country)) {
    n_cn <- sum(stacked$country == cn)
    prev <- mean(stacked$Y_binary[stacked$country == cn], na.rm = TRUE)
    cat(sprintf("  %s: n = %d, prevalence = %.1f%%\n", cn, n_cn, prev * 100))
  }

  list(
    data         = stacked,
    Xvars_common = Xvars_common,
    outcome_tag  = outcome_tag
  )
}


#' Run leave-one-country-out cross-validation
#'
#' For each country, trains a SuperLearner on all other countries using only
#' common proxy predictors, then evaluates on the held-out country.
#'
#' @param pooled Output from build_pooled_dataset()
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @return data.frame with one row per held-out country
run_loco_cv <- function(pooled, sl_learners, params) {

  d       <- pooled$data
  Xvars   <- pooled$Xvars_common
  countries <- unique(d$country)

  # Use mlr3 library from sl_learners
  mlr3_lib <- sl_learners$library

  results <- list()

  for (held_out in countries) {
    cat(sprintf("\n[LOCO] Held out: %s\n", held_out))

    train_data <- d[d$country != held_out, ]
    test_data  <- d[d$country == held_out, ]

    cat(sprintf("  Train: %d rows (%s)\n", nrow(train_data),
                paste(setdiff(countries, held_out), collapse = " + ")))
    cat(sprintf("  Test:  %d rows (%s)\n", nrow(test_data), held_out))

    # Fit mlr3 SL on training countries
    fit <- tryCatch({
      t0 <- proc.time()
      res <- mlr3_SL_clustered(
        d          = train_data,
        Xvars      = Xvars,
        outcome    = "Y_binary",
        population = "pooled",
        id         = "pooled_cluster",
        folds      = params$K,
        mlr3_library = mlr3_lib,
        outcome_type = "binomial",
        prescreen  = TRUE
      )
      elapsed <- (proc.time() - t0)["elapsed"]
      cat(sprintf("  SL fit in %.1f min\n", elapsed / 60))
      res
    }, error = function(e) {
      cat(sprintf("  *** LOCO SL fit FAILED for held_out=%s ***\n", held_out))
      cat(sprintf("  Error: %s\n", e$message))
      NULL
    })

    cat(sprintf("  Fit result: %s\n", if (is.null(fit)) "NULL (FAILED)" else "OK"))

    if (is.null(fit)) {
      results[[held_out]] <- data.frame(
        held_out = held_out,
        n_train = nrow(train_data), n_test = nrow(test_data),
        auc = NA_real_, brier = NA_real_,
        calib_intercept = NA_real_, calib_slope = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Predict on held-out country using saved preprocessing recipe
    cat(sprintf("  Model Xvars: %d, pre_recipe_cols: %d\n",
                length(fit$Xvars), length(fit$pre_recipe_cols)))
    preds <- tryCatch({
      p <- predict_on_new_data(fit, test_data, Xvars)
      cat(sprintf("  Prediction succeeded: %d values\n", length(p)))
      p
    }, error = function(e) {
      cat(sprintf("  LOCO prediction FAILED for %s: %s\n", held_out, e$message))
      cat(sprintf("  Traceback: %s\n", paste(deparse(sys.calls()), collapse="\n")))
      NULL
    })

    # If prediction failed, try manual preprocessing fallback
    if (is.null(preds) || length(preds) != nrow(test_data)) {
      cat("  Trying direct prediction with manual preprocessing...\n")
      preds <- tryCatch({
        final_covars <- fit$Xvars
        pre_cols <- fit$pre_recipe_cols

        X0 <- test_data %>% dplyr::select(dplyr::any_of(Xvars)) %>% as.data.frame()
        cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
        cov0 <- cov0[, !sapply(cov0, function(x) all(is.na(x))), drop = FALSE]

        # Impute missing
        cov0 <- as.data.frame(
          ck37r::impute_missing_values(cov0, type = "standard",
                                       add_indicators = TRUE,
                                       prefix = "missing_")$data
        )

        # Align to pre-recipe columns
        for (col in setdiff(pre_cols, colnames(cov0))) cov0[[col]] <- 0
        cov0 <- cov0[, pre_cols, drop = FALSE]

        # Apply recipe
        cov0 <- recipes::bake(fit$auto_recipe, new_data = cov0)
        cov0 <- data.frame(cov0)

        # Final alignment
        for (col in setdiff(final_covars, colnames(cov0))) cov0[[col]] <- 0
        cov0 <- cov0[, final_covars, drop = FALSE]

        # Use the mlr3 wrapper's predict method
        as.numeric(fit$sl_fit$predict(cov0))
      }, error = function(e) {
        cat(sprintf("  Manual prediction also failed for %s: %s\n", held_out, e$message))
        NULL
      })
    }

    if (is.null(preds) || length(preds) != nrow(test_data)) {
      results[[held_out]] <- data.frame(
        held_out = held_out,
        n_train = nrow(train_data), n_test = nrow(test_data),
        auc = NA_real_, brier = NA_real_,
        calib_intercept = NA_real_, calib_slope = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Evaluate
    Y_test <- test_data$Y_binary
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y_test, preds, quiet = TRUE))),
      error = function(e) NA_real_
    )
    brier_val <- mean((Y_test - preds)^2, na.rm = TRUE)

    # Calibration: logit(pred) ~ Y
    calib <- tryCatch({
      preds_clamp <- pmin(pmax(preds, 1e-6), 1 - 1e-6)
      glm_cal <- glm(Y_test ~ qlogis(preds_clamp), family = binomial())
      coef(glm_cal)
    }, error = function(e) c(NA_real_, NA_real_))

    # Admin1 predictions on held-out country
    admin1_preds <- NULL
    if ("Admin1" %in% colnames(test_data)) {
      admin1_preds <- data.frame(
        Admin1 = test_data$Admin1,
        pred   = preds,
        obs    = Y_test,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(!is.na(Admin1)) %>%
        dplyr::group_by(Admin1) %>%
        dplyr::summarise(
          n          = dplyr::n(),
          obs_prev   = mean(obs, na.rm = TRUE),
          pred_prev  = mean(pred, na.rm = TRUE),
          .groups    = "drop"
        )
    }

    results[[held_out]] <- data.frame(
      held_out        = held_out,
      n_train         = nrow(train_data),
      n_test          = nrow(test_data),
      prev_test       = mean(Y_test, na.rm = TRUE),
      auc             = auc_val,
      brier           = brier_val,
      null_brier      = mean(Y_test) * (1 - mean(Y_test)),
      calib_intercept = calib[1],
      calib_slope     = calib[2],
      stringsAsFactors = FALSE
    )

    cat(sprintf("  AUC = %.3f | Brier = %.4f | Prevalence = %.1f%%\n",
                auc_val, brier_val, mean(Y_test) * 100))

    # Store admin1 predictions as attribute
    attr(results[[held_out]], "admin1_preds") <- admin1_preds
  }

  dplyr::bind_rows(results)
}


#' Predict on new data using a fitted DHS_SL_clustered model
#'
#' Applies the saved preprocessing recipe from the training fit to new data
#' and generates predictions. Mirrors the prediction logic in one_bootstrap_rep().
#'
#' @param fit Output from DHS_SL_clustered()
#' @param newdata data.frame of new observations
#' @param Xvars Character vector of predictor column names
#' @return numeric vector of predictions
predict_on_new_data <- function(fit, newdata, Xvars) {

  pre_cols     <- fit$pre_recipe_cols
  final_covars <- fit$Xvars

  cat(sprintf("    [predict] newdata: %d rows, %d Xvars requested\n",
              nrow(newdata), length(Xvars)))

  # Select raw predictor columns available in newdata
  available <- intersect(Xvars, colnames(newdata))
  cat(sprintf("    [predict] %d/%d Xvars found in newdata\n",
              length(available), length(Xvars)))

  X0 <- newdata[, available, drop = FALSE] %>% as.data.frame()
  cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)

  # Drop all-NA columns
  all_na <- sapply(cov0, function(x) all(is.na(x)))
  if (any(all_na)) {
    cat(sprintf("    [predict] Dropping %d all-NA columns\n", sum(all_na)))
    cov0 <- cov0[, !all_na, drop = FALSE]
  }

  cat(sprintf("    [predict] Before impute: %d cols\n", ncol(cov0)))

  # Impute missing values (same method as training)
  cov0 <- as.data.frame(
    ck37r::impute_missing_values(cov0, type = "standard",
                                 add_indicators = TRUE,
                                 prefix = "missing_")$data
  )
  cat(sprintf("    [predict] After impute: %d cols\n", ncol(cov0)))

  # Align to the pre-recipe columns the model expects.
  missing_pre <- setdiff(pre_cols, colnames(cov0))
  cat(sprintf("    [predict] pre_recipe_cols: %d, missing %d (padding with 0)\n",
              length(pre_cols), length(missing_pre)))
  for (col in missing_pre) cov0[[col]] <- 0
  cov0 <- cov0[, pre_cols, drop = FALSE]

  # Apply the saved recipe (step_corr + step_normalize)
  cov0 <- recipes::bake(fit$auto_recipe, new_data = cov0)
  cov0 <- data.frame(cov0)
  cat(sprintf("    [predict] After bake: %d cols\n", ncol(cov0)))

  # Final alignment to the exact post-recipe columns
  missing_final <- setdiff(final_covars, colnames(cov0))
  cat(sprintf("    [predict] final_covars: %d, missing %d (padding with 0)\n",
              length(final_covars), length(missing_final)))
  for (col in missing_final) cov0[[col]] <- 0
  cov0 <- cov0[, final_covars, drop = FALSE]

  # Predict using the model wrapper (works for both sl3 and mlr3)
  result <- as.numeric(fit$sl_fit$predict(cov0))
  cat(sprintf("    [predict] Got %d predictions, range: [%.4f, %.4f]\n",
              length(result), min(result, na.rm=TRUE), max(result, na.rm=TRUE)))
  result
}


# =============================================================================
# GEE-Only LOCO Analysis
#
# Uses only Google Earth Engine (remotely sensed) variables for LOCO CV.
# Tests whether freely available satellite-derived data alone can predict
# micronutrient deficiency prevalence across countries.
# =============================================================================

#' Build a pooled dataset using only GEE variables
#'
#' Same structure as build_pooled_dataset() but restricts predictors to
#' variables matching "gee_" prefix (Google Earth Engine remote sensing).
#' This tests whether freely available satellite-derived data alone can
#' predict micronutrient deficiency across countries.
#'
#' @param all_merged Named list of merged data frames (one per country)
#' @param all_configs Output of get_country_configs()
#' @param outcome_tag Character, e.g. "child_vitA"
#' @return list with data, Xvars_common, outcome_tag
build_pooled_gee_only <- function(all_merged, all_configs, outcome_tag) {

  pooled <- build_pooled_dataset(all_merged, all_configs, outcome_tag)

  if (is.null(pooled) || nrow(pooled$data) == 0) return(pooled)

  # Restrict to GEE variables only
  gee_vars <- grep("^gee_", pooled$Xvars_common, value = TRUE)
  cat(sprintf("[pool_gee] %d GEE variables out of %d total common predictors\n",
              length(gee_vars), length(pooled$Xvars_common)))

  if (length(gee_vars) == 0) {
    warning("No GEE variables found in common predictors")
    return(NULL)
  }

  pooled$Xvars_common <- gee_vars
  pooled
}


#' Run leave-one-country-out CV using only GEE (remote sensing) predictors
#'
#' Identical to run_loco_cv() but filters the pooled dataset to only
#' use variables with "gee_" prefix. This tests whether freely available
#' satellite-derived data can predict micronutrient deficiency across borders.
#'
#' @param pooled_gee Output from build_pooled_gee_only()
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @return data.frame with one row per held-out country
run_loco_gee_only <- function(pooled_gee, sl_learners, params) {

  if (is.null(pooled_gee) || nrow(pooled_gee$data) == 0) {
    warning("No pooled GEE data available for LOCO")
    return(data.frame())
  }

  cat(sprintf("\n[LOCO-GEE] Using %d GEE-only predictors\n",
              length(pooled_gee$Xvars_common)))

  # Delegate to the standard LOCO function — it uses pooled$Xvars_common
  # which we've already restricted to gee_ variables
  run_loco_cv(pooled_gee, sl_learners, params)
}


#' Summarize transportability results across outcomes
#'
#' @param loco_list Named list of run_loco_cv() outputs (one per outcome)
#' @return data.frame with all LOCO results, one row per held-out country x outcome
summarize_transportability <- function(loco_list) {
  out <- dplyr::bind_rows(lapply(names(loco_list), function(tag) {
    res <- loco_list[[tag]]
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$outcome <- tag
    res
  }))

  if (nrow(out) > 0) {
    # Save CSV
    out_dir <- here::here("results", "tables")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(out, file.path(out_dir, "transportability_loco_results.csv"))
    cat(sprintf("[transportability] Saved %d rows to results/tables/transportability_loco_results.csv\n",
                nrow(out)))
  }

  out
}
