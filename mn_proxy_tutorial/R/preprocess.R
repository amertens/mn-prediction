# =============================================================================
# R/preprocess.R
#
# Preprocessing functions: missingness introduction, imputation,
# standardisation, domain grouping, and leakage prevention.
# =============================================================================

# ── Domain definitions ────────────────────────────────────────────────────────

#' Return a named list grouping covariate columns by thematic domain.
#'
#' @param df Data frame (column names are inspected).
#' @return Named list: domain_name → character vector of column names.
make_domain_groups <- function(df) {
  all_cols <- colnames(df)
  domains  <- list(
    wash    = grep("^wash_",          all_cols, value = TRUE),
    ses     = grep("^ses_",           all_cols, value = TRUE),
    climate = grep("^clim_",          all_cols, value = TRUE),
    disease = grep("^dis_",           all_cols, value = TRUE),
    indiv   = grep("^indiv_",         all_cols, value = TRUE),
    spatial = grep("^[xy]_coord$",    all_cols, value = TRUE)
  )
  domains[sapply(domains, length) > 0]   # drop empty groups
}

#' Return a flat character vector of all predictor columns.
#' Excludes outcome columns, identifiers, and weights.
#'
#' @param df         Data frame.
#' @param exclude    Column names to exclude (outcomes, IDs, weights).
#' @return Character vector of predictor column names.
get_predictor_cols <- function(df,
                                exclude = c("Y", "D", "svy_weight",
                                            "indiv_id", "cluster_id",
                                            "admin1_id", "admin1_name")) {
  domains <- make_domain_groups(df)
  all_pred <- unlist(domains, use.names = FALSE)
  setdiff(all_pred, exclude)
}


# ── Missingness introduction ──────────────────────────────────────────────────

#' Introduce MCAR and MAR missingness into selected columns.
#'
#' MCAR: each value is independently deleted with probability `mcar_rate`.
#' MAR : missingness in `mar_vars` depends on the value of `ses_wealth`
#'       (lower wealth → higher missingness probability).
#'
#' @param df        Data frame (not modified in place; a copy is returned).
#' @param mcar_vars Character vector of column names for MCAR deletion.
#' @param mar_vars  Character vector of column names for MAR deletion.
#' @param mcar_rate MCAR deletion probability.
#' @param mar_rate  Maximum MAR deletion probability.
#' @param seed      RNG seed.
#' @return Modified data frame with NAs introduced.
introduce_missingness <- function(df,
                                   mcar_vars = c("ses_educ", "clim_temp"),
                                   mar_vars  = c("dis_malaria", "wash_sanitation"),
                                   mcar_rate = 0.05,
                                   mar_rate  = 0.12,
                                   seed      = 7) {
  set.seed(seed)
  df_miss <- df

  # MCAR: delete at fixed rate regardless of covariate values
  for (v in mcar_vars) {
    if (v %in% colnames(df_miss)) {
      mask <- runif(nrow(df_miss)) < mcar_rate
      df_miss[[v]][mask] <- NA
    }
  }

  # MAR: missingness probability increases for lower-SES individuals
  if ("ses_wealth" %in% colnames(df_miss)) {
    # Normalise ses_wealth to [0,1] for probability calculation
    w_norm <- (df_miss$ses_wealth - min(df_miss$ses_wealth, na.rm = TRUE)) /
              diff(range(df_miss$ses_wealth, na.rm = TRUE))
    p_miss <- mar_rate * (1 - w_norm)   # poorer → more likely missing
    for (v in mar_vars) {
      if (v %in% colnames(df_miss)) {
        mask <- runif(nrow(df_miss)) < p_miss
        df_miss[[v]][mask] <- NA
      }
    }
  }

  df_miss
}


# ── Imputation ────────────────────────────────────────────────────────────────

#' Impute missing values in predictor columns.
#'
#' Supports two strategies:
#'   "median"    — replace NA with column median (numeric) or mode (factor).
#'   "indicator" — as above, but also add binary missing-indicator columns.
#'
#' Imputation statistics are computed on `df_train` and applied to `df_apply`
#' to prevent leakage when used within cross-validation.
#'
#' @param df_train  Data frame used to compute imputation statistics.
#' @param df_apply  Data frame to impute (defaults to df_train).
#' @param pred_cols Character vector of predictor columns to impute.
#' @param method    "median" or "indicator".
#' @return A list: $data (imputed df_apply), $stats (named list of impute values).
impute_data <- function(df_train,
                         df_apply  = df_train,
                         pred_cols = NULL,
                         method    = c("median", "indicator")) {
  method <- match.arg(method)
  if (is.null(pred_cols)) pred_cols <- get_predictor_cols(df_train)

  # Compute imputation statistics on training data
  stats <- lapply(pred_cols, function(v) {
    x <- df_train[[v]]
    if (is.numeric(x)) {
      median(x, na.rm = TRUE)
    } else {
      # Mode
      tbl <- table(x, useNA = "no")
      names(tbl)[which.max(tbl)]
    }
  })
  names(stats) <- pred_cols

  # Apply to df_apply
  df_out <- df_apply
  new_cols <- character(0)

  for (v in pred_cols) {
    if (!(v %in% colnames(df_out))) next
    miss <- is.na(df_out[[v]])
    if (any(miss)) {
      if (method == "indicator") {
        ind_name <- paste0("miss_", v)
        df_out[[ind_name]] <- as.integer(miss)
        new_cols <- c(new_cols, ind_name)
      }
      df_out[[v]][miss] <- stats[[v]]
    }
  }

  list(data = df_out, stats = stats, indicator_cols = new_cols)
}


# ── Standardisation ───────────────────────────────────────────────────────────

#' Z-score standardise numeric predictor columns.
#'
#' Mean and SD are computed from `df_train` and applied to `df_apply`
#' to prevent leakage.
#'
#' @param df_train  Training data frame.
#' @param df_apply  Data frame to standardise (defaults to df_train).
#' @param pred_cols Columns to standardise (defaults to all numeric predictors).
#' @return List: $data (standardised df_apply), $params (list of mean, sd per col).
standardise_data <- function(df_train,
                              df_apply  = df_train,
                              pred_cols = NULL) {
  if (is.null(pred_cols)) pred_cols <- get_predictor_cols(df_train)

  # Only standardise numeric columns
  num_cols <- pred_cols[sapply(df_train[, pred_cols, drop = FALSE], is.numeric)]

  params <- lapply(num_cols, function(v) {
    list(mean = mean(df_train[[v]], na.rm = TRUE),
         sd   =   sd(df_train[[v]], na.rm = TRUE))
  })
  names(params) <- num_cols

  df_out <- df_apply
  for (v in num_cols) {
    if (!(v %in% colnames(df_out))) next
    mu <- params[[v]]$mean
    sg <- params[[v]]$sd
    if (!is.na(sg) && sg > 0)
      df_out[[v]] <- (df_out[[v]] - mu) / sg
  }

  list(data = df_out, params = params)
}


# ── Master preprocessing pipeline ─────────────────────────────────────────────

#' Full preprocessing pipeline: impute → standardise → return.
#'
#' Designed to be called separately on training and held-out folds so that
#' imputation and standardisation statistics are always computed from training
#' data only (leakage prevention).
#'
#' @param df_train  Training data frame.
#' @param df_apply  Data frame to preprocess (default df_train; pass fold
#'                  validation data to apply training statistics).
#' @param pred_cols Predictor columns (NULL → auto-detected via domain groups).
#' @param method    Imputation method: "median" or "indicator".
#' @return List: $data (preprocessed df_apply), $impute_stats, $scale_params.
preprocess_survey <- function(df_train,
                               df_apply  = df_train,
                               pred_cols = NULL,
                               method    = "indicator") {
  if (is.null(pred_cols)) pred_cols <- get_predictor_cols(df_train)

  # Step 1: impute
  imp   <- impute_data(df_train, df_apply, pred_cols, method)

  # Step 2: standardise (using updated pred_cols that include indicator cols)
  all_pred <- c(pred_cols, imp$indicator_cols)
  scl      <- standardise_data(df_train  = imp$data,    # for stats
                                df_apply  = imp$data,
                                pred_cols = all_pred)

  list(
    data          = scl$data,
    impute_stats  = imp$stats,
    scale_params  = scl$params,
    pred_cols     = all_pred
  )
}


#' Apply previously computed preprocessing statistics to new data.
#'
#' Use this to preprocess validation/test data using training-derived stats
#' (the proper within-CV approach).
#'
#' @param df_new       New data frame (validation or test fold).
#' @param impute_stats Named list of imputation values (from preprocess_survey).
#' @param scale_params Named list of (mean, sd) per column.
#' @param method       Same method as used during training.
#' @return Preprocessed data frame.
apply_preprocessing <- function(df_new, impute_stats, scale_params,
                                 method = "indicator") {
  df_out <- df_new

  # Impute
  for (v in names(impute_stats)) {
    if (!(v %in% colnames(df_out))) next
    miss <- is.na(df_out[[v]])
    if (any(miss)) {
      if (method == "indicator") {
        ind_name <- paste0("miss_", v)
        df_out[[ind_name]] <- as.integer(miss)
      }
      df_out[[v]][miss] <- impute_stats[[v]]
    }
  }

  # Standardise
  for (v in names(scale_params)) {
    if (!(v %in% colnames(df_out))) next
    mu <- scale_params[[v]]$mean
    sg <- scale_params[[v]]$sd
    if (!is.na(sg) && sg > 0)
      df_out[[v]] <- (df_out[[v]] - mu) / sg
  }

  df_out
}
