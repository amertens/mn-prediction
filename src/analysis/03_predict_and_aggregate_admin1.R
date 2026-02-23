# =============================================================================
# src/analysis/03_predict_and_aggregate_admin1.R
#
# Extracts individual-level predictions from fitted SL models, converts to
# deficiency indicators, aggregates to Admin1 prevalence, and saves:
#   • results/tables/gambia_admin1_prevalence_{tag}.csv
#   • results/figures/admin1_scatter_{tag}.png
#
# Also reports within-sample CV performance:
#   • Continuous : RMSE, MAE, R²   (from cross-validated predictions)
#   • Binary     : AUC, Brier, calibration slope/intercept
#
# Strategy:
#   - If a binary SL fit is available (sl_results_bin), use its predicted
#     probabilities directly as prevalence estimates (preferred).
#   - Otherwise, threshold the continuous SL's cross-validated predictions at
#     cfg$outcomes[[tag]]$cutoff to produce 0/1 deficiency calls, then average
#     by Admin1.  This is a deterministic transformation; no probability
#     calibration or threshold tuning is applied.
#
# Inputs  : cfg, gw_data_list, sl_results_cont, sl_results_bin
# Outputs : CSVs in results/tables/, PNGs in results/figures/
#           admin1_prev_list  - named list of admin1 prevalence tables
#           cv_perf_table     - data.frame of CV performance metrics
# =============================================================================

cat("\n[03] Predicting and aggregating to Admin1...\n")

library(ggplot2)

dir.create(cfg$out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_figures, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: apply deficiency threshold ------------------------------------
# Returns 1 (deficient) or 0 (not deficient) for a vector of continuous values.
apply_threshold <- function(x, cutoff, direction = "less") {
  if (direction == "less")    as.integer(x < cutoff)
  else                        as.integer(x > cutoff)
}


# ---- Helper: calibration slope/intercept -----------------------------------
# Regresses observed binary outcome on logit(predicted probability).
# Returns slope and intercept; perfect calibration → slope=1, intercept=0.
calibration_stats <- function(obs, pred_prob) {
  pred_prob <- pmin(pmax(pred_prob, 1e-6), 1 - 1e-6)  # clamp from (0,1)
  logit_p   <- log(pred_prob / (1 - pred_prob))
  fit       <- suppressWarnings(glm(obs ~ logit_p, family = binomial()))
  coef_vals <- coef(fit)
  data.frame(
    calib_intercept = unname(coef_vals[1]),
    calib_slope     = unname(coef_vals[2])
  )
}


# ---- Helper: CV performance for one outcome --------------------------------
cv_performance <- function(tag, res_cont, res_bin, outcome_cfg, d_orig) {

  admin1_col <- cfg$admin1
  perf_rows  <- list()

  # -- Continuous performance --
  if (!is.null(res_cont)) {
    r    <- res_cont$res
    Y    <- r$Y
    Yhat <- r$yhat_full

    # Continuous metrics
    mse_model  <- mean((Y - Yhat)^2, na.rm = TRUE)
    mse_mean   <- mean((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE)
    rmse       <- sqrt(mse_model)
    mae        <- mean(abs(Y - Yhat), na.rm = TRUE)
    R2         <- 1 - mse_model / var(Y, na.rm = TRUE)
    R2_vs_mean <- 1 - mse_model / mse_mean
    corr       <- cor(Y, Yhat, use = "complete.obs")

    # Binary metrics from continuous predictions (threshold at cutoff)
    # We use the CONTINUOUS predicted value as the "score" for AUC;
    # direction is reversed (lower score → more deficient → positive class).
    cutoff  <- outcome_cfg$cutoff
    dir     <- outcome_cfg$cutoff_dir
    Y_bin   <- apply_threshold(Y, cutoff, dir)
    # For AUC: positive class = deficient (Y_bin==1).
    # For Vitamin A and iron, a LOWER predicted value → more deficient,
    # so we negate Yhat as the score (higher score = more deficient).
    score  <- if (dir == "less") -Yhat else Yhat
    ok     <- !is.na(Y_bin) & !is.na(score)
    auc_cont <- tryCatch({
      roc_obj <- pROC::roc(response = Y_bin[ok], predictor = score[ok],
                           quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA_real_)

    Yhat_bin <- apply_threshold(Yhat, cutoff, dir)
    brier_cont <- mean((Yhat_bin - Y_bin)^2, na.rm = TRUE)

    perf_rows[["continuous"]] <- data.frame(
      outcome    = tag,
      model_type = "continuous",
      rmse       = rmse,
      mae        = mae,
      R2         = R2,
      R2_vs_mean = R2_vs_mean,
      correlation= corr,
      auc        = auc_cont,
      brier      = brier_cont,
      calib_intercept = NA_real_,
      calib_slope     = NA_real_
    )
  }

  # -- Binary performance --
  if (!is.null(res_bin)) {
    r      <- res_bin$res
    Y_bin  <- r$Y           # already 0/1
    Yhat_p <- r$yhat_full   # predicted probability

    auc_bin <- tryCatch({
      roc_obj <- pROC::roc(response = Y_bin, predictor = Yhat_p,
                           quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA_real_)

    brier_bin  <- mean((Yhat_p - Y_bin)^2, na.rm = TRUE)
    calib      <- tryCatch(
      calibration_stats(Y_bin, Yhat_p),
      error = function(e) data.frame(calib_intercept = NA_real_, calib_slope = NA_real_)
    )

    perf_rows[["binary"]] <- data.frame(
      outcome    = tag,
      model_type = "binary",
      rmse       = NA_real_,
      mae        = NA_real_,
      R2         = NA_real_,
      R2_vs_mean = NA_real_,
      correlation= NA_real_,
      auc        = auc_bin,
      brier      = brier_bin,
      calib_intercept = calib$calib_intercept,
      calib_slope     = calib$calib_slope
    )
  }

  do.call(rbind, perf_rows)
}


# ---- Helper: Admin1 prevalence aggregation ---------------------------------
admin1_prevalence <- function(tag, res_cont, res_bin, d_orig, outcome_cfg) {

  admin1_col <- cfg$admin1
  cutoff     <- outcome_cfg$cutoff
  dir        <- outcome_cfg$cutoff_dir

  # Merge predictions with Admin1 column from the original dataset
  # dataid was created in step 01 and exists in both res$res and d_orig
  base_df <- d_orig %>%
    dplyr::select(dataid, dplyr::all_of(admin1_col)) %>%
    dplyr::distinct()

  # --- Observed deficiency (from continuous or pre-computed binary) ---------
  # Prefer pre-computed binary column; fall back to threshold on continuous.
  bin_col <- outcome_cfg$binary
  if (!is.null(bin_col) && bin_col %in% colnames(d_orig)) {
    obs_src <- d_orig %>%
      dplyr::select(dataid, dplyr::all_of(admin1_col),
                    Y_obs_bin = dplyr::all_of(bin_col)) %>%
      dplyr::mutate(Y_obs_bin = as.integer(.data$Y_obs_bin))
  } else {
    # Threshold continuous observed values
    cont_col <- outcome_cfg$continuous
    obs_src <- d_orig %>%
      dplyr::select(dataid, dplyr::all_of(admin1_col),
                    Y_cont = dplyr::all_of(cont_col)) %>%
      dplyr::mutate(Y_obs_bin = apply_threshold(.data$Y_cont, cutoff, dir)) %>%
      dplyr::select(-Y_cont)
  }

  # --- Predicted deficiency --------------------------------------------------
  # Prefer binary SL predicted probabilities; fall back to thresholded
  # continuous predictions.  For binary SL, predicted prevalence =
  # mean(predicted probability).  For thresholded continuous,
  # predicted prevalence = mean(predicted_biomarker < cutoff).

  if (!is.null(res_bin)) {
    pred_src <- res_bin$res %>%
      dplyr::select(dataid, yhat_pred = yhat_full)
    pred_type <- "binary_prob"
  } else if (!is.null(res_cont)) {
    pred_src <- res_cont$res %>%
      dplyr::mutate(yhat_pred = apply_threshold(yhat_full, cutoff, dir)) %>%
      dplyr::select(dataid, yhat_pred)
    pred_type <- "continuous_threshold"
  } else {
    cat(sprintf("  [skip] no fitted model for %s\n", tag))
    return(NULL)
  }

  # Join observed + predicted + Admin1
  prev_df <- obs_src %>%
    dplyr::left_join(pred_src, by = "dataid") %>%
    dplyr::filter(!is.na(.data[[admin1_col]]),
                  !is.na(Y_obs_bin),
                  !is.na(yhat_pred))

  # Aggregate to Admin1
  admin1_tbl <- prev_df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(
      n                  = dplyr::n(),
      obs_prevalence     = mean(Y_obs_bin,  na.rm = TRUE),
      pred_prevalence    = mean(yhat_pred,  na.rm = TRUE),
      .groups            = "drop"
    ) %>%
    dplyr::mutate(
      outcome    = tag,
      label      = outcome_cfg$label,
      pred_type  = pred_type
    ) %>%
    dplyr::rename(Admin1 = dplyr::all_of(admin1_col))

  admin1_tbl
}


# ---- Helper: scatter plot ---------------------------------------------------
make_scatter <- function(tbl, outcome_cfg, outfile) {

  if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))

  r2_val  <- cor(tbl$obs_prevalence, tbl$pred_prevalence,
                 use = "complete.obs")^2
  n_admin <- nrow(tbl)

  p <- ggplot(tbl, aes(x = obs_prevalence, y = pred_prevalence,
                       label = Admin1)) +
    geom_point(size = 2.5, colour = "#2166AC", alpha = 0.8) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 15,
                             colour = "grey40") +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey30") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    labs(
      title    = paste0("Admin1 prevalence: ", outcome_cfg$label),
      subtitle = sprintf("n = %d Admin1 units  |  R\u00b2 = %.3f", n_admin, r2_val),
      x        = "Observed prevalence",
      y        = "Predicted prevalence",
      caption  = paste0("Prediction type: ", tbl$pred_type[1])
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # ggrepel is optional; fall back to plain geom_text if absent
  tryCatch(
    ggsave(outfile, plot = p, width = 7, height = 6, dpi = 180),
    error = function(e) {
      # Try without ggrepel labels
      p2 <- p + ggplot2::geom_text(size = 2.8, colour = "grey40",
                                   vjust = -0.5, hjust = 0.5)
      ggsave(outfile, plot = p2, width = 7, height = 6, dpi = 180)
    }
  )
  invisible(p)
}


# ============================================================================
# Main loop over outcomes
# ============================================================================
admin1_prev_list <- list()
cv_perf_rows     <- list()

for (tag in names(cfg$outcomes)) {

  cat(sprintf("  Processing: %s\n", tag))
  outcome_cfg <- cfg$outcomes[[tag]]
  d_orig      <- gw_data_list[[tag]]
  res_cont    <- sl_results_cont[[tag]]
  res_bin     <- sl_results_bin[[tag]]

  # CV performance
  perf <- tryCatch(
    cv_performance(tag, res_cont, res_bin, outcome_cfg, d_orig),
    error = function(e) {
      cat(sprintf("    [warn] CV perf error: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(perf)) cv_perf_rows[[tag]] <- perf

  # Admin1 prevalence
  a1_tbl <- tryCatch(
    admin1_prevalence(tag, res_cont, res_bin, d_orig, outcome_cfg),
    error = function(e) {
      cat(sprintf("    [warn] Admin1 error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(a1_tbl)) {
    admin1_prev_list[[tag]] <- a1_tbl

    # Save CSV
    csv_file <- file.path(
      cfg$out_tables,
      sprintf("gambia_admin1_prevalence_%s.csv", outcome_cfg$table_tag)
    )
    readr::write_csv(a1_tbl, csv_file)
    cat(sprintf("    Saved: %s\n", basename(csv_file)))

    # Save scatter plot
    png_file <- file.path(
      cfg$out_figures,
      sprintf("admin1_scatter_%s.png", outcome_cfg$scatter_tag)
    )
    tryCatch(
      make_scatter(a1_tbl, outcome_cfg, png_file),
      error = function(e) cat(sprintf("    [warn] plot error: %s\n",
                                      conditionMessage(e)))
    )
    cat(sprintf("    Saved: %s\n", basename(png_file)))
  }
}

# ---- Save CV performance table ---------------------------------------------
cv_perf_table <- do.call(rbind, cv_perf_rows)
if (!is.null(cv_perf_table) && nrow(cv_perf_table) > 0) {
  cv_csv <- file.path(cfg$out_tables, "gambia_cv_performance.csv")
  readr::write_csv(cv_perf_table, cv_csv)
  cat(sprintf("\n  CV performance saved: %s\n", basename(cv_csv)))

  print(cv_perf_table[, c("outcome", "model_type", "rmse", "auc", "brier")],
        row.names = FALSE)
}

cat("[03] Done.\n\n")
