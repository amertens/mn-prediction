# =============================================================================
# src/analysis/05_domain_ablation.R
#
# Domain (predictor group) ablation study.
#
# For each outcome × each predictor domain:
#   1. Fit a "full" SL model with ALL predictor domains
#   2. Fit a "reduced" SL model with the focal domain REMOVED
#   3. Compare CV performance: RMSE for continuous, AUC + Brier for binary
#   4. Compute Δ (full − reduced) to measure domain contribution
#
# Larger |Δ| = the domain contributes more predictive information.
# For RMSE: positive Δ means dropping the domain *increases* RMSE (harmful).
# For AUC : positive Δ means dropping the domain *decreases* AUC (harmful).
#
# The full-model fit is reused from sl_results_cont / sl_results_bin where
# already computed; only the reduced fits require additional computation.
#
# Outputs:
#   results/tables/gambia_domain_ablation_summary.csv
#
# Inputs  : cfg, gw_data_list, Xvars_full, domain_vars,
#           sl_results_cont, sl_results_bin,
#           DHS_SL_clustered (defined in 02_fit_sl_models.R)
# =============================================================================

cat("\n[05] Domain ablation...\n")

dir.create(cfg$out_tables, showWarnings = FALSE, recursive = TRUE)


# ---- Helper: extract CV metrics from a fitted DHS_SL_clustered result ------
extract_metrics <- function(res_obj, cutoff, cutoff_dir, model_type) {

  if (is.null(res_obj)) return(NULL)

  r    <- res_obj$res
  Y    <- r$Y
  Yhat <- r$yhat_full

  if (model_type == "continuous") {
    rmse_val <- sqrt(mean((Y - Yhat)^2, na.rm = TRUE))
    mae_val  <- mean(abs(Y - Yhat),    na.rm = TRUE)
    R2_val   <- 1 - mean((Y - Yhat)^2, na.rm = TRUE) / var(Y, na.rm = TRUE)

    # AUC treating continuous as a score for binary deficiency
    Y_bin   <- apply_threshold(Y, cutoff, cutoff_dir)
    score   <- if (cutoff_dir == "less") -Yhat else Yhat
    ok      <- !is.na(Y_bin) & is.finite(score)
    auc_val <- tryCatch({
      roc_obj <- pROC::roc(response = Y_bin[ok], predictor = score[ok],
                           quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA_real_)

    Yhat_bin    <- apply_threshold(Yhat, cutoff, cutoff_dir)
    brier_val   <- mean((Yhat_bin - Y_bin)^2, na.rm = TRUE)

    data.frame(rmse = rmse_val, mae = mae_val, R2 = R2_val,
               auc = auc_val,  brier = brier_val)

  } else {   # binary
    auc_val <- tryCatch({
      roc_obj <- pROC::roc(response = Y, predictor = Yhat,
                           quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA_real_)

    brier_val <- mean((Yhat - Y)^2, na.rm = TRUE)
    data.frame(rmse = NA_real_, mae = NA_real_, R2 = NA_real_,
               auc = auc_val, brier = brier_val)
  }
}


# ---- Helper: fit a reduced model (one domain removed) ----------------------
fit_reduced <- function(d, Xvars_full_local, remove_domain_cols,
                        outcome_col, population, id, K, sl) {

  Xvars_reduced <- setdiff(Xvars_full_local, remove_domain_cols)
  Xvars_reduced <- Xvars_reduced[Xvars_reduced %in% colnames(d)]

  if (length(Xvars_reduced) == 0) return(NULL)

  tryCatch(
    DHS_SL_clustered(
      d          = d,
      Xvars      = Xvars_reduced,
      outcome    = outcome_col,
      population = population,
      id         = id,
      folds      = K,
      CV         = FALSE,
      prescreen  = TRUE,
      sl         = sl
    ),
    error = function(e) {
      cat(sprintf("      [error] %s\n", conditionMessage(e)))
      NULL
    }
  )
}


# ============================================================================
# Main ablation loop
# ============================================================================
ablation_rows <- list()

for (tag in names(cfg$outcomes)) {

  outcome_cfg <- cfg$outcomes[[tag]]
  d_orig      <- gw_data_list[[tag]]
  res_cont    <- sl_results_cont[[tag]]
  res_bin     <- sl_results_bin[[tag]]

  # Decide primary model type
  use_binary  <- !is.null(res_bin)
  model_type  <- if (use_binary) "binary" else "continuous"
  outcome_col <- if (use_binary) outcome_cfg$binary else outcome_cfg$continuous
  sl_obj      <- if (use_binary) slmod2_bin else slmod
  res_full    <- if (use_binary) res_bin else res_cont

  if (is.null(res_full)) {
    cat(sprintf("  [skip] no model for %s\n", tag))
    next
  }

  cat(sprintf("\n  Ablation for: %s  (model: %s)\n", tag, model_type))

  d_fit    <- d_orig[!is.na(d_orig[[outcome_col]]), ]
  Xvars_b  <- Xvars_full[Xvars_full %in% colnames(d_fit)]

  # Full-model metrics (reuse existing fit)
  full_metrics <- extract_metrics(res_full, outcome_cfg$cutoff,
                                  outcome_cfg$cutoff_dir, model_type)
  if (is.null(full_metrics)) next

  cat(sprintf("    Full model: RMSE=%.4f  AUC=%.4f  Brier=%.4f\n",
              ifelse(is.na(full_metrics$rmse), NA, full_metrics$rmse),
              full_metrics$auc, full_metrics$brier))

  # Loop over domains
  for (dom_name in names(cfg$domains)) {

    dom_cols <- domain_vars[[dom_name]]  # columns in this domain
    # Remove only columns that actually appear in Xvars_b
    remove_cols <- intersect(dom_cols, Xvars_b)

    if (length(remove_cols) == 0) {
      cat(sprintf("    Domain %-10s: [no columns, skip]\n", dom_name))
      next
    }

    cat(sprintf("    Domain %-10s (removing %d cols): ",
                dom_name, length(remove_cols)))

    res_red <- fit_reduced(
      d                  = d_fit,
      Xvars_full_local   = Xvars_b,
      remove_domain_cols = remove_cols,
      outcome_col        = outcome_col,
      population         = outcome_cfg$population,
      id                 = cfg$cluster_id,
      K                  = cfg$K,
      sl                 = sl_obj
    )

    red_metrics <- extract_metrics(res_red, outcome_cfg$cutoff,
                                   outcome_cfg$cutoff_dir, model_type)

    if (is.null(red_metrics)) {
      cat("failed\n")
      row_i <- data.frame(
        outcome      = tag,
        model_type   = model_type,
        domain       = dom_name,
        n_dom_cols   = length(remove_cols),
        rmse_full    = full_metrics$rmse,
        rmse_reduced = NA_real_,
        delta_rmse   = NA_real_,
        auc_full     = full_metrics$auc,
        auc_reduced  = NA_real_,
        delta_auc    = NA_real_,
        brier_full   = full_metrics$brier,
        brier_reduced= NA_real_,
        delta_brier  = NA_real_
      )
    } else {
      cat(sprintf("AUC=%.4f  Brier=%.4f\n",
                  red_metrics$auc, red_metrics$brier))

      row_i <- data.frame(
        outcome       = tag,
        model_type    = model_type,
        domain        = dom_name,
        n_dom_cols    = length(remove_cols),
        rmse_full     = full_metrics$rmse,
        rmse_reduced  = red_metrics$rmse,
        delta_rmse    = red_metrics$rmse - full_metrics$rmse,   # + = domain helps
        auc_full      = full_metrics$auc,
        auc_reduced   = red_metrics$auc,
        delta_auc     = full_metrics$auc - red_metrics$auc,     # + = domain helps
        brier_full    = full_metrics$brier,
        brier_reduced = red_metrics$brier,
        delta_brier   = red_metrics$brier - full_metrics$brier  # + = domain helps
      )
    }

    ablation_rows[[paste(tag, dom_name, sep = "_")]] <- row_i
  }
}

# ---- Compile and save ablation table ----------------------------------------
ablation_table <- do.call(rbind, ablation_rows)

if (!is.null(ablation_table) && nrow(ablation_table) > 0) {
  # Sort by outcome then |delta_auc| descending
  ablation_table <- ablation_table %>%
    dplyr::arrange(outcome,
                   dplyr::desc(abs(ifelse(is.na(delta_auc), 0, delta_auc))))

  abl_csv <- file.path(cfg$out_tables, "gambia_domain_ablation_summary.csv")
  readr::write_csv(ablation_table, abl_csv)
  cat(sprintf("\n  Ablation table saved: %s\n", basename(abl_csv)))

  # Print top contributions per outcome
  cat("\n  Top domain contributions (by |delta_auc|):\n")
  top <- ablation_table %>%
    dplyr::filter(!is.na(delta_auc)) %>%
    dplyr::group_by(outcome) %>%
    dplyr::slice_max(order_by = abs(delta_auc), n = 3) %>%
    dplyr::ungroup() %>%
    dplyr::select(outcome, domain, delta_auc, delta_brier)
  print(top, digits = 4, row.names = FALSE)
} else {
  cat("  [warn] No ablation results produced.\n")
}

cat("[05] Done.\n\n")
