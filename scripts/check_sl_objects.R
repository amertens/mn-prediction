# =============================================================================
# scripts/check_sl_objects.R
#
# Diagnostic check of SL fit objects:
# 1. Are models non-empty (not just 131 bytes)?
# 2. Does predict() work on original train_data?
# 3. Does predict() produce non-degenerate output?
# 4. Does predict() change when features are permuted?
# 5. Do learner_models (serialization-safe fallback) exist?
# 6. What are the SL weights?
#
# Run after pipeline completes:
#   source("scripts/check_sl_objects.R")
# =============================================================================

library(here)
library(targets)
library(pROC)

cat("\n=== SuperLearner Object Diagnostics ===\n\n")

countries <- c("gambia", "ghana", "sierraleone", "malawi")
outcomes  <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

results <- list()

for (ctry in countries) {
  for (oc in outcomes) {
    suffix <- paste0(ctry, "_", oc)
    tname <- paste0("sl_fit_", suffix)

    sl <- tryCatch(tar_read_raw(tname, store = "_targets"), error = function(e) NULL)
    if (is.null(sl)) {
      cat(sprintf("  [SKIP] %s: target not found\n", suffix))
      next
    }

    for (mtype in c("bin_fit", "cont_fit")) {
      fit_obj <- sl[[mtype]]
      label <- paste0(suffix, " (", mtype, ")")

      if (is.null(fit_obj) || is.null(fit_obj$res)) {
        cat(sprintf("  [EMPTY] %s: no model\n", label))
        next
      }

      Y <- as.numeric(fit_obj$res$Y)
      yhat_cv <- as.numeric(fit_obj$res$yhat_full)
      n <- length(Y)
      n_events <- sum(Y == 1, na.rm = TRUE)
      covars <- fit_obj$Xvars
      td <- fit_obj$train_data

      # Check 1: Basic structure
      has_td <- !is.null(td) && nrow(td) > 0
      has_lm <- !is.null(fit_obj$sl_fit$learner_models)
      has_wt <- !is.null(fit_obj$sl_fit$sl_weights)

      # Check 2: Predict on original data
      pred_ok <- FALSE
      pred_degen <- TRUE
      if (has_td) {
        pred <- tryCatch(
          as.numeric(fit_obj$sl_fit$predict(data.frame(td))),
          error = function(e) NULL
        )
        if (!is.null(pred) && length(pred) == n) {
          pred_ok <- TRUE
          pred_degen <- length(unique(round(pred, 6))) <= 1
        }
      }

      # Check 3: Permutation changes predictions
      perm_changes <- NA
      if (pred_ok && !pred_degen && length(covars) > 0) {
        td_perm <- data.frame(td)
        first_var <- covars[1]
        td_perm[[first_var]] <- td_perm[[first_var]][sample(nrow(td_perm))]
        pred_perm <- tryCatch(
          as.numeric(fit_obj$sl_fit$predict(td_perm)),
          error = function(e) NULL
        )
        if (!is.null(pred_perm)) {
          perm_changes <- !isTRUE(all.equal(pred, pred_perm))
        }
      }

      # Check 4: SL weights
      wt_str <- "N/A"
      if (has_wt) {
        wts <- fit_obj$sl_fit$sl_weights
        nonzero <- names(wts[wts > 0.01])
        if (is.null(nonzero)) nonzero <- which(wts > 0.01)
        wt_str <- paste(sprintf("%s=%.2f", names(wts)[wts > 0.01], wts[wts > 0.01]),
                        collapse = ", ")
      }

      # Check 5: AUC (binary only)
      auc_val <- NA
      if (mtype == "bin_fit" && length(unique(Y)) == 2) {
        auc_val <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(Y, yhat_cv, quiet = TRUE))),
          error = function(e) NA
        )
      }

      status <- if (!has_td) "NO_TRAINDATA"
                else if (!pred_ok) "PREDICT_FAIL"
                else if (pred_degen) "DEGENERATE"
                else if (isFALSE(perm_changes)) "PERM_UNCHANGED"
                else "OK"

      cat(sprintf("  [%s] %s: n=%d, p=%d, events=%d, AUC=%.3f, predict=%s, perm_changes=%s, fallback=%s, weights=[%s]\n",
                  status, label, n, length(covars), n_events,
                  ifelse(is.na(auc_val), NA, auc_val),
                  ifelse(pred_ok, ifelse(pred_degen, "DEGEN", "OK"), "FAIL"),
                  ifelse(is.na(perm_changes), "N/A", ifelse(perm_changes, "YES", "NO")),
                  ifelse(has_lm, "YES", "NO"),
                  wt_str))

      results[[label]] <- data.frame(
        suffix = suffix, model = mtype, n = n, p = length(covars),
        n_events = n_events, auc = auc_val,
        has_train_data = has_td, predict_ok = pred_ok,
        predict_degenerate = pred_degen, perm_changes = perm_changes,
        has_learner_models = has_lm, has_weights = has_wt,
        status = status, stringsAsFactors = FALSE
      )
    }
    rm(sl); gc(verbose = FALSE)
  }
}

cat("\n=== Summary ===\n")
if (length(results) > 0) {
  summary_df <- dplyr::bind_rows(results)
  cat(sprintf("Total models: %d\n", nrow(summary_df)))
  cat(sprintf("OK: %d\n", sum(summary_df$status == "OK")))
  cat(sprintf("DEGENERATE: %d\n", sum(summary_df$status == "DEGENERATE")))
  cat(sprintf("PREDICT_FAIL: %d\n", sum(summary_df$status == "PREDICT_FAIL")))
  cat(sprintf("EMPTY: %d\n", sum(summary_df$status == "NO_TRAINDATA")))
  cat(sprintf("Has fallback models: %d\n", sum(summary_df$has_learner_models)))

  # Show problem cases
  problems <- summary_df[summary_df$status != "OK" & summary_df$status != "NO_TRAINDATA", ]
  if (nrow(problems) > 0) {
    cat("\nProblem cases:\n")
    print(problems[, c("suffix", "model", "n", "auc", "status")])
  }
}

cat("\n=== Done ===\n")
