# =============================================================================
# R/sensitivity/gp_sensitivity.R
#
# Does dropping the Gaussian process from the full-mode stack cost anything?
#
# gaussianprocess was removed from the default full stack in 2026-08 because
# kernlab::gausspr is O(n^3) in training rows: measured at p = 140 on this
# machine it takes 0.2 s at n = 400, 1.0 s at n = 800 and 9.4 s at n = 1600.
# The individual-level LOCO targets pool n = 10,011 and 13,107, where that curve
# implies 8-18 hours per target -- three of them burned 4.9 CPU-hours each
# without finishing.
#
# Dropping a learner because it is slow is a decision about compute, not about
# accuracy, so it needs checking rather than assuming. This refits ONE small
# country x outcome with the GP included and compares held-out performance
# against the same fit without it. Small n is the point: it keeps the GP
# tractable while still telling us whether it contributes signal.
# =============================================================================

#' Held-out performance summary from a fitted SL object.
#'
#' Uses the cross-validated predictions (`yhat_full`), NOT `yhat_insample`.
#' The in-sample column is optimistic by construction and comparing on it would
#' reward the more flexible stack automatically.
#'
#' @param fit One element of fit_mlr3_models() output (`bin_fit` or `cont_fit`)
#' @return one-row data.frame, or NULL if the fit is unusable
gp_perf_summary <- function(fit) {
  if (is.null(fit) || is.null(fit$res)) return(NULL)
  d <- fit$res
  if (!all(c("Y", "yhat_full") %in% names(d))) return(NULL)
  ok <- is.finite(d$Y) & is.finite(d$yhat_full)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 10) return(NULL)

  binary <- length(unique(d$Y)) == 2
  auc <- NA_real_
  if (binary) {
    auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(d$Y, d$yhat_full, quiet = TRUE))),
                    error = function(e) NA_real_)
  }
  data.frame(
    n        = nrow(d),
    binary   = binary,
    auc      = auc,
    rmse     = sqrt(mean((d$Y - d$yhat_full)^2)),
    cor      = suppressWarnings(stats::cor(d$Y, d$yhat_full)),
    mean_yhat = mean(d$yhat_full),
    stringsAsFactors = FALSE
  )
}

#' Compare a with-GP fit against the production (no-GP) fit.
#'
#' @param fit_no_gp  sl_fit_* target built with the default stack
#' @param fit_with_gp same country x outcome refit with with_gp = TRUE
#' @param label      country x outcome being compared
#' @return data.frame, one row per (arm x model type), plus a delta attribute
compare_gp_sensitivity <- function(fit_no_gp, fit_with_gp, label = "") {
  rows <- list()
  for (arm in c("no_gp", "with_gp")) {
    f <- if (arm == "no_gp") fit_no_gp else fit_with_gp
    for (mt in c("bin_fit", "cont_fit")) {
      s <- gp_perf_summary(f[[mt]])
      if (is.null(s)) next
      s$arm <- arm; s$model <- sub("_fit$", "", mt); s$label <- label
      rows[[length(rows) + 1L]] <- s
    }
  }
  if (!length(rows)) {
    warning("gp_sensitivity: no comparable fits for ", label)
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  out <- out[, c("label", "model", "arm", "n", "binary", "auc", "rmse",
                 "cor", "mean_yhat")]

  # Report the difference explicitly so the answer does not depend on the
  # reader eyeballing two rows.
  for (mt in unique(out$model)) {
    a <- out[out$model == mt & out$arm == "no_gp", ]
    b <- out[out$model == mt & out$arm == "with_gp", ]
    if (nrow(a) == 1 && nrow(b) == 1) {
      metric <- if (isTRUE(a$binary)) "AUC" else "correlation"
      va <- if (isTRUE(a$binary)) a$auc else a$cor
      vb <- if (isTRUE(b$binary)) b$auc else b$cor
      cat(sprintf("[gp_sensitivity] %s / %s: %s %.4f (no GP) vs %.4f (with GP), delta = %+.4f\n",
                  label, mt, metric, va, vb, vb - va))
    }
  }

  dir.create(here::here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(out, here::here("results", "tables", "gp_sensitivity.csv"))
  out
}
