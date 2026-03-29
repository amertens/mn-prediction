# =============================================================================
# scripts/variable_importance.R
#
# Variable importance for top-performing models.
# Fits a standalone ranger model on the same training data used by the SL,
# then extracts impurity-based variable importance. This avoids the BART
# serialization issues that make direct extraction from SL objects infeasible.
# Saves results to results/tables/variable_importance_*.csv
# =============================================================================

library(targets)
library(data.table)
library(dplyr)
library(ranger)
library(here)

TARGETS_STORE <- here("_targets_full")

#' Compute variable importance by fitting a ranger model on the SL training data
#'
#' @param sl_fit The sl_fit object from the targets store
#' @param use_binary Whether this is a binary outcome
#' @return data.frame with variable and importance columns
compute_ranger_importance <- function(sl_fit, use_binary = TRUE) {

  fit_obj <- if (use_binary) sl_fit$bin_fit else sl_fit$cont_fit
  if (is.null(fit_obj)) stop("No fitted model found")

  covars <- fit_obj$Xvars
  train_data <- fit_obj$train_data

  if (is.null(train_data)) stop("No train_data stored")

  Y <- as.numeric(fit_obj$res$Y)
  X <- as.data.frame(train_data)[, covars, drop = FALSE]

  cat(sprintf("  Fitting ranger on %d obs x %d vars...\n", nrow(X), ncol(X)))

  # Fit ranger with impurity importance
  if (use_binary) {
    rf <- ranger::ranger(
      x = X, y = factor(Y),
      num.trees = 1000,
      importance = "impurity",
      min.node.size = 10,
      probability = TRUE
    )
  } else {
    rf <- ranger::ranger(
      x = X, y = Y,
      num.trees = 1000,
      importance = "impurity",
      min.node.size = 10
    )
  }

  imp <- rf$variable.importance
  vi <- data.frame(
    variable = names(imp),
    importance = as.numeric(imp),
    stringsAsFactors = FALSE
  ) |>
    arrange(desc(importance)) |>
    mutate(importance_norm = importance / max(importance))

  cat(sprintf("  Ranger OOB error: %.4f\n", rf$prediction.error))
  vi
}


# ── Run for top-2 models ────────────────────────────────────────────────────

# Top-2 by AUC with reasonable event counts:
#   1. Ghana child_iron  (AUC=0.828, 208 events)
#   2. Ghana women_b12   (AUC=0.855, 40 events)
top_models <- list(
  list(target = "sl_fit_ghana_child_iron",   country = "Ghana", outcome = "child_iron"),
  list(target = "sl_fit_ghana_women_b12",    country = "Ghana", outcome = "women_b12")
)

out_dir <- here("results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (m in top_models) {
  cat(sprintf("\n=== %s — %s ===\n", m$country, m$outcome))

  sl_fit <- tryCatch(
    targets::tar_read_raw(m$target, store = TARGETS_STORE),
    error = function(e) {
      cat(sprintf("  Failed to load %s: %s\n", m$target, e$message))
      NULL
    }
  )
  if (is.null(sl_fit)) next

  vi <- tryCatch(
    compute_ranger_importance(sl_fit, use_binary = TRUE),
    error = function(e) {
      cat(sprintf("  VI computation failed: %s\n", e$message))
      NULL
    }
  )
  if (is.null(vi)) next

  vi$country <- m$country
  vi$outcome <- m$outcome

  outfile <- file.path(out_dir,
                        sprintf("variable_importance_%s_%s.csv",
                                tolower(m$country), m$outcome))
  write.csv(vi, outfile, row.names = FALSE)
  cat(sprintf("  Saved %d variables to %s\n", nrow(vi), outfile))
  cat("  Top 10:\n")
  print(head(vi[, c("variable", "importance", "importance_norm")], 10))
}

cat("\nDone.\n")
