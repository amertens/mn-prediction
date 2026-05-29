# =============================================================================
# sandbox/00_setup.R
#
# Common data-loading + LOCO evaluator for every sandbox experiment.
# Source this first from each script.
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(tidyr); library(here); library(sf)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

STORE <- here::here("_targets_full")
SANDBOX_CACHE <- here::here("sandbox", "cache")
SANDBOX_RESULTS <- here::here("sandbox", "results")
dir.create(SANDBOX_CACHE,   showWarnings = FALSE, recursive = TRUE)
dir.create(SANDBOX_RESULTS, showWarnings = FALSE, recursive = TRUE)

COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi")
COUNTRY_LABELS <- c(Gambia = "gambia", Ghana = "ghana",
                     SierraLeone = "sierraleone", Malawi = "malawi")
OUTCOMES <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

# --------------------------------------------------------------------------- #
# Cached pooled data loader (per outcome). Cache makes subsequent scripts ~10x
# faster than going back to tar_read every time.
# --------------------------------------------------------------------------- #
load_pooled <- function(outcome, refresh = FALSE) {
  cache_path <- file.path(SANDBOX_CACHE, paste0("pooled_", outcome, ".rds"))
  if (!refresh && file.exists(cache_path)) return(readRDS(cache_path))
  source(here::here("R", "area_level_comparison.R"))
  gee_list <- setNames(lapply(COUNTRIES, function(c)
    tryCatch(targets::tar_read_raw(paste0("gee_admin2_", c), store = STORE),
              error = function(e) NULL)),
    names(COUNTRY_LABELS))
  svy_list <- setNames(lapply(COUNTRIES, function(c)
    tryCatch(targets::tar_read_raw(paste0("svy_admin2_", c, "_", outcome),
                                     store = STORE),
              error = function(e) NULL)),
    names(COUNTRY_LABELS))
  # Drop countries without this outcome OR without GEE data (Côte d'Ivoire etc.).
  ok <- !vapply(svy_list, is.null, logical(1)) & !vapply(gee_list, is.null, logical(1))
  if (sum(ok) < 2) {
    warning(sprintf("Only %d countries have data for outcome '%s'", sum(ok), outcome))
    return(NULL)
  }
  svy_list <- svy_list[ok]; gee_list <- gee_list[ok]
  pooled <- build_area_loco_dataset(svy_list, gee_list)
  pooled$svy_list <- svy_list
  saveRDS(pooled, cache_path)
  pooled
}

# Convenience: load all 4 outcomes' pooled data once.
load_all_pooled <- function(refresh = FALSE) {
  out <- list()
  for (oc in OUTCOMES) out[[oc]] <- load_pooled(oc, refresh = refresh)
  out
}

# --------------------------------------------------------------------------- #
# Unified LOCO evaluator.
#
# Takes a prediction function `predict_fn(train, test, vars, ...)` that returns
# numeric predictions for `test$svy_prev`, runs it across all 4 LOCO holdouts,
# returns a tidy 4-row data.frame with held_out, Pearson r, Spearman ρ,
# MAE (pp), RMSE (pp), bias (pp), n_train, n_test.
# --------------------------------------------------------------------------- #
loco_evaluate <- function(pooled, predict_fn, outcome = NA_character_,
                            method_name = "method",
                            country_names = c("Gambia","Ghana","SierraLeone","Malawi"),
                            ...) {
  pd <- pooled$pooled_data
  vars <- pooled$common_gee_vars
  res <- list()
  for (held_out in country_names) {
    train <- pd[pd$country != held_out, , drop = FALSE]
    test  <- pd[pd$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    p <- tryCatch(predict_fn(train, test, vars, ...),
                  error = function(e) { cat("  err:", conditionMessage(e), "\n"); NULL })
    if (is.null(p) || length(p) != nrow(test)) next
    obs <- test$svy_prev
    ok <- is.finite(obs) & is.finite(p)
    if (sum(ok) < 3) next
    obs <- obs[ok]; p <- p[ok]
    res[[held_out]] <- data.frame(
      method   = method_name, outcome  = outcome,
      held_out = held_out,    n_train  = nrow(train), n_test = length(obs),
      pearson_r  = suppressWarnings(round(cor(obs, p),                    3)),
      spearman_r = suppressWarnings(round(cor(obs, p, method = "spearman"), 3)),
      mae_pp     = round(mean(abs(obs - p)) * 100, 2),
      rmse_pp    = round(sqrt(mean((obs - p)^2)) * 100, 2),
      mean_bias_pp = round(mean(p - obs) * 100, 2),
      stringsAsFactors = FALSE
    )
  }
  if (length(res) == 0) return(data.frame())
  dplyr::bind_rows(res)
}

# Run a predict_fn across all 4 outcomes and tally a summary.
loco_evaluate_all <- function(pooled_all, predict_fn, method_name, ...) {
  rows <- list()
  for (oc in names(pooled_all)) {
    r <- loco_evaluate(pooled_all[[oc]], predict_fn, outcome = oc,
                        method_name = method_name, ...)
    if (nrow(r) > 0) rows[[oc]] <- r
  }
  dplyr::bind_rows(rows)
}

# Pretty print a results table grouped by outcome.
print_loco <- function(res) {
  for (oc in unique(res$outcome)) {
    cat(sprintf("\n--- %s ---\n", oc))
    print(res[res$outcome == oc,
              c("held_out","pearson_r","spearman_r","mae_pp","mean_bias_pp")],
          row.names = FALSE)
  }
  cat(sprintf("\nmean Pearson r = %.3f | median = %.3f | n = %d\n",
              mean(res$pearson_r, na.rm = TRUE),
              median(res$pearson_r, na.rm = TRUE), nrow(res)))
}

# --------------------------------------------------------------------------- #
# Variable screening: drop near-constant and high-NA covariates on TRAINING
# data, return a clean character vector of usable variables.
# --------------------------------------------------------------------------- #
screen_vars <- function(train, vars, max_na_frac = 0.5) {
  vars[vapply(vars, function(v) {
    x <- train[[v]]; n <- length(x)
    na_frac <- sum(!is.finite(x)) / n
    valid <- x[is.finite(x)]
    na_frac < max_na_frac && length(unique(valid)) > 1
  }, logical(1))]
}

# Median-impute NAs in (train, test) using train medians.
impute_median <- function(train, test, vars) {
  Xtr <- as.data.frame(train[, vars, drop = FALSE])
  Xte <- as.data.frame(test [, vars, drop = FALSE])
  med <- vapply(Xtr, median, numeric(1), na.rm = TRUE)
  for (v in vars) {
    btr <- !is.finite(Xtr[[v]]); if (any(btr)) Xtr[[v]][btr] <- med[v]
    bte <- !is.finite(Xte[[v]]); if (any(bte)) Xte[[v]][bte] <- med[v]
  }
  list(train = Xtr, test = Xte, medians = med)
}

cat("[sandbox] setup loaded: load_pooled(), load_all_pooled(), loco_evaluate(),",
    "loco_evaluate_all(), screen_vars(), impute_median()\n")
