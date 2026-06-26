# =============================================================================
# _tmp_verify_oof.R — VERIFICATION ONLY (no writes to results/)
#
# Confirms that fit_area_sl()'s reported loo_preds (== predict(sl_fit, fit_df))
# are IN-SAMPLE refits, not out-of-fold, by comparing against:
#   (a) genuine 10-fold OOF (train-fold-only median imputation)
#   (b) genuine LOO
# for Ghana + Malawi child_iron (MSE / continuous loss — the headline variant).
# =============================================================================
suppressWarnings(suppressMessages({
  library(targets); library(dplyr)
  library(mlr3superlearner); library(mlr3extralearners)
}))

STORE   <- "_targets_full"
OUTCOME <- "child_iron"
SEED    <- 20260615

SL_LIB <- list(
  list("glmnet", alpha = 1,   id = "lasso"),
  list("glmnet", alpha = 0,   id = "ridge"),
  list("glmnet", alpha = 0.5, id = "elastic_net"),
  list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
  list("xgboost", max_depth = 2, eta = 0.05, nrounds = 200,
       min_child_weight = 3, subsample = 0.8, id = "xgb")
)

source("R/area_level_comparison.R")

load_country <- function(cn) {
  svy <- tar_read_raw(paste0("svy_admin2_", cn, "_", OUTCOME), store = STORE)
  gee <- tar_read_raw(paste0("gee_admin2_", cn), store = STORE)
  df <- dplyr::inner_join(
    svy |> dplyr::select(Admin2, svy_prev, dplyr::any_of("n_svy")),
    gee, by = "Admin2") |>
    dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))
  gee_vars <- setdiff(names(gee), c("Admin1", "Admin2"))
  valid <- gee_vars[vapply(gee_vars, function(v) {
    x <- df[[v]]; sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  list(df = df, vars = valid, svy = svy, gee = gee)
}

impute_with <- function(target, train, vars) {
  for (v in vars) {
    x <- train[[v]]; if (inherits(x, "haven_labelled")) x <- as.double(unclass(x))
    med <- median(x[is.finite(x)], na.rm = TRUE)
    tx <- target[[v]]; if (inherits(tx, "haven_labelled")) tx <- as.double(unclass(tx))
    tx[!is.finite(tx)] <- med
    target[[v]] <- tx
  }
  target
}

fit_one <- function(tr_i, vars) {
  tryCatch(mlr3superlearner::mlr3superlearner(
    data = data.frame(Y = tr_i$svy_prev, tr_i[, vars, drop = FALSE]),
    target = "Y", library = SL_LIB, outcome_type = "continuous",
    folds = min(5L, nrow(tr_i) - 1L)), error = function(e) {cat("   fit err:", conditionMessage(e), "\n"); NULL})
}
pred_one <- function(fit, te_i, vars) {
  p <- tryCatch(as.numeric(predict(fit, data.frame(te_i[, vars, drop = FALSE]))),
                error = function(e) NULL)
  if (is.null(p)) return(rep(NA_real_, nrow(te_i)))
  pmin(pmax(p, 0), 1)
}
rr <- function(y, p) { ok <- is.finite(y) & is.finite(p); round(cor(y[ok], p[ok]), 3) }
mae <- function(y, p) { ok <- is.finite(y) & is.finite(p); round(mean(abs(y[ok]-p[ok]))*100, 2) }

run <- function(cn, label) {
  cat(sprintf("\n========== %s child_iron (continuous/MSE) ==========\n", label))
  d <- load_country(cn); df <- d$df; vars <- d$vars
  n <- nrow(df)
  cat(sprintf("  n=%d areas, p=%d GEE vars\n", n, length(vars)))

  ## (0) What fit_area_sl reports today (in-sample refit predict)
  fa <- fit_area_sl(d$svy, d$gee, outcome_type = "continuous")
  insample_r   <- fa$metrics$pearson_r
  insample_mae <- fa$metrics$mae_pp
  cat(sprintf("  [fit_area_sl reported]  r=%.3f  MAE=%.2f pp  (uses predict on training)\n",
              insample_r, insample_mae))

  ## (1) Direct demonstration: full-data fit, predict on training == in-sample
  df_i_full <- impute_with(df, df, vars)
  fit_full  <- fit_one(df_i_full, vars)
  p_train   <- pred_one(fit_full, df_i_full, vars)
  cat(sprintf("  [predict() on training] r=%.3f  MAE=%.2f pp\n",
              rr(df$svy_prev, p_train), mae(df$svy_prev, p_train)))

  ## (2) Genuine 10-fold OOF (train-fold-only imputation)
  set.seed(SEED); K <- 10L
  fold <- sample(rep_len(seq_len(K), n)); oof <- rep(NA_real_, n)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    tr_i <- impute_with(df[tr,,drop=FALSE], df[tr,,drop=FALSE], vars)
    te_i <- impute_with(df[te,,drop=FALSE], df[tr,,drop=FALSE], vars)
    fit <- fit_one(tr_i, vars); if (is.null(fit)) next
    oof[te] <- pred_one(fit, te_i, vars)
  }
  cat(sprintf("  [GENUINE 10-fold OOF]   r=%.3f  MAE=%.2f pp\n",
              rr(df$svy_prev, oof), mae(df$svy_prev, oof)))

  ## (3) Genuine LOO
  loo <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    tr_i <- impute_with(df[tr,,drop=FALSE], df[tr,,drop=FALSE], vars)
    te_i <- impute_with(df[i,,drop=FALSE],  df[tr,,drop=FALSE], vars)
    fit <- fit_one(tr_i, vars); if (is.null(fit)) next
    loo[i] <- pred_one(fit, te_i, vars)
  }
  cat(sprintf("  [GENUINE LOO]           r=%.3f  MAE=%.2f pp\n",
              rr(df$svy_prev, loo), mae(df$svy_prev, loo)))
}

run("ghana",  "Ghana")
run("malawi", "Malawi")
cat("\nDONE\n")
