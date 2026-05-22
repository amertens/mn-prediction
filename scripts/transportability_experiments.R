# =============================================================================
# scripts/transportability_experiments.R
#
# Develop universal, parsimonious cross-country (LOCO) transportability models
# for mapping micronutrient deficiency from NON-SURVEY proxy predictors only
# (GEE remote sensing + IHME + Malaria Atlas + food security).
#
# Target = Admin-2 survey prevalence (the quantity we map). For each held-out
# country we train on the remaining countries and predict that country's
# Admin-2 prevalences. We score WITHIN-COUNTRY spatial transfer (Pearson /
# Spearman r), error (RMSE/MAE in pp), calibration and national bias, across a
# grid of {prescreen x model x scale x weighting} recipes.
#
# Usage:  Rscript scripts/transportability_experiments.R [outcome1 outcome2 ...]
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(glmnet); library(ranger)
})
set.seed(20260521)
STORE <- "_targets"
COUNTRIES <- c("gambia", "ghana", "malawi", "sierraleone")
OUTCOMES  <- c("child_vitA", "child_iron", "women_vitA", "women_iron")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) OUTCOMES <- args

# ---- harmonization: collapse year-suffixed GEE names to stems -------------
harmonize_gee <- function(x) {
  x <- gsub("_(19|20)[0-9]{2}([0-9]{4})?", "", x)  # strip YYYY and YYYYMMDD
  x <- gsub("__+", "_", x)
  sub("_$", "", x)
}

as_num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  if (is.factor(v)) return(suppressWarnings(as.numeric(as.character(v))))
  suppressWarnings(as.numeric(v))
}

# ---- build enriched Admin-2 covariates for one country --------------------
build_country_admin2 <- function(cn) {
  # 1) GEE from the clean Admin-2 polygon extraction, year-harmonized & averaged
  g <- tar_read_raw(paste0("gee_admin2_", cn), store = STORE)
  gee_cols <- grep("^gee_", names(g), value = TRUE)
  stems <- harmonize_gee(gee_cols)
  gee_mat <- sapply(gee_cols, function(c) as_num(g[[c]]))
  gee_df <- as.data.frame(
    sapply(unique(stems), function(s) {
      cols <- gee_cols[stems == s]
      if (length(cols) == 1) gee_mat[, cols]
      else rowMeans(gee_mat[, cols, drop = FALSE], na.rm = TRUE)
    })
  )
  names(gee_df) <- unique(stems)
  gee_df$Admin2 <- as.character(g$Admin2)

  # 2) IHME / MAP / food-security: aggregate cluster-level merged to Admin2 mean
  m <- tar_read_raw(paste0("merged_", cn), store = STORE)
  dom_cols <- grep("^(ihme_|MAP_|fsec_)", names(m), value = TRUE)
  agg <- data.frame(Admin2 = as.character(m$Admin2))
  for (c in dom_cols) agg[[c]] <- as_num(m[[c]])
  agg <- agg %>%
    group_by(Admin2) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()

  out <- dplyr::full_join(gee_df, agg, by = "Admin2")
  out
}

cat("Building enriched Admin-2 covariates per country...\n")
cov_list <- lapply(COUNTRIES, function(cn) {
  d <- build_country_admin2(cn)
  cat(sprintf("  %s: %d areas x %d covariates\n", cn, nrow(d), ncol(d) - 1))
  d
})
names(cov_list) <- COUNTRIES

# ---- assemble pooled LOCO dataset for one outcome -------------------------
build_pooled <- function(outcome) {
  frames <- list()
  for (cn in COUNTRIES) {
    s <- tryCatch(tar_read_raw(paste0("svy_admin2_", cn, "_", outcome), store = STORE),
                  error = function(e) NULL)
    if (is.null(s) || !any(is.finite(s$svy_prev))) next
    sv <- s %>% transmute(Admin2 = as.character(Admin2), svy_prev,
                          n_svy = if ("n_svy" %in% names(.)) n_svy else NA_real_) %>%
      filter(is.finite(svy_prev))
    d <- dplyr::inner_join(sv, cov_list[[cn]], by = "Admin2")
    if (nrow(d) == 0) next
    d$country <- cn
    frames[[cn]] <- d
  }
  if (length(frames) < 2) return(NULL)
  # common covariates = present (not all-NA) in every country present
  covnames <- function(df) {
    cv <- setdiff(names(df), c("Admin2", "svy_prev", "n_svy", "country"))
    cv[sapply(cv, function(c) any(is.finite(df[[c]])))]
  }
  common <- Reduce(intersect, lapply(frames, covnames))
  pooled <- bind_rows(lapply(frames, function(df)
    df[, c("country", "Admin2", "svy_prev", "n_svy", common), drop = FALSE]))
  list(data = pooled, predictors = common, countries = names(frames))
}

# ---- helpers --------------------------------------------------------------
logit  <- function(p) qlogis(pmin(pmax(p, 0.005), 0.995))
ilogit <- plogis

prep_X <- function(train, test, preds) {
  Xtr <- as.matrix(train[, preds, drop = FALSE])
  Xte <- as.matrix(test[,  preds, drop = FALSE])
  med <- apply(Xtr, 2, function(z) median(z[is.finite(z)]))
  med[!is.finite(med)] <- 0
  for (j in seq_along(preds)) {
    Xtr[!is.finite(Xtr[, j]), j] <- med[j]
    Xte[!is.finite(Xte[, j]), j] <- med[j]
  }
  # drop zero-variance on train
  sdv <- apply(Xtr, 2, sd); keep <- is.finite(sdv) & sdv > 1e-8
  list(Xtr = Xtr[, keep, drop = FALSE], Xte = Xte[, keep, drop = FALSE],
       preds = preds[keep])
}

# correlation prescreen (|cor| with logit prevalence on train), keep top K
screen_corr <- function(Xtr, ytr, K) {
  r <- abs(suppressWarnings(apply(Xtr, 2, function(z) cor(z, ytr))))
  r[!is.finite(r)] <- 0
  sel <- order(r, decreasing = TRUE)[seq_len(min(K, ncol(Xtr)))]
  sort(sel)
}

# ---- model fitters: return prediction on the LOGIT (linear) scale ---------
fit_glmnet <- function(Xtr, ytr, Xte, alpha, w = NULL) {
  if (ncol(Xtr) < 2) return(rep(mean(ytr), nrow(Xte)))
  cv <- tryCatch(cv.glmnet(Xtr, ytr, alpha = alpha, weights = w, nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(ytr), nrow(Xte)))
  as.numeric(predict(cv, newx = Xte, s = "lambda.min"))
}
fit_rf <- function(Xtr, ytr, Xte, w = NULL) {
  df <- data.frame(y = ytr, Xtr)
  rf <- ranger(y ~ ., data = df, num.trees = 800, min.node.size = 5,
               case.weights = w, respect.unordered.factors = TRUE)
  predict(rf, data.frame(Xte))$predictions
}
fit_ols <- function(Xtr, ytr, Xte, w = NULL) {
  df <- data.frame(y = ytr, Xtr); dft <- data.frame(Xte)
  m <- lm(y ~ ., data = df, weights = w)
  as.numeric(predict(m, dft))
}

n_nonzero_glmnet <- function(Xtr, ytr, alpha, w = NULL) {
  cv <- tryCatch(cv.glmnet(Xtr, ytr, alpha = alpha, weights = w, nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv)) return(NA_integer_)
  b <- as.matrix(coef(cv, s = "lambda.min"))
  sum(b[-1, 1] != 0)
}

# within-country centering of a matrix by a country grouping vector
center_by <- function(X, grp) {
  for (g in unique(grp)) {
    idx <- grp == g
    X[idx, ] <- sweep(X[idx, , drop = FALSE], 2,
                      colMeans(X[idx, , drop = FALSE]), "-")
  }
  X
}

# ---- one LOCO run for a given recipe --------------------------------------
run_recipe <- function(pooled, recipe) {
  d <- pooled$data; preds0 <- pooled$predictors
  rows <- list()
  for (ho in pooled$countries) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 10 || nrow(te) < 4) next
    ytr <- logit(tr$svy_prev); yte <- te$svy_prev
    w <- if (recipe$weight) pmax(tr$n_svy, 1) else NULL

    pp <- prep_X(tr, te, preds0)
    Xtr <- pp$Xtr; Xte <- pp$Xte
    # prescreen (on the scale we will fit)
    ytr_fit <- ytr
    if (isTRUE(recipe$center)) {
      Xtr <- center_by(Xtr, tr$country)
      Xte <- sweep(Xte, 2, colMeans(Xte), "-")          # center test by own means
      for (g in unique(tr$country)) {
        idx <- tr$country == g
        ytr_fit[idx] <- ytr[idx] - mean(ytr[idx])
      }
    }
    if (recipe$screen == "corr") {
      sel <- screen_corr(Xtr, ytr_fit, recipe$K)
      Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
    }
    np <- ncol(Xtr)
    level <- mean(ytr)                                   # train grand mean (logit)
    if (recipe$model == "mean") {
      pred <- rep(mean(tr$svy_prev), nrow(te))
    } else {
      plogit <- switch(recipe$model,
        lasso = fit_glmnet(Xtr, ytr_fit, Xte, 1, w),
        enet  = fit_glmnet(Xtr, ytr_fit, Xte, 0.5, w),
        ridge = fit_glmnet(Xtr, ytr_fit, Xte, 0, w),
        rf    = fit_rf(Xtr, ytr_fit, Xte, w),
        ols   = fit_ols(Xtr, ytr_fit, Xte, w)
      )
      if (isTRUE(recipe$center)) plogit <- plogit + level # add back level guess
      pred <- ilogit(plogit)
    }
    pred <- pmin(pmax(pred, 0), 1)
    ok <- is.finite(pred) & is.finite(yte)
    if (sum(ok) < 4) next
    pr <- suppressWarnings(cor(yte[ok], pred[ok]))
    sr <- suppressWarnings(cor(yte[ok], pred[ok], method = "spearman"))
    cal <- tryCatch(coef(lm(yte[ok] ~ pred[ok]))[2], error = function(e) NA)
    rows[[ho]] <- data.frame(
      outcome = pooled$outcome, recipe = recipe$id, held_out = ho,
      n_train = nrow(tr), n_test = sum(ok), n_pred = np,
      pearson_r = round(pr, 3), spearman_r = round(sr, 3),
      rmse_pp = round(sqrt(mean((yte[ok] - pred[ok])^2)) * 100, 2),
      mae_pp  = round(mean(abs(yte[ok] - pred[ok])) * 100, 2),
      bias_pp = round(mean(pred[ok] - yte[ok]) * 100, 2),
      nat_bias_pp = round((mean(pred[ok]) - mean(yte[ok])) * 100, 2),
      calib_slope = round(as.numeric(cal), 2),
      stringsAsFactors = FALSE)
  }
  if (length(rows)) bind_rows(rows) else NULL
}

# ---- recipe grid ----------------------------------------------------------
recipes <- list()
add <- function(id, model, screen = "none", K = NA, weight = TRUE, center = FALSE)
  recipes[[length(recipes) + 1]] <<- list(id = id, model = model, screen = screen,
                                           K = K, weight = weight, center = center)
# --- pooled (no country FE) -------------------------------------------------
add("baseline_mean",    "mean")
add("lasso_all",        "lasso")
add("enet_all",         "enet")
add("ridge_all",        "ridge")
add("rf_all",           "rf")
add("enet_corr30",      "enet",  "corr", 30)
add("enet_corr15",      "enet",  "corr", 15)
add("lasso_corr15",     "lasso", "corr", 15)
add("rf_corr30",        "rf",    "corr", 30)
add("ols_corr8",        "ols",   "corr", 8)
# --- within-country centered (country fixed effects on train) ---------------
add("C_enet_all",       "enet",  center = TRUE)
add("C_ridge_all",      "ridge", center = TRUE)
add("C_enet_corr30",    "enet",  "corr", 30, center = TRUE)
add("C_enet_corr15",    "enet",  "corr", 15, center = TRUE)
add("C_lasso_corr15",   "lasso", "corr", 15, center = TRUE)
add("C_rf_corr30",      "rf",    "corr", 30, center = TRUE)
add("C_ols_corr8",      "ols",   "corr", 8,  center = TRUE)

# ---- run ------------------------------------------------------------------
all_fold <- list()
for (oc in OUTCOMES) {
  cat(sprintf("\n================= %s =================\n", oc))
  pooled <- build_pooled(oc)
  if (is.null(pooled)) { cat("  insufficient data\n"); next }
  pooled$outcome <- oc
  cat(sprintf("  pooled: %d areas, %d common predictors, countries: %s\n",
              nrow(pooled$data), length(pooled$predictors),
              paste(pooled$countries, collapse = ", ")))
  for (rc in recipes) {
    res <- tryCatch(run_recipe(pooled, rc), error = function(e) {
      cat(sprintf("    [%s] ERROR: %s\n", rc$id, e$message)); NULL })
    if (!is.null(res)) all_fold[[paste(oc, rc$id)]] <- res
  }
}

# ---- within-country CV ceiling (same-country training available) ----------
# 5-fold CV within each country (enet on top-15 corr-screened, logit scale).
# Shows the achievable spatial r when NOT transferring across borders.
ceiling_rows <- list()
for (oc in OUTCOMES) {
  pooled <- build_pooled(oc); if (is.null(pooled)) next
  d <- pooled$data; preds0 <- pooled$predictors
  for (cn in pooled$countries) {
    dc <- d[d$country == cn, ]
    if (nrow(dc) < 20) next
    folds <- sample(rep(1:5, length.out = nrow(dc)))
    oof <- rep(NA_real_, nrow(dc))
    for (k in 1:5) {
      tr <- dc[folds != k, ]; te <- dc[folds == k, ]
      ytr <- logit(tr$svy_prev)
      pp <- prep_X(tr, te, preds0)
      sel <- screen_corr(pp$Xtr, ytr, 15)
      oof[folds == k] <- ilogit(fit_glmnet(pp$Xtr[, sel, drop = FALSE], ytr,
                                            pp$Xte[, sel, drop = FALSE], 0.5,
                                            pmax(tr$n_svy, 1)))
    }
    ok <- is.finite(oof) & is.finite(dc$svy_prev)
    if (sum(ok) < 8) next
    ceiling_rows[[paste(oc, cn)]] <- data.frame(
      outcome = oc, country = cn, n = sum(ok),
      cv_pearson = round(cor(dc$svy_prev[ok], oof[ok]), 3),
      cv_spearman = round(cor(dc$svy_prev[ok], oof[ok], method = "spearman"), 3),
      cv_rmse_pp = round(sqrt(mean((dc$svy_prev[ok] - oof[ok])^2)) * 100, 2),
      stringsAsFactors = FALSE)
  }
}
ceiling_df <- bind_rows(ceiling_rows)
if (nrow(ceiling_df)) {
  write.csv(ceiling_df, "results/tables/transportability_within_country_ceiling.csv",
            row.names = FALSE)
}

fold_df <- bind_rows(all_fold)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(fold_df, "results/tables/transportability_experiments_folds.csv", row.names = FALSE)

# ---- summary: average across held-out countries ---------------------------
summary_df <- fold_df %>%
  group_by(outcome, recipe) %>%
  summarise(n_folds = n(),
            mean_pearson = round(mean(pearson_r, na.rm = TRUE), 3),
            mean_spearman = round(mean(spearman_r, na.rm = TRUE), 3),
            mean_rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
            mean_mae_pp  = round(mean(mae_pp, na.rm = TRUE), 2),
            mean_abs_natbias = round(mean(abs(nat_bias_pp), na.rm = TRUE), 2),
            med_npred = round(median(n_pred), 0),
            .groups = "drop") %>%
  arrange(outcome, desc(mean_spearman))
write.csv(summary_df, "results/tables/transportability_experiments_summary.csv", row.names = FALSE)

cat("\n\n############ SUMMARY (sorted by mean Spearman within outcome) ############\n")
options(width = 200)
print(as.data.frame(summary_df))

cat("\n\n############ Best recipe per outcome (by mean Spearman) ############\n")
best <- summary_df %>% group_by(outcome) %>% slice_max(mean_spearman, n = 1) %>% ungroup()
print(as.data.frame(best))

cat("\n\n############ Recipe ranking AVERAGED across outcomes ############\n")
overall <- summary_df %>% group_by(recipe) %>%
  summarise(mean_spearman = round(mean(mean_spearman), 3),
            mean_pearson = round(mean(mean_pearson), 3),
            mean_rmse_pp = round(mean(mean_rmse_pp), 2),
            med_npred = round(median(med_npred), 0), .groups = "drop") %>%
  arrange(desc(mean_spearman))
print(as.data.frame(overall))

cat("\n\n############ Within-country CV ceiling (same-country training) ############\n")
if (exists("ceiling_df") && nrow(ceiling_df)) {
  print(as.data.frame(ceiling_df))
  cat(sprintf("\n  Mean within-country CV: Pearson=%.3f Spearman=%.3f (vs LOCO best Spearman=%.3f)\n",
              mean(ceiling_df$cv_pearson, na.rm = TRUE),
              mean(ceiling_df$cv_spearman, na.rm = TRUE),
              max(overall$mean_spearman, na.rm = TRUE)))
}
cat("\nDONE\n")
