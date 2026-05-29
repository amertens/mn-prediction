# =============================================================================
# sandbox/16_higher_k_and_sl_comparison.R
#
# Two follow-ups to the spatial_plus_top8 winner (mean LOCO r = 0.330):
#
#   Part A. Extend the K sweep. We tested K in {3..8} and K=8 was best.
#           Does K = 10, 12, 15, 20, 30, 40, 60 continue to climb, or
#           does it plateau / degrade? Tests whether the soil-feature
#           pool itself is the bottleneck.
#
#   Part B. With the WINNING feature set (lon, lat + top-K soil), does a
#           SuperLearner (or any non-linear learner: ranger, xgboost)
#           outperform the simple linear-additions GAM? If even SL on
#           the restricted features can't beat OLS-on-restricted, that
#           is direct evidence that the relationships are genuinely
#           linear/smooth and the prior 17-hour SL was failing because
#           of feature noise, not modelling capacity.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "config.R"))
cc_list <- get_country_configs()

pooled_all <- load_all_pooled()
example_svy_list <- pooled_all[[1]]$svy_list
for (oc in names(pooled_all))
  pooled_all[[oc]]$pooled_data <-
    add_admin2_centroids(pooled_all[[oc]]$pooled_data, cc_list, example_svy_list)

# Pull the all-source univariate ranking. (No need to refit script 13 - cached.)
all_uni <- read.csv(file.path(SANDBOX_RESULTS, "13_univariate_loco.csv"))
all_uni <- all_uni[order(-all_uni$mean_r), ]

# ── Part A: extend K sweep to large K ───────────────────────────────────────
make_spatial_plus <- function(extras) function(train, test, vars) {
  o <- fit_predict_spatial_coords(train, test, vars, extra_covars = extras)
  if (is.null(o)) NULL else o$pred
}

ks <- c(8, 10, 12, 15, 20, 30, 40, 60)
sweep_rows <- list()
for (K in ks) {
  extras <- head(all_uni$variable, K)
  name <- sprintf("spatial_top%d", K)
  cat(sprintf("\n[K-sweep] %s (%d extras) ...\n", name, length(extras)))
  r <- loco_evaluate_all(pooled_all, make_spatial_plus(extras), method_name = name)
  sweep_rows[[name]] <- r
}
sweep <- dplyr::bind_rows(sweep_rows)
summ_A <- sweep |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
                    median_r = round(median(pearson_r, na.rm = TRUE), 3),
                    wins = sum(pearson_r > 0.1, na.rm = TRUE),
                    n_splits = sum(!is.na(pearson_r)),
                    .groups = "drop")
# Preserve K order
summ_A$K <- as.integer(sub("spatial_top", "", summ_A$method))
summ_A <- summ_A[order(summ_A$K), ]
cat("\n===== Part A: extended K sweep =====\n")
print(as.data.frame(summ_A), row.names = FALSE)

# ── Part B: compare learners on the WINNER's feature set (top-8 soil) ───────
top_extras <- head(all_uni$variable, 8)
cat(sprintf("\n===== Part B: learner comparison on (lon, lat + top-8 soil) =====\n"))
cat("features:", paste(top_extras, collapse = "\n         "), "\n\n")

# Learner 1: GAM with linear extras (= our current winning method)
fn_gam <- function(train, test, vars) {
  o <- fit_predict_spatial_coords(train, test, vars, extra_covars = top_extras)
  if (is.null(o)) NULL else o$pred
}

# Learner 2: ranger (random forest) on (lon, lat, top-8)
fn_ranger <- function(train, test, vars) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  imp <- impute_median(train, test, top_extras)
  fit_df <- data.frame(Y = train$svy_prev, lon = train$lon, lat = train$lat,
                        imp$train)
  test_df <- data.frame(            lon = test$lon, lat = test$lat,
                          imp$test)
  fit <- tryCatch(ranger::ranger(Y ~ ., data = fit_df, num.trees = 500,
                                    min.node.size = 5,
                                    case.weights = pmax(train$n_svy %||% rep(1, nrow(fit_df)), 1)),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pmin(pmax(predict(fit, data = test_df)$predictions, 0), 1)
}

# Learner 3: xgboost
fn_xgb <- function(train, test, vars) {
  if (!requireNamespace("xgboost", quietly = TRUE)) return(NULL)
  imp <- impute_median(train, test, top_extras)
  Xtr <- as.matrix(data.frame(lon = train$lon, lat = train$lat, imp$train))
  Xte <- as.matrix(data.frame(lon = test$lon,  lat = test$lat,  imp$test))
  dtr <- xgboost::xgb.DMatrix(Xtr, label = train$svy_prev,
                                weight = pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1))
  fit <- tryCatch(xgboost::xgb.train(
                     data = dtr, nrounds = 300, max_depth = 4, eta = 0.05,
                     subsample = 0.8, colsample_bytree = 0.8,
                     objective = "reg:squarederror", verbose = 0),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pmin(pmax(predict(fit, Xte), 0), 1)
}

# Learner 4: SuperLearner (mlr3superlearner) on the small feature set
fn_sl <- function(train, test, vars) {
  if (!requireNamespace("mlr3superlearner", quietly = TRUE)) return(NULL)
  imp <- impute_median(train, test, top_extras)
  fit_df <- data.frame(Y = train$svy_prev, lon = train$lon, lat = train$lat,
                        imp$train)
  test_df <- data.frame(Y = 0, lon = test$lon, lat = test$lat, imp$test)
  lib <- list(
    list("glmnet", alpha = 0.5, id = "elastic_net"),
    list("glmnet", alpha = 1,   id = "lasso"),
    list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
    list("xgboost", max_depth = 4, eta = 0.05, nrounds = 300,
          subsample = 0.8, colsample_bytree = 0.8, id = "xgb"),
    list("mean", id = "mean")
  )
  fit <- tryCatch(suppressWarnings(suppressMessages(
    mlr3superlearner::mlr3superlearner(data = fit_df, target = "Y",
                                          library = lib,
                                          outcome_type = "continuous",
                                          folds = min(nrow(fit_df), 5L)))),
    error = function(e) { cat("    SL fail:", conditionMessage(e),"\n"); NULL })
  if (is.null(fit)) return(NULL)
  p <- tryCatch(as.numeric(stats::predict(fit, test_df)),
                error = function(e) NULL)
  if (is.null(p)) return(NULL)
  pmin(pmax(p, 0), 1)
}

learners <- list(
  spatial_plus_linear_OLS = fn_gam,
  spatial_plus_ranger     = fn_ranger,
  spatial_plus_xgb        = fn_xgb,
  spatial_plus_SL         = fn_sl
)

res_B <- list()
for (m in names(learners)) {
  cat(sprintf("\n[B] %s ...\n", m))
  r <- loco_evaluate_all(pooled_all, learners[[m]], method_name = m)
  res_B[[m]] <- r
}
res_B_df <- dplyr::bind_rows(res_B)
summ_B <- res_B_df |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
                    median_r = round(median(pearson_r, na.rm = TRUE), 3),
                    wins = sum(pearson_r > 0.1, na.rm = TRUE),
                    n_splits = sum(!is.na(pearson_r)),
                    mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_r))
cat("\n===== Part B: learner comparison on identical 10-feature input =====\n")
print(as.data.frame(summ_B), row.names = FALSE)

readr::write_csv(sweep, file.path(SANDBOX_RESULTS, "16a_higher_k_sweep.csv"))
readr::write_csv(res_B_df,
                  file.path(SANDBOX_RESULTS, "16b_learner_comparison.csv"))
cat(sprintf("\nDONE -- see sandbox/results/16{a,b}_*.csv\n"))
