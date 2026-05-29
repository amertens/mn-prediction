# =============================================================================
# sandbox/14_spatial_plus_topgee.R
#
# Stack spatial GAM (the new winner) with top-K univariate-LOCO GEE
# predictors from script 13. Two flavours:
#
#   (a) "spatial_plus_features"
#       Single GAM with thin-plate spline on (lon, lat) + linear top-K
#       GEE features.
#
#   (b) "spatial_then_residual_en"
#       Stage 1: spatial GAM on coords only.
#       Stage 2: elastic-net on top-K GEE features predicting residuals.
#       Final prediction = stage1 + stage2.
#
# Compare to the spatial-only baseline (0.285) and the H6 combined-filter
# (0.224). The question: can we recover ANY useful signal from GEE on
# top of geography?
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

# Add centroids to pooled data.
source(here::here("R", "config.R"))
cc_list <- get_country_configs()

pooled_all <- load_all_pooled()
example_svy_list <- pooled_all[[1]]$svy_list
for (oc in names(pooled_all)) {
  pooled_all[[oc]]$pooled_data <- add_admin2_centroids(
    pooled_all[[oc]]$pooled_data, cc_list, example_svy_list)
}

# Load top-K from univariate ranking (script 13). Re-rank by mean_r > 0.
top_uni <- tryCatch(
  read.csv(file.path(SANDBOX_RESULTS, "13_univariate_loco.csv")),
  error = function(e) NULL)
if (is.null(top_uni)) {
  stop("Need sandbox/results/13_univariate_loco.csv — run script 13 first")
}
top_uni <- top_uni[order(-top_uni$mean_r), ]
cat("top 10 single-variable predictors:\n")
print(head(top_uni[, c("variable","mean_r","wins","losses")], 10), row.names = FALSE)

top5  <- head(top_uni$variable, 5)
top10 <- head(top_uni$variable, 10)
top20 <- head(top_uni$variable, 20)

# (a) Spatial + linear top features.
fit_spatial_plus_features <- function(top_vars) function(train, test, vars) {
  out <- fit_predict_spatial_coords(train, test, vars,
                                       extra_covars = top_vars)
  if (is.null(out)) NULL else out$pred
}

# (b) Spatial baseline + elastic-net on residuals over top features.
fit_spatial_then_resid <- function(top_vars) function(train, test, vars) {
  s <- fit_predict_spatial_coords(train, test, vars)
  if (is.null(s)) return(NULL)
  resid_tr <- train$svy_prev - s$train_pred
  imp <- impute_median(train, test, top_vars)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  wt  <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  if (ncol(Xtr) < 2) return(s$pred)
  fit <- tryCatch(glmnet::cv.glmnet(Xtr, resid_tr, alpha = 0.5,
                                       nfolds = min(nrow(Xtr), 10L),
                                       weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(s$pred)
  r_te <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.min"))
  pmin(pmax(s$pred + r_te, 0), 1)
}

res <- list()
configs <- list(
  list(name = "spatial_only",                    fn = function(train, test, vars) {
    o <- fit_predict_spatial_coords(train, test, vars); if (is.null(o)) NULL else o$pred }),
  list(name = "spatial_plus_top5_linear",        fn = fit_spatial_plus_features(top5)),
  list(name = "spatial_plus_top10_linear",       fn = fit_spatial_plus_features(top10)),
  list(name = "spatial_then_resid_en_top5",      fn = fit_spatial_then_resid(top5)),
  list(name = "spatial_then_resid_en_top10",     fn = fit_spatial_then_resid(top10)),
  list(name = "spatial_then_resid_en_top20",     fn = fit_spatial_then_resid(top20))
)
for (cfg in configs) {
  cat(sprintf("\n[H7] %s ...\n", cfg$name))
  r <- loco_evaluate_all(pooled_all, cfg$fn, method_name = cfg$name)
  res[[cfg$name]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "14_spatial_plus_topgee.csv"))
cat(sprintf("\nwritten: %s\n",
            file.path(SANDBOX_RESULTS, "14_spatial_plus_topgee.csv")))
