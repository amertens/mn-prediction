# =============================================================================
# sandbox/13_univariate_loco_ranking.R
#
# For each of the 152 common GEE variables, fit a single-variable LOCO
# regression model (per outcome) and record the resulting LOCO Pearson r.
# Then ask: does ANY single GEE variable beat the spatial GAM (r = 0.285)
# on average?
#
# This is the most direct test of "is there ANY feature in the GEE set
# carrying transportable signal beyond geography?" If no single variable
# beats spatial GAM, the answer is no; if some do, those are the candidates
# to fold into a spatial-plus-feature stacked model.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

pooled_all <- load_all_pooled()
vars <- pooled_all[[1]]$common_gee_vars
countries <- pooled_all[[1]]$country_names

cat(sprintf("evaluating %d variables x %d outcomes x %d holdouts = %d univariate fits\n",
            length(vars), length(pooled_all), length(countries),
            length(vars) * length(pooled_all) * length(countries)))

uni_predict <- function(train, test, vars) {
  # vars here is a single variable.
  imp <- impute_median(train, test, vars)
  Y   <- train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  df  <- data.frame(Y = Y, x = imp$train[[vars]])
  fit <- tryCatch(stats::lm(Y ~ x, data = df, weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  as.numeric(stats::predict(fit, newdata = data.frame(x = imp$test[[vars]])))
}

all_rows <- list()
t0 <- Sys.time()
for (oc in names(pooled_all)) {
  cat(sprintf("\n[uni] outcome %s ...\n", oc))
  p <- pooled_all[[oc]]
  per_var <- list()
  for (i in seq_along(vars)) {
    v <- vars[i]
    # Skip vars that are NA in the pooled data for this outcome.
    if (sum(!is.finite(p$pooled_data[[v]])) > nrow(p$pooled_data) * 0.5) next
    res <- loco_evaluate(
      pooled = list(pooled_data = p$pooled_data,
                     common_gee_vars = v,
                     country_names = countries),
      predict_fn = function(train, test, vars) uni_predict(train, test, v),
      outcome = oc, method_name = v)
    if (nrow(res) > 0) {
      res$variable <- v
      per_var[[v]] <- res
    }
    if (i %% 50 == 0) cat(sprintf("  ... %d/%d (elapsed %.1fs)\n", i, length(vars),
                                    as.numeric(Sys.time()-t0, units = "secs")))
  }
  all_rows[[oc]] <- dplyr::bind_rows(per_var)
}
res <- dplyr::bind_rows(all_rows)

# Per-variable summary across all (outcome × holdout) splits.
summary_per_var <- res |>
  dplyr::group_by(variable) |>
  dplyr::summarise(
    n_splits = dplyr::n(),
    mean_r   = round(mean(pearson_r, na.rm = TRUE),   3),
    median_r = round(median(pearson_r, na.rm = TRUE), 3),
    wins     = sum(pearson_r > 0.1, na.rm = TRUE),
    losses   = sum(pearson_r < -0.1, na.rm = TRUE),
    mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_r))

cat("\n===== top 20 single-variable LOCO predictors =====\n")
print(as.data.frame(head(summary_per_var, 20)), row.names = FALSE)

cat("\n===== bottom 5 (anti-predictive) =====\n")
print(as.data.frame(tail(summary_per_var, 5)), row.names = FALSE)

cat(sprintf("\nbest single variable mean_r = %.3f (target to beat: spatial GAM = 0.285)\n",
            max(summary_per_var$mean_r, na.rm = TRUE)))
cat(sprintf("variables with mean_r > 0.10: %d / %d\n",
            sum(summary_per_var$mean_r > 0.10, na.rm = TRUE), nrow(summary_per_var)))
cat(sprintf("variables with mean_r > 0.20: %d / %d\n",
            sum(summary_per_var$mean_r > 0.20, na.rm = TRUE), nrow(summary_per_var)))

readr::write_csv(summary_per_var,
                  file.path(SANDBOX_RESULTS, "13_univariate_loco.csv"))
cat(sprintf("\nwritten: %s\n",
            file.path(SANDBOX_RESULTS, "13_univariate_loco.csv")))
