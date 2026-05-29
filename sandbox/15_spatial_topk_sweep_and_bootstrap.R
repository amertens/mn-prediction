# =============================================================================
# sandbox/15_spatial_topk_sweep_and_bootstrap.R
#
# Three-part follow-up to the spatial_plus_top5 winner from script 14:
#
#   Part A. K sweep on "all-source" top-K (top by univariate LOCO mean r).
#           Test K = 3, 4, 5, 6, 7, 8 to find the optimum.
#
#   Part B. K sweep on "SoilGrids-only" top-K to strengthen the mechanism
#           story (and avoid noise from low-signal non-soil vars).
#
#   Part C. Bootstrap CIs (B = 500, polygon-level resampling within
#           held-out country) on the WINNING configuration per outcome.
#           Honest reporting for the manuscript.
#
#   Part D. Per-outcome univariate-LOCO ranking: which SoilGrids features
#           are consistently top across iron / vitA / B12 / folate / zinc?
#           Identifies the "robust soil feature set" to promote.
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

# Pull the all-source ranking from script 13.
all_uni <- read.csv(file.path(SANDBOX_RESULTS, "13_univariate_loco.csv"))
all_uni <- all_uni[order(-all_uni$mean_r), ]
soil_uni <- all_uni[grepl("^gee_soil", all_uni$variable), ]
cat("all-source top 8:\n");  print(head(all_uni[, c("variable","mean_r")], 8), row.names = FALSE)
cat("\nsoil-only top 8:\n"); print(head(soil_uni[, c("variable","mean_r")], 8), row.names = FALSE)

# Predict helper: spatial GAM + linear extras.
make_spatial_plus <- function(extras) function(train, test, vars) {
  o <- fit_predict_spatial_coords(train, test, vars, extra_covars = extras)
  if (is.null(o)) NULL else o$pred
}

# ── Part A + B: K sweeps ────────────────────────────────────────────────────
results <- list()
for (K in 3:8) {
  for (kind in c("all", "soil")) {
    pool <- if (kind == "all") all_uni$variable else soil_uni$variable
    extras <- head(pool, K)
    name <- sprintf("spatial_%s_top%d", kind, K)
    cat(sprintf("\n[K-sweep] %s ... (vars: %s)\n", name,
                paste(substr(extras, 1, 30), collapse=", ")))
    r <- loco_evaluate_all(pooled_all, make_spatial_plus(extras),
                            method_name = name)
    results[[name]] <- r
  }
}
sweep <- dplyr::bind_rows(results)

# Summary
summ <- sweep |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r   = round(mean(pearson_r,  na.rm = TRUE), 3),
                    median_r = round(median(pearson_r,na.rm = TRUE), 3),
                    wins     = sum(pearson_r > 0.1, na.rm = TRUE),
                    n_splits = sum(!is.na(pearson_r)),
                    mean_mae = round(mean(mae_pp,     na.rm = TRUE), 2),
                    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_r))
cat("\n===== K sweep summary =====\n")
print(as.data.frame(summ), row.names = FALSE)
readr::write_csv(sweep, file.path(SANDBOX_RESULTS, "15a_spatial_topk_sweep.csv"))

best_config <- summ$method[1]
cat(sprintf("\n>>> best config: %s (mean r = %.3f)\n", best_config, summ$mean_r[1]))

# Resolve best K + kind for bootstrap.
best_parts <- strsplit(best_config, "_")[[1]]
best_kind  <- best_parts[2]               # all / soil
best_K     <- as.integer(sub("top","", best_parts[3]))
best_extras <- head(if (best_kind == "all") all_uni$variable else soil_uni$variable,
                     best_K)
cat(sprintf("best_extras (%d):\n  - %s\n", length(best_extras),
            paste(best_extras, collapse = "\n  - ")))

# ── Part C: bootstrap CIs on best config ────────────────────────────────────
B <- 500L
boot_rows <- list()
for (oc in OUTCOMES) {
  cat(sprintf("\n[boot] %s ...\n", oc))
  p <- pooled_all[[oc]]
  for (held_out in p$country_names) {
    train <- p$pooled_data[p$pooled_data$country != held_out, , drop = FALSE]
    test  <- p$pooled_data[p$pooled_data$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    out <- fit_predict_spatial_coords(train, test, p$common_gee_vars,
                                         extra_covars = best_extras)
    if (is.null(out)) next
    obs <- test$svy_prev; pred <- out$pred
    ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
    if (length(obs) < 3) next
    r_point <- suppressWarnings(cor(obs, pred))
    set.seed(2026)
    rs <- numeric(B); n <- length(obs)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      rs[b] <- suppressWarnings(cor(obs[idx], pred[idx]))
    }
    rs <- rs[is.finite(rs)]
    boot_rows[[paste(oc, held_out)]] <- data.frame(
      outcome = oc, held_out = held_out, n_test = n,
      r_point = round(r_point, 3),
      r_lo_80 = round(stats::quantile(rs, 0.1, na.rm = TRUE), 3),
      r_hi_80 = round(stats::quantile(rs, 0.9, na.rm = TRUE), 3),
      r_lo_95 = round(stats::quantile(rs, 0.025, na.rm = TRUE), 3),
      r_hi_95 = round(stats::quantile(rs, 0.975, na.rm = TRUE), 3),
      ci_excludes_zero = stats::quantile(rs, 0.025, na.rm = TRUE) > 0,
      stringsAsFactors = FALSE)
  }
}
boot <- dplyr::bind_rows(boot_rows)
cat("\n===== bootstrap CIs on best config =====\n")
print(boot, row.names = FALSE)
readr::write_csv(boot, file.path(SANDBOX_RESULTS, "15b_bootstrap_best.csv"))
cat(sprintf("\nholdouts with 95%% CI strictly > 0: %d / %d\n",
            sum(boot$ci_excludes_zero, na.rm = TRUE), nrow(boot)))

# ── Part D: per-outcome univariate stability ────────────────────────────────
# Compute univariate LOCO mean r per OUTCOME for each variable, then ask
# which SoilGrids variables are in the top-10 for *all four* outcomes.
cat("\n===== per-outcome univariate stability =====\n")
uni_predict <- function(train, test, vars) {
  imp <- impute_median(train, test, vars); Y <- train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  fit <- tryCatch(stats::lm(Y ~ x, data = data.frame(Y=Y, x=imp$train[[vars]]),
                              weights = wt), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  as.numeric(stats::predict(fit, newdata = data.frame(x = imp$test[[vars]])))
}
soil_vars <- soil_uni$variable
per_oc <- list()
for (oc in OUTCOMES) {
  p <- pooled_all[[oc]]
  v_means <- vapply(soil_vars, function(v) {
    r <- loco_evaluate(list(pooled_data = p$pooled_data,
                             common_gee_vars = v,
                             country_names = p$country_names),
                        predict_fn = function(train, test, vars) uni_predict(train, test, v),
                        outcome = oc, method_name = v)
    if (nrow(r) == 0) NA_real_ else mean(r$pearson_r, na.rm = TRUE)
  }, numeric(1))
  per_oc[[oc]] <- data.frame(variable = soil_vars, mean_r = round(v_means, 3),
                              stringsAsFactors = FALSE)
  per_oc[[oc]]$rank <- rank(-per_oc[[oc]]$mean_r, na.last = "keep", ties.method = "min")
}
# Wide table: variable × outcome's rank
ranks_wide <- per_oc[[1]][, c("variable")]
for (oc in OUTCOMES) ranks_wide[[paste0(oc,"_rank")]] <- per_oc[[oc]]$rank
ranks_wide$max_rank <- apply(ranks_wide[, -1], 1, max)
ranks_wide <- ranks_wide[order(ranks_wide$max_rank), ]
cat("\nSoilGrids variables ranked by WORST-outcome rank (smaller = more stable):\n")
print(head(ranks_wide, 12), row.names = FALSE)
readr::write_csv(ranks_wide, file.path(SANDBOX_RESULTS, "15c_soil_stability.csv"))

cat("\nDONE — see sandbox/results/15{a,b,c}_*.csv\n")
