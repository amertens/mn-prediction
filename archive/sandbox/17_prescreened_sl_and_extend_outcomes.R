# =============================================================================
# sandbox/17_prescreened_sl_and_extend_outcomes.R
#
# Two follow-ups:
#
# Part A. Test the new fit_predict_sl_prescreened() SuperLearner (METHOD 20)
#         against the leaders (spatial_plus_soil, spatial_coords, H6 combined,
#         baseline penalised) on the original 4 LOCO outcomes. Does the
#         prescreening repair the 17-hour SL's failures?
#
# Part B. Extend the leaders to outcomes 5 & 6 (women_b12 and women_folate),
#         which have data in 3 of the 4 countries (Ghana, Sierra Leone,
#         Malawi). 3-country LOCO = train on 2, test on 1.
#         Outcomes 7 & 8 (child_zinc, women_zinc) are Malawi-only, so we
#         report a within-country 5-fold CV mean r instead.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "config.R"))
cc_list <- get_country_configs()

# Override OUTCOMES + extend.
ALL_OUTCOMES <- c("child_vitA", "women_vitA", "child_iron", "women_iron",
                   "women_b12", "women_folate", "child_zinc", "women_zinc")
LOCO_OUTCOMES_3PLUS <- c(ALL_OUTCOMES[1:4],  # 4-country (existing)
                         "women_b12", "women_folate")  # 3-country
WITHIN_ONLY <- c("child_zinc", "women_zinc")  # Malawi only

# Load pooled per outcome. For 3-country outcomes load_pooled() will return
# pooled data missing Gambia — that's fine, build_area_loco_dataset handles it.
pooled_all <- list()
for (oc in LOCO_OUTCOMES_3PLUS) {
  cat(sprintf("[load] %s ...\n", oc))
  pooled_all[[oc]] <- load_pooled(oc)
}
# Augment all with admin-2 centroids.
example_svy_list <- pooled_all[[1]]$svy_list
for (oc in names(pooled_all))
  pooled_all[[oc]]$pooled_data <-
    add_admin2_centroids(pooled_all[[oc]]$pooled_data, cc_list, example_svy_list)

# Cross-outcome pooled list — required by combined_filter.
cross_pooled <- setNames(lapply(pooled_all, function(x) x$pooled_data),
                          names(pooled_all))

# ── Predict wrappers ────────────────────────────────────────────────────────
wrap <- function(fn, ...) function(train, test, vars) {
  o <- fn(train, test, vars, ...); if (is.null(o)) NULL else o$pred
}

methods <- list(
  spatial_plus_soil = wrap(fit_predict_spatial_plus_soil),
  sl_prescreened    = wrap(fit_predict_sl_prescreened, K_max = 40L,
                              min_within_ratio = 0.30, min_uni_abs_r = 0.05),
  penalized_baseline = wrap(fit_predict_penalized, model_type = "continuous")
)

# Run LOCO on each outcome where it's possible (3-country needs >= 3 country
# names in the country_names entry of pooled).
res_all <- list()
for (oc in names(pooled_all)) {
  cat(sprintf("\n===== outcome: %s =====\n", oc))
  p <- pooled_all[[oc]]
  cat(sprintf("  countries: %s; rows: %d\n",
              paste(p$country_names, collapse=","), nrow(p$pooled_data)))
  if (length(p$country_names) < 3) {
    cat("  skipping LOCO (need >= 3 countries)\n"); next
  }
  for (m in names(methods)) {
    t_m <- Sys.time()
    cat(sprintf("  [%s] ...", m))
    r <- loco_evaluate(p, methods[[m]], outcome = oc, method_name = m)
    cat(sprintf(" mean_r=%+.3f  (%.0fs)\n",
                if (nrow(r) > 0) mean(r$pearson_r, na.rm = TRUE) else NA,
                as.numeric(Sys.time() - t_m, units = "secs")))
    flush.console()
    if (nrow(r) > 0) res_all[[paste(oc, m)]] <- r
  }
}
res_df <- dplyr::bind_rows(res_all)

# Summary across outcomes
summ <- res_df |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
                    median_r = round(median(pearson_r, na.rm = TRUE), 3),
                    wins = sum(pearson_r > 0.1, na.rm = TRUE),
                    n_splits = sum(!is.na(pearson_r)),
                    mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_r))

cat("\n===== summary across 6 LOCO-able outcomes (4 + women_b12 + women_folate) =====\n")
print(as.data.frame(summ), row.names = FALSE)

# Per-outcome × method table
per_oc <- res_df |>
  dplyr::group_by(method, outcome) |>
  dplyr::summarise(mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
                    wins = sum(pearson_r > 0.1, na.rm = TRUE),
                    n_splits = sum(!is.na(pearson_r)),
                    .groups = "drop") |>
  tidyr::pivot_wider(names_from = outcome, values_from = c(mean_r, wins))
cat("\n===== per outcome mean r =====\n")
print(as.data.frame(per_oc), row.names = FALSE)

readr::write_csv(res_df,
                  file.path(SANDBOX_RESULTS,
                              "17_prescreened_sl_and_extend.csv"))
cat(sprintf("\nwritten: %s\n",
            file.path(SANDBOX_RESULTS,
                        "17_prescreened_sl_and_extend.csv")))
