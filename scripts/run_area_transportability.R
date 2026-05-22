# =============================================================================
# scripts/run_area_transportability.R
#
# Generate the final area-level transportability deliverables using the
# selected universal/parsimonious recipe (R/transportability_area.R):
#   * results/tables/transportability_area_loco_metrics.csv
#   * results/tables/transportability_area_loco_predictions.csv  (Admin-2 level)
#   * results/transportability/area_loco_predictions.rds         (for dashboard)
#   * results/tables/transportability_area_selected_vars.csv     (interpretability)
#
# Reads inputs from the cached `_targets` store.
# =============================================================================
suppressMessages({ library(targets); library(dplyr) })
source("R/transportability_area.R")
set.seed(20260521)

STORE <- "_targets"
COUNTRIES <- c("gambia", "ghana", "malawi", "sierraleone")
OUTCOMES  <- c("child_vitA", "child_iron", "women_vitA", "women_iron")

cat("Building enriched Admin-2 covariates...\n")
cov_list <- lapply(COUNTRIES, function(cn) {
  g <- tar_read_raw(paste0("gee_admin2_", cn), store = STORE)
  m <- tar_read_raw(paste0("merged_", cn), store = STORE)
  d <- build_admin2_covariates(g, m)
  cat(sprintf("  %s: %d areas x %d covariates\n", cn, nrow(d), ncol(d) - 1))
  d
})
names(cov_list) <- COUNTRIES

all_metrics <- list(); all_preds <- list(); all_sel <- list()
for (oc in OUTCOMES) {
  svy_list <- lapply(COUNTRIES, function(cn)
    tryCatch(tar_read_raw(paste0("svy_admin2_", cn, "_", oc), store = STORE),
             error = function(e) NULL))
  names(svy_list) <- COUNTRIES
  pooled <- assemble_area_transport(svy_list, cov_list, outcome = oc)
  if (is.null(pooled)) { cat(sprintf("[%s] insufficient data\n", oc)); next }
  cat(sprintf("\n[%s] %d areas, %d common predictors, countries: %s\n",
              oc, nrow(pooled$data), length(pooled$predictors),
              paste(pooled$countries, collapse = ", ")))
  res <- run_area_transport_loco(pooled, AREA_TRANSPORT_RECIPE)
  if (!is.null(res$metrics)) {
    print(res$metrics[, c("held_out","n_test","n_selected","pearson_r",
                          "spearman_r","rmse_pp","mae_pp","nat_bias_pp")])
    all_metrics[[oc]] <- res$metrics
  }
  if (!is.null(res$predictions)) all_preds[[oc]] <- res$predictions
  sel <- unlist(res$selected_vars)
  if (length(sel)) all_sel[[oc]] <- data.frame(outcome = oc, var = sel)
}

metrics_df <- bind_rows(all_metrics)
preds_df   <- bind_rows(all_preds)
sel_df     <- bind_rows(all_sel)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/transportability", showWarnings = FALSE, recursive = TRUE)
write.csv(metrics_df, "results/tables/transportability_area_loco_metrics.csv", row.names = FALSE)
write.csv(preds_df,   "results/tables/transportability_area_loco_predictions.csv", row.names = FALSE)
saveRDS(preds_df,     "results/transportability/area_loco_predictions.rds")

# Variable-selection frequency (how often each predictor is retained) ------
if (nrow(sel_df)) {
  freq <- sel_df |> count(var, name = "n_folds_selected") |>
    arrange(desc(n_folds_selected))
  write.csv(freq, "results/tables/transportability_area_selected_vars.csv", row.names = FALSE)
}

cat("\n\n############ LOCO summary by outcome (selected recipe) ############\n")
options(width = 200)
summ <- metrics_df |> group_by(outcome) |>
  summarise(n_countries = n(),
            mean_pearson = round(mean(pearson_r, na.rm = TRUE), 3),
            mean_spearman = round(mean(spearman_r, na.rm = TRUE), 3),
            mean_rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
            mean_mae_pp  = round(mean(mae_pp, na.rm = TRUE), 2),
            med_n_selected = round(median(n_selected), 0), .groups = "drop")
print(as.data.frame(summ))
cat(sprintf("\nOVERALL (simple country mean) Spearman = %.3f | Pearson = %.3f | median selected vars = %d\n",
            mean(metrics_df$spearman_r, na.rm = TRUE),
            mean(metrics_df$pearson_r, na.rm = TRUE),
            round(median(metrics_df$n_selected), 0)))

# n_test-weighted means (down-weights tiny test countries e.g. Sierra Leone n=14)
wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[is.finite(x)], na.rm = TRUE)
cat(sprintf("OVERALL (n-weighted)        Spearman = %.3f | Pearson = %.3f\n",
            wm(metrics_df$spearman_r, metrics_df$n_test),
            wm(metrics_df$pearson_r,  metrics_df$n_test)))

# Pooled within-country-demeaned correlation across ALL held-out Admin-2 units:
# the single cleanest measure of spatial transfer (level removed per country).
pooled_r <- preds_df |>
  group_by(outcome, country) |>
  mutate(s_c = survey_prev - mean(survey_prev),
         m_c = modeled_prev - mean(modeled_prev)) |>
  ungroup()
cat(sprintf("POOLED within-country-demeaned spatial r (all areas): Pearson = %.3f | Spearman = %.3f\n",
            cor(pooled_r$s_c, pooled_r$m_c, use = "complete.obs"),
            cor(rank(pooled_r$s_c), rank(pooled_r$m_c), use = "complete.obs")))

if (exists("freq")) {
  cat("\nTop 20 most frequently selected predictors:\n")
  print(head(as.data.frame(freq), 20))
}
cat("\nDONE\n")
