# =============================================================================
# sandbox/22_all_methods_comparison.R
#
# Sensitivity / supplementary analysis: comparison of every method and
# SuperLearner library tested in this work, applied to ONE outcome
# (child_iron — strongest signal, all 4 countries with data).
#
# Two evaluation tracks per method:
#   - WITHIN-COUNTRY  : 5-fold CV on each country's admin-2 polygons
#                        (separately per country, averaged).
#   - LEAVE-ONE-COUNTRY-OUT : 4 holdouts, train on 3, predict on 1.
#
# Outputs:
#   sandbox/results/22_all_methods_child_iron.csv
#   sandbox/results/22_all_methods_child_iron.png  (2-panel figure)
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here); library(sf)
  library(ggplot2); library(patchwork)
})
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))
source(here::here("sandbox", "00_setup.R"))

OUTCOME <- "child_iron"
cat(sprintf("===== outcome: %s =====\n\n", OUTCOME))

# ── Load pooled area-level data ─────────────────────────────────────────
pooled <- load_pooled(OUTCOME)
cc_list <- get_country_configs()
pooled$pooled_data <- add_admin2_centroids(pooled$pooled_data, cc_list,
                                              pooled$svy_list)
pooled$pooled_data$svy_prev <- as.numeric(pooled$pooled_data$svy_prev)
cat(sprintf("pooled: %d admin-2 polygons across %d countries\n",
            nrow(pooled$pooled_data), length(pooled$country_names)))

# ── Cross-outcome pooled (for combined_filter) ──────────────────────────
cross_pool <- list()
for (oc in c("child_vitA","women_vitA","child_iron","women_iron")) {
  p <- tryCatch(load_pooled(oc), error = function(e) NULL)
  if (!is.null(p)) {
    p$pooled_data <- add_admin2_centroids(p$pooled_data, cc_list, p$svy_list)
    cross_pool[[oc]] <- p$pooled_data
  }
}

# ── Define methods to evaluate ──────────────────────────────────────────
wrap_loco <- function(fn, ...) function(train, test, vars) {
  o <- fn(train, test, vars, ...); if (is.null(o)) NULL else o$pred
}

METHODS <- list(
  baseline             = wrap_loco(fit_predict_baseline,    model_type = "continuous"),
  glm                  = wrap_loco(fit_predict_glm,         model_type = "continuous"),
  penalized            = wrap_loco(fit_predict_penalized,   model_type = "continuous"),
  fh                   = wrap_loco(fit_predict_fh),
  mixed                = wrap_loco(fit_predict_mixed,       model_type = "continuous"),
  two_stage            = wrap_loco(fit_predict_two_stage),
  forward              = wrap_loco(fit_predict_forward_loco, model_type = "continuous"),
  gam                  = wrap_loco(fit_predict_gam,         model_type = "continuous"),
  grouplasso           = wrap_loco(fit_predict_grouplasso),
  quasibinomial        = wrap_loco(fit_predict_quasibinomial),
  betareg              = wrap_loco(fit_predict_betareg),
  dag                  = wrap_loco(fit_predict_dag,         outcome_tag = OUTCOME,
                                     model_type = "continuous"),
  invariance_filter    = wrap_loco(fit_predict_invariance_filter,
                                     model_type = "continuous"),
  combined_filter      = function(train, test, vars)
                          {o <- fit_predict_combined_filter(train, test, vars,
                              cross_outcome_pooled = cross_pool,
                              this_outcome = OUTCOME);
                           if (is.null(o)) NULL else o$pred},
  spatial_coords       = wrap_loco(fit_predict_spatial_coords),
  spatial_plus_soil    = wrap_loco(fit_predict_spatial_plus_soil),
  sl_prescreened       = wrap_loco(fit_predict_sl_prescreened,
                                      K_max = 40L, min_within_ratio = 0.30,
                                      min_uni_abs_r = 0.05),
  # Standalone fast learners on the (lon, lat + top-8 soil) feature set
  spatial_ranger       = function(train, test, vars) {
    if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
    soil <- intersect(.SPATIAL_PLUS_SOIL_DEFAULT_VARS, colnames(train))
    if (!all(c("lon","lat") %in% colnames(train))) return(NULL)
    Xtr <- as.data.frame(train[, c("lon","lat", soil), drop = FALSE])
    Xte <- as.data.frame(test [, c("lon","lat", soil), drop = FALSE])
    med <- vapply(Xtr, median, numeric(1), na.rm = TRUE)
    for (v in names(Xtr)) {
      Xtr[[v]][!is.finite(Xtr[[v]])] <- med[v]
      Xte[[v]][!is.finite(Xte[[v]])] <- med[v] }
    fit <- tryCatch(ranger::ranger(Y ~ ., data = data.frame(Y = train$svy_prev, Xtr),
                                      num.trees = 500, min.node.size = 5,
                                      case.weights = pmax(train$n_svy, 1)),
                     error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pmin(pmax(predict(fit, data = Xte)$predictions, 0), 1)
  },
  spatial_xgb         = function(train, test, vars) {
    if (!requireNamespace("xgboost", quietly = TRUE)) return(NULL)
    soil <- intersect(.SPATIAL_PLUS_SOIL_DEFAULT_VARS, colnames(train))
    if (!all(c("lon","lat") %in% colnames(train))) return(NULL)
    Xtr <- as.matrix(as.data.frame(train[, c("lon","lat", soil), drop = FALSE]))
    Xte <- as.matrix(as.data.frame(test [, c("lon","lat", soil), drop = FALSE]))
    Xtr[!is.finite(Xtr)] <- 0; Xte[!is.finite(Xte)] <- 0
    dtr <- xgboost::xgb.DMatrix(Xtr, label = train$svy_prev,
                                  weight = pmax(train$n_svy, 1))
    fit <- tryCatch(xgboost::xgb.train(dtr, nrounds = 300, max_depth = 4,
                                          eta = 0.05, subsample = 0.8,
                                          colsample_bytree = 0.8,
                                          objective = "reg:squarederror",
                                          verbose = 0),
                     error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pmin(pmax(predict(fit, Xte), 0), 1)
  }
)

# ── LOCO across all methods ─────────────────────────────────────────────
cat(sprintf("\n===== running LOCO for %d methods =====\n", length(METHODS)))
pd <- pooled$pooled_data
loco_rows <- list()
for (m in names(METHODS)) {
  cat(sprintf("[loco %-25s] ", m))
  t0 <- Sys.time()
  per_holdout <- list()
  for (held in pooled$country_names) {
    train <- pd[pd$country != held, , drop = FALSE]
    test  <- pd[pd$country == held, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    p <- tryCatch(METHODS[[m]](train, test, pooled$common_gee_vars),
                  error = function(e) NULL)
    if (is.null(p) || length(p) != nrow(test)) next
    ok <- is.finite(p) & is.finite(test$svy_prev)
    if (sum(ok) < 3) next
    per_holdout[[held]] <- data.frame(
      track = "LOCO", method = m, country = held,
      n = sum(ok),
      pearson_r  = round(cor(test$svy_prev[ok], p[ok]), 3),
      spearman_r = round(cor(test$svy_prev[ok], p[ok], method = "spearman"), 3),
      mae_pp = round(mean(abs(test$svy_prev[ok] - p[ok])) * 100, 2))
  }
  cat(sprintf("%.1fs\n", as.numeric(Sys.time()-t0, units = "secs")))
  loco_rows[[m]] <- dplyr::bind_rows(per_holdout)
}
loco_results <- dplyr::bind_rows(loco_rows)

# ── Within-country 5-fold CV (subset of methods that work at small n) ──
within_methods <- c("baseline","glm","penalized","gam","spatial_coords",
                    "spatial_plus_soil","spatial_ranger","spatial_xgb",
                    "sl_prescreened")
cat(sprintf("\n===== running within-country 5-fold CV for %d methods =====\n",
            length(within_methods)))
within_rows <- list()
for (ctry in pooled$country_names) {
  d <- pd[pd$country == ctry, , drop = FALSE]
  if (nrow(d) < 8) next
  set.seed(42)
  K <- min(5, nrow(d) - 1)
  folds <- sample(rep(seq_len(K), length.out = nrow(d)))
  for (m in within_methods) {
    t0 <- Sys.time()
    preds <- rep(NA_real_, nrow(d))
    for (k in seq_len(K)) {
      tr <- d[folds != k, , drop = FALSE]; te <- d[folds == k, , drop = FALSE]
      if (nrow(tr) < 3) next
      p <- tryCatch(METHODS[[m]](tr, te, pooled$common_gee_vars),
                    error = function(e) NULL)
      if (!is.null(p) && length(p) == nrow(te)) preds[folds == k] <- p
    }
    ok <- is.finite(preds) & is.finite(d$svy_prev)
    if (sum(ok) < 3) next
    within_rows[[paste(ctry, m)]] <- data.frame(
      track = "within", method = m, country = ctry,
      n = sum(ok),
      pearson_r  = round(cor(d$svy_prev[ok], preds[ok]), 3),
      spearman_r = round(cor(d$svy_prev[ok], preds[ok], method = "spearman"), 3),
      mae_pp = round(mean(abs(d$svy_prev[ok] - preds[ok])) * 100, 2))
    cat(sprintf("[within %-12s / %-25s] r = %+.2f (%.1fs)\n", ctry, m,
                cor(d$svy_prev[ok], preds[ok]),
                as.numeric(Sys.time()-t0, units="secs")))
  }
}
within_results <- dplyr::bind_rows(within_rows)

all <- dplyr::bind_rows(within_results, loco_results)
csv_path <- file.path(SANDBOX_RESULTS,
                        paste0("22_all_methods_", OUTCOME, ".csv"))
readr::write_csv(all, csv_path)

# ── Summary table ───────────────────────────────────────────────────────
cat("\n===== summary (mean Pearson r per method) =====\n")
summary_tbl <- all |>
  dplyr::group_by(track, method) |>
  dplyr::summarise(
    mean_r   = round(mean(pearson_r,  na.rm = TRUE), 3),
    median_r = round(median(pearson_r,na.rm = TRUE), 3),
    n_folds  = dplyr::n(),
    mean_mae = round(mean(mae_pp,     na.rm = TRUE), 2),
    .groups = "drop") |>
  dplyr::arrange(track, dplyr::desc(mean_r))
print(as.data.frame(summary_tbl), row.names = FALSE)

# ── Build the figure ────────────────────────────────────────────────────
cat("\n===== building figure =====\n")
plot_df <- all
# Method ordering: by mean LOCO r across methods (so the figure is sortable)
method_order <- loco_results |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r = mean(pearson_r, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(mean_r)
plot_df$method <- factor(plot_df$method, levels = method_order$method)

p_loco <- plot_df |> dplyr::filter(track == "LOCO") |>
  ggplot(aes(x = method, y = pearson_r, colour = country)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_summary(aes(group = method), fun = mean, geom = "point",
                shape = 4, size = 4, colour = "black", stroke = 1.2) +
  coord_flip() +
  labs(title = sprintf("(b) LOCO transportability — %s", OUTCOME),
        subtitle = "Dot = per-holdout Pearson r; X = mean across 4 holdouts",
        x = NULL, y = "Pearson r (predicted vs observed admin-2 prevalence)",
        colour = "Held-out") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_within <- plot_df |> dplyr::filter(track == "within") |>
  ggplot(aes(x = method, y = pearson_r, colour = country)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_summary(aes(group = method), fun = mean, geom = "point",
                shape = 4, size = 4, colour = "black", stroke = 1.2) +
  coord_flip() +
  labs(title = sprintf("(a) Within-country 5-fold CV — %s", OUTCOME),
        subtitle = "Dot = per-country Pearson r; X = mean across 4 countries",
        x = NULL, y = "Pearson r",
        colour = "Country") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# Same method levels on both panels for visual alignment.
p_within$data$method <- factor(p_within$data$method,
                                 levels = method_order$method)

fig <- p_within + p_loco + plot_layout(ncol = 2)
fig_path <- file.path(SANDBOX_RESULTS,
                        paste0("22_all_methods_", OUTCOME, ".png"))
ggsave(fig_path, fig, width = 13, height = 7, dpi = 150)
cat(sprintf("written:\n  %s\n  %s\n", csv_path, fig_path))
