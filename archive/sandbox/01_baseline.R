# =============================================================================
# sandbox/01_baseline.R
#
# Reproduce the current best of the pipeline's parsimonious LOCO comparators
# (penalized GLM, mixed-effects, two-stage, GAM, quasi-binomial) using the
# sandbox LOCO evaluator. Establishes the bar that every hypothesis script
# below has to beat.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

pooled_all <- load_all_pooled()

# Wrap each existing fit_predict_* into a predict_fn(train, test, vars) closure
# returning a numeric prediction for nrow(test).
wrap <- function(fn, ...) function(train, test, vars) {
  out <- fn(train, test, vars, ...); if (is.null(out)) NULL else out$pred
}

methods <- list(
  penalized      = wrap(fit_predict_penalized,    model_type = "continuous"),
  mixed          = wrap(fit_predict_mixed,        model_type = "continuous"),
  two_stage      = wrap(fit_predict_two_stage),
  gam            = wrap(fit_predict_gam,          model_type = "continuous"),
  forward        = wrap(fit_predict_forward_loco, model_type = "continuous"),
  quasibinomial  = wrap(fit_predict_quasibinomial)
)

all_res <- list()
for (m in names(methods)) {
  cat(sprintf("[baseline] %s ...\n", m))
  r <- loco_evaluate_all(pooled_all, methods[[m]], method_name = m)
  all_res[[m]] <- r
}
baseline <- dplyr::bind_rows(all_res)

# Per-outcome best-method tally
summary <- baseline |>
  dplyr::group_by(method) |>
  dplyr::summarise(mean_r   = round(mean(pearson_r,  na.rm = TRUE), 3),
                    median_r = round(median(pearson_r,na.rm = TRUE), 3),
                    mean_mae = round(mean(mae_pp,     na.rm = TRUE), 2),
                    n        = sum(!is.na(pearson_r))) |>
  dplyr::arrange(dplyr::desc(mean_r))

cat("\n===== baseline LOCO summary (continuous, mean Pearson r) =====\n")
print(as.data.frame(summary), row.names = FALSE)

readr::write_csv(baseline, file.path(SANDBOX_RESULTS, "01_baseline.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "01_baseline.csv")))
