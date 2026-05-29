# =============================================================================
# scripts/run_sl_prescreened_main.R
#
# Standalone driver for fit_predict_sl_prescreened() (R/benchmark_models.R
# METHOD 20). Reads pooled admin-2 data directly from the _targets_full
# store, runs the prescreened SuperLearner LOCO across each outcome where
# at least three countries have data, writes a tidy CSV that the dashboard
# data-prep script can consume.
#
# This script does the same thing as the `sl_prescreened_<outcome>` +
# `sl_prescreened_all` targets in _targets.R but doesn't require
# tar_make() (so it doesn't trigger the package-version cascade we hit
# earlier). Useful when the user wants results NOW.
#
# Usage:
#   source(here::here("scripts", "run_sl_prescreened_main.R"))
#
# Output:
#   results/tables/sl_prescreened_main.csv
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here); library(sf)
})

source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))

STORE <- here::here("_targets_full")

COUNTRIES_LOWER <- c("gambia","ghana","sierraleone","malawi")
COUNTRY_LABELS  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")

# Outcomes with their country availability.
OUTCOME_COUNTRIES <- list(
  child_vitA   = COUNTRIES_LOWER,
  women_vitA   = COUNTRIES_LOWER,
  child_iron   = COUNTRIES_LOWER,
  women_iron   = COUNTRIES_LOWER,
  women_b12    = setdiff(COUNTRIES_LOWER, "gambia"),
  women_folate = setdiff(COUNTRIES_LOWER, "gambia")
)

cc_list <- get_country_configs()

t_all <- Sys.time()
all_results <- list()

for (otag in names(OUTCOME_COUNTRIES)) {
  cc <- OUTCOME_COUNTRIES[[otag]]
  if (length(cc) < 2) next
  cat(sprintf("\n===== %s (countries: %s) =====\n", otag,
              paste(cc, collapse = ", ")))
  labels <- names(COUNTRY_LABELS)[match(cc, COUNTRY_LABELS)]

  # Load svy + gee.
  svy_list <- setNames(lapply(cc, function(c)
    tryCatch(targets::tar_read_raw(paste0("svy_admin2_", c, "_", otag),
                                     store = STORE),
              error = function(e) NULL)), labels)
  gee_list <- setNames(lapply(cc, function(c)
    tryCatch(targets::tar_read_raw(paste0("gee_admin2_", c), store = STORE),
              error = function(e) NULL)), labels)
  ok <- !vapply(svy_list, is.null, logical(1)) &
        !vapply(gee_list, is.null, logical(1))
  if (sum(ok) < 2) {
    cat("  skip — too few countries with both svy + gee\n"); next
  }
  svy_list <- svy_list[ok]; gee_list <- gee_list[ok]

  pooled <- build_area_loco_dataset(svy_list, gee_list)
  if (is.null(pooled) || nrow(pooled$pooled_data) < 10) {
    cat("  skip — pooled data too small\n"); next
  }
  pooled$pooled_data <- add_admin2_centroids(pooled$pooled_data,
                                                cc_list, svy_list)

  t_oc <- Sys.time()
  out <- run_area_benchmarks_loco(
    pooled_data    = pooled$pooled_data,
    gee_vars       = pooled$common_gee_vars,
    country_names  = pooled$country_names,
    adjacency_list = NULL,
    outcome_label  = otag,
    methods        = "sl_prescreened",
    augment_features = FALSE)
  cat(sprintf("  elapsed: %.0fs\n",
              as.numeric(Sys.time() - t_oc, units = "secs")))
  if (nrow(out) > 0) {
    out$eval_type <- "loco_sl_only"
    out$outcome   <- otag
    cat(sprintf("  mean LOCO Pearson r = %+.3f (%d / %d holdouts)\n",
                mean(out$pearson_r, na.rm = TRUE),
                sum(out$pearson_r > 0.1, na.rm = TRUE),
                sum(!is.na(out$pearson_r))))
    all_results[[otag]] <- out
  }
}

if (length(all_results) > 0) {
  result <- dplyr::bind_rows(all_results)
  out_dir <- here::here("results", "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  csv_path <- file.path(out_dir, "sl_prescreened_main.csv")
  readr::write_csv(result, csv_path)
  cat(sprintf("\n===== summary =====\n"))
  print(result |> dplyr::group_by(outcome) |>
          dplyr::summarise(mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
                            wins   = sum(pearson_r > 0.1, na.rm = TRUE),
                            n      = dplyr::n(), .groups = "drop"))
  cat(sprintf("\nwritten: %s\n", csv_path))
} else {
  warning("No SL prescreened results produced")
}

cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(Sys.time() - t_all, units = "mins")))
