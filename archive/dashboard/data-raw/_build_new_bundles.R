# Lightweight builder for the two new dashboard bundles (no GADM/geodata).
# Mirrors sections 7d & 7e of 01_prepare_dashboard_data.R so the new tabs have
# data without re-running the full (heavy) prep.
`%||%` <- function(a, b) if (is.null(a)) b else a
root <- getwd()
if (!dir.exists(file.path(root, "results"))) {
  root <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
}
tdir <- file.path(root, "results", "tables")
ddir <- file.path(root, "dashboard", "data")
dir.create(ddir, showWarnings = FALSE, recursive = TRUE)

read_csv_if <- function(rel) {
  p <- file.path(tdir, rel)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}

model_diagnostics <- list(
  binary      = read_csv_if("diagnostics_binary.csv"),
  continuous  = read_csv_if("diagnostics_continuous.csv"),
  calibrated  = read_csv_if("diagnostics_binary_calibrated.csv"),
  roc         = read_csv_if("roc_curves.csv"),
  pr          = read_csv_if("pr_curves.csv"),
  calibration = read_csv_if("calibration_tables.csv"),
  build_time  = Sys.time()
)
saveRDS(model_diagnostics, file.path(ddir, "model_diagnostics.rds"))

benchmarks <- list(
  benchmarks      = read_csv_if("benchmarks_all.csv"),
  area_comparison = read_csv_if("area_comparison_all.csv"),
  admin2_error    = read_csv_if("admin2_error_all.csv"),
  sl_prescreened  = read_csv_if("sl_prescreened_main.csv"),
  area_transport  = read_csv_if("transportability_area_loco_metrics.csv"),
  build_time      = Sys.time()
)
saveRDS(benchmarks, file.path(ddir, "benchmarks.rds"))

cat(sprintf("model_diagnostics: binary=%d continuous=%d calibrated=%d roc=%d pr=%d calib=%d\n",
            nrow(model_diagnostics$binary %||% data.frame()),
            nrow(model_diagnostics$continuous %||% data.frame()),
            nrow(model_diagnostics$calibrated %||% data.frame()),
            nrow(model_diagnostics$roc %||% data.frame()),
            nrow(model_diagnostics$pr %||% data.frame()),
            nrow(model_diagnostics$calibration %||% data.frame())))
cat(sprintf("benchmarks: benchmarks=%d area_comparison=%d admin2_error=%d sl_prescreened=%d area_transport=%d\n",
            nrow(benchmarks$benchmarks %||% data.frame()),
            nrow(benchmarks$area_comparison %||% data.frame()),
            nrow(benchmarks$admin2_error %||% data.frame()),
            nrow(benchmarks$sl_prescreened %||% data.frame()),
            nrow(benchmarks$area_transport %||% data.frame())))
