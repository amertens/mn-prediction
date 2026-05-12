# Validate deploy.R prerequisites without actually deploying
dashboard_dir <- here::here("dashboard")
data_dir <- file.path(dashboard_dir, "data")

cat("Checking required data files...\n")
required_files <- c(
  "admin2_predictions.rds", "admin2_population.rds",
  "admin2_boundaries.rds", "admin1_boundaries.rds",
  "national_estimates.rds", "cv_performance.rds",
  "metadata.rds", "gbd_estimates.rds"
)
existing <- list.files(data_dir)
missing <- setdiff(required_files, existing)
if (length(missing) == 0) {
  cat("  ✓ All required data files present\n")
} else {
  cat("  ✗ Missing:", missing, "\n")
}

cat("\nChecking deployable files...\n")
deploy_files <- c(
  file.path(dashboard_dir, c("app.R", "global.R")),
  list.files(file.path(dashboard_dir, "R"), full.names = TRUE),
  list.files(file.path(dashboard_dir, "data"), full.names = TRUE)
)
total_size <- sum(file.size(deploy_files[file.exists(deploy_files)]))
cat(sprintf("  Total bundle size: %.1f MB across %d files\n",
            total_size / 1024 / 1024,
            sum(file.exists(deploy_files))))

cat("\nChecking shinyapps.io account...\n")
accts <- rsconnect::accounts()
print(accts)

cat("\nReady to deploy. Run:\n")
cat("  Rscript dashboard/deploy.R\n")
