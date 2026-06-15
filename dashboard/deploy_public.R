# =============================================================================
# deploy_public.R — deploy the lean, policymaker-facing public app
#
# Deploys app_public.R as a SEPARATE shinyapps.io app (distinct URL) from the
# full internal app. Run AFTER refreshing dashboard data:
#   Rscript dashboard/deploy_public.R
# =============================================================================
if (!requireNamespace("rsconnect", quietly = TRUE))
  stop("Install 'rsconnect': install.packages('rsconnect')")

APP_NAME  <- "micronutrient-burden-public"
APP_TITLE <- "Micronutrient Burden (public)"

dashboard_dir <- here::here("dashboard")
data_dir <- file.path(dashboard_dir, "data")

required_files <- c("admin2_predictions.rds", "admin2_population.rds",
                    "admin2_boundaries.rds", "admin1_boundaries.rds",
                    "national_estimates.rds", "cv_performance.rds",
                    "metadata.rds")
missing <- setdiff(required_files, list.files(data_dir))
if (length(missing) > 0)
  stop("Missing data files: ", paste(missing, collapse = ", "),
       "\nRun dashboard/data-raw/01_prepare_dashboard_data.R first.")

accts <- rsconnect::accounts()
if (nrow(accts) == 0)
  stop("No shinyapps.io account configured (rsconnect::setAccountInfo).")
cat(sprintf("Deploying public app as: %s\n", accts$name[1]))

# Bundle the public entry point + shared module/data files (exclude app.R so
# shinyapps does not treat the full app as the entry point).
deploy_files <- c(
  "app_public.R", "global.R",
  file.path("R",    list.files(file.path(dashboard_dir, "R"))),
  file.path("data", list.files(file.path(dashboard_dir, "data")))
)
if (dir.exists(file.path(dashboard_dir, "www")))
  deploy_files <- c(deploy_files,
                    file.path("www", list.files(file.path(dashboard_dir, "www"))))

exists_check <- file.exists(file.path(dashboard_dir, deploy_files))
if (!all(exists_check))
  stop("Missing appFiles: ", paste(deploy_files[!exists_check], collapse = ", "))

cat(sprintf("\nDeploying %d files as '%s'...\n", length(deploy_files), APP_NAME))
rsconnect::deployApp(
  appDir         = dashboard_dir,
  appFiles       = deploy_files,
  appPrimaryDoc  = "app_public.R",
  appName        = APP_NAME,
  appTitle       = APP_TITLE,
  account        = "amertens",
  server         = "shinyapps.io",
  forceUpdate    = TRUE,
  launch.browser = FALSE
)
cat(sprintf("\nDone. Visit: https://%s.shinyapps.io/%s/\n", accts$name[1], APP_NAME))
