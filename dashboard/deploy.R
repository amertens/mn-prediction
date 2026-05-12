# =============================================================================
# deploy.R — Deploy the dashboard to shinyapps.io
#
# Prerequisites (one-time setup):
#
# 1. Sign up for a free shinyapps.io account at https://www.shinyapps.io/
#
# 2. From the shinyapps.io dashboard, get your account name and tokens
#    (Account → Tokens → Show → Copy to clipboard).
#
# 3. Run the rsconnect setup command pasted from the dashboard. It looks like:
#
#      rsconnect::setAccountInfo(name = 'YOUR_ACCOUNT',
#                                token = 'XXXX',
#                                secret = 'YYYY')
#
#    This stores credentials in ~/.config/R/rsconnect/accounts/.
#
# 4. Install required packages locally:
#
#      install.packages(c("rsconnect", "bslib", "bsicons", "leaflet",
#                         "plotly", "reactable", "sf", "shiny",
#                         "dplyr", "tidyr", "htmltools"))
#
# 5. Refresh the dashboard data (so /data/ is up to date):
#
#      Rscript dashboard/data-raw/01_prepare_dashboard_data.R
#      Rscript dashboard/data-raw/02_gbd_placeholder.R
#
# 6. Deploy with:
#
#      Rscript dashboard/deploy.R
# =============================================================================

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install the 'rsconnect' package: install.packages('rsconnect')")
}

# ── Configuration ──────────────────────────────────────────────────────────
APP_NAME <- "micronutrient-burden"  # change if you want a different URL
APP_TITLE <- "Micronutrient Burden Dashboard"

# ── Sanity checks ──────────────────────────────────────────────────────────
dashboard_dir <- here::here("dashboard")
data_dir <- file.path(dashboard_dir, "data")

required_files <- c(
  "admin2_predictions.rds",
  "admin2_population.rds",
  "admin2_boundaries.rds",
  "admin1_boundaries.rds",
  "national_estimates.rds",
  "cv_performance.rds",
  "metadata.rds",
  "gbd_estimates.rds"
)

missing <- setdiff(required_files, list.files(data_dir))
if (length(missing) > 0) {
  cat("\nMissing data files:\n")
  cat("  ", missing, sep = "\n  ")
  cat("\nRun the following first:\n")
  cat("  Rscript dashboard/data-raw/01_prepare_dashboard_data.R\n")
  cat("  Rscript dashboard/data-raw/02_gbd_placeholder.R\n\n")
  stop("Cannot deploy with missing data files.")
}

# ── Verify rsconnect account is configured ────────────────────────────────
accts <- rsconnect::accounts()
if (nrow(accts) == 0) {
  stop(paste(
    "No shinyapps.io account configured.",
    "Run rsconnect::setAccountInfo() with credentials from",
    "https://www.shinyapps.io/admin/#/tokens — see the comments at the",
    "top of this script for instructions."
  ))
}
cat(sprintf("Deploying as: %s\n", accts$name[1]))

# ── Files to include in the bundle ────────────────────────────────────────
# rsconnect::deployApp() expects appFiles as paths RELATIVE to appDir,
# without any leading directory prefix. We work inside dashboard/ to build
# the list, then deploy.
deploy_files <- c(
  "app.R",
  "global.R",
  file.path("R",    list.files(file.path(dashboard_dir, "R"))),
  file.path("data", list.files(file.path(dashboard_dir, "data")))
)
if (dir.exists(file.path(dashboard_dir, "www"))) {
  deploy_files <- c(deploy_files,
                    file.path("www",
                              list.files(file.path(dashboard_dir, "www"))))
}

# Verify each file actually exists relative to the dashboard directory
exists_check <- file.exists(file.path(dashboard_dir, deploy_files))
if (!all(exists_check)) {
  cat("Missing files:\n")
  cat(paste0("  ", deploy_files[!exists_check]), sep = "\n")
  stop("Some appFiles do not exist")
}

cat(sprintf("\nDeploying %d files:\n", length(deploy_files)))
for (f in deploy_files) cat("  ", f, "\n")

# ── Deploy ─────────────────────────────────────────────────────────────────
rsconnect::deployApp(
  appDir         = dashboard_dir,
  appFiles       = deploy_files,
  appName        = APP_NAME,
  appTitle       = APP_TITLE,
  account        = "amertens",
  server         = "shinyapps.io",
  forceUpdate    = TRUE,
  launch.browser = FALSE
)

cat("\n✓ Deployment complete.\n")
cat(sprintf("  Visit: https://%s.shinyapps.io/%s/\n",
            accts$name[1], APP_NAME))
