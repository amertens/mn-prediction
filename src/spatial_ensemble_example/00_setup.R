# 00_setup.R
# Spatial ensemble learning workflow (simulated data)
# This script:
# - installs/loads packages
# - sets global options and seeds
# - defines a simple project directory structure

options(stringsAsFactors = FALSE)

# Reproducible seed (set once here; downstream scripts rely on it)
set.seed(12345)

# ---- Packages ----
pkgs <- c(
  "data.table",
  "sf",
  "terra",
  "mgcv",
  "glmnet",
  "ranger",
  "xgboost",
  "blockCV",
  "ggplot2"
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Paths ----
dir.create("data", showWarnings = FALSE)
dir.create("data/sim", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

message("Setup complete.")
