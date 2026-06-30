#!/usr/bin/env Rscript
# Build full-coverage area-level RECIPE district predictions for the dashboard
# (docs/AREA_LEVEL_RECIPE_SPEC.md). Trains the recipe on surveyed districts and
# predicts EVERY admin-2 polygon using the universal (gee, polygon-zonal) feature
# set — the only one available at unsurveyed polygons. Output schema matches
# dashboard/data/admin2_area_predictions.rds so the app can consume it as a layer.
#
# Run from the repo root:  Rscript dashboard/data-raw/build_recipe_predictions.R
suppressMessages({library(targets); library(glmnet)})
tar_config_set(store = "_targets_full")
source("R/corrected/00_corrected_utils.R"); source("R/corrected/area_recipe.R")
cfgs <- local({ source("R/config.R", local = TRUE); get_country_configs() })
meta <- tryCatch(readRDS("dashboard/data/metadata.rds"), error = function(e) NULL)
COUNTRIES <- c("Gambia", "Ghana", "Sierra Leone", "Malawi")

classify_who <- function(prev, oc_tag) {
  th <- if (!is.null(meta$who_thresholds)) meta$who_thresholds[[oc_tag]] else NULL
  if (is.null(th)) return(NA_character_)
  vapply(prev, function(p) {
    if (!is.finite(p)) return(NA_character_)
    if (p < th["none"]) "Low" else if (p < th["mild"]) "Mild"
    else if (p < th["moderate"]) "Moderate" else "Severe"
  }, character(1))
}

rows <- list()
for (co in COUNTRIES) { cc <- cfgs[[co]]; low <- tolower(gsub(" ", "", co))
  gee <- tryCatch(tar_read_raw(paste0("area_covariates_", low))$gee_admin2, error = function(e) NULL)
  if (is.null(gee)) next
  for (ocn in names(cc$outcomes)) { oc <- cc$outcomes[[ocn]]
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", low, "_", ocn)), error = function(e) NULL)
    sv <- tryCatch(tar_read_raw(paste0("svy_admin2_", low, "_", ocn)), error = function(e) NULL)
    if (is.null(od) || is.null(sv)) next
    fc <- tryCatch(ar_full_coverage(od, sv, gee, cc, oc), error = function(e) NULL)
    if (is.null(fc)) next
    rows[[length(rows) + 1]] <- data.frame(
      country = co, outcome = ocn, Admin2 = fc$Admin2,
      pred_prev = fc$pred_prev, obs_prev = fc$direct_prev,
      ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_,
      n_survey = fc$n_svy, who_class = classify_who(fc$pred_prev, ocn),
      stringsAsFactors = FALSE)
  }
  cat("done:", co, "\n")
}
out <- do.call(rbind, rows)
saveRDS(out, "dashboard/data/admin2_recipe_predictions.rds")
cat(sprintf("Wrote dashboard/data/admin2_recipe_predictions.rds: %d rows, %d country-outcomes, %d polygons\n",
            nrow(out), length(unique(paste(out$country, out$outcome))), length(unique(out$Admin2))))
