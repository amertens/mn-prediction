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
# Registry KEYS, not display labels: get_country_configs() keys Sierra Leone as
# "SierraLeone". Looking it up by the display label "Sierra Leone" silently
# yields NULL, whose zero-length $country later collides with the 14-row Admin2
# column ("arguments imply differing number of rows: 0, 14") and drops the whole
# country. Same pitfall documented at _targets.R:1429.
COUNTRIES <- c("Gambia", "Ghana", "SierraLeone", "Malawi")

classify_who <- function(prev, oc_tag) {
  th <- if (!is.null(meta$who_thresholds)) meta$who_thresholds[[oc_tag]] else NULL
  if (is.null(th)) return(NA_character_)
  vapply(prev, function(p) {
    if (!is.finite(p)) return(NA_character_)
    if (p < th["none"]) "Low" else if (p < th["mild"]) "Mild"
    else if (p < th["moderate"]) "Moderate" else "Severe"
  }, character(1))
}

rows   <- list()
failed <- character(0)
for (co in COUNTRIES) {
  cc <- cfgs[[co]]
  if (is.null(cc)) stop("No country config registered under key '", co,
                        "'. Known keys: ", paste(names(cfgs), collapse = ", "))
  low   <- tolower(co)
  label <- cc$country                      # display label for the dashboard join
  gee <- tryCatch(tar_read_raw(paste0("area_covariates_", low))$gee_admin2,
                  error = function(e) NULL)
  if (is.null(gee)) {
    failed <- c(failed, sprintf("%s: no area_covariates_%s target", label, low)); next
  }
  n_ok <- 0L
  for (ocn in names(cc$outcomes)) { oc <- cc$outcomes[[ocn]]
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", low, "_", ocn)), error = function(e) NULL)
    sv <- tryCatch(tar_read_raw(paste0("svy_admin2_", low, "_", ocn)), error = function(e) NULL)
    if (is.null(od) || is.null(sv)) {
      failed <- c(failed, sprintf("%s/%s: missing outcome_data/svy_admin2 target", label, ocn)); next
    }
    # Record WHY a cell drops out instead of swallowing it: a silent NULL here is
    # what hid Sierra Leone's absence for seven weeks.
    fc <- tryCatch(ar_full_coverage(od, sv, gee, cc, oc),
                   error = function(e) structure(list(), err = conditionMessage(e)))
    if (!is.data.frame(fc)) {
      failed <- c(failed, sprintf("%s/%s: %s", label, ocn,
                  attr(fc, "err") %||% "returned NULL (too few surveyed areas?)")); next
    }
    rows[[length(rows) + 1]] <- data.frame(
      country = label, outcome = ocn, Admin2 = fc$Admin2,
      pred_prev = fc$pred_prev, obs_prev = fc$direct_prev,
      ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_,
      n_survey = fc$n_svy, who_class = classify_who(fc$pred_prev, ocn),
      stringsAsFactors = FALSE)
    n_ok <- n_ok + 1L
  }
  cat(sprintf("  %-14s %d/%d outcomes\n", label, n_ok, length(cc$outcomes)))
}
if (length(failed)) {
  cat("\nCells not produced:\n"); cat(paste0("  - ", failed, collapse = "\n"), "\n\n")
}
out <- do.call(rbind, rows)
saveRDS(out, "dashboard/data/admin2_recipe_predictions.rds")
cat(sprintf("Wrote dashboard/data/admin2_recipe_predictions.rds: %d rows, %d country-outcomes, %d polygons\n",
            nrow(out), length(unique(paste(out$country, out$outcome))), length(unique(out$Admin2))))
