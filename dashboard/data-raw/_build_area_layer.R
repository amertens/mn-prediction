# Build the full-coverage area-level SAE map layer from area_model_* targets.
# Extracts area_preds (one row per Admin-2 polygon, surveyed + unsurveyed) into
# admin2_area_predictions.rds, in the same column schema the map helper expects
# (get_country_admin2 / get_country_admin1 in dashboard/R/fct_helpers.R).
# Mirrors section logic added to 01_prepare_dashboard_data.R.
suppressWarnings({
root <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
if (!dir.exists(file.path(root, "_targets_full"))) root <- getwd()
S <- file.path(root, "_targets_full", "objects")
DASH <- file.path(root, "dashboard", "data")
source(file.path(root, "R", "admin2_key_hygiene.R"))

countries <- c(gambia = "Gambia", ghana = "Ghana",
               sierraleone = "Sierra Leone", malawi = "Malawi")
outcome_tags <- c("child_vitA","women_vitA","child_iron","women_iron",
                  "women_folate","women_b12","child_zinc","women_zinc")

# WHO public-health classification (same thresholds as 01_prepare_dashboard_data.R)
who_thresholds <- list(
  child_vitA = c(none=0.02, mild=0.10, moderate=0.20),
  women_vitA = c(none=0.02, mild=0.10, moderate=0.20),
  child_iron = c(none=0.05, mild=0.20, moderate=0.40),
  women_iron = c(none=0.05, mild=0.20, moderate=0.40),
  women_folate = c(none=0.05, mild=0.20, moderate=0.40),
  women_b12 = c(none=0.05, mild=0.20, moderate=0.40),
  child_zinc = c(none=0.10, mild=0.20, moderate=0.30),
  women_zinc = c(none=0.10, mild=0.20, moderate=0.30))
classify_who <- function(prev, oc) {
  th <- who_thresholds[[oc]]
  if (is.null(th) || is.na(prev)) return(NA_character_)
  if (prev < th["none"]) return("Low")
  if (prev < th["mild"]) return("Mild")
  if (prev < th["moderate"]) return("Moderate")
  "Severe"
}

rows <- list()
for (low in names(countries)) {
  for (oc in outcome_tags) {
    obj <- file.path(S, paste0("area_model_", low, "_", oc))
    if (!file.exists(obj)) next
    am <- tryCatch(readRDS(obj), error = function(e) NULL)
    ap <- am$area_preds
    if (is.null(ap) || !"area_pred_prev" %in% colnames(ap)) next
    pred <- as.numeric(ap$area_pred_prev)
    obs  <- if ("svy_prev" %in% colnames(ap)) {
      hs <- if ("has_survey" %in% colnames(ap)) as.logical(ap$has_survey) else !is.na(ap$svy_prev)
      ifelse(hs, ap$svy_prev, NA_real_)
    } else NA_real_
    # Same water/duplicate-name treatment as section 1c of
    # 01_prepare_dashboard_data.R — area_preds is one row per GADM polygon and
    # so still carries the lake polygons and the repeated Admin-2 names.
    layer <- clean_admin2_keys(
      data.frame(
        country   = unname(countries[low]),
        outcome   = oc,
        Admin2    = ap$Admin2,
        pred_prev = pred,
        obs_prev  = obs,
        ci_lo     = NA_real_, ci_hi = NA_real_, ci_width = NA_real_,
        n_survey  = if ("n_svy" %in% colnames(ap)) ap$n_svy else NA_integer_,
        stringsAsFactors = FALSE),
      sprintf("area layer %s/%s", low, oc))
    layer$who_class <- vapply(layer$pred_prev, classify_who, character(1), oc = oc)
    rows[[paste(low, oc)]] <- layer
  }
}
out <- do.call(rbind, rows)
saveRDS(out, file.path(DASH, "admin2_area_predictions.rds"))
cat(sprintf("admin2_area_predictions.rds: %d rows, %d country x outcome combos\n",
            nrow(out), length(rows)))
print(tapply(out$Admin2, out$country, function(x) length(unique(x))))
cat("share of district-rows that are unsurveyed (area model fills these):",
    round(mean(is.na(out$obs_prev)), 2), "\n")
})
