# =============================================================================
# 01_prepare_dashboard_data.R
#
# Extracts pipeline outputs from _targets_full/objects/ and prepares
# dashboard-ready RDS files in dashboard/data/.
#
# Run this script to refresh the dashboard data after pipeline updates.
# The dashboard never reads from _targets_full directly — only from
# the curated files this script produces.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(sf)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Configuration ───────────────────────────────────────────────────────────
TARGETS_DIR  <- here::here("_targets_full", "objects")
DASHBOARD_DATA <- here::here("dashboard", "data")
dir.create(DASHBOARD_DATA, showWarnings = FALSE, recursive = TRUE)

countries <- c("gambia", "ghana", "sierraleone", "malawi")
country_labels <- c(
  gambia      = "Gambia",
  ghana       = "Ghana",
  sierraleone = "Sierra Leone",
  malawi      = "Malawi"
)
gadm_codes <- c(gambia = "GMB", ghana = "GHA", sierraleone = "SLE", malawi = "MWI")
survey_years <- c(gambia = 2018, ghana = 2017, sierraleone = 2013, malawi = 2016)

outcome_labels <- c(
  child_vitA   = "Vitamin A deficiency (children)",
  women_vitA   = "Vitamin A deficiency (women)",
  child_iron   = "Iron deficiency (children)",
  women_iron   = "Iron deficiency (women)",
  women_folate = "Folate deficiency (women)",
  women_b12    = "Vitamin B12 deficiency (women)",
  child_zinc   = "Zinc deficiency (children)",
  women_zinc   = "Zinc deficiency (women)"
)

outcome_short <- c(
  child_vitA   = "Vit A (child)",
  women_vitA   = "Vit A (women)",
  child_iron   = "Iron (child)",
  women_iron   = "Iron (women)",
  women_folate = "Folate (women)",
  women_b12    = "B12 (women)",
  child_zinc   = "Zinc (child)",
  women_zinc   = "Zinc (women)"
)

# WHO public health classification thresholds (prevalence as proportion)
# Used for vitamin A deficiency (RBP < 0.70 µmol/L); other outcomes use
# similar tiered thresholds adapted from WHO/IZiNCG/Lancet guidance.
who_thresholds <- list(
  child_vitA   = c(none = 0.02, mild = 0.10, moderate = 0.20),
  women_vitA   = c(none = 0.02, mild = 0.10, moderate = 0.20),
  child_iron   = c(none = 0.05, mild = 0.20, moderate = 0.40),
  women_iron   = c(none = 0.05, mild = 0.20, moderate = 0.40),
  women_folate = c(none = 0.05, mild = 0.20, moderate = 0.40),
  women_b12    = c(none = 0.05, mild = 0.20, moderate = 0.40),
  child_zinc   = c(none = 0.10, mild = 0.20, moderate = 0.30),
  women_zinc   = c(none = 0.10, mild = 0.20, moderate = 0.30)
)

# ── Helper: safe RDS load ───────────────────────────────────────────────────
safe_read <- function(target_name) {
  path <- file.path(TARGETS_DIR, target_name)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) {
    warning(sprintf("Failed to read %s: %s", target_name, e$message))
    NULL
  })
}

# ── Helper: classify prevalence into WHO public health categories ──────────
classify_who <- function(prev, outcome) {
  th <- who_thresholds[[outcome]]
  if (is.null(th) || is.na(prev)) return(NA_character_)
  if (prev < th["none"])     return("Low")
  if (prev < th["mild"])     return("Mild")
  if (prev < th["moderate"]) return("Moderate")
  return("Severe")
}

# =============================================================================
# 1. ADMIN-2 PREDICTIONS (the core data for the map explorer)
# =============================================================================
# For each country × outcome, we need:
#   - Admin2 name
#   - SuperLearner predicted prevalence (continuous probability average)
#   - Survey-weighted observed prevalence (NA where unsurveyed)
#   - Conformal CI lower and upper bounds
#   - Number of survey observations per district
#   - WHO classification

cat("\n── Building Admin-2 prediction tables ──\n")

admin2_predictions <- list()

for (ctry in countries) {
  for (oc in names(outcome_labels)) {
    suffix <- paste0(ctry, "_", oc)

    sl_pred <- safe_read(paste0("admin2_sl_", suffix))
    svy_obs <- safe_read(paste0("svy_admin2_", suffix))
    conf_ci <- safe_read(paste0("conformal_ci_", suffix))

    if (is.null(sl_pred) || nrow(sl_pred) == 0) next

    # Standardize column names; admin2_sl_* has Admin2 + pred_prev
    if (!"Admin2" %in% colnames(sl_pred)) next

    # Find prediction column (admin2_sl_* uses 'sl_prev')
    pred_col <- intersect(c("sl_prev", "pred_prev", "predicted_prev", "yhat", "fit"),
                          colnames(sl_pred))[1]
    if (is.na(pred_col) || is.null(pred_col)) next

    df <- data.frame(
      country = country_labels[ctry],
      outcome = oc,
      Admin1  = if ("Admin1" %in% colnames(sl_pred)) sl_pred$Admin1 else NA_character_,
      Admin2  = sl_pred$Admin2,
      pred_prev = sl_pred[[pred_col]],
      stringsAsFactors = FALSE
    )

    # Merge survey-weighted observed prevalence (uses 'svy_prev', 'n_svy')
    if (!is.null(svy_obs) && nrow(svy_obs) > 0) {
      svy_pred_col <- intersect(c("svy_prev", "obs_prev", "prev", "prevalence"),
                                colnames(svy_obs))[1]
      svy_n_col <- intersect(c("n_svy", "n", "n_obs", "n_indiv"),
                              colnames(svy_obs))[1]
      svy_keep <- data.frame(
        Admin2 = svy_obs$Admin2,
        obs_prev = if (!is.na(svy_pred_col)) svy_obs[[svy_pred_col]] else NA_real_,
        n_survey = if (!is.na(svy_n_col)) svy_obs[[svy_n_col]] else NA_integer_,
        stringsAsFactors = FALSE
      )
      df <- left_join(df, svy_keep, by = "Admin2")
    } else {
      df$obs_prev <- NA_real_
      df$n_survey <- NA_integer_
    }

    # Merge Admin-1 conformal CIs (these are Admin-1, not Admin-2 — broadcast
    # to constituent districts via Admin1 join when possible)
    df$ci_lo <- NA_real_
    df$ci_hi <- NA_real_
    df$ci_width <- NA_real_

    if (!is.null(conf_ci) && !is.null(conf_ci$admin1_ci) &&
        nrow(conf_ci$admin1_ci) > 0) {
      a1_ci <- conf_ci$admin1_ci
      # Standardize admin1_ci column names
      ci_lo_col <- intersect(c("ci_lo", "lower", "lo"), colnames(a1_ci))[1]
      ci_hi_col <- intersect(c("ci_hi", "upper", "hi"), colnames(a1_ci))[1]
      a1_keep <- data.frame(
        Admin1 = a1_ci$Admin1,
        ci_lo = if (!is.na(ci_lo_col)) a1_ci[[ci_lo_col]] else NA_real_,
        ci_hi = if (!is.na(ci_hi_col)) a1_ci[[ci_hi_col]] else NA_real_,
        stringsAsFactors = FALSE
      )
      a1_keep$ci_width <- a1_keep$ci_hi - a1_keep$ci_lo

      # Join by Admin1 — broadcasts CI to all admin-2 in same admin-1
      df <- df |>
        select(-ci_lo, -ci_hi, -ci_width) |>
        left_join(a1_keep, by = "Admin1")
    }

    # WHO classification on point estimate
    df$who_class <- vapply(df$pred_prev, classify_who, character(1),
                            outcome = oc)

    admin2_predictions[[suffix]] <- df
    cat(sprintf("  %s: %d districts\n", suffix, nrow(df)))
  }
}

admin2_all <- bind_rows(admin2_predictions)
cat(sprintf("\nTotal Admin-2 prediction rows: %d\n", nrow(admin2_all)))

saveRDS(admin2_all, file.path(DASHBOARD_DATA, "admin2_predictions.rds"))


# =============================================================================
# 1b. TRANSPORTABILITY (leave-one-country-out) Admin-2 predictions
# =============================================================================
# Out-of-sample modeled prevalence for each country, produced by training the
# universal/parsimonious area-level model (R/transportability_area.R) on the
# OTHER countries only. Powers the "transportability error" difference map:
# how far a model transported from other countries lands from the local survey.
cat("\n── Building transportability (LOCO) Admin-2 predictions ──\n")
tryCatch({
  source(here::here("R", "transportability_area.R"))
  loco_outcomes <- c("child_vitA", "child_iron", "women_vitA", "women_iron")

  cov_list <- lapply(countries, function(cn) {
    g <- safe_read(paste0("gee_admin2_", cn))
    m <- safe_read(paste0("merged_", cn))
    if (is.null(g) || is.null(m)) return(NULL)
    build_admin2_covariates(g, m)
  })
  names(cov_list) <- countries
  cov_list <- cov_list[!vapply(cov_list, is.null, logical(1))]

  loco_rows <- list()
  for (oc in loco_outcomes) {
    sl <- lapply(names(cov_list), function(cn) safe_read(paste0("svy_admin2_", cn, "_", oc)))
    names(sl) <- names(cov_list)
    pooled <- assemble_area_transport(sl, cov_list, oc)
    if (is.null(pooled)) next
    res <- run_area_transport_loco(pooled, AREA_TRANSPORT_RECIPE)
    if (!is.null(res$predictions)) loco_rows[[oc]] <- res$predictions
  }
  loco_all <- bind_rows(loco_rows)
  if (nrow(loco_all) > 0) {
    loco_all <- loco_all |>
      transmute(country, outcome, Admin2,
                loco_modeled_prev = modeled_prev, survey_prev, n_svy)
    saveRDS(loco_all, file.path(DASHBOARD_DATA, "transportability_loco.rds"))
    cat(sprintf("  Saved %d LOCO predictions across %d outcomes\n",
                nrow(loco_all), length(loco_rows)))
  } else {
    cat("  No LOCO predictions produced (insufficient inputs)\n")
  }
}, error = function(e) cat(sprintf("  LOCO prep skipped: %s\n", e$message)))


# =============================================================================
# 2. NATIONAL ESTIMATES (pre-computed from pipeline)
# =============================================================================
# Source: results/tables/national_estimates_all.csv
# Columns: country, outcome, n, obs_prev, model_n, ci_lo, ci_hi, pred_prev, ...

cat("\n── Building national estimates table ──\n")

natl_path <- here::here("results", "tables", "national_estimates_all.csv")
if (file.exists(natl_path)) {
  natl <- read.csv(natl_path, stringsAsFactors = FALSE)
  saveRDS(natl, file.path(DASHBOARD_DATA, "national_estimates.rds"))
  cat(sprintf("  %d national estimate rows saved\n", nrow(natl)))
} else {
  warning("national_estimates_all.csv not found")
}


# =============================================================================
# 3. CV PERFORMANCE (for methods explainer / model card)
# =============================================================================

cat("\n── Building CV performance table ──\n")

cv_path <- here::here("results", "tables", "cv_performance_all.csv")
if (file.exists(cv_path)) {
  cv_perf <- read.csv(cv_path, stringsAsFactors = FALSE)
  saveRDS(cv_perf, file.path(DASHBOARD_DATA, "cv_performance.rds"))
  cat(sprintf("  %d CV performance rows saved\n", nrow(cv_perf)))
}


# =============================================================================
# 4. ADMIN-2 GADM BOUNDARIES (for choropleth maps)
# =============================================================================
# Cached in data/external_cache/gadm/ — load and save as simplified geometries

cat("\n── Building Admin-2 boundary geopackages ──\n")

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' required", pkg))
  }
}
require_pkg("geodata")

admin2_boundaries <- list()
admin1_boundaries <- list()

gadm_cache <- here::here("data", "external_cache", "gadm")
dir.create(gadm_cache, showWarnings = FALSE, recursive = TRUE)

for (ctry in countries) {
  iso3 <- gadm_codes[ctry]
  cat(sprintf("  Loading GADM for %s (%s)...\n", country_labels[ctry], iso3))

  # Admin-2
  gadm2 <- tryCatch(
    geodata::gadm(iso3, level = 2, path = gadm_cache),
    error = function(e) NULL
  )
  if (!is.null(gadm2)) {
    sf2 <- sf::st_as_sf(gadm2)
    if (!"Admin2" %in% colnames(sf2)) sf2$Admin2 <- sf2$NAME_2
    if (!"Admin1" %in% colnames(sf2)) sf2$Admin1 <- sf2$NAME_1
    sf2 <- sf2[, c("Admin1", "Admin2", "geometry")]
    sf2$country <- country_labels[ctry]
    # Simplify geometry to reduce file size (tolerance in degrees ≈ 100m at equator)
    sf2 <- sf::st_simplify(sf2, dTolerance = 0.001, preserveTopology = TRUE)
    admin2_boundaries[[ctry]] <- sf2
    cat(sprintf("    Admin-2: %d polygons\n", nrow(sf2)))
  }

  # Admin-1
  gadm1 <- tryCatch(
    geodata::gadm(iso3, level = 1, path = gadm_cache),
    error = function(e) NULL
  )
  if (!is.null(gadm1)) {
    sf1 <- sf::st_as_sf(gadm1)
    if (!"Admin1" %in% colnames(sf1)) sf1$Admin1 <- sf1$NAME_1
    sf1 <- sf1[, c("Admin1", "geometry")]
    sf1$country <- country_labels[ctry]
    sf1 <- sf::st_simplify(sf1, dTolerance = 0.001, preserveTopology = TRUE)
    admin1_boundaries[[ctry]] <- sf1
    cat(sprintf("    Admin-1: %d polygons\n", nrow(sf1)))
  }
}

saveRDS(admin2_boundaries, file.path(DASHBOARD_DATA, "admin2_boundaries.rds"))
saveRDS(admin1_boundaries, file.path(DASHBOARD_DATA, "admin1_boundaries.rds"))


# =============================================================================
# 5. POPULATION DATA (for "people at risk" counts)
# =============================================================================
# WorldPop population per Admin-2 district. We use the pre-extracted values
# from the external predictor cache when available, falling back to NA.

cat("\n── Building Admin-2 population table ──\n")

pop_rows <- list()
for (ctry in countries) {
  ext_file <- file.path(here::here("data", "external_cache"),
                        paste0(ctry, "_external_predictors.rds"))
  if (!file.exists(ext_file)) next

  ext <- tryCatch(readRDS(ext_file), error = function(e) NULL)
  if (is.null(ext)) next

  # WorldPop columns are typically named worldpop_pop_total or pop_total
  pop_col <- intersect(c("worldpop_pop_total", "pop_total", "worldpop_total"),
                        colnames(ext))[1]
  if (is.na(pop_col)) next

  # CRITICAL: deduplicate by Admin2 — Malawi cache has cartesian-join inflation
  # (16.7M rows from duplicate "Lake Malawi" / "Lake Chilwa" polygon names).
  # We take the first value per unique Admin2 (population is constant across
  # duplicates by definition).
  ext_dedup <- ext[!duplicated(ext$Admin2),
                    intersect(colnames(ext), c("Admin1", "Admin2", pop_col)),
                    drop = FALSE]
  cat(sprintf("    %s: deduplicated %d → %d rows\n",
              country_labels[ctry], nrow(ext), nrow(ext_dedup)))

  pop_rows[[ctry]] <- data.frame(
    country = country_labels[ctry],
    Admin2 = ext_dedup$Admin2,
    Admin1 = if ("Admin1" %in% colnames(ext_dedup)) ext_dedup$Admin1 else NA_character_,
    population = ext_dedup[[pop_col]],
    pop_year = survey_years[ctry],  # WorldPop pulled at survey year
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s: %d districts with population data\n", country_labels[ctry],
              sum(!is.na(pop_rows[[ctry]]$population))))
}

if (length(pop_rows) > 0) {
  pop_all <- bind_rows(pop_rows)

  # ── 2023 projection ─────────────────────────────────────────────────────
  # Apply country-specific annual growth rates from World Bank Population
  # Estimates (2018–2023 averages). This is a uniform district-level
  # projection and does not account for differential urbanization or
  # internal migration — a flag is set in metadata so the dashboard can
  # surface this caveat to users.
  growth_rates <- c(
    Gambia       = 0.029,   # 2.9%/yr
    Ghana        = 0.020,   # 2.0%/yr
    `Sierra Leone` = 0.021, # 2.1%/yr
    Malawi       = 0.026    # 2.6%/yr
  )
  TARGET_YEAR <- 2023L
  pop_all$population_2023 <- mapply(function(pop, ctry, yr) {
    g <- growth_rates[[ctry]]
    if (is.null(g) || is.na(pop)) return(NA_real_)
    pop * (1 + g)^(TARGET_YEAR - yr)
  }, pop_all$population, pop_all$country, pop_all$pop_year)

  # ── Age/sex demographic scaling ────────────────────────────────────────
  # The pipeline predicts prevalence among specific subgroups (children
  # 6–59 months, women 15–49 years), but our WorldPop denominator is total
  # population. Use country-specific demographic shares from UN World
  # Population Prospects 2022 to compute subgroup-specific populations.
  # These are national-level shares applied uniformly to all districts —
  # a simplification that ignores within-country demographic heterogeneity
  # but is much better than using total population.
  demog_shares <- list(
    # share = fraction of total population in each subgroup
    Gambia       = list(child_6_59m = 0.13, women_15_49 = 0.245),
    Ghana        = list(child_6_59m = 0.10, women_15_49 = 0.255),
    `Sierra Leone` = list(child_6_59m = 0.12, women_15_49 = 0.245),
    Malawi       = list(child_6_59m = 0.13, women_15_49 = 0.240)
  )
  # Sources: UN WPP 2022 single-year-age estimates for relevant year.
  # 6–59 months ≈ ages 0–4 minus first 6 months of births ≈ ~88% of under-5.

  pop_all$pop_child <- NA_real_
  pop_all$pop_women <- NA_real_
  pop_all$pop_child_2023 <- NA_real_
  pop_all$pop_women_2023 <- NA_real_
  for (i in seq_len(nrow(pop_all))) {
    s <- demog_shares[[as.character(pop_all$country[i])]]
    if (is.null(s)) next
    pop_all$pop_child[i] <- pop_all$population[i] * s$child_6_59m
    pop_all$pop_women[i] <- pop_all$population[i] * s$women_15_49
    pop_all$pop_child_2023[i] <- pop_all$population_2023[i] * s$child_6_59m
    pop_all$pop_women_2023[i] <- pop_all$population_2023[i] * s$women_15_49
  }

  saveRDS(pop_all, file.path(DASHBOARD_DATA, "admin2_population.rds"))
  cat(sprintf("\n  Added 2023 projection (growth rates: %s)\n",
              paste(sprintf("%s=%.1f%%", names(growth_rates),
                            growth_rates * 100), collapse = ", ")))
  cat(sprintf("  Added age/sex-specific populations (child 6-59m, women 15-49)\n"))
} else {
  warning("No population data found — population-weighted counts will be unavailable")
}


# =============================================================================
# 6. CÔTE D'IVOIRE OUT-OF-SAMPLE PREDICTIONS
# =============================================================================
# The pipeline produces oos_civ_<outcome> targets — predictions for an
# unsurveyed country (Côte d'Ivoire) using the pooled multi-country model.
# Here we extract them for the dashboard.

cat("\n── Building Côte d'Ivoire OOS predictions ──\n")

oos_rows <- list()
oos_outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

# Try to load Côte d'Ivoire boundaries
civ_bnd <- NULL
civ_admin2_sf <- NULL
civ_gadm <- tryCatch(
  geodata::gadm("CIV", level = 2, path = gadm_cache),
  error = function(e) NULL
)
if (!is.null(civ_gadm)) {
  civ_admin2_sf <- sf::st_as_sf(civ_gadm)
  civ_admin2_sf$Admin2 <- civ_admin2_sf$NAME_2
  civ_admin2_sf$Admin1 <- civ_admin2_sf$NAME_1
  civ_admin2_sf <- civ_admin2_sf[, c("Admin1", "Admin2", "geometry")]
  civ_admin2_sf <- sf::st_simplify(civ_admin2_sf, dTolerance = 0.001,
                                     preserveTopology = TRUE)
  cat(sprintf("  Côte d'Ivoire: %d Admin-2 polygons\n", nrow(civ_admin2_sf)))
}

for (oc in oos_outcomes) {
  oos <- safe_read(paste0("oos_civ_", oc))
  if (is.null(oos)) next

  # OOS targets are typically a list with $predictions data.frame
  preds <- if (is.data.frame(oos)) oos
            else if (!is.null(oos$predictions)) oos$predictions
            else if (!is.null(oos$admin2)) oos$admin2
            else NULL
  if (is.null(preds) || nrow(preds) == 0) next

  pred_col <- intersect(c("pred_prev", "predicted_prev", "yhat", "fit"),
                          colnames(preds))[1]
  ci_lo_col <- intersect(c("ci_lo", "lower"), colnames(preds))[1]
  ci_hi_col <- intersect(c("ci_hi", "upper"), colnames(preds))[1]

  if (is.na(pred_col)) next

  oos_rows[[oc]] <- data.frame(
    outcome = oc,
    Admin2  = preds$Admin2,
    Admin1  = if ("Admin1" %in% colnames(preds)) preds$Admin1 else NA_character_,
    pred_prev = preds[[pred_col]],
    ci_lo   = if (!is.na(ci_lo_col)) preds[[ci_lo_col]] else NA_real_,
    ci_hi   = if (!is.na(ci_hi_col)) preds[[ci_hi_col]] else NA_real_,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s: %d districts\n", oc, nrow(oos_rows[[oc]])))
}

if (length(oos_rows) > 0) {
  oos_all <- bind_rows(oos_rows)
  saveRDS(list(predictions = oos_all, boundaries = civ_admin2_sf),
          file.path(DASHBOARD_DATA, "oos_cote_divoire.rds"))
  cat(sprintf("  Saved %d total OOS rows for Côte d'Ivoire\n", nrow(oos_all)))
} else {
  saveRDS(list(predictions = data.frame(), boundaries = civ_admin2_sf),
          file.path(DASHBOARD_DATA, "oos_cote_divoire.rds"))
  cat("  No OOS prediction targets found — saving empty placeholder\n")
}


# =============================================================================
# 7. SHAP DISTRICT FACTORS + DOMAIN ABLATION
# =============================================================================
# These come from the new pipeline outputs. If not yet produced, save NULL.

cat("\n── Loading SHAP and domain importance ──\n")

shap_path <- here::here("results", "tables", "shap_district_factors.csv")
varimp_path <- here::here("results", "tables", "single_var_importance.csv")
ablation_path <- here::here("results", "tables", "domain_ablation_all.csv")

shap_factors <- if (file.exists(shap_path))
  read.csv(shap_path, stringsAsFactors = FALSE) else NULL
varimp <- if (file.exists(varimp_path))
  read.csv(varimp_path, stringsAsFactors = FALSE) else NULL
ablation <- if (file.exists(ablation_path))
  read.csv(ablation_path, stringsAsFactors = FALSE) else NULL

saveRDS(list(shap = shap_factors, varimp = varimp, ablation = ablation),
        file.path(DASHBOARD_DATA, "importance.rds"))
cat(sprintf("  SHAP rows: %d, varimp rows: %d, ablation rows: %d\n",
            nrow(shap_factors %||% data.frame()),
            nrow(varimp %||% data.frame()),
            nrow(ablation %||% data.frame())))


# =============================================================================
# 8. METADATA (labels, lookups, version info)
# =============================================================================

metadata <- list(
  countries        = country_labels,
  gadm_codes       = gadm_codes,
  survey_years     = survey_years,
  outcome_labels   = outcome_labels,
  outcome_short    = outcome_short,
  who_thresholds   = who_thresholds,
  fieldwork_dates  = list(
    Gambia       = "13 March – 4 May 2018",
    Ghana        = "27 April – 9 June 2017",
    `Sierra Leone` = "11 November – 2 December 2013",
    Malawi       = "8 December 2015 – 16 February 2016"
  ),
  build_timestamp  = Sys.time(),
  build_targets_store = "_targets_full"
)

saveRDS(metadata, file.path(DASHBOARD_DATA, "metadata.rds"))

cat("\n=== Dashboard data preparation complete ===\n")
cat(sprintf("Output directory: %s\n", DASHBOARD_DATA))
cat("Files created:\n")
for (f in list.files(DASHBOARD_DATA)) {
  size_kb <- file.size(file.path(DASHBOARD_DATA, f)) / 1024
  cat(sprintf("  %-40s  %6.1f KB\n", f, size_kb))
}
