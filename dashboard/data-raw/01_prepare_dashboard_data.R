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

# Admin-2 key hygiene: drop GADM inland-water polygons and collapse repeated
# Admin-2 names. extract_gee_admin2() already applies this, which is why the
# FH/BYM2 layers are clean, but extract_area_covariates() deliberately keeps one
# row per polygon ("consumers should filter with is_water_admin2() before
# fitting") — and the area/recipe layers built here are those consumers. Without
# this the default map paints deficiency prevalence on Lake Malawi.
source(here::here("R", "admin2_key_hygiene.R"))

# Apply the standard treatment to a dashboard prediction layer, then recompute
# the WHO class from the collapsed prevalence. dedupe_admin2_key() averages
# numeric columns and takes the first value of everything else — taking the
# first who_class would leave a label that disagrees with the averaged
# pred_prev, so it is derived again here.
clean_pred_layer <- function(d, oc, what) {
  if (is.null(d) || !nrow(d)) return(d)
  d <- clean_admin2_keys(d, what)
  d$who_class <- vapply(as.numeric(d$pred_prev), classify_who, character(1),
                        outcome = oc)
  d
}

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
        Admin1 = if ("Admin1" %in% colnames(svy_obs)) as.character(svy_obs$Admin1)
                 else NA_character_,
        Admin2 = svy_obs$Admin2,
        obs_prev = if (!is.na(svy_pred_col)) svy_obs[[svy_pred_col]] else NA_real_,
        n_survey = if (!is.na(svy_n_col)) svy_obs[[svy_n_col]] else NA_integer_,
        stringsAsFactors = FALSE
      )
      # Pair-join where both sides carry Admin1; see the note on the area layer.
      by_svy <- admin2_join_by(df, svy_keep)
      if (identical(by_svy, "Admin2")) svy_keep$Admin1 <- NULL
      df <- left_join(df, svy_keep, by = by_svy)
    } else {
      df$obs_prev <- NA_real_
      df$n_survey <- NA_integer_
    }

    # Conformal prediction intervals. Prefer the DISTRICT-level intervals
    # (conformal_ci$admin2_ci, keyed on Admin2 — no Admin-1 -> Admin-2 broadcast).
    # Fall back to the Admin-1 intervals broadcast to constituent districts only
    # if district-level intervals are unavailable AND Admin1 is known.
    df$ci_lo <- NA_real_
    df$ci_hi <- NA_real_
    df$ci_width <- NA_real_

    if (!is.null(conf_ci) && !is.null(conf_ci$admin2_ci) &&
        nrow(conf_ci$admin2_ci) > 0 && "Admin2" %in% colnames(conf_ci$admin2_ci)) {
      a2_ci <- conf_ci$admin2_ci
      ci_lo_col <- intersect(c("ci_lo", "lower", "lo"), colnames(a2_ci))[1]
      ci_hi_col <- intersect(c("ci_hi", "upper", "hi"), colnames(a2_ci))[1]
      a2_keep <- data.frame(
        Admin1 = if ("Admin1" %in% colnames(a2_ci)) as.character(a2_ci$Admin1)
                 else NA_character_,
        Admin2 = a2_ci$Admin2,
        ci_lo = if (!is.na(ci_lo_col)) a2_ci[[ci_lo_col]] else NA_real_,
        ci_hi = if (!is.na(ci_hi_col)) a2_ci[[ci_hi_col]] else NA_real_,
        stringsAsFactors = FALSE
      )
      a2_keep$ci_width <- a2_keep$ci_hi - a2_keep$ci_lo
      # Resolve the key and prune BEFORE the pipe: computing it inside the
      # left_join() argument list left the mutation of a2_keep racing lazy
      # argument evaluation.
      by_ci <- admin2_join_by(df, a2_keep)
      if (identical(by_ci, "Admin2")) a2_keep$Admin1 <- NULL
      df <- df |>
        select(-ci_lo, -ci_hi, -ci_width) |>
        left_join(a2_keep, by = by_ci)        # district-level, no broadcast
    } else if (!is.null(conf_ci) && !is.null(conf_ci$admin1_ci) &&
               nrow(conf_ci$admin1_ci) > 0 && !all(is.na(df$Admin1))) {
      a1_ci <- conf_ci$admin1_ci
      ci_lo_col <- intersect(c("ci_lo", "lower", "lo"), colnames(a1_ci))[1]
      ci_hi_col <- intersect(c("ci_hi", "upper", "hi"), colnames(a1_ci))[1]
      a1_keep <- data.frame(
        Admin1 = a1_ci$Admin1,
        ci_lo = if (!is.na(ci_lo_col)) a1_ci[[ci_lo_col]] else NA_real_,
        ci_hi = if (!is.na(ci_hi_col)) a1_ci[[ci_hi_col]] else NA_real_,
        stringsAsFactors = FALSE
      )
      a1_keep$ci_width <- a1_keep$ci_hi - a1_keep$ci_lo
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
# 1c. AREA-LEVEL SAE PREDICTIONS (full district coverage, surveyed + unsurveyed)
# =============================================================================
# From area_model_<country>_<outcome>$area_preds — predicts EVERY Admin-2
# polygon (not just surveyed districts), powering the map explorer's
# "Area-level SAE (all districts)" layer. Same column schema as
# admin2_predictions so get_country_admin2()/get_country_admin1() can consume it.

cat("\n── Building area-level SAE map layer (all districts) ──\n")
area_rows <- list()
for (ctry in countries) {
  for (oc in names(outcome_labels)) {
    am <- safe_read(paste0("area_model_", ctry, "_", oc))
    ap <- am$area_preds
    if (is.null(ap) || !"area_pred_prev" %in% colnames(ap)) next
    hs <- if ("has_survey" %in% colnames(ap)) as.logical(ap$has_survey)
          else !is.na(ap$svy_prev)
    area_rows[[paste(ctry, oc)]] <- clean_pred_layer(
      data.frame(
        country = country_labels[ctry], outcome = oc,
        # CARRY Admin1. Without it clean_admin2_keys() falls back to the
        # name-only key and AVERAGES same-named districts in different regions
        # -- exactly what dedupe_admin2_key()'s own comment calls "the bug this
        # migration fixed". Malawi has 243 polygons under 239 names, so the
        # dashboard layer collapsed to 239 rows and the map explorer painted
        # TA Lundu / TA Malemia / TA Ngabu / TA Pemba (four genuinely distinct
        # district pairs) with one shared prediction each. admin2_boundaries.rds
        # already carries Admin1, so the map can join on the pair once this does.
        Admin1 = if ("Admin1" %in% colnames(ap)) as.character(ap$Admin1)
                 else NA_character_,
        Admin2 = ap$Admin2,
        pred_prev = as.numeric(ap$area_pred_prev),
        obs_prev = if ("svy_prev" %in% colnames(ap))
                     ifelse(hs, ap$svy_prev, NA_real_) else NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_,
        n_survey = if ("n_svy" %in% colnames(ap)) ap$n_svy else NA_integer_,
        who_class = NA_character_,
        stringsAsFactors = FALSE),
      oc, sprintf("area layer %s/%s", ctry, oc))
  }
}
if (length(area_rows) > 0) {
  area_all <- bind_rows(area_rows)
  saveRDS(area_all, file.path(DASHBOARD_DATA, "admin2_area_predictions.rds"))
  cat(sprintf("  %d area-level rows across %d country × outcome combos\n",
              nrow(area_all), length(area_rows)))
} else {
  cat("  No area_model_* objects found — skipping area layer.\n")
}


# =============================================================================
# 1d. FAY-HERRIOT / EMPIRICAL-BAYES SAE LAYER (full coverage, with intervals)
# =============================================================================
# Full-coverage area-level model that blends each district's direct survey
# estimate with an indicator-based synthetic prediction, with design-based 95%
# intervals (narrow where surveyed, wide where not). Self-contained builder;
# sourced in an isolated environment so it cannot clobber variables here.

cat("\n── Building Fay-Herriot SAE map layer (all districts, with intervals) ──\n")
tryCatch(
  source(here::here("dashboard", "data-raw", "_build_fh_layer.R"),
         local = new.env()),
  error = function(e) cat("  FH layer build skipped:", conditionMessage(e), "\n"))


# =============================================================================
# 1e. SL -> BYM2 SAE LAYER (full coverage, spatially-smoothed credible intervals)
# =============================================================================
# Feeds the area-level SuperLearner mean into a Bayesian spatial (BYM2) model so
# every district gets a 95% credible interval that borrows strength from
# neighbours — tighter and better-calibrated than the Fay-Herriot intervals
# (see archive/sandbox/sae_sl_hybrid_prototype.R). Self-contained builder; sourced in an
# isolated environment. Requires INLA + spdep; skipped cleanly if unavailable.
cat("\n── Building SL → BYM2 SAE map layer (all districts, spatial intervals) ──\n")
tryCatch(
  source(here::here("dashboard", "data-raw", "_build_bym2_layer.R"),
         local = new.env()),
  error = function(e) cat("  BYM2 layer build skipped:", conditionMessage(e), "\n"))


# =============================================================================
# 1b. TRANSPORTABILITY (leave-one-country-out) Admin-2 predictions
# =============================================================================
# Out-of-sample modeled prevalence for each country, produced by training the
# universal/parsimonious area-level model (R/transportability_area.R) on the
# OTHER countries only. Powers the "transportability error" difference map:
# how far a model transported from other countries lands from the local survey.
cat("\n── Building transportability (LOCO) Admin-2 predictions ──\n")
tryCatch({
  # Prefer the canonical results file produced by
  # scripts/run_area_transportability.R, which uses the HARMONIZED outcome
  # (uniform WHO cutoffs) and the same recipe. This keeps the dashboard's
  # transport-error layer consistent with the analysis deliverables and
  # carries each area's harmonized survey prevalence (loco_survey_prev).
  loco_results <- here::here("results", "transportability", "area_loco_predictions.rds")
  if (file.exists(loco_results)) {
    lr <- readRDS(loco_results)
    loco_all <- lr |>
      transmute(country, outcome, Admin2,
                loco_modeled_prev = modeled_prev,
                loco_survey_prev  = survey_prev, n_svy)
    saveRDS(loco_all, file.path(DASHBOARD_DATA, "transportability_loco.rds"))
    cat(sprintf("  Loaded %d harmonized LOCO predictions from results/\n", nrow(loco_all)))
  } else {
    cat("  results/transportability/area_loco_predictions.rds not found —\n")
    cat("  run scripts/run_area_transportability.R first to generate it.\n")
  }
}, error = function(e) cat(sprintf("  LOCO prep skipped: %s\n", e$message)))


# ── Optimized primary SL: prescreened SuperLearner LOCO results ───────────
# Source: results/tables/sl_prescreened_main.csv (from
#   scripts/run_sl_prescreened_main.R or `sl_prescreened_all` target).
# This is the OPTIMIZED primary SL workflow — fast (~30s per LOCO fold)
# with five-stage feature prescreening. Picks up automatically if the
# CSV exists.
cat("\n── Loading optimized SL prescreened LOCO results ──\n")
sl_pre_path <- here::here("results", "tables", "sl_prescreened_main.csv")
if (file.exists(sl_pre_path)) {
  sl_pre <- read.csv(sl_pre_path, stringsAsFactors = FALSE)
  saveRDS(sl_pre, file.path(DASHBOARD_DATA, "sl_prescreened_loco.rds"))
  cat(sprintf("  %d SL-prescreened LOCO rows saved\n", nrow(sl_pre)))
} else {
  cat("  results/tables/sl_prescreened_main.csv not found - skipping.\n")
}


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
    # Drop GADM inland-water polygons (Malawi ships Lake Malawi as 8 separate
    # Admin-2 features, one per bordering district, plus Lake Chilwa x3, Lake
    # Chiuta and Lake Malombe). Removing the polygon rather than leaving it
    # unpredicted is deliberate: an unpredicted lake renders as a grey "no data"
    # district, which reads as a place we failed to model. Dropped, the basemap
    # shows through and the lake looks like a lake.
    n_before <- nrow(sf2)
    sf2 <- sf2[!is_water_admin2(as.character(sf2$Admin2)), , drop = FALSE]
    if (nrow(sf2) < n_before)
      cat(sprintf("    dropped %d water-body polygon(s)\n", n_before - nrow(sf2)))
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
  # WorldPop bleeds across the shoreline, so the lake polygons carry a nonzero
  # count (Lake Malawi ~1,959 people). Dropping them keeps the denominator on
  # land; it moves Malawi's total by 0.05%.
  ext_dedup <- drop_water_admin2(ext_dedup,
                                 sprintf("population %s", country_labels[ctry]))

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
# 7b. CLUSTER-RESOLUTION COMPARISON (admin-2 vs survey-cluster)
# =============================================================================
# Sensitivity analysis comparing predictive skill when prevalence is modeled at
# administrative-2 resolution vs. the finer survey-cluster GPS buffers. Two
# inputs, both optional (the panel degrades to an empty state if absent):
#   results/cluster_vs_admin2_comparison.csv  — within-country CV: per
#       country x outcome, Pearson r and MAE (pp) at each resolution.
#   results/cluster_vs_admin2_LOCO.csv        — pooled leave-one-country-out
#       transportability, per resolution x method x held-out country.

cat("\n── Loading cluster-resolution comparison ──\n")

cmp_path  <- here::here("results", "cluster_vs_admin2_comparison.csv")
loco_path <- here::here("results", "cluster_vs_admin2_LOCO.csv")

cluster_comparison <- if (file.exists(cmp_path)) {
  cc_df <- read.csv(cmp_path, stringsAsFactors = FALSE)
  # Derived deltas (cluster − admin2): positive r-delta = cluster more skillful;
  # positive MAE-delta = cluster less accurate (MAE: lower is better).
  cc_df$delta_r   <- cc_df$r_cluster   - cc_df$r_admin2
  cc_df$delta_mae <- cc_df$mae_cluster - cc_df$mae_admin2
  cc_df
} else NULL

cluster_loco <- if (file.exists(loco_path)) {
  read.csv(loco_path, stringsAsFactors = FALSE)
} else NULL

saveRDS(
  list(comparison = cluster_comparison, loco = cluster_loco,
       build_time = Sys.time()),
  file.path(DASHBOARD_DATA, "cluster_resolution.rds"))
cat(sprintf("  comparison rows: %d, LOCO rows: %d\n",
            nrow(cluster_comparison %||% data.frame()),
            nrow(cluster_loco %||% data.frame())))


# =============================================================================
# 7b-ii. WHERE THE LEVEL COMES FROM — anchoring, and the choice of unit
# =============================================================================
# The single design choice that most affects the district map, and until now the
# app showed its OUTPUT without ever showing the choice. Anchoring each district
# prediction to its region's design-based survey total more than doubles rank
# correlation (0.170 -> 0.406) and cuts mean absolute bias from 3.6 to 1.5 pp,
# better in 20 of 24 country x outcome cells; anchoring to the NATIONAL total
# instead buys almost nothing (+0.006). Read together with the resolution
# comparison, these say the model supplies the pattern and the survey supplies
# the level -- which is also why a transported map for a country with no survey
# of its own can be read as a ranking but not as a set of prevalences.
anchor_path <- here::here("results", "tables", "admin1_arms.csv")
resolution_path <- here::here("results", "tables", "resolution_comparison.csv")
anchor_arms <- if (file.exists(anchor_path))
  read.csv(anchor_path, stringsAsFactors = FALSE) else NULL
resolution_levels <- if (file.exists(resolution_path))
  read.csv(resolution_path, stringsAsFactors = FALSE) else NULL

saveRDS(list(arms = anchor_arms, levels = resolution_levels,
             build_time = Sys.time()),
        file.path(DASHBOARD_DATA, "anchoring.rds"))
cat(sprintf("\n-- Anchoring / level --\n  anchor arms: %d rows | resolution rows: %d\n",
            nrow(anchor_arms %||% data.frame()),
            nrow(resolution_levels %||% data.frame())))


# =============================================================================
# 7c. TRANSPORT CALIBRATION (predicted vs. true prevalence, multi-level)
# =============================================================================
# Powers the interactive Transportability panel: predicted-vs-observed
# prevalence under leave-one-country-out (LOCO) holdout, selectable across
# outcome, aggregation level, and model/approach. Three optional sources:
#   national, individual-level SL  -> transportability_loco_results.csv
#   national, area-level SL        -> area_loco_comparison.csv (model_type)
#   admin-2, area-level SL         -> results/transportability/area_loco_predictions.rds
# Admin-1 is derived in-app by population-weighting the admin-2 predictions.

cat("\n── Building transport-calibration data ──\n")

tc_nat_indiv <- {
  p <- here::here("results", "tables", "transportability_loco_results.csv")
  if (file.exists(p)) {
    x <- read.csv(p, stringsAsFactors = FALSE)
    data.frame(outcome = x$outcome, held_out = x$held_out,
               true_prev = x$prev_test, pred_prev = x$pred_prev,
               auc = x$auc, calib_slope = x$calib_slope,
               stringsAsFactors = FALSE)
  } else NULL
}

tc_nat_area <- {
  p <- here::here("results", "tables", "area_loco_comparison.csv")
  if (file.exists(p)) {
    x <- read.csv(p, stringsAsFactors = FALSE)
    data.frame(outcome = x$outcome, held_out = x$held_out,
               model_type = x$model_type, true_prev = x$obs_prev,
               pred_prev = x$pred_prev, pearson_r = x$pearson_r,
               mae_pp = x$mae_pp, stringsAsFactors = FALSE)
  } else NULL
}

tc_adm2_area <- {
  p <- here::here("results", "transportability", "area_loco_predictions.rds")
  if (file.exists(p)) {
    x <- readRDS(p)
    data.frame(outcome = x$outcome, country = x$country, Admin2 = x$Admin2,
               true_prev = x$survey_prev, pred_prev = x$modeled_prev,
               n_svy = x$n_svy, stringsAsFactors = FALSE)
  } else NULL
}

saveRDS(list(nat_indiv = tc_nat_indiv, nat_area = tc_nat_area,
             adm2_area = tc_adm2_area, build_time = Sys.time()),
        file.path(DASHBOARD_DATA, "transport_calibration.rds"))
cat(sprintf("  national/indiv: %d rows, national/area: %d rows, admin-2/area: %d rows\n",
            nrow(tc_nat_indiv %||% data.frame()),
            nrow(tc_nat_area %||% data.frame()),
            nrow(tc_adm2_area %||% data.frame())))


# =============================================================================
# 7d. MODEL DIAGNOSTICS (ROC / PR / calibration curves + metrics)
# =============================================================================
# Surfaces the per-model discrimination & calibration diagnostics produced by
# the `diagnostics_all` / `diagnostics_calibrated_all` pipeline targets, which
# write the following CSVs (one row per country × outcome, or per curve point):
#   diagnostics_binary.csv            — ROC-AUC, PR-AUC, Brier skill, ECE/MCE,
#                                       calibration intercept/slope (binary)
#   diagnostics_continuous.csv        — RMSE, MAE, R^2, r (continuous)
#   diagnostics_binary_calibrated.csv — same as binary after Platt recalibration
#   roc_curves.csv                    — fpr, tpr per country × outcome
#   pr_curves.csv                     — recall, precision per country × outcome
#   calibration_tables.csv            — binned mean_pred vs mean_obs
# The dashboard previously only showed ROC-AUC / Brier (via cv_performance);
# this bundle exposes the full diagnostic set.

cat("\n── Building model diagnostics bundle ──\n")

read_csv_if <- function(rel) {
  p <- here::here("results", "tables", rel)
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
saveRDS(model_diagnostics, file.path(DASHBOARD_DATA, "model_diagnostics.rds"))
cat(sprintf(paste0("  binary: %d, continuous: %d, calibrated: %d, ",
                   "roc pts: %d, pr pts: %d, calib bins: %d\n"),
            nrow(model_diagnostics$binary %||% data.frame()),
            nrow(model_diagnostics$continuous %||% data.frame()),
            nrow(model_diagnostics$calibrated %||% data.frame()),
            nrow(model_diagnostics$roc %||% data.frame()),
            nrow(model_diagnostics$pr %||% data.frame()),
            nrow(model_diagnostics$calibration %||% data.frame())))


# =============================================================================
# 7e. METHOD BENCHMARKS (SuperLearner vs small-area-estimation methods)
# =============================================================================
# Exposes the area-level model-comparison tables that were previously computed
# by the pipeline but never displayed:
#   benchmarks_all.csv                     — SL vs baseline/GLM/elastic-net/
#                                            Fay-Herriot/BYM2 under LOCO
#   area_comparison_all.csv                — individual-level SL (aggregated to
#                                            Admin-2) vs area-level SL, within-
#                                            country
#   admin2_error_all.csv                   — per country × outcome Admin-2
#                                            prediction accuracy (MAE/RMSE/r)
#   sl_prescreened_main.csv                — optimized prescreened SL, LOCO
#   transportability_area_loco_metrics.csv — enriched area-transport LOCO metrics

cat("\n── Building method-benchmarks bundle ──\n")

benchmarks <- list(
  benchmarks      = read_csv_if("benchmarks_all.csv"),
  area_comparison = read_csv_if("area_comparison_all.csv"),
  admin2_error    = read_csv_if("admin2_error_all.csv"),
  sl_prescreened  = read_csv_if("sl_prescreened_main.csv"),
  area_transport  = read_csv_if("transportability_area_loco_metrics.csv"),
  build_time      = Sys.time()
)
saveRDS(benchmarks, file.path(DASHBOARD_DATA, "benchmarks.rds"))
cat(sprintf(paste0("  benchmarks: %d, area_comparison: %d, admin2_error: %d, ",
                   "sl_prescreened: %d, area_transport: %d\n"),
            nrow(benchmarks$benchmarks %||% data.frame()),
            nrow(benchmarks$area_comparison %||% data.frame()),
            nrow(benchmarks$admin2_error %||% data.frame()),
            nrow(benchmarks$sl_prescreened %||% data.frame()),
            nrow(benchmarks$area_transport %||% data.frame())))


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
