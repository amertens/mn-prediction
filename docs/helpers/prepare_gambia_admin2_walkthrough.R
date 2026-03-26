# =============================================================================
# docs/helpers/prepare_gambia_admin2_walkthrough.R
#
# Lightweight helper script that prepares data for the stakeholder walkthrough
# QMD without requiring the full SL modeling pipeline.
#
# What it does:
#   1. Loads the Gambia merged dataset
#   2. Computes survey-weighted Admin-2 VAD prevalence (children)
#   3. Downloads GADM Admin-2 boundaries
#   4. Checks for full pipeline outputs; if missing, creates synthetic
#      predictions matched to the known performance metrics (r = 0.966,
#      MAE = 2.55 pp) so the walkthrough renders with realistic visuals.
#   5. Saves all outputs to docs/helpers/walkthrough_data/
#
# Run this script ONCE from RStudio before rendering the QMD:
#   source("docs/helpers/prepare_gambia_admin2_walkthrough.R")
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(srvyr)
  library(survey)
  library(sf)
  library(geodata)
  library(readr)
  library(labelled)
})

options(survey.lonely.psu = "adjust")

out_dir <- here::here("docs", "helpers", "walkthrough_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Preparing Gambia Admin-2 walkthrough data ===\n\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("[1] Loading Gambia merged dataset...\n")
d_raw <- readRDS(here::here("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"))
if (inherits(d_raw, "sf")) d_raw <- sf::st_drop_geometry(d_raw)

# Children with non-missing VAD outcome
d_child <- d_raw[d_raw$gw_child_flag == 1L & !is.na(d_raw$gw_cRBPAdjThurn), ]
cat(sprintf("  Children with RBP data: %d\n", nrow(d_child)))
cat(sprintf("  Admin-2 units: %d\n", length(unique(d_child$Admin2))))

# ── 2. Survey-weighted Admin-2 VAD prevalence ───────────────────────────────
cat("\n[2] Computing survey-weighted Admin-2 prevalence...\n")

svy_df <- d_child %>%
  dplyr::select(gw_cnum, gw_svy_weight, Admin1, Admin2, gw_cVAD_Thurn) %>%
  dplyr::filter(!is.na(gw_cVAD_Thurn), !is.na(Admin2), !is.na(gw_svy_weight)) %>%
  dplyr::mutate(
    gw_cnum = as.factor(gw_cnum),
    gw_cVAD_Thurn = as.numeric(gw_cVAD_Thurn),
    gw_svy_weight = as.numeric(gw_svy_weight)
  )

svy_des <- srvyr::as_survey_design(
  svy_df,
  ids     = gw_cnum,
  weights = gw_svy_weight
)

svy_admin2 <- svy_des %>%
  group_by(Admin2) %>%
  summarise(
    svy_prev = survey_mean(gw_cVAD_Thurn, vartype = c("se", "ci"), na.rm = TRUE),
    n_svy = n(),
    .groups = "drop"
  )

# National
svy_national <- svy_des %>%
  summarise(
    svy_prev = survey_mean(gw_cVAD_Thurn, vartype = c("se", "ci"), na.rm = TRUE),
    n = n()
  )

cat(sprintf("  National VAD prevalence: %.1f%% [%.1f%%, %.1f%%]\n",
            svy_national$svy_prev * 100,
            svy_national$svy_prev_low * 100,
            svy_national$svy_prev_upp * 100))
cat(sprintf("  Admin-2 units with survey estimate: %d\n", nrow(svy_admin2)))

# ── 3. Download GADM boundaries ─────────────────────────────────────────────
cat("\n[3] Downloading GADM Admin-2 boundaries...\n")

poly_a2 <- geodata::gadm(country = "GMB", level = 2, path = tempdir())
poly_a2 <- sf::st_as_sf(poly_a2) %>%
  dplyr::select(NAME_1, NAME_2) %>%
  dplyr::rename(Admin1_gadm = NAME_1, Admin2_gadm = NAME_2) %>%
  sf::st_transform(crs = 4326)

cat(sprintf("  GADM polygons: %d\n", nrow(poly_a2)))

# Fuzzy-match Admin2 names
survey_names <- sort(unique(svy_admin2$Admin2))
gadm_names   <- sort(unique(poly_a2$Admin2_gadm))

# Try direct match first, then fuzzy
match_table <- data.frame(
  Admin2 = survey_names,
  Admin2_gadm = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(survey_names)) {
  sn <- survey_names[i]
  # Direct match
  if (sn %in% gadm_names) {
    match_table$Admin2_gadm[i] <- sn
  } else {
    # Fuzzy match
    dists <- adist(tolower(sn), tolower(gadm_names))[1, ]
    best  <- which.min(dists)
    if (dists[best] <= 3) {
      match_table$Admin2_gadm[i] <- gadm_names[best]
    }
  }
}

n_matched <- sum(!is.na(match_table$Admin2_gadm))
cat(sprintf("  Matched %d / %d survey names to GADM polygons\n",
            n_matched, length(survey_names)))

# ── 4. Check for full pipeline outputs; create synthetic if missing ─────────
cat("\n[4] Checking for full pipeline outputs...\n")

full_table_path <- here::here("results", "admin2", "tables",
                               "gambia_child_vitA_admin2_full_table.csv")
boot_ci_path    <- here::here("results", "admin2", "tables",
                               "gambia_child_vitA_admin2_boot_ci.csv")

has_full_table <- file.exists(full_table_path)
has_boot_ci    <- file.exists(boot_ci_path)

if (has_full_table) {
  cat("  Found full pipeline table -- using real predictions.\n")
  admin2_full <- read_csv(full_table_path, show_col_types = FALSE)
} else {
  cat("  Full pipeline table NOT found.\n")
  cat("  Creating synthetic predictions matched to known performance\n")
  cat("  (r = 0.966, MAE = 2.55 pp) for walkthrough rendering.\n")
  cat("  >>> Run scripts/06_admin2_predictions_map.R for real predictions. <<<\n")

  set.seed(2024)
  # Generate synthetic predictions with known correlation structure
  # Target: r = 0.966 with svy_prev, MAE ~ 2.55 pp
  obs <- svy_admin2$svy_prev
  n_a2 <- length(obs)

  # Use a correlated noise model: pred = obs + noise
  # To achieve r ~ 0.966 with MAE ~ 0.0255
  noise_sd <- 0.03  # tuned to give ~2.5 pp MAE
  noise <- rnorm(n_a2, mean = -0.003, sd = noise_sd)  # slight negative bias (-0.31 pp)
  sl_prev_synth <- pmin(pmax(obs + noise, 0.01), 0.99)

  admin2_full <- svy_admin2 %>%
    dplyr::mutate(
      sl_prev    = sl_prev_synth,
      abs_error  = abs(sl_prev - svy_prev),
      bias       = sl_prev - svy_prev,
      synthetic  = TRUE
    )
}

if (has_boot_ci) {
  cat("  Found bootstrap CI table.\n")
  admin2_ci <- read_csv(boot_ci_path, show_col_types = FALSE)
} else {
  cat("  Bootstrap CIs NOT found.\n")
  cat("  Creating illustrative uncertainty intervals for rendering.\n")

  set.seed(2025)
  admin2_ci <- admin2_full %>%
    dplyr::mutate(
      boot_mean = sl_prev,
      ci_width  = pmax(0.04, abs(0.15 - 0.08 * sl_prev) + rnorm(n(), 0, 0.02)),
      ci_lo     = pmax(0, sl_prev - ci_width / 2),
      ci_hi     = pmin(1, sl_prev + ci_width / 2)
    ) %>%
    dplyr::select(Admin2, boot_mean, ci_lo, ci_hi, ci_width)
}

# ── 5. Save outputs ─────────────────────────────────────────────────────────
cat("\n[5] Saving walkthrough data...\n")

saveRDS(svy_admin2,     file.path(out_dir, "svy_admin2.rds"))
saveRDS(svy_national,   file.path(out_dir, "svy_national.rds"))
saveRDS(poly_a2,        file.path(out_dir, "gadm_admin2.rds"))
saveRDS(match_table,    file.path(out_dir, "name_match_table.rds"))
saveRDS(admin2_full,    file.path(out_dir, "admin2_full.rds"))
saveRDS(admin2_ci,      file.path(out_dir, "admin2_ci.rds"))

is_synthetic <- !has_full_table
saveRDS(is_synthetic,   file.path(out_dir, "is_synthetic.rds"))

cat("\n  Files saved to docs/helpers/walkthrough_data/:\n")
for (f in list.files(out_dir)) cat(sprintf("    %s\n", f))

cat("\n=== Done. You can now render the walkthrough QMD. ===\n")
if (is_synthetic) {
  cat("\n  NOTE: Predictions are SYNTHETIC (matched to known performance).\n")
  cat("  Run the full pipeline (scripts/06_admin2_predictions_map.R)\n")
  cat("  to generate real model predictions.\n")
}
