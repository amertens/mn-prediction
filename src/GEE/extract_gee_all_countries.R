# =============================================================================
# src/GEE/extract_gee_all_countries.R
#
# Re-extract admin-2 GEE zonal means via the EE API (extract_gee_ee_api.R) for
# ALL countries, so every country uses the SAME gee_a2_ column names — required
# for the cross-country LOCO common-vocabulary set (the legacy terra export used
# different names than the new EE extraction).
#
# Writes data/GEE/<ConfigKey>_<year>_admin2_gee.csv, backing up any existing
# legacy file to <name>.legacy first. Per-country per-layer cache -> resumable.
#
#   Rscript src/GEE/extract_gee_all_countries.R            # all countries
#   Rscript src/GEE/extract_gee_all_countries.R Ghana      # one
# =============================================================================
suppressPackageStartupMessages({ library(sf); library(dplyr); library(readr); library(here) })
source(here("R", "data_prep.R"))
source(here("src", "GEE", "extract_gee_ee_api.R"))

# ── EE init (same as the Tanzania runner) ────────────────────────────────────
ee_project <- Sys.getenv("EE_PROJECT", unset = "mn-prediction-420517")
rgee_py <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (file.exists(rgee_py)) Sys.setenv(RETICULATE_PYTHON = rgee_py)
reticulate::py_run_string(sprintf("import ee; ee.Initialize(project='%s')", ee_project))
invisible(rgee::ee$Number(1)$getInfo())

# ConfigKey -> (gadm_code, survey_year). All 5 (incl. Tanzania) re-extract under
# the enhanced manifest (yoy + soil 20-50cm); caches were cleared so it rebuilds.
COUNTRIES <- list(
  Gambia      = list(gadm = "GMB", year = 2018L),
  Ghana       = list(gadm = "GHA", year = 2017L),
  SierraLeone = list(gadm = "SLE", year = 2013L),
  Malawi      = list(gadm = "MWI", year = 2016L),
  Tanzania    = list(gadm = "TZA", year = 2010L)
)

args <- commandArgs(trailingOnly = TRUE)
sel  <- if (length(args)) args else names(COUNTRIES)

for (ck in sel) {
  cc <- COUNTRIES[[ck]]; if (is.null(cc)) { message("unknown: ", ck); next }
  message("\n===== ", ck, " (", cc$year, ") =====")
  poly <- load_gadm_cached(cc$gadm, level = 2) |>
    sf::st_as_sf() |> dplyr::select(NAME_1, NAME_2) |>
    dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)
  poly$Admin1 <- trimws(poly$Admin1); poly$Admin2 <- trimws(poly$Admin2)

  adm2 <- extract_gee_layers(poly, cc$year, id_col = "Admin2",
                             col_prefix = "gee_a2_", verbose = TRUE,
                             cache_dir = here("data", "GEE",
                                              paste0(".cache_", tolower(ck), cc$year)))
  out <- here("data", "GEE", sprintf("%s_%d_admin2_gee.csv", ck, cc$year))
  if (file.exists(out) && !file.exists(paste0(out, ".legacy")))
    file.copy(out, paste0(out, ".legacy"))          # preserve the legacy export
  readr::write_csv(adm2, out)
  message(sprintf("  wrote %s (%d districts x %d cols)", out, nrow(adm2), ncol(adm2) - 2L))
}
message("\nDone. Re-run each country's 2_GW_*_data_merge.R to pick up aligned gee_a2_ names.")
