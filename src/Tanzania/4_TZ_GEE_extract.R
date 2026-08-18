# =============================================================================
# src/Tanzania/4_TZ_GEE_extract.R
#
# Auto-extract the Tanzania (2010) GEE proxy layers from the Earth Engine API
# (no manual Code-Editor export). Produces the two CSVs the merge step expects:
#
#   data/GEE/Tanzania_2010_admin2_gee.csv   (admin-2 zonal means -> gee_a2_)
#   data/GEE/TZ2010_buffers_<date>.csv      (cluster buffers      -> gee_)
#
# PREREQUISITES (see src/GEE/README_GEE_API_SETUP.md):
#   - EE API access + rgee installed; run rgee::ee_Initialize() successfully.
#   - Tanzania_GMS_cleaned.rds exists (src/Tanzania/1_GW_Tanzania_data_clean.R)
#     so cluster GPS are available for the buffer product.
#
# Run interactively (recommended first time, to watch per-layer logs):
#   source(here::here("src/Tanzania/4_TZ_GEE_extract.R"))
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(here)
})
source(here("R", "data_prep.R"))                 # load_gadm_cached()
source(here("src", "GEE", "extract_gee_ee_api.R"))

YEAR <- 2010L

# Per-layer cache: successful layers are cached here, so a re-run only pulls the
# missing/failed ones. Set GEE_FORCE=true to re-extract everything, or delete the
# dir. (Failed layers are never cached, so they always retry.)
CACHE_DIR <- here("data", "GEE", ".cache_tza2010")
FORCE     <- tolower(Sys.getenv("GEE_FORCE", "false")) %in% c("1", "true", "yes")

# ── Initialise Earth Engine ────────────────────────────────────────────────
# Account amertens@berkeley.edu, Cloud project mn-prediction-420517. The Python
# env with earthengine-api is the reticulate venv (the old conda rgee env was
# removed). Env broken? See src/GEE/README_GEE_API_SETUP.md.
#
# NOTE: no drive = TRUE — extract_gee_layers() uses rgee::ee_extract(), which
# pulls reductions straight into R. Google Drive is only needed for async
# ee_image_to_drive exports (and enabling it requires ticking the Drive
# permission box at auth time). We don't use those.
if (!requireNamespace("rgee", quietly = TRUE))
  stop("Install rgee first: install.packages('rgee')")

ee_user    <- Sys.getenv("EE_USER",    unset = "amertens@berkeley.edu")
ee_project <- Sys.getenv("EE_PROJECT", unset = "mn-prediction-420517")  # Cloud project

# Point reticulate at the venv that has earthengine-api (if .Renviron's
# EARTHENGINE_PYTHON isn't already doing so).
rgee_py <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (file.exists(rgee_py)) Sys.setenv(RETICULATE_PYTHON = rgee_py)

# Initialize via the Python EE client directly (reads cached credentials from
# ~/.config/earthengine/). This deliberately AVOIDS rgee::ee_Initialize(), whose
# interactive "asset root folder" setup loops forever when the home already
# exists (a known rgee 1.1.7 bug). We never touch the asset home — the
# extraction is read-only — so this is the clean, non-interactive path.
# One-time auth (if the token is stale): run rgee::ee_Authenticate() first.
stopifnot(nzchar(ee_project))
reticulate::py_run_string(sprintf("import ee; ee.Initialize(project='%s')", ee_project))
ee <- rgee::ee
invisible(ee$Number(1)$getInfo())  # cheap sanity call; errors early if not connected

# ── Admin-2 product ────────────────────────────────────────────────────────
poly <- load_gadm_cached("TZA", level = 2) |>
  sf::st_as_sf() |>
  dplyr::select(NAME_1, NAME_2) |>
  dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)
poly$Admin1 <- trimws(poly$Admin1); poly$Admin2 <- trimws(poly$Admin2)

message("\n=== Tanzania admin-2 GEE extraction (", YEAR, ") ===")
adm2 <- extract_gee_layers(poly, YEAR, id_col = "Admin2",
                           col_prefix = "gee_a2_", verbose = TRUE,
                           cache_dir = file.path(CACHE_DIR, "admin2"), force = FORCE)

# Parity report vs Sierra Leone's legacy admin-2 file (naming overlap).
ref <- here("data/GEE/SierraLeone_2013_admin2_gee.csv")
if (file.exists(ref)) {
  ref_cols <- names(readr::read_csv(ref, n_max = 0, show_col_types = FALSE))
  shared   <- intersect(names(adm2), ref_cols)
  message(sprintf("[parity] %d/%d admin-2 cols match SL legacy names (%.0f%%)",
                  length(shared), length(names(adm2)),
                  100 * length(shared) / max(1, length(names(adm2)))))
}
out_adm2 <- here("data/GEE/Tanzania_2010_admin2_gee.csv")
readr::write_csv(adm2, out_adm2)
message(sprintf("  wrote %s (%d districts x %d cols)", out_adm2,
                nrow(adm2), ncol(adm2) - 2L))

# ── Cluster-buffer product ─────────────────────────────────────────────────
clean_rds <- here("data", "IPD", "Tanzania 2010", "Tanzania_GMS_cleaned.rds")
if (file.exists(clean_rds)) {
  d <- readRDS(clean_rds)
  clusters <- d |>
    dplyr::distinct(gw_cnum, latitude, longitude) |>
    dplyr::filter(!is.na(latitude), !is.na(longitude)) |>
    dplyr::rename(cnum = gw_cnum)
  clusters_sf <- sf::st_as_sf(clusters, coords = c("longitude", "latitude"),
                              crs = 4326)
  message("\n=== Tanzania cluster-buffer GEE extraction (", YEAR, ") ===")
  buffers <- extract_gee_cluster_buffers(clusters_sf, YEAR, id_col = "cnum",
                                         radii_km = c(10, 25, 50), verbose = TRUE,
                                         cache_dir = file.path(CACHE_DIR, "buffers"),
                                         force = FORCE)
  out_buf <- here(sprintf("data/GEE/TZ2010_buffers_%s.csv",
                          format(Sys.Date(), "%m.%d.%Y")))
  readr::write_csv(buffers, out_buf)
  message(sprintf("  wrote %s (%d clusters x %d cols)", out_buf,
                  nrow(buffers), ncol(buffers) - 1L))
} else {
  warning("Tanzania_GMS_cleaned.rds not found — run 1_GW_Tanzania_data_clean.R ",
          "first to get cluster GPS for the buffer product.")
}

message("\nDone. Tanzania GEE CSVs written to data/GEE/.")
