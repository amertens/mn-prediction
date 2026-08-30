# =============================================================================
# scripts/covariates/01_extract_alphaearth.R                        [STAGE 1]
#
# EXPENSIVE DOWNLOAD — RUN ONCE. Admin-2 zonal means of the AlphaEarth
# Foundations satellite embedding (GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL):
# 64 dimensions (A00..A63), 10 m native, annual, 2017-2025.
#
# Output: data/AlphaEarth/<CODE>_alphaearth_admin2.csv
#         columns Admin1, Admin2, gee_alphaearth_A00 .. _A63
#
# WHY ONE FIXED YEAR (2017) FOR EVERY COUNTRY
# -------------------------------------------
# AEF starts in 2017; the surveys run 2010-2018. Using each country's nearest
# available year would make embedding vintage collinear with country -- a
# country fixed effect wearing a covariate label, which is exactly the failure
# mode GEE_PARITY_EXCLUDED exists to prevent. One shared vintage keeps the
# cross-country comparison clean and pushes the temporal mismatch into a single
# documented caveat instead of a hidden confound.
#
# Measured drift over Sierra Leone's 14 districts, 2017 vs 2024 (a 7-year gap,
# the Tanzania worst case): between-district Pearson r 0.70-0.97 by dimension.
# Structural signal is largely stable; treat AEF as a time-invariant land-
# surface descriptor, not as a covariate measured at survey time.
#
# NORMALISATION: AEF pixel vectors are unit-norm. A polygon mean is not, so the
# aggregate is re-normalised to the unit hypersphere (see AEF_L2_RENORM). Per-
# dimension z-scoring is deliberately NOT applied -- it destroys the angular
# geometry the embedding was trained on.
#
#   Rscript scripts/covariates/01_extract_alphaearth.R            # all countries
#   Rscript scripts/covariates/01_extract_alphaearth.R Ghana      # one
#   AEF_FORCE=1 Rscript ...                                       # ignore cache
# =============================================================================
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(here)
})
source(here("R", "data_prep.R"))            # load_gadm_cached()
source(here("R", "admin2_key_hygiene.R"))   # clean_admin2_keys()
source(here("src", "GEE", "extract_gee_ee_api.R"))  # extract_gee_layers()

AEF_YEAR      <- as.integer(Sys.getenv("AEF_YEAR", "2017"))
AEF_L2_RENORM <- !identical(tolower(Sys.getenv("AEF_L2_RENORM", "true")), "false")
FORCE         <- nzchar(Sys.getenv("AEF_FORCE", ""))

COUNTRIES <- list(
  Gambia      = list(gadm = "GMB"),
  Ghana       = list(gadm = "GHA"),
  SierraLeone = list(gadm = "SLE"),
  Malawi      = list(gadm = "MWI"),
  Tanzania    = list(gadm = "TZA")
)

# One-entry manifest. yoy = FALSE on purpose: a year-over-year embedding delta
# would double the block to 128 columns for a change signal we cannot interpret.
AEF_MANIFEST <- list(list(
  family = "AlphaEarth", asset = "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL",
  kind = "ImageCollection", bands = NULL, temporal = "filter_year",
  yoy = FALSE, scale = 10, avail = c(2017L, 2025L),
  note = "AlphaEarth Foundations annual embedding field; 64 dims, unit-norm per pixel."
))

ee_project <- Sys.getenv("EE_PROJECT", unset = "mn-prediction-420517")
rgee_py <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (file.exists(rgee_py)) Sys.setenv(RETICULATE_PYTHON = rgee_py)
reticulate::py_run_string(sprintf("import ee; ee.Initialize(project='%s')", ee_project))
invisible(rgee::ee$Number(1)$getInfo())
message("Earth Engine initialised (project ", ee_project, ")")

out_dir <- here("data", "AlphaEarth"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
sel  <- if (length(args)) args else names(COUNTRIES)

for (cn in sel) {
  cc <- COUNTRIES[[cn]]; if (is.null(cc)) { message("unknown country: ", cn); next }
  out <- file.path(out_dir, sprintf("%s_alphaearth_admin2.csv", cc$gadm))
  if (file.exists(out) && !FORCE) { message(cn, ": cached -> ", basename(out)); next }

  message("\n=== ", cn, " (", cc$gadm, ", year ", AEF_YEAR, ") ===")
  poly <- load_gadm_cached(cc$gadm, level = 2) |>
    sf::st_as_sf() |> dplyr::select(NAME_1, NAME_2) |>
    dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)
  poly$Admin1 <- trimws(as.character(poly$Admin1))
  poly$Admin2 <- trimws(as.character(poly$Admin2))

  res <- tryCatch(
    extract_gee_layers(poly, AEF_YEAR, manifest = AEF_MANIFEST, id_col = "Admin2",
                       col_prefix = "gee_", verbose = TRUE, scale_floor = 250,
                       cache_dir = file.path(out_dir, paste0(".cache_", cc$gadm)),
                       force = FORCE),
    error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  if (is.null(res)) next

  emb <- grep("^gee_AlphaEarth_A[0-9]+$", names(res), value = TRUE)
  if (length(emb) != 64L) {
    message("  expected 64 embedding columns, got ", length(emb), " -- skipping")
    next
  }
  # Legacy-vocabulary naming: lower-case stem, no year token. AEF is used as a
  # time-invariant descriptor, so a year stamp would only re-create the
  # cross-country name mismatch this whole exercise exists to remove.
  names(res)[match(emb, names(res))] <- paste0("gee_alphaearth_", sub("^gee_AlphaEarth_", "", emb))
  emb <- grep("^gee_alphaearth_A[0-9]+$", names(res), value = TRUE)
  emb <- emb[order(emb)]

  if (AEF_L2_RENORM) {
    M <- as.matrix(res[, emb]); nrm <- sqrt(rowSums(M^2, na.rm = TRUE))
    nrm[!is.finite(nrm) | nrm == 0] <- NA_real_
    res[, emb] <- M / nrm
    message(sprintf("  L2-renormalised (pre-norm range %.3f-%.3f)",
                    min(nrm, na.rm = TRUE), max(nrm, na.rm = TRUE)))
  }

  res <- clean_admin2_keys(res, what = paste(cn, "AlphaEarth"))
  keep <- c(intersect(c("Admin1", "Admin2"), names(res)), emb)
  readr::write_csv(res[, keep, drop = FALSE], out)
  message(sprintf("  wrote %s (%d areas x 64 dims)", basename(out), nrow(res)))
}
message("\nStage 1 complete.")
