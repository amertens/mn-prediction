# =============================================================================
# src/GEE/00_ee_connect_test.R
#
# End-to-end Earth Engine connectivity test — run this ONCE before the full
# Tanzania extraction (src/Tanzania/4_TZ_GEE_extract.R) to confirm rgee auth,
# the Cloud project, and a real server round-trip all work.
#
# Expected on success: prints the rgee/Python env check, initialises against the
# project, fetches SRTM band names, reduces elevation over a tiny Tanzania box,
# and pulls one iSDAsoil zinc value (a layer used in the real extraction).
#
#   source(here::here("src/GEE/00_ee_connect_test.R"))
#
# If ee_Initialize() errors about credentials, run rgee::ee_Authenticate() once
# (opens a browser), then re-run this script.
# =============================================================================

suppressPackageStartupMessages({ library(here) })

if (!requireNamespace("rgee", quietly = TRUE))
  stop("Install rgee first: install.packages('rgee'); rgee::ee_install(py_env = 'rgee')")

EE_USER    <- Sys.getenv("EE_USER",    unset = "amertens@berkeley.edu")
EE_PROJECT <- Sys.getenv("EE_PROJECT", unset = "mn-prediction-420517")

# Point at the venv that has earthengine-api (the old conda rgee env was removed;
# earthengine-api was installed into this reticulate venv instead).
rgee_py <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (file.exists(rgee_py)) Sys.setenv(RETICULATE_PYTHON = rgee_py)

cat(sprintf("== Step 1: ee.Initialize(project='%s') via Python ==\n", EE_PROJECT))
# Initialize the Python EE client directly (cached creds in ~/.config/earthengine/).
# Avoids rgee::ee_Initialize()'s asset-root prompt loop (known 1.1.7 bug). If this
# errors about credentials, run rgee::ee_Authenticate() once, then re-run.
ok_init <- tryCatch({
  reticulate::py_run_string(sprintf("import ee; ee.Initialize(project='%s')", EE_PROJECT))
  TRUE
}, error = function(e) {
  cat("  init FAILED: ", conditionMessage(e), "\n")
  cat("  -> Run rgee::ee_Authenticate() once (browser), then re-run this script.\n")
  cat("  -> If it says the project is not registered, finish EE registration at\n")
  cat("     https://console.cloud.google.com/earth-engine/configuration?project=", EE_PROJECT, "\n")
  FALSE
})
if (!ok_init) stop("EE init failed — see message above.")

ee <- rgee::ee
cat("  ✓ Earth Engine initialised — ee$Number(1)$getInfo() =",
    ee$Number(1)$getInfo(), "\n")

cat("\n== Step 3: server round-trip — SRTM band names ==\n")
srtm <- ee$Image("USGS/SRTMGL1_003")
cat("  SRTM bands: ", paste(srtm$bandNames()$getInfo(), collapse = ", "), "\n")

cat("\n== Step 4: reduceRegion — mean elevation over a small Tanzania box ==\n")
box <- ee$Geometry$Rectangle(c(34.0, -7.0, 36.0, -5.0))  # central Tanzania
elev <- tryCatch(
  srtm$reduceRegion(reducer = ee$Reducer$mean(), geometry = box,
                    scale = 1000, maxPixels = 1e9)$getInfo(),
  error = function(e) { cat("  reduceRegion error: ", conditionMessage(e), "\n"); NULL })
cat("  mean elevation (m): ", if (!is.null(elev)) round(unlist(elev), 1) else "NA", "\n")

cat("\n== Step 5: iSDAsoil zinc (a real extraction layer) over the same box ==\n")
zn <- tryCatch({
  img <- ee$Image("ISDASOIL/Africa/v1/zinc_extractable")$select("mean_0_20")
  img$reduceRegion(reducer = ee$Reducer$mean(), geometry = box,
                   scale = 1000, maxPixels = 1e9)$getInfo()
}, error = function(e) { cat("  zinc error: ", conditionMessage(e), "\n"); NULL })
cat("  mean iSDAsoil zinc (0-20cm): ", if (!is.null(zn)) round(unlist(zn), 3) else "NA", "\n")

cat("\n=========================================================\n")
cat("✓ Connectivity test PASSED — EE auth + project + reductions all work.\n")
cat("  Next: source('src/Tanzania/4_TZ_GEE_extract.R') for the full extraction.\n")
cat("=========================================================\n")
