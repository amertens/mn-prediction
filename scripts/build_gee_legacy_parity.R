# =============================================================================
# scripts/build_gee_legacy_parity.R
#
# Extract a country's admin-2 GEE covariates from the Earth Engine API under the
# LEGACY column names, so a country without local .tif exports in
# data/<Country>_GEE_rasters/ can still join the shared cross-country predictor
# vocabulary used by every pooled / LOCO analysis.
#
# Writes data/GEE/<ISO3>_legacy_parity_admin2_gee.csv, which
# extract_gee_admin2() / extract_area_covariates() pick up automatically when the
# country has no raster directory (see .append_legacy_parity_cols in
# R/admin2_analysis.R).
#
# Run:
#   Rscript scripts/build_gee_legacy_parity.R Tanzania
#   Rscript scripts/build_gee_legacy_parity.R SierraLeone --validate
#
# --validate runs the extraction for a country that DOES have legacy rasters and
# reports the per-variable correlation between the EE values and the stored
# raster-derived values, which is how we check that a new country's covariates
# are on the same scale as the existing ones. Sierra Leone (14 districts) is the
# cheap choice.
#
# Prerequisite: Earth Engine access (src/GEE/README_GEE_API_SETUP.md). Verify
# with src/GEE/00_ee_connect_test.R first.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(sf); library(readr)
})

args     <- commandArgs(trailingOnly = TRUE)
country  <- if (length(args) >= 1) args[1] else "Tanzania"
validate <- "--validate" %in% args
force    <- "--force" %in% args

# Tuning knobs. Do NOT raise --scale-floor above the 250 m default without
# re-validating: measured on Sierra Leone, a 1000 m floor still reproduces every
# depth-MEAN column but biases the soil stdev bands and WSF 10-20% high (135/135
# columns pass at 250 m, only 117/135 at 1000 m), because variance and mask
# quantities do not survive coarse aggregation. --batch above ~40 features
# exceeds Earth Engine's 10 MB request limit on detailed admin-2 geometry.
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) default else as.numeric(args[i + 1L])
}
scale_floor <- argval("--scale-floor", 250)
batch       <- argval("--batch", 40)

source(here("R", "config.R"))
source(here("R", "data_prep.R"))                 # load_gadm_cached()
source(here("R", "admin2_analysis.R"))           # legacy_parity_csv_path()
source(here("src", "GEE", "extract_gee_legacy_parity.R"))

cfg <- get_country_configs()[[country]]
if (is.null(cfg)) stop("Unknown country: ", country,
                       " (known: ", paste(names(get_country_configs()), collapse = ", "), ")")

# ---- Earth Engine init (same path as src/Tanzania/4_TZ_GEE_extract.R) --------
rgee_py <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (file.exists(rgee_py)) Sys.setenv(RETICULATE_PYTHON = rgee_py)
ee_project <- Sys.getenv("EE_PROJECT", unset = "mn-prediction-420517")
reticulate::py_run_string(sprintf("import ee; ee.Initialize(project='%s')", ee_project))
invisible(rgee::ee$Number(1)$getInfo())
message("Earth Engine initialised (project ", ee_project, ")")

cache_dir <- here("data", "GEE", paste0(".cache_parity_", cfg$gadm_code))
message("\n=== legacy-parity admin-2 GEE extraction: ", country,
        " (", cfg$gadm_code, ") ===")
t0  <- Sys.time()
message(sprintf("(scale_floor = %d m, batch = %d features/request)",
                scale_floor, batch))
df  <- extract_legacy_parity_admin2(cfg$gadm_code, verbose = TRUE,
                                    cache_dir = cache_dir, force = force,
                                    scale_floor = scale_floor, batch = batch)
message(sprintf("elapsed: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

expected <- gee_legacy_parity_colnames()
gotcols  <- setdiff(names(df), c("Admin1", "Admin2"))
missing  <- setdiff(expected, gotcols)
allna    <- gotcols[vapply(df[gotcols], function(x) all(is.na(x)), logical(1))]
message(sprintf("\n%d/%d expected columns present; %d all-NA",
                length(intersect(expected, gotcols)), length(expected), length(allna)))
if (length(missing)) message("MISSING: ", paste(missing, collapse = ", "))
if (length(allna))   message("ALL-NA : ", paste(allna, collapse = ", "))

out_csv <- legacy_parity_csv_path(cfg)   # keyed on ISO3, see R/admin2_analysis.R
readr::write_csv(df, out_csv)
message(sprintf("wrote %s (%d areas x %d covariates)", out_csv, nrow(df), length(gotcols)))

# ---- Optional validation against the stored raster-derived values ------------
if (validate) {
  message("\n=== validation vs raster-derived gee_admin2_", tolower(country), " ===")
  legacy <- tryCatch(targets::tar_read_raw(paste0("gee_admin2_", tolower(country))),
                     error = function(e) NULL)
  if (is.null(legacy)) {
    message("no stored gee_admin2_ target for ", country, " -- skipping validation")
  } else {
    shared <- intersect(gotcols, names(legacy))
    m <- dplyr::inner_join(
      df[, c("Admin2", shared)],
      legacy[, c("Admin2", shared)],
      by = "Admin2", suffix = c(".ee", ".ras"))
    # Pass rule: the EE column must both RANK areas like the raster one
    # (r >= 0.9) and sit on the same SCALE (mean within 10%). Rank agreement
    # alone is not enough — a constant scale offset on one country is read by
    # the pooled model as a real country difference.
    res <- do.call(rbind, lapply(shared, function(v) {
      a <- m[[paste0(v, ".ee")]]; b <- m[[paste0(v, ".ras")]]
      ok <- is.finite(a) & is.finite(b)
      r <- if (sum(ok) > 2 && stats::sd(a[ok]) > 0 && stats::sd(b[ok]) > 0)
             round(stats::cor(a[ok], b[ok]), 3) else NA_real_
      ratio <- round(mean(a[ok]) / mean(b[ok]), 3)
      data.frame(variable = v, n = sum(ok), r = r,
                 mean_ee = round(mean(a[ok]), 4), mean_ras = round(mean(b[ok]), 4),
                 ratio = ratio,
                 pass = isTRUE(r >= 0.9) && isTRUE(ratio > 0.9) && isTRUE(ratio < 1.1))
    }))
    out_val <- here("results", "sensitivity",
                    paste0("gee_legacy_parity_validation_", tolower(country), ".csv"))
    dir.create(dirname(out_val), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(res, out_val)
    message(sprintf("%d shared variables; median r = %.3f; %d PASS (r>=0.9 and mean within 10%%)",
                    nrow(res), stats::median(res$r, na.rm = TRUE), sum(res$pass)))
    if (any(!res$pass)) {
      message("\nFAILING variables (do not ship these for a new country):")
      print(res[!res$pass, c("variable", "n", "r", "mean_ee", "mean_ras", "ratio")],
            row.names = FALSE)
    }
    message("wrote ", out_val)
  }
}

message("\nDone.")
