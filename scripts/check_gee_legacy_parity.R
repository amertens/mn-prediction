# =============================================================================
# scripts/check_gee_legacy_parity.R
#
# Reports the shared admin-2 GEE predictor vocabulary across countries and how
# much of it each country actually supplies.
#
# The pooled and LOCO analyses take a strict INTERSECTION of covariate names
# across countries, so one country that names its covariates differently (or is
# missing a domain) silently deletes those covariates for everyone. This script
# makes that visible before a pipeline run rather than after.
#
# Run (reads the {targets} store, so the gee_admin2_* targets must exist):
#   Rscript scripts/check_gee_legacy_parity.R
#
# Writes metadata/gee_legacy_common_vars.txt — the reference vocabulary that
# src/GEE/extract_gee_legacy_parity.R is built to reproduce for a new country.
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("R", "config.R"))
source(here("src", "GEE", "extract_gee_legacy_parity.R"))

countries <- names(get_country_configs())
vars <- list()
for (cn in countries) {
  tname <- paste0("gee_admin2_", tolower(cn))
  x <- tryCatch(targets::tar_read_raw(tname), error = function(e) NULL)
  if (is.null(x)) { message(sprintf("%-14s -- target %s not in store", cn, tname)); next }
  vars[[cn]] <- setdiff(names(x), c("Admin1", "Admin2"))
  message(sprintf("%-14s %4d gee_ variables, %d areas", cn, length(vars[[cn]]), nrow(x)))
}
if (length(vars) < 2) stop("Need at least two countries in the store.")

common <- sort(Reduce(intersect, vars))
message(sprintf("\nIntersection across %d countries: %d variables",
                length(vars), length(common)))

# What each country costs the shared vocabulary by being included.
if (length(vars) > 2) {
  message("\nCost of including each country (variables lost from the intersection):")
  for (cn in names(vars)) {
    without <- length(Reduce(intersect, vars[setdiff(names(vars), cn)]))
    message(sprintf("  %-14s %4d -> %4d  (%+d)", cn, without, length(common),
                    length(common) - without))
  }
}

# Reference vocabulary the legacy-parity extractor targets.
ref_path <- here("metadata", "gee_legacy_common_vars.txt")
dir.create(dirname(ref_path), showWarnings = FALSE, recursive = TRUE)
writeLines(common, ref_path)
message("\nwrote ", ref_path)

spec_cols <- gee_legacy_parity_colnames()
message(sprintf("\nParity extractor covers %d of the %d shared variables.",
                length(intersect(spec_cols, common)), length(common)))
gap <- setdiff(common, spec_cols)
if (length(gap)) {
  message("Not reproducible from Earth Engine (excluded after validation):")
  for (g in gap) message("  ", g)
  message("\nReasons:")
  for (fam in names(GEE_PARITY_EXCLUDED))
    message(sprintf("  %-16s %s", fam, GEE_PARITY_EXCLUDED[[fam]]))
}
