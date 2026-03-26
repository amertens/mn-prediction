# =============================================================================
# scripts/download_gadm.R
#
# Download GADM admin boundaries for all countries to a local cache.
# Run this ONCE with internet access. All pipeline code then uses the cache.
# =============================================================================

library(geodata)

gadm_dir <- here::here("data", "gadm")
dir.create(gadm_dir, recursive = TRUE, showWarnings = FALSE)

countries <- list(
  GMB = "Gambia",
  GHA = "Ghana",
  SLE = "Sierra Leone",
  MWI = "Malawi"
)

for (code in names(countries)) {
  for (level in 1:2) {
    cat(sprintf("Downloading GADM %s level %d (%s)...\n", code, level, countries[[code]]))
    tryCatch({
      geodata::gadm(code, level = level, path = gadm_dir)
      cat("  OK\n")
    }, error = function(e) {
      cat(sprintf("  FAILED: %s\n", e$message))
    })
  }
}

cat("\nCached files:\n")
print(list.files(gadm_dir, pattern = "\\.rds$"))
cat(sprintf("\nAll GADM data cached in: %s\n", gadm_dir))
cat("Pipeline will use these without internet access.\n")
