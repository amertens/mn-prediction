# =============================================================================
# scripts/validate_assay_lineage.R
#
# WS3, deliverable 2. Fail if any modeled outcome lacks an assay-lineage row.
#
# The point is that a country-outcome cannot enter a model without a recorded
# provenance row. UNKNOWN is an acceptable VALUE in a lineage row; a missing
# ROW is not, because it means an outcome is being modelled that nobody has
# written down the measurement basis for.
#
# Exits 1 on failure so it can gate a pipeline run, and is also wired as a
# pipeline target with an md5 stamp on the CSV so an edited lineage file
# re-triggers validation.
#
# Run:
#   Rscript scripts/validate_assay_lineage.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
targets::tar_source(here("R"))

# The validator itself lives in R/assay_lineage.R because a pipeline target
# calls it and tar_source() loads R/ but not scripts/.

res <- validate_assay_lineage()

cat(sprintf("assay lineage: %d rows for %s required country-outcomes\n",
            res$n_rows, res$n_required))
if (length(res$extra))
  cat("rows not matching a configured outcome (retained, not an error):\n  ",
      paste(res$extra, collapse = "\n  "), "\n", sep = "")

if (!res$ok) {
  if (length(res$missing))
    cat("\nMISSING lineage rows for:\n  ", paste(res$missing, collapse = "\n  "), "\n", sep = "")
  if (length(res$blank_sources))
    cat("\nBLANK source cells in: ",
        paste(sprintf("%s (%d)", names(res$blank_sources), res$blank_sources),
              collapse = ", "), "\n", sep = "")
  cat("\nFAIL: every modeled outcome needs a lineage row with a source for every asserted cell.\n")
  quit(status = 1L)
}

cat("PASS: every configured country-outcome has a lineage row, and every source cell is populated.\n")
