# =============================================================================
# scripts/covariates/run_all.R
#
# Run the whole covariate pipeline, stage 1 -> 4.
#
#   Rscript scripts/covariates/run_all.R              # skips cached downloads
#   Rscript scripts/covariates/run_all.R --from 2     # start at stage 2
#
# Stage 1 is the only expensive step and it is cached per country: re-running
# is cheap and safe. Stages 2-4 are pure transformations of what is already on
# disk, so they can be re-run freely while harmonisation rules are being
# revised -- which is the point of keeping download and harmonisation apart.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
from <- if ("--from" %in% args) as.integer(args[which(args == "--from") + 1L]) else 1L

stages <- c(
  "01_extract_alphaearth.R",
  "02_build_country_predictors.R",
  "03_harmonize.R",
  "04_document_and_qc.R")

root <- here::here("scripts", "covariates")
for (i in seq_along(stages)) {
  if (i < from) { message("skipping stage ", i, " (", stages[i], ")"); next }
  message("\n############ STAGE ", i, ": ", stages[i], " ############")
  src <- file.path(root, stages[i])
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(src))
  if (!identical(status, 0L)) stop("stage ", i, " failed (exit ", status, ")")
}
message("\nAll stages complete. Harmonised outputs in data/covariates/harmonized/")
