#!/usr/bin/env Rscript
# Examine shared proxy column names across countries for a given outcome,
# and domain composition of the shared set.
suppressMessages(library(here))
cache_dir <- here::here("sandbox_fe", "cache")
otag <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(otag)) otag <- "child_vitA"
countries <- c("Gambia","Ghana","SierraLeone","Malawi")
objs <- list()
for (cn in countries) {
  fn <- file.path(cache_dir, sprintf("%s__%s.rds", cn, otag))
  if (file.exists(fn)) objs[[cn]] <- readRDS(fn)
}
cat("Outcome:", otag, "| countries:", paste(names(objs), collapse=", "), "\n\n")
colsets <- lapply(objs, function(o) colnames(o$X))
for (cn in names(objs)) cat(sprintf("  %-12s %d cols\n", cn, length(colsets[[cn]])))
shared <- Reduce(intersect, colsets)
cat(sprintf("\nShared across ALL %d countries: %d cols\n", length(objs), length(shared)))

dom_of <- function(cn) {
  pref <- c(DHS="^dhs", MICS="^mics_", IHME="^ihme_", LSMS="^lsms_",
            MAP="^MAP_", WFP="^wfp_", FLUNET="^flunet_", GEE="^gee_",
            FSEC="^fsec_", SOIL="^soil_")
  out <- rep("OTHER", length(cn))
  for (nm in names(pref)) out[grepl(pref[[nm]], cn)] <- nm
  out
}
cat("\nDomain composition of SHARED set:\n")
print(table(dom_of(shared)))

# Also report missingness of shared cols per country (mean fraction NA)
cat("\nMean NA fraction (shared cols) per country:\n")
for (cn in names(objs)) {
  X <- objs[[cn]]$X[, shared, drop=FALSE]
  cat(sprintf("  %-12s %.3f\n", cn, mean(is.na(as.matrix(sapply(X, as.numeric))))))
}
saveRDS(shared, file.path(cache_dir, sprintf("_shared_%s.rds", otag)))
