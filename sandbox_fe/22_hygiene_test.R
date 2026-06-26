#!/usr/bin/env Rscript
suppressMessages({library(here);library(dplyr)})
source(here::here("R","config.R")); source(here::here("R","data_prep.R"))
cfgs <- get_country_configs()

for (cn in c("Ghana","Gambia")) {
  cc <- cfgs[[cn]]; merged <- load_merged_data(cc$data_path); oc <- cc$outcomes[["child_iron"]]
  if (is.null(oc)) oc <- cc$outcomes[["women_iron"]]
  on  <- build_outcome_dataset(merged, cc, oc, clean_predictors = TRUE)
  off <- build_outcome_dataset(merged, cc, oc, clean_predictors = FALSE)
  dropped <- setdiff(off$Xvars, on$Xvars)
  cat(sprintf("\n=== %s %s ===\n", cn, oc$tag))
  cat(sprintf("  Xvars: clean=%d  raw=%d  dropped=%d\n", length(on$Xvars), length(off$Xvars), length(dropped)))
  cat(sprintf("  dropped MAP _Count: %d | MAP dup-snapshot: %d | IHME counts: %d\n",
              sum(grepl("^MAP_.*_Count$", dropped)),
              sum(grepl("^MAP_[0-9]{6}_", dropped) & !grepl("_Count$", dropped)),
              sum(grepl("^ihme_", dropped))))
  # confirm a Rate sibling survived for each dropped MAP _Count
  mc <- grep("^MAP_.*_Count$", dropped, value=TRUE)
  ok_rate <- all(sapply(mc, function(x) sub("_Count$","_Rate",x) %in% on$Xvars))
  cat(sprintf("  every dropped MAP _Count has surviving _Rate sibling: %s\n", ok_rate))
  # MAP snapshot years remaining
  maprem <- grep("^MAP_[0-9]{6}_", on$Xvars, value=TRUE)
  cat("  MAP snapshot years remaining:", paste(sort(unique(substr(sub("^MAP_","",maprem),1,6))),collapse=", "),
      sprintf(" (survey_year=%d)\n", cc$survey_year))
  cat(sprintf("  bundle still valid: Xvars_bundle=%d (>=2: %s)\n",
              length(on$Xvars_bundle), length(on$Xvars_bundle)>=2))
}
