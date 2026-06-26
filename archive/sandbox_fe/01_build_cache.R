#!/usr/bin/env Rscript
# =============================================================================
# 01_build_cache.R
# Build compact per-(country x outcome) caches of PROXY-ONLY predictors plus
# outcome / cluster / Admin2, so feature-engineering experiments load fast.
#
# Focus outcomes (present across multiple countries -> usable for LOCO transfer):
#   women_iron, child_vitA, child_iron, women_vitA
# =============================================================================
suppressMessages({
  library(here)
  library(dplyr)
})
source(here::here("R", "config.R"))
source(here::here("R", "data_prep.R"))

cfgs <- get_country_configs()
focus_outcomes <- c("women_iron", "child_vitA", "child_iron", "women_vitA",
                    "women_b12", "women_folate")

out_dir <- here::here("sandbox_fe", "cache")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

manifest <- list()

for (country in names(cfgs)) {
  cc <- cfgs[[country]]
  if (!file.exists(cc$data_path)) { cat("skip (no data):", country, "\n"); next }
  cat("\n##### Loading", country, "#####\n")
  merged <- tryCatch(load_merged_data(cc$data_path),
                     error = function(e) { cat("  LOAD ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(merged)) next

  for (otag in focus_outcomes) {
    oc <- cc$outcomes[[otag]]
    if (is.null(oc)) next
    od <- tryCatch(build_outcome_dataset(merged, cc, oc),
                   error = function(e) { cat("  build error", otag, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(od)) next
    d <- od$data
    if (nrow(d) == 0 || length(od$Xvars) == 0) next

    # Assemble compact object: outcome cols, cluster, admin2, proxy X
    y_cont <- if (!is.null(oc$continuous) && oc$continuous %in% colnames(d)) d[[oc$continuous]] else NA
    y_bin  <- if (!is.null(oc$binary) && oc$binary %in% colnames(d)) d[[oc$binary]] else NA
    admin2 <- if (cc$admin2_col %in% colnames(d)) d[[cc$admin2_col]] else NA
    clust  <- if (cc$cluster_id %in% colnames(d)) d[[cc$cluster_id]] else seq_len(nrow(d))

    X <- d[, intersect(od$Xvars, colnames(d)), drop = FALSE]
    # tag each proxy column by domain prefix
    dom_of <- function(cn) {
      pref <- c(DHS="^dhs", MICS="^mics_", IHME="^ihme_", LSMS="^lsms_",
                MAP="^MAP_", WFP="^wfp_", FLUNET="^flunet_", GEE="^gee_",
                FSEC="^fsec_", SOIL="^soil_", CHIRPS="^chirps_", WORLDPOP="^worldpop_",
                MAP2="^map2_", NTL="^ntl_", GDL="^gdl_", CROP="^crop_", IPC="^ipc_", ACLED="^acled_")
      out <- rep("OTHER", length(cn))
      for (nm in names(pref)) out[grepl(pref[[nm]], cn)] <- nm
      out
    }
    domains <- dom_of(colnames(X))

    obj <- list(
      country = country, outcome = otag,
      y_cont = as.numeric(y_cont), y_bin = as.numeric(y_bin),
      admin2 = as.character(admin2), cluster = as.character(clust),
      X = X, domains = domains,
      cutoff = oc$cutoff, cutoff_scale = oc$cutoff_scale
    )
    fn <- file.path(out_dir, sprintf("%s__%s.rds", country, otag))
    saveRDS(obj, fn)
    n_bin_pos <- sum(obj$y_bin == 1, na.rm = TRUE)
    n_bin_tot <- sum(!is.na(obj$y_bin))
    manifest[[length(manifest)+1]] <- data.frame(
      country=country, outcome=otag, n=nrow(X), p=ncol(X),
      n_bin_pos=n_bin_pos, n_bin_tot=n_bin_tot,
      prev=ifelse(n_bin_tot>0, round(n_bin_pos/n_bin_tot,3), NA))
    cat(sprintf("  saved %s__%s: n=%d p=%d prev=%.3f\n", country, otag,
                nrow(X), ncol(X), ifelse(n_bin_tot>0, n_bin_pos/n_bin_tot, NA)))
  }
  rm(merged); gc()
}

man <- do.call(rbind, manifest)
write.csv(man, file.path(out_dir, "_manifest.csv"), row.names = FALSE)
cat("\n=== MANIFEST ===\n"); print(man)
