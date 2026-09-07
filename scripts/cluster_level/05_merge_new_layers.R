# =============================================================================
# scripts/cluster_level/05_merge_new_layers.R
#
# Add the 2026-09-07 layers to the cluster predictor table so the cluster
# track sees the same vocabulary additions as the district track:
#   data/covariates/cluster/predictors_cluster_ihme_raster.csv    IHME 5 km surfaces at the buffer (17 validated columns)
#   data/covariates/cluster/predictors_cluster_livestock.csv      GLW4 head/km2, TLU/km2, TLU per person, ruminant share
#   data/covariates/cluster/predictors_cluster_water_distance.csv distance to permanent / any surface water and to the coast
#   data/covariates/cluster/predictors_cluster_espen.csv          ESPEN endemicity class and MDA history (district value)
# Columns already present (the cluster extraction's own ihme_* CGF / anaemia
# columns) are replaced by the new value; metadata rows are added with the
# same domain labels the Admin-2 set uses. Idempotent: re-running replaces.
#
#   Rscript scripts/cluster_level/05_merge_new_layers.R
# -> data/covariates/cluster/predictors_cluster.csv (+ _metadata.csv), backups *.pre_new_layers
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
CDIR <- "data/covariates/cluster"; OT <- Sys.getenv("CL_OUT_TAG", "")
P  <- read.csv(file.path(CDIR, paste0("predictors_cluster", OT, ".csv")), check.names = FALSE)
MD <- read.csv(file.path(CDIR, paste0("predictors_cluster", OT, "_metadata.csv")), check.names = FALSE)
for (f in c(paste0("predictors_cluster", OT, ".csv"), paste0("predictors_cluster", OT, "_metadata.csv")))
  if (!file.exists(file.path(CDIR, paste0(f, ".pre_new_layers")))) file.copy(file.path(CDIR, f), file.path(CDIR, paste0(f, ".pre_new_layers")))
key <- c("country", "cluster")
blocks <- list(
  list(file = "predictors_cluster_ihme_raster.csv", domain_fn = function(v) ifelse(grepl("anemia|stunting|wasting|underweight", v), "Nutrition status (MODELLED SURFACE)", "Infant and child morbidity/mortality"), role = "slow", source = "IHME (5 km surfaces, buffer mean)", validated_only = TRUE),
  list(file = "predictors_cluster_livestock.csv", domain_fn = function(v) "Livestock density", role = "static", source = "GLW4 2020"),
  list(file = "predictors_cluster_water_distance.csv", domain_fn = function(v) "Water and coast proximity", role = "static", source = "JRC GSW / LSIB (Earth Engine)"),
  list(file = "predictors_cluster_espen.csv", domain_fn = function(v) "Helminth burden and control", role = "slow", source = "WHO ESPEN (district value)"))
for (b in blocks) {
  f <- file.path(CDIR, sub("[.]csv$", paste0(OT, ".csv"), b$file)); if (!file.exists(f)) f <- file.path(CDIR, b$file)   # ESPEN is a district value: untagged
  if (!file.exists(f)) { cat("  missing", b$file, "\n"); next }
  X <- read.csv(f, check.names = FALSE); X <- X[, intersect(names(X), c(key, setdiff(names(X), c("Admin1", "Admin2", "radius_km")))), drop = FALSE]
  if (isTRUE(b$validated_only)) { rm <- read.csv("data/covariates/harmonized/predictors_admin2_ihme_raster_metadata.csv"); ok <- unique(rm$column[is.finite(rm$rho_all) & rm$rho_all >= 0.5]); X <- X[, c(key, intersect(ok, names(X))), drop = FALSE] }
  vals <- setdiff(names(X), key); X <- X |> distinct(across(all_of(key)), .keep_all = TRUE)
  P <- P |> select(-any_of(vals)) |> left_join(X, by = key)
  MD <- MD[!MD$column %in% vals, ]
  cov <- sapply(vals, function(v) mean(is.finite(P[[v]]))); ctry <- sapply(vals, function(v) paste(sort(unique(P$country[is.finite(P[[v]])])), collapse = "|"))
  MD <- bind_rows(MD, data.frame(column = vals, domain = b$domain_fn(vals), role = b$role, source = b$source, completeness = round(cov, 3), countries = ctry, n_countries = sapply(strsplit(ctry, "[|]"), length), stringsAsFactors = FALSE))
  cat(sprintf("  %-40s %2d columns | mean completeness %.2f\n", b$file, length(vals), mean(cov)))
}
write.csv(P, file.path(CDIR, paste0("predictors_cluster", OT, ".csv")), row.names = FALSE); write.csv(MD, file.path(CDIR, paste0("predictors_cluster", OT, "_metadata.csv")), row.names = FALSE)
cat(sprintf("cluster predictors: %d clusters x %d columns; metadata %d rows; domains: %s\n", nrow(P), ncol(P), nrow(MD), paste(sort(unique(MD$domain)), collapse = " | ")))
cat("DONE\n")
