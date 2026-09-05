# =============================================================================
# scripts/protocol_v2/42_alphaearth_as_domain.R   [AE-01, prep]
#
# ALPHAEARTH EMBEDDINGS ARE ALREADY IN THE VOCABULARY -- FILED UNDER AGRICULTURE
#
# The 64 AlphaEarth satellite-embedding dimensions (Google's 2025 foundation
# embedding, one 64-vector per district) sit in predictors_admin2_shared.csv
# as aef_A00..aef_A63 with domain "Agricultural production, land use" and
# source "MapSPAM / AEF". They are two thirds of that 93-column domain. So the
# domain ablation (DA-01) that found agriculture dead weight under transport
# was mostly ablating a satellite embedding, and the agriculture PCs are
# mostly embedding PCs. This prep script copies the 64 columns to a new file
# under new names (sat_A00..) so script 39 can score them as their OWN domain
# (ADDON_DROP='^aef_' gives the 'replace' set: embedding out of agriculture and
# into a block of its own; 'addon_only' is the embedding alone).
#
#   Rscript scripts/protocol_v2/42_alphaearth_as_domain.R
# -> data/covariates/harmonized/predictors_admin2_alphaearth.csv
# =============================================================================
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
MD <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"))
aef <- grep("^aef_A[0-9]{2}$", names(S), value = TRUE)
cat("aef columns in shared:", length(aef), "| metadata domain:", paste(unique(MD$domain[MD$column %in% aef]), collapse = "; "), "| source:", paste(unique(MD$source[MD$column %in% aef]), collapse = "; "), "\n")
cat("agriculture domain size:", sum(MD$domain == "Agricultural production, land use"), "of which aef:", sum(MD$column[MD$domain == "Agricultural production, land use"] %in% aef), "\n")
out <- S[, c("country", "Admin1", "Admin2", aef)]; names(out)[-(1:3)] <- sub("^aef_", "sat_", aef)
cat("coverage by country:", paste(sprintf("%s %.2f", unique(out$country), tapply(is.finite(out$sat_A00), out$country, mean)[unique(out$country)]), collapse = " | "), "\n")
write.csv(out, file.path(HDIR, "predictors_admin2_alphaearth.csv"), row.names = FALSE)
cat("DONE\n")
