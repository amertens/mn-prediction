# =============================================================================
# sandbox_parsimony/R/07_join_diagnostics.R
#
# Is the Admin-2 key joining survey prevalence to the right polygon?
#
# Two spatial-autocorrelation checks per country:
#   moran_cov  -- Moran I of the GEE covariates themselves. Raster aggregates
#                 over adjacent districts MUST be spatially smooth. A value
#                 near 0 means the covariate table is attached to the wrong
#                 Admin-2 rows (or the centroids are).
#   moran_out  -- Moran I of the survey prevalence. Low is not by itself
#                 damning (deficiency need not be spatially smooth) but low
#                 outcome + high covariate autocorrelation + negative CV skill
#                 across every model is the signature of a broken outcome join.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")

#' Moran I with k-nearest-neighbour weights built from Admin-2 centroids
moran_knn <- function(x, lon, lat, k = 5L) {
  ok <- is.finite(x) & is.finite(lon) & is.finite(lat)
  if (sum(ok) < 15) return(NA_real_)
  x <- x[ok]; co <- cbind(lon[ok], lat[ok]); n <- length(x)
  D <- as.matrix(stats::dist(co))
  diag(D) <- Inf
  W <- matrix(0, n, n)
  kk <- min(k, n - 1)
  for (i in seq_len(n)) W[i, order(D[i, ])[seq_len(kk)]] <- 1 / kk
  z <- x - mean(x)
  num <- sum(W * outer(z, z))
  den <- sum(z^2)
  if (den <= 0) return(NA_real_)
  (n / sum(W)) * num / den
}

rows <- list()
for (oc in names(pooled_all)) {
  P <- pooled_all[[oc]]
  d <- P$data
  # a handful of covariates that are unambiguously smooth in space
  probe <- intersect(c("gee_elevation", "gee_ndvi", "gee_popdensity", "gee_trmm",
                       "gee_soiltotalcarbon_mean_0_20"), P$predictors)
  for (ctry in unique(d$country)) {
    s <- d[d$country == ctry, ]
    s <- s[is.finite(s$svy_prev), ]
    if (nrow(s) < 15) next
    mcov <- if (length(probe))
      mean(vapply(probe, function(v) moran_knn(s[[v]], s$lon, s$lat), numeric(1)),
           na.rm = TRUE) else NA_real_
    rows[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, n = nrow(s),
      pct_coords_missing = round(100 * mean(!is.finite(s$lon)), 0),
      moran_cov = round(mcov, 3),
      moran_out = round(moran_knn(s$svy_prev, s$lon, s$lat), 3),
      stringsAsFactors = FALSE)
  }
}
res <- bind_rows(rows) |> arrange(country, outcome)
cat("\n=== Spatial-autocorrelation join check (Moran I, 5-NN) ===\n")
cat("moran_cov near 0 => the covariate table is not aligned to the polygons.\n\n")
print(as.data.frame(res), row.names = FALSE)
write.csv(res, "sandbox_parsimony/out/join_diagnostics.csv", row.names = FALSE)

# ---- name-match audit: how many survey Admin2 labels find a covariate row? --
cat("\n=== Admin-2 key match rates (survey rows vs covariate table) ===\n")
cov_list <- readRDS("sandbox_parsimony/out/cov_list.rds")
STORE <- "_targets_full/objects"
mr <- list()
for (oc in names(pooled_all)) {
  for (ctry in names(cov_list)) {
    f <- file.path(STORE, paste0("svy_admin2_", ctry, "_", oc))
    if (!file.exists(f)) next
    s <- readRDS(f)
    a2 <- as.character(s$Admin2)
    cv <- as.character(cov_list[[ctry]]$Admin2)
    mr[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry,
      n_survey_areas = length(a2), n_cov_areas = length(cv),
      pct_survey_matched = round(100 * mean(a2 %in% cv), 0),
      n_unmatched = sum(!(a2 %in% cv)),
      example_unmatched = paste(head(a2[!(a2 %in% cv)], 3), collapse = " | "),
      stringsAsFactors = FALSE)
  }
}
mr <- bind_rows(mr) |> arrange(country, outcome)
print(as.data.frame(mr), row.names = FALSE)
write.csv(mr, "sandbox_parsimony/out/key_match_rates.csv", row.names = FALSE)
