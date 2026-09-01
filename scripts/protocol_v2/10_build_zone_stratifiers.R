# =============================================================================
# scripts/protocol_v2/10_build_zone_stratifiers.R
#
# Categorical ZONE stratifiers at Admin-2: Koppen-Geiger climate class and the
# IFPRI/HarvestChoice Agro-Ecological Zones for Sub-Saharan Africa.
#
# WHY, GIVEN 41 CLIMATE AND 93 AGRICULTURE COLUMNS ALREADY
# --------------------------------------------------------
# The existing climate and land-use predictors are all CONTINUOUS surfaces -
# monthly temperature, rainfall, vegetation indices, crop shares. A zone class
# is a different kind of variable: it is the standard stratifier that
# agricultural and nutrition work in this region actually uses, it is stable
# rather than year-specific, and it encodes interactions (hot AND dry AND
# bimodal rainfall) that a linear model over separate continuous columns cannot
# reach. Two derived forms are produced per source:
#
#   *_class      the modal (most common) zone class in the district
#   *_purity     the share of the district covered by that modal class, which
#                distinguishes a district that sits squarely inside one zone
#                from one straddling a boundary
#   *_n_classes  how many distinct classes occur, a simple heterogeneity count
#
# The class itself is an unordered CATEGORY. It is emitted as an integer code
# plus, for the models, a set of one-hot shares for the classes that are common
# enough to be usable (>= 5 percent of districts in at least one country), so
# no learner is misled into treating "class 12" as greater than "class 4".
#
#   Rscript scripts/protocol_v2/10_build_zone_stratifiers.R
# -> data/covariates/harmonized/zone_stratifiers_admin2.csv
# =============================================================================
suppressPackageStartupMessages({library(sf); library(terra); library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
sf::sf_use_s2(FALSE)

HDIR <- "data/covariates/harmonized"
LC <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
        sierraleone = "SierraLeone")
B <- readRDS("dashboard/data/admin2_boundaries.rds")

polys <- do.call(rbind, lapply(names(LC), function(lc) {
  b <- sf::st_transform(sf::st_make_valid(B[[lc]]), 4326)
  d <- sf::st_drop_geometry(b)[, c("Admin1", "Admin2")]
  s <- sf::st_sf(country = LC[[lc]], Admin1 = as.character(d$Admin1),
                 Admin2 = as.character(d$Admin2), geometry = sf::st_geometry(b))
  s[!sf::st_is_empty(s), ]
}))
cat("polygons:", nrow(polys), "\n")

#' Modal class, its purity, and the class count, per polygon
zone_summary <- function(rast_path, prefix) {
  if (!file.exists(rast_path)) {
    cat("  missing:", rast_path, "\n"); return(NULL)
  }
  r <- terra::rast(rast_path)
  v <- terra::vect(polys)
  ex <- terra::extract(r, v)
  names(ex)[2] <- "z"
  # categorical rasters (AEZ ships one) come back as a factor; the class CODE
  # is what we want, and comparing a factor with > is meaningless
  if (is.factor(ex$z)) ex$z <- as.integer(ex$z)
  ex$z <- suppressWarnings(as.numeric(ex$z))
  s <- ex |> group_by(ID) |>
    summarise(
      cls = { t <- table(z[is.finite(z) & z > 0])
              if (!length(t)) NA_real_ else as.numeric(names(t)[which.max(t)]) },
      pur = { t <- table(z[is.finite(z) & z > 0])
              if (!length(t)) NA_real_ else max(t) / sum(t) },
      ncl = { t <- table(z[is.finite(z) & z > 0]); length(t) },
      .groups = "drop")
  out <- data.frame(cls = s$cls, pur = round(s$pur, 4), ncl = s$ncl)
  names(out) <- paste0(prefix, c("_class", "_purity", "_n_classes"))
  out
}

res <- polys |> sf::st_drop_geometry()

kp <- "data/Koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
cat("[Koppen-Geiger 1991-2020]\n")
z1 <- zone_summary(kp, "koppen")
if (!is.null(z1)) res <- cbind(res, z1)

# AEZ ships zipped; extract once into a cache directory
aezdir <- file.path("data", "Agro-Ecological Zones for Africa South of the Sahara")
cache <- file.path(aezdir, "_extracted")
dir.create(cache, showWarnings = FALSE)
zipf <- file.path(aezdir, "AEZ16 r2.0 - TIF.zip")
if (file.exists(zipf) && !length(list.files(cache, "\\.tif$")))
  utils::unzip(zipf, exdir = cache)
aeztif <- list.files(cache, "\\.tif$", full.names = TRUE, recursive = TRUE)[1]
cat("[AEZ16 Sub-Saharan Africa]\n")
if (!is.na(aeztif)) {
  z2 <- zone_summary(aeztif, "aez16")
  if (!is.null(z2)) res <- cbind(res, z2)
} else cat("  no AEZ tif found\n")

# one-hot shares for classes common enough to be usable
for (pfx in c("koppen", "aez16")) {
  cc <- paste0(pfx, "_class")
  if (!cc %in% names(res)) next
  tb <- table(res[[cc]])
  common <- as.numeric(names(tb)[tb >= max(5, 0.05 * nrow(res))])
  for (k in common)
    res[[sprintf("%s_is_%d", pfx, k)]] <- as.integer(res[[cc]] == k)
  cat(sprintf("  %s: %d distinct classes, %d common enough to one-hot\n",
              pfx, length(tb), length(common)))
}

write.csv(res, file.path(HDIR, "zone_stratifiers_admin2.csv"), row.names = FALSE)
cat("\nwrote", file.path(HDIR, "zone_stratifiers_admin2.csv"), "with",
    ncol(res) - 3, "columns\n")
newc <- setdiff(names(res), c("country", "Admin1", "Admin2"))
for (v in head(newc, 12)) {
  s <- tapply(res[[v]], res$country, function(z) mean(z, na.rm = TRUE))
  cat(sprintf("  %-20s complete=%.2f | %s\n", v, mean(is.finite(res[[v]])),
              paste(sprintf("%s=%.2f", names(s), s), collapse = " ")))
}
cat("\nDONE\n")
