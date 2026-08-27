# =============================================================================
# sandbox_parsimony/R/23_travel_time.R
#
# TRAVEL TIME TO CITIES for every country, from one consistent global product.
#
# Why this variable: of everything missing from the five-country intersection
# (FINDINGS.md section 4), travel time to cities has the strongest prior claim on
# dietary diversity -- market access governs whether a household can buy animal-
# source foods, fortified staples and vegetables. It is also the variable whose
# absence is least defensible, since four of the five countries already have a
# version of it sitting in data/<Country>_GEE_rasters/.
#
# WHY NOT REUSE THOSE FILES: they are per-country clips of the Weiss et al.
# (2018) MAP surface, and Tanzania has none. Mixing a MAP-derived column for
# four countries with a different product for the fifth would give the pooled
# and LOCO models a covariate that means slightly different things in different
# countries -- exactly the failure mode FINDINGS.md section 4 documents. So this
# uses ONE product for all five.
#
# SOURCE:   Nelson et al. (2019), travel time to cities, via
#           geodata::travel_time(to = "city", size = 5, up = TRUE)
#           https://geodata.ucdavis.edu/geodata/travel/travel_time_to_cities_u5.tif
# LICENCE:  CC BY 4.0
# SIZE:     406 MB global GeoTIFF, 30 arcsec, cached under out/_geodata/
#           (gitignored; delete it when done)
# UNITS:    minutes of travel to the nearest city of >= 50,000 people
#
# Writes sandbox_parsimony/out/travel_time_admin2.csv.
# =============================================================================
suppressPackageStartupMessages({library(terra); library(sf); library(dplyr)})

GEO <- "sandbox_parsimony/out/_geodata"
OUT <- "sandbox_parsimony/out/travel_time_admin2.csv"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi", "tanzania")

r <- tryCatch(geodata::travel_time(to = "city", size = 5, up = TRUE, path = GEO),
              error = function(e) NULL)
if (is.null(r)) stop("travel_time raster unavailable; run once with network access")

rows <- list()
for (ctry in COUNTRIES) {
  ac <- tryCatch(readRDS(file.path("_targets_full/objects",
                                   paste0("area_covariates_", ctry))),
                 error = function(e) NULL)
  if (is.null(ac) || is.null(ac$polygons)) { message("no polygons: ", ctry); next }
  pg <- ac$polygons
  pv <- terra::vect(sf::st_transform(pg, terra::crs(r)))
  rc <- terra::crop(r, terra::ext(pv))          # crop before extract: 406 MB global
  ex <- terra::extract(rc, pv, fun = NULL, ID = TRUE)
  names(ex)[2] <- "tt"
  # Travel time is heavily right-skewed and a district mean is dominated by its
  # remotest pixel, so keep the median and the log-mean as well. log1p because
  # a district containing a city has pixels at 0 minutes.
  agg <- ex |> filter(is.finite(tt)) |> group_by(ID) |>
    summarise(gee_traveltime_city_mean   = mean(tt),
              gee_traveltime_city_median = stats::median(tt),
              gee_traveltime_city_logmean = mean(log1p(tt)),
              gee_traveltime_city_p90    = stats::quantile(tt, .9, names = FALSE),
              .groups = "drop")
  out <- data.frame(country = ctry,
                    Admin2 = as.character(pg$Admin2)[agg$ID],
                    stringsAsFactors = FALSE)
  for (v in setdiff(names(agg), "ID")) out[[v]] <- agg[[v]]
  out <- out[!duplicated(out$Admin2), , drop = FALSE]
  rows[[ctry]] <- out
  message(sprintf("  %-12s %4d districts, median travel time %.0f min (IQR %.0f-%.0f)",
                  ctry, nrow(out), stats::median(out$gee_traveltime_city_median),
                  stats::quantile(out$gee_traveltime_city_median, .25),
                  stats::quantile(out$gee_traveltime_city_median, .75)))
}

tt <- bind_rows(rows)
write.csv(tt, OUT, row.names = FALSE)
message(sprintf("\nSaved %s (%d rows)", OUT, nrow(tt)))

# sanity: does it agree with the per-country MAP clips where those exist?
old <- readRDS("sandbox_parsimony/out/cov_list.rds")
for (ctry in intersect(names(old), COUNTRIES)) {
  if (!"gee_accessibility" %in% names(old[[ctry]])) next
  cmp <- dplyr::inner_join(
    tt[tt$country == ctry, c("Admin2", "gee_traveltime_city_median")],
    old[[ctry]][, c("Admin2", "gee_accessibility")], by = "Admin2")
  if (nrow(cmp) > 10)
    message(sprintf("  check %-12s r(new travel time, existing MAP accessibility) = %.2f  (n=%d)",
                    ctry, suppressWarnings(cor(cmp[[2]], cmp[[3]],
                                               use = "complete.obs")), nrow(cmp)))
}
