# =============================================================================
# scripts/covariates/build_rtfp_price_block.R   [RT-01, 2026-09-15]
#
# LOCAL FOOD PRICES AROUND THE SURVEY, FROM THE WORLD BANK REAL-TIME FOOD
# PRICES (RTFP) MARKET PANEL
#
# RTFP (Andree 2021; vintage 2026-02-10, the RA's download) gives a monthly,
# market-level food price index and staple prices from 2007 on, built from
# WFP / national price monitoring and filled by the RTFP model where the
# monitoring has gaps. Of the four countries only The Gambia (28 markets, 22
# districts) and Malawi (129 markets, 31 districts) are in the panel; Ghana
# and Sierra Leone are not (the world file was checked), so the block is
# NA there and the metadata says so. The index is on a common national base,
# so a market's level can be compared with the national mean at the same
# date (checked: the Malawi index tracks the maize price across markets).
#
# Exposure window: the 12 months ending in the last fieldwork month of the
# micronutrient survey, so the prices are the ones the sampled households
# faced in the year before their blood draw:
#   Malawi      MNS 2015-16, fieldwork Dec 2015 - Feb 2016 -> Mar 2015 - Feb 2016
#               (the El Nino lean season: maize prices peaked in Jan-Feb 2016)
#   The Gambia  GMNS 2018,   fieldwork Jan - Apr 2018      -> May 2017 - Apr 2018
# (fieldwork dates from metadata/survey_years.csv, the FW-01 cluster dates).
# The "previous year" is the 12 months before the window.
#
# Market indicators (log points unless stated):
#   fpi_rel_national         mean log food price index in the window, minus the
#                            cross-market mean (how expensive food is here
#                            relative to the country)
#   fpi_inflation_12m        mean log index in the window minus the previous year
#   fpi_volatility           sd of month-on-month log index changes over the 24
#                            months ending at the window end
#   fpi_seasonal_range       max minus min log index within the window
#   staple_rel_national      the same relative level for the main staple
#                            (maize in Malawi, rice in The Gambia)
#   staple_inflation_12m     staple inflation over the window
# Admin-2 values are inverse-distance-weighted means of the 3 nearest markets
# to the unit centroid (weight 1 / max(distance, 2 km)), plus
#   rtfp_dist_market_km      distance from the unit centroid to the nearest market
#   rtfp_n_markets           markets inside the unit polygon (audit column,
#                            not appended to the shared set)
#
#   Rscript -e "source('scripts/covariates/build_rtfp_price_block.R')"
# -> data/covariates/harmonized/predictors_admin2_rtfp.csv (+ _metadata.csv)
# -> metadata/rtfp_market_indicators.csv   (one row per market, for audit)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
RTFP <- "data/RA_2026-09/extracted/RTFP market food price"
HDIR <- "data/covariates/harmonized"
K_NEAREST <- 3L; MIN_KM <- 2
num <- function(x) suppressWarnings(as.numeric(x))

SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
spine <- SH[, c("country", "Admin1", "Admin2")]
BND <- readRDS("dashboard/data/admin2_boundaries.rds")

# window end = the month of the survey's last fieldwork date, from the single
# source of survey dates (metadata/survey_years.csv, FW-01); staple = the main
# cereal in the panel
SY <- read.csv("metadata/survey_years.csv", stringsAsFactors = FALSE)
win_end <- function(cn) { e <- as.Date(SY$fieldwork_end[SY$country == cn]); as.Date(format(e, "%Y-%m-01")) }
fw_lab  <- function(cn) sprintf("%s to %s", format(as.Date(SY$fieldwork_start[SY$country == cn]), "%b %Y"), format(as.Date(SY$fieldwork_end[SY$country == cn]), "%b %Y"))
WIN <- list(
  Gambia = list(file = "GMB_RTFP_mkt_2007_2026-02-10.csv", bnd = "gambia", end = win_end("Gambia"), staple = "c_rice",  fieldwork = fw_lab("Gambia")),
  Malawi = list(file = "MWI_RTFP_mkt_2007_2026-02-10.csv", bnd = "malawi", end = win_end("Malawi"), staple = "c_maize", fieldwork = fw_lab("Malawi")))
months_back <- function(end, n) seq(end, by = "-1 month", length.out = n)   # end and the n-1 months before it

mk_rows <- list(); out_rows <- list()
for (cn in names(WIN)) {
  w <- WIN[[cn]]; cat(sprintf("\n[%s] window ending %s (fieldwork %s), staple %s\n", cn, format(w$end, "%b %Y"), w$fieldwork, sub("^c_", "", w$staple)))
  d <- read.csv(file.path(RTFP, w$file), check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$mkt_name != "Market Average" & is.finite(num(d$lat)) & is.finite(num(d$lon)), ]
  d$date <- as.Date(d$price_date); d$fpi <- log(num(d$c_food_price_index)); d$stp <- log(num(d[[w$staple]]))
  cur <- months_back(w$end, 12); prev <- months_back(cur[12], 13)[-1]; vol <- months_back(w$end, 25)   # 25 months -> 24 differences
  M <- d |> group_by(mkt_name, adm1_name, adm2_name, lat = num(lat), lon = num(lon)) |> arrange(date, .by_group = TRUE) |>
    summarise(n_months = sum(date %in% cur),
              fpi_cur = mean(fpi[date %in% cur]), fpi_prev = mean(fpi[date %in% prev]),
              fpi_volatility = sd(diff(fpi[date %in% vol])), fpi_seasonal_range = diff(range(fpi[date %in% cur])),
              stp_cur = mean(stp[date %in% cur]), stp_prev = mean(stp[date %in% prev]),
              coverage = first(num(data_coverage)), interpolated = first(num(spatially_interpolated)), .groups = "drop") |>
    filter(n_months == 12L) |>
    mutate(fpi_rel_national = fpi_cur - mean(fpi_cur), fpi_inflation_12m = fpi_cur - fpi_prev,
           staple_rel_national = stp_cur - mean(stp_cur), staple_inflation_12m = stp_cur - stp_prev, country = cn)
  cat(sprintf("  %d markets with the full window | food price level spread (sd of log) %.3f | national inflation %.1f%% | staple inflation %.1f%%\n",
              nrow(M), sd(M$fpi_rel_national), 100 * (exp(mean(M$fpi_inflation_12m)) - 1), 100 * (exp(mean(M$staple_inflation_12m)) - 1)))
  mk_rows[[cn]] <- M

  # Admin-2 units: centroid distances to every market, IDW over the k nearest
  b <- BND[[w$bnd]]; b$Admin1 <- as.character(b$Admin1); b$Admin2 <- as.character(b$Admin2)
  cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(b)))
  pts  <- sf::st_as_sf(M, coords = c("lon", "lat"), crs = 4326)
  D <- matrix(as.numeric(sf::st_distance(cent, pts)) / 1000, nrow = nrow(b))            # km
  inside <- lengths(sf::st_intersects(b, pts))
  ind <- c("fpi_rel_national", "fpi_inflation_12m", "fpi_volatility", "fpi_seasonal_range", "staple_rel_national", "staple_inflation_12m")
  U <- data.frame(country = cn, Admin1 = b$Admin1, Admin2 = b$Admin2, stringsAsFactors = FALSE)
  for (v in ind) U[[paste0("rtfp_", v)]] <- vapply(seq_len(nrow(b)), function(i) {
    o <- order(D[i, ])[seq_len(min(K_NEAREST, ncol(D)))]; wt <- 1 / pmax(D[i, o], MIN_KM); sum(wt * M[[v]][o]) / sum(wt) }, 0)
  U$rtfp_dist_market_km <- apply(D, 1, min); U$rtfp_n_markets <- inside
  cat(sprintf("  %d Admin-2 units | %d with a market inside | median distance to nearest market %.1f km (max %.1f)\n",
              nrow(U), sum(inside > 0), median(U$rtfp_dist_market_km), max(U$rtfp_dist_market_km)))
  out_rows[[cn]] <- U
}
MK <- bind_rows(mk_rows); write.csv(MK, "metadata/rtfp_market_indicators.csv", row.names = FALSE)
OUT <- bind_rows(out_rows)
OUT <- left_join(spine, OUT, by = c("country", "Admin1", "Admin2"))
stopifnot(nrow(OUT) == nrow(spine), sum(is.finite(OUT$rtfp_fpi_rel_national)) == sum(spine$country %in% names(WIN)))
write.csv(OUT, file.path(HDIR, "predictors_admin2_rtfp.csv"), row.names = FALSE)

cols <- grep("^rtfp_", names(OUT), value = TRUE)
desc <- c(rtfp_fpi_rel_national = "Mean log RTFP food price index over the 12 months ending in the survey's last fieldwork month, minus the national cross-market mean (local price level; index on a common national base).",
          rtfp_fpi_inflation_12m = "Log change in the mean food price index between the 12-month exposure window and the 12 months before it.",
          rtfp_fpi_volatility = "SD of month-on-month log changes in the food price index over the 24 months ending at the window end.",
          rtfp_fpi_seasonal_range = "Max minus min log food price index within the 12-month window (seasonal price swing).",
          rtfp_staple_rel_national = "Local staple price level relative to the national cross-market mean (maize in Malawi, rice in The Gambia), log points, 12-month window.",
          rtfp_staple_inflation_12m = "Log change in the staple price between the window and the previous 12 months.",
          rtfp_dist_market_km = "Great-circle distance (km) from the Admin-2 centroid to the nearest RTFP market.",
          rtfp_n_markets = "RTFP markets inside the Admin-2 polygon (audit column; not appended to the shared set).")
md <- data.frame(column = cols, source = "World Bank Real-Time Food Prices (RTFP) market panel, vintage 2026-02-10", domain = "Market prices (RTFP)", subnational = TRUE,
                 assumption = paste(unname(desc[cols]), "Admin-2 value = inverse-distance-weighted mean of the 3 nearest markets to the unit centroid (weight 1/max(km, 2)). Windows end in the survey's last fieldwork month (metadata/survey_years.csv): Malawi Mar 2015 - Feb 2016, The Gambia May 2017 - Apr 2018. Ghana and Sierra Leone are not in the RTFP panel (NA)."),
                 stringsAsFactors = FALSE)
md$n_countries <- vapply(md$column, function(v) sum(tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
md$countries <- vapply(md$column, function(v) { s <- tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))); paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
md$completeness <- round(vapply(md$column, function(v) mean(is.finite(OUT[[v]])), 0), 3)
write.csv(md, file.path(HDIR, "predictors_admin2_rtfp_metadata.csv"), row.names = FALSE)

cat("\n=== RTFP block ===\n"); print(md[, c("column", "n_countries", "countries", "completeness")], row.names = FALSE)
cat("\ncountry summaries (unit means):\n")
print(OUT |> filter(is.finite(rtfp_fpi_rel_national)) |> group_by(country) |>
        summarise(n = n(), across(c(rtfp_fpi_rel_national, rtfp_fpi_inflation_12m, rtfp_fpi_volatility, rtfp_fpi_seasonal_range, rtfp_staple_inflation_12m, rtfp_dist_market_km), ~ round(mean(.x), 3))), width = 200)
cat("\nDONE\n")
