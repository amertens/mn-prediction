# =============================================================================
# scripts/protocol_v2/38_food_price_timematched.R   [FP-02]
#
# FOOD PRICES MATCHED TO THE FIELDWORK WINDOW, WITH SEASONALITY
#
# Script 07 put nine WFP-derived price columns into the vocabulary, all
# multi-year averages (survey year +/- 2). Two things were wrong with that for
# a seasonal quantity. First, prices at the time of the blood draw are what
# bear on status, and the surveys fell in very different parts of the price
# cycle (Gambia Jan-Apr, Ghana Apr-Jun, Malawi Dec-Feb, the lean season).
# Second, script 07 took Gambia's survey year as 2021; the fieldwork dates
# (FW-01) put it in January-April 2018, so its window was three years late.
#
# This script builds market-level indicators the price literature uses and
# matches them to each district's fieldwork window from FW-01:
#   seasonal profile   per market x commodity series, log price minus that
#                      year's median (removes inflation and exchange-rate
#                      trend), median by calendar month over all years with
#                      >= 6 months; series with < 3 years borrow the
#                      country x commodity profile
#   fpt_*_seas_amp     max - min of the profile (log points ~ % seasonal gap;
#                      Gilbert, Christiaensen & Kaminski 2017)
#   fpt_*_seas_pos     where in the cycle the fieldwork fell: profile value
#                      averaged over the window months (+ = high-price season)
#   fpt_*_anom_z       ALPS-style anomaly: (deviation from the same-month norm)
#                      / series SD, averaged over the window (WFP ALPS; FAO IPA)
#   fpt_*_rel_win      window price relative to the national median for the
#                      same commodity+unit in the same months (spatial index)
#   fpt_staple_usd_kg  absolute staple price, USD per kg, window median
#                      (only units parseable to kg)
#   fpt_*_to_staple    nutritious-group relative price against the staple
#                      (relative-price-of-nutritious-foods, Headey & Alderman)
#   fpt_months_to_peak circular months from the fieldwork month to the staple
#                      price peak (0 = surveyed at the peak, 6 = at the trough)
# Market values reach districts by inverse-distance weighting over the three
# nearest markets (as in 07) and, separately, cluster GPS points, so the same
# columns exist for a future cluster-level model.
#
# Nothing here touches predictors_admin2_shared.csv; the add-on is tested with
# script 39 before anyone decides whether it enters the vocabulary.
#
#   Rscript scripts/protocol_v2/38_food_price_timematched.R
# -> data/covariates/harmonized/predictors_admin2_food_tm.csv (+ _metadata.csv)
# -> data/covariates/harmonized/predictors_cluster_food_tm.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"; HDIR <- "data/covariates/harmonized"
ISO <- c(Gambia = "gmb", Ghana = "gha", Malawi = "mwi", SierraLeone = "sle")
LC  <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
CATS <- c(staple = "cereals and tubers", animal = "meat, fish and eggs", pulses = "pulses and nuts", vegfruit = "vegetables and fruits", oils = "oil and fats")
GPS <- list(Gambia = c(f = "data/IPD/Gambia/Gambia_GMS_GPS_cleaned.csv", id = "MICS_Cluster_Number"),
            Ghana = c(f = "data/IPD/Ghana/Ghana_GMS_GPS_cleaned.csv", id = "cnum"),
            Malawi = c(f = "data/IPD/Malawi/Malawi_GMS_GPS_cleaned.csv", id = "gw_cnum"),
            SierraLeone = c(f = "data/IPD/Sierra Leone/Sierra Leone_GMS_GPS_cleaned.csv", id = "cnum"))
FW  <- read.csv(file.path(OUTDIR, "fieldwork_windows_admin2.csv"), stringsAsFactors = FALSE)
FWC <- read.csv(file.path(OUTDIR, "fieldwork_windows_cluster.csv"), stringsAsFactors = FALSE)
FW$date_first <- as.Date(FW$date_first); FW$date_last <- as.Date(FW$date_last); FW$date_med <- as.Date(FW$date_med)
FWC$date_med <- as.Date(FWC$date_med)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
cent <- do.call(rbind, lapply(names(LC), function(cn) { b <- BND[[LC[[cn]]]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = cn, Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
cent <- cent[is.finite(cent$lon) & is.finite(cent$lat), ]
haversine_km <- function(lon1, lat1, lon2, lat2) { r <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2; 2 * r * asin(pmin(1, sqrt(a))) }
kg_factor <- function(unit) { u <- toupper(trimws(unit)); out <- rep(NA_real_, length(u)); out[u == "KG"] <- 1; out[u == "L"] <- 1
  k <- suppressWarnings(as.numeric(sub("^([0-9.]+) ?KG$", "\\1", u))); k[!grepl("^[0-9.]+ ?KG$", u)] <- NA; out[is.finite(k)] <- k[is.finite(k)]
  g <- suppressWarnings(as.numeric(sub("^([0-9.]+) ?G$", "\\1", u))); g[!grepl("^[0-9.]+ ?G$", u)] <- NA; out[is.finite(g)] <- g[is.finite(g)] / 1000; out }
dist_mat <- function(px, py, mx, my) outer(seq_along(px), seq_along(mx), Vectorize(function(i, j) haversine_km(px[i], py[i], mx[j], my[j])))
idw3 <- function(D, val) vapply(seq_len(nrow(D)), function(i) { d <- D[i, ]; d[!is.finite(val)] <- Inf; o <- order(d)[seq_len(min(3, sum(is.finite(val))))]
  if (!length(o) || !is.finite(d[o[1]])) return(NA_real_); wt <- 1 / pmax(d[o], 1); sum(wt * val[o]) / sum(wt) }, 0)
nearest <- function(D, val) vapply(seq_len(nrow(D)), function(i) { d <- D[i, ]; d[!is.finite(val)] <- Inf; o <- which.min(d); if (!length(o) || !is.finite(d[o])) NA_real_ else val[o] }, 0)
circ_months <- function(a, b) { d <- abs(a - b) %% 12; pmin(d, 12 - d) }
month_seq <- function(d1, d2) { y1 <- as.integer(format(d1, "%Y")); m1 <- as.integer(format(d1, "%m")); y2 <- as.integer(format(d2, "%Y")); m2 <- as.integer(format(d2, "%m")); seq(y1 * 12 + m1 - 1, y2 * 12 + m2 - 1) }

dist_rows <- list(); clus_rows <- list(); notes <- list()
for (cn in names(ISO)) {
  f <- file.path("data", "food_price", sprintf("wfp_food_prices_%s.csv", ISO[[cn]])); if (!file.exists(f)) { notes[[cn]] <- "no WFP file"; next }
  w <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE, progress = FALSE)) |> as.data.frame()
  w <- w[!is.na(w$date) & !startsWith(as.character(w$date), "#"), ]
  w$year <- as.integer(substr(as.character(w$date), 1, 4)); w$month <- as.integer(substr(as.character(w$date), 6, 7)); w$ym <- w$year * 12 + w$month - 1
  w$price <- suppressWarnings(as.numeric(w$price)); w$usdprice <- suppressWarnings(as.numeric(w$usdprice))
  w$latitude <- suppressWarnings(as.numeric(w$latitude)); w$longitude <- suppressWarnings(as.numeric(w$longitude))
  w$cat <- names(CATS)[match(w$category, CATS)]
  w <- w[is.finite(w$price) & w$price > 0 & is.finite(w$latitude) & is.finite(w$longitude) & !is.na(w$cat) & is.finite(w$ym), ]
  fw <- FW[FW$country == cn, ]; fwc <- FWC[FWC$country == cn & !is.na(FWC$date_med), ]
  if (!nrow(fw)) { notes[[cn]] <- "no fieldwork dates (FW-01); time-matched features not built"; cat(cn, notes[[cn]], "\n"); next }
  span <- month_seq(min(fw$date_first), max(fw$date_last)); yr_med <- as.integer(stats::median(fw$year_med))
  # price type: retail if it exists around the survey, else wholesale
  near <- w[abs(w$year - yr_med) <= 1, ]
  ptype <- if (sum(near$pricetype == "Retail", na.rm = TRUE) > 200) "Retail" else "Wholesale"; w <- w[w$pricetype == ptype, ]
  w <- w[w$year <= yr_med + 2, ]                       # reference period ends two years after the survey
  w$series <- paste(w$market, w$commodity, w$unit, sep = "||"); w$lp <- log(w$price)
  # year level per series; seasonal profile from series-years with >= 6 months
  w <- w |> group_by(series, year) |> mutate(r = lp - stats::median(lp), n_in_year = dplyr::n()) |> ungroup()
  prof_s <- w |> filter(n_in_year >= 6) |> group_by(series, month) |> summarise(s = stats::median(r), .groups = "drop")
  nyrs <- w |> filter(n_in_year >= 6) |> group_by(series) |> summarise(n_years = dplyr::n_distinct(year), .groups = "drop")
  prof_c <- w |> filter(n_in_year >= 6) |> group_by(commodity, month) |> summarise(s_c = stats::median(r), .groups = "drop")
  w <- w |> left_join(nyrs, by = "series") |> left_join(prof_s, by = c("series", "month")) |> left_join(prof_c, by = c("commodity", "month"))
  w$s_use <- ifelse(is.finite(w$n_years) & w$n_years >= 3 & is.finite(w$s), w$s, w$s_c)
  w$e <- w$r - w$s_use
  w <- w |> group_by(series) |> mutate(sd_e = if (sum(is.finite(e)) >= 12) stats::sd(e, na.rm = TRUE) else NA_real_) |> ungroup()
  w$z <- ifelse(is.finite(w$sd_e) & w$sd_e > 0, w$e / w$sd_e, NA_real_)
  w$kg <- kg_factor(w$unit); w$usd_kg <- ifelse(is.finite(w$kg) & is.finite(w$usdprice), w$usdprice / w$kg, NA_real_)
  # national same-month median by commodity+unit, for the spatial relative price
  w <- w |> group_by(commodity, unit, ym) |> mutate(rel_nat = lp - stats::median(lp)) |> ungroup()
  # seasonal profile summary per series -> per market x cat (static)
  amp <- w |> filter(is.finite(s_use)) |> group_by(series, market, longitude, latitude, cat) |>
    summarise(amp = { p <- tapply(s_use, month, stats::median); if (length(p) >= 8) max(p) - min(p) else NA_real_ },
              peak = { p <- tapply(s_use, month, stats::median); if (length(p) >= 8) as.integer(names(p)[which.max(p)]) else NA_integer_ }, .groups = "drop") |>
    group_by(market, longitude, latitude, cat) |> summarise(seas_amp = stats::median(amp, na.rm = TRUE), peak_month = { p <- peak[is.finite(peak)]; if (length(p)) as.integer(round(stats::median(p))) else NA_integer_ }, .groups = "drop")
  # window-month values per market x cat x ym (over the country's fieldwork span)
  win <- w[w$ym %in% span, ] |> group_by(market, longitude, latitude, cat, ym) |>
    summarise(anom_z = stats::median(z, na.rm = TRUE), seas_pos = stats::median(s_use, na.rm = TRUE), rel_win = stats::median(rel_nat, na.rm = TRUE), usd_kg = stats::median(usd_kg, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
  mk <- w |> distinct(market, longitude, latitude)
  build_points <- function(px, py, keys, wins) {   # wins: list of month vectors per point
    D <- dist_mat(px, py, mk$longitude, mk$latitude); out <- keys
    for (k in c("staple", "pulses", "animal", "vegfruit")) {
      for (v in c("anom_z", "seas_pos", "rel_win", "usd_kg")) {
        if (v == "usd_kg" && k != "staple") next
        vals <- vapply(seq_along(px), function(i) { ms <- wins[[i]]; if (!length(ms)) return(NA_real_)
          x <- vapply(ms, function(m) { sub <- win[win$cat == k & win$ym == m, ]; if (!nrow(sub)) return(NA_real_)
            val <- sub[[v]][match(mk$market, sub$market)]; idw3(D[i, , drop = FALSE], val) }, 0); if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) }, 0)
        out[[paste0("fpt_", k, "_", v)]] <- round(vals, 4) }
      sub <- amp[amp$cat == k, ]; val <- sub$seas_amp[match(mk$market, sub$market)]
      out[[paste0("fpt_", k, "_seas_amp")]] <- round(idw3(D, val), 4)
      if (k == "staple") { pk <- sub$peak_month[match(mk$market, sub$market)]; out$fpt_staple_peak_month <- nearest(D, pk) }
    }
    for (k in c("animal", "pulses", "vegfruit")) out[[paste0("fpt_", k, "_to_staple")]] <- round(out[[paste0("fpt_", k, "_rel_win")]] - out$fpt_staple_rel_win, 4)
    out
  }
  # districts
  cc <- cent[cent$country == cn, ]; key <- paste(cc$Admin1, cc$Admin2); fk <- paste(fw$Admin1, fw$Admin2)
  wins <- lapply(seq_len(nrow(cc)), function(i) { j <- match(key[i], fk); if (is.na(j)) integer(0) else month_seq(fw$date_first[j], fw$date_last[j]) })
  dd <- build_points(cc$lon, cc$lat, cc[, c("country", "Admin1", "Admin2")], wins)
  dd$fpt_fieldwork_month <- fw$month_med[match(key, fk)]
  dd$fpt_months_to_peak <- circ_months(dd$fpt_fieldwork_month, dd$fpt_staple_peak_month)
  dist_rows[[cn]] <- dd
  # clusters
  g <- GPS[[cn]]; gp <- tryCatch(read.csv(g[["f"]], stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(gp) && all(c("latitude", "longitude", g[["id"]]) %in% names(gp))) {
    gp <- gp[is.finite(gp$latitude) & is.finite(gp$longitude), ]; cid <- as.character(gp[[g[["id"]]]])
    j <- match(cid, as.character(fwc$cluster))
    winc <- lapply(seq_along(cid), function(i) if (is.na(j[i])) integer(0) else { d <- fwc$date_med[j[i]]; month_seq(d, d) })
    cd <- build_points(gp$longitude, gp$latitude, data.frame(country = cn, cluster = cid, lat = gp$latitude, lon = gp$longitude, stringsAsFactors = FALSE), winc)
    cd$fpt_fieldwork_month <- fwc$month_med[j]; cd$fpt_months_to_peak <- circ_months(cd$fpt_fieldwork_month, cd$fpt_staple_peak_month)
    clus_rows[[cn]] <- cd
  }
  notes[[cn]] <- sprintf("%s prices, %d markets, fieldwork %s..%s (%d window months), %d districts dated of %d, %d clusters", ptype, nrow(mk), format(min(fw$date_first)), format(max(fw$date_last)), length(span), sum(!is.na(match(key, fk))), nrow(cc), length(clus_rows[[cn]][["cluster"]]))
  cat(sprintf("%-12s %s\n", cn, notes[[cn]]))
}
FT <- bind_rows(dist_rows); CT <- bind_rows(clus_rows)
write.csv(FT, file.path(HDIR, "predictors_admin2_food_tm.csv"), row.names = FALSE)
if (nrow(CT)) write.csv(CT, file.path(HDIR, "predictors_cluster_food_tm.csv"), row.names = FALSE)
newcols <- setdiff(names(FT), c("country", "Admin1", "Admin2"))
cov_by <- sapply(newcols, function(cl) paste(sprintf("%s=%.2f", names(ISO), vapply(names(ISO), function(cn) { z <- FT[[cl]][FT$country == cn]; if (!length(z)) 0 else mean(is.finite(z)) }, 0)), collapse = ";"))
MDT <- data.frame(column = newcols, domain = "Food prices, time-matched", source = "WFP / HDX market prices x FW-01 fieldwork windows",
  n_countries = vapply(newcols, function(cl) sum(vapply(names(ISO), function(cn) { z <- FT[[cl]][FT$country == cn]; length(z) > 0 && mean(is.finite(z)) > 0.5 }, TRUE)), 0L),
  completeness = round(vapply(newcols, function(cl) mean(is.finite(FT[[cl]])), 0), 3), subnational = TRUE, coverage_by_country = unname(cov_by), stringsAsFactors = FALSE)
write.csv(MDT, file.path(HDIR, "predictors_admin2_food_tm_metadata.csv"), row.names = FALSE)
cat("\n===== FP-02: time-matched food-price add-on =====\n"); print(MDT[, c("column", "n_countries", "completeness", "coverage_by_country")], row.names = FALSE)
cat("\nwithin-country summaries (median over districts):\n")
print(as.data.frame(FT |> group_by(country) |> summarise(across(all_of(c("fpt_staple_anom_z", "fpt_staple_seas_pos", "fpt_staple_seas_amp", "fpt_staple_usd_kg", "fpt_staple_rel_win", "fpt_months_to_peak", "fpt_fieldwork_month")), ~ round(stats::median(.x, na.rm = TRUE), 3)), .groups = "drop")), row.names = FALSE)
cat("\nDONE\n")
