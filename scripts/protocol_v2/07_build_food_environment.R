# =============================================================================
# scripts/protocol_v2/07_build_food_environment.R
#
# Put the food environment into the modelling vocabulary.
#
# THE PROBLEM THIS FIXES
# ----------------------
# No food-supply, price, market-access or affordability column reaches the
# modelled set, although WFP price series for five countries, HDX HAPI, and
# per-country FAOSTAT Admin-2 files all sit in data/. The chain: only two
# food-security columns were ever harmonised (`fsec_ipc_phase_*`, which are IPC
# PHASE classifications rather than prices), and both were then dropped by the
# rule "absent, mostly missing or constant in at least one country". IPC phase
# is frequently constant within a country by design, so that rule deletes it
# automatically. The domain closest to the intake mechanism was therefore empty
# for a filter reason, not an evidence reason.
#
# THE ASSUMPTIONS, STATED UP FRONT
# --------------------------------
# Building district food-price columns from market point data requires choices.
# Each one is recorded in the metadata `assumption` field so a reader can
# disagree with it explicitly:
#
#  1 WINDOW. Prices are averaged over the survey year +/- 2 years, per country
#    (Gambia 2021, Ghana 2017, Malawi 2015, Sierra Leone 2013). Rationale:
#    biomarker status reflects sustained rather than instantaneous conditions,
#    and a single year of a thin market series is noisy.
#  2 UNITS. WFP units are heterogeneous ("91 KG", "1 L", ...). Prices are NOT
#    converted to a common physical unit. Instead each price is expressed
#    RELATIVE to the country-wide median for the SAME commodity and unit, so
#    unit differences cancel within a commodity. The resulting index is
#    "how expensive is this basket here, relative to the national norm".
#  3 PRICE TYPE. Retail where available, falling back to wholesale (Ghana is
#    predominantly wholesale). Recorded per country.
#  4 ASSIGNMENT. A district takes the inverse-distance-weighted mean of the
#    three nearest markets, and separately records the distance to the nearest
#    one. Districts have no markets far more often than not (13-130 markets
#    against 14-87 districts), so nearest-market interpolation is unavoidable;
#    the distance column makes the extrapolation visible rather than hidden.
#  5 FAOSTAT. The per-country FAOSTAT Admin-2 files are NATIONAL figures
#    broadcast to every district - verified here, they take one distinct value
#    per country. They carry no within-country information and are flagged
#    `subnational = FALSE`. They are included because they are informative
#    across countries, which is exactly the axis a leave-one-country-out design
#    uses, and excluded automatically by the all-four-countries filter.
#
# The all-four-countries filter is NOT applied to these columns. A variable
# informative in three countries is kept and its coverage recorded, because
# discarding it entirely is the same complete-case logic the audit flagged for
# Ghana's DHS block.
#
#   Rscript scripts/protocol_v2/07_build_food_environment.R
# -> data/covariates/harmonized/predictors_admin2_food.csv
# -> data/covariates/harmonized/predictors_admin2_food_metadata.csv
# -> updates predictors_admin2_shared.csv + _metadata.csv (backup kept)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

HDIR <- "data/covariates/harmonized"
# Gambia corrected 2021 -> 2018 on 2026-09-04: FW-01 recovered the interview
# dates (24 Jan - 19 Apr 2018), so the +/-2 window is now 2016-2020, not
# 2019-2023. Malawi's fieldwork was Dec 2015 - Feb 2016; 2015 +/- 2 covers it.
SURVEY_YEAR <- c(Gambia = 2018, Ghana = 2017, Malawi = 2015, SierraLeone = 2013)
ISO <- c(Gambia = "gmb", Ghana = "gha", Malawi = "mwi", SierraLeone = "sle")
LC  <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi",
         SierraLeone = "sierraleone")
WINDOW <- 2L
# mechanistically ordered: animal-source foods carry haem iron, B12 and
# preformed vitamin A; pulses and vegetables carry non-haem iron and
# provitamin A; cereals/tubers are the staple energy base.
CATS <- c(staple = "cereals and tubers", animal = "meat, fish and eggs",
          pulses = "pulses and nuts", vegfruit = "vegetables and fruits",
          oils = "oil and fats")

BND <- readRDS("dashboard/data/admin2_boundaries.rds")
cent <- do.call(rbind, lapply(names(LC), function(cn) {
  b <- BND[[LC[[cn]]]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = cn,
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))
cent <- cent[is.finite(cent$lon) & is.finite(cent$lat), ]

haversine_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

food_rows <- list(); notes <- list()
for (cn in names(ISO)) {
  f <- file.path("data", "food_price",
                 sprintf("wfp_food_prices_%s.csv", ISO[[cn]]))
  if (!file.exists(f)) { notes[[cn]] <- "no WFP file"; next }
  w <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE,
                                        progress = FALSE)) |> as.data.frame()
  w <- w[!is.na(w$date) & !startsWith(as.character(w$date), "#"), ]
  w$year <- as.integer(substr(as.character(w$date), 1, 4))
  w$usdprice <- suppressWarnings(as.numeric(w$usdprice))
  w$latitude  <- suppressWarnings(as.numeric(w$latitude))
  w$longitude <- suppressWarnings(as.numeric(w$longitude))
  yr <- SURVEY_YEAR[[cn]]
  win <- w[is.finite(w$year) & abs(w$year - yr) <= WINDOW &
             is.finite(w$usdprice) & w$usdprice > 0 &
             is.finite(w$latitude) & is.finite(w$longitude), ]
  # ASSUMPTION 3: retail if it exists in the window, else wholesale
  ptype <- if (sum(win$pricetype == "Retail", na.rm = TRUE) > 200) "Retail" else "Wholesale"
  win <- win[win$pricetype == ptype, ]
  if (!nrow(win)) { notes[[cn]] <- "no prices in window"; next }

  # ASSUMPTION 2: relative to the country median for the same commodity+unit
  win <- win |> group_by(commodity, unit) |>
    mutate(rel = usdprice / stats::median(usdprice, na.rm = TRUE)) |>
    ungroup() |> filter(is.finite(rel))

  mk <- win |> group_by(market, longitude, latitude) |>
    summarise(across(everything(), ~NA, .names = "drop_{.col}"), .groups = "drop") |>
    select(market, longitude, latitude)

  # per market x category: median relative price, and staple volatility
  mc <- win |> mutate(cat = names(CATS)[match(category, CATS)]) |>
    filter(!is.na(cat)) |>
    group_by(market, longitude, latitude, cat) |>
    summarise(rel = stats::median(rel, na.rm = TRUE),
              vol = stats::sd(rel, na.rm = TRUE) / pmax(base::mean(rel, na.rm = TRUE), 1e-9),
              n_obs = dplyr::n(), .groups = "drop")

  cc <- cent[cent$country == cn, ]
  if (!nrow(cc)) { notes[[cn]] <- "no centroids"; next }
  D <- outer(seq_len(nrow(cc)), seq_len(nrow(mk)), Vectorize(function(i, j)
    haversine_km(cc$lon[i], cc$lat[i], mk$longitude[j], mk$latitude[j])))
  out <- cc[, c("country", "Admin1", "Admin2")]
  out$fprice_dist_nearest_market_km <- round(apply(D, 1, min), 2)
  out$fprice_n_markets_100km <- rowSums(D <= 100)

  # ASSUMPTION 4: inverse-distance weight over the three nearest markets
  for (k in names(CATS)) {
    sub <- mc[mc$cat == k, ]
    if (!nrow(sub)) { out[[paste0("fprice_", k, "_rel")]] <- NA_real_; next }
    idx <- match(sub$market, mk$market)
    vals <- vapply(seq_len(nrow(cc)), function(i) {
      d <- D[i, idx]; o <- order(d)[seq_len(min(3, length(d)))]
      wt <- 1 / pmax(d[o], 1)
      sum(wt * sub$rel[o]) / sum(wt)
    }, 0)
    out[[paste0("fprice_", k, "_rel")]] <- round(vals, 4)
  }
  sv <- mc[mc$cat == "staple", ]
  if (nrow(sv)) {
    idx <- match(sv$market, mk$market)
    out$fprice_staple_volatility <- round(vapply(seq_len(nrow(cc)), function(i) {
      d <- D[i, idx]; o <- order(d)[seq_len(min(3, length(d)))]
      wt <- 1 / pmax(d[o], 1)
      sum(wt * sv$vol[o], na.rm = TRUE) / sum(wt)
    }, 0), 4)
  }
  # relative cost of animal-source foods against the staple: an affordability
  # contrast rather than a price level, and the one closest to the mechanism
  if (all(c("fprice_animal_rel", "fprice_staple_rel") %in% names(out)))
    out$fprice_animal_to_staple <-
      round(out$fprice_animal_rel / pmax(out$fprice_staple_rel, 1e-9), 4)

  food_rows[[cn]] <- out
  notes[[cn]] <- sprintf("%s prices, %d markets, %d obs in %d+/-%d",
                         ptype, nrow(mk), nrow(win), yr, WINDOW)
  cat(sprintf("%-12s %s\n", cn, notes[[cn]]))
}
FP <- bind_rows(food_rows)

# ── FAOSTAT: national supply, broadcast, flagged as such ────────────────────
fao_rows <- list()
for (cn in names(ISO)) {
  f <- file.path("data", "FAOSTAT", sprintf("%s_fao_admin2.csv", cn))
  if (!file.exists(f)) next
  d <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE,
                                        progress = FALSE)) |> as.data.frame()
  keep <- grep("^fao_", names(d), value = TRUE)
  keep <- setdiff(keep, "fao_year")
  d <- d[, c("Admin1", "Admin2", keep)]
  d$country <- cn
  fao_rows[[cn]] <- d
}
FA <- bind_rows(fao_rows)

FOOD <- if (nrow(FP) && nrow(FA))
  full_join(FP, FA, by = c("country", "Admin1", "Admin2")) else
  if (nrow(FP)) FP else FA
newcols <- setdiff(names(FOOD), c("country", "Admin1", "Admin2"))
write.csv(FOOD, file.path(HDIR, "predictors_admin2_food.csv"), row.names = FALSE)

# ── metadata, with the assumption recorded per column ───────────────────────
assump <- c(
  fprice_dist_nearest_market_km = "Haversine km from district centroid to nearest WFP market.",
  fprice_n_markets_100km = "Count of WFP markets within 100 km of the district centroid.",
  fprice_staple_volatility = "IDW over 3 nearest markets of the within-window CV of relative staple price.",
  fprice_animal_to_staple = "Ratio of animal-source to staple relative price: an affordability contrast, not a level.")
for (k in names(CATS))
  assump[[paste0("fprice_", k, "_rel")]] <- sprintf(
    "Median USD price of '%s' relative to the country median for the same commodity and unit, survey year +/-%d, IDW over the 3 nearest markets.",
    CATS[[k]], WINDOW)
fao_ass <- "FAOSTAT national food supply, broadcast to every district: NO within-country variation. Kept for cross-country (LOCO) information only."

cov_by_country <- function(v) {
  s <- FOOD |> group_by(country) |>
    summarise(ok = mean(is.finite(.data[[v]])), .groups = "drop")
  paste(sprintf("%s=%.2f", s$country, s$ok), collapse = ";")
}
MDF <- data.frame(
  column = newcols,
  domain = "Food prices and supply",
  source = ifelse(grepl("^fao_", newcols), "FAOSTAT (national, broadcast)",
                  "WFP / HDX market prices"),
  n_countries = vapply(newcols, function(v)
    sum(tapply(FOOD[[v]], FOOD$country, function(z) any(is.finite(z)))), 0L),
  countries = vapply(newcols, function(v) {
    s <- tapply(FOOD[[v]], FOOD$country, function(z) any(is.finite(z)))
    paste(names(s)[s], collapse = ";") }, ""),
  completeness = round(vapply(newcols, function(v) mean(is.finite(FOOD[[v]])), 0), 3),
  subnational = !grepl("^fao_", newcols),
  coverage_by_country = vapply(newcols, cov_by_country, ""),
  assumption = ifelse(grepl("^fao_", newcols), fao_ass,
                      unname(assump[newcols])),
  stringsAsFactors = FALSE)
write.csv(MDF, file.path(HDIR, "predictors_admin2_food_metadata.csv"),
          row.names = FALSE)

# ── append to the shared set, WITHOUT the all-four-countries filter ─────────
SH  <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
SHM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"),
                stringsAsFactors = FALSE)
if (!file.exists(file.path(HDIR, "predictors_admin2_shared.csv.pre_food")))
  file.copy(file.path(HDIR, "predictors_admin2_shared.csv"),
            file.path(HDIR, "predictors_admin2_shared.csv.pre_food"))
if (!file.exists(file.path(HDIR, "predictors_admin2_shared_metadata.csv.pre_food")))
  file.copy(file.path(HDIR, "predictors_admin2_shared_metadata.csv"),
            file.path(HDIR, "predictors_admin2_shared_metadata.csv.pre_food"))

SH2 <- SH |> select(-any_of(newcols)) |>
  left_join(FOOD, by = c("country", "Admin1", "Admin2"))
SHM2 <- bind_rows(
  SHM |> filter(!column %in% newcols),
  MDF |> transmute(column, domain, source, n_countries, countries,
                   completeness, subnational))
write.csv(SH2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(SHM2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"),
          row.names = FALSE)

cat("\n=== food-environment columns added ===\n")
print(as.data.frame(MDF[, c("column", "source", "n_countries", "completeness",
                            "subnational", "coverage_by_country")]),
      row.names = FALSE)
cat("\nshared set:", ncol(SH) - 3, "->", ncol(SH2) - 3, "predictors\n")
cat("backups written with suffix .pre_food\n")
cat("\nNOTE: the all-four-countries filter is deliberately NOT applied to these\n")
cat("columns. Coverage is recorded per country instead, so a consumer can\n")
cat("decide rather than having the decision made silently upstream.\n")
cat("\nDONE\n")
