# =============================================================================
# scripts/protocol_v2/60_build_acled_conflict.R   [AB-02, 2026-09-15]
#
# Conflict exposure per Admin-2 from an ACLED data export: the one domain in
# the West Africa data-landscape review (Table 1, ACLED) that the shared set
# had no column for at all.
#
# INPUT (a one-time manual step, ACLED requires a registered account and
# does not permit redistribution):
#   1. https://acleddata.com -> Data Export Tool, logged in.
#   2. Countries: Gambia, Ghana, Malawi, Sierra Leone (and any country you
#      plan to add). Dates: 2007-01-01 to the latest. All event types.
#   3. Save the CSV(s) under data/ACLED/ (any file name, *.csv). Several files
#      are read and stacked; duplicates on event_id_cnty are removed.
#
# WHAT IS BUILT, per district, over the 36 months ending at the survey's
# fieldwork end (metadata/survey_years.csv):
#   acled_events_36m            events per 100,000 people (WorldPop)
#   acled_fatalities_36m        fatalities per 100,000 people
#   acled_violent_events_36m    battles + explosions/remote violence + violence
#                               against civilians, per 100,000
#   acled_civilian_targeting_36m events with civilian targeting, per 100,000
#                               (ACLED codes this from 2023 exports; NA when the
#                               column is absent)
#   acled_protest_riot_36m      protests + riots per 100,000
#   acled_months_with_event_36m share of the 36 months with at least one event
#   acled_any_event_36m         1 if any event in the window
#   acled_dist_nearest_event_km distance from the district centroid to the
#                               nearest event in the window (km), any type
# Point-in-polygon on the GADM Admin-2 polygons (data/admin_boundaries); events
# with geo_precision 3 (only the admin-1 is known) are kept for the national
# count but dropped from the district assignment.
#
# For these four surveys (2013-2018) most districts have zero events, so the
# block mainly matters for the next countries (Nigeria, Burkina Faso, Mali,
# Niger); it is built now so the vocabulary exists when they arrive.
#
#   Rscript -e "source('scripts/protocol_v2/60_build_acled_conflict.R')"
# -> data/covariates/harmonized/predictors_admin2_conflict.csv (+ _metadata.csv)
# -> updates predictors_admin2_shared.csv + _metadata.csv (.pre_conflict backup)
# Order in a rebuild: builder -> 07 -> 08 -> 59 -> 60 (this) -> 53.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
HDIR <- "data/covariates/harmonized"
SY <- read.csv("metadata/survey_years.csv", stringsAsFactors = FALSE)
SY <- SY[SY$in_protocol %in% c(TRUE, "TRUE"), ]
GADM <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
ACLED_NAME <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
WINDOW_MONTHS <- 36L
num <- function(x) suppressWarnings(as.numeric(x))

files <- list.files("data/ACLED", pattern = "[.]csv$", full.names = TRUE)
if (!length(files)) stop("no ACLED export under data/ACLED/ - see the header of this script for the manual download step")
A <- bind_rows(lapply(files, function(f) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)))
if ("event_id_cnty" %in% names(A)) A <- A[!duplicated(A$event_id_cnty), ]
need <- c("event_date", "country", "latitude", "longitude", "event_type", "fatalities")
stopifnot(all(need %in% names(A)))
A$event_date <- as.Date(A$event_date, tryFormats = c("%Y-%m-%d", "%d %B %Y", "%d-%b-%y", "%m/%d/%Y"))
A$geo_precision <- if ("geo_precision" %in% names(A)) num(A$geo_precision) else 1
A$civ <- if ("civilian_targeting" %in% names(A)) nzchar(as.character(A$civilian_targeting)) & !is.na(A$civilian_targeting) else NA
cat(sprintf("[acled] %d events in %d file(s), %s to %s\n", nrow(A), length(files), min(A$event_date, na.rm = TRUE), max(A$event_date, na.rm = TRUE)))

SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
spine <- SH[, c("country", "Admin1", "Admin2")]
POP <- tryCatch({ p <- readRDS("dashboard/data/admin2_population.rds"); p$country <- gsub(" ", "", p$country); p }, error = function(e) NULL)

rows <- list()
for (cn in names(GADM)) {
  s <- SY[SY$country == cn, ]; if (!nrow(s)) next
  end <- as.Date(s$fieldwork_end[1]); start <- seq(end, by = sprintf("-%d months", WINDOW_MONTHS), length.out = 2)[2]
  a <- A[A$country == ACLED_NAME[[cn]] & !is.na(A$event_date) & A$event_date > start & A$event_date <= end, ]
  poly <- sf::st_transform(load_gadm_cached(GADM[[cn]], level = 2), 4326)
  a1 <- as.character(poly$NAME_1); a2 <- as.character(poly$NAME_2)
  sp <- spine[spine$country == cn, ]
  cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(poly)))
  ok <- is.finite(num(a$latitude)) & is.finite(num(a$longitude)) & a$geo_precision < 3
  pts <- if (any(ok)) sf::st_as_sf(a[ok, ], coords = c("longitude", "latitude"), crs = 4326) else NULL
  ix <- if (!is.null(pts)) vapply(sf::st_within(pts, poly, sparse = TRUE), function(z) if (length(z)) z[1] else NA_integer_, integer(1)) else integer(0)
  ev <- if (!is.null(pts)) data.frame(i = ix, month = format(pts$event_date, "%Y-%m"), type = pts$event_type,
                                      fat = num(pts$fatalities), civ = pts$civ, stringsAsFactors = FALSE) else data.frame()
  ev <- ev[!is.na(ev$i), , drop = FALSE]
  agg <- function(i, f) if (nrow(ev)) f(ev[ev$i == i, , drop = FALSE]) else f(ev)
  out <- data.frame(country = cn, Admin1 = a1, Admin2 = a2, stringsAsFactors = FALSE)
  out$n_events <- vapply(seq_len(nrow(poly)), function(i) agg(i, nrow), 0)
  out$n_fat    <- vapply(seq_len(nrow(poly)), function(i) agg(i, function(d) sum(d$fat, na.rm = TRUE)), 0)
  out$n_viol   <- vapply(seq_len(nrow(poly)), function(i) agg(i, function(d) sum(d$type %in% c("Battles", "Explosions/Remote violence", "Violence against civilians"))), 0)
  out$n_civ    <- vapply(seq_len(nrow(poly)), function(i) agg(i, function(d) if (all(is.na(d$civ))) NA_real_ else sum(d$civ, na.rm = TRUE)), 0)
  out$n_prot   <- vapply(seq_len(nrow(poly)), function(i) agg(i, function(d) sum(d$type %in% c("Protests", "Riots"))), 0)
  out$acled_months_with_event_36m <- vapply(seq_len(nrow(poly)), function(i) agg(i, function(d) length(unique(d$month))), 0) / WINDOW_MONTHS
  out$acled_any_event_36m <- as.numeric(out$n_events > 0)
  out$acled_dist_nearest_event_km <- if (!is.null(pts) && nrow(pts)) {
    d <- sf::st_distance(cent, sf::st_geometry(pts)); as.numeric(apply(d, 1, min)) / 1000 } else NA_real_
  # per 100,000 people where a population is on file, else raw counts (flagged in the metadata)
  pop <- if (!is.null(POP)) POP$population[match(paste(cn, a1, a2), paste(POP$country, POP$Admin1, POP$Admin2))] else NA
  per <- function(x) if (all(is.na(pop))) x else x / pmax(pop, 1) * 1e5
  out$acled_events_36m <- per(out$n_events); out$acled_fatalities_36m <- per(out$n_fat)
  out$acled_violent_events_36m <- per(out$n_viol); out$acled_civilian_targeting_36m <- per(out$n_civ)
  out$acled_protest_riot_36m <- per(out$n_prot)
  out <- out[, c("country", "Admin1", "Admin2", grep("^acled_", names(out), value = TRUE))]
  cat(sprintf("  %-12s window %s..%s: %d events in country, %d geolocated to a district (%d districts with any), pop weighting: %s\n",
              cn, start, end, nrow(a), nrow(ev), sum(out$acled_any_event_36m), if (all(is.na(pop))) "NONE (raw counts)" else "WorldPop"))
  rows[[cn]] <- out
}
C <- bind_rows(rows)
cc <- grep("^acled_", names(C), value = TRUE)
write.csv(C, file.path(HDIR, "predictors_admin2_conflict.csv"), row.names = FALSE)
MD <- data.frame(column = cc, source = "ACLED (Armed Conflict Location & Event Data Project) export",
                 domain = "Conflict and insecurity", subnational = TRUE,
                 assumption = "36-month window ending at the survey's fieldwork end (metadata/survey_years.csv); events assigned to GADM Admin-2 by point-in-polygon, geo_precision 3 excluded; rates per 100,000 people from dashboard/data/admin2_population.rds where present. ACLED data are licensed, not redistributable: the export must be downloaded by the user.",
                 stringsAsFactors = FALSE)
MD$n_countries <- vapply(MD$column, function(v) sum(tapply(C[[v]], C$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
MD$countries <- vapply(MD$column, function(v) { s <- tapply(C[[v]], C$country, function(z) any(is.finite(z))); paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
MD$completeness <- round(vapply(MD$column, function(v) mean(is.finite(C[[v]])), 0), 3)
write.csv(MD, file.path(HDIR, "predictors_admin2_conflict_metadata.csv"), row.names = FALSE)

SHM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
SHM$subnational <- as.logical(toupper(as.character(SHM$subnational)))
for (f in c("predictors_admin2_shared.csv", "predictors_admin2_shared_metadata.csv"))
  if (!file.exists(file.path(HDIR, paste0(f, ".pre_conflict")))) file.copy(file.path(HDIR, f), file.path(HDIR, paste0(f, ".pre_conflict")))
SH2 <- SH |> select(-any_of(cc)) |> left_join(C, by = c("country", "Admin1", "Admin2"))
SHM2 <- bind_rows(SHM |> filter(!column %in% cc), MD |> transmute(column, domain, source, n_countries, countries, completeness, subnational))
stopifnot(nrow(SH2) == nrow(SH), ncol(SH2) - 3L == nrow(SHM2))
write.csv(SH2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(SHM2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"), row.names = FALSE)
cat("\nshared set:", ncol(SH) - 3, "->", ncol(SH2) - 3, "predictors\nDONE\n")
