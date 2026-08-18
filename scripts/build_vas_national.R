#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_vas_national.R
#
# Vitamin A supplementation (VAS) coverage as a NATIONAL predictor, broadcast to
# every GADM Admin-2 in each country -> data/VAS/<Country>_vas_admin2.csv
# with `vas_` columns (constant within country).
#
# WHY: two-dose VAS coverage in children 6-59 mo is the proximate driver of
# population vitamin A status and the covariate GBD leans on for its modelled VAD
# decline. Fortification (GFDx) is a different channel; VAS is the supplementation
# channel and nothing in the current predictor set captures it. See
# docs/ra_new_data_sources_specs.md §5.
#
# ACCESS (public, no key): World Bank WDI mirror of the UNICEF estimate,
# indicator SN.ITK.VITA.ZS ("Vitamin A supplementation coverage rate, % of
# children 6-59 months"), via the WDI package. National-only -> helps the
# cross-country pooled / transfer model, NOT within-country Admin-2 ranking.
#
# USAGE:  Rscript scripts/build_vas_national.R              # all countries
#         Rscript scripts/build_vas_national.R Ghana Tanzania
# =============================================================================
suppressMessages({library(here); library(dplyr)})
if (!requireNamespace("WDI", quietly = TRUE))
  stop("Install WDI: install.packages('WDI')")

ISO2 <- c(Gambia="GM", Ghana="GH", SierraLeone="SL", Malawi="MW", Tanzania="TZ")
ISO3 <- c(Gambia="GMB",Ghana="GHA",SierraLeone="SLE",Malawi="MWI",Tanzania="TZA")
SURVEY_YEAR <- c(Gambia=2018, Ghana=2017, SierraLeone=2013, Malawi=2016, Tanzania=2010)
OUT_DIR <- here::here("data","VAS"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

args <- commandArgs(trailingOnly=TRUE)
countries <- if (length(args)) args else names(ISO2)
bad <- setdiff(countries, names(ISO2)); if (length(bad)) stop("Unknown country: ", paste(bad, collapse=", "))

cat("## Vitamin A supplementation coverage (WDI SN.ITK.VITA.ZS)\n")
raw <- WDI::WDI(country = unname(ISO2[countries]),
                indicator = c(vas_vita_coverage_pct = "SN.ITK.VITA.ZS"),
                start = 2000, end = 2024, extra = FALSE)

source(here::here("R","data_prep.R"))  # load_gadm_cached()

pick_nearest <- function(df_c, yr) {
  d <- df_c[is.finite(df_c$vas_vita_coverage_pct), ]
  if (!nrow(d)) return(list(val=NA_real_, yr=NA_integer_))
  i <- which.min(abs(d$year - yr))
  list(val = d$vas_vita_coverage_pct[i], yr = d$year[i])
}

for (cn in countries) {
  cat("\n#####", cn, "#####\n")
  df_c <- raw[raw$iso2c == ISO2[[cn]], ]
  hit  <- pick_nearest(df_c, SURVEY_YEAR[[cn]])
  if (is.na(hit$val)) { cat("  no VAS value available; skipping\n"); next }
  cat(sprintf("  VAS coverage = %.1f%% (year %d; survey year %d)\n",
              hit$val, hit$yr, SURVEY_YEAR[[cn]]))
  gadm <- tryCatch(load_gadm_cached(ISO3[[cn]], level=2), error=function(e) NULL)
  if (is.null(gadm)) { cat("  GADM load failed; skipping\n"); next }
  g <- unique(sf::st_drop_geometry(gadm)[, intersect(c("NAME_1","NAME_2"), names(gadm)), drop=FALSE])
  names(g)[names(g)=="NAME_1"] <- "Admin1"; names(g)[names(g)=="NAME_2"] <- "Admin2"
  g$Admin1 <- trimws(g$Admin1); g$Admin2 <- trimws(g$Admin2)
  g$vas_vita_coverage_pct <- hit$val
  g$vas_year <- hit$yr
  fn <- file.path(OUT_DIR, sprintf("%s_vas_admin2.csv", cn))
  write.csv(g, fn, row.names=FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows)\n", basename(fn), nrow(g)))
}
cat("\nDONE. Domain VAS=list(prefix=\"vas_\") + guarded merge join are wired in config/merges.\n")
