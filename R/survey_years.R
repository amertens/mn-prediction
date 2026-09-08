# =============================================================================
# R/survey_years.R
#
# ONE SOURCE OF TRUTH FOR THE SURVEY YEAR OF EACH COUNTRY
#
# Every extractor that matches a time-varying layer to a survey used to carry
# its own year table, and they drifted (AU-01 finding 7: Gambia 2020 in the
# Earth Engine scripts and 2021 in one harmoniser against 2018 elsewhere;
# Malawi 2015 in the protocol scripts against 2016 in R/config.R). The table
# metadata/survey_years.csv is built from the interview dates of every dated
# cluster (FW-01): survey_year is the calendar year of the respondent-weighted
# median interview date, so Malawi (fieldwork 8 Dec 2015 to 15 Feb 2016,
# median 22 Jan 2016, 78% of respondents in 2016) is 2016.
#
#   survey_years()                    c(Gambia = 2018, Ghana = 2017, Malawi = 2016, SierraLeone = 2013)
#   survey_years(keys = "lower")      the same with lower-case names
#   survey_years(protocol_only = FALSE)   adds Tanzania (2010), used by the national-supply builders
#   survey_years_table()              the full table (fieldwork window, median date, provenance)
#
# Python extractors read the same file through scripts/protocol_v2/survey_years.py.
# Auto-sourced by tar_source("R/"); scripts that do not source R/ call
# source("R/survey_years.R") (or here::here("R", "survey_years.R")) first.
# =============================================================================
.survey_years_file <- function() {
  d <- getwd()
  for (i in 1:6) {
    f <- file.path(d, "metadata", "survey_years.csv")
    if (file.exists(f)) return(f)
    d <- dirname(d)
  }
  stop("metadata/survey_years.csv not found in or above ", getwd())
}

survey_years_table <- function() {
  t <- read.csv(.survey_years_file(), stringsAsFactors = FALSE)
  t$survey_year <- as.integer(t$survey_year)
  t$in_protocol <- as.logical(toupper(as.character(t$in_protocol)))
  t
}

survey_years <- function(countries = NULL, keys = c("title", "lower"), protocol_only = TRUE) {
  keys <- match.arg(keys)
  t <- survey_years_table()
  if (protocol_only) t <- t[t$in_protocol, ]
  v <- stats::setNames(t$survey_year, t$country)
  if (keys == "lower") names(v) <- tolower(names(v))
  if (!is.null(countries)) {
    miss <- setdiff(countries, names(v))
    if (length(miss)) stop("survey_years(): no row for ", paste(miss, collapse = ", "))
    v <- v[countries]
  }
  v
}
