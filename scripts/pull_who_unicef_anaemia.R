# =============================================================================
# scripts/pull_who_unicef_anaemia.R
#
# Pull WHO GHO modelled-prevalence estimates for anaemia (children 6-59m and
# women 15-49) — these are the WHO/UNICEF joint Bayesian estimates also
# disseminated through data.unicef.org. The GHO OData API is fully public, no
# auth needed.
#
# Output:
#   results/external/who_anaemia_country.csv
#     country, year, indicator, prevalence, value_low, value_high, n_pop
# =============================================================================
suppressPackageStartupMessages({
  library(here); library(httr); library(jsonlite); library(dplyr); library(readr)
})

OUT <- here::here("results", "external")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

# Restrict to nutrition_anaemia_* indicators that report prevalence (%)
INDICATORS <- c(
  women_repro_age = "NUTRITION_ANAEMIA_REPRODUCTIVEAGE_PREV",
  women_pregnant  = "NUTRITION_ANAEMIA_PREGNANT_PREV",
  women_nonpreg   = "NUTRITION_ANAEMIA_NONPREGNANT_PREV",
  children_6_59m  = "NUTRITION_ANAEMIA_CHILDREN_PREV"
)

TARGET_ISO3 <- c(
  "GMB", "GHA", "SLE", "MWI",     # training cohort
  "CIV",                           # external (in-flight)
  "KEN", "TZA", "ETH"              # future-extension candidates
)

fetch <- function(code) {
  url <- sprintf("https://ghoapi.azureedge.net/api/%s", code)
  cat(sprintf("[gho] fetching %s\n", code))
  resp <- httr::GET(url)
  httr::stop_for_status(resp)
  j <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                           simplifyDataFrame = TRUE)
  j$value
}

all_rows <- list()
for (i in seq_along(INDICATORS)) {
  code  <- INDICATORS[i]
  label <- names(INDICATORS)[i]
  df    <- fetch(code)
  df$indicator_label <- label
  df$indicator_code  <- code
  all_rows[[label]]  <- df
}

raw <- dplyr::bind_rows(all_rows)
cat(sprintf("[gho] raw rows: %d, cols: %d\n", nrow(raw), ncol(raw)))
cat("columns sample:\n")
print(colnames(raw))

# Filter to our target countries (SpatialDim = ISO3 for COUNTRY rows)
# Drop Dim1 strata that aren't the all-population total. Anaemia indicators are
# stratified by ANAEMIA severity (mild/moderate/severe/total). We want "Total".
target_rows <- raw |>
  dplyr::filter(SpatialDimType == "COUNTRY",
                SpatialDim %in% TARGET_ISO3,
                # Dim1 stratifies by anaemia severity (MILD/MODERATE/SEVERE/TOTAL)
                # and (for children) by sex (SEX_BTSX = both sexes). The
                # "any anaemia" total is reported under both SEX_BTSX and
                # SEVERITY_TOTAL — keep only those.
                is.na(Dim1) | Dim1 %in% c("SEX_BTSX", "SEVERITY_TOTAL")) |>
  dplyr::transmute(
    iso3       = SpatialDim,
    year       = TimeDim,
    indicator  = indicator_label,
    indicator_code,
    pop_dim    = Dim1,
    value_num  = NumericValue,
    value_str  = Value,
    value_low  = Low,
    value_high = High
  ) |>
  dplyr::arrange(iso3, indicator, dplyr::desc(year))
cat(sprintf("[gho] target-country rows: %d\n", nrow(target_rows)))

readr::write_csv(target_rows, file.path(OUT, "who_anaemia_country.csv"))
cat(sprintf("[gho] wrote %s\n", file.path(OUT, "who_anaemia_country.csv")))

# Latest year per country × indicator: quick coverage check
latest <- target_rows |>
  dplyr::group_by(iso3, indicator) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(iso3, indicator, year, value_num, value_low, value_high)
cat("\n[gho] latest-year prevalence per country×indicator:\n")
print(as.data.frame(latest))
readr::write_csv(latest, file.path(OUT, "who_anaemia_country_latest.csv"))
cat(sprintf("[gho] wrote %s\n", file.path(OUT, "who_anaemia_country_latest.csv")))
