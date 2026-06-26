# =============================================================================
# scripts/pull_vmnis_validation.R
#
# Build a focused VMNIS national-prevalence table for manuscript validation.
# Pulls from the local WHO VMNIS long-format dump at
#   C:/Users/andre/OneDrive/Documents/mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta
# (no internet/scraping needed; the file is already in the mn-proxies cache).
#
# Output:
#   results/external/vmnis_national.csv
#     country, year, indicator, mn_group, Population, prevalence, mean, n
#     for: Gambia, Ghana, Sierra Leone, Malawi, Cote d'Ivoire,
#          Kenya, Tanzania, Ethiopia
#     restricted to nationally-representative, both urban+rural surveys
#
# This is held-out ground truth for the LOCO / external-validation comparison.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(haven); library(tidyr); library(readr)
})

VMNIS_PATH <- "C:/Users/andre/OneDrive/Documents/mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta"
stopifnot(file.exists(VMNIS_PATH))

OUT_DIR <- here::here("results", "external")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TARGET_COUNTRIES <- c(
  "Gambia", "Ghana", "Sierra Leone", "Malawi",         # Training cohort
  "Cote d'Ivoire", "Côte d'Ivoire",                    # External (in-flight)
  "Kenya", "United Republic of Tanzania",
  "Tanzania (United Republic of)", "Tanzania",
  "Ethiopia"                                            # Future-extension candidates
)

# Same mn_group mapping as national_prediction/R/data_cleaning.R::clean_vmnis_for_heatmap
mn_map <- function(x) dplyr::case_when(
  x %in% c("Ferritin", "Serum transferrin receptor",
           "Iron (serum or plasma)", "Total iron binding capacity (TIBC)",
           "Transferrin saturation", "Mean corpuscular volume (MCV)") ~ "Iron",
  x %in% c("Retinol (plasma or serum)", "Beta-carotene")              ~ "Vitamin A",
  x == "Zinc (plasma or serum)"                                       ~ "Zinc",
  x %in% c("Goitre by ultrasound", "Iodine")                          ~ "Iodine",
  x %in% c("Folate (plasma or serum)", "Folate (red blood cell)")     ~ "Folate",
  x == "Vitamin B12"                                                  ~ "Vitamin B12",
  x == "25-Hydroxyvitamin D"                                          ~ "Vitamin D",
  x == "Haemoglobin"                                                  ~ "Haemoglobin",
  TRUE                                                                ~ NA_character_
)

cat("[vmnis] reading long-format file (~520 MB)...\n")
d <- haven::read_dta(VMNIS_PATH)
cat(sprintf("[vmnis] rows: %d, countries: %d\n",
            nrow(d), length(unique(d$Country))))

# Confirm which of our target-name spellings actually appear in the file
seen <- intersect(TARGET_COUNTRIES, unique(d$Country))
cat(sprintf("[vmnis] target countries present: %s\n", paste(seen, collapse = "; ")))

cleaned <- d |>
  dplyr::filter(
    Country %in% TARGET_COUNTRIES,
    Representativeness == "national",
    Areacovered == "both urban and rural"
  ) |>
  dplyr::mutate(mn_group = mn_map(Indicator)) |>
  tidyr::drop_na(mn_group) |>
  dplyr::select(
    country         = Country,
    year            = Beginyear,
    Indicator,
    mn_group,
    Population,
    Mean,
    Geomean,
    Median,
    dplyr::starts_with("Prevalence"),
    Dataadjustedfor,
    Surveymethodology,
    Agefrom, Ageto
  )

cat(sprintf("[vmnis] rows after filter (national, both U+R, target countries, mapped mn_group): %d\n",
            nrow(cleaned)))

# Per-(country, mn_group, Population) keep the most recent survey
latest <- cleaned |>
  dplyr::group_by(country, mn_group, Population) |>
  dplyr::arrange(dplyr::desc(year)) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

cat(sprintf("[vmnis] latest-per-(country,group,pop) rows: %d\n", nrow(latest)))

# Wide summary of latest prevalence by country × mn_group × Population
summary_table <- latest |>
  dplyr::select(country, year, mn_group, Population,
                Indicator, Mean, Geomean, Median,
                dplyr::starts_with("Prevalence"))

out_full   <- file.path(OUT_DIR, "vmnis_national_full.csv")
out_latest <- file.path(OUT_DIR, "vmnis_national.csv")
readr::write_csv(cleaned,       out_full)
readr::write_csv(summary_table, out_latest)
cat(sprintf("[vmnis] wrote %s (%d rows)\n", out_full,   nrow(cleaned)))
cat(sprintf("[vmnis] wrote %s (%d rows)\n", out_latest, nrow(summary_table)))

# Quick console summary by country × mn_group: how much do we have?
coverage <- summary_table |>
  dplyr::count(country, mn_group, name = "n_pop_groups") |>
  tidyr::pivot_wider(names_from = mn_group, values_from = n_pop_groups,
                     values_fill = 0)
cat("\n[vmnis] coverage matrix (count of Population groups per country×mn_group):\n")
print(as.data.frame(coverage))
