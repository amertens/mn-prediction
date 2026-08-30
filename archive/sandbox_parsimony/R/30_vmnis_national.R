# =============================================================================
# sandbox_parsimony/R/30_vmnis_national.R
#
# A NATIONAL-level micronutrient model, built from national-level data.
#
# Motivation. FINDINGS.md sections 19-20: the surveys are designed and weighted
# for national estimates, not Admin-2, and the Admin-2 LOCO models recover
# national prevalence badly (best model ~5 pp mean absolute error, 2-3x out on
# iron). If the goal is an honest national prevalence for a country with no
# survey, the natural design is not "aggregate a district model" but "fit a
# model whose unit of observation IS the country-survey".
#
# WHO's Vitamin and Mineral Nutrition Information System (VMNIS) is the dataset
# for that: one row per survey x indicator x population group, across ~190
# countries and four decades. results/external/vmnis_national_full.csv holds
# only the 8 countries scripts/pull_vmnis_validation.R was pointed at; this
# reads the full long-format dump.
#
# SOURCE: WHO VMNIS long-format export, already on disk at
#   C:/Users/andre/OneDrive/Documents/mn-proxies/data/VMNIS/
#   VMNISIndicator_long_format.dta  (522 MB)
#
# Output: sandbox_parsimony/out/vmnis_national_slim.rds -- one row per
#   country x year x micronutrient x population, restricted to nationally
#   representative surveys, with the deficiency prevalence and sample size.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(haven)})

VMNIS <- "C:/Users/andre/OneDrive/Documents/mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta"
OUT   <- "sandbox_parsimony/out/vmnis_national_slim.rds"

KEEP <- c("Indicator", "CountryCode", "Country", "Representativeness",
          "Representativenessname", "Areacovered", "Beginyear", "Endyear",
          "Population", "Gender", "Agefrom", "Ageto", "Ageunit", "Samplesize",
          "Indicatorunit", "Mean", "Geomean", "Deficiencycutoff",
          "Prevalenceofdeficiency", "Dataadjustedfor", "Surveymethodology",
          "SurveyId", "Surveystatus")

# Same indicator -> micronutrient mapping the project already uses
# (scripts/pull_vmnis_validation.R).
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
  TRUE                                                                ~ NA_character_)

if (file.exists(OUT)) {
  slim <- readRDS(OUT)
  message("using cached ", OUT)
} else {
  stopifnot(file.exists(VMNIS))
  message("reading the 522 MB VMNIS dump (slow, once) ...")
  raw <- haven::read_dta(VMNIS, col_select = dplyr::all_of(KEEP))
  message(sprintf("  %d rows, %d countries", nrow(raw), dplyr::n_distinct(raw$Country)))

  slim <- raw |>
    mutate(across(where(haven::is.labelled), haven::as_factor),
           across(c(Beginyear, Endyear, Samplesize, Mean, Geomean,
                    Deficiencycutoff, Prevalenceofdeficiency),
                  ~ suppressWarnings(as.numeric(.x))),
           mn_group = mn_map(as.character(Indicator)),
           year = dplyr::coalesce(Endyear, Beginyear),
           country = as.character(Country),
           prev = Prevalenceofdeficiency / 100) |>
    filter(!is.na(mn_group), is.finite(prev), prev >= 0, prev <= 1,
           is.finite(year))
  saveRDS(slim, OUT)
  message("saved ", OUT)
}

cat("\n=== VMNIS, all rows with a deficiency prevalence ===\n")
cat(sprintf("rows %d | countries %d | years %d-%d\n", nrow(slim),
            dplyr::n_distinct(slim$country), min(slim$year), max(slim$year)))

# --- how nationally representative is each row? ----------------------------
cat("\n-- Representativeness --\n")
print(as.data.frame(slim |> count(rep = as.character(Representativenessname),
                                  sort = TRUE) |> head(8)))
cat("\n-- Areacovered --\n")
print(as.data.frame(slim |> count(area = as.character(Areacovered), sort = TRUE) |> head(6)))

# National = the `Representativeness` CODE says national (the companion
# `Representativenessname` column names the SUBNATIONAL unit and is blank
# exactly when the survey is national -- filtering on the name is backwards),
# and the survey covers both urban and rural. Anything else is a subnational or
# special-population survey and cannot serve as a national outcome.
nat <- slim |>
  filter(grepl("^national", trimws(as.character(Representativeness)), ignore.case = TRUE),
         grepl("both", as.character(Areacovered), ignore.case = TRUE))

cat(sprintf("\nnationally representative rows: %d across %d countries\n",
            nrow(nat), dplyr::n_distinct(nat$country)))

cat("\n=== usable country-years by micronutrient x population ===\n")
tab <- nat |>
  mutate(pop = as.character(Population)) |>
  filter(pop %in% c("Preschool-age children", "Non-pregnant women (NPW)",
                    "Women of reproductive age", "School-age children",
                    "Pregnant women")) |>
  group_by(mn_group, pop) |>
  summarise(rows = n(), countries = n_distinct(country),
            median_year = round(median(year)), .groups = "drop") |>
  arrange(desc(countries))
print(as.data.frame(tab |> filter(countries >= 8)), row.names = FALSE)

saveRDS(nat, "sandbox_parsimony/out/vmnis_national_rep.rds")
message("\nSaved -> sandbox_parsimony/out/vmnis_national_rep.rds")
