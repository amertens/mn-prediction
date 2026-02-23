
# ===== 0) Packages =====
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "janitor", "countrycode")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ===== 1) Read the Excel (country sheet) =====
# If you're running in the same environment where you uploaded the file, use the /mnt/data path.
# Otherwise, put the file in your working directory and use the short name.

path <- "C:/Users/andre/Downloads/UNICEF-WHO-Birthweight-data-coverage-2023.xlsx"
if (!file.exists(path)) path <- "UNICEF-WHO-Birthweight-data-coverage-2023.xlsx"

raw <- readxl::read_excel(
  path,
  sheet = "country birthweight coverage",
  col_names = FALSE
)

# ===== 2) Promote the embedded header row and drop empty columns =====
# Find the row whose first cell says "ISO3 code" (that is the header row)
hdr_row <- which(raw[[1]] == "ISO3 code")[1]
stopifnot(!is.na(hdr_row))

names(raw) <- as.character(unlist(raw[hdr_row, ]))
dat <- raw |>
  dplyr::slice((hdr_row + 1):dplyr::n()) |>
  dplyr::select(where(~ !all(is.na(.x))))  # drop completely empty columns

# Tidy column names and make duplicates unique (e.g., "UNICEF", "UNICEF_2")
names(dat) <- janitor::make_clean_names(names(dat), unique = TRUE)

# ===== 3) Rename the important variables (robust to minor header text changes) =====
# After clean_names(), your columns will look like:
# iso3_code, unicef, unicef_2, un_major, un_minor, sdg_region, country, ...
dat <- dat |>
  dplyr::rename(
    iso3               = iso3_code,
    unicef_region      = unicef,
    unicef_subregion   = unicef_2,
    un_major           = un_major,
    un_minor           = un_minor,
    sdg_region         = sdg_region,
    country            = country,
    # long/awkward headers (use regex-based matching to be safe)
    pct_without_weight = tidyselect::matches("^percentage_of_births_without_a_birthweight"),
    method_flag        = tidyselect::matches("^blank_as_reported"),
    survey_publication_year       = tidyselect::matches("^survey_publication_year1$"),
    survey_publication_year_short = tidyselect::matches("^survey_publication_year_short$"),
    outside_2015_2021             = tidyselect::matches("^outside_of_2015_2021"),
    year                          = tidyselect::matches("^year_assigned"),
    source_information            = tidyselect::matches("^source_information$"),
    source_short                  = tidyselect::matches("^short_source4$")
  )

# ===== 4) Type conversions & derived fields =====
dat <- dat |>
  dplyr::mutate(
    iso3   = as.character(iso3),
    country = as.character(country),
    year   = suppressWarnings(as.integer(year)),
    pct_without_weight = suppressWarnings(as.numeric(pct_without_weight)),
    # Coverage (% with birthweight recorded) = 100 - % without
    coverage_pct = dplyr::if_else(
      is.na(pct_without_weight), NA_real_,
      pmax(0, pmin(100, 100 - pct_without_weight))
    ),
    # Recode method flags
    method_flag = dplyr::case_when(
      is.na(method_flag) | method_flag == "" ~ "as_reported",
      method_flag == "R" ~ "reanalysis_survey",
      method_flag == "C" ~ "calculated_admin",
      TRUE ~ as.character(method_flag)
    ),
    # Treat any non-blank marker as TRUE; blank as FALSE
    outside_2015_2021 = !(is.na(outside_2015_2021) | trimws(outside_2015_2021) == "")
  ) |>
  dplyr::arrange(iso3)

# Optional: quick QA checks
# dplyr::count(dat, iso3) |> dplyr::filter(n > 1)     # duplicates?
# summary(dat$pct_without_weight); summary(dat$coverage_pct)

# ===== 5) Outputs for merging =====

# (A) Long (tidy) form: one row per country-year
birthweight_long <- dat |>
  dplyr::select(
    iso3, country, year,
    pct_without_weight,           # % births without a recorded weight
    coverage_pct,                 # % births WITH a recorded weight
    method_flag, source_short
  ) |>
  dplyr::filter(!is.na(year))     # drop rows without an assigned year

# (B) Wide (year columns) — choose the measure you want to merge.
# Here we use 'coverage_pct' because it reads naturally as "coverage".
# If you prefer the raw "% without birthweight", change values_from to pct_without_weight.
birthweight_wide <- birthweight_long |>
  dplyr::select(iso3, country, year, coverage_pct) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = coverage_pct
  ) |>
  # Make year column names syntactic (helpful in base R & some IDEs)
  dplyr::rename_with(~ paste0("Y", .x), tidyselect::matches("^\\d{4}$")) |>
  dplyr::arrange(iso3)

# ===== 6) Examples: merging into your wide country×years dataset =====

# Example A: your dataset has an ISO3 column already -------------------
# your_wide_df <- read.csv("your_wide_table.csv")
# merged <- your_wide_df |>
#   dplyr::left_join(birthweight_wide, by = "iso3")

# Example B: your dataset has country names (not ISO3) -----------------
# Build an iso3 column, then join by iso3. Tweak special cases if needed.
# your_wide_df <- read.csv("your_wide_table.csv")
# your_wide_df <- your_wide_df |>
#   dplyr::mutate(
#     iso3 = countrycode::countrycode(country, "country.name", "iso3c"),
#     iso3 = dplyr::case_when(
#       country %in% c("Kosovo", "Republic of Kosovo", "Kosovo (under UNSCR 1244)") ~ "XKX",
#       country %in% c("State of Palestine", "Palestine") ~ "PSE",
#       TRUE ~ iso3
#     )
#   )
# merged <- your_wide_df |>
#   dplyr::left_join(birthweight_wide, by = "iso3")

# ===== 7) (Optional) Save outputs =====
# write.csv(birthweight_long, "birthweight_coverage_long.csv", row.names = FALSE)
# write.csv(birthweight_wide, "birthweight_coverage_wide.csv", row.names = FALSE)
