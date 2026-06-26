# =============================================================================
# scripts/build_kenya_benchmark.R
#
# Build a Kenya NATIONAL benchmark table combining all currently-available
# external reference estimates. This is the "hypothetical-Kenya sensitivity"
# referenced in the manuscript: we can compare our model's national-mean
# prediction (when Kenya is added in v2) against KNMS 2011, WHO 2019, VMNIS.
#
# KNMS 2011 does NOT release admin-1 (county) estimates because the sample
# was not powered for them (Kenya has 47 counties; KNMS sampled by 4 broad
# ecological zones). Subnational Kenya validation therefore requires either
# (a) the full KNMS microdata via DUA, or (b) a different survey (KDHS 2022).
#
# Input files:
#   results/external/knms2011_national_summary.csv  (extracted by hand from PDF)
#   results/external/who_anaemia_country_latest.csv (WHO GHO API)
#   results/external/vmnis_national.csv             (WHO VMNIS .dta)
#
# Output:
#   results/external/kenya_national_benchmark.csv
# =============================================================================
suppressPackageStartupMessages({library(here); library(dplyr); library(readr)})
OUT <- here::here("results", "external")

knms <- readr::read_csv(file.path(OUT, "knms2011_national_summary.csv"),
                        show_col_types = FALSE)
who  <- readr::read_csv(file.path(OUT, "who_anaemia_country_latest.csv"),
                        show_col_types = FALSE)
vmn  <- readr::read_csv(file.path(OUT, "vmnis_national.csv"),
                        show_col_types = FALSE)

# Restrict each source to Kenya
who_kenya <- who |> dplyr::filter(iso3 == "KEN") |>
  dplyr::transmute(
    population = dplyr::recode(indicator,
      children_6_59m  = "preschool_6_59m",
      women_pregnant  = "pregnant_women",
      women_nonpreg   = "non_pregnant_women_15_49y"),
    outcome        = "anaemia",
    prevalence_pct = value_num,
    ci_low         = value_low,
    ci_high        = value_high,
    year           = year,
    source         = "WHO GHO 2019 (UNICEF/WHO joint estimate)"
  )

vmn_kenya <- vmn |> dplyr::filter(country == "Kenya") |>
  dplyr::transmute(
    population = dplyr::recode(Population,
      `Children (12-23 months)`               = "preschool_6_59m",
      `Children (24-59 months)`               = "preschool_6_59m",
      `Preschool-age children (6-59 months)`  = "preschool_6_59m",
      `Pregnant women (15-49 years)`          = "pregnant_women",
      `Pregnant women`                        = "pregnant_women",
      `Non-pregnant women (15-49 years)`      = "non_pregnant_women_15_49y",
      `Non-pregnant women of reproductive age`= "non_pregnant_women_15_49y",
      `Non-pregnant women`                    = "non_pregnant_women_15_49y",
      `School-age children (5-14 years)`      = "school_age_5_14y",
      `Adolescent girls (10-19 years)`        = "non_pregnant_women_15_49y",
      `Women of reproductive age (15-49 years)` = "non_pregnant_women_15_49y",
      .default = as.character(Population)),
    outcome        = dplyr::recode(mn_group,
      Haemoglobin    = "anaemia",
      Iron           = "iron_deficiency",
      `Vitamin A`    = "vitamin_a_deficiency",
      Zinc           = "zinc_deficiency",
      Folate         = "folate_deficiency",
      `Vitamin B12`  = "b12_deficiency",
      Iodine         = "iodine_deficiency",
      .default       = mn_group),
    # VMNIS reports several deficiency-prevalence columns alongside mean/median.
    # Coalesce the deficiency-type prevalences (rows that carry only a mean/
    # median biomarker concentration legitimately stay NA). 2026 fix: this was
    # previously hardcoded NA, discarding usable prevalences (e.g. Kenya folate
    # deficiency 31.5%/32.1% from KNMS 2011 via VMNIS).
    prevalence_pct = dplyr::coalesce(
      suppressWarnings(as.numeric(Prevalenceofdeficiency)),
      suppressWarnings(as.numeric(Prevalenceoflowvalues)),
      suppressWarnings(as.numeric(Prevalenceofinsufficiency)),
      suppressWarnings(as.numeric(Prevalenceofinadequacy))),
    ci_low         = NA_real_,
    ci_high        = NA_real_,
    year           = year,
    source         = sprintf("WHO VMNIS (%s, year=%d)", Indicator, year),
    Mean = Mean
  ) |> dplyr::select(-Mean)

knms_kenya <- knms |>
  dplyr::mutate(
    country = "Kenya",
    ci_low = NA_real_, ci_high = NA_real_, year = 2011L,
    source = "KNMS 2011 (national, no admin-1 release)"
  ) |>
  dplyr::select(population, outcome, prevalence_pct, ci_low, ci_high, year, source)

bench <- dplyr::bind_rows(knms_kenya, who_kenya, vmn_kenya) |>
  dplyr::arrange(population, outcome, dplyr::desc(year))

readr::write_csv(bench, file.path(OUT, "kenya_national_benchmark.csv"))
cat(sprintf("[bench] wrote %s (%d rows)\n",
            file.path(OUT, "kenya_national_benchmark.csv"), nrow(bench)))

# Pretty summary
cat("\n[bench] Kenya benchmarks by population × outcome (latest year per source):\n")
print(bench |>
        dplyr::filter(!is.na(prevalence_pct)) |>
        dplyr::group_by(population, outcome) |>
        dplyr::summarise(
          n_sources = dplyr::n(),
          values = paste(sprintf("%.1f%% (%s)", prevalence_pct,
                                  substr(source, 1, 20)),
                         collapse = "; "),
          .groups = "drop") |>
        as.data.frame())
