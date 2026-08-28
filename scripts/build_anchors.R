# =============================================================================
# scripts/build_anchors.R
#
# WS8c. Consolidate every national deficiency-prevalence anchor available for
# this project's countries and outcomes into one table, with provenance and an
# explicit statement of what is NOT available.
#
# WS7 needs a national anchor to constrain a transported map's
# population-weighted mean. Two arms:
#
#   own_survey  each held-out country's own design-based national value. Needs
#               no external data and measures the VALUE of anchoring, since it
#               is the best anchor that could exist.
#   external    a published national estimate. Measures anchor QUALITY, which
#               is the question that matters for a country with no survey.
#
# This script builds the external arm, and its main finding is about coverage.
#
# SOURCES ALREADY IN THE REPOSITORY
# ---------------------------------
# VMNIS was pulled by scripts/pull_vmnis_validation.R from a local WHO
# long-format extract held outside this repository
# (mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta), not from a machine
# endpoint. That is recorded as the reframe branch in the provenance table:
# there is no reachable API to re-pull from, so the loader validates a schema
# rather than fetching.
#
# Stevens et al. 2022 (scripts/pull_stevens2022.R) reports a COMPOSITE "any
# micronutrient deficiency" prevalence, not per-nutrient. It cannot anchor a
# single-nutrient map and is recorded with anchor_type "composite" so it cannot
# be joined into a per-nutrient anchor by accident.
#
# Run:
#   Rscript scripts/build_anchors.R
#
# Writes
#   metadata/anchors/national_anchors.csv
#   metadata/external_provenance.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
targets::tar_source(here("R"))

configs <- get_country_configs()
survey_year <- vapply(configs, function(cc) as.integer(cc$survey_year), integer(1))
names(survey_year) <- vapply(configs, function(cc) cc$country, character(1))
TARGET_COUNTRIES <- names(survey_year)

WOMEN_POPS <- c("Non-pregnant women (NPW)", "Women of reproductive age",
                "Non-pregnant, non-lactating women (NPNLW)")
CHILD_POPS <- c("Preschool-age children")

#' Map a VMNIS (mn_group, Population) pair to a project outcome tag.
#' Returns NA when the pair is not one this project models; pregnant women,
#' lactating women, men, adolescents and school-age children all fall here
#' because none of them is a modelled population.
.anchor_outcome <- function(mn_group, population) {
  nut <- dplyr::case_when(
    mn_group == "Vitamin A"   ~ "vitA",
    mn_group == "Iron"        ~ "iron",
    mn_group == "Folate"      ~ "folate",
    mn_group == "Vitamin B12" ~ "b12",
    mn_group == "Zinc"        ~ "zinc",
    TRUE ~ NA_character_)
  pop <- dplyr::case_when(
    population %in% CHILD_POPS ~ "child",
    population %in% WOMEN_POPS ~ "women",
    TRUE ~ NA_character_)
  ifelse(is.na(nut) | is.na(pop), NA_character_, paste0(pop, "_", nut))
}

# ── VMNIS ────────────────────────────────────────────────────────────────────
vmnis_path <- here("results", "external", "vmnis_national.csv")
if (!file.exists(vmnis_path))
  stop("VMNIS extract missing. Run scripts/pull_vmnis_validation.R first.", call. = FALSE)
v <- readr::read_csv(vmnis_path, show_col_types = FALSE)

vm <- v |>
  dplyr::filter(country %in% TARGET_COUNTRIES) |>
  dplyr::mutate(outcome = .anchor_outcome(mn_group, Population)) |>
  dplyr::filter(!is.na(outcome), !is.na(Prevalenceofdeficiency)) |>
  dplyr::transmute(
    country, outcome,
    anchor_prevalence = Prevalenceofdeficiency / 100,
    anchor_year = as.integer(year),
    survey_year = survey_year[country],
    year_gap = abs(as.integer(year) - survey_year[country]),
    anchor_type = "single_nutrient",
    population_label = Population,
    indicator = Indicator,
    source = "WHO VMNIS",
    source_url = "https://www.who.int/data/nutrition/nlis/info/vitamin-and-mineral-nutrition-information-system",
    source_note = paste("Derived by scripts/pull_vmnis_validation.R from a local WHO",
                        "long-format extract outside this repository; no machine",
                        "endpoint was used."),
    license = "WHO VMNIS terms of use apply; see the WHO nutrition data platform.",
    usable_as_anchor = TRUE
  )

# ── Stevens et al. 2022, composite ───────────────────────────────────────────
stevens <- list()
for (f in c(psc = "stevens2022_psc.csv", npw = "stevens2022_npw.csv")) {
  p <- here("results", "external", f)
  if (!file.exists(p)) next
  s <- readr::read_csv(p, show_col_types = FALSE)
  if (!all(c("whoname", "year", "any_def") %in% names(s))) next
  stevens[[f]] <- s |>
    dplyr::filter(whoname %in% TARGET_COUNTRIES) |>
    dplyr::transmute(
      country = whoname,
      outcome = if (grepl("psc", f)) "child_any_deficiency" else "women_any_deficiency",
      anchor_prevalence = any_def,
      anchor_year = as.integer(year),
      survey_year = survey_year[whoname],
      year_gap = abs(as.integer(year) - survey_year[whoname]),
      # Flagged so it can never be joined onto a single-nutrient outcome.
      anchor_type = "composite",
      population_label = if (grepl("psc", f)) "Preschool-age children" else "Non-pregnant women",
      indicator = "any micronutrient deficiency (composite)",
      source = "Stevens et al. 2022",
      source_url = "https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(22)00367-9/fulltext",
      source_note = "Composite across nutrients; cannot anchor a single-nutrient map.",
      license = "See the publication's terms.",
      usable_as_anchor = FALSE
    )
}

anchors <- dplyr::bind_rows(vm, dplyr::bind_rows(stevens)) |>
  dplyr::arrange(country, outcome, anchor_year)

dir.create(here("metadata", "anchors"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(anchors, here("metadata", "anchors", "national_anchors.csv"))

# ── Provenance ───────────────────────────────────────────────────────────────
prov <- dplyr::bind_rows(
  data.frame(
    artifact = "results/external/vmnis_national.csv",
    source = "WHO VMNIS",
    access_method = "local long-format extract, no machine endpoint reachable",
    url = "https://www.who.int/data/nutrition/nlis/info/vitamin-and-mineral-nutrition-information-system",
    retrieved_by = "scripts/pull_vmnis_validation.R",
    version_or_date = format(file.info(vmnis_path)$mtime, "%Y-%m-%d"),
    license = "WHO VMNIS terms of use",
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "results/external/stevens2022_psc.csv, results/external/stevens2022_npw.csv",
    source = "Stevens et al. 2022, Lancet Global Health",
    access_method = "publication supplementary data",
    url = "https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(22)00367-9/fulltext",
    retrieved_by = "scripts/pull_stevens2022.R",
    version_or_date = format(file.info(here("results", "external", "stevens2022_psc.csv"))$mtime, "%Y-%m-%d"),
    license = "See publication",
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "metadata/anchors/national_anchors.csv",
    source = "derived",
    access_method = "scripts/build_anchors.R",
    url = NA_character_,
    retrieved_by = "scripts/build_anchors.R",
    version_or_date = format(Sys.Date(), "%Y-%m-%d"),
    license = "inherits from the sources above",
    stringsAsFactors = FALSE)
)
readr::write_csv(prov, here("metadata", "external_provenance.csv"))

# ── Coverage report: what is missing is the point ────────────────────────────
PRIMARY <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
grid <- expand.grid(country = TARGET_COUNTRIES, outcome = PRIMARY,
                    stringsAsFactors = FALSE)
have <- anchors |> dplyr::filter(usable_as_anchor)
grid$has_anchor <- paste(grid$country, grid$outcome) %in%
  paste(have$country, have$outcome)
grid$year_gap <- vapply(seq_len(nrow(grid)), function(i) {
  h <- have[have$country == grid$country[i] & have$outcome == grid$outcome[i], ]
  if (!nrow(h)) NA_integer_ else min(h$year_gap)
}, numeric(1))

cat(sprintf("\nanchors written: %d rows (%d usable single-nutrient)\n",
            nrow(anchors), sum(anchors$usable_as_anchor)))
cat("\n=== coverage of the four primary outcomes ===\n")
print(as.data.frame(grid), row.names = FALSE)
cat(sprintf("\ncovered: %d of %d primary country-outcome cells\n",
            sum(grid$has_anchor), nrow(grid)))
cat("\nwrote metadata/anchors/national_anchors.csv and metadata/external_provenance.csv\n")
