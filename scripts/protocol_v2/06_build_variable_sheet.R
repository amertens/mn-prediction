# =============================================================================
# scripts/protocol_v2/06_build_variable_sheet.R
#
# One metadata sheet for the modelling vocabulary, built to be an RA WORKSHEET.
#
# Every field is either DOCUMENTED (traced to a committed source) or
# AI_PROPOSED (derived here by rule or by inference and NOT verified). The
# `*_source` columns say which, for every row, so a reader never has to guess
# whether a definition came from a data dictionary or from a name pattern.
# That distinction is the point of the sheet: the taxonomy is load-bearing
# (the domain index is the best-transporting arm) and most of it has never
# been checked against source documentation.
#
# SOURCES, in priority order
#   1 data/covariates/harmonized/data_dictionary.csv   294 rows, rich and
#     genuinely documented: provider, pathway, canonical unit, temporal kind,
#     licence, source columns, collapse policy.
#   2 results/tables/predictor_inventory.csv           8383 rows, carries
#     conceptual_domain / conceptual_subdomain / description for 5978.
#   3 data/covariates/harmonized/predictors_admin2_shared_metadata.csv
#     the 18-domain taxonomy the domain scores actually use.
#   4 the data itself, for type, levels, IQR and per-country coverage.
#
#   Rscript scripts/protocol_v2/06_build_variable_sheet.R
# -> results/tables/protocol_v2/variable_sheet.csv
# -> results/tables/protocol_v2/variable_sheet_gaps.csv   (what an RA must fill)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SH  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
                check.names = FALSE)
SHM <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv",
                stringsAsFactors = FALSE)
HA  <- read.csv("data/covariates/harmonized/predictors_admin2_harmonized.csv",
                check.names = FALSE)
DD  <- read.csv("data/covariates/harmonized/data_dictionary.csv",
                stringsAsFactors = FALSE)
INV <- read.csv("results/tables/predictor_inventory.csv", stringsAsFactors = FALSE)

KEYS <- c("country", "Admin1", "Admin2")
shared_cols <- setdiff(names(SH), KEYS)
harm_cols   <- setdiff(names(HA), KEYS)
vars <- sort(union(shared_cols, harm_cols))
cat("vocabulary:", length(vars), "variables (",
    length(shared_cols), "shared /", length(harm_cols), "harmonized )\n")

# ── AI-PROPOSED sub-domain, by transparent rule on the name ─────────────────
# Rule: strip trailing measurement qualifiers (depth bands, months, years,
# statistic suffixes), then keep the leading two tokens. This groups
# soil_calcium_mean_0_20 and soil_calcium_stdev_20_50 under "soil_calcium",
# and lst_night_m03_t0 with lst_night_m08_t0 under "lst_night". It is a
# NAMING regularity, not a semantic judgement - hence ai_proposed.
propose_subdomain <- function(v) {
  x <- sub("_t[0-9]+$", "", v)
  x <- sub("_(m[0-9]{2}|y[0-9]{4}|[0-9]+_[0-9]+)$", "", x)
  x <- sub("_(mean|stdev|sd|min|max|median|sum|frac|share|prod|count|pct)$", "", x)
  tok <- strsplit(x, "_")[[1]]
  if (length(tok) >= 2) paste(tok[1], tok[2], sep = "_") else x
}

# ── empirical profile from the data ─────────────────────────────────────────
profile_var <- function(v) {
  src <- if (v %in% shared_cols) SH else HA
  x <- suppressWarnings(as.numeric(src[[v]]))
  ctry <- src$country
  ok <- is.finite(x)
  nd <- dplyr::n_distinct(x[ok])
  is_int <- all(abs(x[ok] - round(x[ok])) < 1e-9)
  vtype <- if (!any(ok)) "empty"
  else if (nd <= 2) "binary"
  else if (nd <= 10 && is_int) "categorical"
  else "continuous"
  q <- if (any(ok)) stats::quantile(x[ok], c(.25, .5, .75), names = FALSE)
  else rep(NA_real_, 3)
  # per-country presence and worst-country missingness
  pc <- tapply(ok, ctry, mean)
  present <- names(pc)[!is.na(pc) & pc > 0]
  const <- tapply(x, ctry, function(z) {
    z <- z[is.finite(z)]; length(z) > 0 && stats::sd(z) == 0
  })
  data.frame(
    var_type = vtype,
    n_levels = if (vtype %in% c("binary", "categorical")) nd else NA_integer_,
    median = if (vtype == "continuous") signif(q[2], 5) else NA_real_,
    IQR = if (vtype == "continuous") signif(q[3] - q[1], 5) else NA_real_,
    min = if (any(ok)) signif(min(x[ok]), 5) else NA_real_,
    max = if (any(ok)) signif(max(x[ok]), 5) else NA_real_,
    pct_missing_overall = round(100 * (1 - mean(ok)), 1),
    worst_country_pct_missing =
      if (length(pc)) round(100 * (1 - min(pc, na.rm = TRUE)), 1) else NA_real_,
    countries_present = paste(present, collapse = ";"),
    n_countries_present = length(present),
    constant_in_some_country = any(const, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

prof <- bind_rows(lapply(vars, profile_var))
prof$variable <- vars

# ── join documented sources ─────────────────────────────────────────────────
dd <- DD |> transmute(variable = canonical,
                      dd_provider = provider, dd_family = family,
                      dd_domain = domain, dd_definition = pathway,
                      dd_unit = canonical_unit, dd_temporal = temporal_kind,
                      dd_source_columns = source_columns,
                      dd_collapse = collapse_policy, dd_note = source_note)
inv <- INV |> filter(variable %in% vars) |>
  group_by(variable) |>
  slice(1) |> ungroup() |>
  transmute(variable, inv_source = source, inv_source_kind = source_kind,
            inv_domain = conceptual_domain, inv_subdomain = conceptual_subdomain,
            inv_description = description, inv_note = note)
shm <- SHM |> transmute(variable = column, assigned_domain = domain,
                        shm_source = source, shm_completeness = completeness,
                        shm_subnational = subnational)

sheet <- prof |>
  left_join(dd, by = "variable") |>
  left_join(inv, by = "variable") |>
  left_join(shm, by = "variable") |>
  mutate(
    in_shared_373    = variable %in% shared_cols,
    in_harmonized_294 = variable %in% harm_cols,
    country_specific = n_countries_present < 4,
    proposed_subdomain = vapply(variable, propose_subdomain, ""),
    # definition, with provenance
    definition = dplyr::coalesce(
      dplyr::na_if(dd_definition, "NA"),
      dplyr::na_if(inv_description, "NA")),
    definition_source = dplyr::case_when(
      !is.na(dplyr::na_if(dd_definition, "NA")) ~ "data_dictionary",
      !is.na(dplyr::na_if(inv_description, "NA")) ~ "predictor_inventory",
      TRUE ~ "MISSING_needs_RA"),
    data_origin = dplyr::coalesce(
      dplyr::na_if(dd_provider, "NA"), dplyr::na_if(inv_source, "NA"),
      dplyr::na_if(shm_source, "NA")),
    data_origin_source = dplyr::case_when(
      !is.na(dplyr::na_if(dd_provider, "NA")) ~ "data_dictionary",
      !is.na(dplyr::na_if(inv_source, "NA")) ~ "predictor_inventory",
      !is.na(dplyr::na_if(shm_source, "NA")) ~ "shared_metadata",
      TRUE ~ "MISSING_needs_RA"),
    unit = dplyr::na_if(dd_unit, "NA"),
    unit_source = ifelse(is.na(unit), "MISSING_needs_RA", "data_dictionary"),
    subdomain_source = "AI_PROPOSED_unverified",
    # blank columns for the RA to complete
    ra_verified_definition = NA_character_,
    ra_verified_subdomain = NA_character_,
    ra_mechanism_tag = NA_character_,
    ra_distal_proximal = NA_character_,
    ra_measurement_year = NA_character_,
    ra_notes = NA_character_)

sheet <- sheet |>
  select(variable, data_origin, data_origin_source,
         assigned_domain, dd_domain, inv_domain,
         proposed_subdomain, inv_subdomain, subdomain_source,
         definition, definition_source,
         unit, unit_source, dd_temporal, dd_collapse, dd_source_columns,
         var_type, n_levels, median, IQR, min, max,
         pct_missing_overall, worst_country_pct_missing,
         countries_present, n_countries_present, country_specific,
         constant_in_some_country, in_shared_373, in_harmonized_294,
         dd_note, inv_note,
         ra_verified_definition, ra_verified_subdomain, ra_mechanism_tag,
         ra_distal_proximal, ra_measurement_year, ra_notes) |>
  arrange(assigned_domain, proposed_subdomain, variable)

write.csv(sheet, file.path(OUTDIR, "variable_sheet.csv"), row.names = FALSE)

gaps <- sheet |>
  filter(in_shared_373,
         definition_source == "MISSING_needs_RA" | unit_source == "MISSING_needs_RA") |>
  select(variable, assigned_domain, proposed_subdomain, data_origin,
         definition_source, unit_source, var_type, countries_present)
write.csv(gaps, file.path(OUTDIR, "variable_sheet_gaps.csv"), row.names = FALSE)

cat("\n=== sheet:", nrow(sheet), "variables x", ncol(sheet), "fields ===\n")
cat("\ndefinition provenance (all rows):\n")
print(as.data.frame(sheet |> count(definition_source)), row.names = FALSE)
cat("\ndefinition provenance (the 373 actually modelled):\n")
print(as.data.frame(sheet |> filter(in_shared_373) |> count(definition_source)),
      row.names = FALSE)
cat("\nunit provenance (the 373):\n")
print(as.data.frame(sheet |> filter(in_shared_373) |> count(unit_source)),
      row.names = FALSE)
cat("\nvariable types (the 373):\n")
print(as.data.frame(sheet |> filter(in_shared_373) |> count(var_type)), row.names = FALSE)
cat("\nproposed sub-domains per assigned domain (the 373):\n")
print(as.data.frame(sheet |> filter(in_shared_373) |> group_by(assigned_domain) |>
  summarise(vars = n(), proposed_subdomains = n_distinct(proposed_subdomain),
            .groups = "drop") |> arrange(desc(vars))), row.names = FALSE)
cat("\nrows an RA must fill (in the modelled set):", nrow(gaps), "\n")
cat("\nDONE\n")
