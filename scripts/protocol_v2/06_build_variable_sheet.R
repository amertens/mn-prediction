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
# ONE SHEET, NOT TWO. Everything an RA needs is in a single table: the rows
# needing work are flagged, prioritised and sorted to the top, rather than
# split into a separate gaps file that can drift out of sync with the main one.
# The .xlsx carries real cell highlighting; the .csv carries the same
# information in FLAG_* columns for anything that reads it programmatically.
#
#   Rscript scripts/protocol_v2/06_build_variable_sheet.R
# -> results/tables/protocol_v2/variable_sheet.csv
# -> results/tables/protocol_v2/variable_sheet.xlsx   (highlighted worksheet)
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

# ── flags: what needs an RA's attention, and why ────────────────────────────
# A domain whose variables all propose the SAME sub-domain has no internal
# structure; one where every variable proposes its OWN has no grouping. Both
# mean the proposal is useless and the sub-domain must be authored by hand.
dom_stats <- sheet |> filter(in_shared_373) |> group_by(assigned_domain) |>
  summarise(dom_n = n(), dom_subdomains = n_distinct(proposed_subdomain),
            .groups = "drop") |>
  mutate(subdomain_proposal_useless = dom_n >= 6 &
           (dom_subdomains == 1 | dom_subdomains == dom_n))

sheet <- sheet |>
  left_join(dom_stats, by = "assigned_domain") |>
  mutate(
    FLAG_definition_missing   = definition_source == "MISSING_needs_RA",
    FLAG_unit_missing         = unit_source == "MISSING_needs_RA",
    FLAG_subdomain_useless    = tidyr::replace_na(subdomain_proposal_useless, FALSE),
    FLAG_domain_too_large     = tidyr::replace_na(dom_n >= 30, FALSE),
    FLAG_absent_in_a_country  = worst_country_pct_missing >= 100,
    FLAG_constant_in_country  = constant_in_some_country,
    FLAG_high_missingness     = pct_missing_overall >= 25,
    FLAG_empty                = var_type == "empty",
    n_flags = 0)
flagcols <- grep("^FLAG_", names(sheet), value = TRUE)
sheet$n_flags <- rowSums(as.matrix(sheet[, flagcols]), na.rm = TRUE)

sheet <- sheet |>
  mutate(
    ra_priority = dplyr::case_when(
      FLAG_empty ~ 1L,
      FLAG_definition_missing | FLAG_unit_missing ~ 1L,
      FLAG_subdomain_useless | FLAG_domain_too_large ~ 2L,
      FLAG_absent_in_a_country | FLAG_constant_in_country ~ 2L,
      FLAG_high_missingness ~ 3L,
      TRUE ~ 3L),
    ra_action = dplyr::case_when(
      FLAG_empty ~
        "EMPTY: no finite values anywhere. Confirm whether extraction failed, then drop or repair.",
      FLAG_definition_missing & FLAG_subdomain_useless ~
        "Find the provider's definition AND author a real sub-domain (the proposal is degenerate).",
      FLAG_definition_missing ~
        "Find the provider's definition, unit, measurement year and native resolution; cite the source.",
      FLAG_subdomain_useless ~
        "Author a real sub-domain: the name-rule proposal gives no usable grouping for this domain.",
      FLAG_domain_too_large ~
        "Large domain: split into coherent sub-domains of roughly 3-15 members.",
      FLAG_absent_in_a_country | FLAG_constant_in_country ~
        "Absent or constant in at least one country: confirm this is real, not an extraction failure. It is what the all-four-countries filter deletes on.",
      FLAG_high_missingness ~
        "High missingness: confirm the coverage is genuine before it is imputed.",
      TRUE ~
        "Verify the proposed sub-domain, then add mechanism tags and the distal-proximal rating."))

sheet <- sheet |>
  select(ra_priority, ra_action, n_flags,
         variable, data_origin, data_origin_source,
         assigned_domain, dom_n, dd_domain, inv_domain,
         proposed_subdomain, inv_subdomain, subdomain_source,
         definition, definition_source,
         unit, unit_source, dd_temporal, dd_collapse, dd_source_columns,
         var_type, n_levels, median, IQR, min, max,
         pct_missing_overall, worst_country_pct_missing,
         countries_present, n_countries_present, country_specific,
         constant_in_some_country, in_shared_373, in_harmonized_294,
         all_of(flagcols), dd_note, inv_note,
         ra_verified_definition, ra_verified_subdomain, ra_mechanism_tag,
         ra_distal_proximal, ra_measurement_year, ra_notes) |>
  arrange(ra_priority, desc(n_flags), assigned_domain, proposed_subdomain,
          variable)

write.csv(sheet, file.path(OUTDIR, "variable_sheet.csv"), row.names = FALSE)

# ── highlighted worksheet ───────────────────────────────────────────────────
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "variables")
  openxlsx::writeData(wb, "variables", sheet, withFilter = TRUE)
  nr <- nrow(sheet)
  cn <- function(x) which(names(sheet) == x)

  red    <- openxlsx::createStyle(bgFill = "#F4CCCC")  # must be filled
  amber  <- openxlsx::createStyle(bgFill = "#FCE5CD")  # must be checked
  blue   <- openxlsx::createStyle(bgFill = "#CFE2F3")  # RA writes here
  grey   <- openxlsx::createStyle(fgFill = "#EFEFEF", textDecoration = "italic")
  header <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9D9D9",
                                  halign = "left", border = "bottom")
  openxlsx::addStyle(wb, "variables", header, rows = 1,
                     cols = seq_along(sheet), gridExpand = TRUE)

  # red: the cells that are actually missing
  for (v in c("definition", "unit")) {
    src <- paste0(sub("^(.*)$", "\\1", v), "_source")
    openxlsx::conditionalFormatting(
      wb, "variables", cols = cn(v), rows = 2:(nr + 1),
      rule = sprintf('INDIRECT(ADDRESS(ROW(),%d))="MISSING_needs_RA"', cn(src)),
      style = red)
  }
  # amber: unverified proposals and coverage oddities the RA must confirm
  openxlsx::conditionalFormatting(
    wb, "variables", cols = cn("proposed_subdomain"), rows = 2:(nr + 1),
    rule = sprintf('INDIRECT(ADDRESS(ROW(),%d))=TRUE', cn("FLAG_subdomain_useless")),
    style = amber)
  for (v in c("countries_present", "constant_in_some_country",
              "pct_missing_overall", "worst_country_pct_missing")) {
    openxlsx::conditionalFormatting(
      wb, "variables", cols = cn(v), rows = 2:(nr + 1),
      rule = sprintf('OR(INDIRECT(ADDRESS(ROW(),%d))=TRUE,INDIRECT(ADDRESS(ROW(),%d))=TRUE)',
                     cn("FLAG_absent_in_a_country"), cn("FLAG_constant_in_country")),
      style = amber)
  }
  # blue: the columns the RA fills in
  openxlsx::addStyle(wb, "variables", blue,
                     rows = 2:(nr + 1),
                     cols = cn("ra_verified_definition"):cn("ra_notes"),
                     gridExpand = TRUE, stack = TRUE)
  # grey: machine-computed context, not to be edited
  openxlsx::addStyle(wb, "variables", grey, rows = 2:(nr + 1),
                     cols = cn("var_type"):cn("in_harmonized_294"),
                     gridExpand = TRUE, stack = TRUE)

  openxlsx::freezePane(wb, "variables", firstActiveRow = 2, firstActiveCol = 5)
  openxlsx::setColWidths(wb, "variables", cols = seq_along(sheet),
                         widths = "auto")
  openxlsx::setColWidths(wb, "variables",
                         cols = c(cn("ra_action"), cn("definition")), widths = 60)

  openxlsx::addWorksheet(wb, "legend")
  legend <- data.frame(
    colour = c("red", "amber", "blue", "grey", "", "priority 1", "priority 2",
               "priority 3"),
    meaning = c(
      "Missing. No documented value exists; find it from the provider and record it.",
      "Unverified or odd. A machine proposal with no usable structure, or a coverage pattern that decides whether the variable survives the all-four-countries filter. Confirm it is real.",
      "Yours to fill. The ra_* columns are the deliverable.",
      "Computed from the data. Context only - do not edit.",
      "",
      "No documented definition or unit, or the variable is empty.",
      "The sub-domain proposal is unusable, the domain is too large to be one score, or coverage needs confirming.",
      "Routine verification: check the proposed sub-domain, add mechanism tags and the distal-proximal rating."),
    stringsAsFactors = FALSE)
  openxlsx::writeData(wb, "legend", legend)
  openxlsx::addStyle(wb, "legend", header, rows = 1, cols = 1:2, gridExpand = TRUE)
  openxlsx::addStyle(wb, "legend", red,   rows = 2, cols = 1)
  openxlsx::addStyle(wb, "legend", amber, rows = 3, cols = 1)
  openxlsx::addStyle(wb, "legend", blue,  rows = 4, cols = 1)
  openxlsx::addStyle(wb, "legend", grey,  rows = 5, cols = 1)
  openxlsx::setColWidths(wb, "legend", cols = 1:2, widths = c(14, 110))

  openxlsx::saveWorkbook(wb, file.path(OUTDIR, "variable_sheet.xlsx"),
                         overwrite = TRUE)
  cat("wrote highlighted worksheet: variable_sheet.xlsx\n")
}
# the separate gaps file is retired: everything is in the one sheet now
unlink(file.path(OUTDIR, "variable_sheet_gaps.csv"))

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
cat("\nRA workload, by priority (the 373 modelled):\n")
print(as.data.frame(sheet |> filter(in_shared_373) |>
  count(ra_priority, name = "variables")), row.names = FALSE)
cat("\nflag counts (the 373):\n")
for (f in flagcols)
  cat(sprintf("  %-26s %d\n", f, sum(sheet[[f]][sheet$in_shared_373], na.rm = TRUE)))
cat("\nDONE\n")
