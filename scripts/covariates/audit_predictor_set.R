# =============================================================================
# scripts/covariates/audit_predictor_set.R
#
# One-table audit of the shared Admin-2 predictor set, for sign-off before a
# model re-run. Every column of predictors_admin2_shared.csv gets a row with
# its provenance (source, domain), coverage (which countries, finite share per
# country), whether it varies within a country, its scale, and the flags a
# reviewer needs: national constant, partial coverage, all-missing in some
# country, value-identical duplicate, external survey-derived, modelled
# outcome-adjacent surface, dropped-by-policy.
#
#   Rscript -e "source('scripts/covariates/audit_predictor_set.R')"
# -> results/tables/predictor_audit_<date>.csv          one row per column
# -> results/tables/predictor_audit_<date>_by_source.csv one row per source
# -> results/tables/predictor_audit_<date>_summary.md
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")   # assign_tier_v2()
HDIR <- "data/covariates/harmonized"; ODIR <- "results/tables"
STAMP <- format(Sys.Date(), "%Y-%m-%d")

S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
M <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
M$subnational <- as.logical(toupper(as.character(M$subnational)))
cols <- setdiff(names(S), c("country", "Admin1", "Admin2"))
COUNTRIES <- sort(unique(S$country))
stopifnot(setequal(cols, M$column))
cat(sprintf("[audit] %d Admin-2 rows, %d predictors, %d countries\n", nrow(S), length(cols), length(COUNTRIES)))

# per-column x country finite share and within-country sd
fin <- sapply(COUNTRIES, function(cc) vapply(cols, function(v) mean(is.finite(S[[v]][S$country == cc])), 0))
sdc <- sapply(COUNTRIES, function(cc) vapply(cols, function(v) { z <- S[[v]][S$country == cc]; z <- z[is.finite(z)]
  if (length(z) > 1) stats::sd(z) else NA_real_ }, 0))
colnames(fin) <- paste0("finite_", COUNTRIES); colnames(sdc) <- paste0("sd_", COUNTRIES)

# value-identical duplicates (same values everywhere both are finite)
key <- vapply(cols, function(v) paste(signif(S[[v]], 7), collapse = "|"), "")
dup_of <- vapply(seq_along(cols), function(i) { j <- which(key == key[i]); j <- j[j < i]; if (length(j)) cols[j[1]] else "" }, "")

A <- data.frame(column = cols, stringsAsFactors = FALSE) |>
  left_join(M, by = "column") |>
  mutate(
    finite_all   = round(vapply(cols, function(v) mean(is.finite(S[[v]])), 0), 3),
    n_countries_finite = rowSums(fin > 0.5),
    within_country_var = rowSums(!is.na(sdc) & sdc > 0) > 0,
    min = vapply(cols, function(v) suppressWarnings(min(S[[v]], na.rm = TRUE)), 0),
    median = vapply(cols, function(v) suppressWarnings(stats::median(S[[v]], na.rm = TRUE)), 0),
    max = vapply(cols, function(v) suppressWarnings(max(S[[v]], na.rm = TRUE)), 0),
    duplicate_of = dup_of,
    flag_national_constant = !within_country_var,
    flag_partial_coverage  = n_countries_finite < length(COUNTRIES),
    flag_low_completeness  = finite_all < 0.7,
    flag_duplicate         = nzchar(dup_of),
    flag_survey_derived    = grepl("^dhs_", column),
    flag_modelled_surface  = if ("modelled_surface" %in% names(M)) as.logical(toupper(as.character(modelled_surface))) %in% TRUE else grepl("MODELLED SURFACE", domain, fixed = TRUE),
    flag_survey_year_matched = if ("year_offset_max_abs" %in% names(M)) !is.na(year_offset_max_abs) & year_offset_max_abs == 0 else grepl("_sy$|_t0$|_sy_", column) | grepl("^dhs_|^map_sy_|^fprice_|^fsec_|^ihme_", column),
    flag_year_offset_3plus   = if ("year_offset_max_abs" %in% names(M)) !is.na(year_offset_max_abs) & year_offset_max_abs >= 3 else NA,
    flag_declared_subnational_mismatch = if ("subnational" %in% names(M)) !is.na(subnational) & (subnational != within_country_var) else NA
  )
if (!"tier" %in% names(A)) A$tier <- assign_tier_v2(A)
A <- cbind(A, round(fin, 3), round(sdc, 4))
A$min[!is.finite(A$min)] <- NA; A$max[!is.finite(A$max)] <- NA; A$median[!is.finite(A$median)] <- NA

out_col <- file.path(ODIR, sprintf("predictor_audit_%s.csv", STAMP))
write.csv(A, out_col, row.names = FALSE)

BS <- A |> group_by(source) |>
  summarise(n_columns = n(), n_domains = n_distinct(domain), tier = paste(sort(unique(tier)), collapse = "|"),
            n_all_countries = sum(!flag_partial_coverage), n_subnational = sum(within_country_var),
            n_national_constant = sum(flag_national_constant), n_low_completeness = sum(flag_low_completeness),
            n_duplicates = sum(flag_duplicate), n_modelled_surface = sum(flag_modelled_surface),
            n_year_offset_3plus = if ("year_offset_max_abs" %in% names(A)) sum(flag_year_offset_3plus, na.rm = TRUE) else NA_integer_,
            max_year_offset = if ("year_offset_max_abs" %in% names(A)) suppressWarnings(max(year_offset_max_abs, na.rm = TRUE)) else NA_real_,
            across(all_of(paste0("finite_", COUNTRIES)), ~ round(mean(.x), 2)), .groups = "drop") |>
  arrange(desc(n_columns))
BS$max_year_offset[!is.finite(BS$max_year_offset)] <- NA
BT <- A |> group_by(tier) |> summarise(n_columns = n(), n_sources = n_distinct(source), n_subnational = sum(within_country_var),
                                       n_national_constant = sum(flag_national_constant), .groups = "drop")
out_src <- file.path(ODIR, sprintf("predictor_audit_%s_by_source.csv", STAMP))
write.csv(BS, out_src, row.names = FALSE)

md <- c(sprintf("# Predictor set audit, %s", STAMP), "",
        sprintf("`%s`: %d Admin-2 rows x %d predictors, %d countries (%s).", "predictors_admin2_shared.csv",
                nrow(S), length(cols), length(COUNTRIES), paste(COUNTRIES, collapse = ", ")), "",
        "| flag | columns |", "|---|---:|",
        sprintf("| present in all %d countries | %d |", length(COUNTRIES), sum(!A$flag_partial_coverage)),
        sprintf("| partial coverage (< %d countries) | %d |", length(COUNTRIES), sum(A$flag_partial_coverage)),
        sprintf("| national constants (no within-country variation) | %d |", sum(A$flag_national_constant)),
        sprintf("| completeness < 70%% | %d |", sum(A$flag_low_completeness)),
        sprintf("| value-identical duplicates | %d |", sum(A$flag_duplicate)),
        sprintf("| survey-derived (DHS) | %d |", sum(A$flag_survey_derived)),
        sprintf("| modelled outcome-adjacent surfaces (V2_DROP_MODELLED sensitivity) | %d |", sum(A$flag_modelled_surface)),
        sprintf("| data year >= 3 years from the survey in some country (TA-01) | %s |", if ("year_offset_max_abs" %in% names(A)) sum(A$flag_year_offset_3plus, na.rm = TRUE) else "not stamped"),
        sprintf("| declared subnational flag disagrees with the data | %s |", if ("subnational" %in% names(A)) sum(A$flag_declared_subnational_mismatch, na.rm = TRUE) else "not stamped"), "",
        "## By tier (TP-01: V2_PREDICTOR_TIERS; national constants dropped at fit time unless V2_KEEP_NATIONAL=1)", "",
        "| tier | columns | sources | subnational | national const. |", "|---|---:|---:|---:|---:|",
        sprintf("| %s | %d | %d | %d | %d |", BT$tier, BT$n_columns, BT$n_sources, BT$n_subnational, BT$n_national_constant), "",
        "## By source", "", "| source | tier | columns | domains | all countries | subnational | national const. | <70% complete | duplicates | modelled | offset >= 3 y | max offset |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        sprintf("| %s | %s | %d | %d | %d | %d | %d | %d | %d | %d | %s | %s |", BS$source, BS$tier, BS$n_columns, BS$n_domains, BS$n_all_countries,
                BS$n_subnational, BS$n_national_constant, BS$n_low_completeness, BS$n_duplicates, BS$n_modelled_surface,
                ifelse(is.na(BS$n_year_offset_3plus), "-", BS$n_year_offset_3plus), ifelse(is.na(BS$max_year_offset), "-", BS$max_year_offset)), "",
        "## Temporal alignment (TA-01): columns whose data year is >= 3 years from the survey", "",
        if ("year_offset_max_abs" %in% names(A) && any(A$flag_year_offset_3plus, na.rm = TRUE)) {
          d <- A[A$flag_year_offset_3plus %in% TRUE, ] |> group_by(alignment_rule) |> summarise(n = n(), example = column[1], year_used = year_used[1], max_offset = max(year_offset_max_abs), .groups = "drop") |> arrange(desc(max_offset))
          c("| rule | columns | example | year used | max offset |", "|---|---:|---|---|---:|", sprintf("| %s | %d | `%s` | %s | %d |", d$alignment_rule, d$n, d$example, d$year_used, d$max_offset))
        } else "- none", "",
        "## Partial-coverage columns", "",
        "| column | source | countries | completeness |", "|---|---|---|---:|",
        with(A[A$flag_partial_coverage, ], sprintf("| %s | %s | %s | %.2f |", column, source, countries, finite_all)), "",
        "## Duplicates", "", if (any(A$flag_duplicate)) with(A[A$flag_duplicate, ], sprintf("- `%s` duplicates `%s`", column, duplicate_of)) else "- none")
writeLines(md, file.path(ODIR, sprintf("predictor_audit_%s_summary.md", STAMP)))
cat(sprintf("-> %s\n-> %s\n-> %s\n", out_col, out_src, file.path(ODIR, sprintf("predictor_audit_%s_summary.md", STAMP))))
print(as.data.frame(BT), row.names = FALSE); print(as.data.frame(BS[, c("source", "tier", "n_columns", "n_national_constant", "n_year_offset_3plus", "max_year_offset")]), row.names = FALSE)
cat("DONE\n")
