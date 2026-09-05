# =============================================================================
# scripts/protocol_v2/37_fieldwork_windows.R   [FW-01]
#
# WHEN WAS EACH CLUSTER SURVEYED? Fieldwork dates for time-matching covariates.
#
# Biomarkers are seasonal (ferritin and retinol fall with infection; serum
# zinc with fasting state and time of day) and food prices are seasonal, yet
# every covariate in the modelling vocabulary is an annual or multi-year
# summary. Matching climate and price covariates to the fieldwork window
# needs the date each cluster was visited. This script extracts it from
# whatever each survey carries, and reports what it could not find rather
# than assuming a window:
#   Gambia        gw_cIntDate / gw_wIntDate (Stata day counts, origin 1960-01-01)
#   Ghana         gw_month + gw_year (month resolution)
#   Malawi        date_interview in clean_malawi_mn_data.RDS, linked to the
#                 store's clusters by cluster number (34% of respondents lack
#                 a date; cluster medians are taken over those who have one)
#   Sierra Leone  no interview date. Reconstructed from child date of birth
#                 (month/year) + age in months where both exist; otherwise
#                 reported as missing.
#
#   Rscript scripts/protocol_v2/37_fieldwork_windows.R
# -> results/tables/protocol_v2/fieldwork_windows_cluster.csv
# -> results/tables/protocol_v2/fieldwork_windows_admin2.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; STORE <- "_targets_full"
rd <- function(name) tryCatch(tar_read_raw(name, store = STORE)$data, error = function(e) NULL)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
parts <- list()

# ── Gambia ───────────────────────────────────────────────────────────────────
for (on in c("child_vitA", "women_vitA", "child_iron", "women_iron")) { d <- rd(paste0("outcome_data_gambia_", on)); if (is.null(d)) next
  col <- if (grepl("^child", on)) "gw_cIntDate" else "gw_wIntDate"; if (!col %in% names(d)) next
  parts[[length(parts) + 1L]] <- data.frame(country = "Gambia", cluster = as.character(d$gw_cnum), Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2),
    date = as.Date(num(d[[col]]), origin = "1960-01-01"), src = on, stringsAsFactors = FALSE) }
# ── Ghana ────────────────────────────────────────────────────────────────────
for (on in c("child_vitA", "women_vitA")) { d <- rd(paste0("outcome_data_ghana_", on)); if (is.null(d)) next
  yr <- num(d$gw_year); if (all(is.na(yr))) yr <- rep(2017, nrow(d)); mo <- num(d$gw_month)
  cat("Ghana", on, "years:", paste(names(table(yr)), table(yr), collapse = " "), "| months:", paste(names(table(mo)), table(mo), collapse = " "), "\n")
  parts[[length(parts) + 1L]] <- data.frame(country = "Ghana", cluster = as.character(d$gw_cnum), Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2),
    date = as.Date(ifelse(is.finite(mo) & is.finite(yr), sprintf("%d-%02d-15", as.integer(yr), as.integer(mo)), NA)), src = on, stringsAsFactors = FALSE) }
# ── Malawi ───────────────────────────────────────────────────────────────────
m <- readRDS("data/IPD/Malawi/clean_malawi_mn_data.RDS")
st <- rd("outcome_data_malawi_child_vitA")
if (!is.null(st)) { a <- unique(as.character(m$cluster)); b <- unique(as.character(st$gw_cnum))
  cat(sprintf("Malawi cluster ids: clean %d | store %d | overlap %d\n", length(a), length(b), length(intersect(a, b)))) }
parts[[length(parts) + 1L]] <- data.frame(country = "Malawi", cluster = as.character(m$cluster), Admin1 = as.character(m$Admin1), Admin2 = as.character(m$Admin2),
  date = as.Date(m$date_interview), src = "clean_mn", stringsAsFactors = FALSE)
# ── Sierra Leone: reconstruct if possible ────────────────────────────────────
d <- rd("outcome_data_sierraleone_child_vitA")
if (!is.null(d)) {
  agecols <- names(d)[grepl("age", names(d), ignore.case = TRUE) & !grepl("^dhs|^ihme|^gee|cat|percent|village|average|coverage|usage|manage|storage|image", names(d), ignore.case = TRUE)]
  cat("Sierra Leone age-like columns:", paste(head(agecols, 20), collapse = ", "), "\n")
  for (h in head(agecols, 20)) { v <- num(d[[h]]); if (sum(is.finite(v)) > 0) cat(sprintf("   %-22s finite %4d  range %s..%s\n", h, sum(is.finite(v)), signif(min(v, na.rm = TRUE), 4), signif(max(v, na.rm = TRUE), 4))) }
  mo_col <- intersect(c("gw_cAgeCalcDays", "gw_cAgeMonths", "gw_cAgeMonth", "gw_age_months", "gw_cAgeMonths_GMNS"), names(d))
  dob_m <- num(d$gw_cDOBMonth); dob_y <- num(d$gw_cDOBYear)
  if (length(mo_col)) {
    # interview date = date of birth (mid-month) + age; gw_cAgeCalcDays is age in
    # days (370 of 486 children), the others age in months
    dob <- as.Date(ifelse(is.finite(dob_y) & is.finite(dob_m) & dob_m >= 1 & dob_m <= 12, sprintf("%d-%02d-15", as.integer(dob_y), as.integer(dob_m)), NA))
    age_days <- if (mo_col[1] == "gw_cAgeCalcDays") num(d[[mo_col[1]]]) else round(num(d[[mo_col[1]]]) * 30.44)
    dt <- dob + age_days
    dt[!is.na(dt) & (format(dt, "%Y") < "2012" | format(dt, "%Y") > "2014")] <- NA
    cat("Sierra Leone reconstructed dates from DOB +", mo_col[1], ": finite", sum(!is.na(dt)), "of", nrow(d), "| range", format(range(dt, na.rm = TRUE)), "\n")
    parts[[length(parts) + 1L]] <- data.frame(country = "SierraLeone", cluster = as.character(d$gw_cnum), Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2),
      date = dt, src = "reconstructed_dob_age", stringsAsFactors = FALSE)
  } else {
    cat("Sierra Leone: no age-in-months column; fieldwork dates NOT available\n")
    parts[[length(parts) + 1L]] <- data.frame(country = "SierraLeone", cluster = as.character(d$gw_cnum), Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2),
      date = as.Date(NA), src = "none", stringsAsFactors = FALSE)
  }
}
ALL <- bind_rows(parts); ALL <- ALL[!is.na(ALL$cluster), ]
CL <- ALL |> group_by(country, cluster) |> summarise(Admin1 = names(sort(table(Admin1), decreasing = TRUE))[1], Admin2 = names(sort(table(Admin2), decreasing = TRUE))[1],
  n_resp = dplyr::n(), n_dated = sum(!is.na(date)), date_med = as.Date(stats::median(date, na.rm = TRUE)), date_min = suppressWarnings(min(date, na.rm = TRUE)), date_max = suppressWarnings(max(date, na.rm = TRUE)), .groups = "drop")
CL$date_min[!is.finite(CL$date_min)] <- NA; CL$date_max[!is.finite(CL$date_max)] <- NA
CL$month_med <- as.integer(format(CL$date_med, "%m")); CL$year_med <- as.integer(format(CL$date_med, "%Y"))
write.csv(CL, file.path(OUTDIR, "fieldwork_windows_cluster.csv"), row.names = FALSE)
A2 <- CL |> filter(!is.na(date_med)) |> group_by(country, Admin1, Admin2) |> summarise(n_clusters = dplyr::n(), date_med = as.Date(stats::median(date_med)),
  date_first = min(date_med), date_last = max(date_med), .groups = "drop") |>
  mutate(month_med = as.integer(format(date_med, "%m")), year_med = as.integer(format(date_med, "%Y")),
         month_first = as.integer(format(date_first, "%m")), month_last = as.integer(format(date_last, "%m")),
         n_months = pmax(1L, as.integer(round(as.numeric(date_last - date_first) / 30.4)) + 1L))
write.csv(A2, file.path(OUTDIR, "fieldwork_windows_admin2.csv"), row.names = FALSE)
cat("\n===== FW-01: fieldwork windows =====\n")
print(as.data.frame(CL |> group_by(country) |> summarise(clusters = dplyr::n(), dated = sum(n_dated > 0), first = min(date_med, na.rm = TRUE), last = max(date_med, na.rm = TRUE),
  span_days = as.integer(max(date_med, na.rm = TRUE) - min(date_med, na.rm = TRUE)), .groups = "drop")), row.names = FALSE)
cat("\nclusters by calendar month:\n"); print(table(CL$country, CL$month_med, useNA = "ifany"))
cat("\ndistricts with dates:", nrow(A2), "| by country:", paste(names(table(A2$country)), table(A2$country), collapse = " "), "\n")
cat("\nDONE\n")
