# =============================================================================
# sandbox/12_audit_gee_data.R
#
# Quality audit of the existing GEE raster extracts. Looks for:
#   1. NA rates per variable per country (extraction failures, polygon edge
#      effects).
#   2. Year-suffix distribution (is the data actually multi-year, or is
#      everything from 2019 regardless of survey year?).
#   3. Constant / near-constant variables (zero discriminative power).
#   4. Cross-country availability: how many of the 152 "common" GEE
#      variables actually exist consistently across all four countries.
#   5. Outliers / Inf / impossible values.
#   6. Per-variable correlation with the outcomes — is anything in here
#      a strong single predictor that we'd want as a smoothed addition
#      to the spatial GAM?
#
# Output: sandbox/results/12_gee_audit.csv with per-variable diagnostics.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

countries <- c("Gambia","Ghana","SierraLeone","Malawi")
gee_list <- setNames(lapply(c("gambia","ghana","sierraleone","malawi"), function(c)
  targets::tar_read_raw(paste0("gee_admin2_", c), store = STORE)),
  countries)

cat("===== GEE extract shapes =====\n")
for (n in names(gee_list)) cat(sprintf("  %s: %d polygons x %d cols\n",
                                         n, nrow(gee_list[[n]]),
                                         ncol(gee_list[[n]])))

# Common variables across all four
all_vars <- Reduce(intersect, lapply(gee_list, colnames))
cat("\ncommon vars across all 4 countries:", length(all_vars), "\n")

# Variable categorization by data-source prefix
src_prefix <- sub("^(gee_)?", "", all_vars)
src <- sub("_.*$", "", src_prefix)
cat("\n===== variables by data-source prefix =====\n")
print(table(src))

# Year distribution
year_suffix <- regmatches(all_vars, regexpr("_(20[01][0-9]|2020|2021|2022)", all_vars))
year_suffix <- sub("^_", "", year_suffix)
cat("\n===== year-suffixed variables =====\n")
cat("variables with year suffix:", sum(nchar(year_suffix) > 0), "/", length(all_vars), "\n")
print(table(year_suffix[nchar(year_suffix) > 0]))

# ── Per-variable diagnostics across countries ────────────────────────────
audit <- data.frame(
  variable = all_vars, na_frac = NA_real_, sd_value = NA_real_,
  constant_in_some_country = FALSE, has_inf = FALSE,
  range_min = NA_real_, range_max = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(all_vars)) {
  v <- all_vars[i]
  vals_all <- unlist(lapply(gee_list, function(g) g[[v]]))
  audit$na_frac[i] <- mean(!is.finite(vals_all))
  finite_vals <- vals_all[is.finite(vals_all)]
  audit$sd_value[i] <- if (length(finite_vals) > 1) stats::sd(finite_vals) else NA_real_
  audit$range_min[i] <- if (length(finite_vals)) min(finite_vals) else NA_real_
  audit$range_max[i] <- if (length(finite_vals)) max(finite_vals) else NA_real_
  audit$has_inf[i]   <- any(is.infinite(vals_all))
  audit$constant_in_some_country[i] <-
    any(vapply(gee_list, function(g) {
      x <- g[[v]]; xf <- x[is.finite(x)]
      length(xf) == 0 || length(unique(xf)) <= 1
    }, logical(1)))
}
cat("\n===== distribution of NA rates =====\n")
print(quantile(audit$na_frac, c(0, .1, .25, .5, .75, .9, .95, 1)))
cat("\nvars with NA frac > 0.50:", sum(audit$na_frac > 0.5), "\n")
cat("vars constant in some country:", sum(audit$constant_in_some_country), "\n")
cat("vars with Inf values:        ", sum(audit$has_inf), "\n")

readr::write_csv(audit, file.path(SANDBOX_RESULTS, "12_gee_audit.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "12_gee_audit.csv")))
