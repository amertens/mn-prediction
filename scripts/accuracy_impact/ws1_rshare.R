# =============================================================================
# scripts/accuracy_impact/ws1_rshare.R
#
# WS1 acceptance. Recompute r_share for every arm under the validated ceiling.
#
# Section 9 of docs/SESSION_FINDINGS_FOR_REVIEW.md reports mean r_share above 1
# for every anchoring arm, reaching 2.06, and leaves the cause open. WS1a, WS1b
# and WS1c establish that the analytic ceiling is biased low. This script
# restates the same measured correlations against the empirical ceiling so the
# size of the correction is visible, and reports the ratio as a MEDIAN as well
# as a mean because a ratio with a near-zero denominator has no usable mean.
#
# Nothing is overwritten: the output is a new file carrying a `ceiling` column.
#
#   Rscript scripts/accuracy_impact/ws1_rshare.R
# -> results/tables/r_share_revised.csv
# =============================================================================
suppressPackageStartupMessages({library(here)})
TDIR <- here("results", "tables")
arms <- read.csv(file.path(TDIR, "admin1_arms.csv"), stringsAsFactors = FALSE)
rel  <- read.csv(file.path(TDIR, "reliability_empirical.csv"), stringsAsFactors = FALSE)
rel  <- rel[rel$scheme == "within", ]

key <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
arms$.k <- paste(key(arms$country), arms$outcome)
rel$.k  <- paste(key(rel$country),  rel$outcome)
arms$r_max_emp <- rel$r_max_emp[match(arms$.k, rel$.k)]

# The same guard the production code applies: a ratio is only reported where the
# denominator is far enough from zero to carry meaning. add_reliability_columns()
# uses 0.05; the same cut is kept so the two are comparable.
ratio <- function(r, rm) ifelse(is.finite(rm) & rm > 0.05, r / rm, NA_real_)
arms$r_share_analytic  <- round(ratio(arms$pearson_r, arms$r_max), 3)
arms$r_share_empirical <- round(ratio(arms$pearson_r, arms$r_max_emp), 3)

out <- arms[, c("country","outcome","arm","n_admin2","pearson_r",
                "r_max","r_max_emp","r_share_analytic","r_share_empirical")]
readr::write_csv(out, file.path(TDIR, "r_share_revised.csv"))

summ <- do.call(rbind, lapply(split(out, out$arm), function(d) data.frame(
  arm = d$arm[1], cells = nrow(d),
  mean_r = round(mean(d$pearson_r, na.rm = TRUE), 3),
  n_analytic  = sum(is.finite(d$r_share_analytic)),
  mean_share_analytic = round(mean(d$r_share_analytic, na.rm = TRUE), 2),
  med_share_analytic  = round(stats::median(d$r_share_analytic, na.rm = TRUE), 2),
  n_empirical = sum(is.finite(d$r_share_empirical)),
  mean_share_empirical = round(mean(d$r_share_empirical, na.rm = TRUE), 2),
  med_share_empirical  = round(stats::median(d$r_share_empirical, na.rm = TRUE), 2),
  n_over_1_analytic  = sum(d$r_share_analytic  > 1, na.rm = TRUE),
  n_over_1_empirical = sum(d$r_share_empirical > 1, na.rm = TRUE),
  stringsAsFactors = FALSE)))
readr::write_csv(summ, file.path(TDIR, "r_share_revised_summary.csv"))

cat("=== r_share under the analytic ceiling and under the empirical one ===\n")
print(as.data.frame(summ), row.names = FALSE)
cat(sprintf("\ncells with a usable analytic denominator: %d of %d\n",
            sum(is.finite(out$r_share_analytic)), nrow(out)))
cat(sprintf("cells with a usable empirical denominator: %d of %d\n",
            sum(is.finite(out$r_share_empirical)), nrow(out)))
cat(sprintf("r_share > 1 under analytic: %d; under empirical: %d\n",
            sum(out$r_share_analytic > 1, na.rm = TRUE),
            sum(out$r_share_empirical > 1, na.rm = TRUE)))
