# =============================================================================
# scripts/protocol_v2/24_sign_convention_check.R   [QA-01]
#
# REGRESSION TEST: THE SIGNAL-SCAN SIGN CONVENTION
#
# On 2026-09-02 the direction of the project's two strongest replicated
# associations was reported backwards in the manuscript, the slide outline,
# the dashboard and a published summary: legume consumption and cattle
# ownership were described as tracking LESS deficiency. In both scan outputs a
# POSITIVE meta_z means MORE deficiency (the continuous scan negates the
# biomarker, y = -wmean(t, w), precisely to match the binary scan), and both
# indicators are positive. This script asserts the convention against the
# result files so a future rerun that flips it fails loudly instead of
# propagating.
#
# It checks three things:
#   1. p4 header still documents the negation.
#   2. In p1 AND p4 the anchor indicators have the expected sign:
#        dhs_c_fg_legumes  ~ child_vitA   positive
#        dhs_hh_cattle     ~ child_iron   positive
#        spam_share_roots  ~ child_iron   negative
#        dhs_c_wasted      ~ women_vitA   positive   (p1 only; p4 lacks it)
#   3. The dashboard headline bundle carries "MORE deficiency" for those rows.
#
# Exit status is non-zero on any failure, so it can gate a rebuild.
#
#   Rscript scripts/protocol_v2/24_sign_convention_check.R
# =============================================================================
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
fails <- character(0)
note <- function(ok, msg) { cat(if (ok) "  PASS " else "  FAIL ", msg, "\n"); if (!ok) fails <<- c(fails, msg) }

hdr <- readLines("scripts/signal_probes/p4_admin1_continuous_scan.R", n = 40)
note(any(grepl("NEGATED before correlation", hdr)), "p4 header documents the outcome negation")
src <- readLines("scripts/signal_probes/p4_admin1_continuous_scan.R")
note(any(grepl("y = -wmean\\(t, w\\)", src)), "p4 code negates the biomarker (y = -wmean(t, w))")

expect <- data.frame(
  predictor = c("dhs_c_fg_legumes", "dhs_hh_cattle", "spam_share_roots", "dhs_c_wasted"),
  group     = c("child_vitA",       "child_iron",    "child_iron",       "women_vitA"),
  sign      = c(+1, +1, -1, +1), stringsAsFactors = FALSE)
for (f in c("p1_admin1_scan_predictors.csv", "p4_admin1_continuous_predictors.csv")) {
  d <- read.csv(file.path("results/tables/signal_probes", f), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(expect))) {
    r <- d[d$predictor == expect$predictor[i] & d$group == expect$group[i], ]
    if (!nrow(r)) { cat("  skip ", f, expect$predictor[i], "(absent)\n"); next }
    note(sign(r$meta_z[1]) == expect$sign[i],
         sprintf("%s: %s ~ %s meta_z = %+.2f has expected sign %+d", f, expect$predictor[i],
                 expect$group[i], r$meta_z[1], expect$sign[i]))
  }
}

b <- tryCatch(readRDS("dashboard/data/nutrient_signal.rds"), error = function(e) NULL)
if (!is.null(b)) {
  h <- b$headline
  note(grepl("MORE deficiency", h$direction[grepl("children", h$outcome) & grepl("Vitamin A", h$outcome)]),
       "dashboard headline: child vitamin A row says MORE deficiency")
  note(grepl("MORE deficiency", h$direction[grepl("children", h$outcome) & grepl("Iron", h$outcome)]),
       "dashboard headline: child iron row says MORE deficiency")
} else cat("  skip  dashboard bundle not built\n")

cat(sprintf("\n%d failure(s)\n", length(fails)))
if (length(fails)) quit(status = 1L)
