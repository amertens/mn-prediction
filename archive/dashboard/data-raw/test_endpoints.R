# Smoke test endpoints by simulating reactive flows.
# This catches errors in module servers without launching the full UI.
options(shiny.testmode = TRUE)
setwd(here::here("dashboard"))
source("global.R")

# 1. Map module — exercise data flow
cat("--- Map explorer reactive flow ---\n")
df <- get_country_admin2("ghana", "women_iron",
                          admin2_bnds, admin2_pred, admin2_pop)
stopifnot(nrow(df) > 0)
stopifnot("pred_prev" %in% colnames(df))
stopifnot("who_class" %in% colnames(df))
cat(sprintf("  OK: %d rows, predictions in [%.3f, %.3f]\n",
            nrow(df), min(df$pred_prev, na.rm=TRUE),
            max(df$pred_prev, na.rm=TRUE)))

# 2. National burden — verify aggregate outputs
cat("--- National burden ---\n")
for (ctry in names(meta$countries)) {
  for (oc in names(meta$outcome_labels)) {
    d <- get_country_admin2(ctry, oc, admin2_bnds, admin2_pred, admin2_pop)
    if (is.null(d) || all(is.na(d$pred_prev))) next
    natl <- national_aggregate(d)
    if (!is.na(natl$pred_prev_natl)) {
      stopifnot(natl$pred_prev_natl >= 0 && natl$pred_prev_natl <= 1)
    }
  }
}
cat("  OK: all national aggregates in valid [0,1] range\n")

# 3. Confidence badges
cat("--- Confidence badges ---\n")
stopifnot(confidence_badge(0.05) == "High confidence")
stopifnot(confidence_badge(0.15) == "Moderate confidence")
stopifnot(confidence_badge(0.25) == "Low confidence")
stopifnot(confidence_badge(NA) == "Unknown")
cat("  OK\n")

# 4. WHO classifications
cat("--- WHO classifications ---\n")
classes <- table(admin2_pred$who_class, useNA = "ifany")
print(classes)

cat("\nAll endpoint tests PASSED\n")
