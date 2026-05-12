# Smoke test (run from dashboard/ working dir): exercise core helper paths
setwd(here::here("dashboard"))
source("global.R")

stopifnot(nrow(admin2_pred) > 0)
stopifnot(nrow(admin2_pop) > 0)
stopifnot(length(admin2_bnds) == 4)

# Each country × outcome combination should produce a valid sf
for (ctry in names(meta$countries)) {
  outcomes_for_ctry <- unique(admin2_pred$outcome[
    admin2_pred$country == meta$countries[ctry]
  ])
  for (oc in outcomes_for_ctry) {
    df <- get_country_admin2(ctry, oc, admin2_bnds, admin2_pred, admin2_pop)
    stopifnot(!is.null(df))
    stopifnot(nrow(df) > 0)
    natl <- national_aggregate(df)
    cat(sprintf("  %s × %s: %d districts, %d with predictions, national=%.2f%%\n",
                ctry, oc, nrow(df), sum(!is.na(df$pred_prev)),
                100 * natl$pred_prev_natl))
  }
}

cat("\nAll combinations OK\n")
