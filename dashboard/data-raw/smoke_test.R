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
    # Difference layers must exist
    stopifnot(all(c("diff_prev", "loco_pred_prev", "loco_diff") %in% names(df)))
    natl <- national_aggregate(df)
    df1 <- get_country_admin1(ctry, oc, admin1_bnds, admin2_bnds,
                              admin2_pred, admin2_pop)
    stopifnot(all(c("diff_prev", "loco_diff") %in% names(df1)))
    cat(sprintf("  %s x %s: %d districts, %d preds, %d diff, %d LOCO, national=%.2f%%\n",
                ctry, oc, nrow(df), sum(!is.na(df$pred_prev)),
                sum(!is.na(df$diff_prev)), sum(!is.na(df$loco_diff)),
                100 * natl$pred_prev_natl))
  }
}

cat("\nAll combinations OK\n")
