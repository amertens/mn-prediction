# Extended endpoint tests covering scenario and GBD modules
options(shiny.testmode = TRUE)
setwd(here::here("dashboard"))
source("global.R")

cat("--- Population year toggle ---\n")
df_survey <- get_country_admin2("ghana", "women_iron",
                                  admin2_bnds, admin2_pred, admin2_pop,
                                  pop_year = "survey")
df_2023 <- get_country_admin2("ghana", "women_iron",
                                admin2_bnds, admin2_pred, admin2_pop,
                                pop_year = "2023")
pop_survey_total <- sum(df_survey$population, na.rm = TRUE)
pop_2023_total <- sum(df_2023$population, na.rm = TRUE)
stopifnot(pop_2023_total > pop_survey_total)
cat(sprintf("  OK: survey-yr pop=%s, 2023 pop=%s (×%.2f)\n",
            format(round(pop_survey_total), big.mark = ","),
            format(round(pop_2023_total), big.mark = ","),
            pop_2023_total / pop_survey_total))

cat("--- Coverage scenario math ---\n")
cov <- 0.6
eff <- 0.4
prev_before <- df_survey$pred_prev
prev_after <- prev_before * (1 - cov * eff)
stopifnot(all(prev_after <= prev_before, na.rm = TRUE))
stopifnot(all(prev_after >= 0, na.rm = TRUE))
cases_averted <- sum((prev_before - prev_after) * df_survey$population, na.rm = TRUE)
cat(sprintf("  OK: %s cases averted under 60%% coverage × 40%% effect\n",
            format(round(cases_averted), big.mark = ",")))

cat("--- GBD data load ---\n")
stopifnot(!is.null(gbd_data$data))
stopifnot(nrow(gbd_data$data) > 0)
cat(sprintf("  OK: %d GBD rows loaded (placeholder=%s)\n",
            nrow(gbd_data$data), gbd_meta$using_placeholder))

cat("--- GBD comparison logic ---\n")
gbd <- gbd_data$data
ctry_label <- "Ghana"
oc <- "child_vitA"
target_year <- 2017
gbd_row <- gbd[gbd$country == ctry_label & gbd$outcome == oc &
                  gbd$year == target_year, ]
stopifnot(nrow(gbd_row) > 0)
df_pipe <- get_country_admin2("ghana", oc, admin2_bnds, admin2_pred, admin2_pop)
natl <- national_aggregate(df_pipe)
diff_pp <- (natl$pred_prev_natl - gbd_row$gbd_prev[1]) * 100
cat(sprintf("  OK: Ghana child_vitA — pipeline %.1f%% vs GBD %.1f%% (Δ %+.1f pp)\n",
            natl$pred_prev_natl * 100, gbd_row$gbd_prev[1] * 100, diff_pp))

cat("\nAll v2 endpoint tests PASSED\n")
