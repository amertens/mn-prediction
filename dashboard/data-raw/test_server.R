# Server-side checks — run from the repo root:
#   Rscript dashboard/data-raw/test_server.R
# smoke_test.R proves the UI constructs and the data joins; this exercises the
# reactives that decide which numbers a user actually sees. Added when the
# prediction-model default was unified, because a UI that builds fine can still
# open on the wrong estimator.

owd <- setwd(here::here("dashboard"))
on.exit(setwd(owd), add = TRUE)
source("global.R")
library(shiny)

fails <- character(0)
note  <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- c(fails, msg)
}

cat("Default estimator\n")
note(DEFAULT_PRED_MODEL == "recipe",
     "DEFAULT_PRED_MODEL is the recipe layer (docs/AREA_LEVEL_RECIPE_SPEC.md)")
note(!is.null(admin2_recipe_pred) && nrow(admin2_recipe_pred) > 0,
     "the recipe layer is loaded and non-empty")
note(DEFAULT_PRED_MODEL %in% unname(pred_model_choices()),
     "the default appears in the selector choices")
note(identical(pred_model_data(DEFAULT_PRED_MODEL), admin2_recipe_pred),
     "pred_model_data() resolves the default to the recipe table")

cat("\nMap explorer and district profiles agree\n")
shiny::testServer(mod_map_explorer_server, args = list(id = "map"), {
  session$setInputs(country = "ghana", outcome = "women_iron",
                    admin_level = "admin2", layer = "pred_prev",
                    pred_model = DEFAULT_PRED_MODEL)
  note(identical(pred_df(), admin2_recipe_pred),
       "map explorer opens on the recipe layer")
  map_val <<- as.data.frame(map_data())
})
shiny::testServer(mod_district_profile_server, args = list(id = "district"), {
  session$setInputs(country = "ghana", scope = "district",
                    dp_model = DEFAULT_PRED_MODEL)
  note(identical(dp_pred(), admin2_recipe_pred),
       "district profiles open on the same layer")
})

# The same district must not read differently on the two tabs.
d_map <- map_val[is.finite(map_val$pred_prev), c("Admin2", "pred_prev")]
d_dp  <- pred_model_data(DEFAULT_PRED_MODEL)
d_dp  <- d_dp[d_dp$country == "Ghana" & d_dp$outcome == "women_iron",
              c("Admin2", "pred_prev")]
m <- merge(d_map, d_dp, by = "Admin2", suffixes = c("_map", "_dp"))
note(nrow(m) > 0 && max(abs(m$pred_prev_map - m$pred_prev_dp), na.rm = TRUE) < 1e-9,
     sprintf("%d Ghana districts match between the two tabs", nrow(m)))

cat("\nScenario ranges\n")
shiny::testServer(mod_scenarios_server, args = list(id = "scenarios"), {
  session$setInputs(country = "ghana", outcome = "women_iron",
                    targeting = "all", coverage = 50, effect = 30, top_n = 10)
  b <- as.data.frame(base_data())
  note(all(c("prev_lo", "prev_hi") %in% names(b)),
       "base data carries a 95% prevalence range")
  note(sum(is.finite(b$prev_lo)) > 0,
       sprintf("%d of %d districts got a range", sum(is.finite(b$prev_lo)), nrow(b)))
  note(all(b$prev_lo <= b$pred_prev + 1e-9, na.rm = TRUE) &&
       all(b$prev_hi >= b$pred_prev - 1e-9, na.rm = TRUE),
       "the range brackets the point estimate")

  s <- as.data.frame(scenario_data())
  lo <- sum(s$cases_averted_lo, na.rm = TRUE)
  pt <- sum(s$cases_averted,    na.rm = TRUE)
  hi <- sum(s$cases_averted_hi, na.rm = TRUE)
  note(lo <= pt + 1e-6 && pt <= hi + 1e-6,
       sprintf("cases averted %.0f <= %.0f <= %.0f", lo, pt, hi))
  note(is.finite(pt) && pt > 0, "cases averted is a positive, finite number")
})

cat("\nInterval-bearing layers are sane\n")
# A small-area fit that collapses toward zero still produces a tidy-looking
# table, so compare each layer's median district against the national survey
# figure. This is how the BYM2 collapse of 2026-08-27 was found, after it had
# already been deployed.
.ne <- read.csv(here::here("results", "tables", "national_estimates_all.csv"),
                stringsAsFactors = FALSE)
layer_health <- function(d) {
  s <- do.call(rbind, lapply(split(d, list(d$country, d$outcome), drop = TRUE), function(x)
    data.frame(country = x$country[1], outcome = x$outcome[1],
               med = median(x$pred_prev, na.rm = TRUE))))
  s$survey <- vapply(seq_len(nrow(s)), function(i) {
    r <- .ne[.ne$country == s$country[i] & .ne$outcome == s$outcome[i], ]
    if (nrow(r)) r$obs_prev[1] else NA_real_
  }, numeric(1))
  s$ratio <- s$med / s$survey
  sum(s$ratio < 0.25 | s$ratio > 4, na.rm = TRUE)
}
fh_bad <- layer_health(admin2_fh_pred)
note(fh_bad <= 2, sprintf("Fay-Herriot: %d of 24 cells off by >4x from the survey", fh_bad))
if (!is.null(admin2_bym2_pred)) {
  # Gates at 3, the count before the degenerate-SE guard was briefly applied
  # here and took it to 10. See the note in _build_bym2_layer.R.
  by_bad <- layer_health(admin2_bym2_pred)
  note(by_bad <= 3, sprintf("BYM2: %d of 24 cells off by >4x from the survey", by_bad))
}
# Whatever supplies scenario ranges must be one of the sane ones.
rng <- attach_prevalence_range(
  as.data.frame(get_country_admin2("ghana", "women_iron", admin2_bnds,
                                   pred_model_data(DEFAULT_PRED_MODEL), admin2_pop)),
  "Ghana", "women_iron")
w <- median(rng$prev_hi - rng$prev_lo, na.rm = TRUE)
note(is.finite(w) && w > 0.02,
     sprintf("scenario range is a usable width (median %.1f pp, not a collapsed one)", 100 * w))

cat("\nStart here\n")
shiny::testServer(mod_start_here_server,
                  args = list(id = "start", go_to = function(...) invisible(NULL)), {
  h <- output$hero
  note(!is.null(h) && nzchar(as.character(h$html %||% h)),
       "hero stat renders")
  g <- output$ghana_example
  note(!is.null(g) && nzchar(as.character(g$html %||% g)),
       "worked example renders")
})

if (length(fails)) {
  cat("\nFAILURES:\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("\nAll server checks passed.\n")
