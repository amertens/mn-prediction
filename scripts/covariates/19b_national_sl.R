# =============================================================================
# scripts/covariates/19b_national_sl.R   [NAT-SL, 2026-09-28]
#
# The national (WHO VMNIS) model of script 19 with a SuperLearner added, for
# the MNF15 slide "Can proxy models recover national-level prevalence?".
#
#   1. Leave-one-country-out on every VMNIS panel, as script 19 section 1, with
#      the SuperLearner (mean, ridge, lasso, elastic net, random forest; NNLS;
#      inner folds grouped by country) beside the null, ridge and forest arms,
#      all on the same outer folds, imputation and screen.
#   2. Our four countries' national vitamin A prevalence predicted with the
#      country held out of VMNIS, by ridge (script 19's learner) and by the
#      SuperLearner, at the survey years of metadata/survey_years.csv (script
#      19's table predates the 7 September year fix: The Gambia 2021, Malawi
#      2015). The survey value is the respondent-weighted mean of the district
#      prevalences in the current targets_v2.csv, as script 19 took it.
#
# Does not touch script 19's tables.
#
#   Rscript scripts/covariates/19b_national_sl.R
# -> results/tables/national_vmnis_loco_sl.csv         metrics, every arm
#    results/tables/national_vmnis_loco_sl_pred.csv    out-of-country predictions
#    results/tables/national_levels_sl.csv             our four countries, vitamin A
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(here); library(glmnet); library(ranger)
})
source(here("R", "national_covariates.R"))
source(here("R", "national_vmnis.R"))
source(here("R", "survey_years.R")); SURVEY_YEAR <- survey_years(keys = "lower")
set.seed(20260928L)

SCREEN_K <- as.integer(Sys.getenv("NAT_SCREEN", "150"))
KNN_K    <- as.integer(Sys.getenv("NAT_KNN", "5"))
COV_SRC  <- Sys.getenv("NAT_COV", "wdi")
ISO3 <- c(gambia = "GMB", ghana = "GHA", malawi = "MWI", sierraleone = "SLE")
TGC  <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")

nat <- vmnis_national()
cov <- load_national_covariates(COV_SRC)
V <- cov$vars; LOGV <- cov$log_vars
dat <- nat |> inner_join(cov$df |> select(iso3c, year, all_of(V)), by = c("iso3c", "year"))
message(sprintf("VMNIS joined to %s covariates: %d rows / %d countries", cov$source, nrow(dat), n_distinct(dat$iso3c)))

# ── 1. leave-one-country-out, every panel, SuperLearner beside the old arms ──
panels <- unique(do.call(rbind, lapply(NAT_PANEL_MAP, function(x) data.frame(mn = x[1], pop = x[2], stringsAsFactors = FALSE))))
loco <- list(); pred <- list()
for (i in seq_len(nrow(panels))) {
  lab <- paste(panels$mn[i], "|", panels$pop[i])
  d <- dat |> filter(mn_group == panels$mn[i], pop == panels$pop[i])
  r <- tryCatch(fit_national_loco(d, V, LOGV, lab, cov$source, SCREEN_K, KNN_K, sl = TRUE),
                error = function(e) { message("  ", lab, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  loco[[lab]] <- r$metrics; pred[[lab]] <- r$predictions
  message(sprintf("  %-42s %3d surveys / %2d countries", lab, r$metrics$n_surveys[1], r$metrics$n_countries[1]))
}
LT <- bind_rows(loco)
readr::write_csv(LT, here("results", "tables", "national_vmnis_loco_sl.csv"))
readr::write_csv(bind_rows(pred), here("results", "tables", "national_vmnis_loco_sl_pred.csv"))
print(as.data.frame(LT |> select(panel, model, n_countries, mae_pp, spearman) |> arrange(panel, mae_pp)), row.names = FALSE)

# ── 2. our four countries, vitamin A, held out of VMNIS ──
TG <- read.csv(here("results", "tables", "protocol_v2", "targets_v2.csv"), stringsAsFactors = FALSE)
rows <- list()
for (oc in c("child_vitA", "women_vitA")) {
  pan <- national_panel_for(oc); lab <- paste(pan[1], "|", pan[2])
  d <- dat |> filter(mn_group == pan[1], pop == pan[2])
  for (cn in names(ISO3)) {
    t <- TG[TG$country == TGC[[cn]] & TG$outcome == oc & is.finite(TG$y_prev) & is.finite(TG$n_raw), ]
    if (!nrow(t)) next
    yr <- SURVEY_YEAR[[cn]]
    lv <- sapply(c("ridge", "sl"), function(l) tryCatch(
      predict_national_level(d, V, LOGV, ISO3[[cn]], yr, cov$df, SCREEN_K, KNN_K, learner = l), error = function(e) NA_real_))
    tr_prev <- d$prev[is.finite(d$prev) & d$iso3c != ISO3[[cn]]]
    rows[[length(rows) + 1L]] <- data.frame(
      outcome = oc, country = cn, panel = lab, year = yr,
      survey_pp = round(100 * sum(t$n_raw * t$y_prev) / sum(t$n_raw), 2),
      ridge_pp = round(100 * lv[["ridge"]], 2), sl_pp = round(100 * lv[["sl"]], 2),
      null_pp = round(100 * stats::plogis(mean(stats::qlogis(pmin(pmax(tr_prev, .005), .995)))), 2),
      stringsAsFactors = FALSE)
    message(sprintf("  %-11s %-10s %d  survey %5.1f  ridge %5.1f  SL %5.1f", cn, oc, yr,
                    rows[[length(rows)]]$survey_pp, rows[[length(rows)]]$ridge_pp, rows[[length(rows)]]$sl_pp))
  }
}
LV <- bind_rows(rows)
readr::write_csv(LV, here("results", "tables", "national_levels_sl.csv"))
print(LV, row.names = FALSE)
