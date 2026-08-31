# =============================================================================
# scripts/covariates/19_national_composition.R
#
# LEVEL FROM VMNIS, PATTERN FROM THE SUBNATIONAL MODEL.
#
# The subnational work ends at an awkward place. Anchoring a district map to a
# region's design-based survey total more than doubles its rank correlation and
# cuts absolute bias from 3.60 to 1.49 pp -- but a country with no survey has no
# anchor to use, so its transported map is a ranking and not a set of
# prevalences. That is the gap this closes: predict the national LEVEL from the
# WHO VMNIS panel (~70 countries, national covariates), and shift the
# transported district PATTERN onto it.
#
# THE ARMS, all scored on the same districts against the same survey truth:
#   transported (no level)       the LOCO map as it comes out of the model
#   + VMNIS level                the composition -- level predicted with the
#                                target country held out of VMNIS entirely
#   + national null level        the level a VMNIS model WITHOUT covariates
#                                gives (the training-country mean). Separates
#                                "having any level at all" from "having a
#                                predicted one", which the composition needs or
#                                it takes credit for the wrong thing.
#   + true national (oracle)     the survey's own national total. Not
#                                deployable; it bounds what any level fix can
#                                buy.
#
# The shift is monotone on the logit scale, so the district RANKING is identical
# in all four arms. Spearman is therefore constant by construction and is
# reported only as a check. What moves is MAE and bias, which is the point.
# Correlation is undefined in the cells where the transported map is flat
# across a country's districts, so it is averaged over the cells that have one
# and the count is printed alongside.
#
# SCOPE, stated up front. VMNIS carries Folate, Vitamin A, Vitamin B12,
# Vitamin D and Zinc and has NO iron panel; the LOCO transport predictions cover
# iron and vitamin A. The intersection is vitamin A, so this runs on 2 outcomes
# x 4 countries. Vitamin A / preschool children is also VMNIS's strongest panel,
# so read the result as an upper bound, not an average over nutrients.
#
# WHAT IT FOUND, 2026-08-31: THE COMPOSITION DOES NOT WORK, AND CANNOT.
#
#   arm                        MAE pp   |bias| pp   better than no level
#   transported (no level)       5.81       3.35    --
#   + national null level        8.22       4.69    1/8
#   + VMNIS level               12.70      10.03    1/8
#   + true national (oracle)     5.58       0.73    4/8
#
# Two separate failures, and the second is the one that settles it.
#
# 1. The VMNIS LEVEL IS ON A DIFFERENT SCALE. The model predicts 40.9 pp for
#    Sierra Leone against a true 12.2, and 31.2 for Malawi against 9.2; the
#    VMNIS vitamin A panel averages 24.9 pp while our four surveys run 9-20.
#    The no-covariate null beats the covariate model in 6 of 8 cells, so this
#    is not a weak model -- it is a model predicting a different quantity. This
#    is the same cross-survey biomarker LEVEL offset that defeats LOCO
#    transport: VMNIS reports what national surveys report, with their own
#    assays and inflammation adjustments, not what ours would have measured.
#
# 2. EVEN A PERFECT NATIONAL LEVEL BUYS ALMOST NOTHING: 5.81 -> 5.58 pp, a
#    0.23 pp gain, and it wins only 4 of 8 cells. So the ceiling on this whole
#    approach is a quarter of a percentage point, and no amount of further
#    VMNIS modelling can raise it. The reason is the same one the Admin-1
#    anchoring result already showed: a national number moves every district by
#    the same amount, including the ones that were already right. What pays is
#    the ADMIN-1 anchor, which corrects region by region (r 0.170 -> 0.406) --
#    and that needs a survey, which is exactly what the no-survey country does
#    not have.
#
# The VMNIS model is kept and promoted anyway, because it is good at the thing
# it is actually good at: RANKING countries. Vitamin A / preschool children
# reaches r = 0.655 against a noise ceiling of 0.818, i.e. 80% of what is
# attainable. That supports "which countries are worst", not "what is the
# prevalence in this district".
#
#   Rscript scripts/covariates/19_national_composition.R
# -> results/tables/national_vmnis_loco.csv        the VMNIS model itself
# -> results/tables/national_vmnis_ceiling.csv     variance-component ceilings
# -> results/tables/national_composition.csv       the four arms
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(here); library(glmnet); library(ranger)
})
source(here("R", "national_covariates.R"))
source(here("R", "national_vmnis.R"))

SCREEN_K <- as.integer(Sys.getenv("NAT_SCREEN", "150"))
KNN_K    <- as.integer(Sys.getenv("NAT_KNN", "5"))
COV_SRC  <- Sys.getenv("NAT_COV", "wdi")

# Survey year per country, for pointing the national model at the right row.
SURVEY_YEAR <- c(gambia = 2021, ghana = 2017, malawi = 2015, sierraleone = 2013)
ISO3 <- c(gambia = "GMB", ghana = "GHA", malawi = "MWI", sierraleone = "SLE")

nat <- vmnis_national()
cov <- load_national_covariates(COV_SRC)
V <- cov$vars; LOGV <- cov$log_vars
message(sprintf("VMNIS: %d rows / %d countries | covariates: %s (%d)",
                nrow(nat), n_distinct(nat$iso3c), cov$source, length(V)))

dat <- nat |>
  inner_join(cov$df |> select(iso3c, year, all_of(V)), by = c("iso3c", "year"))
message(sprintf("joined to covariates: %d rows / %d countries",
                nrow(dat), n_distinct(dat$iso3c)))

# ---------------------------------------------------------------------------
# 1. The VMNIS model, scored leave-one-country-out
# ---------------------------------------------------------------------------
cat("\n===================== VMNIS national model (LOCO) =====================\n")
panels <- unique(do.call(rbind, lapply(NAT_PANEL_MAP, function(x)
  data.frame(mn = x[1], pop = x[2], stringsAsFactors = FALSE))))

loco <- list(); loco_pred <- list()
for (i in seq_len(nrow(panels))) {
  mn <- panels$mn[i]; pp <- panels$pop[i]
  lab <- paste(mn, "|", pp)
  d <- dat |> filter(mn_group == mn, pop == pp)
  r <- tryCatch(fit_national_loco(d, V, LOGV, lab, cov$source, SCREEN_K, KNN_K),
                error = function(e) { message("  ", lab, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) { message("  skipped (too few countries): ", lab); next }
  loco[[lab]] <- r$metrics; loco_pred[[lab]] <- r$predictions
  message(sprintf("  %-42s %3d surveys / %2d countries",
                  lab, r$metrics$n_surveys[1], r$metrics$n_countries[1]))
}
loco_tab <- bind_rows(loco)
readr::write_csv(loco_tab, here("results", "tables", "national_vmnis_loco.csv"))
print(as.data.frame(loco_tab |> arrange(panel, mae_pp)), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. The ceiling, per panel -- what fraction of it the model reaches
# ---------------------------------------------------------------------------
cat("\n======================= noise ceiling by panel ========================\n")
adj_level <- function(x) {
  x <- trimws(as.character(x))
  case_when(grepl("inflam", x, ignore.case = TRUE) ~ "inflammation_adjusted",
            grepl("^none$", x, ignore.case = TRUE) ~ "unadjusted",
            x == "" | grepl("not specified", x, ignore.case = TRUE) ~ "unspecified",
            TRUE ~ "other")
}
ceil <- list()
for (i in seq_len(nrow(panels))) {
  mn <- panels$mn[i]; pp <- panels$pop[i]; lab <- paste(mn, "|", pp)
  d <- nat |> filter(mn_group == mn, pop == pp) |>
    mutate(method = paste(adj_level(Dataadjustedfor),
                          trimws(as.character(Indicator))))
  cc <- tryCatch(national_noise_ceiling(d), error = function(e) NULL)
  if (is.null(cc)) next
  ceil[[lab]] <- cbind(panel = lab, cc)
}
ceil_tab <- bind_rows(ceil)
if (nrow(ceil_tab)) {
  best <- loco_tab |> filter(model != "null") |> group_by(panel) |>
    summarise(best_r = max(pearson, na.rm = TRUE), .groups = "drop")
  ceil_tab <- ceil_tab |> left_join(best, by = "panel") |>
    mutate(r_share = round(best_r / r_max_report, 2))
  readr::write_csv(ceil_tab, here("results", "tables", "national_vmnis_ceiling.csv"))
  print(as.data.frame(ceil_tab |> select(panel, surveys, countries,
    sd_country, sd_method, sd_resid, sd_sampling,
    r_max_report, r_max_standardised, best_r, r_share)), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 3. The composition
# ---------------------------------------------------------------------------
cat("\n============ composition: level from VMNIS, pattern from LOCO =========\n")
tp <- read.csv(here("results", "tables",
                    "transportability_area_loco_predictions.csv"),
               stringsAsFactors = FALSE)

rows <- list(); lvl_rows <- list()
for (oc in intersect(unique(tp$outcome), names(NAT_PANEL_MAP))) {
  pan <- national_panel_for(oc); lab <- paste(pan[1], "|", pan[2])
  d <- dat |> filter(mn_group == pan[1], pop == pan[2])
  if (!nrow(d)) { message("  no VMNIS panel for ", oc); next }

  for (cn in unique(tp$country[tp$outcome == oc])) {
    iso <- ISO3[[cn]]; yr <- SURVEY_YEAR[[cn]]
    if (is.null(iso)) next
    pat <- tp[tp$outcome == oc & tp$country == cn, , drop = FALSE]
    pat <- pat[is.finite(pat$modeled_prev) & is.finite(pat$survey_prev), , drop = FALSE]
    if (nrow(pat) < 8) next

    w <- pat$n_svy; w[!is.finite(w)] <- 0
    if (sum(w) <= 0) w <- rep(1, nrow(pat))
    true_nat <- sum(w * pat$survey_prev) / sum(w)

    # The predicted level, with this country held out of VMNIS entirely.
    lvl <- tryCatch(predict_national_level(d, V, LOGV, iso, yr, cov$df,
                                           SCREEN_K, KNN_K),
                    error = function(e) NA_real_)
    # The no-covariate comparator: what a level costs when it carries no
    # information beyond the training countries' average.
    tr_prev <- d$prev[is.finite(d$prev) & d$iso3c != iso]
    lvl_null <- if (length(tr_prev) >= 10)
      stats::plogis(mean(stats::qlogis(pmin(pmax(tr_prev, .005), .995)))) else NA_real_

    lvl_rows[[length(lvl_rows) + 1L]] <- data.frame(
      outcome = oc, country = cn, panel = lab, year = yr,
      true_national_pp = round(100 * true_nat, 2),
      vmnis_level_pp   = round(100 * lvl, 2),
      null_level_pp    = round(100 * lvl_null, 2),
      vmnis_err_pp     = round(100 * abs(lvl - true_nat), 2),
      null_err_pp      = round(100 * abs(lvl_null - true_nat), 2),
      stringsAsFactors = FALSE)

    arms <- list(
      `transported (no level)`  = pat$modeled_prev,
      `+ national null level`   = shift_to_level(pat$modeled_prev, lvl_null, w),
      `+ VMNIS level`           = shift_to_level(pat$modeled_prev, lvl, w),
      `+ true national (oracle)`= shift_to_level(pat$modeled_prev, true_nat, w))

    for (a in names(arms)) {
      p <- arms[[a]]; ok <- is.finite(p)
      if (sum(ok) < 8) next
      rows[[length(rows) + 1L]] <- data.frame(
        outcome = oc, country = cn, panel = lab, arm = a,
        n_admin2 = sum(ok),
        mae_pp  = round(100 * mean(abs(p[ok] - pat$survey_prev[ok])), 2),
        bias_pp = round(100 * mean(p[ok] - pat$survey_prev[ok]), 2),
        pearson  = round(suppressWarnings(stats::cor(p[ok], pat$survey_prev[ok])), 3),
        spearman = round(suppressWarnings(stats::cor(
          p[ok], pat$survey_prev[ok], method = "spearman")), 3),
        stringsAsFactors = FALSE)
    }
    message(sprintf("  %-11s %-11s true=%5.1f pp  vmnis=%5.1f pp  null=%5.1f pp",
                    cn, oc, 100 * true_nat, 100 * lvl, 100 * lvl_null))
  }
}

comp <- bind_rows(rows)
lvls <- bind_rows(lvl_rows)
readr::write_csv(comp, here("results", "tables", "national_composition.csv"))
readr::write_csv(lvls, here("results", "tables", "national_composition_levels.csv"))

cat("\n--- how well is the national LEVEL itself predicted? (pp error) ---\n")
print(as.data.frame(lvls), row.names = FALSE)

cat("\n--- district accuracy by arm (mean over cells) ---\n")
ARM_ORDER <- c("transported (no level)", "+ national null level",
               "+ VMNIS level", "+ true national (oracle)")
print(as.data.frame(comp |>
  mutate(arm = factor(arm, levels = ARM_ORDER)) |> arrange(arm) |>
  group_by(arm) |>
  summarise(cells = n(), mae_pp = round(mean(mae_pp), 2),
            absbias_pp = round(mean(abs(bias_pp)), 2),
            # Correlation is undefined where the transported map is flat across
            # a country's districts (Gambia women_vitA, Sierra Leone after a
            # large shift saturates the logit). Averaging over the cells that
            # do have one, and saying how many, beats printing NA.
            r_cells = sum(is.finite(pearson)),
            pearson = round(mean(pearson, na.rm = TRUE), 3),
            spearman = round(mean(spearman, na.rm = TRUE), 3), .groups = "drop")),
  row.names = FALSE)

cat("\n--- per cell, MAE by arm ---\n")
print(as.data.frame(comp |> select(country, outcome, arm, mae_pp) |>
  tidyr::pivot_wider(names_from = arm, values_from = mae_pp)), row.names = FALSE)

base <- comp |> filter(arm == "transported (no level)") |>
  select(country, outcome, base_mae = mae_pp)
win <- comp |> filter(arm == "+ VMNIS level") |>
  left_join(base, by = c("country", "outcome")) |>
  mutate(better = mae_pp < base_mae)
cat(sprintf("\nVMNIS level beats no level in %d of %d cells.\n",
            sum(win$better, na.rm = TRUE), nrow(win)))
cat("\n-> results/tables/national_composition.csv\n")
