# =============================================================================
# global.R — loaded once at app startup, shared across all sessions
# =============================================================================
# REBUILT 2026-09-13 on the corrected protocol (protocol v2, RR-10 tables).
# Every number the app shows now comes from dashboard/data-raw/05_build_
# protocol_v2_bundles.R, which reads the committed protocol result tables and
# fits one thing: the deployment ranking (the zero-tuning index fitted on all
# surveyed districts of a country and applied to every district). The
# person-level SuperLearner, the area-level recipe, the Fay-Herriot and BYM2
# layers, the old leaderboard, the corrected-methods P1-P8 comparison, the GBD
# placeholder and the Sierra Leone chiefdom layer are gone: each rested on an
# evaluation the audit withdrew, or on a model the protocol found no better
# than chance. The decks and the manuscript say the same things this app does.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(leaflet)
  library(plotly)
  library(reactable)
  library(htmltools)
})

DATA_DIR <- "data"
BRIEF_DIR <- "briefs"
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
.rds <- function(f) { p <- file.path(DATA_DIR, f); if (file.exists(p)) readRDS(p) else NULL }

# ── Data ───────────────────────────────────────────────────────────────────
IDX <- .rds("admin2_index.rds")
if (is.null(IDX)) stop("dashboard/data/admin2_index.rds is missing: run dashboard/data-raw/05_build_protocol_v2_bundles.R")
idx_districts <- IDX$districts     # one row per district x outcome, every district
idx_national  <- IDX$national      # the survey's national prevalence per cell (the anchor)
idx_fits      <- IDX$fits          # per cell: exact per-predictor weights (logit scale)
idx_xr        <- IDX$xr            # per country: the rank-normal predictor matrix the model sees

admin2_pop  <- readRDS(file.path(DATA_DIR, "admin2_population.rds"))
admin2_bnds <- readRDS(file.path(DATA_DIR, "admin2_boundaries.rds"))
admin1_bnds <- readRDS(file.path(DATA_DIR, "admin1_boundaries.rds"))
meta        <- readRDS(file.path(DATA_DIR, "metadata.rds"))

CIV <- .rds("civ_index.rds")                 # Cote d'Ivoire from climate and soil
EV  <- .rds("protocol_evidence.rds") %||% list()   # the tables behind the trust tabs
CAT <- .rds("predictor_catalogue.rds")       # every predictor, described
UE  <- .rds("uncertainty_ensembles.rds") %||% list()   # stability ensembles (UE-01): rank ranges, WHO exceedance, weight ranges

data_build_time <- format(IDX$build_time, "%Y-%m-%d %H:%M")
PROTOCOL_LABEL  <- IDX$protocol %||% "protocol v2"

# ── Source helpers and modules ────────────────────────────────────────────
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f, local = FALSE)

# ── Lookups ────────────────────────────────────────────────────────────────
country_choices  <- setNames(names(meta$countries), meta$countries)
outcomes_present <- intersect(names(meta$outcome_labels), unique(idx_districts$outcome))
outcome_choices  <- setNames(outcomes_present, meta$outcome_labels[outcomes_present])
outcomes_for <- function(ck) {
  ocs <- intersect(names(meta$outcome_labels), unique(idx_districts$outcome[idx_districts$country_key == ck]))
  setNames(ocs, meta$outcome_labels[ocs])
}
# short outcome names for axes and tables
outcome_short <- c(child_vitA = "Vitamin A, children", women_vitA = "Vitamin A, women",
                   child_iron = "Iron, children", women_iron = "Iron, women",
                   women_folate = "Folate, women", women_b12 = "B12, women",
                   child_zinc = "Zinc, children", women_zinc = "Zinc, women",
                   child_selenium = "Selenium, children", women_selenium = "Selenium, women",
                   women_iodine = "Iodine, women")
# Method names as a first-time reader sees them. One name per method,
# everywhere in the app.
arm_label <- c(domain_index = "This model (public data)", spatial_plus_domain = "Neighbouring districts + public data",
               spatial = "Neighbouring districts' average", domain_enet = "Regression on the same data groups",
               raw_enet = "Regression on all data layers", region_mean_jk = "Survey's regional averages",
               null_train_mean = "No information (same value everywhere)")
# Domain names as displayed (18 Sep terminology; the data-side names are fixed
# by the DP-01 prefix rule, so the mapping is display-only)
domain_display <- c("Infection and inflammation burden" = "Infectious disease burden",
                    "Anaemia and haemoglobin" = "Anaemia (modelled)",
                    "Dietary inadequacy (MODELLED SURFACE" = "Dietary inadequacy (modelled)",
                    "Nutrition status (MODELLED SURFACE)" = "Nutrition status (modelled)")
dom_disp <- function(x) { y <- unname(domain_display[x]); ifelse(is.na(y), x, y) }

# ── Headline numbers, read from the evidence tables ───────────────────────
# The same quantities the decks and the manuscript quote, computed here so no
# tab types a number by hand.
g1 <- function(x) if (length(x) && is.finite(x[1])) x[1] else NA_real_
.tbl <- function(name) EV[[name]]
bm <- function(est, tg, arm, col = "mean_spearman") {
  B <- .tbl("benchmarks_summary"); if (is.null(B)) return(NA_real_)
  g1(B[[col]][B$estimand == est & B$target == tg & B$arm == arm])
}
Q <- list()
Q$infill      <- bm("infill", "level", "domain_index");   Q$infill_prev <- bm("infill", "prev", "domain_index")
Q$infill_jk   <- bm("infill", "level", "region_mean_jk"); Q$infill_jk_prev <- bm("infill", "prev", "region_mean_jk")
Q$infill_sp   <- bm("infill", "level", "spatial")
Q$region      <- bm("region", "level", "domain_index")
Q$tr          <- bm("country", "level", "domain_index");  Q$tr_prev <- bm("country", "prev", "domain_index")
Q$tr_pos      <- bm("country", "level", "domain_index", "cells_positive")
Q$tr_n        <- bm("country", "level", "domain_index", "cells")
Q$topk        <- bm("infill", "level", "domain_index", "mean_topk")
Q$topk_tr     <- bm("country", "level", "domain_index", "mean_topk")
local({
  ND <- .tbl("nested_domains"); CS1 <- .tbl("climate_soil_admin1"); A1 <- .tbl("admin1_transport"); NC <- .tbl("transport_null")
  Q$cs     <<- if (is.null(ND)) NA else g1(mean(ND$spearman[ND$arm == "fixed_cs" & ND$target == "level"], na.rm = TRUE))
  Q$cs_pos <<- if (is.null(ND)) NA else sum(ND$spearman[ND$arm == "fixed_cs" & ND$target == "level"] > 0, na.rm = TRUE)
  Q$cs_a1  <<- if (is.null(CS1)) NA else g1(mean(CS1$spearman[CS1$set == "climate_soil" & CS1$target == "level"], na.rm = TRUE))
  Q$tr_a1  <<- if (is.null(A1)) NA else g1(mean(A1$spearman[A1$arm == "domain_index" & A1$target == "level"], na.rm = TRUE))
  Q$null_d <<- if (is.null(NC)) NA else g1(NC$null_mean_q95[NC$tier == "admin2"])
  Q$null_a1 <<- if (is.null(NC)) NA else g1(NC$null_mean_q95[NC$tier == "admin1"])
  NT <- .tbl("targeting_summary")
  cap <- function(a) if (is.null(NT)) NA else g1(NT$mean_capture[NT$estimand == "infill" & NT$arm == a])
  Q$cap_index <<- cap("domain_index"); Q$cap_jk <<- cap("region_mean_jk")
  Q$cap_null  <<- cap("null_train_mean"); Q$cap_oracle <<- cap("oracle_ceiling")
  RC <- .tbl("risk_summary")
  Q$band_exact <<- if (is.null(RC)) NA else g1(RC$exact_admin2[RC$scheme == "who_vitA" & RC$arm == "domain_index" & RC$estimand == "infill"])
  Q$band_w1    <<- if (is.null(RC)) NA else g1(RC$within1_admin2[RC$scheme == "who_vitA" & RC$arm == "domain_index" & RC$estimand == "infill"])
  AR <- .tbl("design_summary")
  ar <- function(design) if (is.null(AR)) NA else g1(AR$mae[AR$fraction == 0.05 & AR$rank_from == "prev" & AR$set == "climate_soil" & AR$design == design])
  Q$ar_a1 <<- ar("A1_anchor_rank"); Q$ar_c <<- ar("C_regional_survey"); Q$ar_b <<- ar("B_district_survey")
  Q$ar_b_match <<- if (is.null(AR) || !is.finite(Q$ar_a1)) NA else {
    b <- AR[AR$rank_from == "prev" & AR$set == "climate_soil" & AR$design == "B_district_survey", ]
    g1(min(b$fraction[b$mae <= Q$ar_a1])) }
  VC <- .tbl("ceiling")
  Q$ceiling_prev  <<- if (is.null(VC)) NA else g1(mean(VC$ceiling_vc[VC$rung == "admin2" & VC$target == "prev"], na.rm = TRUE))
  Q$ceiling_level <<- if (is.null(VC)) NA else g1(mean(VC$ceiling_vc[VC$rung == "admin2" & VC$target == "level"], na.rm = TRUE))
  Q$vc_at <<- if (is.null(VC)) NA else { v <- VC[VC$rung == "admin2" & VC$target == "prev", ]; sum(v$achieved_spearman >= v$ceiling_vc, na.rm = TRUE) }
  Q$vc_n  <<- if (is.null(VC)) NA else sum(is.finite(VC$ceiling_vc[VC$rung == "admin2" & VC$target == "prev"]))
  CAL <- .tbl("worst_fifth_calibration")
  Q$cal_top <<- if (is.null(CAL)) NA else g1(CAL$share_in_survey_worst_fifth[CAL$band == "80 to 100%"])
  Q$cal_low <<- if (is.null(CAL)) NA else g1(CAL$share_in_survey_worst_fifth[CAL$band == "0 to 20%"])
  TC <- .tbl("training_curve")
  if (!is.null(TC)) { tc <- TC[TC$arm == "domain_index" & TC$target == "level", ]
    lc <- tapply(tc$spearman, tc$n_train_countries, mean, na.rm = TRUE)
    Q$lc1 <<- g1(lc[["1"]]); Q$lc_max <<- g1(lc[[as.character(max(tc$n_train_countries))]])
    Q$lc_step <<- (Q$lc_max - Q$lc1) / (max(tc$n_train_countries) - 1) }
  ID <- .tbl("importance_domains")
  if (!is.null(ID)) { e <- ID[ID$scope == "pooled" & ID$target == "level" & ID$domain %in% c("Satellite embedding", "Climate and weather", "Soil characteristics"), ]
    s <- tapply(e$share, e$outcome, sum); Q$env_lo <<- g1(min(s)); Q$env_hi <<- g1(max(s)) }
  WS <- .tbl("weight_sources")
  Q$sparse20_tr <<- if (is.null(WS)) NA else g1(WS$mean_spearman[WS$estimand == "country" & WS$arm == "sparse20" & WS$target == "level"])
  MB <- .tbl("geostat_cells")
  if (!is.null(MB)) { m <- MB[MB$estimand == "infill", ]
    gm <- function(tg, a, col = "rho_agg") g1(mean(m[[col]][m$target == tg & m$arm == a], na.rm = TRUE))
    Q$mbg_rank <<- gm("level", "mbg"); Q$mbg_rank_index <<- gm("level", "domain_index")
    Q$mbg_err <<- gm("prev", "mbg", "wmae_agg"); Q$mbg_err_index <<- gm("prev", "domain_index", "wmae_agg") }
  IL <- .tbl("individual_level")
  Q$il_auc <<- if (is.null(IL)) NA else g1(mean(IL$auc, na.rm = TRUE))
})
# 2026-09-27: quantities for the external-validation, roadmap and planning panels
local({
  XP <- .tbl("xv_pooled")
  if (!is.null(XP)) {
    XP$set <- trimws(XP$set)
    gx <- function(s, col) g1(XP[[col]][XP$set == s])
    Q$xv_level       <<- gx("Africa / iSDA / level", "mean_rho")
    Q$xv_level_cells <<- gx("Africa / iSDA / level", "cells")
    Q$xv_level_pos   <<- gx("Africa / iSDA / level", "positive")
    Q$xv_level_null  <<- gx("Africa / iSDA / level", "block_null_p95")
    Q$xv_prev        <<- gx("Africa / iSDA / prev", "mean_rho")
    Q$xv_sg_level    <<- gx("Africa / SoilGrids / level", "mean_rho")
    Q$xv_off         <<- gx("Off-continent / SoilGrids / prev", "mean_rho")
    Q$xv_off_pos     <<- gx("Off-continent / SoilGrids / prev", "positive")
    Q$xv_off_cells   <<- gx("Off-continent / SoilGrids / prev", "cells")
    Q$xv_off_null    <<- gx("Off-continent / SoilGrids / prev", "block_null_p95")
  }
  XC <- .tbl("xv_cells")
  Q$xv_countries <<- if (!is.null(XC)) length(unique(XC$country)) else NA
  BC <- .tbl("benchmarks_cells")
  if (!is.null(BC)) {
    s <- BC$spearman[BC$estimand == "infill" & BC$target == "level" & BC$arm == "domain_index"]
    s <- s[is.finite(s)]
    Q$strong_n <<- sum(s >= 0.5); Q$strong_mean <<- mean(s[s >= 0.5]); Q$infill_cells <<- length(s)
  }
  HR <- .tbl("headroom")
  if (!is.null(HR)) { h <- HR[is.finite(HR$r_max_emp), ]; Q$rmax_med <<- stats::median(h$r_max_emp) }
  XO <- .tbl("cross_outcome")
  if (!is.null(XO)) {
    xo <- XO[XO$domain_set == "cs" & XO$target == "level" & XO$model == "domain_index" & XO$outcome_class == "iron / vitA", ]
    Q$xo_base <<- g1(xo$base_mean[xo$arm == "block_only_same_nutrient"])
    Q$xo_same <<- g1(xo$arm_mean[xo$arm == "block_only_same_nutrient"])
    Q$xo_added <<- g1(xo$mean_delta[xo$arm == "same_nutrient"])
  }
  RC2 <- .tbl("rank_coverage")
  Q$stab_cov <<- if (!is.null(RC2)) mean(RC2$coverage_90, na.rm = TRUE) else NA
  PS <- .tbl("planner_summary")
  if (!is.null(PS)) {
    p5 <- PS[PS$fraction == 0.5, ]
    Q$plan_random <<- g1(p5$spearman[p5$arm == "random"])
    Q$plan_spread <<- g1(p5$spearman[p5$arm == "spread_model"])
    Q$plan_delta  <<- g1(p5$mean_delta[p5$arm == "spread_model"])
    Q$plan_better <<- g1(p5$cells_better[p5$arm == "spread_model"])
    Q$plan_cells  <<- g1(p5$cells[p5$arm == "spread_model"])
  }
  SEp <- .tbl("selenium_protocol")
  if (!is.null(SEp)) {
    pr <- SEp[SEp$analysis == "protocol" & SEp$arm == "domain_index", ]
    Q$se_mean <<- if (nrow(pr)) mean(pr$spearman, na.rm = TRUE) else NA
  }
  CPc <- .tbl("conformal_prev")
  if (!is.null(CPc)) {
    Q$cal_cov      <<- mean(CPc$loo_coverage_90, na.rm = TRUE)
    Q$cal_half_med <<- stats::median(CPc$half_width_pp, na.rm = TRUE)
    Q$stabprev_cov <<- mean(CPc$stability_band_coverage, na.rm = TRUE)
  }
  # quantities the reader-facing text quotes (text review, 2026-09-27)
  if (!is.null(BC)) {
    s <- BC$spearman[BC$estimand == "infill" & BC$target == "level" & BC$arm == "domain_index"]
    s <- s[is.finite(s)]; Q$weak_mean <<- mean(s[s < 0.5]); Q$weak_n <<- sum(s < 0.5)
  }
  if (!is.null(SEp)) {
    sl <- SEp[SEp$analysis == "protocol" & SEp$estimand == "infill" & SEp$target == "level" &
                SEp$outcome %in% c("child_selenium", "women_selenium"), ]
    Q$se_idx_lo <<- min(sl$spearman[sl$arm == "domain_index"]); Q$se_idx_hi <<- max(sl$spearman[sl$arm == "domain_index"])
    Q$se_jk_lo  <<- min(sl$spearman[sl$arm == "region_mean_jk"]); Q$se_jk_hi <<- max(sl$spearman[sl$arm == "region_mean_jk"])
  }
  RCs <- .tbl("risk_summary")
  if (!is.null(RCs)) {
    r <- RCs[RCs$scheme == "who_vitA" & RCs$arm == "region_mean_jk" & RCs$estimand == "infill", ]
    Q$band_exact_jk <<- g1(r$exact_admin2); Q$band_w1_jk <<- g1(r$within1_admin2)
  }
  PVt <- .tbl("planner_validation")
  Q$plan_transport <<- if (!is.null(PVt)) mean(PVt$spearman_transport_only, na.rm = TRUE) else NA
  Q$plan_pps <<- if (!is.null(PS)) g1(PS$spearman[PS$fraction == 0.5 & PS$arm == "pps"]) else NA
  # Start here, the three questions (national level, district prevalence, ranking)
  if (!is.null(PVt)) {
    at <- function(a, col) mean(abs(PVt[[col]][PVt$arm == a & PVt$fraction == 0.5]), na.rm = TRUE)
    Q$nat_bias_random <<- at("random", "nat_bias_pp");  Q$nat_ci_random <<- at("random", "nat_ci_pp")
    Q$nat_bias_strat  <<- at("spread_model", "nat_strat_bias_pp"); Q$nat_ci_strat <<- at("spread_model", "nat_strat_ci_pp")
    Q$nat_bias_ext    <<- at("extremes_model", "nat_bias_pp")
  }
  NV <- .tbl("national_vmnis")
  if (!is.null(NV)) {
    NV$model[is.na(NV$model) | NV$model == ""] <- "null"
    best <- tapply(NV$mae_pp[NV$model != "null"], NV$panel[NV$model != "null"], min)
    Q$natpred_lo <<- min(best); Q$natpred_hi <<- max(best)
    Q$natpred_ctry_lo <<- min(NV$n_countries); Q$natpred_ctry_hi <<- max(NV$n_countries)
    nul <- tapply(NV$mae_pp[NV$model == "null"], NV$panel[NV$model == "null"], min)[names(best)]
    Q$natpred_panels <<- length(best); Q$natpred_nobetter <<- sum(best > nul - 0.5, na.rm = TRUE)
  }
  NL <- .tbl("national_levels")
  Q$natpred_own_hi <<- if (!is.null(NL)) max(NL$vmnis_err_pp, na.rm = TRUE) else NA
  if (!is.null(CPc)) Q$prev_err <<- mean(CPc$median_abs_err_pp, na.rm = TRUE)
  SMq <- .tbl("survey_design_meta")
  if (!is.null(SMq)) { Q$anchor_n_lo <<- 0.05 * min(SMq$n_raw, na.rm = TRUE); Q$anchor_n_hi <<- 0.05 * max(SMq$n_raw, na.rm = TRUE) }
})
Q$n_predictors <- if (!is.null(CAT)) nrow(CAT$variables) else 570
Q$n_domains    <- if (!is.null(CAT)) nrow(CAT$domains) else 28
Q$n_in_model   <- if (!is.null(CAT) && "in_model" %in% names(CAT$variables)) sum(CAT$variables$in_model, na.rm = TRUE) else NA
Q$n_dhs        <- if (!is.null(CAT) && "tier" %in% names(CAT$variables)) sum(grepl("^DHS", CAT$variables$tier)) else NA
Q$n_districts  <- length(unique(paste(idx_districts$country, idx_districts$Admin1, idx_districts$Admin2)))
Q$n_surveyed   <- sum(idx_national$n_surveyed[!duplicated(idx_national$country)])
f2 <- function(x, d = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")
pc <- function(x, d = 0) ifelse(is.finite(x), paste0(formatC(100 * x, format = "f", digits = d), "%"), "NA")

# ── Colours ────────────────────────────────────────────────────────────────
who_colors <- c("Low" = "#2c7bb6", "Mild" = "#abd9e9", "Moderate" = "#fdae61",
                "Severe" = "#d7191c", "No data" = "#cccccc")
PROXY_COL <- "#0F7B8A"; SURVEY_COL <- "#C8641E"; GREY_COL <- "#8c8c8c"

# ── Caveats ────────────────────────────────────────────────────────────────
selenium_caveat <- paste("Selenium was measured in Malawi only, so this ranking cannot be checked in other countries.",
                         "With each district hidden in turn, the model ranks selenium better than the survey's own",
                         "regional averages do; selenium in food follows soil geology, which the data layers capture.",
                         "WHO has no severity bands for selenium, so no severity class is shown.")
zinc_caveat <- paste("Zinc was measured in Malawi only, so this ranking cannot be checked in other countries.",
                     "Most of the difference between districts reflects the time of day blood was drawn rather than",
                     "place, and the model ranks zinc no better than chance.")
biomarker_caveats <- list(
  women_vitA   = paste("Vitamin A deficiency in women is below 3% in every survey, so there is little difference",
                       "between districts to rank. The marker used (retinol-binding protein) is also affected by",
                       "inflammation. Treat this ranking as weak."),
  child_vitA   = paste("Vitamin A is measured with retinol-binding protein, adjusted for inflammation and converted",
                       "to a retinol value with each survey's own conversion, then compared with the 0.70 cut-off."),
  women_b12    = paste("B12 was measured in three of the four countries (not The Gambia) with serum B12, an imperfect",
                       "marker. Treat district differences as indicative."),
  women_folate = paste("Folate was measured in three of the four countries. National prevalence ranges from 19% to",
                       "79% between them, so the percentages are less certain than the ranking."),
  child_zinc   = zinc_caveat,
  women_zinc   = zinc_caveat,
  child_selenium = selenium_caveat,
  women_selenium = selenium_caveat,
  women_iodine = paste("Iodine (urinary iodine below 100 micrograms per litre) was measured in Malawi only, so this",
                       "ranking cannot be checked in other countries. Salt iodisation can change iodine status quickly,",
                       "so treat this as a single-survey result.")
)
GENERAL_CAVEAT <- paste(
  "The four surveys were carried out between 2013 and 2018 by different teams, so prevalence levels are not",
  "directly comparable between countries. The ranking of districts carries over to a new country; prevalence",
  "levels do not, and need at least a national survey figure. Rankings for a country without a survey are less",
  "accurate than for a surveyed country (see How well it works).")

# ── Site banner ────────────────────────────────────────────────────────────
# The headline is the talks' thesis line (18 Sep terminology decisions).
SITE_SCOPE_HEADLINE <- "Modelling extends a survey's reach. It does not replace the survey."
SITE_SCOPE_BODY <- "The model ranks districts; percentages are planning estimates based on the national survey figure."
SITE_SCOPE_POINTER <- ""
site_banner <- div(
  class = "alert alert-info",
  style = "margin:0 0 10px;border-radius:0;text-align:center;font-size:0.88em;padding:6px 12px;",
  bsicons::bs_icon("info-circle"), " ",
  strong(SITE_SCOPE_HEADLINE), " ", SITE_SCOPE_BODY,
  tags$span(style = "color:#4a6b7c;", " ", SITE_SCOPE_POINTER)
)

# ── About and glossary ─────────────────────────────────────────────────────
about_content <- div(
  h5("About this dashboard", style = "margin-top: 0;"),
  p("This dashboard ranks the districts of The Gambia, Ghana, Sierra Leone and Malawi by how likely they are to",
    " have high levels of micronutrient deficiency. It uses public data and is checked against each country's",
    " national biomarker survey. It also ranks the districts of Cote d'Ivoire, which has not had such a survey.",
    " It is meant for ministries, funders and researchers deciding where to focus and where a future survey",
    " should collect samples."),
  h6("Method"),
  p(sprintf(paste("Each district is described by %d public data layers in %d groups, including satellite imagery,",
                  "climate, soil, crops, livestock, malaria and other diseases, food prices, and summaries of public",
                  "household surveys. The model uses %s of these layers. It leaves out summaries of the DHS surveys, so",
                  "that it relies only on data a country without a recent DHS survey could also collect. The layers in",
                  "each group are condensed into a few scores, each score is weighted by how closely it followed",
                  "deficiency in the surveyed districts, and the weighted scores are added up. No setting was adjusted",
                  "to improve the results. Every accuracy figure comes from districts, regions or whole countries that",
                  "were hidden from the model while it was built."),
            Q$n_predictors, Q$n_domains, if (is.finite(Q$n_in_model)) Q$n_in_model else "383")),
  h6("Accuracy in brief"),
  tags$ul(
    tags$li(sprintf("Inside a surveyed country, the model's ranking matches the survey's at %s, compared with %s for the survey's own regional averages and at most %s for a random ranking.",
                    f2(Q$infill), f2(Q$infill_jk), f2(Q$null_d))),
    tags$li(sprintf("In a country left out of the model, it scores %s with all data layers and %s with climate and soil layers only. The climate-and-soil version was chosen after seeing these results, so it is now being tested on new surveys.",
                    f2(Q$tr), f2(Q$cs))),
    tags$li(sprintf("Most districts have only one or two survey clusters, so even a perfect model would reach only about %s against the survey's district figures. This model reaches about two thirds of that.",
                    f2(Q$ceiling_level))),
    tags$li(sprintf("Compared with survey results held by the WHO for six more countries, the ranking scored %s in four African countries and %s in Pakistan and India (regional results).",
                    f2(Q$xv_level), f2(Q$xv_off)))
  ),
  h6("Citation"),
  p(em("Mertens et al. (in preparation). What geospatial covariates can and cannot do for sub-national",
       " micronutrient estimation: a protocol-first re-analysis of four national biomarker surveys.")),
  h6("Surveys"),
  p("The Gambia 2018, Ghana 2017, Sierra Leone 2013 and Malawi 2015 to 2016. Vitamin A and iron in children",
    " and women; folate and B12 in women in three countries; zinc, selenium and iodine in Malawi only."),
  p(em(sprintf("Data built %s. Analysis version for the project team: %s.", data_build_time, PROTOCOL_LABEL)), style = "color: #888; font-size: 0.85em;")
)

glossary_content <- div(
  h5("Glossary", style = "margin-top: 0;"),
  tags$dl(
    tags$dt("Priority score"),
    tags$dd("A district's place in its country's ranking, from 100 (ranked worst) down to near 0 (ranked best).",
            " It shows which districts are likely to be worse off, not how severe the problem is."),
    tags$dt("Ranking accuracy"),
    tags$dd("How closely the model's order of districts matches the survey's order, from 0 (no better than a random",
            " order) to 1 (the same order). It is a Spearman correlation. Differences smaller than 0.03 are ties."),
    tags$dt("Random ranking"),
    tags$dd(sprintf("The score a random order of districts reaches in 95%% of tries: %s for districts and %s for regions.",
                    f2(Q$null_d), f2(Q$null_a1))),
    tags$dt("Best achievable score"),
    tags$dd("The highest ranking accuracy any model could reach against the survey's district figures. It is below 1",
            " because most districts have only one or two survey clusters, so the survey figures are themselves uncertain."),
    tags$dt("Estimated prevalence"),
    tags$dd("The model's ranking converted into a percentage using the country's national survey figure. The order of",
            " districts comes from the model and the overall level from the survey. Where the model's district",
            " percentages are not informative, every district is shown close to the national figure."),
    tags$dt("Reliability of the district percentages"),
    tags$dd("How well the model's district percentages matched the survey's district figures in districts the model had",
            " not seen, as a correlation from 0 to 1: below 0.10 not informative, 0.10 to 0.30 weak, 0.30 to 0.50",
            " moderate, above 0.50 good. It affects the percentages, not the ranking."),
    tags$dt("Checked 90% range"),
    tags$dd(sprintf(paste("A range around each estimated prevalence, built from how far the model's estimates missed the",
                          "survey figures in districts the model had not seen. In that check, 90%% ranges contained the",
                          "survey figure %s of the time. Districts without survey data get the same width as the",
                          "surveyed districts in their country."), pc(Q$cal_cov))),
    tags$dt("Rank range when re-estimated"),
    tags$dd(sprintf(paste("How far a district's rank moves when the model is re-estimated on different samples of the",
                          "surveyed districts. It shows how firmly the model places a district. It is not a confidence",
                          "interval: in countries left out of the model, the survey's rank fell inside this range %s of the time."),
                    pc(Q$stab_cov))),
    tags$dt("Placed in the worst fifth (share of runs)"),
    tags$dd(sprintf(paste("How often a surveyed district was placed in its country's worst fifth when the model was",
                          "re-estimated 40 times with that district hidden. It points in the right direction but overstates",
                          "certainty: districts placed there in at least 80%% of runs were in the survey's worst fifth %s of",
                          "the time, against 20%% by chance."), pc(Q$cal_top))),
    tags$dt("The three tests"),
    tags$dd("District hidden: one district is hidden and predicted from the rest of its country. Region hidden: a whole",
            " region is hidden. Country left out: the model is built without that country, which is the situation of a",
            " country with no survey."),
    tags$dt("Percentage points"),
    tags$dd("The difference between two percentages. From 20% to 23% is 3 percentage points."),
    tags$dt("District and region"),
    tags$dd("District is the second administrative level and region the first. Rankings are made for districts.")
  )
)


# ── Concise version: shorter caveats, About and Glossary ──────────────────
# The full caveats stay available to the Technical notes page.
biomarker_caveats_long <- biomarker_caveats
biomarker_caveats <- list(
  women_vitA = "Below 3% in every survey, so there is little to rank. Treat as weak.",
  child_vitA = "Measured with retinol-binding protein, adjusted for inflammation.",
  women_b12 = "Measured in three countries (not The Gambia). Treat as indicative.",
  women_folate = "Measured in three countries. The percentages are less certain than the ranking.",
  child_zinc = "Malawi only. The model ranks zinc no better than chance.",
  women_zinc = "Malawi only. The model ranks zinc no better than chance.",
  child_selenium = "Malawi only, so it cannot be checked in another country. No WHO bands exist.",
  women_selenium = "Malawi only, so it cannot be checked in another country. No WHO bands exist.",
  women_iodine = "Malawi only. Salt iodisation can change iodine status quickly.")
about_content <- div(
  h5("About this dashboard", style = "margin-top: 0;"),
  p("District rankings of micronutrient deficiency for The Gambia, Ghana, Sierra Leone, Malawi and Cote d'Ivoire,",
    " made from public data and checked against national biomarker surveys. For ministries, funders and researchers",
    " deciding where to focus and where to survey next."),
  p("Method, tests and limits: see Technical notes."),
  h6("Citation"),
  p(em("Mertens et al. (in preparation). What geospatial covariates can and cannot do for sub-national",
       " micronutrient estimation: a protocol-first re-analysis of four national biomarker surveys.")),
  p(em(sprintf("Data built %s.", data_build_time)), style = "color: #888; font-size: 0.85em;"))
glossary_content <- div(
  h5("Glossary", style = "margin-top: 0;"),
  tags$dl(
    tags$dt("Priority score"), tags$dd("Place in the country's ranking: 100 = ranked worst."),
    tags$dt("Ranking accuracy"), tags$dd("Match between the model's order of districts and the survey's: 0 = random, 1 = the same."),
    tags$dt("Estimated prevalence"), tags$dd("The ranking converted to a percentage using the national survey figure."),
    tags$dt("Checked range"), tags$dd(sprintf("A range that held the survey's own figure %s of the time in tests.", pc(Q$cal_cov))),
    tags$dt("Rank range"), tags$dd("How far a rank moves when the model is re-estimated. Not a confidence interval."),
    tags$dt("Country left out"), tags$dd("A test where the model is built without that country, as for a country with no survey."),
    tags$dt("More"), tags$dd("Technical notes has full definitions.")))
