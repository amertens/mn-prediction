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
                   child_zinc = "Zinc, children", women_zinc = "Zinc, women")
arm_label <- c(domain_index = "Proxy index (this dashboard)", spatial_plus_domain = "Neighbour smoother + proxies",
               spatial = "Neighbour smoother alone", domain_enet = "Penalised regression, domain components",
               raw_enet = "Penalised regression, all columns", region_mean_jk = "Survey's own regional average",
               null_train_mean = "No information (training mean)")

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
Q$n_predictors <- if (!is.null(CAT)) nrow(CAT$variables) else 454
Q$n_domains    <- if (!is.null(CAT)) nrow(CAT$domains) else 24
Q$n_districts  <- length(unique(paste(idx_districts$country, idx_districts$Admin1, idx_districts$Admin2)))
Q$n_surveyed   <- sum(idx_national$n_surveyed[!duplicated(idx_national$country)])
f2 <- function(x, d = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")
pc <- function(x, d = 0) ifelse(is.finite(x), paste0(formatC(100 * x, format = "f", digits = d), "%"), "NA")

# ── Colours ────────────────────────────────────────────────────────────────
who_colors <- c("Low" = "#2c7bb6", "Mild" = "#abd9e9", "Moderate" = "#fdae61",
                "Severe" = "#d7191c", "No data" = "#cccccc")
PROXY_COL <- "#0F7B8A"; SURVEY_COL <- "#C8641E"; GREY_COL <- "#8c8c8c"

# ── Caveats ────────────────────────────────────────────────────────────────
biomarker_caveats <- list(
  women_vitA   = paste("Vitamin A in women rests on retinol-binding protein, a weaker marker in women",
                       "and one moved by inflammation. Prevalence is under 3 percent in every survey, so the",
                       "district ranking has little to work with; the reliability ceiling for this outcome is low."),
  child_vitA   = paste("Vitamin A is measured by retinol-binding protein, adjusted for inflammation and converted",
                       "to retinol with each survey's own calibration line before the 0.70 cut-off."),
  women_b12    = paste("B12 is measured by serum B12, a marker of limited specificity, and in three countries only.",
                       "Read district differences as indicative."),
  women_folate = paste("Folate is measured in three countries only. National prevalence ranges from 19 to 79 percent",
                       "across them, so the level, more than the ranking, is the fragile part."),
  child_zinc   = paste("Zinc is measured in Malawi only, so nothing about it can be checked across borders. The",
                       "district differences are largely the time of the blood draw, not geography."),
  women_zinc   = paste("Zinc is measured in Malawi only, so nothing about it can be checked across borders. The",
                       "district differences are largely the time of the blood draw, not geography.")
)
GENERAL_CAVEAT <- paste(
  "The four surveys were run between 2013 and 2018 by different teams, so levels are not comparable",
  "across countries. The model carries rankings across borders; it does not carry levels, and a",
  "prevalence figure needs a national survey number to anchor it. Rankings in a country with no",
  "survey are rougher than inside a surveyed country; How well it works has the numbers.")

# ── Site banner ────────────────────────────────────────────────────────────
SITE_SCOPE_HEADLINE <- "Working estimates from public data, scored against four biomarker surveys."
SITE_SCOPE_BODY <- paste("The model ranks districts; it does not measure prevalence. Use the ranking, and read",
                         "any percentage as a planning figure anchored to the national survey.")
SITE_SCOPE_POINTER <- "How well it works shows what the ranking reaches and where it stops."
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
  p("District rankings of micronutrient deficiency for The Gambia, Ghana, Sierra Leone and Malawi,",
    " built from public data and scored against the four national biomarker surveys, with a ranking",
    " for Cote d'Ivoire, which has no survey. For ministries, funders and researchers deciding where",
    " to look first and where the next survey should sample."),
  h6("Method"),
  p(sprintf(paste("To every district we attach %d public data layers in %d groups: satellite imagery, climate,",
                  "soil, crops, livestock, malaria and other morbidity, prices and household-survey aggregates.",
                  "Each group is summarised into a few axes, weighted by how well it tracked deficiency where",
                  "blood was drawn, and summed. Nothing is tuned. Every district is scored with itself hidden",
                  "from the model, every country with the whole country hidden, and every method is compared",
                  "with the survey's own regional average, a neighbour smoother, and a permutation null."),
            Q$n_predictors, Q$n_domains)),
  h6("What it reaches"),
  tags$ul(
    tags$li(sprintf("Inside a surveyed country: ranking accuracy %s against %s for the survey's own regional averages and %s for chance.",
                    f2(Q$infill), f2(Q$infill_jk), f2(Q$null_d))),
    tags$li(sprintf("In a country never used in training: %s with everything, %s from climate and soil alone, positive in %s of %s country-outcome pairs.",
                    f2(Q$tr), f2(Q$cs), Q$tr_pos, Q$tr_n)),
    tags$li(sprintf("A perfect predictor could reach about %s given the survey's own noise; the model is about two-thirds of the way.",
                    f2(Q$ceiling_level)))
  ),
  h6("Citation"),
  p(em("Mertens et al. (in preparation). What geospatial covariates can and cannot do for sub-national",
       " micronutrient estimation: a protocol-first re-analysis of four national biomarker surveys.")),
  h6("Surveys"),
  p("The Gambia 2018, Ghana 2017, Sierra Leone 2013, Malawi 2015 to 2016. Vitamin A and iron in children",
    " and women; folate and B12 in women in three countries; zinc in Malawi."),
  p(em(sprintf("Data build: %s (%s)", data_build_time, PROTOCOL_LABEL)), style = "color: #888; font-size: 0.85em;")
)

glossary_content <- div(
  h5("Glossary", style = "margin-top: 0;"),
  tags$dl(
    tags$dt("Priority score"),
    tags$dd("Where a district sits in its country's ranking, from 100 (ranked worst) to near 0 (ranked best).",
            " It is the model's output. It says which districts are likely worst, not how bad."),
    tags$dt("Ranking accuracy"),
    tags$dd("How closely the model's order of districts matches the survey's, from 0 to 1 (Spearman correlation).",
            " Differences under 0.03 are ties over a few dozen districts."),
    tags$dt("Chance level"),
    tags$dd(sprintf("What a ranking reaches with no information, measured by shuffling the outcome: %s across districts, %s across regions.",
                    f2(Q$null_d), f2(Q$null_a1))),
    tags$dt("Reliability ceiling"),
    tags$dd("The best ranking accuracy a perfect predictor could reach, because the survey's own district values rest",
            " on one or two clusters and are themselves noisy."),
    tags$dt("Chance of being in the worst fifth"),
    tags$dd("The model was refitted on 40 random splits with the district hidden each time; this is how often the",
            " district landed in the worst fifth. Computed for surveyed districts, where there is a survey to check against."),
    tags$dt("Planning prevalence"),
    tags$dd("The ranking turned into a percentage by anchoring it to the country's national survey prevalence.",
            " The order comes from the model; the level comes from the survey."),
    tags$dt("In-fill, region, transport"),
    tags$dd("The three tests: a district hidden inside a surveyed country; a whole region hidden; a whole country hidden."),
    tags$dt("Percentage points (pp)"),
    tags$dd("The plain gap between two percentages. From 20 to 23 percent is 3 points."),
    tags$dt("District and region"),
    tags$dd("District is the second administrative level; region the first. Rankings are made for districts.")
  )
)
