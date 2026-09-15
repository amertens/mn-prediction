# =============================================================================
# scripts/protocol_v2/59_build_addback_sources.R   [AB-01, 2026-09-15]
#
# Add-back sources: everything from the West Africa data-landscape review
# (Mertens et al., submission draft) that is on disk or openly fetchable and was
# NOT in the shared Admin-2 predictor set, plus two on-disk subnational layers
# the review does not catalogue. Nothing here is filtered for coverage; the
# consumer filters on the metadata (n_countries, subnational), as for scripts
# 07 and 08.
#
#   A. UNICEF vitamin A supplementation coverage (WDI mirror)   national, 4/4
#   B. FluNet influenza surveillance at the survey year           national, 2/4
#   C. GFDx fortification PROGRAMME fields (legislation in force at the survey
#      year, intake, industrial processing, standards, potential nutrient
#      delivery) - the part of GFDx script 08 did not take   national, 4/4
#   D. WHO/UNICEF Global Anaemia Estimates at the survey year     national, 4/4
#   E. Tang et al. 2026 (WFP / MIMI) HCES-based nutrient inadequacy, Ghana
#      Admin-1 broadcast to districts                             Ghana only
#   F. Global Data Lab subnational HDI at the survey year         Admin-1, 4/4
#   G. MICS immunisation indicators by subnational region from the WHO Health
#      Inequality Data Repository (HEAT), nearest MICS round   Admin-1, 4/4
#      (Malawi MICS 2014 at district level; the only MICS block that covers
#      Sierra Leone and Malawi without microdata)
#
# LEAKAGE POLICY (LK-02). A predictor is leakage only when it is measured on the
# outcome individuals (same survey instance). Every block here is an external
# source, so anaemia (D) and per-nutrient inadequacy (E) are admissible; both
# are flagged in `domain` so a sensitivity run can drop them
# (V2_DROP_MODELLED=1 removes every "MODELLED SURFACE" domain at fit time).
#
# National blocks are constant within a country (subnational = FALSE): they can
# move a country's level in a pooled or LOCO model and can never change a
# district ranking. Same convention as FAOSTAT and GFDx in scripts 07/08.
#
# Order in a rebuild: builder -> 07 -> 08 -> 59 (this) -> 53.
#
#   Rscript -e "source('scripts/protocol_v2/59_build_addback_sources.R')"
# -> data/covariates/harmonized/predictors_admin2_addback.csv
# -> data/covariates/harmonized/predictors_admin2_addback_metadata.csv
# -> updates predictors_admin2_shared.csv + _metadata.csv (.pre_addback backup)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

HDIR <- "data/covariates/harmonized"
source("R/survey_years.R"); source("R/admin2_keys.R"); SURVEY_YEAR <- survey_years()
COUNTRIES <- names(SURVEY_YEAR)
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(as.character(x))))
num <- function(x) suppressWarnings(as.numeric(x))

SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE,
               stringsAsFactors = FALSE)
PREV_ADDBACK_COLS <- if (file.exists(file.path(HDIR, "predictors_admin2_addback.csv")))
  setdiff(names(read.csv(file.path(HDIR, "predictors_admin2_addback.csv"), nrows = 1, check.names = FALSE)),
          c("country", "Admin1", "Admin2")) else character()
spine <- SH[, c("country", "Admin1", "Admin2")]
cat(sprintf("[addback] spine: %d Admin-2 rows, %d countries; shared set %d predictors\n",
            nrow(spine), length(unique(spine$country)), ncol(SH) - 3L))

blocks <- list(); meta <- list()
add_meta <- function(cols, source, domain, subnational, assumption) {
  if (!length(cols)) return(invisible())
  meta[[length(meta) + 1]] <<- data.frame(
    column = cols, source = source, domain = domain,
    subnational = subnational, assumption = assumption, stringsAsFactors = FALSE)
}
# nearest available year to the survey year, ties to the earlier year
pick_year <- function(years, target) {
  years <- years[is.finite(years)]; if (!length(years)) return(NA_real_)
  years[order(abs(years - target), years)][1]
}
match_names <- function(src, tgt, label) {
  # JK-01 (2026-09-15): one matcher for every name-keyed block - exact on the
  # normalised key, then the aliases file, then Jaro-Winkler within 0.15 - and
  # every decision is written to metadata/crosswalks/review_59_<label>.csv so
  # the fuzzy step can be checked line by line (the first version kept the
  # source spelling on exact hits and silently emptied the HEAT block; caught
  # by the audit).
  admin2_match_v2(src, tgt, aliases_csv = "metadata/crosswalks/aliases_59_addback.csv",
                  review_csv = file.path("metadata", "crosswalks", sprintf("review_59_%s.csv", gsub("[^A-Za-z0-9]+", "_", label))), label = label)
}
national_block <- function(df_country, cols) {
  # df_country: one row per country with `country` + cols -> broadcast to spine
  spine |> left_join(df_country[, c("country", cols)], by = "country")
}
# The Gambia: the spine's Admin1 are the six old divisions; DHS, MICS, GDL and
# most national statistics report the eight LGAs. Central River (Maccarthy
# Island) splits by river bank: Kuntaur LGA is the north bank (Saloums, Nianija,
# Niani, Sami), Janjanbureh LGA the south bank (Niaminas, Fulladu West,
# Janjanbureh); Banjul division holds both Banjul and Kanifing LGAs.
gambia_lga <- function(sp) {
  ifelse(sp$Admin2 == "Kanifing", "Kanifing",
  ifelse(sp$Admin1 == "Banjul", "Banjul",
  ifelse(sp$Admin1 == "Western", "Brikama",
  ifelse(sp$Admin1 == "Lower River", "Mansakonko",
  ifelse(sp$Admin1 == "North Bank", "Kerewan",
  ifelse(sp$Admin1 == "Upper River", "Basse",
  ifelse(sp$Admin2 %in% c("Lower Saloum", "Upper Saloum", "Nianija", "Niani", "Sami"), "Kuntaur",
  ifelse(sp$Admin2 %in% c("Fulladu West", "Janjanbureh", "Niamina Dankunku", "Niamina East", "Niamina West"),
         "Janjanbureh", NA_character_))))))))
}

# ── A. UNICEF vitamin A supplementation coverage (national, broadcast) ───────
cat("\n[A. UNICEF VAS]\n")
vas <- lapply(COUNTRIES, function(cn) {
  f <- file.path("data", "VAS", sprintf("%s_vas_admin2.csv", cn))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  # one national value per country (the builder wrote it to every district)
  data.frame(country = cn,
             vas_vita_coverage_pct_sy = num(d$vas_vita_coverage_pct[1]),
             vas_year_used = num(d$vas_year[1]))
}) |> bind_rows()
if (nrow(vas)) {
  cat(sprintf("  %s\n", paste(sprintf("%s=%.0f%% (%d)", vas$country, vas$vas_vita_coverage_pct_sy, vas$vas_year_used), collapse = ", ")))
  blocks$vas <- national_block(vas, "vas_vita_coverage_pct_sy")
  add_meta("vas_vita_coverage_pct_sy", "UNICEF VAS (World Bank WDI mirror SN.ITK.VITA.ZS)",
           "Food fortification and supplementation", FALSE,
           "NATIONAL two-dose vitamin A supplementation coverage of children 6-59 months at the year nearest the survey (scripts/build_vas_national.R), broadcast to every district. External programme coverage, not measured on the biomarker survey's respondents (LK-02).")
}

# ── B. FluNet at the survey year (national) ──────────────────────────────────
cat("\n[B. FluNet]\n")
fl_path <- "data/FluNet/VIW_FNT.csv"
if (file.exists(fl_path)) {
  F <- read.csv(fl_path, stringsAsFactors = FALSE)
  FL <- bind_rows(lapply(COUNTRIES, function(cn) {
    f <- F[F$COUNTRY_CODE == ISO[[cn]] & F$ISO_YEAR == SURVEY_YEAR[[cn]], ]
    if (!nrow(f)) return(NULL)
    sp <- sum(num(f$SPEC_PROCESSED_NB), na.rm = TRUE)
    pa <- sum(num(f$INF_A), na.rm = TRUE); pall <- sum(num(f$INF_ALL), na.rm = TRUE)
    data.frame(country = cn,
               flunet_specimens_per_week_sy = sp / nrow(f),
               flunet_share_positive_sy     = if (sp > 0) pall / sp else NA_real_,
               flunet_influenza_a_share_sy  = if (pall > 0) pa / pall else NA_real_)
  }))
  if (nrow(FL)) {
    print(FL, row.names = FALSE)
    fc <- setdiff(names(FL), "country")
    blocks$flunet <- national_block(FL, fc)
    add_meta(fc, "WHO FluNet (VIW_FNT export)", "Infectious disease surveillance", FALSE,
             "NATIONAL weekly virological surveillance summed over the survey year, broadcast to every district; reported for Ghana and Sierra Leone only (no Gambia or Malawi rows in FluNet). Constant within a country.")
  }
} else cat("  FluNet export absent\n")

# ── C. GFDx fortification programme fields (national, broadcast) ─────────────
cat("\n[C. GFDx programme]\n")
gf_path <- "data/GFDx/GFDxDataSet.csv"
if (file.exists(gf_path)) {
  gf <- suppressWarnings(readr::read_csv(gf_path, show_col_types = FALSE, progress = FALSE,
                                         guess_max = 100000)) |> as.data.frame()
  # REDCap export: country_name is written on the first row of each country only
  cc <- num(gf$country_code); ok <- !is.na(cc)
  pos <- cummax(ifelse(ok, seq_along(cc), 0L))                 # last non-missing position
  gf$cc <- ifelse(pos > 0L, cc[pmax(pos, 1L)], NA_real_)        # last observation carried forward
  gf <- gf |> group_by(cc) |>
    mutate(cn = { v <- country_name; v <- v[!is.na(v)]; if (length(v)) v[1] else NA_character_ }) |>
    ungroup() |> as.data.frame()
  GF_NAME <- c(Gambia = "Gambia, Republic of The", Ghana = "Ghana", Malawi = "Malawi",
               SierraLeone = "Sierra Leone")
  VEH <- c(wheat = "Wheat flour", maize = "Maize flour", oil = "Oil", salt = "Salt", rice = "Rice")
  # standard_nutrient codes from the GFDx indicator compendium (13 Apr 2022,
  # data/GFDx/GFDx-indicator-compendium-13-Apr-2022.docx): 1 B6, 2 B12,
  # 3 calcium, 4 fluoride, 5 folate, 6 iodine, 7 iron, 8 niacin, 9 riboflavin,
  # 10 selenium, 11 thiamin, 12 vitamin A, 13 vitamin D, 14 vitamin E, 15 zinc.
  # (An earlier version of this script inferred 8 = zinc from the target levels;
  # 8 is niacin. Corrected 2026-09-15 from the compendium the RA collected.)
  NUT <- c(`7` = "iron", `12` = "vita", `5` = "folic", `2` = "b12", `15` = "zinc", `6` = "iodine")
  rows <- list()
  for (cn in COUNTRIES) {
    Y <- SURVEY_YEAR[[cn]]
    s <- gf[!is.na(gf$cn) & gf$cn == GF_NAME[[cn]], , drop = FALSE]
    if (!nrow(s)) { cat("  ", cn, "not in GFDx\n"); next }
    rec <- list(country = cn)
    deliver <- setNames(rep(0, length(NUT)), NUT)
    for (v in names(VEH)) {
      e <- s[s$event_name == VEH[[v]], , drop = FALSE]
      base <- e[is.na(e$repeat_instrument), , drop = FALSE]
      status <- if (nrow(base)) num(base$status_food[1]) else NA_real_       # 1 mandatory, 2 voluntary, 3 none
      eff <- if (nrow(base)) num(base$effective_year[1]) else NA_real_
      if (!is.finite(eff) && nrow(base)) eff <- num(base$fortification_year[1])
      mand <- is.finite(status) && status == 1 && is.finite(eff) && eff <= Y
      rec[[sprintf("gfdx_%s_mandatory_sy", v)]] <- as.numeric(mand)
      rec[[sprintf("gfdx_%s_years_mandatory_sy", v)]] <- if (mand) Y - eff else 0
      it <- e[e$repeat_instrument %in% "intake" & is.finite(num(e$food_intake)), , drop = FALSE]
      yi <- pick_year(num(it$food_intake_year), Y)
      intake <- if (is.finite(yi)) num(it$food_intake[num(it$food_intake_year) == yi][1]) else NA_real_
      rec[[sprintf("gfdx_%s_intake_g_sy", v)]] <- intake
      ip <- e[e$repeat_instrument %in% "industrially_processed" & is.finite(num(e$industrially_processed_pc)), , drop = FALSE]
      yp <- pick_year(num(ip$ip_year), Y)
      ip_pc <- if (is.finite(yp)) num(ip$industrially_processed_pc[num(ip$ip_year) == yp][1]) else NA_real_
      rec[[sprintf("gfdx_%s_ip_pc_sy", v)]] <- ip_pc
      # standards -> potential delivery (mg/capita/day) = intake g/1000 kg x level mg/kg x industrially processed share
      nc <- e[e$repeat_instrument %in% "nutrients_compounds" & is.finite(num(e$standard_nutrient)), , drop = FALSE]
      for (code in names(NUT)) {
        lvl <- num(nc$nutrient_level[num(nc$standard_nutrient) == num(code)])
        lvl <- if (length(lvl) && any(is.finite(lvl))) max(lvl, na.rm = TRUE) else 0
        if (!mand) lvl <- 0
        if (is.finite(intake) && is.finite(ip_pc))
          deliver[[NUT[[code]]]] <- deliver[[NUT[[code]]]] + intake / 1000 * lvl * ip_pc / 100
      }
    }
    for (n in names(deliver)) rec[[sprintf("gfdx_%s_delivery_mg_sy", n)]] <- deliver[[n]]
    rec$gfdx_n_vehicles_mandatory_sy <- sum(unlist(rec[grep("^gfdx_[a-z]+_mandatory_sy$", names(rec))]))
    rows[[cn]] <- as.data.frame(rec, stringsAsFactors = FALSE)
  }
  G <- bind_rows(rows)
  if (nrow(G)) {
    gc <- setdiff(names(G), "country")
    # drop columns that are zero / NA everywhere (a vehicle no country fortifies)
    gc <- gc[vapply(gc, function(v) any(is.finite(G[[v]])) && length(unique(stats::na.omit(G[[v]]))) > 1, TRUE)]
    print(G[, c("country", intersect(c("gfdx_n_vehicles_mandatory_sy", "gfdx_iron_delivery_mg_sy",
                                       "gfdx_vita_delivery_mg_sy", "gfdx_folic_delivery_mg_sy"), gc))], row.names = FALSE)
    blocks$gfdx_prog <- national_block(G, gc)
    add_meta(gc, "GFDx (Global Fortification Data Exchange) programme fields",
             "Food fortification and supplementation", FALSE,
             "NATIONAL, broadcast to every district. mandatory_sy = mandatory legislation in force by the survey year; years_mandatory_sy = years since it took effect; intake_g_sy = g/capita/day of the vehicle at the nearest year; ip_pc_sy = share industrially processed at the nearest year; *_delivery_mg_sy = sum over vehicles of intake x standard level x industrially-processed share for vehicles with mandatory legislation in force (potential, not measured, delivery). Nutrient codes 7/12/5/2/15/6 = iron/vitamin A/folic acid/B12/zinc/iodine per the GFDx indicator compendium (13 Apr 2022, data/GFDx/).")
  }
} else cat("  GFDx dataset absent\n")

# ── D. WHO/UNICEF Global Anaemia Estimates at the survey year (national) ─────
cat("\n[D. WHO anaemia]\n")
wa_path <- "results/external/who_anaemia_country.csv"
if (!file.exists(wa_path)) {
  cat("  cache absent; run scripts/pull_who_unicef_anaemia.R first\n")
} else {
  W <- read.csv(wa_path, stringsAsFactors = FALSE)
  W <- W[W$pop_dim %in% c("SEVERITY_TOTAL", "SEX_BTSX", "") | is.na(W$pop_dim), ]
  W <- W[!duplicated(W[, c("iso3", "year", "indicator")]), ]
  IND <- c(children_6_59m = "who_anaemia_child_pct_sy", women_nonpreg = "who_anaemia_nonpreg_pct_sy",
           women_pregnant = "who_anaemia_pregnant_pct_sy")
  WA <- bind_rows(lapply(COUNTRIES, function(cn) {
    Y <- SURVEY_YEAR[[cn]]; w <- W[W$iso3 == ISO[[cn]], ]
    rec <- list(country = cn)
    for (i in names(IND)) {
      wi <- w[w$indicator == i, ]
      yy <- pick_year(num(wi$year), Y)
      rec[[IND[[i]]]] <- if (is.finite(yy)) num(wi$value_num[num(wi$year) == yy][1]) else NA_real_
      if (i %in% c("children_6_59m", "women_nonpreg")) {
        y5 <- pick_year(num(wi$year), Y - 5)
        v5 <- if (is.finite(y5)) num(wi$value_num[num(wi$year) == y5][1]) else NA_real_
        rec[[sub("_pct_sy$", "_trend5_pp", IND[[i]])]] <- rec[[IND[[i]]]] - v5
      }
    }
    as.data.frame(rec, stringsAsFactors = FALSE)
  }))
  wc <- setdiff(names(WA), "country")
  print(WA, row.names = FALSE)
  blocks$who_anaemia <- national_block(WA, wc)
  add_meta(wc, "WHO Global Anaemia Estimates (GHO, WHO/UNICEF joint estimates)",
           "Nutrition status (MODELLED SURFACE)", FALSE,
           "NATIONAL modelled anaemia prevalence at the survey year (children 6-59 months, women 15-49, pregnant, non-pregnant) and the 5-year change, broadcast to every district. External modelled estimate (LK-02 admissible); for Malawi the WHO model ingests the MDHS 2015-16 haemoglobin sample of which the MNS is a subsample - a second-order overlap that cannot be removed. Dropped under V2_DROP_MODELLED=1.")
}

# ── E. Tang et al. 2026 (WFP/MIMI) nutrient inadequacy, Ghana Admin-1 ────────
cat("\n[E. Tang 2026 / MIMI, Ghana]\n")
tg_path <- "data/WFP_LSFF_2026/tangS1_adm1_vulnerability.csv"
xw_path <- "data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv"
if (file.exists(tg_path) && file.exists(xw_path)) {
  TG <- read.csv(tg_path, stringsAsFactors = FALSE); TG <- TG[TG$country == "Ghana", ]
  XW <- read.csv(xw_path, stringsAsFactors = FALSE)
  gh <- spine[spine$country == "Ghana", ]
  reg10 <- XW$adm1_paper_10[match(kk(gh$Admin1), kk(XW$admin1_16))]
  reg10[is.na(reg10)] <- gh$Admin1[is.na(reg10)]
  j <- match(kk(reg10), kk(TG$adm1_paper))
  cat(sprintf("  Ghana districts mapped to a Tang region: %d of %d\n", sum(!is.na(j)), nrow(gh)))
  tc <- c(mimi_mpi = "mpi", mimi_vita_inadequate_pct = "vitA_pct", mimi_folate_inadequate_pct = "folate_pct",
          mimi_b12_inadequate_pct = "b12_pct", mimi_iron_inadequate_pct = "iron_pct", mimi_zinc_inadequate_pct = "zinc_pct")
  E <- gh
  for (n in names(tc)) E[[n]] <- num(TG[[tc[[n]]]][j])
  blocks$tang <- E
  add_meta(names(tc), "Tang et al. 2026 Nature Food (WFP / MIMI), Supplementary Table S1",
           "Dietary inadequacy (MODELLED SURFACE, HCES)", TRUE,
           "Ghana ONLY. Modelled probability of inadequate intake from GLSS 2016/17 apparent household consumption (Harmonized Average Requirement; iron full-probability, others cut-point), reported for the 10 pre-2018 regions and broadcast to districts through the 16-to-10 crosswalk. External HCES-based estimate (LK-02 admissible) and outcome-adjacent by construction: flagged MODELLED SURFACE. Using it as a predictor removes it as an independent Ghana check (docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md); the Cote d'Ivoire check is unaffected.")
} else cat("  Tang table or crosswalk absent\n")

# ── F. Global Data Lab subnational HDI at the survey year ────────────────────
cat("\n[F. GDL subnational HDI]\n")
gdl_path <- "data/external_cache/GDL-Subnational-HDI-data.csv"
if (file.exists(gdl_path)) {
  D <- read.csv(gdl_path, stringsAsFactors = FALSE)
  GDL_NAME <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
  XW <- if (file.exists(xw_path)) read.csv(xw_path, stringsAsFactors = FALSE) else NULL
  gdl_rows <- list()
  for (cn in COUNTRIES) {
    Y <- SURVEY_YEAR[[cn]]
    d <- D[D$Country == GDL_NAME[[cn]] & D$Level == "Subnat", ]
    if (!nrow(d)) { cat("  ", cn, "no subnational rows\n"); next }
    yy <- pick_year(unique(num(d$Year)), Y); y5 <- pick_year(unique(num(d$Year)), Y - 5)
    reg <- unique(d$Region)
    val <- vapply(reg, function(r) num(d$shdi[d$Region == r & num(d$Year) == yy][1]), 0)
    val5 <- vapply(reg, function(r) num(d$shdi[d$Region == r & num(d$Year) == y5][1]), 0)
    sp <- spine[spine$country == cn, ]
    # the level GDL reports at, per country: Gambia LGAs and Malawi districts are
    # the spine's Admin1 (Malawi groups some districts), Ghana's 10 old regions map
    # through the crosswalk, Sierra Leone's 14 districts are the spine's Admin2
    if (cn == "SierraLeone") {
      key <- match_names(sp$Admin2, reg, "Sierra Leone Admin2 -> GDL")
    } else if (cn == "Ghana" && !is.null(XW)) {
      r10 <- XW$adm1_paper_10[match(kk(sp$Admin1), kk(XW$admin1_16))]
      r10[is.na(r10)] <- sp$Admin1[is.na(r10)]
      key <- match_names(r10, reg, "Ghana Admin1(10) -> GDL")
    } else if (cn == "Malawi") {
      # grouped units list their districts in parentheses
      key <- vapply(sp$Admin1, function(a) {
        hit <- reg[vapply(reg, function(r) kk(a) == kk(sub(" \\(.*$", "", r)) ||
                                            grepl(kk(a), kk(r), fixed = TRUE), TRUE)]
        if (length(hit)) hit[1] else NA_character_ }, "")
      cat(sprintf("    %-28s %3d matched, %3d unmatched of %3d\n", "Malawi Admin1 -> GDL",
                  sum(!is.na(key)), sum(is.na(key)), length(key)))
    } else if (cn == "Gambia") {
      key <- match_names(gambia_lga(sp), reg, "Gambia district -> LGA -> GDL")
    } else {
      key <- match_names(sp$Admin1, reg, paste(cn, "Admin1 -> GDL"))
    }
    sp$gdl_shdi_sy <- unname(val[key]); sp$gdl_shdi_trend5 <- unname(val[key] - val5[key])
    gdl_rows[[cn]] <- sp
  }
  if (length(gdl_rows)) {
    blocks$gdl <- bind_rows(gdl_rows)
    add_meta(c("gdl_shdi_sy", "gdl_shdi_trend5"), "Global Data Lab Subnational HDI (SHDI)",
             "Education, employment, SES", TRUE,
             "Subnational HDI at the year nearest the survey and its 5-year change, at the level GDL reports (Gambia LGAs and Malawi districts = Admin1, Ghana 10 old regions via the 16-to-10 crosswalk, Sierra Leone districts = Admin2), broadcast to districts within the unit. Not in the landscape review; on disk from an earlier download.")
  }
} else cat("  GDL cache absent\n")

# ── G. MICS by subnational region, WHO Health Inequality Data Repository ─────
# The HEAT export (data/WHO_HEAT/*.xlsx, filtered to the four countries and the
# "Subnational region" dimension in heat_subnational_4countries.csv) carries
# MICS estimates for every country and round, which is the only MICS source on
# disk for Sierra Leone (2010, 2017) and Malawi (2014, at district level). Only
# the immunisation dataset holds MICS rows for these countries, so the block is
# 18 immunisation indicators; the full MICS harmonisation still needs microdata.
cat("\n[G. HEAT MICS by region]\n")
heat_path <- "data/WHO_HEAT/heat_subnational_4countries.csv"
if (file.exists(heat_path)) {
  H <- read.csv(heat_path, stringsAsFactors = FALSE)
  H <- H[grepl("MICS", H$source) & H$dimension == "Subnational region", ]
  heat_rows <- list(); heat_cols <- character()
  for (cn in COUNTRIES) {
    Y <- SURVEY_YEAR[[cn]]; h <- H[H$iso3 == ISO[[cn]], ]
    if (!nrow(h)) { cat("  ", cn, "no MICS rows\n"); next }
    yr <- pick_year(unique(num(h$date)), Y); h <- h[num(h$date) == yr, ]
    w <- tidyr::pivot_wider(h[, c("subgroup", "indicator_abbr", "estimate")],
                            names_from = "indicator_abbr", values_from = "estimate",
                            values_fn = function(z) mean(z, na.rm = TRUE)) |> as.data.frame()
    ind <- setdiff(names(w), "subgroup"); names(w)[match(ind, names(w))] <- paste0("mics_heat_", ind, "_sy")
    reg <- w$subgroup; sp <- spine[spine$country == cn, ]
    key <- if (cn == "Gambia") match_names(gambia_lga(sp), reg, "Gambia LGA -> MICS region") else
      if (cn == "Ghana" && !is.null(XW)) {
        r10 <- XW$adm1_paper_10[match(kk(sp$Admin1), kk(XW$admin1_16))]; r10[is.na(r10)] <- sp$Admin1[is.na(r10)]
        match_names(r10, reg, "Ghana Admin1(10) -> MICS region")
      } else if (cn == "SierraLeone") {
        # HEAT labels the four provinces north / east / south / west
        match_names(sub("ern$", "", tolower(sp$Admin1)), reg, "Sierra Leone province -> MICS region")
      } else if (cn == "Malawi") {
        # district-level MICS 2014; the three regions in the same file are skipped
        match_names(sp$Admin1, setdiff(reg, c("Central", "Northern", "Southern")), "Malawi district -> MICS region")
      } else match_names(sp$Admin1, reg, paste(cn, "Admin1 -> MICS region"))
    j <- match(key, reg)
    for (v in setdiff(names(w), "subgroup")) sp[[v]] <- num(w[[v]][j])
    sp$mics_heat_round_year <- yr
    cat(sprintf("  %-12s MICS %d, %d indicators, %d regions\n", cn, yr, length(ind), length(reg)))
    heat_rows[[cn]] <- sp; heat_cols <- union(heat_cols, setdiff(names(w), "subgroup"))
  }
  if (length(heat_rows)) {
    HB <- bind_rows(heat_rows); HB$mics_heat_round_year <- NULL
    blocks$heat_mics <- HB
    add_meta(heat_cols, "WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region",
             "Infant and child morbidity/mortality", TRUE,
             "MICS estimate for the subnational region at the round nearest the survey (Gambia 2018 LGAs, Ghana 2017 old regions, Sierra Leone 2010 provinces, Malawi 2014 districts), broadcast to the districts inside the region. Only the HEAT immunisation dataset carries MICS rows for these countries, so the block is immunisation coverage; a different survey programme from DHS, measured on different people.")
  }
} else cat("  HEAT extract absent (see docs/findings/SANDBOX_LOG_2026-09.md AB-01 for the download URLs)\n")

# ── assemble and append, without the four-country filter ────────────────────
if (!length(blocks)) stop("no blocks built")
AB <- spine
for (b in names(blocks)) {
  by <- intersect(c("country", "Admin1", "Admin2"), names(blocks[[b]]))
  j <- tryCatch(left_join(AB, blocks[[b]], by = by), error = function(e) NULL)
  if (is.null(j) || nrow(j) != nrow(AB)) { cat("  [join]", b, "changed row count or failed, skipped\n"); next }
  AB <- j
}
newcols <- setdiff(names(AB), c("country", "Admin1", "Admin2"))
empty <- newcols[vapply(newcols, function(v) !any(is.finite(AB[[v]])), TRUE)]
if (length(empty)) { cat("  dropping all-missing:", paste(empty, collapse = ", "), "\n")
  AB <- AB[, setdiff(names(AB), empty), drop = FALSE]; newcols <- setdiff(newcols, empty) }
# value-identical columns carry nothing twice (e.g. a vehicle mandatory for one
# year in one country makes years_mandatory == mandatory); keep the first
key <- vapply(newcols, function(v) paste(signif(AB[[v]], 7), collapse = "|"), "")
dup <- newcols[duplicated(key)]
if (length(dup)) { cat("  dropping value-identical duplicates:", paste(dup, collapse = ", "), "\n")
  AB <- AB[, setdiff(names(AB), dup), drop = FALSE]; newcols <- setdiff(newcols, dup) }
write.csv(AB, file.path(HDIR, "predictors_admin2_addback.csv"), row.names = FALSE)

MD <- bind_rows(meta) |> filter(column %in% newcols) |> distinct(column, .keep_all = TRUE)
MD$n_countries <- vapply(MD$column, function(v)
  sum(tapply(AB[[v]], AB$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
MD$countries <- vapply(MD$column, function(v) {
  s <- tapply(AB[[v]], AB$country, function(z) any(is.finite(z)))
  paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
MD$completeness <- round(vapply(MD$column, function(v) mean(is.finite(AB[[v]])), 0), 3)
MD$coverage_by_country <- vapply(MD$column, function(v) {
  s <- AB |> group_by(country) |> summarise(ok = mean(is.finite(.data[[v]])), .groups = "drop")
  paste(sprintf("%s=%.2f", s$country, s$ok), collapse = ";") }, "")
write.csv(MD, file.path(HDIR, "predictors_admin2_addback_metadata.csv"), row.names = FALSE)

SHM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
SHM$subnational <- as.logical(toupper(as.character(SHM$subnational)))
for (f in c("predictors_admin2_shared.csv", "predictors_admin2_shared_metadata.csv"))
  if (!file.exists(file.path(HDIR, paste0(f, ".pre_addback"))))
    file.copy(file.path(HDIR, f), file.path(HDIR, paste0(f, ".pre_addback")))
# columns of an EARLIER add-back run that this run no longer produces (a
# duplicate dropped, a block withdrawn) must leave the shared set too
stale <- setdiff(PREV_ADDBACK_COLS, newcols)
if (length(stale)) cat("  removing stale add-back columns:", paste(stale, collapse = ", "), "\n")
SH2 <- SH |> select(-any_of(c(newcols, stale))) |> left_join(AB, by = c("country", "Admin1", "Admin2"))
SHM2 <- bind_rows(SHM |> filter(!column %in% c(newcols, stale)),
                  MD |> transmute(column, domain, source, n_countries, countries, completeness, subnational))
stopifnot(nrow(SH2) == nrow(SH), ncol(SH2) - 3L == nrow(SHM2))
write.csv(SH2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(SHM2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"), row.names = FALSE)

cat("\n=== added ===\n")
print(as.data.frame(MD[, c("column", "source", "n_countries", "completeness", "subnational")]), row.names = FALSE)
cat("\nshared set:", ncol(SH) - 3, "->", ncol(SH2) - 3, "predictors\n")
cat("DONE\n")
