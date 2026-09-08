# =============================================================================
# scripts/protocol_v2/08_build_extra_sources.R
#
# Wire the unused data sources into the modelling vocabulary.
#
# THE PROBLEM
# -----------
# Several substantial sources sit in data/ and reach ZERO modelled predictors.
# The largest is IHME: 93 GB across 13 topic folders (anaemia, child growth
# failure, double burden, exclusive breastfeeding, education, WASH, u5
# mortality, diarrhoea treatment, measles) whose harmonised columns were logged
# as "kept" and then removed by the four-country filter, because their Admin-2
# completeness fell below 50 percent in at least one country. HFID is
# misleadingly named: it is not health facilities but humanitarian
# food-insecurity data, carrying Food Consumption Score and reduced Coping
# Strategies Index at Admin-2 - direct dietary-quality measures - of which only
# two IPC PHASE columns were ever taken, and those were dropped as well.
#
# A DEFECT THIS FIXES, NOT JUST A GAP
# -----------------------------------
# IHME's adm2_name is NOT Admin-2 in every country. For Malawi it is the 27
# DISTRICTS, which are this project's Admin-1; the project's Admin-2 is the
# ~239 Traditional Authorities. Measured: 27 of 30 IHME Malawi names match
# Admin-1 exactly and 0 match Admin-2. The existing harmonisation fuzzy-matched
# those district names onto TA names at Jaro-Winkler 0.15 and "recovered" 23 of
# 30 - every one a silent mis-linkage of a district value onto an unrelated TA.
# This script matches each country at the level its names actually live on and
# broadcasts down where the source is coarser, recording that as an assumption
# instead of hiding it inside a fuzzy join.
#
# ASSUMPTIONS, RECORDED PER COLUMN IN THE METADATA
#  1 YEAR. IHME and HFID values are taken at the year nearest each country's
#    survey (Gambia 2021, Ghana 2017, Malawi 2015, Sierra Leone 2013) rather
#    than one fixed vintage, so the predictor is contemporaneous with the
#    outcome it is asked to explain.
#  2 STRATA. IHME rows are population-weighted across age and sex strata where
#    pop is present, and simple-averaged where it is not.
#  3 LEVEL. Gambia, Ghana and Sierra Leone join at Admin-2. Malawi joins at
#    Admin-1 and is broadcast to its Traditional Authorities: within-district
#    variation is not available from IHME and is not invented.
#  4 NAMES. Exact match on a normalised name first, Jaro-Winkler <= 0.15 for
#    the remainder, and only within the correct administrative level. The match
#    rate is reported per country rather than assumed.
#  5 NATIONAL SOURCES. GFDx fortification is country-level and broadcast to
#    every district, flagged subnational = FALSE, exactly as FAOSTAT is.
#
# The four-country filter is NOT applied. Per-country coverage is recorded so a
# consumer decides - the same correction made for the food environment.
#
#   Rscript scripts/protocol_v2/08_build_extra_sources.R
# -> data/covariates/harmonized/predictors_admin2_extrasrc.csv
# -> data/covariates/harmonized/predictors_admin2_extrasrc_metadata.csv
# -> updates predictors_admin2_shared.csv + _metadata.csv (.pre_extrasrc backup)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

HDIR <- "data/covariates/harmonized"
SURVEY_YEAR <- c(Gambia = 2018, Ghana = 2017, Malawi = 2015, SierraLeone = 2013)   # Gambia fieldwork Jan-Apr 2018 (FW-01); was 2021 until 2026-09-04
COUNTRIES <- names(SURVEY_YEAR)
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))

SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
spine <- SH[, c("country", "Admin1", "Admin2")]

match_names <- function(src, tgt, label) {
  out <- as.character(src)
  hit <- kk(out) %in% kk(tgt)
  n_exact <- sum(hit)
  if (any(!hit) && requireNamespace("stringdist", quietly = TRUE) && length(tgt)) {
    dm <- stringdist::stringdistmatrix(kk(out[!hit]), kk(tgt), method = "jw", p = 0.1)
    best <- apply(dm, 1, which.min); dist <- apply(dm, 1, min)
    out[!hit] <- ifelse(dist <= 0.15, tgt[best], NA_character_)
  } else if (any(!hit)) out[!hit] <- NA_character_
  cat(sprintf("    %-26s %3d exact, %3d after fuzzy, %3d unmatched of %3d\n",
              label, n_exact, sum(!is.na(out)), sum(is.na(out)), length(src)))
  out
}

blocks <- list(); meta <- list()
add_meta <- function(cols, source, domain, subnational, assumption) {
  if (!length(cols)) return(invisible())
  meta[[length(meta) + 1]] <<- data.frame(
    column = cols, source = source, domain = domain,
    subnational = subnational, assumption = assumption,
    stringsAsFactors = FALSE)
}

# ── 1. IHME ────────────────────────────────────────────────────────────────
cat("\n[IHME]\n")
ihme_files <- c(Gambia = "IHME_Gambia_data.dta", Ghana = "IHME_Ghana_data.dta",
                Malawi = "IHME_Malawi_data.dta",
                SierraLeone = "IHME_Sierra_Leone_data.dta")
# the level at which each country's IHME adm2_name actually lives, measured
IHME_LEVEL <- c(Gambia = "Admin2", Ghana = "Admin2",
                Malawi = "Admin1", SierraLeone = "Admin2")
ihme_blocks <- list()
for (cn in COUNTRIES) {
  p <- file.path("data", "IHME", ihme_files[[cn]])
  if (!file.exists(p)) { cat("  ", cn, "missing\n"); next }
  d <- haven::read_dta(p)
  d <- d[d$adm_level == 2, ]
  d$measure <- as.character(d$measure)
  d$metric  <- as.character(d$metric)
  d <- d[d$metric %in% c("Percent", "Prevalence", "Rate") |
           grepl("prevalence|proportion|coverage|attainment", d$measure,
                 ignore.case = TRUE), ]
  d <- d[is.finite(d$mean), ]
  if (!nrow(d)) { cat("  ", cn, "no usable rows\n"); next }
  # NEAREST YEAR PER INDICATOR, not per country. The series end in different
  # years (child growth failure, education, WASH, ORS stop at 2017; anaemia,
  # EBF, MCV run to 2019). One year per country kept only the indicators whose
  # series reached that year: with 2021 as Gambia's year the pick was 2019 and
  # 18 of Gambia's 29 indicators went all-NA (found 2026-09-04).
  d$.key <- paste0("ihme_", kk(d$measure))
  d$year <- suppressWarnings(as.numeric(d$year)); d <- d[is.finite(d$year), ]
  yp <- d |> group_by(.key) |>
    summarise(.pick = year[which.min(abs(year - SURVEY_YEAR[[cn]]))], .groups = "drop")
  d <- d |> inner_join(yp, by = ".key") |> filter(year == .pick)
  pick <- paste(range(d$.pick), collapse = "-")
  d$.w <- if ("pop" %in% names(d)) suppressWarnings(as.numeric(d$pop)) else 1
  d$.w[!is.finite(d$.w) | d$.w <= 0] <- 1
  agg <- d |> group_by(.a = as.character(adm2_name), .key) |>
    summarise(v = sum(mean * .w) / sum(.w), .groups = "drop")
  w <- tidyr::pivot_wider(agg, names_from = ".key", values_from = "v") |>
    as.data.frame()
  lvl <- IHME_LEVEL[[cn]]
  tgt <- unique(spine[[lvl]][spine$country == cn])
  cat(sprintf("  %-12s years %s, %3d indicators, join at %s\n",
              cn, pick, ncol(w) - 1, lvl))
  w$.match <- match_names(w$.a, tgt, paste(cn, lvl))
  w <- w[!is.na(w$.match), , drop = FALSE]
  if (!nrow(w)) next
  w[[lvl]] <- w$.match
  w$.a <- NULL; w$.match <- NULL
  # COLLAPSE BEFORE JOINING. Several IHME names can fuzzy-match onto the same
  # GADM unit (Ghana: 280 matched names onto 260 targets), and joining without
  # collapsing fans the spine - the defect class this project has already fixed
  # a dozen times. Average the duplicates and say how many there were.
  ndup <- sum(duplicated(w[[lvl]]))
  if (ndup) {
    w <- w |> group_by(.data[[lvl]]) |>
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    cat(sprintf("    collapsed %d duplicate %s matches by averaging\n", ndup, lvl))
  }
  w$country <- cn
  ihme_blocks[[cn]] <- spine[spine$country == cn, ] |>
    left_join(w, by = c("country", lvl))
}
if (length(ihme_blocks)) {
  IH <- bind_rows(ihme_blocks)
  ihcols <- setdiff(names(IH), c("country", "Admin1", "Admin2"))
  keep <- ihcols[vapply(ihcols, function(v)
    sum(tapply(IH[[v]], IH$country, function(z) any(is.finite(z))), na.rm = TRUE) >= 2,
    TRUE)]
  IH <- IH[, c("country", "Admin1", "Admin2", keep), drop = FALSE]
  # RASTER SWAP (IH-02, script 48). Where a 5 km IHME surface exists AND its
  # zonal mean agrees with the name-joined rollup where both exist (Spearman
  # >= 0.5 over matched districts), use the WorldPop-weighted zonal mean: no
  # name matching (Ghana's 2019 splits were leaving 26% of surveyed districts
  # empty), and true Admin-2 resolution in Malawi instead of a district value
  # broadcast to its Traditional Authorities. The diarrhoea rate surfaces do
  # not reproduce the tabular incidence / deaths columns (rho ~ 0) and stay
  # tabular, as do the six indicators with no surface.
  swap <- character(0)
  rf <- file.path(HDIR, "predictors_admin2_ihme_raster.csv")
  rm <- file.path(HDIR, "predictors_admin2_ihme_raster_metadata.csv")
  if (file.exists(rf) && file.exists(rm)) {
    IR <- read.csv(rf, check.names = FALSE); RM <- read.csv(rm, stringsAsFactors = FALSE)
    okcol <- unique(RM$column[is.finite(RM$rho_all) & RM$rho_all >= 0.5])
    swap <- intersect(intersect(okcol, names(IR)), keep)
    if (length(swap)) {
      IH <- IH |> select(-all_of(swap)) |>
        left_join(IR[, c("country", "Admin1", "Admin2", swap)], by = c("country", "Admin1", "Admin2"))
      IH <- IH[, c("country", "Admin1", "Admin2", keep), drop = FALSE]
    }
    writeLines(sprintf("  -> %d IHME columns taken from raster zonal means (script 48); tabular kept for: %s",
                       length(swap), paste(setdiff(keep, swap), collapse = ", ")))
  }
  blocks$ihme <- IH
  # anaemia and growth-failure surfaces are MODEL OUTPUTS estimating something
  # close to the outcome; flagged so they can be excluded from primary models.
  near <- grepl("anemia|anaemia|stunting|wasting|underweight", keep)
  add_meta(keep, "IHME (modelled surfaces)",
           ifelse(near, "Nutrition status (MODELLED SURFACE)",
                  "Infant and child morbidity/mortality"),
           TRUE,
           paste0(ifelse(keep %in% swap, "ZONAL MEAN of the 5 km IHME GeoTIFF surface over the spine polygon, WorldPop-weighted, nearest year to the survey (script 48). ", ""),
           ifelse(near,
                  "IHME modelled Admin-2 surface, year nearest survey, population-weighted over age and sex. MODELLED SURFACE: an estimate produced by someone else's model, not a direct measurement. Provenance checked: not fitted to the surveys scored here. Included by default; V2_DROP_MODELLED=1 excludes for sensitivity.",
                  "IHME modelled Admin-2 surface, year nearest survey, population-weighted over age and sex. Malawi joins at Admin-1 (its adm2_name is the district) and is broadcast to Traditional Authorities.")))
  cat(sprintf("  -> %d IHME columns kept (%d flagged NEAR-OUTCOME)\n",
              length(keep), sum(near)))
}

# ── 2. HFID food security (FCS, rCSI, IPC phase) ────────────────────────────
cat("\n[HFID food security]\n")
hf <- tryCatch(suppressWarnings(readr::read_csv("data/HFID/hfid_hv1.csv",
                 show_col_types = FALSE, progress = FALSE)) |> as.data.frame(),
               error = function(e) NULL)
if (!is.null(hf)) {
  nm <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi",
          SierraLeone = "Sierra Leone")
  hf$year <- suppressWarnings(as.integer(substr(as.character(hf$year_month), 1, 4)))
  vars <- intersect(c("ipc_phase_fews", "ipc_phase_ipcch", "fcs_lit", "rcsi_lit",
                      "fcs_rt mean", "rcsi_rt mean"), names(hf))
  hb <- list()
  for (cn in COUNTRIES) {
    s <- hf[!is.na(hf$ADMIN0) & hf$ADMIN0 == nm[[cn]] & !is.na(hf$ADMIN2), ]
    if (!nrow(s)) next
    s <- s[is.finite(s$year) & abs(s$year - SURVEY_YEAR[[cn]]) <= 2, ]
    if (!nrow(s)) { cat("  ", cn, "no rows in the survey window\n"); next }
    for (v in vars) s[[v]] <- suppressWarnings(as.numeric(s[[v]]))
    a <- s |> group_by(.a = as.character(ADMIN2)) |>
      summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    tgt <- unique(spine$Admin2[spine$country == cn])
    a$.match <- match_names(a$.a, tgt, paste(cn, "Admin2"))
    a <- a[!is.na(a$.match), , drop = FALSE]
    if (!nrow(a)) next
    a$Admin2 <- a$.match; a$.a <- NULL; a$.match <- NULL
    if (sum(duplicated(a$Admin2))) {
      a <- a |> group_by(Admin2) |>
        summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
        as.data.frame()
    }
    a$country <- cn
    hb[[cn]] <- a
  }
  if (length(hb)) {
    HF <- bind_rows(hb)
    vcols <- setdiff(names(HF), c("country", "Admin2"))
    names(HF)[match(vcols, names(HF))] <-
      paste0("fsec_", gsub("[ ]+", "_", vcols))
    blocks$hfid <- HF
    hfc <- setdiff(names(HF), c("country", "Admin2"))
    add_meta(hfc, "HFID (FEWS NET / IPC / WFP mVAM)", "Food prices and supply",
             TRUE,
             "Admin-2 food-security indicator averaged over the survey year +/-2 years. fcs = Food Consumption Score, rcsi = reduced Coping Strategies Index, ipc_phase = IPC/CH phase 1-5. Higher fcs is better; higher rcsi and ipc_phase are worse.")
    cat(sprintf("  -> %d HFID columns\n", length(hfc)))
  }
}

# ── 3. GFDx fortification (national, broadcast) ─────────────────────────────
cat("\n[GFDx fortification]\n")
gf <- tryCatch(suppressWarnings(readr::read_csv("data/GFDx/GFDxDataSet.csv",
                 show_col_types = FALSE, progress = FALSE)) |> as.data.frame(),
               error = function(e) NULL)
if (!is.null(gf) && "country_name" %in% names(gf)) {
  nm <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi",
          SierraLeone = "Sierra Leone")
  sub <- gf[gf$country_name %in% nm, , drop = FALSE]
  num <- names(sub)[vapply(sub, is.numeric, TRUE)]
  num <- setdiff(num, grep("year|code|instance|_complete$|population",
                           num, value = TRUE))
  keepn <- num[vapply(num, function(v) {
    s <- tapply(sub[[v]], sub$country_name, function(z) any(is.finite(z)))
    sum(s, na.rm = TRUE) >= 3 &&
      length(unique(stats::na.omit(sub[[v]]))) > 1
  }, TRUE)]
  keepn <- head(keepn, 12)
  if (length(keepn)) {
    g <- sub[, c("country_name", keepn), drop = FALSE] |>
      group_by(country_name) |>
      summarise(across(all_of(keepn), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    g$country <- names(nm)[match(g$country_name, nm)]
    g$country_name <- NULL
    vn <- setdiff(names(g), "country")
    names(g)[match(vn, names(g))] <- paste0("gfdx_", vn)
    blocks$gfdx <- spine |> left_join(g, by = "country")
    gnear <- grepl("anemia|anaemia|zinc_def|deficien", vn)
    add_meta(paste0("gfdx_", vn), "GFDx",
             ifelse(gnear, "Nutrition status (MODELLED SURFACE)",
                    "Food fortification and supplementation"), FALSE,
             ifelse(gnear,
               "GFDx NATIONAL prevalence broadcast to every district; no within-country variation. Anaemia is WHO 2011 (pre-dates every survey here); zinc is Wessells & Brown 2012 from FAO food-balance-sheet availability, so neither is fitted to these surveys.",
               "Global Fortification Data Exchange, NATIONAL value broadcast to every district: no within-country variation. Kept for cross-country (LOCO) information only."))
    cat(sprintf("  -> %d GFDx columns (national, broadcast)\n", length(vn)))
  } else cat("  no usable numeric GFDx columns for these countries\n")
} else cat("  GFDx unreadable or unexpected schema\n")


# ── 4. GEE demography and urbanisation (WorldPop age-sex, GHS-SMOD) ─────────
# Built by scripts/protocol_v2/09_extract_gee_demography.py. The vocabulary had
# 428 predictors and not one described WHO LIVES THERE, although age and sex
# composition is the most direct area-level determinant of nutritional
# REQUIREMENT. WorldPop age-sex had been used only for post-stratification
# weights, never as a predictor.
cat("
[GEE demography]
")
dem <- tryCatch(read.csv("data/covariates/harmonized/gee_demography_admin2.csv",
                         stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(dem)) {
  blocks$demog <- dem
  dcols <- setdiff(names(dem), c("country", "Admin1", "Admin2"))
  add_meta(dcols, "WorldPop age-sex / JRC GHS-SMOD (Earth Engine)",
           ifelse(grepl("smod", dcols),
                  "Ruralness, population density, built environment",
                  "Household assets and characteristics"),
           TRUE,
           "Earth Engine zonal sum at 1 km over the GADM Admin-2 polygon, expressed as a SHARE of total population so it is comparable across countries. TEMPORAL CAVEAT: the WorldPop age-sex collection carries a single 2020 vintage, so this is 2020 composition matched to surveys from 2013-2021 - a lag of up to 7 years for Sierra Leone. GHS-SMOD is the degree-of-urbanisation class, 10 rural to 30 urban centre.")
  cat(sprintf("  -> %d demography columns
", length(dcols)))
} else cat("  demography CSV absent; run 09_extract_gee_demography.py first
")


# ── 5. Zone stratifiers (Koppen-Geiger, AEZ16) ─────────────────────────────
# Built by scripts/protocol_v2/10_build_zone_stratifiers.R. The 41 climate and
# 93 agriculture columns are all CONTINUOUS surfaces; a zone class is the
# stratifier this literature actually uses, is stable rather than year-specific,
# and encodes interactions a linear model over separate continuous columns
# cannot reach. Classes are emitted one-hot so no learner reads class 12 as
# greater than class 4.
cat("
[zone stratifiers]
")
zs <- tryCatch(read.csv("data/covariates/harmonized/zone_stratifiers_admin2.csv",
                        stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(zs)) {
  # the raw class CODE is an unordered category: keep purity, heterogeneity and
  # the one-hot indicators, drop the code itself so nothing treats it as ordinal
  zs <- zs[, !grepl("_class$", names(zs)), drop = FALSE]
  blocks$zones <- zs
  zcols <- setdiff(names(zs), c("country", "Admin1", "Admin2"))
  add_meta(zcols, "Koppen-Geiger 1991-2020 / IFPRI-HarvestChoice AEZ16",
           "Climate and weather", TRUE,
           "Modal zone class over the Admin-2 polygon. *_purity is the share of the district in its modal class (1 = squarely inside one zone); *_n_classes counts distinct zones present; *_is_N are one-hot indicators for classes common enough to be usable. The raw ordinal class code is deliberately NOT carried, because the classes are unordered.")
  cat(sprintf("  -> %d zone columns
", length(zcols)))
} else cat("  zone CSV absent; run 10_build_zone_stratifiers.R first
")


# ── 6. Relative Wealth Index, year-matched density, GPW 2010 age-sex ────────
# Built by scripts/protocol_v2/11_extract_gee_rwi_density.py.
# RWI is the highest-resolution SES surface publicly available for these
# countries (2.4 km, validated against DHS wealth), against a current SES block
# of nine prior-round DHS aggregates. Density is taken at each country's OWN
# survey year, which supersedes the 2020-only density from step 09. GPW 2010
# age-sex is carried as the second temporal bracket: checked against the
# catalogue, WorldPop age-sex exists only for 2020 and GPW only for 2010, there
# is no annual age-sex product, and those two vintages bracket the 2013-2021
# survey window.
cat("
[RWI / density / GPW]
")
rw <- tryCatch(read.csv("data/covariates/harmonized/gee_rwi_density_admin2.csv",
                        stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(rw)) {
  rw$wpop_density_year <- NULL   # provenance, not a predictor
  blocks$rwi <- rw
  rcols <- setdiff(names(rw), c("country", "Admin1", "Admin2"))
  add_meta(rcols,
           "Meta/Data for Good RWI; WorldPop; CIESIN GPW v4.11 (Earth Engine)",
           ifelse(grepl("^gpw2010", rcols), "Household assets and characteristics",
           ifelse(grepl("^rwi", rcols), "Education, employment, SES",
                  "Ruralness, population density, built environment")),
           TRUE,
           ifelse(grepl("^rwi", rcols),
             "Relative Wealth Index, mean/SD/count of Meta 2.4 km points falling inside the Admin-2 polygon. rwi_n_points says how many points the district estimate rests on, so a district resting on three is distinguishable from one resting on three hundred.",
           ifelse(grepl("^gpw2010", rcols),
             "GPW v4.11 2010 age-sex share. SECOND TEMPORAL BRACKET: WorldPop age-sex exists only for 2020 and GPW only for 2010; there is no annual age-sex product, and these two vintages bracket the 2013-2021 survey window. Composition moves slowly, so the bracket is an honest representation of what is knowable.",
             "WorldPop population per km2 (log1p) at each country's OWN survey year (Gambia 2018, Ghana 2017, Malawi 2015, Sierra Leone 2013) - the time-matched replacement for the 2020-only density. Until 2026-09-07 this was the log of the polygon COUNT (AU-01 finding 3).")))
  cat(sprintf("  -> %d RWI/density/GPW columns
", length(rcols)))
} else cat("  RWI CSV absent; run 11_extract_gee_rwi_density.py first
")

# ── 7. Livestock density (GLW4 2015; script 48) ─────────────────────────────
writeLines("\n[GLW4 livestock]")
lf <- file.path(HDIR, "predictors_admin2_livestock.csv")
if (file.exists(lf)) {
  LV <- read.csv(lf, check.names = FALSE); blocks$livestock <- LV
  lcols <- setdiff(names(LV), c("country", "Admin1", "Admin2"))
  add_meta(lcols, "Gridded Livestock of the World 4 (2020, 5 arc-min, dasymetric; FAO catalog)",
           "Livestock density", TRUE,
           "Head per km2 (cattle, sheep, goats, pigs, chickens), tropical livestock units per km2 (0.7/0.1/0.1/0.2/0.01) and per person (WorldPop), ruminant share of TLU; area-weighted zonal mean of the GLW4-2020 dasymetric density surface (script 48; the 2015 Dataverse release ranks districts identically, Spearman 0.99-1.00). Proxy for animal-source food availability (iron, zinc, B12, preformed vitamin A).")
  writeLines(sprintf("  -> %d livestock columns", length(lcols)))
} else writeLines("  livestock CSV absent; run scripts/protocol_v2/48_ihme_raster_and_livestock.R")

# ── 8. Distance to surface water and to the coast (Earth Engine) ─────────────
writeLines("\n[water and coast distance]")
wf <- file.path(HDIR, "predictors_admin2_water_distance.csv")
if (file.exists(wf)) {
  WD <- read.csv(wf, check.names = FALSE); blocks$wdist <- WD
  wcols <- setdiff(names(WD), c("country", "Admin1", "Admin2"))
  add_meta(wcols, "JRC Global Surface Water 1.4; USDOS LSIB 2017 (Earth Engine)",
           "Water and coast proximity", TRUE,
           "Distance (km) from each pixel to the nearest permanent (occurrence >= 50%) or any (>= 10%) surface water and to the ocean, distance transform at 500 m / 2 km in Web Mercator corrected by cos(latitude); polygon mean and minimum (scripts/covariates/gee_water_coast_distance.py). Fish and aquatic-food access; the iodine geography (inland, elevated soils are iodine-poor).")
  writeLines(sprintf("  -> %d water/coast distance columns", length(wcols)))
} else writeLines("  water-distance CSV absent; run scripts/covariates/gee_water_coast_distance.py")

# ── 9. Helminth burden and control (ESPEN; script 49) ────────────────────────
writeLines("\n[ESPEN helminths]")
ef <- file.path(HDIR, "predictors_admin2_espen.csv")
if (file.exists(ef)) {
  ES <- read.csv(ef, check.names = FALSE); blocks$espen <- ES
  ecols <- setdiff(names(ES), c("country", "Admin1", "Admin2"))
  add_meta(ecols, "WHO ESPEN implementation-unit database 2014-2025 (portal export, no key)",
           "Helminth burden and control", TRUE,
           "Programme endemicity class midpoint (%) for soil-transmitted helminths and schistosomiasis at the year nearest the survey and at the earliest reported year; share of 2014-2018 with mass drug administration delivered; mean reported epidemiological coverage over delivered years. IU = ESPEN ADM2, name-matched to the spine with aliases for post-split districts; Malawi at district, broadcast to Traditional Authorities (script 49).")
  writeLines(sprintf("  -> %d helminth columns", length(ecols)))
} else writeLines("  ESPEN CSV absent; run scripts/protocol_v2/49_espen_admin2_block.R")

# ── 10. Malaria Atlas at the survey year (script 52) ─────────────────────────
# Replaces the release-year map_malaria* / map_interventions* columns (AU-01
# finding 4: rasters requested at the dataset's release date, 4-12 years after
# the surveys, three releases of one product kept as separate columns). The
# static blood-disorder surfaces (2012 release) are untouched.
writeLines("\n[Malaria Atlas, survey year]")
mf <- file.path(HDIR, "predictors_admin2_map_sy.csv")
MAP_DROP <- character(0)
if (file.exists(mf)) {
  MS <- read.csv(mf, check.names = FALSE); blocks$map_sy <- MS
  mcols <- setdiff(names(MS), c("country", "Admin1", "Admin2"))
  MAP_DROP <- grep("^map_(malaria|interventions)", names(SH), value = TRUE)
  add_meta(mcols, "Malaria Atlas Project (latest release, layer at each survey year)",
           "Malaria incidence and treatment", TRUE,
           "Pf parasite rate, incidence and mortality rates, reproductive number, ITN access / use / use rate, IRS coverage and effective-treatment coverage from the latest MAP release, the annual layer at each country's survey year (Gambia 2018, Ghana 2017, Malawi 2015, Sierra Leone 2013), clipped to the country and averaged over the spine polygon (area-weighted, as the block it replaces). Script 52.")
  writeLines(sprintf("  -> %d survey-year columns; %d release-year map_ columns dropped", length(mcols), length(MAP_DROP)))
} else writeLines("  survey-year MAP CSV absent; run scripts/protocol_v2/52_map_survey_year.R")

# ── assemble and append, without the four-country filter ────────────────────
if (!length(blocks)) stop("no blocks built")
EX <- spine
for (b in names(blocks)) {
  by <- intersect(c("country", "Admin1", "Admin2"), names(blocks[[b]]))
  j <- tryCatch(left_join(EX, blocks[[b]], by = by), error = function(e) NULL)
  if (is.null(j) || nrow(j) != nrow(EX)) {
    cat("  [join]", b, "changed row count or failed, skipped\n"); next
  }
  EX <- j
}
newcols <- setdiff(names(EX), c("country", "Admin1", "Admin2"))
# a column with no finite value anywhere carries nothing and would only show up
# later as an "empty" flag in the variable sheet
empty <- newcols[vapply(newcols, function(v) !any(is.finite(EX[[v]])), TRUE)]
if (length(empty)) {
  cat("  dropping", length(empty), "all-missing columns:",
      paste(empty, collapse = ", "), "
")
  EX <- EX[, setdiff(names(EX), empty), drop = FALSE]
  newcols <- setdiff(newcols, empty)
}
write.csv(EX, file.path(HDIR, "predictors_admin2_extrasrc.csv"), row.names = FALSE)

MD <- bind_rows(meta) |> filter(column %in% newcols) |> distinct(column, .keep_all = TRUE)
MD$n_countries <- vapply(MD$column, function(v)
  sum(tapply(EX[[v]], EX$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
MD$countries <- vapply(MD$column, function(v) {
  s <- tapply(EX[[v]], EX$country, function(z) any(is.finite(z)))
  paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
MD$completeness <- round(vapply(MD$column, function(v) mean(is.finite(EX[[v]])), 0), 3)
MD$coverage_by_country <- vapply(MD$column, function(v) {
  s <- EX |> group_by(country) |>
    summarise(ok = mean(is.finite(.data[[v]])), .groups = "drop")
  paste(sprintf("%s=%.2f", s$country, s$ok), collapse = ";") }, "")
write.csv(MD, file.path(HDIR, "predictors_admin2_extrasrc_metadata.csv"),
          row.names = FALSE)

SHM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"),
                stringsAsFactors = FALSE)
# the Python builder writes True/False, which read.csv keeps as character; bind_rows below needs logical
SHM$subnational <- as.logical(toupper(as.character(SHM$subnational)))
for (f in c("predictors_admin2_shared.csv", "predictors_admin2_shared_metadata.csv"))
  if (!file.exists(file.path(HDIR, paste0(f, ".pre_extrasrc"))))
    file.copy(file.path(HDIR, f), file.path(HDIR, paste0(f, ".pre_extrasrc")))

SH2 <- SH |> select(-any_of(c(newcols, MAP_DROP))) |>
  left_join(EX, by = c("country", "Admin1", "Admin2"))
SHM2 <- bind_rows(SHM |> filter(!column %in% c(newcols, MAP_DROP)),
                  MD |> transmute(column, domain, source, n_countries, countries,
                                  completeness, subnational))
write.csv(SH2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(SHM2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"),
          row.names = FALSE)

cat("\n=== added ===\n")
print(as.data.frame(MD[, c("column", "source", "n_countries", "completeness",
                           "subnational")]), row.names = FALSE)
cat("\nshared set:", ncol(SH) - 3, "->", ncol(SH2) - 3, "predictors\n")
cat("NEAR-OUTCOME flagged:", sum(grepl("NEAR-OUTCOME", MD$domain)), "\n")
cat("\nDONE\n")
