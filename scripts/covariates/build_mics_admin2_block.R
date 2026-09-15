# =============================================================================
# scripts/covariates/build_mics_admin2_block.R   [MC-01, 2026-09-15]
#
# MICS MICRODATA -> ADMIN-2 PREDICTORS FOR THE FOUR COUNTRIES
#
# The one household-survey programme besides DHS that measures, on the same
# households, what the biomarkers respond to: the household salt iodine test,
# infant and young child feeding, WASH, wealth, child anthropometry and recent
# illness, women's education, malaria prophylaxis in pregnancy. Four rounds,
# all public microdata (data/MICS/, UNICEF MICS):
#
#   The Gambia    MICS6 2018   390 clusters, GPS in hand -> GADM district by
#                              point-in-polygon; LGA (HH7) mean where a district
#                              has fewer than MIN_CL clusters. THE GMNS 2018 WAS
#                              FIELDED INSIDE THIS SAMPLE (70 of its clusters carry
#                              the biomarker respondents), so those clusters are
#                              dropped first (LK-02; metadata/mns_dhs_overlap_clusters.csv,
#                              programme MICS).
#   Ghana         MICS6 2017-18  10 regions (HH7) -> the 16 GADM regions through
#                              the 16-to-10 crosswalk -> every district (no GPS yet)
#   Malawi        MICS5 2013-14  31 districts and cities (HH7) -> the district's
#                              TAs; the four city estimates go to the city TAs
#                              (GADM ENGTYPE City / Urban) and the district
#                              estimate to the rest (no GPS released for this round)
#   Sierra Leone  MICS6 2017   14 districts (HH7A) = the spine's Admin-2
#
# Indicators (weighted means; hhweight / chweight / wmweight; a unit needs
# MIN_N respondents or it takes its parent's value, recorded in mics_level):
#   household  salt_iodised_15ppm, salt_any_iodine (among tested), water_improved,
#              water_piped, sanitation_improved, open_defecation, handwash_soap,
#              electricity (MICS6 only), wealth_poorest40, wealth_score_mean,
#              head_no_education
#   children   c_stunted, c_wasted, c_underweight, c_mean_haz (WHO z-scores,
#              flags clean), c_diarrhoea_2wk, c_fever_2wk, c_cough_2wk (MICS6),
#              c_vas_6mo (MICS5 only), c_mdd (6-23 months, >= 5 of 8 groups incl.
#              breastmilk, WHO 2021), c_fg_dairy, c_fg_flesh, c_fg_eggs,
#              c_fg_legumes, c_fg_vita_fv, c_fg_n (mean groups)
#   women      w_no_education, w_secondary_plus, w_literate, w_iptp_sp (SP in
#              last pregnancy; Malawi: any malaria prophylaxis)
# Code lists follow the JMP and WHO definitions and are written out as
# metadata/mics_codes.csv for review.
#
#   Rscript -e "source('scripts/covariates/build_mics_admin2_block.R')"
# -> data/covariates/harmonized/predictors_admin2_mics.csv (+ _metadata.csv)
# -> data/covariates/harmonized/mics_unit_estimates.csv   (per country x unit, with n)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(haven); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/admin2_key_hygiene.R"); source("R/admin2_keys.R"); source("R/mns_dhs_overlap.R")
MDIR <- "data/MICS"; HDIR <- "data/covariates/harmonized"
MIN_N <- 30L; MIN_CL <- 5L
num <- function(x) suppressWarnings(as.numeric(unclass(x)))
wm  <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (sum(ok) == 0) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
yes <- function(x) ifelse(is.na(num(x)) | num(x) %in% c(8, 9, 98, 99), NA, as.numeric(num(x) == 1))
SP <- admin2_spine()
spine <- SP[, c("country", "Admin1", "Admin2", "engtype_2")]

SURV <- list(
  Gambia      = list(dir = "Gambia 2018",                      round = 6L, year = 2018L),
  Ghana       = list(dir = "Ghana 2017",                       round = 6L, year = 2017L),
  Malawi      = list(dir = "Malawi MICS 2013-14 SPSS Datasets", round = 5L, year = 2014L),
  SierraLeone = list(dir = "Sierra Leone MICS6 Datasets",      round = 6L, year = 2017L))

# ── code lists (JMP / WHO) ───────────────────────────────────────────────────
WATER_IMPROVED <- c(11, 12, 13, 14, 21, 31, 41, 51, 61, 71, 72, 91, 92)
WATER_PIPED    <- c(11, 12, 13, 14)
SAN_IMPROVED   <- c(11, 12, 13, 15, 18, 21, 22, 24, 31)   # flush to sewer/septic/pit/unknown, VIP, pit with slab, composting
SAN_OPEN       <- 95
codes <- data.frame(item = c("water_improved", "water_piped", "sanitation_improved", "open_defecation", "salt_tested", "salt_15ppm", "literate", "iycf_groups"),
                    codes = c(paste(WATER_IMPROVED, collapse = ","), paste(WATER_PIPED, collapse = ","), paste(SAN_IMPROVED, collapse = ","), "95", "SA1/SI1 in 1,2,3", "SA1/SI1 == 3", "WB14/WB7 == 3 of 1,2,3",
                              "grains BD8B|BD8C|BD8E; legumes BD8M; dairy BD7D|BD7E|BD8A|BD8N; flesh BD8I|BD8J|BD8L; eggs BD8K; vitA BD8D|BD8F|BD8G; other BD8H; breastmilk BD3"),
                    note = c("JMP improved incl. delivered (tanker, cart, kiosk) and packaged water", "piped into dwelling/yard/neighbour/public tap", "JMP improved; Ghana 14 (open drain) and Malawi 14 (elsewhere) unimproved", "no facility / bush / field", "tested households only", "at least 15 ppm", "reads the whole sentence", "WHO 2021 MDD: >= 5 of 8 groups, children 6-23 months"))
dir.create("metadata", showWarnings = FALSE); write.csv(codes, "metadata/mics_codes.csv", row.names = FALSE)

# ── per-survey household / child / woman tables ─────────────────────────────
read_mics <- function(cn) {
  s <- SURV[[cn]]; d <- file.path(MDIR, s$dir); r6 <- s$round == 6L
  hh <- read_sav(file.path(d, "hh.sav")) |> as.data.frame(); ch <- read_sav(file.path(d, "ch.sav")) |> as.data.frame(); wm_ <- read_sav(file.path(d, "wm.sav")) |> as.data.frame()
  hh$HH1 <- as.integer(num(hh$HH1)); ch$HH1 <- as.integer(num(ch$HH1)); wm_$HH1 <- as.integer(num(wm_$HH1))
  # geography labels carried on every file
  geo <- function(x) { g <- data.frame(HH1 = x$HH1, stringsAsFactors = FALSE)
    g$adm1 <- as.character(haven::as_factor(x$HH7)); if ("HH7A" %in% names(x)) g$adm2 <- as.character(haven::as_factor(x$HH7A)); g }
  # household
  salt <- num(if (r6) hh$SA1 else hh$SI1); tested <- salt %in% 1:3
  toilet <- num(if (r6) hh$WS11 else hh$WS8); water <- num(hh$WS1)
  H <- cbind(geo(hh), data.frame(
    w = num(hh$hhweight),
    salt_iodised_15ppm = ifelse(tested, as.numeric(salt == 3), NA), salt_any_iodine = ifelse(tested, as.numeric(salt %in% 2:3), NA),
    water_improved = ifelse(is.na(water) | water == 99, NA, as.numeric(water %in% WATER_IMPROVED)), water_piped = ifelse(is.na(water) | water == 99, NA, as.numeric(water %in% WATER_PIPED)),
    sanitation_improved = ifelse(is.na(toilet) | toilet == 99, NA, as.numeric(toilet %in% SAN_IMPROVED)), open_defecation = ifelse(is.na(toilet) | toilet == 99, NA, as.numeric(toilet == SAN_OPEN)),
    handwash_soap = if (r6) as.numeric(num(hh$HW1) %in% 1:3 & num(hh$HW3) == 1) else as.numeric(num(hh$HW1) == 1 & num(hh$HW3A) == 1),
    electricity = if ("HC8" %in% names(hh)) { v <- num(hh$HC8); ifelse(is.na(v) | v == 9, NA, as.numeric(v %in% c(1, 2))) } else NA_real_,
    wealth_poorest40 = as.numeric(num(hh$windex5) <= 2), wealth_score_mean = num(hh$wscore),
    head_no_education = ifelse(num(hh$helevel) == 9, NA, as.numeric(num(hh$helevel) == if (r6) 0 else 1)), stringsAsFactors = FALSE))
  # children under 5
  zs <- function(z, fl) { z <- num(z); f <- num(fl); z[!is.finite(z) | z > 90 | (is.finite(f) & f != 0)] <- NA; z }
  haz <- zs(ch$HAZ2, ch$HAZFLAG); whz <- zs(ch$WHZ2, ch$WHZFLAG); waz <- zs(ch$WAZ2, ch$WAZFLAG)
  age <- num(ch$CAGE); iy <- is.finite(age) & age >= 6 & age <= 23
  g1 <- function(...) { m <- do.call(cbind, lapply(list(...), function(v) if (v %in% names(ch)) as.numeric(num(ch[[v]]) == 1) else rep(0, nrow(ch)))); as.numeric(rowSums(m, na.rm = TRUE) > 0) }
  fg <- data.frame(grains = g1("BD8B", "BD8C", "BD8E"), legumes = g1("BD8M"), dairy = g1("BD7D", "BD7E", "BD8A", "BD8N"), flesh = g1("BD8I", "BD8J", "BD8L"),
                   eggs = g1("BD8K"), vita = g1("BD8D", "BD8F", "BD8G"), other = g1("BD8H"), bm = as.numeric(num(ch$BD3) == 1))
  fg[!iy, ] <- NA
  ngrp <- rowSums(fg)
  vas <- if (r6) rep(NA_real_, nrow(ch)) else { v <- num(ch$IM3VD); ifelse(is.na(v) | v %in% c(97, 98, 99), NA, as.numeric(v %in% c(1:31, 44, 66))) }
  C <- cbind(geo(ch), data.frame(
    w = num(ch$chweight),
    c_stunted = as.numeric(haz < -2), c_wasted = as.numeric(whz < -2), c_underweight = as.numeric(waz < -2), c_mean_haz = haz,
    c_diarrhoea_2wk = yes(ch$CA1), c_fever_2wk = yes(if (r6) ch$CA14 else ch$CA6AA), c_cough_2wk = if (r6) yes(ch$CA16) else NA_real_,
    c_vas_6mo = vas, c_mdd = as.numeric(ngrp >= 5), c_fg_dairy = fg$dairy, c_fg_flesh = fg$flesh, c_fg_eggs = fg$eggs, c_fg_legumes = fg$legumes,
    c_fg_vita_fv = fg$vita, c_fg_n = ngrp, stringsAsFactors = FALSE))
  # women 15-49
  wl <- num(wm_$welevel); wl[wl == 9] <- NA
  lit <- num(if (r6) wm_$WB14 else wm_$WB7)
  W <- cbind(geo(wm_), data.frame(
    w = num(wm_$wmweight),
    w_no_education = as.numeric(wl == if (r6) 0 else 1),
    w_secondary_plus = as.numeric(wl >= if (cn == "Malawi") 3 else 2),
    w_literate = ifelse(lit %in% 1:3, as.numeric(lit == 3), NA),
    w_iptp_sp = if (cn == "Malawi") yes(wm_$MN13) else yes(wm_$MN16), stringsAsFactors = FALSE))
  list(H = H, C = C, W = W)
}

agg_unit <- function(D, keys) {   # weighted means per unit for every indicator column, with n
  ind <- setdiff(names(D), c("HH1", "adm1", "adm2", "w", "unit", keys))
  D |> group_by(across(all_of(keys))) |> summarise(n = n(), n_cl = n_distinct(HH1), across(all_of(ind), ~ wm(.x, w)), .groups = "drop") |> as.data.frame()
}

out_rows <- list(); unit_rows <- list(); IND <- NULL
for (cn in names(SURV)) {
  cat(sprintf("\n[%s] MICS%d %d\n", cn, SURV[[cn]]$round, SURV[[cn]]$year))
  M <- read_mics(cn)
  if (cn == "Gambia") for (k in names(M)) M[[k]] <- drop_mns_overlap(M[[k]], "Gambia", "HH1", programme = "MICS")
  ind_h <- setdiff(names(M$H), c("HH1", "adm1", "adm2", "w")); ind_c <- setdiff(names(M$C), c("HH1", "adm1", "adm2", "w")); ind_w <- setdiff(names(M$W), c("HH1", "adm1", "adm2", "w"))
  IND <- c(ind_h, ind_c, ind_w)
  sp <- spine[spine$country == cn, ]

  if (cn == "Gambia") {
    # clusters -> GADM district by point-in-polygon on the displaced GPS
    gps <- sf::st_read("data/MICS/GPS/GMB_2018_MICS_v01_M/GPS Datasets/GambiaMICS2018GPS/GambiaMICS2018GPS.shp", quiet = TRUE) |> sf::st_transform(4326)
    poly <- readRDS("data/admin_boundaries/gadm41_GMB_2.rds") |> sf::st_transform(4326)
    hit <- sf::st_intersects(gps, poly); idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, 1L)
    cl <- data.frame(HH1 = as.integer(gps$HH1), Admin1 = poly$NAME_1[idx], Admin2 = poly$NAME_2[idx], stringsAsFactors = FALSE)
    miss <- cl[is.na(cl$Admin2), ]; if (nrow(miss)) { nn <- sf::st_nearest_feature(gps[is.na(idx), ], poly); cl$Admin1[is.na(cl$Admin2)] <- poly$NAME_1[nn]; cl$Admin2[is.na(cl$Admin2)] <- poly$NAME_2[nn] }
    cat(sprintf("  %d GPS clusters -> districts (%d by nearest polygon); clusters after dropping the GMNS overlap: %d\n", nrow(cl), nrow(miss), length(unique(M$H$HH1))))
    for (k in names(M)) M[[k]] <- M[[k]] |> left_join(cl, by = "HH1")
    A2 <- lapply(M, agg_unit, keys = c("Admin1", "Admin2")); A1 <- lapply(M, agg_unit, keys = "Admin1")
    U <- Reduce(function(x, y) full_join(x, y, by = c("Admin1", "Admin2"), suffix = c("", ".y")), lapply(A2, function(a) a[, setdiff(names(a), c("n", "n_cl"))]))
    # the household table decides the unit's level (n households and clusters)
    U$n <- A2$H$n[match(paste(U$Admin1, U$Admin2), paste(A2$H$Admin1, A2$H$Admin2))]; U$n_cl <- A2$H$n_cl[match(paste(U$Admin1, U$Admin2), paste(A2$H$Admin1, A2$H$Admin2))]
    P <- Reduce(function(x, y) full_join(x, y, by = "Admin1", suffix = c("", ".y")), lapply(A1, function(a) a[, setdiff(names(a), c("n", "n_cl"))]))
    U <- sp[, c("Admin1", "Admin2")] |> left_join(U[, c("Admin1", "Admin2", "n", "n_cl", IND)], by = c("Admin1", "Admin2"))
    U$n[is.na(U$n)] <- 0L; U$n_cl[is.na(U$n_cl)] <- 0L
    U$mics_level <- ifelse(U$n >= MIN_N & U$n_cl >= MIN_CL, "admin2", "admin1_broadcast")
    low <- which(U$mics_level != "admin2"); j <- match(U$Admin1[low], P$Admin1); for (v in IND) U[[v]][low] <- P[[v]][j]
    cat(sprintf("  %d districts from their own clusters, %d from the LGA mean\n", sum(U$mics_level == "admin2"), length(low)))
  } else if (cn == "SierraLeone") {
    for (k in names(M)) { M[[k]]$Admin2 <- admin2_match_v2(M[[k]]$adm2, sp$Admin2, review_csv = if (k == "H") "metadata/crosswalks/review_mics_SierraLeone_district.csv" else NULL, label = "SL MICS district")
      M[[k]]$Admin1 <- sp$Admin1[match(M[[k]]$Admin2, sp$Admin2)] }   # pair key throughout (JK-01)
    A2 <- lapply(M, agg_unit, keys = c("Admin1", "Admin2"))
    U <- Reduce(function(x, y) full_join(x, y, by = c("Admin1", "Admin2"), suffix = c("", ".y")), lapply(A2, function(a) a[, setdiff(names(a), c("n", "n_cl"))]))
    U$n <- A2$H$n[match(paste(U$Admin1, U$Admin2), paste(A2$H$Admin1, A2$H$Admin2))]
    U <- sp[, c("Admin1", "Admin2")] |> left_join(U[, c("Admin1", "Admin2", "n", IND)], by = c("Admin1", "Admin2")); U$mics_level <- "admin2"
    cat(sprintf("  %d of %d districts matched\n", sum(!is.na(U$n)), nrow(sp)))
  } else if (cn == "Malawi") {
    # HH7 = 27 districts + 4 cities. District -> Admin1; city -> the city TAs (ENGTYPE City / Urban) of that district
    lab <- unique(M$H$adm1); is_city <- grepl("city$", lab, ignore.case = TRUE)
    base <- sub("\\s*city$", "", lab, ignore.case = TRUE)
    a1 <- admin2_match_v2(base, unique(sp$Admin1), aliases_csv = "metadata/crosswalks/aliases_mics_malawi.csv", review_csv = "metadata/crosswalks/review_mics_Malawi_district.csv", label = "Malawi MICS district")
    map <- data.frame(adm1 = lab, Admin1 = a1, city = is_city, stringsAsFactors = FALSE)
    for (k in names(M)) M[[k]] <- M[[k]] |> left_join(map, by = "adm1")
    A <- lapply(M, agg_unit, keys = c("Admin1", "city"))
    U0 <- Reduce(function(x, y) full_join(x, y, by = c("Admin1", "city"), suffix = c("", ".y")), lapply(A, function(a) a[, setdiff(names(a), c("n", "n_cl"))]))
    U0$n <- A$H$n[match(paste(U0$Admin1, U0$city), paste(A$H$Admin1, A$H$city))]
    sp$city <- sp$engtype_2 %in% c("City", "Urban") & sp$Admin1 %in% U0$Admin1[U0$city]
    U <- sp[, c("Admin1", "Admin2", "city")] |> left_join(U0[, c("Admin1", "city", "n", IND)], by = c("Admin1", "city"))
    # a city TA in a district without a separate city estimate takes the district estimate
    miss <- which(is.na(U$n)); if (length(miss)) { j <- match(paste(U$Admin1[miss], FALSE), paste(U0$Admin1, U0$city)); for (v in c("n", IND)) U[[v]][miss] <- U0[[v]][j] }
    U$mics_level <- ifelse(U$city, "admin1_city", "admin1_broadcast"); U$city <- NULL
    cat(sprintf("  %d MICS districts/cities -> %d TAs (%d city TAs with their own city estimate)\n", length(lab), nrow(U), sum(U$mics_level == "admin1_city")))
  } else if (cn == "Ghana") {
    XW <- read.csv("data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv", stringsAsFactors = FALSE)
    r10 <- admin2_match_v2(unique(M$H$adm1), unique(XW$adm1_paper_10), review_csv = "metadata/crosswalks/review_mics_Ghana_region.csv", label = "Ghana MICS region")
    map <- data.frame(adm1 = unique(M$H$adm1), r10 = r10, stringsAsFactors = FALSE)
    for (k in names(M)) M[[k]] <- M[[k]] |> left_join(map, by = "adm1")
    A <- lapply(M, agg_unit, keys = "r10")
    U0 <- Reduce(function(x, y) full_join(x, y, by = "r10", suffix = c("", ".y")), lapply(A, function(a) a[, setdiff(names(a), c("n", "n_cl"))]))
    U0$n <- A$H$n[match(U0$r10, A$H$r10)]
    sp$r10 <- XW$adm1_paper_10[match(admin2_kk(sp$Admin1), admin2_kk(XW$admin1_16))]; sp$r10[is.na(sp$r10)] <- sp$Admin1[is.na(sp$r10)]
    U <- sp[, c("Admin1", "Admin2", "r10")] |> left_join(U0[, c("r10", "n", IND)], by = "r10"); U$r10 <- NULL; U$mics_level <- "admin1_broadcast"
    cat(sprintf("  10 regions -> %d districts (%d matched)\n", nrow(U), sum(!is.na(U$n))))
  }
  U$country <- cn
  unit_rows[[cn]] <- U
  cat(sprintf("  weighted means: salt >= 15 ppm %.2f | improved water %.2f | improved sanitation %.2f | stunted %.2f | MDD %.2f | women secondary+ %.2f\n",
              mean(U$salt_iodised_15ppm, na.rm = TRUE), mean(U$water_improved, na.rm = TRUE), mean(U$sanitation_improved, na.rm = TRUE), mean(U$c_stunted, na.rm = TRUE), mean(U$c_mdd, na.rm = TRUE), mean(U$w_secondary_plus, na.rm = TRUE)))
}
ALL <- bind_rows(unit_rows)
OUT <- spine[, c("country", "Admin1", "Admin2")] |> left_join(ALL[, c("country", "Admin1", "Admin2", "n", "mics_level", IND)], by = c("country", "Admin1", "Admin2"))
stopifnot(nrow(OUT) == nrow(spine))
names(OUT)[names(OUT) %in% IND] <- paste0("mics_", IND); names(OUT)[names(OUT) == "n"] <- "mics_n_hh"
write.csv(OUT, file.path(HDIR, "predictors_admin2_mics.csv"), row.names = FALSE)
write.csv(ALL, file.path(HDIR, "mics_unit_estimates.csv"), row.names = FALSE)

# ── metadata ─────────────────────────────────────────────────────────────────
dom <- c(salt_iodised_15ppm = "Food fortification and supplementation", salt_any_iodine = "Food fortification and supplementation",
         water_improved = "Water and sanitation", water_piped = "Water and sanitation", sanitation_improved = "Water and sanitation", open_defecation = "Water and sanitation", handwash_soap = "Water and sanitation",
         electricity = "Household assets and characteristics", wealth_poorest40 = "Household assets and characteristics", wealth_score_mean = "Household assets and characteristics", head_no_education = "Education, employment, SES",
         c_stunted = "Child anthropometry", c_wasted = "Child anthropometry", c_underweight = "Child anthropometry", c_mean_haz = "Child anthropometry",
         c_diarrhoea_2wk = "Infection and inflammation burden", c_fever_2wk = "Infection and inflammation burden", c_cough_2wk = "Infection and inflammation burden",
         c_vas_6mo = "Food fortification and supplementation", c_mdd = "Infant and young child feeding", c_fg_dairy = "Infant and young child feeding", c_fg_flesh = "Infant and young child feeding",
         c_fg_eggs = "Infant and young child feeding", c_fg_legumes = "Infant and young child feeding", c_fg_vita_fv = "Infant and young child feeding", c_fg_n = "Infant and young child feeding",
         w_no_education = "Education, employment, SES", w_secondary_plus = "Education, employment, SES", w_literate = "Education, employment, SES", w_iptp_sp = "Infection and inflammation burden")
desc <- c(salt_iodised_15ppm = "Share of households with tested salt at >= 15 ppm iodine (rapid test kit; tested households)", salt_any_iodine = "Share of tested households with any iodine in salt (> 0 ppm)",
          water_improved = "Share of households using an improved drinking-water source (JMP, incl. delivered and packaged water)", water_piped = "Share of households with piped drinking water (dwelling, yard, neighbour or public tap)",
          sanitation_improved = "Share of households using an improved sanitation facility (JMP)", open_defecation = "Share of households with no facility (open defecation)",
          handwash_soap = "Share of households with an observed handwashing place with soap", electricity = "Share of households with electricity (grid or off-grid; MICS6 only)",
          wealth_poorest40 = "Share of households in the two poorest national wealth quintiles", wealth_score_mean = "Mean MICS wealth score (national standardisation, within-country ranking only)",
          head_no_education = "Share of household heads with no or pre-primary education", c_stunted = "Share of children under 5 with height-for-age z < -2 (WHO, flagged values excluded)",
          c_wasted = "Share of children under 5 with weight-for-height z < -2", c_underweight = "Share of children under 5 with weight-for-age z < -2", c_mean_haz = "Mean height-for-age z-score, children under 5",
          c_diarrhoea_2wk = "Share of children under 5 with diarrhoea in the last two weeks", c_fever_2wk = "Share of children under 5 with fever in the last two weeks", c_cough_2wk = "Share of children under 5 with cough in the last two weeks (MICS6 only)",
          c_vas_6mo = "Share of children 6-59 months with a vitamin A dose recorded (MICS5 only; MICS6 dropped the item)", c_mdd = "Share of children 6-23 months meeting minimum dietary diversity (>= 5 of 8 food groups incl. breastmilk, WHO 2021)",
          c_fg_dairy = "Share of children 6-23 months given milk, formula, yogurt or cheese yesterday", c_fg_flesh = "Share of children 6-23 months given meat, organ meat or fish yesterday", c_fg_eggs = "Share of children 6-23 months given eggs yesterday",
          c_fg_legumes = "Share of children 6-23 months given beans, lentils or nuts yesterday", c_fg_vita_fv = "Share of children 6-23 months given vitamin-A-rich fruit or vegetables yesterday", c_fg_n = "Mean number of the 8 food groups consumed yesterday, children 6-23 months",
          w_no_education = "Share of women 15-49 with no or pre-primary education", w_secondary_plus = "Share of women 15-49 with secondary education or higher", w_literate = "Share of women 15-49 who read a whole sentence",
          w_iptp_sp = "Share of women with a live birth in the last two years who took SP/Fansidar (Malawi MICS5: any malaria prophylaxis) during the pregnancy")
basis <- "MICS microdata, weighted (hhweight / chweight / wmweight). The Gambia MICS6 2018: clusters to GADM districts by GPS point-in-polygon after dropping the 70 clusters in which the GMNS 2018 was fielded (LK-02), LGA mean where a district has < 30 households or < 5 clusters; Ghana MICS6 2017-18: 10 regions broadcast through the 16-to-10 crosswalk; Malawi MICS5 2013-14: 27 districts + 4 cities, city estimates to the city TAs, district estimates to the other TAs; Sierra Leone MICS6 2017: 14 districts. Salt, water, sanitation and literacy code lists in metadata/mics_codes.csv."
cols <- paste0("mics_", IND)
md <- data.frame(column = c(cols, "mics_n_hh"), source = "MICS microdata (UNICEF MICS6 2017/2018, MICS5 2013-14)", domain = c(unname(dom[IND]), "Household assets and characteristics"), subnational = TRUE,
                 assumption = paste0(c(unname(desc[IND]), "Households behind the unit estimate (audit column)"), ". ", basis), stringsAsFactors = FALSE)
md$n_countries <- vapply(md$column, function(v) sum(tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
md$countries <- vapply(md$column, function(v) { s <- tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))); paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
md$completeness <- round(vapply(md$column, function(v) mean(is.finite(OUT[[v]])), 0), 3)
write.csv(md, file.path(HDIR, "predictors_admin2_mics_metadata.csv"), row.names = FALSE)
cat("\n=== MICS block ===\n"); print(md[, c("column", "n_countries", "countries", "completeness")], row.names = FALSE)
cat("\nunit levels:\n"); print(table(OUT$country, OUT$mics_level, useNA = "ifany"))
cat("DONE\n")
