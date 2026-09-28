# =============================================================================
# explore/scripts/17_geonutrition.R   [probe GN-01]
#
# MEASURED GRAIN AND SOIL CHEMISTRY AS PREDICTORS (Malawi)
#
# WHY THIS AND NOT MORE OF THE SAME. MX-01 built nutrient-specific features by
# multiplying crop PRODUCTION (MapSPAM) by a generic food-composition table -
# i.e. it assumed maize everywhere has the same zinc. It failed, and showed no
# nutrient specificity. GeoNutrition measures what that assumption got wrong:
# the actual mineral concentration of the grain grown at 1,812 georeferenced
# Malawi sites, plus the soil chemistry underneath it.
#
# It also fills a gap this session verified: the 575-column store contains NO
# environmental selenium at all (iSDA carries Al, Ca, CEC, Fe, Mg, P, K, S, C,
# Zn - no Se) and no environmental iodine, only two fortification-PROGRAMME
# fields. Malawi selenium is the most reliably measured outcome in the study
# (RL-01: district reliability 0.94) and had nothing mechanistic to predict it
# with.
#
# SOURCE. Gashu et al. 2022, Scientific Data, "Cereal grain mineral
# micronutrient and soil chemistry data from GeoNutrition surveys in Ethiopia
# and Malawi". figshare 10.6084/m9.figshare.15911973, CC BY 4.0. Malawi
# national sampling April-June 2018; 820 of the 1,900 sites were drawn from the
# 2015/16 Malawi DHS frame - the same frame as this project's biomarker survey.
#
# TIMING CAVEAT, STATED UP FRONT. The biomarker survey is Dec 2015 - Feb 2016;
# GeoNutrition sampled the 2018 harvest. Soil chemistry is slow-varying and the
# mismatch is acceptable for it; GRAIN concentration reflects the 2018 season
# and is a proxy for the district's typical grain, not for what was eaten in
# 2015/16.
#
#   Rscript explore/scripts/17_geonutrition.R
# -> explore/out/17_geonutrition_admin2.csv   the district block
#    explore/out/17_geonutrition_scores.csv   the test
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
EXP_ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
setwd(EXP_ROOT)
source("explore/R/harness.R")
source("R/config.R")
source("R/data_prep.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
GN_CSV <- "explore/data/geonutrition/MWI_CropSoilChemData_CSV/MWI_CropSoilData_NA.csv"
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# ── 1. aggregate the 1,812 points to Admin-2 ────────────────────────────────
G <- read.csv(GN_CSV, fileEncoding = "latin1", check.names = FALSE,
              stringsAsFactors = FALSE)
message("GeoNutrition Malawi: ", nrow(G), " sites, ", ncol(G), " columns")

# the variables worth carrying: measured grain minerals for the nutrients this
# project models, and the soil chemistry that drives them
GRAIN <- c("Se_grain", "Zn_grain", "Fe_grain", "Ca_grain", "S_grain",
           "Mg_grain", "P_grain", "Cu_grain", "Mn_grain")
SOIL  <- c("Se_Tot_Aqu", "Se_Org_Seq", "Se_Sol_Seq", "Se_Ads_Seq",
           "Zn_DTPA", "Zn_E", "Zn_LogKd", "Zn_Tot_Aqu",
           "Fe_DTPA", "Cu_DTPA", "Mn_DTPA",
           "S_Tot_Aqu", "S_Org_Seq", "S_Sol_Seq",
           "pH_w", "C_org", "N_pct", "eCEC", "P_Olsen")
VARS <- intersect(c(GRAIN, SOIL), names(G))
message("carrying ", length(VARS), " variables (",
        sum(GRAIN %in% VARS), " grain, ", sum(SOIL %in% VARS), " soil)")

G <- G[is.finite(num(G$Latitude)) & is.finite(num(G$Longitude)), ]
gadm <- sf::st_as_sf(load_gadm_cached("MWI", level = 2))[, c("NAME_1", "NAME_2")]
pts <- sf::st_as_sf(G, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
j <- sf::st_drop_geometry(sf::st_join(pts, gadm, join = sf::st_intersects))
j <- j[!is.na(j$NAME_2), ]
message("joined ", nrow(j), " of ", nrow(G), " sites to GADM districts")

# Grain minerals are right-skewed concentrations; aggregate on the log scale and
# take the district MEDIAN of the sites, which is robust to the few very high
# values that ICP-MS produces.
agg <- list()
for (v in VARS) {
  x <- num(j[[v]])
  x[!is.finite(x) | x <= 0] <- NA
  agg[[paste0("gn_", v)]] <- tapply(log(x), paste(j$NAME_1, j$NAME_2),
                                    function(z) stats::median(z, na.rm = TRUE))
}
keys <- names(agg[[1]])
GN <- data.frame(country = "Malawi",
                 Admin1 = sub("\\|.*$", "", gsub(" \\|", "|", keys)),
                 stringsAsFactors = FALSE)
sp <- strsplit(keys, " ")
GN <- data.frame(country = "Malawi", key = keys, stringsAsFactors = FALSE)
for (nm in names(agg)) GN[[nm]] <- as.numeric(agg[[nm]][keys])
# recover Admin1 / Admin2 from the key by matching against gadm
gd <- unique(sf::st_drop_geometry(gadm))
gd$key <- paste(gd$NAME_1, gd$NAME_2)
GN <- merge(GN, gd, by = "key", all.x = TRUE)
names(GN)[names(GN) == "NAME_1"] <- "Admin1"
names(GN)[names(GN) == "NAME_2"] <- "Admin2"
GN$Admin1 <- trimws(GN$Admin1); GN$Admin2 <- trimws(GN$Admin2)
GN$n_sites <- as.numeric(table(paste(j$NAME_1, j$NAME_2))[GN$key])
GN <- GN[, c("country", "Admin1", "Admin2", "n_sites",
             grep("^gn_", names(GN), value = TRUE))]
exp_write(GN, "17_geonutrition_admin2")
message("district block: ", nrow(GN), " districts, ",
        sum(grepl("^gn_", names(GN))), " variables")

# ── 2. how well does it cover the analysis districts? ───────────────────────
E <- exp_load()
tgt <- unique(E$TG[E$TG$country == "Malawi", c("Admin1", "Admin2")])
hit <- sum(paste(tgt$Admin1, tgt$Admin2) %in% paste(GN$Admin1, GN$Admin2))
hit2 <- sum(tgt$Admin2 %in% GN$Admin2)
cat(sprintf("\ncoverage of Malawi's %d target districts: pair-key %d, Admin2-only %d\n",
            nrow(tgt), hit, hit2))
cat(sprintf("median sites per covered district: %.0f (range %.0f-%.0f)\n",
            stats::median(GN$n_sites), min(GN$n_sites), max(GN$n_sites)))

# ── 3. Malawi selenium and iodine district targets (not in targets_v2) ──────
CFG <- get_country_configs()
cc <- CFG[[grep("malawi", names(CFG), ignore.case = TRUE)[1]]]
dat <- load_merged_data(cc$data_path)
extra_targets <- list()
for (on in names(cc$outcomes)) {
  if (!grepl("selenium|iodine", on)) next
  oc <- cc$outcomes[[on]]
  if (is.null(oc$continuous) || !oc$continuous %in% names(dat)) next
  keep <- outcome_population_mask(dat, cc, oc, label = "[GN-01]")
  d <- dat[keep, ]
  y <- num(d[[oc$continuous]])
  w <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d))
    num(d[[cc$weight_col]]) else rep(1, nrow(d))
  a1col <- if (!is.null(cc$admin1_col) && cc$admin1_col %in% names(d)) cc$admin1_col else "Admin1"
  ok <- is.finite(y) & y > 0 & !is.na(d[[cc$admin2_col]]) & !is.na(d[[a1col]])
  w[!is.finite(w) | w <= 0] <- NA
  ly <- -log(y)                              # negated log: higher = worse
  g <- paste(d[[a1col]][ok], d[[cc$admin2_col]][ok], sep = "")
  mu <- tapply(seq_len(sum(ok)), g, function(i) {
    yy <- ly[ok][i]; ww <- w[ok][i]; ww[!is.finite(ww)] <- stats::median(ww, na.rm = TRUE)
    if (all(!is.finite(ww))) mean(yy) else stats::weighted.mean(yy, ww)
  })
  nn <- tapply(ly[ok], g, length)
  extra_targets[[on]] <- data.frame(
    country = "Malawi", outcome = on,
    Admin1 = sub(".*$", "", names(mu)), Admin2 = sub("^.*", "", names(mu)),
    y_level = as.numeric(mu), n_eff_cont = as.numeric(nn),
    stringsAsFactors = FALSE)
  message("built district targets for Malawi ", on, ": ", length(mu), " districts")
}
XT <- dplyr::bind_rows(extra_targets)
if (nrow(XT)) exp_write(XT, "17_malawi_se_io_targets")

# ── 4. score ────────────────────────────────────────────────────────────────
GNCOLS <- grep("^gn_", names(GN), value = TRUE)
ridge_on <- function(pick) function(tr, te, y, X, D, aux) {
  cc2 <- intersect(pick, colnames(X))
  if (length(cc2) < 2) return(rep(mean(y[tr]), length(te)))
  .v2_enet(X[tr, cc2, drop = FALSE], y[tr], X[te, cc2, drop = FALSE], alpha = 0)
}
index_plus <- function(pick) function(tr, te, y, X, D, aux) {
  cc2 <- intersect(pick, colnames(X))
  D2 <- if (length(cc2)) cbind(D, X[, cc2, drop = FALSE]) else D
  arm_domain_index_v2(tr, te, y, X, D2, aux)
}
ARMS <- c(exp_baseline_arms(),
          list(gn_only = ridge_on(GNCOLS), index_plus_gn = index_plus(GNCOLS)))

rows <- list()
ix <- exp_cell_index(E)
ix <- ix[ix$country == "Malawi", ]
for (i in seq_len(nrow(ix))) {
  for (tgt2 in c("level", "prev")) {
    cell <- tryCatch(exp_cell(E, "Malawi", ix$outcome[i], tgt2, extra = GN,
                              extra_domain = "GeoNutrition grain and soil"),
                     error = function(e) NULL)
    if (is.null(cell)) next
    r <- exp_infill(cell, ARMS, reps = REPS)
    rows[[paste(ix$outcome[i], tgt2)]] <- r
  }
  message("  scored Malawi ", ix$outcome[i])
}

# selenium and iodine, built above, scored on the level only
if (nrow(XT)) {
  S <- E$S[E$S$country == "Malawi", ]
  CENT <- E$CENT[E$CENT$country == "Malawi", ]
  for (on in unique(XT$outcome)) {
    t <- XT[XT$outcome == on, ]
    m <- merge(t, S, by = c("Admin1", "Admin2"))
    m <- merge(m, GN[, c("Admin1", "Admin2", GNCOLS)], by = c("Admin1", "Admin2"), all.x = TRUE)
    m <- merge(m, CENT[, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
    stopifnot(!anyDuplicated(paste(m$Admin1, m$Admin2)))
    if (nrow(m) < 15) { message("  skip ", on, ": ", nrow(m), " districts"); next }
    preds <- c(intersect(E$PREDS, names(m)), GNCOLS)
    Xr <- prep_predictors_v2(as.matrix(m[, preds, drop = FALSE]))
    dom <- c(E$domain_of, stats::setNames(rep("GeoNutrition grain and soil",
                                              length(GNCOLS)), GNCOLS))
    D <- domain_representation_v2(Xr, dom)
    cell <- list(country = "Malawi", outcome = on, target = "level",
                 y_nat = m$y_level, y_mod = m$y_level, X = Xr, D = D,
                 aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1,
                            y_nat = m$y_level, target = "level", w = m$n_eff_cont),
                 w = m$n_eff_cont, Admin1 = m$Admin1, Admin2 = m$Admin2,
                 n = nrow(m))
    rows[[paste(on, "level")]] <- exp_infill(cell, ARMS, reps = REPS)
    message("  scored Malawi ", on, " (", nrow(m), " districts)")
  }
}

SM <- exp_summarise(dplyr::bind_rows(rows)); exp_write(SM, "17_geonutrition_scores")

cat("\n== Malawi in-fill, level target: median Spearman ==\n")
a <- SM[SM$estimand == "infill" & SM$target == "level", ]
print(aggregate(spearman ~ arm, data = a,
                FUN = function(z) round(median(z, na.rm = TRUE), 3)), row.names = FALSE)

cat("\n== paired: does GeoNutrition add to the index? (level) ==\n")
w <- reshape(a[, c("outcome", "arm", "spearman")], idvar = "outcome",
             timevar = "arm", direction = "wide")
names(w) <- sub("^spearman[.]", "", names(w))
w$gain <- round(w$index_plus_gn - w$domain_index, 3)
print(w[order(-w$gain), c("outcome", "domain_index", "index_plus_gn", "gn_only",
                          "spatial", "gain")], row.names = FALSE, digits = 3)
g <- w$gain[is.finite(w$gain)]
cat(sprintf("\nmean gain %+.4f | better in %d of %d Malawi cells | sign p = %.3f\n",
            mean(g), sum(g > 0), length(g),
            stats::binom.test(sum(g > 0), length(g), 0.5)$p.value))
