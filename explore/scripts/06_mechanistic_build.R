# =============================================================================
# explore/scripts/06_mechanistic_build.R   [probe MX-01, build step]
#
# NUTRIENT-SPECIFIC MECHANISTIC FEATURES FROM THE CROP BASKET
#
# The project's predictor set is nutrient-agnostic: the same 383 columns are
# offered to zinc, B12, folate, iron and vitamin A, and are collapsed into one
# index. The nutrition literature is not agnostic. It names specific
# mechanisms, and two of them are buildable from data already on disk:
#
#   ZINC       absorption is governed by the PHYTATE:ZINC MOLAR RATIO of the
#              diet, not by zinc intake alone (the Wessells & Brown 2012 basis
#              for national zinc-deficiency estimates). A cassava district and
#              a sorghum district can have the same zinc intake and very
#              different absorbable zinc.
#   VITAMIN A  in West Africa the dominant plant source is CRUDE RED PALM OIL,
#              then mango/papaya and dark green leafy vegetables. Oil palm is
#              a SPAM crop; "oilcrops" as a group is not the same variable.
#
# The pipeline currently sees only 5 collapsed SPAM group shares. The raw
# MapSPAM files on disk carry all 42 crops, so the mechanism is buildable.
#
# This script keeps all 42 crops, joins the composition table in
# explore/data/food_composition.csv, and derives per-district nutrient density
# of the PRODUCTION basket. Production is not consumption - that leap is the
# main limitation and is stated in the findings entry.
#
#   Rscript explore/scripts/06_mechanistic_build.R
# -> explore/out/06_mechanistic_features.csv (+ _metadata.csv)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(dplyr)
})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("explore/R/harness.R")
source("R/data_prep.R")           # load_gadm_cached()

RAW <- "data/MapSPAM/raw"
ISO3 <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")

FC <- read.csv("explore/data/food_composition.csv", stringsAsFactors = FALSE)
CROPS <- FC$crop
stopifnot(length(CROPS) == 42, !anyDuplicated(CROPS))

# molar masses for the phytate:zinc ratio
MM_PHY <- 660.0   # phytic acid, g/mol
MM_ZN  <- 65.38   # zinc, g/mol

# ── read the 42-crop production grid, four countries only ───────────────────
# The global grid is 286 MB / 833k rows; 06a_spam_prefilter.py streams it down
# to the 4,739 pixels in these four countries so nothing large is held in
# memory here. Run that first if the cached file is absent.
PRE <- "explore/out/06a_spam_prod_4countries.csv"
if (!file.exists(PRE))
  stop("run: python explore/scripts/06a_spam_prefilter.py   (writes ", PRE, ")")
message("reading prefiltered MapSPAM production grid ...")
P <- data.table::fread(PRE)
P[, iso3 := toupper(trimws(iso3))]
have <- intersect(CROPS, names(P))
message("  ", nrow(P), " pixels, ", length(have), " of 42 crop columns present")

# ── aggregate to Admin-2 ────────────────────────────────────────────────────
out <- list()
for (cn in names(ISO3)) {
  g <- tryCatch(load_gadm_cached(ISO3[[cn]], level = 2), error = function(e) NULL)
  if (is.null(g)) { message("  !! GADM failed for ", cn); next }
  g <- sf::st_as_sf(g)[, c("NAME_1", "NAME_2")]
  pc <- P[iso3 == ISO3[[cn]]]
  if (!nrow(pc)) { message("  no pixels for ", cn); next }
  pts <- sf::st_as_sf(as.data.frame(pc), coords = c("x", "y"), crs = 4326, remove = FALSE)
  j <- sf::st_drop_geometry(sf::st_join(pts, g, join = sf::st_intersects))
  j <- as.data.table(j)[!is.na(NAME_2)]
  agg <- j[, lapply(.SD, sum, na.rm = TRUE), by = .(NAME_1, NAME_2), .SDcols = have]
  agg[, country := cn]
  out[[cn]] <- agg
  message(sprintf("  %-12s %3d districts", cn, nrow(agg)))
}
B <- data.table::rbindlist(out, fill = TRUE)
data.table::setnames(B, c("NAME_1", "NAME_2"), c("Admin1", "Admin2"))
B[, Admin1 := trimws(Admin1)][, Admin2 := trimws(Admin2)]

# ── nutrient density of the production basket ───────────────────────────────
M <- as.matrix(B[, ..have]); M[!is.finite(M)] <- 0
tot <- rowSums(M)
W <- M / pmax(tot, 1e-9)                    # production share per crop
fc <- FC[match(have, FC$crop), ]

dens <- function(col) as.numeric(W %*% fc[[col]])
zn   <- dens("zinc_mg")
fe   <- dens("iron_mg")
phy  <- dens("phytate_mg")
va   <- dens("provita_ug_rae")
fol  <- dens("folate_ug")

share_of <- function(codes) {
  cc <- intersect(codes, have)
  if (!length(cc)) return(rep(0, nrow(B)))
  as.numeric(rowSums(W[, cc, drop = FALSE]))
}

# A district whose entire basket is cash crops has zn = 0, and the ratio is not
# 0 there - it is undefined. Median imputation in prep_predictors_v2() is the
# honest default for those, not a fabricated extreme.
na_if_empty <- function(num, den, floor_den = 1e-6)
  ifelse(den > floor_den, num / pmax(den, floor_den), NA_real_)

F <- data.frame(
  country = B$country, Admin1 = B$Admin1, Admin2 = B$Admin2,
  # ZINC: the Wessells & Brown mechanism
  mx_phytate_zn_molar = na_if_empty(phy / MM_PHY, zn / MM_ZN),
  mx_zinc_density     = zn,
  mx_phytate_density  = phy,
  # IRON: non-haem iron and the phytate load that blocks it
  mx_iron_density     = fe,
  mx_phytate_fe_molar = na_if_empty(phy / MM_PHY, fe / 55.85),
  # VITAMIN A: carotenoid density, and oil palm on its own
  mx_provita_density  = va,
  mx_oilpalm_share    = share_of("oilp"),
  mx_provita_crop_share = share_of(c("oilp", "swpo", "trof", "vege", "plnt")),
  # FOLATE: pulses and green vegetables
  mx_folate_density   = fol,
  mx_pulse_share      = share_of(c("bean", "chic", "cowp", "pige", "lent", "opul", "soyb")),
  # composition axes the 5-group version cannot express
  mx_cereal_share     = share_of(c("whea", "rice", "maiz", "barl", "pmil", "smil", "sorg", "ocer")),
  mx_root_share       = share_of(c("pota", "swpo", "yams", "cass", "orts")),
  mx_cash_share       = share_of(c("sugc", "sugb", "cott", "ofib", "acof", "rcof", "coco", "teas", "toba")),
  mx_rice_share       = share_of("rice"),     # milled rice: low phytate, low zinc
  mx_cassava_share    = share_of("cass"),     # the classic low-micronutrient staple
  mx_crop_diversity   = apply(W, 1, function(w) {   # Shannon diversity of the basket
    w <- w[w > 0]; if (!length(w)) return(0); -sum(w * log(w))
  }),
  stringsAsFactors = FALSE
)
F <- F[tot > 0, ]

exp_write(F, "06_mechanistic_features")

META <- data.frame(
  column = setdiff(names(F), c("country", "Admin1", "Admin2")),
  formula = c(
    "(phytate mg/660) / (zinc mg/65.38), production-share weighted over 42 SPAM crops",
    "production-share weighted zinc mg/100g", "production-share weighted phytate mg/100g",
    "production-share weighted iron mg/100g",
    "(phytate mg/660) / (iron mg/55.85), production-share weighted",
    "production-share weighted provitamin A ug RAE/100g",
    "oil palm share of production", "oil palm + sweet potato + tropical fruit + vegetables + plantain share",
    "production-share weighted folate ug/100g", "pulse share of production",
    "cereal share", "root and tuber share", "non-food cash crop share",
    "rice share", "cassava share", "Shannon diversity of the 42-crop production basket"),
  nutrient = c("zinc", "zinc", "zinc", "iron", "iron", "vitA", "vitA", "vitA",
               "folate", "folate", "general", "general", "general", "general",
               "general", "general"),
  source = "MapSPAM 2010 v2r0 42-crop production x explore/data/food_composition.csv",
  stringsAsFactors = FALSE)
exp_write(META, "06_mechanistic_features_metadata")

# ── sanity check against agronomy ───────────────────────────────────────────
cat("\n== sanity check: extremes per feature (expect cereal/sesame districts",
    "high on phytate:Zn, cassava districts low) ==\n")
for (v in c("mx_phytate_zn_molar", "mx_provita_density", "mx_oilpalm_share",
            "mx_cassava_share", "mx_folate_density")) {
  o <- order(F[[v]], decreasing = TRUE)
  cat("\n", v, "\n  high: ",
      paste(sprintf("%s/%s %.2f", F$country[head(o, 4)], F$Admin2[head(o, 4)],
                    F[[v]][head(o, 4)]), collapse = " | "), "\n  low : ",
      paste(sprintf("%s/%s %.2f", F$country[tail(o, 4)], F$Admin2[tail(o, 4)],
                    F[[v]][tail(o, 4)]), collapse = " | "), "\n", sep = "")
}
cat("\nrows:", nrow(F), " countries:", paste(unique(F$country), collapse = ", "), "\n")
