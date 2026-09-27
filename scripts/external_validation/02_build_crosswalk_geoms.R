# =============================================================================
# scripts/external_validation/02_build_crosswalk_geoms.R              [XV-01]
#
# GADM ADMIN-2 -> VMNIS SURVEY UNIT, AND THE GEOMETRY THE GEE PULL READS
#
# Predictors are area-averaged as a SIMPLE MEAN OF DISTRICTS when aggregating
# to a region (PREREGISTRATION_NEW_COUNTRIES_2026-09.md), so the extraction
# stays at admin-2 and carries the survey unit as a label rather than
# dissolving the polygons first.
#
# Boundary vintage is reconciled to each SURVEY's framework, not GADM's.
# GADM 4.1 is current; the surveys are not:
#   Zambia   Muchinga (2011) splits back to Northern, except Chama -> Eastern.
#   Sudan    Central Darfur (2012) -> Western Darfur; East Darfur (2012) ->
#            Southern Darfur; West Kurdufan (2013) splits by pre-2013 parent
#            (En Nuhud, Ghebeish -> Northern Kordofan; As Salam, Lagawa ->
#            Southern Kordofan); Abyei is dropped as disputed territory.
#   Nigeria  37 states -> the 6 geopolitical zones the NFCMS reports on.
# Every non-identity mapping carries a `note` so it is auditable.
#
#   Rscript scripts/external_validation/02_build_crosswalk_geoms.R
# -> data/external_validation/gadm_to_vmnis_crosswalk.csv
#    data/external_validation/admin2_new_countries.geojson
# =============================================================================
suppressPackageStartupMessages({library(sf); library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "data/external_validation"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

NGA_ZONE <- c(
  Benue="North Central", Kogi="North Central", Kwara="North Central",
  Nasarawa="North Central", Niger="North Central", Plateau="North Central",
  `Federal Capital Territory`="North Central",
  Adamawa="North East", Bauchi="North East", Borno="North East",
  Gombe="North East", Taraba="North East", Yobe="North East",
  Jigawa="North West", Kaduna="North West", Kano="North West",
  Katsina="North West", Kebbi="North West", Sokoto="North West", Zamfara="North West",
  Abia="South East", Anambra="South East", Ebonyi="South East",
  Enugu="South East", Imo="South East",
  `Akwa Ibom`="South South", Bayelsa="South South", `Cross River`="South South",
  Delta="South South", Edo="South South", Rivers="South South",
  Ekiti="South West", Lagos="South West", Ogun="South West",
  Ondo="South West", Osun="South West", Oyo="South West")

# ── off-continent arm (XV-02) ───────────────────────────────────────────────
# GADM spelling -> the spelling the VMNIS deposit uses, typos included
# ('Kasmir', 'Glgit', 'Teritory', 'Tribual' are VMNIS's own).
PAK_MAP <- c(`Azad Kashmir`="Azad Jammu & Kasmir (AJK)",
             `Federally Administered Tribal Ar`="Federally Administered Tribual Areas (FATA)",
             `Gilgit-Baltistan`="Glgit-Baltistan",
             Islamabad="Islamabad Capital Teritory",
             `Khyber-Pakhtunkhwa`="Khyber Pakhtunkhwa (previously NWFP)")

IND_MAP <- c(`Andhra Pradesh`="Andhra pradesh", Odisha="Orissa",
             `Jammu and Kashmir`="Jammu & Kashmir", `NCT of Delhi`="Delhi")

# the 30 states CNNS reports; GADM's other UTs (Andaman, Chandigarh, Dadra and
# Nagar Haveli, Daman and Diu, Lakshadweep, Puducherry) have no survey row and
# are dropped rather than folded into a neighbour
IND_UNITS <- c("Andhra pradesh", "Arunachal Pradesh", "Assam", "Bihar",
  "Chhattisgarh", "Delhi", "Goa", "Gujarat", "Haryana", "Himachal Pradesh",
  "Jammu & Kashmir", "Jharkhand", "Karnataka", "Kerala", "Madhya Pradesh",
  "Maharashtra", "Manipur", "Meghalaya", "Mizoram", "Nagaland", "Orissa",
  "Punjab", "Rajasthan", "Sikkim", "Tamil Nadu", "Telangana", "Tripura",
  "Uttar Pradesh", "Uttarakhand", "West Bengal")

ETH_MAP <- c(`Addis Abeba`="Addis Ababa", `Benshangul-Gumaz`="Benishangul-Gumuz",
             `Gambela Peoples`="Gambela", `Harari People`="Harari",
             Oromia="Oromiya", `Southern Nations, Nationalities`="SNNP")

SDN_MAP <- c(`Al Jazirah`="Al Jazeera", `Al Qadarif`="Gadaref", `River Nile`="Nile",
             `North Darfur`="Northern Darfur", `North Kurdufan`="Northern Kordofan",
             `South Darfur`="Southern Darfur", `South Kurdufan`="Southern Kordofan",
             `West Darfur`="Western Darfur",
             `Central Darfur`="Western Darfur", `East Darfur`="Southern Darfur")

SDN_WK <- c(`En Nuhud`="Northern Kordofan", Ghebeish="Northern Kordofan",
            `As Salam`="Southern Kordofan", Lagawa="Southern Kordofan")

map_unit <- function(iso, a1, a2) {
  note <- ""
  u <- a1
  if (iso == "ZMB") {
    if (a1 == "Muchinga") {
      u <- if (a2 == "Chama") "Eastern" else "Northern"
      note <- "Muchinga created 2011; returned to pre-2011 parent"
    } else if (a1 == "North-Western") u <- "North Western"
  } else if (iso == "ETH") {
    if (a1 %in% names(ETH_MAP)) { u <- unname(ETH_MAP[a1]); note <- "spelling" }
  } else if (iso == "SDN") {
    if (a1 == "West Kurdufan") {
      if (a2 == "Abyei") { u <- NA_character_; note <- "Abyei disputed; dropped" }
      else { u <- unname(SDN_WK[a2]); note <- "West Kurdufan created 2013; split by pre-2013 parent" }
    } else if (a1 %in% names(SDN_MAP)) {
      u <- unname(SDN_MAP[a1])
      note <- if (a1 %in% c("Central Darfur","East Darfur")) "created 2012; returned to parent" else "spelling"
    }
  } else if (iso == "NGA") {
    u <- unname(NGA_ZONE[a1]); note <- "state -> NFCMS geopolitical zone"
  } else if (iso == "PAK") {
    if (a1 %in% names(PAK_MAP)) { u <- unname(PAK_MAP[a1]); note <- "spelling" }
  } else if (iso == "IND") {
    if (a1 %in% names(IND_MAP)) { u <- unname(IND_MAP[a1]); note <- "spelling" }
    if (!u %in% IND_UNITS) { u <- NA_character_; note <- "UT not covered by CNNS; dropped" }
  }
  c(u, note)
}

GROUP <- c(ZMB="africa", ETH="africa", SDN="africa", NGA="africa",
           PAK="offcontinent", IND="offcontinent")
CNAME <- c(ZMB="Zambia", ETH="Ethiopia", SDN="Sudan", NGA="Nigeria",
           PAK="Pakistan", IND="India")

geoms <- list(); xw <- list()
for (iso in names(GROUP)) {
  s <- readRDS(sprintf("data/admin_boundaries/gadm41_%s_2.rds", iso))
  s <- sf::st_make_valid(sf::st_transform(s, 4326))
  mm <- t(mapply(map_unit, iso, s$NAME_1, s$NAME_2))
  s$vmnis_unit <- mm[, 1]; s$note <- mm[, 2]
  s$iso3 <- iso
  s$country <- CNAME[[iso]]
  if (any(is.na(s$vmnis_unit)))
    cat(sprintf("[%s] dropping %d district(s) with no survey unit: %s\n", iso,
                sum(is.na(s$vmnis_unit)), paste(s$NAME_2[is.na(s$vmnis_unit)], collapse=", ")))
  s <- s[!is.na(s$vmnis_unit), ]
  xw[[iso]] <- sf::st_drop_geometry(s)[, c("country","iso3","NAME_1","NAME_2","vmnis_unit","note")]
  geoms[[iso]] <- s[, c("country","iso3","NAME_1","NAME_2","vmnis_unit")]
  cat(sprintf("[%s] %d districts -> %d survey units\n", iso, nrow(s), length(unique(s$vmnis_unit))))
}

XW <- do.call(rbind, xw)
names(XW)[names(XW)=="NAME_1"] <- "Admin1"; names(XW)[names(XW)=="NAME_2"] <- "Admin2"
write.csv(XW, file.path(OUT, "gadm_to_vmnis_crosswalk.csv"), row.names = FALSE)

# two geometry files: the Africa arm (iSDA soil available) and the
# off-continent arm (global SoilGrids substitute only)
for (grp in unique(GROUP)) {
  iso_in <- names(GROUP)[GROUP == grp]
  G <- do.call(rbind, geoms[iso_in])
  names(G)[names(G)=="NAME_1"] <- "Admin1"; names(G)[names(G)=="NAME_2"] <- "Admin2"
  G <- sf::st_simplify(G, dTolerance = 0.005, preserveTopology = TRUE)
  G <- G[!sf::st_is_empty(G), ]
  f <- file.path(OUT, if (grp == "africa") "admin2_new_countries.geojson"
                      else "admin2_offcontinent.geojson")
  if (file.exists(f)) file.remove(f)
  sf::st_write(G, f, driver = "GeoJSON", quiet = TRUE)
  cat(sprintf("wrote %s (%d districts)\n", basename(f), nrow(G)))
}
cat(sprintf("wrote gadm_to_vmnis_crosswalk.csv (%d rows)\n", nrow(XW)))
