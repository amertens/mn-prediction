# =============================================================================
# scripts/covariates/build_pcode_crosswalk.R   [JK-01, 2026-09-15]
#
# GADM 4.1 GID_2 -> OCHA COD-AB P-CODE CROSSWALK, BY AREA OF OVERLAP
#
# The spine (metadata/admin2_spine.csv) keys every unit by its GADM 4.1 code.
# Providers that publish subnational data with codes (OCHA, WFP, IPC / Cadre
# Harmonise, HDX HXL tables, ESPEN IUs in some countries) use the COD-AB
# P-codes instead, and the two boundary sets disagree in places: GADM's Gambia
# has 37 districts against the 2013-census 43+ that COD-AB carries, Malawi's
# GADM level 2 is the Traditional Authority while COD-AB ADM2 is the district.
# So the crosswalk is spatial, not by name: each GADM unit takes the P-code of
# the COD-AB ADM2 polygon that covers the largest share of its area, and the
# share is recorded so a many-to-one (or a poor overlap) is visible.
#
# Input: the COD-AB shapefile packages from HDX (downloaded 2026-09-15 with
# the user's approval; metadata/external_provenance.csv), unpacked to
# data/COD_AB/<ISO3>/ as <iso3>_admin<level>.shp with adm<level>_pcode /
# adm<level>_name (https://data.humdata.org/dataset/cod-ab-<iso3>, CC BY-IGO).
# The script skips a country whose folder is absent and says so.
# Result 2026-09-15: 554 of 554 units coded; overlap share < 0.8 for 9 Gambia
# units (COD's 49 districts are finer than GADM's 37), 65 Malawi units (COD
# TAs differ from GADM's) and 4 Sierra Leone units (the 2017 splits); two
# Malawi P-codes are shared by several GADM units.
#
#   Rscript -e "source('scripts/covariates/build_pcode_crosswalk.R')"
# -> metadata/crosswalks/gadm41_gid2_to_pcode.csv   (GADM unit -> P-code, largest overlap)
# -> metadata/crosswalks/pcode_to_gadm41_gid2.csv   (P-code -> GADM unit, for ingesting P-coded tables)
# -> metadata/admin2_spine.csv   (pcode column filled)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
sf::sf_use_s2(FALSE)
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
SP <- read.csv("metadata/admin2_spine.csv", stringsAsFactors = FALSE)

# The COD level that corresponds to the spine's Admin-2: ADM2 for The Gambia
# (49 districts against GADM's 37), Ghana (260, one to one) and Sierra Leone
# (16 against GADM's 14: Karene and Falaba were carved out of Bombali and
# Koinadugu in 2017); ADM3 for Malawi (433 TAs against GADM's 243 units).
LEVEL <- c(Gambia = 2L, Ghana = 2L, Malawi = 3L, SierraLeone = 2L)
find_shp <- function(dir, level) { f <- list.files(dir, pattern = sprintf("_admin%d[.]shp$", level), full.names = TRUE, recursive = TRUE); f[!grepl("_em[.]shp$", f)][1] }
pcode_col <- function(d, level) { nm <- names(d); c(nm[grepl(sprintf("^adm%d_pcode$", level), nm, ignore.case = TRUE)], nm[grepl("pcode", nm, ignore.case = TRUE)])[1] }
name_col  <- function(d, level) { nm <- names(d); c(nm[grepl(sprintf("^adm%d_(name|en)$", level), nm, ignore.case = TRUE)], nm[grepl(sprintf("^adm%d_", level), nm, ignore.case = TRUE)])[1] }

rows <- list(); rev_rows <- list()
for (cn in names(ISO)) {
  dir <- file.path("data", "COD_AB", ISO[[cn]]); lv <- LEVEL[[cn]]
  shp <- if (dir.exists(dir)) find_shp(dir, lv) else NA
  if (is.na(shp)) { cat(sprintf("%s: no COD-AB admin%d shapefile under %s, skipped\n", cn, lv, dir)); next }
  cod <- sf::st_read(shp, quiet = TRUE) |> sf::st_transform(4326) |> sf::st_make_valid()
  pc <- pcode_col(cod, lv); nc <- name_col(cod, lv)
  g <- readRDS(sprintf("data/admin_boundaries/gadm41_%s_2.rds", ISO[[cn]])) |> sf::st_transform(4326) |> sf::st_make_valid()
  g <- g[g$GID_2 %in% SP$gid_2[SP$country == cn], ]
  g$area <- as.numeric(sf::st_area(g))
  inter <- suppressWarnings(sf::st_intersection(g[, c("GID_2", "NAME_1", "NAME_2", "area")], cod[, c(pc, nc)]))
  inter$ov <- as.numeric(sf::st_area(inter)) / inter$area
  best <- sf::st_drop_geometry(inter) |> group_by(GID_2) |> arrange(desc(ov), .by_group = TRUE) |>
    summarise(pcode = first(.data[[pc]]), pcode_name = first(.data[[nc]]), overlap_share = round(first(ov), 3),
              second_share = round(if (n() > 1) ov[2] else 0, 3), n_cod_polys = n(), .groups = "drop")
  out <- SP[SP$country == cn, c("country", "Admin1", "Admin2", "gid_2")] |> left_join(best, by = c("gid_2" = "GID_2"))
  out$cod_level <- lv; out$many_to_one <- out$pcode %in% out$pcode[duplicated(out$pcode)]
  # the reverse direction, used when a P-coded table is ingested: each COD
  # polygon -> the GADM unit that holds the largest share of ITS area
  cod$cod_area <- as.numeric(sf::st_area(cod))
  inter2 <- suppressWarnings(sf::st_intersection(cod[, c(pc, nc, "cod_area")], g[, c("GID_2", "NAME_1", "NAME_2")]))
  inter2$ov <- as.numeric(sf::st_area(inter2)) / inter2$cod_area
  rev <- sf::st_drop_geometry(inter2) |> group_by(pcode = .data[[pc]]) |> arrange(desc(ov), .by_group = TRUE) |>
    summarise(pcode_name = first(.data[[nc]]), gid_2 = first(GID_2), Admin1 = first(NAME_1), Admin2 = first(NAME_2),
              share_of_cod_unit = round(first(ov), 3), n_gadm_polys = n(), .groups = "drop") |> mutate(country = cn, cod_level = lv)
  rev_rows[[cn]] <- rev
  cat(sprintf("%s: %d units, %d matched, %d with overlap < 0.8, %d P-codes shared by several GADM units\n", cn, nrow(out), sum(!is.na(out$pcode)),
              sum(out$overlap_share < 0.8, na.rm = TRUE), length(unique(out$pcode[out$many_to_one]))))
  rows[[cn]] <- out
}
if (!length(rows)) stop("no COD-AB files found; download them to data/COD_AB/<ISO3>/ first")
X <- bind_rows(rows)
dir.create("metadata/crosswalks", showWarnings = FALSE)
write.csv(X, "metadata/crosswalks/gadm41_gid2_to_pcode.csv", row.names = FALSE)
RV <- bind_rows(rev_rows); write.csv(RV, "metadata/crosswalks/pcode_to_gadm41_gid2.csv", row.names = FALSE)
cat(sprintf("reverse crosswalk: %d COD units -> GADM (share of the COD unit inside its GADM unit: median %.2f, %d below 0.6)
", nrow(RV), median(RV$share_of_cod_unit), sum(RV$share_of_cod_unit < 0.6)))
SP$pcode <- X$pcode[match(SP$gid_2, X$gid_2)]
write.csv(SP, "metadata/admin2_spine.csv", row.names = FALSE)
cat(sprintf("crosswalk written for %d units; spine pcode filled for %d\nDONE\n", nrow(X), sum(!is.na(SP$pcode))))
