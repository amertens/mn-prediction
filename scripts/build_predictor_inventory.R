# =============================================================================
# scripts/build_predictor_inventory.R
#
# Enumerate EVERY predictor the pipeline can see, with its source, conceptual
# domain, per-country availability, and which model family actually uses it.
#
# Reads the built targets store (default _targets_full) rather than re-deriving
# anything, so the inventory is what the pipeline really had, not what it is
# supposed to have.
#
#   Rscript scripts/build_predictor_inventory.R [store]
#
# Writes results/tables/predictor_inventory.csv (one row per variable x level)
# and prints a source x level summary.
# =============================================================================
suppressPackageStartupMessages({
  library(targets); library(dplyr); library(readr); library(here)
})

args  <- commandArgs(trailingOnly = TRUE)
STORE <- if (length(args)) args[1] else here("_targets_full")
message("store: ", STORE)

COUNTRIES <- c(Gambia = "gambia", Ghana = "ghana", SierraLeone = "sierraleone",
               Malawi = "malawi", Tanzania = "tanzania")

read_t <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE),
                               error = function(e) { message("  miss: ", nm); NULL })

# ── 1. Column inventories ───────────────────────────────────────────────────
ind_cols <- list(); area_cols <- list()
for (cn in names(COUNTRIES)) {
  key <- COUNTRIES[[cn]]
  d <- read_t(paste0("merged_", key))
  if (!is.null(d)) ind_cols[[cn]] <- colnames(d)
  g <- read_t(paste0("gee_admin2_", key))
  if (!is.null(g)) area_cols[[cn]] <- setdiff(colnames(g), c("Admin1", "Admin2"))
  message(sprintf("%-12s individual %5d cols | area %5d cols", cn,
                  length(ind_cols[[cn]] %||% character(0)),
                  length(area_cols[[cn]] %||% character(0))))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ── 2. Source classification by prefix ──────────────────────────────────────
# Prefixes come from R/config.R (cc$domains) and R/data_prep.R (external_domains).
SOURCES <- tibble::tribble(
  ~prefix,          ~source,                                   ~source_kind,
  "gw_",            "GroundWork/biomarker survey items",        "survey (excluded from models)",
  "dhs2010_",       "DHS 2010 admin-2 aggregate",               "survey aggregate",
  "dhs2013_",       "DHS 2013 admin-2 aggregate",               "survey aggregate",
  "dhs2014_",       "DHS 2014 admin-2 aggregate",               "survey aggregate",
  "dhs2016_",       "DHS 2016 admin-2 aggregate",               "survey aggregate",
  "dhs2017_",       "DHS 2017 admin-2 aggregate",               "survey aggregate",
  "dhs2019_",       "DHS 2019 admin-2 aggregate",               "survey aggregate",
  "mics_",          "UNICEF MICS admin-2 aggregate",            "survey aggregate",
  "ihme_",          "IHME GBD / LBD raster",                    "modelled raster",
  "lsms_",          "LSMS-ISA",                                 "survey aggregate",
  "MAP_",           "Malaria Atlas Project",                    "modelled raster",
  "map2_",          "Malaria Atlas Project (v2 puller)",        "modelled raster",
  "wfp_",           "WFP VAM / market prices",                  "operational data",
  "flunet_",        "WHO FluNet",                               "surveillance",
  "gee_",           "Google Earth Engine zonal means",          "earth observation",
  "chirps_",        "CHIRPS precipitation",                     "earth observation",
  "worldpop_",      "WorldPop",                                 "modelled raster",
  "ntl_",           "Night-time lights",                        "earth observation",
  "soil_",          "SoilGrids / ISRIC",                        "earth observation",
  "gdl_",           "Global Data Lab subnational HDI",          "modelled index",
  "crop_",          "Crop production / MapSPAM",                "modelled raster",
  "ipc_",           "IPC / Cadre Harmonise phase",              "operational data",
  "acled_",         "ACLED conflict events",                    "event data",
  "espen_",         "ESPEN helminth / NTD surveys",             "survey aggregate",
  "mapspam_",       "MapSPAM crop allocation",                  "modelled raster",
  "fsec_",          "Food security composite",                  "derived",
  "faostat_",       "FAOSTAT food balance sheets",              "national statistics",
  "fpn_",           "Food Prices for Nutrition",                "derived index"
)

classify_source <- function(v) {
  # longest matching prefix wins, so dhs2019_ beats a hypothetical dhs_
  ord <- SOURCES[order(-nchar(SOURCES$prefix)), ]
  for (i in seq_len(nrow(ord)))
    if (startsWith(v, ord$prefix[i])) return(ord[i, c("prefix", "source", "source_kind")])
  tibble::tibble(prefix = NA_character_, source = "unclassified / survey-native",
                 source_kind = "unclassified")
}

# ── 3. Conceptual domain from metadata/variable_conceptual_domains.csv ───────
cd <- readr::read_csv(here("metadata", "variable_conceptual_domains.csv"),
                      show_col_types = FALSE)
cd$regex <- paste0("^", gsub("####", "[0-9]{4}", cd$variable_pattern))
concept_of <- function(v) {
  hit <- which(vapply(cd$regex, function(r) grepl(r, v), logical(1)))
  if (!length(hit)) return(c(NA_character_, NA_character_, NA_character_))
  i <- hit[which.max(nchar(cd$variable_pattern[hit]))]   # most specific pattern
  c(cd$level1_domain[i], cd$level2_domain[i], cd$description[i])
}

# ── 4. Assemble ─────────────────────────────────────────────────────────────
build_rows <- function(cols_by_country, level) {
  all_v <- sort(unique(unlist(cols_by_country)))
  if (!length(all_v)) return(NULL)
  present <- vapply(all_v, function(v)
    paste(names(cols_by_country)[vapply(cols_by_country, function(cc) v %in% cc, logical(1))],
          collapse = ";"), character(1))
  n_ctry <- vapply(all_v, function(v)
    sum(vapply(cols_by_country, function(cc) v %in% cc, logical(1))), integer(1))
  src <- do.call(rbind, lapply(all_v, classify_source))
  con <- do.call(rbind, lapply(all_v, concept_of))
  tibble::tibble(
    variable = all_v, level = level,
    prefix = src$prefix, source = src$source, source_kind = src$source_kind,
    conceptual_domain = con[, 1], conceptual_subdomain = con[, 2],
    description = con[, 3],
    countries_present = present, n_countries = n_ctry,
    in_all_countries = n_ctry == length(cols_by_country))
}

inv <- dplyr::bind_rows(
  build_rows(ind_cols,  "individual"),
  build_rows(area_cols, "area (admin-2)"))

# Model reachability: the area-level model sees ONLY the gee_admin2 table; the
# individual-level model sees everything except the gw_ biomarker items.
inv <- inv %>%
  mutate(
    reaches_area_model = level == "area (admin-2)" & in_all_countries,
    reaches_individual_model = level == "individual" &
      !startsWith(variable, "gw_") & !is.na(prefix),
    note = case_when(
      startsWith(variable, "gw_") ~ "excluded: biomarker survey block (proxy-only rule)",
      level == "area (admin-2)" & !in_all_countries ~
        "present in some countries only -> dropped by the LOCO name intersection",
      TRUE ~ NA_character_))

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
out <- here("results", "tables", "predictor_inventory.csv")
readr::write_csv(inv, out)
message("\nwrote ", out, " (", nrow(inv), " rows)\n")

print(inv %>% count(level, source, source_kind, sort = TRUE), n = 60)
message("\n-- area-level variables shared by ALL countries --")
print(inv %>% filter(level == "area (admin-2)") %>% count(in_all_countries))
