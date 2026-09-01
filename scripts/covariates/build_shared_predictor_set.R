# =============================================================================
# scripts/covariates/build_shared_predictor_set.R
#
# Assemble the shared Admin-2 proxy predictor set, with per-column COVERAGE
# FLAGS and DOMAIN METADATA.
#
# WHAT CHANGES HERE, AND WHY
# --------------------------
# Every predictor file this project has built so far applies an all-four-
# countries filter: a column absent or constant in any one country is dropped.
# That is the right rule for a leave-one-country-out model, which needs the same
# vocabulary everywhere, and the wrong rule for a dataset shared with other
# modelling teams. A predictor available in three countries supports
# three-country comparisons and any country added later; dropping it silently
# destroys information other teams could use.
#
# So nothing is dropped for coverage. Every column carries:
#
#   n_countries      how many of the four it is usable in
#   countries        which ones
#   domain           the substantive domain, for domain-level variable importance
#   source           the data source it came from
#   subnational      TRUE if it varies within a country; FALSE if it is a
#                    national constant broadcast down, which cannot contribute
#                    to a district ranking
#
# The consumer filters. The builder reports.
#
# THE SUBNATIONAL FLAG IS NOT COSMETIC. A national covariate broadcast to every
# district is constant within a country. It can shift a country's level and can
# never change a district ranking, so including it in a targeting model without
# the flag would be a silent error of exactly the kind this project keeps
# finding.
#
#   Rscript scripts/covariates/build_shared_predictor_set.R
# -> data/covariates/harmonized/predictors_admin2_shared.csv
# -> data/covariates/harmonized/predictors_admin2_shared_metadata.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here); library(sf)})
targets::tar_source(here("R"))

OUT <- here("data", "covariates", "harmonized", "predictors_admin2_shared.csv")
META <- here("data", "covariates", "harmonized", "predictors_admin2_shared_metadata.csv")
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
cfgs <- get_country_configs()
kk <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(as.character(x))))

# ---------------------------------------------------------------------------
# Domain assignment. Prefix-driven, because every builder in this project
# prefixes its output, and hand-listing 400 columns would rot immediately.
# The order matters: the first pattern that matches wins.
# ---------------------------------------------------------------------------
DOMAIN_RULES <- list(
  c("^dhs_c_fg_|^dhs_c_diet|^dhs_c_mdd",        "Dietary diversity",              "DHS"),
  c("^dhs_hh_(cattle|cows|horses|goats|sheep|chickens|livestock|agland|owns_land)",
                                                "Agricultural production, land use","DHS"),
  c("^dhs_w_(iron|iptp|deworm)|^dhs_c_(vita_capsule|deworm)",
                                                "Food fortification and supplementation","DHS"),
  c("^dhs_w_sib",                               "Adult and maternal mortality",   "DHS"),
  c("^dhs_w_barrier|^dhs_w_health_insurance",   "Healthcare access",              "DHS"),
  c("^dhs_hh_(clean_fuel|solid_fuel|cook_indoors)", "Built environment",          "DHS"),
  c("^dhs_c_",                                  "Infant and child morbidity/mortality","DHS"),
  c("^dhs_w_(edu|literate|no_education|primary|secondary|working|occ|decides|owns_house)",
                                                "Education, employment, SES",     "DHS"),
  c("^dhs_w_(anc|delivery|birth|parity|teen|modern_fp|unmet)", "Fertility, reproductive health","DHS"),
  c("^dhs_w_(bmi|height)",                      "Adult nutrition",                "DHS"),
  c("^dhs_hh_(improved_water|improved_sanit|open_def|water|soap|handwash)",
                                                "Water and sanitation",           "DHS"),
  c("^dhs_hh_",                                 "Household assets and characteristics","DHS"),
  c("^dhs_(AN|CN)_",                            "Adult nutrition",                "DHS"),
  c("^dhs_(CH|CM)_",                            "Infant and child morbidity/mortality","DHS"),
  c("^dhs_(RH|FP)_",                            "Fertility, reproductive health", "DHS"),
  c("^dhs_WS_",                                 "Water and sanitation",           "DHS"),
  c("^dhs_ML_",                                 "Malaria",                        "DHS"),
  c("^dhs_w_",                                  "Fertility, reproductive health", "DHS"),
  c("^dhs_",                                    "Household assets and characteristics","DHS"),
  c("^map_.*blood_disorders",                   "Malaria",                        "Malaria Atlas"),
  c("^map_.*intervention",                      "Malaria incidence and treatment","Malaria Atlas"),
  c("^map_",                                    "Malaria incidence and treatment","Malaria Atlas"),
  c("^ihme_",                                   "Infant and child morbidity/mortality","IHME GBD"),
  c("^fsec_",                                   "Food prices and food security",  "FEWS NET / HFID"),
  c("^wfp_",                                    "Food prices and food security",  "WFP"),
  c("^soilgrids_|^soil_",                       "Soil characteristics",           "SoilGrids / iSDA"),
  c("^spam_|^aef_",                             "Agricultural production, land use","MapSPAM / AEF"),
  c("^ndvi|^evi|^lai|^npp|^wapor|^grassland",   "Ecosystem productivity/greenness","GEE"),
  c("^tclim|^precip|^lst|^aod",                 "Climate and weather",            "GEE"),
  c("^popdens|^ghs|^built|^wsf|^ntl|^human",    "Ruralness, population density, built environment","GEE"),
  c("^lcover|^elevation|^access",               "Built environment",              "GEE")
)
classify <- function(nm) {
  dom <- rep("Unclassified", length(nm)); src <- rep("unknown", length(nm))
  for (r in DOMAIN_RULES) {
    hit <- dom == "Unclassified" & grepl(r[1], nm, ignore.case = TRUE)
    dom[hit] <- r[2]; src[hit] <- r[3]
  }
  data.frame(column = nm, domain = dom, source = src, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 1. The existing harmonized set, and the MAP block built earlier.
# ---------------------------------------------------------------------------
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
res <- H
cat(sprintf("[shared] base harmonized: %d rows x %d predictors\n",
            nrow(H), ncol(H) - 3L))

extra_p <- here("data", "covariates", "harmonized", "predictors_admin2_extra.csv")
if (file.exists(extra_p)) {
  E <- suppressMessages(readr::read_csv(extra_p, show_col_types = FALSE)) |> as.data.frame()
  add <- setdiff(names(E), names(res))
  if (length(add)) {
    res <- dplyr::left_join(res, E[, c("country", "Admin1", "Admin2", add)],
                            by = c("country", "Admin1", "Admin2"))
    cat(sprintf("[shared] + %d columns from the extra-domain build\n", length(add)))
  }
}

# ---------------------------------------------------------------------------
# 2. DHS-derived indicators, aggregated to Admin-2 from the rdhs cache.
#
# Individual records carry a cluster id; the GE file carries that cluster's GPS;
# GADM carries the Admin-2 polygon. Point-in-polygon, then a survey-weighted
# mean per polygon. Done directly rather than through the surveyPrev pipeline,
# which downloads from the DHS API on every call.
# ---------------------------------------------------------------------------
as_num <- function(x) {
  if (inherits(x, "haven_labelled")) return(as.numeric(haven::zap_labels(x)))
  if (is.factor(x)) return(suppressWarnings(as.numeric(as.character(x))))
  suppressWarnings(as.numeric(x))
}
has_vars <- function(df, vars) all(vars %in% names(df))
source(here("src", "DHS", "DHS_extra_indicators_2026-09.R"))

CACHE <- "C:/Users/andre/AppData/Local/andre/rdhs/Cache/datasets"
# The DHS round is chosen to be closest to each country's MICRONUTRIENT survey,
# and its GPS file must be from the SAME round. Getting either wrong is silent:
# cluster numbering overlaps between rounds by coincidence, so a mismatched GPS
# file assigns real values to the wrong districts rather than failing.
#
# Sierra Leone was initially paired as SLIR7ADT with SLGE71FL. That is the 2019
# DHS -- six years after its 2013 micronutrient survey -- joined to a 2016 GPS
# file, giving 58 percent spurious cluster overlap. The correct set is the 61
# series: 2013, 435 clusters, matching the outcome year.
DHS_SETS <- list(
  Gambia      = list(IR="GMIR81DT", KR="GMKR81DT", HR="GMHR81DT", GE="GMGE81FL", year=2019),
  Ghana       = list(IR="GHIR72DT", KR="GHKR72DT", HR="GHHR72DT", GE="GHGE71FL", year=2014),
  Malawi      = list(IR="MWIR7ADT", KR="MWKR7ADT", HR="MWHR7ADT", GE="MWGE7AFL", year=2015),
  SierraLeone = list(IR="SLIR61DT", KR="SLKR61DT", HR="SLHR61DT", GE="SLGE61FL", year=2013))
MN_SURVEY_YEAR <- c(Gambia = 2018, Ghana = 2017, Malawi = 2016, SierraLeone = 2013)

rd <- function(x) { p <- file.path(CACHE, paste0(x, ".rds"))
                    if (file.exists(p)) readRDS(p) else NULL }

dhs_blocks <- list()
for (cn in COUNTRIES) {
  sp <- DHS_SETS[[cn]]
  ge <- rd(sp$GE)
  if (is.null(ge)) { cat(sprintf("[shared] %-12s no GPS file, DHS block skipped\n", cn)); next }
  gg <- sf::st_as_sf(ge)
  cl <- as_num(gg$DHSCLUST); lon <- as_num(gg$LONGNUM); lat <- as_num(gg$LATNUM)
  ok <- is.finite(cl) & is.finite(lon) & is.finite(lat) & !(lon == 0 & lat == 0)
  pts <- sf::st_as_sf(data.frame(clust = cl[ok], lon = lon[ok], lat = lat[ok]),
                      coords = c("lon", "lat"), crs = 4326)
  poly <- tryCatch(load_gadm_cached(cfgs[[cn]]$gadm_code, level = 2),
                   error = function(e) NULL)
  if (is.null(poly)) next
  poly <- sf::st_transform(poly, 4326)
  ix <- suppressMessages(as.integer(sf::st_within(pts, poly, sparse = TRUE) |>
                                      vapply(function(z) if (length(z)) z[1] else NA_integer_, integer(1))))
  a1 <- as.character(poly[[grep("^NAME_1$", names(poly), value = TRUE)[1]]])
  a2 <- as.character(poly[[grep("^NAME_2$", names(poly), value = TRUE)[1]]])
  cmap <- data.frame(clust = pts$clust, Admin1 = a1[ix], Admin2 = a2[ix],
                     stringsAsFactors = FALSE)
  cmap <- cmap[!is.na(cmap$Admin2), , drop = FALSE]

  # GUARD: the GPS file must belong to the same round as the recodes. Cluster
  # numbering repeats across rounds, so a mismatch produces a plausible join
  # onto the wrong districts instead of an error.
  ir_probe <- rd(sp$IR)
  if (!is.null(ir_probe) && "v001" %in% names(ir_probe)) {
    irc <- unique(as_num(ir_probe$v001)); gec <- unique(pts$clust)
    ov <- length(intersect(irc, gec)) / max(1L, length(irc))
    yr_gap <- abs(sp$year - MN_SURVEY_YEAR[[cn]])
    cat(sprintf("[shared] %-12s DHS %d (MN %d, gap %d yr), GPS overlap %.0f%%
",
                cn, sp$year, MN_SURVEY_YEAR[[cn]], yr_gap, 100 * ov))
    if (ov < 0.95) {
      cat(sprintf("[shared] %-12s GPS/recode overlap %.0f%% below 95%%, block REFUSED
",
                  cn, 100 * ov)); next
    }
  }

  acc <- NULL
  for (rc in c("IR", "KR", "HR")) {
    d <- rd(sp[[rc]]); if (is.null(d)) next
    add <- character()
    for (fn in EXTRA_DERIVERS[[rc]]) { r <- fn(d); d <- r$df; add <- c(add, r$added) }
    if (!length(add)) next
    ccol <- if (rc == "HR") "hv001" else "v001"
    wcol <- if (rc == "HR") "hv005" else "v005"
    if (!ccol %in% names(d)) next
    z <- data.frame(clust = as_num(d[[ccol]]),
                    w = if (wcol %in% names(d)) as_num(d[[wcol]]) else 1)
    z$w[!is.finite(z$w) | z$w <= 0] <- 1
    for (v in add) z[[v]] <- as_num(d[[v]])
    z <- dplyr::inner_join(z, cmap, by = "clust")
    if (!nrow(z)) next
    agg <- z |> group_by(Admin1, Admin2) |>
      summarise(across(all_of(add), ~ {
        ok2 <- is.finite(.x) & is.finite(w)
        if (!any(ok2)) NA_real_ else stats::weighted.mean(.x[ok2], w[ok2])
      }), .groups = "drop") |> as.data.frame()
    acc <- if (is.null(acc)) agg else
      dplyr::full_join(acc, agg, by = c("Admin1", "Admin2"))
  }
  if (is.null(acc)) next
  names(acc)[-(1:2)] <- paste0("dhs_", names(acc)[-(1:2)])
  acc$country <- cn
  dhs_blocks[[cn]] <- acc
  cat(sprintf("[shared] %-12s DHS: %d Admin-2 units, %d indicators\n",
              cn, nrow(acc), ncol(acc) - 3L))
}
if (length(dhs_blocks)) {
  DB <- dplyr::bind_rows(dhs_blocks)
  add <- setdiff(names(DB), names(res))
  j <- dplyr::left_join(res, DB[, c("country", "Admin1", "Admin2", add)],
                        by = c("country", "Admin1", "Admin2"))
  if (nrow(j) == nrow(res)) { res <- j
    cat(sprintf("[shared] + %d DHS-derived columns\n", length(add)))
  } else cat("[shared] DHS join changed row count, dropped\n")
}

# ---------------------------------------------------------------------------
# 3. Coverage and domain metadata. NOTHING IS DROPPED FOR COVERAGE.
# ---------------------------------------------------------------------------
cols <- setdiff(names(res), c("country", "Admin1", "Admin2"))
meta <- classify(cols)
cov <- t(vapply(cols, function(v) {
  s <- split(res[[v]], res$country)
  vapply(COUNTRIES, function(cc) {
    z <- s[[cc]]
    if (is.null(z)) return(0L)
    as.integer(mean(is.finite(z)) > 0.5 && stats::sd(z, na.rm = TRUE) > 0)
  }, integer(1))
}, integer(length(COUNTRIES))))
colnames(cov) <- COUNTRIES
meta$n_countries <- rowSums(cov)
meta$countries <- apply(cov, 1, function(z) paste(COUNTRIES[z == 1], collapse = "|"))
meta$completeness <- round(vapply(cols, function(v) mean(is.finite(res[[v]])), numeric(1)), 3)
# Subnational: does it vary WITHIN a country anywhere it is present?
meta$subnational <- vapply(cols, function(v) {
  s <- split(res[[v]], res$country)
  any(vapply(s, function(z) sum(is.finite(z)) > 3 &&
               stats::sd(z, na.rm = TRUE) > 0, logical(1)))
}, logical(1))

readr::write_csv(res, OUT)
readr::write_csv(meta, META)

cat(sprintf("\n[shared] %d Admin-2 rows x %d predictors\n", nrow(res), length(cols)))
cat("\n=== coverage ===\n")
print(as.data.frame(meta |> count(n_countries) |> arrange(desc(n_countries))), row.names = FALSE)
cat("\n=== domains ===\n")
print(as.data.frame(meta |> count(domain, source) |> arrange(desc(n))), row.names = FALSE)
cat(sprintf("\nnot subnational (national constants): %d\n", sum(!meta$subnational)))
cat(sprintf("unclassified domain: %d\n", sum(meta$domain == "Unclassified")))
cat(sprintf("\n-> %s\n-> %s\n", basename(OUT), basename(META)))
