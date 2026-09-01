# =============================================================================
# scripts/accuracy_impact/wsf2_nutrition_proximal.R
#
# WS-F1 and WS-F2. Nutrition-proximal predictors from data already on disk.
#
# WHY THIS IS THE FIRST WS-F STEP AND NOT A DOWNLOAD
# --------------------------------------------------
# Stage 0 established that the harmonised 294 predictors contain NO
# supplementation, iodised-salt, deworming, haemoglobin or sickle-cell columns.
# Section 1.2 of the review document describes the 294 as spanning MAP, IHME,
# WFP and ESPEN; that is true of the individual-level merge and false of the
# harmonised Admin-2 file every area-level result rests on.
#
# The items are not missing from the project, only from that file. A 1.8 GB
# rdhs cache sits at the path in ~/.rdhs.json holding 550 datasets, and the
# five mechanistically relevant items are present in it with confirmed labels:
#
#   hv234   HR  result of salt test for iodine (ppm)
#   h33     KR  received vitamin a1 (most recent)      [h34 preferred if present]
#   h43     KR  drugs for intestinal parasites in last 6 months
#   hw53    KR  haemoglobin level (g/dl - 1 decimal)
#   v453    IR  haemoglobin level (g/dl - 1 decimal)
#   m45_1   IR  during pregnancy, given or bought iron tablets/syrup
#
# So this runs offline, needs no credentials, and supplies the haemoglobin
# inputs the WS-C3 iron-anaemia bridge needs.
#
# WHAT IT WRITES, AND WHERE IT DOES NOT WRITE
# -------------------------------------------
# A SEPARATE table, data/covariates/nutrition_proximal.csv. It is never merged
# into predictors_admin2_harmonized.csv, so no existing model silently changes
# its predictor set. Models opt in by naming a set from WS-F1.
#
# THE JOIN
# --------
# Cluster to district comes from data/DHS/clean/<survey>_cluster_admin_info.rds,
# whose `admin2.name.full` is already the Admin1_Admin2 PAIR key. Splitting it
# gives both columns, so the output carries the pair and downstream joins go
# through admin2_join_by() rather than the district name alone.
#
#   Rscript scripts/accuracy_impact/wsf2_nutrition_proximal.R
# -> data/covariates/nutrition_proximal.csv
# -> data/covariates/nutrition_proximal_coverage.csv
# -> data/metadata/external_provenance.csv   (appended)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

CACHE <- "C:/Users/andre/AppData/Local/andre/rdhs/Cache/datasets"
OUT   <- here("data", "covariates")
META  <- here("data", "metadata")
dir.create(META, showWarnings = FALSE, recursive = TRUE)

# Survey to recode-file map. The phase code identifies the round, and each round
# is the PRIOR DHS for its country, which is what makes these usable as
# predictors for a later micronutrient survey without circularity.
SURVEYS <- list(
  Gambia       = list(stem = "Gambia_2019",       hr = "GMHR81DT", kr = "GMKR81DT", ir = "GMIR81DT", round = "GDHS 2019-20"),
  Ghana        = list(stem = "Ghana_2014",        hr = "GHHR72DT", kr = "GHKR72DT", ir = "GHIR72DT", round = "GDHS 2014"),
  Malawi       = list(stem = "Malawi_2015",       hr = "MWHR7ADT", kr = "MWKR7ADT", ir = "MWIR7ADT", round = "MDHS 2015-16"),
  `Sierra Leone` = list(stem = "Sierra Leone_2013", hr = "SLHR61DT", kr = "SLKR61DT", ir = "SLIR61DT", round = "SLDHS 2013")
)

nz <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
# DHS reserves high codes for "don't know", "missing" and "inconsistent". Coding
# them as values rather than as missing is the classic way to manufacture a
# gradient out of nonresponse, so every reader below names its own codes.
miss <- function(x, codes) { x[x %in% codes] <- NA_real_; x }

get_col <- function(d, nm) if (nm %in% names(d)) nz(d[[nm]]) else rep(NA_real_, nrow(d))

read_cache <- function(f) {
  p <- file.path(CACHE, paste0(f, ".rds"))
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
}

# Weighted district means, keyed on the PAIR.
agg_to_admin2 <- function(vals, wts, clust, cmap, prefix) {
  ok <- is.finite(vals) & is.finite(wts) & wts > 0 & !is.na(clust)
  if (!any(ok)) return(NULL)
  key <- cmap$full[match(clust[ok], cmap$cluster)]
  keep <- !is.na(key)
  if (!any(keep)) return(NULL)
  d <- data.frame(key = key[keep], v = vals[ok][keep], w = wts[ok][keep])
  d |> group_by(key) |>
    summarise(!!paste0(prefix) := stats::weighted.mean(v, w),
              !!paste0(prefix, "_n") := dplyr::n(), .groups = "drop")
}

rows <- list(); cov_rows <- list()
for (cn in names(SURVEYS)) {
  s <- SURVEYS[[cn]]
  ci <- tryCatch(readRDS(here("data", "DHS", "clean",
                              paste0(s$stem, "_cluster_admin_info.rds"))),
                 error = function(e) NULL)
  if (is.null(ci) || is.null(ci$cluster.info$data)) {
    cat("[skip]", cn, "no cluster-admin map\n"); next }
  cmap <- data.frame(cluster = nz(ci$cluster.info$data$cluster),
                     full = as.character(ci$cluster.info$data$admin2.name.full),
                     stringsAsFactors = FALSE)
  cmap <- cmap[!is.na(cmap$cluster) & !is.na(cmap$full), , drop = FALSE]

  parts <- list()
  # ── household recode: iodised salt ────────────────────────────────────────
  hr <- read_cache(s$hr)
  if (!is.null(hr)) {
    w <- get_col(hr, "hv005") / 1e6; cl <- get_col(hr, "hv001")
    salt <- miss(get_col(hr, "hv234"), c(994, 995, 996, 997, 998, 999))
    # Two standard readings: any iodine at all, and adequately iodised at 15 ppm.
    parts$salt_any <- agg_to_admin2(as.numeric(salt > 0), w, cl, cmap, "np_salt_iodised_any")
    parts$salt_ade <- agg_to_admin2(as.numeric(salt >= 15), w, cl, cmap, "np_salt_iodised_15ppm")
  }
  # ── child recode: vitamin A capsule, deworming, haemoglobin ───────────────
  kr <- read_cache(s$kr)
  if (!is.null(kr)) {
    w <- get_col(kr, "v005") / 1e6; cl <- get_col(kr, "v001")
    vas <- if ("h34" %in% names(kr)) miss(get_col(kr, "h34"), c(8, 9)) else
      { h <- miss(get_col(kr, "h33"), c(8, 9)); as.numeric(h > 0) }
    parts$vas   <- agg_to_admin2(as.numeric(vas > 0), w, cl, cmap, "np_vita_capsule")
    dw <- miss(get_col(kr, "h43"), c(8, 9))
    parts$dworm <- agg_to_admin2(as.numeric(dw == 1), w, cl, cmap, "np_deworming")
    hb <- miss(get_col(kr, "hw53"), c(999, 9999)) / 10
    hb[hb < 3 | hb > 22] <- NA_real_
    parts$chb   <- agg_to_admin2(hb, w, cl, cmap, "np_child_hb_gdl")
    parts$canem <- agg_to_admin2(as.numeric(hb < 11), w, cl, cmap, "np_child_anaemia_any")
  }
  # ── women recode: haemoglobin, iron in pregnancy ──────────────────────────
  ir <- read_cache(s$ir)
  if (!is.null(ir)) {
    w <- get_col(ir, "v005") / 1e6; cl <- get_col(ir, "v001")
    hb <- miss(get_col(ir, "v453"), c(999, 9999)) / 10
    hb[hb < 3 | hb > 22] <- NA_real_
    parts$whb   <- agg_to_admin2(hb, w, cl, cmap, "np_women_hb_gdl")
    parts$wanem <- agg_to_admin2(as.numeric(hb < 12), w, cl, cmap, "np_women_anaemia_any")
    m45 <- miss(get_col(ir, "m45_1"), c(8, 9))
    parts$iron  <- agg_to_admin2(as.numeric(m45 == 1), w, cl, cmap, "np_iron_pregnancy")
  }
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) { cat("[skip]", cn, "no indicators built\n"); next }
  tab <- Reduce(function(a, b) dplyr::full_join(a, b, by = "key"), parts)
  sp <- strsplit(tab$key, "_", fixed = TRUE)
  tab$Admin1 <- vapply(sp, function(z) z[1], character(1))
  tab$Admin2 <- vapply(sp, function(z) paste(z[-1], collapse = "_"), character(1))
  tab$country <- cn; tab$dhs_round <- s$round
  rows[[cn]] <- tab
  ind <- grep("^np_", names(tab), value = TRUE); ind <- ind[!grepl("_n$", ind)]
  cov_rows[[cn]] <- data.frame(country = cn, dhs_round = s$round,
    districts = nrow(tab),
    indicator = ind,
    pct_districts_nonmissing = round(100 * vapply(ind,
      function(v) mean(is.finite(tab[[v]])), numeric(1)), 1),
    median_value = round(vapply(ind, function(v)
      stats::median(tab[[v]], na.rm = TRUE), numeric(1)), 4),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-13s %-14s districts=%3d indicators=%d\n", cn, s$round, nrow(tab), length(ind)))
}
if (!length(rows)) stop("No country produced nutrition-proximal columns.")
NP <- dplyr::bind_rows(rows)
front <- c("country", "Admin1", "Admin2", "dhs_round")
NP <- NP[, c(front, setdiff(names(NP), c(front, "key")))]
readr::write_csv(NP, file.path(OUT, "nutrition_proximal.csv"))
COV <- dplyr::bind_rows(cov_rows)
readr::write_csv(COV, file.path(OUT, "nutrition_proximal_coverage.csv"))

cat("\n=== WS-F2: coverage by country and indicator ===\n")
print(as.data.frame(COV |> select(country, indicator, pct_districts_nonmissing, median_value) |>
  tidyr::pivot_wider(names_from = country, values_from =
    c(pct_districts_nonmissing, median_value))), row.names = FALSE)

# ── WS-F1: named per-outcome sets over everything now available ─────────────
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
HC <- setdiff(names(H), c("country", "Admin1", "Admin2"))
npc <- grep("^np_", names(NP), value = TRUE); npc <- npc[!grepl("_n$", npc)]
SETS <- list(
  zinc_set   = c(grep("^soil_zinc|^soilgrids_(ph|cec|clay|organic_carbon)", HC, value = TRUE),
                 grep("^spam_share", HC, value = TRUE)),
  iron_set   = c(grep("^np_(child_hb|women_hb|child_anaemia|women_anaemia|iron_pregnancy)", npc, value = TRUE),
                 grep("^dhs_(c_slept_itn|hh_itn_any)", HC, value = TRUE)),
  vita_set   = c(grep("^np_vita_capsule", npc, value = TRUE),
                 grep("^spam_(share|prod)_(oilcrops|vegetables)", HC, value = TRUE)),
  b12_folate_set = c(grep("^dhs_c_(diet_diversity_score|min_diet_diversity|iron_rich_food|vita_rich_food)", HC, value = TRUE),
                     grep("^spam_(share|prod)_pulses", HC, value = TRUE)),
  program_set = c(grep("^np_(salt_iodised|deworming|vita_capsule|iron_pregnancy)", npc, value = TRUE),
                  grep("^dhs_w_anc(1|4plus)$", HC, value = TRUE))
)
SETS <- lapply(SETS, function(v) unique(v[nzchar(v)]))
setdf <- do.call(rbind, lapply(names(SETS), function(k)
  if (!length(SETS[[k]])) NULL else
    data.frame(set = k, column = SETS[[k]],
               source = ifelse(grepl("^np_", SETS[[k]]), "nutrition_proximal", "harmonized_294"),
               stringsAsFactors = FALSE)))
readr::write_csv(setdf, file.path(OUT, "nutrition_proximal_sets.csv"))
cat("\n=== WS-F1: curated sets ===\n")
print(as.data.frame(setdf |> group_by(set, source) |>
  summarise(n = dplyr::n(), .groups = "drop")), row.names = FALSE)

# ── provenance ─────────────────────────────────────────────────────────────
prov <- data.frame(
  product = "nutrition_proximal.csv",
  source = "DHS standard recodes via the local rdhs cache",
  url = "https://dhsprogram.com (accessed through the rdhs cache; no download performed)",
  retrieved = as.character(Sys.Date()),
  version = paste(vapply(SURVEYS, function(z) z$round, character(1)), collapse = "; "),
  license = "DHS Program data use agreement; microdata not redistributed, only Admin-2 aggregates",
  note = "Built offline from an existing cache. No network request was made.",
  stringsAsFactors = FALSE)
pf <- file.path(META, "external_provenance.csv")
if (file.exists(pf)) {
  old <- utils::read.csv(pf, stringsAsFactors = FALSE)
  prov <- dplyr::bind_rows(old[old$product != prov$product, , drop = FALSE], prov)
}
readr::write_csv(prov, pf)
cat(sprintf("\n-> data/covariates/nutrition_proximal.csv (%d districts, %d indicators)\n",
            nrow(NP), length(npc)))
