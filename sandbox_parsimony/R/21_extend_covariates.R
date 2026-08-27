# =============================================================================
# sandbox_parsimony/R/21_extend_covariates.R
#
# Close the diet-quality gap in the predictor set.
#
# FINDINGS.md section 4 showed the five-country covariate intersection has no
# travel time to cities, no land surface temperature, no climate suite, no
# cropland fraction and no vegetation dynamics -- and section 11 showed the one
# vegetation variable that does survive (NDVI) ranks dead last. This script
# rebuilds the Admin-2 covariate table with those pathways represented.
#
# Four sources, in order of how much they add:
#
#  1. RICH GEE FAMILIES ALREADY ON DISK. area_covariates_<country>$gee_admin2
#     holds accessibility, lst, terraclimate, productivity, dailyevi, lai8days,
#     landcoverlayers and esa worldcereal for Gambia, Ghana, Sierra Leone and
#     Malawi. They vanish from the pooled model only because the name
#     intersection includes Tanzania, which lacks them. Nothing to download.
#
#  2. MapSPAM CROP MIX (data/MapSPAM/*_spam_admin2.csv) -- production and area
#     share of cereals / roots / pulses / oilcrops / vegetables. Present for ALL
#     FIVE countries. This is the most direct measurement of what is grown and
#     therefore plausibly eaten that the project holds, and it is currently
#     unused by the area-level model.
#
#  3. KOPPEN-GEIGER climate class + AGRO-ECOLOGICAL ZONE, extracted from the
#     global rasters already in data/. Categorical agro-ecology, available for
#     every country including Tanzania. Encoded as class-share fractions per
#     district rather than a modal code, since a modal integer is meaningless
#     as a numeric predictor.
#
#  4. TRAVEL TIME TO CITIES for all five countries from one consistent global
#     product (see 23_travel_time.R). This is the variable with the strongest
#     prior claim on dietary diversity and it was absent everywhere. Optional;
#     this script uses it when out/travel_time_admin2.csv is present.
#
# Output: sandbox_parsimony/out/cov_ext.rds, same shape as out/cov_list.rds so
# 00_assemble.R::assemble_outcome() can consume it unchanged.
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
suppressPackageStartupMessages({library(dplyr); library(sf); library(terra)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi", "tanzania")
SPAM_FILE <- c(gambia = "Gambia", ghana = "Ghana", sierraleone = "SierraLeone",
               malawi = "Malawi", tanzania = "Tanzania")

# The families FINDINGS section 4 identified as missing from the intersection.
RICH_FAMILIES <- paste0("^gee_(accessibility|lst|terraclimate|productivity|",
                        "dailyevi|lai8days|landcoverlayers|landcovertype|esa|",
                        "wapor|ghsbuilts|ghspop|fldas)")

# ---------------------------------------------------------------------------
# 0. Admin-2 key hygiene
#
# Two problems that bite any join keyed on the Admin-2 NAME:
#
#  (a) GADM ships inland water bodies as Admin-2 polygons. Malawi has 4 and
#      Tanzania 5, and they are split into multiple pieces so the name repeats.
#      Extracting NDVI, soil chemistry or LST over open water produces
#      meaningless covariates, and these are not populated prediction units.
#      WORSE: "Lake Malawi" (n_svy 21-33) and "Lake Victoria" (n_svy 60-110)
#      appear in the SURVEY tables, so respondents have been geocoded onto a
#      lake. Lake Victoria is one of Tanzania's largest surveyed units.
#      Reported in FINDINGS.md; dropped here.
#
#  (b) Genuine repeated names across Admin-1 regions (Malawi has TA Lundu,
#      TA Malemia, TA Ngabu, TA Pemba). Joining on the name alone then produces
#      a cartesian blow-up -- unfixed, this script inflated Malawi from 256 rows
#      to 29,183. Since the survey side aggregates by NAME too, the consistent
#      resolution is to collapse same-named polygons to one row.
# ---------------------------------------------------------------------------
WATER_PATTERN <- "^Lake |^Water|Reservoir|^Lac |^Sea$"

drop_water <- function(d) {
  if (!"Admin2" %in% names(d)) return(d)
  w <- grepl(WATER_PATTERN, d$Admin2, ignore.case = TRUE)
  if (any(w)) message(sprintf("    dropped %d water-body polygon(s): %s",
                              sum(w), paste(unique(d$Admin2[w]), collapse = ", ")))
  d[!w, , drop = FALSE]
}

#' Collapse repeated Admin-2 names to one row (numeric mean, first non-numeric)
dedupe_admin2 <- function(d) {
  if (!"Admin2" %in% names(d) || !any(duplicated(d$Admin2))) return(d)
  n_dup <- sum(duplicated(d$Admin2))
  out <- d |> group_by(Admin2) |>
    summarise(across(everything(),
                     ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else dplyr::first(.x)),
              .groups = "drop") |> as.data.frame()
  message(sprintf("    collapsed %d duplicate Admin-2 name row(s)", n_dup))
  out
}

# ---------------------------------------------------------------------------
# 1. Categorical agro-ecology as class-share fractions
# ---------------------------------------------------------------------------
#' @param r a categorical SpatRaster
#' @param pg sf polygons with an Admin2 column
#' @param prefix column-name prefix
#' @param min_share drop classes that never exceed this share anywhere
class_shares <- function(r, pg, prefix, min_share = 0.05) {
  pgv <- terra::vect(sf::st_transform(pg, terra::crs(r)))
  ex <- terra::extract(r, pgv)
  names(ex)[2] <- "cls"
  tab <- ex |>
    filter(!is.na(cls)) |>
    count(ID, cls) |>
    group_by(ID) |>
    mutate(share = n / sum(n)) |>
    ungroup()
  if (!nrow(tab)) return(NULL)
  keep <- tab |> group_by(cls) |> summarise(mx = max(share), .groups = "drop") |>
    filter(mx >= min_share) |> pull(cls)
  tab <- tab |> filter(cls %in% keep)
  if (!nrow(tab)) return(NULL)
  w <- tab |> select(ID, cls, share) |>
    tidyr::pivot_wider(names_from = cls, values_from = share, values_fill = 0)
  out <- data.frame(Admin2 = as.character(pg$Admin2)[w$ID])
  for (cl in setdiff(names(w), "ID"))
    out[[paste0(prefix, "_c", cl)]] <- w[[cl]]
  # also a diversity index: how agro-ecologically mixed is the district
  sh <- as.matrix(w[, setdiff(names(w), "ID"), drop = FALSE])
  out[[paste0(prefix, "_shannon")]] <-
    apply(sh, 1, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
  out
}

get_aez <- function() {
  z <- "data/Agro-Ecological Zones for Africa South of the Sahara/AEZ8 r2.0 - TIF.zip"
  if (!file.exists(z)) return(NULL)
  dst <- file.path(tempdir(), "aez8")
  if (!dir.exists(dst)) dir.create(dst, recursive = TRUE)
  tifs <- list.files(dst, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (!length(tifs)) {
    tryCatch(utils::unzip(z, exdir = dst), error = function(e) NULL)
    tifs <- list.files(dst, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  }
  if (!length(tifs)) return(NULL)
  tryCatch(terra::rast(tifs[1]), error = function(e) NULL)
}

koppen <- tryCatch(
  terra::rast("data/Koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"),
  error = function(e) NULL)
aez <- get_aez()
message(sprintf("koppen: %s | aez: %s", !is.null(koppen), !is.null(aez)))

# ---------------------------------------------------------------------------
# 2. MapSPAM crop mix
# ---------------------------------------------------------------------------
read_spam <- function(ctry) {
  f <- file.path("data/MapSPAM", paste0(SPAM_FILE[[ctry]], "_spam_admin2.csv"))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !"Admin2" %in% names(d)) return(NULL)
  num <- setdiff(names(d), c("Admin1", "Admin2"))
  out <- data.frame(Admin2 = as.character(d$Admin2), stringsAsFactors = FALSE)
  for (v in num) out[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  # production totals scale with district size; the SHARES are the dietary
  # signal, so add a per-area intensity term rather than leaving raw tonnage
  # to act as a population proxy.
  if (all(c("spam_prod_total", "spam_parea_total") %in% names(out)))
    out$spam_yield_proxy <- out$spam_prod_total / pmax(out$spam_parea_total, 1e-6)
  out
}

# ---------------------------------------------------------------------------
# Build the extended table per country
# ---------------------------------------------------------------------------
travel <- if (file.exists("sandbox_parsimony/out/travel_time_admin2.csv"))
  read.csv("sandbox_parsimony/out/travel_time_admin2.csv", stringsAsFactors = FALSE) else NULL
if (is.null(travel)) message("NOTE: travel_time_admin2.csv absent - run 23_travel_time.R first")

cov_ext <- list()
for (ctry in COUNTRIES) {
  base <- build_cov(ctry)            # the table 00_assemble.R already produces
  if (is.null(base)) next
  base <- dedupe_admin2(drop_water(base))
  n_base <- ncol(base) - 1

  ac <- .rd(paste0("area_covariates_", ctry))
  gee <- if (!is.null(ac)) ac$gee_admin2 else NULL
  pg  <- if (!is.null(ac)) ac$polygons else NULL

  # (1) rich GEE families, year-harmonised the same way build_cov() does
  if (!is.null(gee)) {
    rich <- grep(RICH_FAMILIES, names(gee), value = TRUE)
    if (length(rich)) {
      stems <- harmonize_names(rich)
      gm <- vapply(rich, function(c) .num(gee[[c]]), numeric(nrow(gee)))
      us <- unique(stems)
      rdf <- as.data.frame(vapply(us, function(s) {
        k <- which(stems == s)
        if (length(k) == 1) gm[, k] else rowMeans(gm[, k, drop = FALSE], na.rm = TRUE)
      }, numeric(nrow(gee))))
      names(rdf) <- us
      rdf$Admin2 <- as.character(gee$Admin2)
      rdf <- dedupe_admin2(drop_water(rdf[, !duplicated(names(rdf))]))
      base <- dplyr::left_join(base, rdf, by = "Admin2")
    }
  }

  # (2) crop mix
  sp <- read_spam(ctry)
  if (!is.null(sp)) base <- dplyr::left_join(base, dedupe_admin2(sp), by = "Admin2")

  # (3) agro-ecology class shares
  if (!is.null(pg) && "Admin2" %in% names(pg)) {
    if (!is.null(koppen)) {
      kk <- tryCatch(class_shares(koppen, pg, "agro_koppen"), error = function(e) NULL)
      if (!is.null(kk)) base <- dplyr::left_join(base, dedupe_admin2(kk), by = "Admin2")
    }
    if (!is.null(aez)) {
      aa <- tryCatch(class_shares(aez, pg, "agro_aez"), error = function(e) NULL)
      if (!is.null(aa)) base <- dplyr::left_join(base, dedupe_admin2(aa), by = "Admin2")
    }
  }

  # (4) travel time to cities, same product for every country
  if (!is.null(travel)) {
    tv <- travel[travel$country == ctry, setdiff(names(travel), "country"), drop = FALSE]
    if (nrow(tv)) base <- dplyr::left_join(base, dedupe_admin2(tv), by = "Admin2")
  }

  # carry Admin1 through so 17_admin1_layer.R can aggregate
  if (!is.null(pg) && !"Admin1" %in% names(base)) {
    a1 <- sf::st_drop_geometry(pg)
    nm1 <- if ("Admin1" %in% names(a1)) "Admin1" else if ("NAME_1" %in% names(a1)) "NAME_1" else NULL
    if (!is.null(nm1))
      base <- dplyr::left_join(
        base, data.frame(Admin2 = as.character(pg$Admin2),
                         Admin1 = as.character(a1[[nm1]]),
                         stringsAsFactors = FALSE) |> distinct(Admin2, .keep_all = TRUE),
        by = "Admin2")
  }
  # Hard guard: every join above is keyed on the Admin-2 name, so a duplicate
  # here means a silent cartesian blow-up downstream. Fail loudly instead.
  stopifnot(!any(duplicated(base$Admin2)))

  base <- base[, !duplicated(names(base))]
  cov_ext[[ctry]] <- base
  message(sprintf("  %-12s %4d areas: %4d -> %4d covariates (+%d)",
                  ctry, nrow(base), n_base, ncol(base) - 1, ncol(base) - 1 - n_base))
}

saveRDS(cov_ext, "sandbox_parsimony/out/cov_ext.rds")

cat("\n=== what the extension adds to the FIVE-country intersection ===\n")
nm <- lapply(cov_ext, function(d) setdiff(names(d), c("Admin2", "Admin1")))
old <- lapply(readRDS("sandbox_parsimony/out/cov_list.rds"),
              function(d) setdiff(names(d), c("Admin2", "Admin1")))
cat(sprintf("before: %d common covariates\nafter : %d common covariates\n",
            length(Reduce(intersect, old)), length(Reduce(intersect, nm))))
newv <- setdiff(Reduce(intersect, nm), Reduce(intersect, old))
cat(sprintf("\n%d newly common variables:\n", length(newv)))
print(newv)
message("\nSaved -> sandbox_parsimony/out/cov_ext.rds")
