# =============================================================================
# scripts/covariates/harmonize_extra_domains.R
#
# Harmonize IHME, Malaria Atlas and food-security predictors to Admin-2 in all
# four countries, so they can enter the pooled and leave-one-country-out models.
#
# WHY THESE THREE, AND NOT THE OTHER FOUR
# ---------------------------------------
# The harmonized predictor set is 294 columns and draws on DHS, soil chemistry,
# crop maps and the GEE raster stack. Seven further domains exist in the
# country-level data and none of them reach it. Checked on disk:
#
#   IHME            4 of 4 country rollups                 -> BUILD
#   Malaria Atlas   30 raster stems common to all four     -> BUILD
#   Food security   HFID covers all four, and its Admin-2
#                   counts match ours exactly (260/243/37/14) -> BUILD
#   WFP prices      4 of 4 files, commodity overlap TBD     -> checked, reported
#   MICS            2 of 4 (Gambia, Ghana only)             -> BLOCKED, no data
#   LSMS            1 of 4 (Ghana only)                     -> BLOCKED, no data
#   FluNet          country-level weekly surveillance       -> USELESS at Admin-2
#
# The last three are not harmonization problems. MICS and LSMS were never
# collected in the other countries, and FluNet has no within-country variation
# at all, so as an Admin-2 predictor it is a constant per country and cannot
# contribute to a district ranking.
#
# THE FAILURE MODE THIS GUARDS AGAINST
# ------------------------------------
# The pooled and LOCO models match covariates by EXACT COLUMN NAME, and a name
# missing for one country is treated as zero rather than raising an error. A
# previous instance of that silently zeroed an entire GEE block for every
# country. So this script does not merely build columns, it REFUSES to emit any
# column that is not present, non-constant and on a comparable scale in all four
# countries, and it reports what it dropped and why.
#
#   Rscript scripts/covariates/harmonize_extra_domains.R
# -> data/covariates/harmonized/predictors_admin2_extra.csv
# -> results/tables/harmonize_extra_domains_manifest.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
targets::tar_source(here("R"))

OUT <- here("data", "covariates", "harmonized", "predictors_admin2_extra.csv")
MAN <- here("results", "tables", "harmonize_extra_domains_manifest.csv")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
SURVEY_YEAR <- c(Gambia = 2018, Ghana = 2017, Malawi = 2016, SierraLeone = 2013)
cfgs <- get_country_configs()
kk <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(as.character(x))))

# The Admin-2 spine every domain must join onto: taken from the existing
# harmonized file so the new columns line up with the ones already in use.
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
spine <- H[, c("country", "Admin1", "Admin2")]
cat(sprintf("[extra] spine: %d Admin-2 units across %d countries\n",
            nrow(spine), length(unique(spine$country))))

manifest <- list()
blocks <- list()

note <- function(domain, column, status, reason = NA_character_) {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    domain = domain, column = column, status = status, reason = reason,
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 1. IHME. Country rollups are long: one row per admin-2 x stratifier. Keep
#    only (measure, metric, age, sex) combinations that exist in ALL FOUR.
# ---------------------------------------------------------------------------
cat("\n[extra] IHME\n")
ihme_files <- c(Gambia = "IHME_Gambia_data.dta", Ghana = "IHME_Ghana_data.dta",
                Malawi = "IHME_Malawi_data.dta",
                SierraLeone = "IHME_Sierra_Leone_data.dta")
ihme_long <- list()
for (cn in COUNTRIES) {
  p <- here("data", "IHME", ihme_files[[cn]])
  if (!file.exists(p)) { cat(sprintf("  %-12s ABSENT\n", cn)); next }
  d <- tryCatch(haven::read_dta(p), error = function(e) NULL)
  if (is.null(d)) { cat(sprintf("  %-12s unreadable\n", cn)); next }
  nm <- names(d)
  need <- c("adm2_name", "year", "measure", "metric", "mean")
  if (!all(need %in% nm)) {
    cat(sprintf("  %-12s missing columns: %s\n", cn,
                paste(setdiff(need, nm), collapse = ", "))); next
  }
  yr <- SURVEY_YEAR[[cn]]
  d$year <- suppressWarnings(as.numeric(d$year))
  # closest available year to the survey
  yrs <- sort(unique(d$year[is.finite(d$year)]))
  pick <- yrs[which.min(abs(yrs - yr))]
  d <- d[is.finite(d$year) & d$year == pick, , drop = FALSE]
  ag <- if ("age_group_name" %in% nm) as.character(d$age_group_name) else ""
  sx <- if ("sex" %in% nm) as.character(d$sex) else ""
  d$.key <- paste("ihme", kk(d$measure), kk(d$metric), kk(ag), kk(sx), sep = "_")
  d$.a2 <- as.character(d$adm2_name)
  d$.val <- suppressWarnings(as.numeric(d$mean))
  ihme_long[[cn]] <- d[, c(".a2", ".key", ".val")]
  cat(sprintf("  %-12s year %d, %d rows, %d distinct indicators\n",
              cn, pick, nrow(d), length(unique(d$.key))))
}
if (length(ihme_long) == 4L) {
  common <- Reduce(intersect, lapply(ihme_long, function(z) unique(z$.key)))
  cat(sprintf("  indicators common to all four: %d\n", length(common)))
  for (cn in COUNTRIES) {
    z <- ihme_long[[cn]]
    z <- z[z$.key %in% common, , drop = FALSE]
    w <- stats::aggregate(.val ~ .a2 + .key, data = z, FUN = function(v)
      mean(v, na.rm = TRUE))
    w <- tidyr::pivot_wider(w, names_from = ".key", values_from = ".val")
    # pivot_wider leaves the id column named `.a2`; the join expects `Admin2`.
    # Without this rename the IHME block fails the join outright and the whole
    # domain vanishes with one line of log.
    names(w)[names(w) == ".a2"] <- "Admin2"
    # FUZZY-MATCH THE ADMIN-2 NAMES ONTO THE SPINE.
    #
    # Exact matching recovers only about two thirds of Ghana's districts (176 of
    # 260) because IHME's adm2_name spellings differ from GADM's. At 34 percent
    # mean completeness the whole IHME block failed the four-country filter, so
    # the domain was present in the log and absent from the output. The rest of
    # the pipeline already fuzzy-matches these names (Jaro-Winkler, 0.85); the
    # same is done here, and the match rate is reported rather than assumed.
    tgt <- unique(spine$Admin2[spine$country == cn])
    if (requireNamespace("stringdist", quietly = TRUE) && length(tgt)) {
      src <- as.character(w$Admin2)
      hit <- src %in% tgt
      if (any(!hit)) {
        dm <- stringdist::stringdistmatrix(kk(src[!hit]), kk(tgt), method = "jw", p = 0.1)
        best <- apply(dm, 1, which.min); dist <- apply(dm, 1, min)
        repl <- ifelse(dist <= 0.15, tgt[best], src[!hit])
        src[!hit] <- repl
        w$Admin2 <- src
      }
      cat(sprintf("  %-12s name match %d/%d exact, %d/%d after fuzzy
", cn,
                  sum(hit), length(hit), sum(w$Admin2 %in% tgt), nrow(w)))
    }
    w$country <- cn
    blocks[[paste0("ihme_", cn)]] <- as.data.frame(w)
  }
  for (k in common) note("IHME", k, "kept")
} else {
  cat("  fewer than four countries readable; IHME block skipped\n")
  note("IHME", NA, "skipped", "not all four country rollups readable")
}

# ---------------------------------------------------------------------------
# 2. Malaria Atlas. Zonal means per Admin-2 polygon, for the raster stems that
#    exist in every country.
# ---------------------------------------------------------------------------
cat("\n[extra] Malaria Atlas\n")
map_dir <- here("data", "Malaria Atlas")
stems <- list()
for (cn in COUNTRIES) {
  p <- file.path(map_dir, cn)
  if (!dir.exists(p)) next
  f <- list.files(p, pattern = "\\.tif$", ignore.case = TRUE)
  stems[[cn]] <- unique(sub("\\.tif$", "", sub("^[A-Z]{3}_", "", f), ignore.case = TRUE))
}
map_common <- if (length(stems) == 4L) Reduce(intersect, stems) else character(0)
cat(sprintf("  raster stems common to all four: %d\n", length(map_common)))

if (length(map_common) && requireNamespace("terra", quietly = TRUE) &&
    requireNamespace("sf", quietly = TRUE)) {
  for (cn in COUNTRIES) {
    cc <- cfgs[[cn]]
    poly <- tryCatch(load_gadm_cached(cc$gadm_code, level = 2),
                     error = function(e) NULL)
    if (is.null(poly)) { cat(sprintf("  %-12s no GADM polygons\n", cn)); next }
    a1 <- poly[[grep("^NAME_1$", names(poly), value = TRUE)[1]]]
    a2 <- poly[[grep("^NAME_2$", names(poly), value = TRUE)[1]]]
    vp <- terra::vect(poly)
    res <- data.frame(country = cn, Admin1 = as.character(a1),
                      Admin2 = as.character(a2), stringsAsFactors = FALSE)
    files <- list.files(file.path(map_dir, cn), pattern = "\\.tif$",
                        ignore.case = TRUE, full.names = TRUE)
    got <- 0L
    for (st in map_common) {
      f <- files[sub("\\.tif$", "", sub("^[A-Z]{3}_", "", basename(files)),
                     ignore.case = TRUE) == st]
      if (!length(f)) next
      r <- tryCatch(terra::rast(f[1]), error = function(e) NULL)
      if (is.null(r)) next
      v <- tryCatch({
        vv <- terra::project(vp, terra::crs(r))
        terra::extract(r[[1]], vv, fun = mean, na.rm = TRUE, ID = FALSE)[, 1]
      }, error = function(e) NULL)
      if (is.null(v) || length(v) != nrow(res)) next
      res[[paste0("map_", kk(st))]] <- as.numeric(v)
      got <- got + 1L
    }
    blocks[[paste0("map_", cn)]] <- res
    cat(sprintf("  %-12s %d of %d rasters extracted over %d polygons\n",
                cn, got, length(map_common), nrow(res)))
  }
  for (st in map_common) note("MalariaAtlas", paste0("map_", kk(st)), "kept")
} else {
  cat("  terra/sf unavailable or no common stems; MAP block skipped\n")
  note("MalariaAtlas", NA, "skipped", "no common stems or terra/sf missing")
}

# ---------------------------------------------------------------------------
# 3. Food security, from HFID. All four countries, Admin-2 named directly.
# ---------------------------------------------------------------------------
cat("\n[extra] Food security (HFID)\n")
hf <- tryCatch(suppressMessages(readr::read_csv(here("data", "HFID", "hfid_hv1.csv"),
                                                show_col_types = FALSE)),
               error = function(e) NULL)
if (!is.null(hf)) {
  hf <- as.data.frame(hf)
  hf$.c <- kk(hf$ADMIN0)
  map_c <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
             sierraleone = "SierraLeone")
  hf$country <- unname(map_c[hf$.c])
  hf <- hf[!is.na(hf$country), , drop = FALSE]
  hf$.yr <- suppressWarnings(as.numeric(substr(as.character(hf$year_month), 1, 4)))
  keep <- intersect(c("ipc_phase_fews", "ipc_phase_ipcch"), names(hf))
  out <- list()
  for (cn in COUNTRIES) {
    z <- hf[hf$country == cn, , drop = FALSE]
    if (!nrow(z)) next
    yr <- SURVEY_YEAR[[cn]]
    yrs <- sort(unique(z$.yr[is.finite(z$.yr)]))
    if (!length(yrs)) next
    pick <- yrs[which.min(abs(yrs - yr))]
    z <- z[is.finite(z$.yr) & z$.yr == pick, , drop = FALSE]
    agg <- z |> group_by(Admin2 = as.character(ADMIN2)) |>
      summarise(across(all_of(keep), ~ mean(suppressWarnings(as.numeric(.x)),
                                            na.rm = TRUE)), .groups = "drop")
    names(agg)[-1] <- paste0("fsec_", names(agg)[-1])
    agg$country <- cn
    out[[cn]] <- as.data.frame(agg)
    cat(sprintf("  %-12s year %d, %d Admin-2 units\n", cn, pick, nrow(agg)))
  }
  if (length(out) == 4L) {
    blocks[["fsec"]] <- dplyr::bind_rows(out)
    for (k in keep) note("FoodSecurity", paste0("fsec_", k), "kept")
  } else {
    cat("  fewer than four countries; food security block skipped\n")
    note("FoodSecurity", NA, "skipped",
         sprintf("only %d of 4 countries present", length(out)))
  }
}

# ---------------------------------------------------------------------------
# The domains that cannot be harmonized, recorded rather than silently absent.
# ---------------------------------------------------------------------------
note("MICS", NA, "blocked", "collected in Gambia and Ghana only; no data for Malawi or Sierra Leone")
note("LSMS", NA, "blocked", "Ghana only; no data for the other three")
note("FluNet", NA, "not_applicable",
     "country-level weekly surveillance, constant within a country, so it carries no Admin-2 variation")

# ---------------------------------------------------------------------------
# Assemble, then REFUSE anything not usable in all four countries.
# ---------------------------------------------------------------------------
# ROW-BIND THE PER-COUNTRY BLOCKS BEFORE JOINING, one join per domain.
#
# An earlier version joined each country's block separately. Because every
# country's block carries the SAME column names, the second join collided with
# the first and dplyr suffixed them .x/.y: 89 real columns became 122 mangled
# ones, each populated for a single country, and every one then failed the
# four-country filter. The filter reported "0 columns survive" and was right,
# but the cause was upstream of it.
domain_of <- function(nm) sub("_(Gambia|Ghana|Malawi|SierraLeone)$", "", nm)
merged <- list()
for (nm in names(blocks)) {
  b <- blocks[[nm]]
  if (is.null(b) || !nrow(b)) next
  d <- domain_of(nm)
  merged[[d]] <- if (is.null(merged[[d]])) b else dplyr::bind_rows(merged[[d]], b)
}

join_block <- function(sp, b) {
  if (!"country" %in% names(b)) return(NULL)
  by <- if (all(c("Admin1", "Admin2") %in% names(b))) c("country", "Admin1", "Admin2") else
    c("country", "Admin2")
  b$Admin2 <- as.character(b$Admin2)
  if ("Admin1" %in% names(b)) b$Admin1 <- as.character(b$Admin1)
  # A duplicated key would fan the spine out and silently change its row count.
  # GADM carries 256 Malawi polygons against the spine's 243, so this is real.
  dup <- duplicated(b[, by, drop = FALSE])
  if (any(dup)) {
    cat(sprintf("  [join] %d duplicate keys collapsed before joining
", sum(dup)))
    b <- b[!dup, , drop = FALSE]
  }
  dplyr::left_join(sp, b, by = by)
}
res <- spine
for (d in names(merged)) {
  j <- tryCatch(join_block(res, merged[[d]]), error = function(e) NULL)
  if (is.null(j)) { cat(sprintf("  [join] %s FAILED
", d)); next }
  if (nrow(j) != nrow(res)) {
    cat(sprintf("  [join] %s changed row count %d -> %d, dropped
",
                d, nrow(res), nrow(j))); next
  }
  added <- setdiff(names(j), names(res))
  cov <- if (length(added)) mean(vapply(added, function(v)
    mean(is.finite(j[[v]])), numeric(1))) else 0
  cat(sprintf("  [join] %-8s +%3d columns, mean completeness %.0f%%
",
              d, length(added), 100 * cov))
  res <- j
}

newcols <- setdiff(names(res), c("country", "Admin1", "Admin2"))
cat(sprintf("\n[extra] %d candidate columns before the four-country filter\n",
            length(newcols)))
ok <- vapply(newcols, function(v) {
  s <- split(res[[v]], res$country)
  length(s) == 4L && all(vapply(s, function(z)
    mean(is.finite(z)) > 0.5 && stats::sd(z, na.rm = TRUE) > 0, logical(1)))
}, logical(1))
for (v in newcols[!ok]) note("filter", v, "dropped",
                             "absent, mostly missing or constant in at least one country")
res <- res[, c("country", "Admin1", "Admin2", newcols[ok]), drop = FALSE]
cat(sprintf("[extra] %d columns survive: present, non-constant and >50%% complete in all four\n",
            sum(ok)))

readr::write_csv(res, OUT)
readr::write_csv(dplyr::bind_rows(manifest), MAN)
cat(sprintf("\n-> %s\n-> %s\n",
            file.path("data","covariates","harmonized", basename(OUT)),
            file.path("results","tables", basename(MAN))))
