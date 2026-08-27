# =============================================================================
# scripts/refresh_gee_coverage.R
#
# Bring every configured country up to the SAME Earth Engine covariate coverage.
#
#   Rscript scripts/refresh_gee_coverage.R                    # PLAN ONLY (default)
#   Rscript scripts/refresh_gee_coverage.R --run              # actually extract
#   Rscript scripts/refresh_gee_coverage.R --run --countries Tanzania,Gambia
#   Rscript scripts/refresh_gee_coverage.R --run --force      # ignore the cache
#   Rscript scripts/refresh_gee_coverage.R --run --rebuild    # ALL families, one
#                                                             # vocabulary (implies --force)
#
# WHY THIS EXISTS
# ---------------
# The pooled / LOCO / out-of-sample models intersect covariates by NAME, so a
# country missing a covariate family does not degrade gracefully -- it removes
# that family for EVERY country. Measured 2026-08-24 against the 40-family
# manifest in src/GEE/gee_layer_manifest.R:
#
#     Ghana / Malawi / SierraLeone   36 / 40 families
#     Gambia                         27 / 40
#     Tanzania                        8 / 40
#
# Tanzania was missing 18 families the other four all had, including ALL ELEVEN
# iSDAsoil micronutrient layers (zinc, iron, calcium, magnesium, potassium,
# phosphorus, sulphur, nitrogen, CEC, total carbon, aluminium) -- the layers most
# directly relevant to micronutrient outcomes -- plus NDVI, DailyEVI, LAI8days,
# Productivity, TRMM, WAPOR and AerosolOptical. That is why predict_oos_pooled()
# reported "0 common GEE variables across 5 countries".
#
# Tanzania was originally extracted with ad-hoc scripts
# (sandbox_parsimony/python/fetch_tanzania_rasters*.py) instead of through the
# manifest, which is how it ended up with a different and much smaller layer set.
# This script exists so that cannot happen again: it reads the country list from
# get_country_configs(), so a country added there is picked up automatically with
# no edit here.
#
# WHAT IT DOES NOT DO
# -------------------
# It does not invent layers. Everything comes from GEE_LAYER_MANIFEST, which was
# reconstructed from the analyst's GEE_export_metadata.xlsx. A family absent from
# the manifest is reported and skipped rather than guessed at: a layer with the
# right NAME and the wrong CONTENT would corrupt the pooled model invisibly.
#
# SAFE TO RE-RUN. extract_gee_layers() caches each family to
# data/GEE/_ee_cache/<Country>/<family>.rds, so only missing or previously-failed
# families are pulled. Existing columns in the country CSV are preserved; new
# ones are joined on Admin2.
# =============================================================================

suppressMessages({
  library(here); library(dplyr); library(readr); library(sf)
})

args     <- commandArgs(trailingOnly = TRUE)
DO_RUN     <- "--run"     %in% args
DO_FORCE   <- "--force"   %in% args
# --rebuild re-extracts EVERY family for every selected country, not just the
# missing ones. Needed because the existing CSVs were written by three different
# generations of the extractor and use three incompatible column vocabularies:
#   Ghana/Malawi/SierraLeone  gee_a2_AerosolOptical_2017   (family + year)
#   Gambia                    gee_a2_NDVI_NDVI             (family + band, no year)
#   Tanzania                  gee_a2_NDVI                  (family only)
# Their pairwise column intersection is ZERO, which is why the pooled model
# cannot use them as-is. Topping up only the missing families would leave that
# heterogeneity in place; --rebuild regenerates all of them through the current
# manifest so every country ends up in ONE vocabulary.
DO_REBUILD <- "--rebuild" %in% args
# --rebuild implies --force. Without it the per-family cache would hand back
# the results of the OLD extractor run, so the "rebuild" would faithfully
# reproduce the very vocabulary it is meant to replace.
if (DO_REBUILD) DO_FORCE <- TRUE
pick <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) NULL else strsplit(args[i + 1L], ",")[[1]]
}
ONLY <- pick("--countries")

say <- function(...) cat(sprintf("[gee-coverage] %s\n", sprintf(...)))

source(here("src", "GEE", "gee_layer_manifest.R"))
for (f in list.files(here("R"), pattern = "[.]R$", full.names = TRUE, recursive = TRUE)) {
  try(source(f), silent = TRUE)
}

cfgs <- get_country_configs()

# PREDICTION-ONLY countries: no survey, so they are absent from
# get_country_configs(), but predict_oos_pooled() still needs covariates for
# them in the SAME vocabulary as the training countries. Cote d'Ivoire is the
# out-of-sample target; its existing rasters are dominated by 2017/2018, so the
# analyst treated it as a ~2018 country and that is the anchor used here.
# Add a row here when a new prediction-only target is introduced.
PREDICTION_ONLY <- list(
  "Cote_dIvoire" = list(survey_year = 2018L, gadm_code = "CIV",
                        raster_dir = here("data", "Cote_dIvoire_GEE_rasters"))
)
cfgs <- c(cfgs, PREDICTION_ONLY[setdiff(names(PREDICTION_ONLY), names(cfgs))])
if (!is.null(ONLY)) {
  unknown <- setdiff(ONLY, names(cfgs))
  if (length(unknown)) stop("unknown country/countries: ", paste(unknown, collapse = ", "))
  cfgs <- cfgs[ONLY]
}

FAMILIES <- vapply(GEE_LAYER_MANIFEST, function(l) l$family, character(1))

csv_path <- function(cn, cc) {
  here("data", "GEE", sprintf("%s_%s_admin2_gee.csv", cn, cc$survey_year))
}

# Which manifest families already have columns in a country's zonal CSV?
families_present <- function(path) {
  if (!file.exists(path)) return(character(0))
  h <- names(suppressMessages(read_csv(path, n_max = 0, show_col_types = FALSE)))
  g <- grep("^gee_a2_", h, value = TRUE)
  FAMILIES[vapply(FAMILIES, function(fa) {
    any(grepl(paste0("^gee_a2_", fa, "_"), g, ignore.case = TRUE))
  }, logical(1))]
}

# ── 1. Plan ─────────────────────────────────────────────────────────────────
plan <- lapply(names(cfgs), function(cn) {
  cc <- cfgs[[cn]]
  p  <- csv_path(cn, cc)
  have <- families_present(p)
  need <- if (DO_REBUILD) FAMILIES else setdiff(FAMILIES, have)
  list(country = cn, cc = cc, path = p, have = have, need = need)
})
names(plan) <- names(cfgs)

say("manifest families: %d", length(FAMILIES))
say("%-14s %-6s %-5s %-5s %s", "COUNTRY", "YEAR", "HAVE", "NEED", "CSV")
for (p in plan) {
  say("%-14s %-6s %-5d %-5d %s", p$country, p$cc$survey_year,
      length(p$have), length(p$need),
      if (file.exists(p$path)) basename(p$path) else "(none yet)")
}
for (p in plan) {
  if (length(p$need)) {
    say("")
    say("%s needs %d family(ies):", p$country, length(p$need))
    say("  %s", paste(sort(p$need), collapse = ", "))
  }
}

total_need <- sum(vapply(plan, function(p) length(p$need), integer(1)))
if (total_need == 0) {
  say("")
  say("nothing to do - every country has every family")
  quit(save = "no")
}

if (!DO_RUN) {
  say("")
  say("PLAN ONLY. %d country-family extraction(s) would run.", total_need)
  say("Re-run with --run to execute; add --countries A,B to restrict.")
  say("")
  say("NOTE: topping up missing families alone will NOT make the countries")
  say("commensurable -- the existing CSVs use three different column")
  say("vocabularies whose intersection is zero. Use --run --rebuild to")
  say("regenerate every family for every country in one vocabulary.")
  quit(save = "no")
}

# ── 2. Extract ──────────────────────────────────────────────────────────────
say("")
say("starting extraction (%d country-family pairs)", total_need)
suppressMessages(library(rgee))
rgee::ee_Initialize(drive = FALSE, quiet = TRUE)
say("Earth Engine initialised")

for (p in plan) {
  if (!length(p$need)) { say("%s: complete, skipping", p$country); next }
  cn <- p$country; cc <- p$cc
  say("")
  say("=== %s (%s): %d family(ies) ===", cn, cc$survey_year, length(p$need))

  gadm <- tryCatch(
    geodata::gadm(cc$gadm_code, level = 2, path = here("data", "gadm")),
    error = function(e) {
      say("  GADM failed for %s: %s", cn, conditionMessage(e)); NULL
    })
  if (is.null(gadm)) next

  g <- sf::st_as_sf(gadm)
  g$Admin2 <- g$NAME_2
  # Same key hygiene the rest of the pipeline applies, so the join below cannot
  # fan out on duplicate names or carry water-body "districts".
  g <- clean_admin2_keys(g, what = paste(cn, "GADM"))

  sub <- Filter(function(l) l$family %in% p$need, GEE_LAYER_MANIFEST)
  new <- tryCatch(
    extract_gee_layers(g, year = as.integer(cc$survey_year), manifest = sub,
                       id_col = "Admin2", col_prefix = "gee_a2_",
                       cache_dir = here("data", "GEE", "_ee_cache", cn),
                       force = DO_FORCE),
    error = function(e) {
      say("  extraction failed for %s: %s", cn, conditionMessage(e)); NULL
    })
  if (is.null(new) || !nrow(new)) { say("  nothing returned for %s", cn); next }

  n_before <- if (file.exists(p$path)) {
    length(names(suppressMessages(read_csv(p$path, n_max = 0, show_col_types = FALSE))))
  } else 0L

  if (DO_REBUILD) {
    # Replace outright. Joining onto the old file would carry the stale
    # vocabulary forward alongside the new one, which is the very problem
    # --rebuild exists to remove. Keep a timestamped copy first.
    if (file.exists(p$path)) {
      bak <- sub("[.]csv$", sprintf("_pre-rebuild_%s.csv", format(Sys.Date(), "%Y%m%d")), p$path)
      file.copy(p$path, bak, overwrite = TRUE)
      say("  kept previous CSV as %s", basename(bak))
    }
    out <- new
  } else if (file.exists(p$path)) {
    old <- suppressMessages(read_csv(p$path, show_col_types = FALSE))
    # Never silently overwrite an existing column: drop the incoming duplicate
    # and report it, so a changed upstream asset is visible rather than absorbed.
    dup <- setdiff(intersect(names(old), names(new)), "Admin2")
    if (length(dup)) {
      say("  %d incoming column(s) already present, keeping existing: %s",
          length(dup), paste(utils::head(dup, 5), collapse = ", "))
      new <- new[, setdiff(names(new), dup), drop = FALSE]
    }
    out <- dplyr::left_join(old, new, by = "Admin2")
  } else {
    out <- new
  }
  write_csv(out, p$path)
  say("  wrote %s: %d cols (was %d)", basename(p$path), ncol(out), n_before)

  still <- setdiff(FAMILIES, families_present(p$path))
  if (length(still)) {
    say("  STILL MISSING (%d): %s", length(still), paste(sort(still), collapse = ", "))
  } else {
    say("  %s now complete", cn)
  }
}

say("")
say("done. Re-run without --run for a fresh coverage report,")
say("then rebuild the pipeline so the new covariates enter the models.")
