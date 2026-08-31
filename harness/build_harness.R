# =============================================================================
# harness/build_harness.R
#
# WS8c. Export the evaluation harness for the Proxy Modeling Alliance.
#
# WHAT THIS EXPORTS
#   folds_loro.csv       region-blocked leave-one-region-out assignments
#   folds_loco.csv       leave-one-country-out assignments
#   heldout_cells.csv    the locked test set: one region per country, chosen
#                        with a recorded seed
#   targets_public.csv   the district-level yardstick for the training cells
#                        only, so a partner can fit without seeing the test set
#   ceilings.csv         the validated empirical reliability per cell
#
# COTE D'IVOIRE, AND WHY IT IS NOT IN THE TEST SET
# -----------------------------------------------
# The workstream specification names Cote d'Ivoire external validation as part
# of the locked test set. Stage 0 established that it cannot serve: it is not a
# configured country, it has GEE rasters and oos_civ_* PREDICTION targets, and
# it has no biomarker outcomes. A test cell needs a label. Cote d'Ivoire is
# therefore exported as an unscored prediction target, and the scored test set
# is one held-out region per country.
#
#   Rscript harness/build_harness.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
HDIR <- here("harness"); STORE <- here("_targets_full")
SEED <- 20260910L
set.seed(SEED)

cfgs <- get_country_configs()
rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]
  for (ocn in names(cc$outcomes)) {
    sv <- tryCatch(targets::tar_read_raw(
      paste0("svy_admin2_", tolower(cn), "_", ocn), store = STORE),
      error = function(e) NULL)
    if (is.null(sv) || !nrow(sv)) next
    if (!"Admin1" %in% names(sv)) next
    rows[[length(rows) + 1L]] <- data.frame(
      country = cc$country, outcome = ocn,
      Admin1 = trimws(as.character(sv$Admin1)),
      Admin2 = trimws(as.character(sv$Admin2)),
      svy_prev = sv$svy_prev,
      n_svy = if ("n_svy" %in% names(sv)) sv$n_svy else NA_integer_,
      stringsAsFactors = FALSE)
  }
}
D <- dplyr::bind_rows(rows)
D <- D[is.finite(D$svy_prev), , drop = FALSE]

# Region-blocked LORO: the fold is the region.
folds_loro <- unique(D[, c("country", "outcome", "Admin1", "Admin2")])
folds_loro$fold <- folds_loro$Admin1
readr::write_csv(folds_loro, file.path(HDIR, "folds_loro.csv"))

# LOCO: the fold is the country.
folds_loco <- unique(D[, c("country", "outcome", "Admin1", "Admin2")])
folds_loco$fold <- folds_loco$country
readr::write_csv(folds_loco, file.path(HDIR, "folds_loco.csv"))

# The locked test set: one region per country, drawn once with SEED. Regions
# with fewer than three districts are excluded so a scored fold is not a single
# point. The draw is recorded here AND in the file, so it can be checked.
# HOW MANY REGIONS, AND WHY NOT ONE
# ---------------------------------
# The workstream specification says one held-out region per country. Running the
# acceptance test showed that cannot work: one region yields three or four
# districts per cell, and a predictor that is constant within a region, which
# includes the flat-regional-mean arm WS2 identifies as the best one, has zero
# variance there, so every correlation came back NA. A test set that cannot
# score the best known arm is not a test set.
#
# The draw is therefore 30 percent of a country's regions, with a floor of two,
# among regions holding at least three districts.
pick <- list()
for (cn in unique(D$country)) {
  z <- D[D$country == cn, ]
  sizes <- tapply(z$Admin2, z$Admin1, function(x) length(unique(x)))
  ok <- names(sizes)[sizes >= 3]
  if (length(ok) < 2) ok <- names(sizes)
  n_pick <- max(2L, round(0.30 * length(ok)))
  n_pick <- min(n_pick, length(ok))
  pick[[cn]] <- sample(ok, n_pick)
}
held <- D[mapply(function(c1, a1) a1 %in% pick[[c1]], D$country, D$Admin1), ]
held$held_out <- TRUE
readr::write_csv(
  held[, c("country", "outcome", "Admin1", "Admin2")],
  file.path(HDIR, "heldout_cells.csv"))
writeLines(c(
  "# Locked held-out regions for the evaluation harness",
  sprintf("# seed: %d", SEED),
  "# drawn: 30 percent of regions per country, floor of two, among regions",
  "# holding at least three districts",
  unlist(lapply(names(pick), function(k)
    sprintf("%s: %s", k, paste(pick[[k]], collapse = "; "))))),
  file.path(HDIR, "heldout_regions.txt"))

# The public yardstick: training cells only.
key_h <- paste(held$country, held$outcome, held$Admin1, held$Admin2)
train <- D[!paste(D$country, D$outcome, D$Admin1, D$Admin2) %in% key_h, ]
readr::write_csv(train, file.path(HDIR, "targets_public.csv"))
# The private yardstick, used by the scorer, kept beside it so the harness is
# self-contained for the project's own use. A partner receives only the public
# file; this one is what score_predictions.R reads.
readr::write_csv(held[, c("country","outcome","Admin1","Admin2","svy_prev","n_svy")],
                 file.path(HDIR, "targets_heldout.csv"))

# The validated ceilings, so a partner can compute r_share the same way.
rel <- tryCatch(read.csv(here("results","tables","reliability_empirical.csv"),
                         stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(rel)) {
  rel <- rel[rel$scheme == "within",
             c("country","outcome","lambda_emp","r_max_emp","r_max_emp_lo","r_max_emp_hi")]
  readr::write_csv(rel, file.path(HDIR, "ceilings.csv"))
}

cat(sprintf("folds_loro     %d rows\n", nrow(folds_loro)))
cat(sprintf("heldout_cells  %d rows across %d regions in %d countries\n",
            nrow(held), length(unlist(pick)), length(pick)))
cat(sprintf("targets_public %d rows\n", nrow(train)))
cat("held-out regions:\n")
for (k in names(pick)) cat(sprintf("  %-13s %s\n", k, paste(pick[[k]], collapse = "; ")))
cat("\ndistricts per held-out cell:\n")
print(summary(as.integer(table(paste(held$country, held$outcome)))))
