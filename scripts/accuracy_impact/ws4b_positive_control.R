# =============================================================================
# scripts/accuracy_impact/ws4b_positive_control.R
#
# WS4b, supplement. The seven named pseudo-targets, including the education
# positive control that Stage 0 could not source.
#
# WHY THIS IS SEPARATE FROM ws4b_skill_curve.R
# --------------------------------------------
# The main WS4b script needs a reliability per target, so it uses the 21 DHS
# indicators that carry a design-based variance in *_dhs_admin2_direct.rds.
# Education, improved water and the wealth index are NOT among those 21. They
# exist only in *_dhs_custom_admin2_wide.rds, as point estimates with no
# standard error, so no reliability can be computed for them.
#
# Their achieved skill can still be measured, and for claim 3.12 that is the
# part that matters. Section 11 claim 2 rests on a positive control stating that
# earth observation predicts district education at r 0.48 to 0.71. Stage 0
# searched results/, sensitivity/, sandbox_parsimony/ and docs/ and found no
# committed table containing it. This script measures it.
#
# A RESOLUTION CAVEAT THAT MATTERS FOR READING THESE NUMBERS
# ---------------------------------------------------------
# The DHS custom Admin-2 files use each country's DHS district system, which is
# NOT the analytic system the micronutrient cells use: Ghana returns 219 units
# here against the 75 analytic districts, and Gambia 35 against 30. More units,
# each measured on a different sample, is a different estimation problem. These
# correlations are therefore evidence about whether the covariates can track a
# well-measured district quantity AT ALL, and are not directly comparable with
# the micronutrient cells' correlations.
#
# SAME EXCLUSION AS THE MAIN SCRIPT. The 97 DHS-derived covariates are removed.
# Predicting dhs2014_w_no_education_adm2 from a covariate set containing DHS
# education aggregates would be the outcome predicting itself, which is what the
# first run of the main script did before the exclusion was added.
#
#   Rscript scripts/accuracy_impact/ws4b_positive_control.R
# -> results/tables/positive_control_targets.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

SEED <- 20260909L; K_SCREEN <- 20L; MIN_D <- 12L
set.seed(SEED)
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- grep("^dhs", setdiff(names(H), c("country","Admin1","Admin2")),
             value = TRUE, invert = TRUE)

# The seven pseudo-targets the workstream names, mapped onto the custom wide
# file's columns. The year prefix differs by country, so the stem is matched.
STEMS <- c(
  education_none      = "w_no_education_adm2",
  education_secondary = "w_secondary_plus_adm2",
  anaemia_child       = "c_anemia_any_adm2",
  haemoglobin_child   = "c_mean_hemoglobin_adm2",
  stunting            = "c_stunted_adm2",
  wealth_index        = "hh_wealth_mean_adm2",
  vaccination_full    = "c_fully_vaccinated_adm2",
  water_improved      = "hh_improved_water_adm2",
  sanitation_improved = "hh_improved_sanitation_adm2"
)
FILES <- list(Gambia = "Gambia_2019", Ghana = "Ghana_2014",
              Malawi = "Malawi_2015", `Sierra Leone` = "Sierra Leone_2013")

rows <- list()
for (cn in names(FILES)) {
  p <- here("data", "DHS", "clean", paste0(FILES[[cn]], "_dhs_custom_admin2_wide.rds"))
  if (!file.exists(p)) { cat("[skip]", cn, "\n"); next }
  w <- readRDS(p)
  if (!"Admin2" %in% names(w)) next
  hc <- H[H$country == gsub(" ", "", cn) | H$country == cn, , drop = FALSE]
  if (!nrow(hc)) next
  for (nm in names(STEMS)) {
    col <- grep(paste0(STEMS[[nm]], "$"), names(w), value = TRUE)
    if (!length(col)) next
    z <- data.frame(Admin2 = trimws(as.character(w$Admin2)),
                    y = suppressWarnings(as.numeric(w[[col[1]]])),
                    stringsAsFactors = FALSE)
    z <- z[is.finite(z$y), , drop = FALSE]
    if (nrow(z) < MIN_D) next
    # Name-only is the only key available: the custom wide file carries no
    # Admin1 column, so admin2_join_by() takes its documented fallback. Recorded
    # rather than silently accepted, because this is the defect class Section 8
    # tracks. The risk is confined to Malawi, whose duplicate district names are
    # water polygons and TA units; the count check below reports the join width.
    m <- dplyr::inner_join(z, hc, by = admin2_join_by(z, hc))
    m <- m[is.finite(m$y), , drop = FALSE]
    if (nrow(m) > nrow(z)) {
      cat(sprintf("    [fan] %s %s: join widened %d -> %d rows; cell skipped\n",
                  cn, nm, nrow(z), nrow(m)))
      next
    }
    if (nrow(m) < MIN_D || dplyr::n_distinct(m$Admin1) < 3) next
    X <- as.matrix(m[, COVS, drop = FALSE])
    oof <- rep(NA_real_, nrow(m))
    for (rg in unique(m$Admin1)) {
      i <- which(m$Admin1 == rg)
      fit <- .ds_fit(X[-i, , drop = FALSE], m$y[-i], k_screen = K_SCREEN)
      pr <- .ds_predict(fit, X[i, , drop = FALSE])
      if (!is.null(pr)) oof[i] <- pr
    }
    fin <- is.finite(oof)
    if (sum(fin) < MIN_D || stats::sd(m$y[fin]) == 0) next
    rows[[length(rows) + 1L]] <- data.frame(
      country = cn, target = nm, column = col[1], n_districts = sum(fin),
      r_oof = round(suppressWarnings(stats::cor(m$y[fin], oof[fin])), 4),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-13s %-20s n=%3d  r=%+.3f\n", cn, nm, sum(fin),
                utils::tail(rows, 1)[[1]]$r_oof))
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
readr::write_csv(res, here("results", "tables", "positive_control_targets.csv"))

cat("\n=== achieved out-of-fold r by target, leave-one-region-out ===\n")
print(as.data.frame(res |> group_by(target) |>
  summarise(countries = dplyr::n(), min_r = round(min(r_oof), 3),
            med_r = round(stats::median(r_oof), 3),
            max_r = round(max(r_oof), 3), .groups = "drop")), row.names = FALSE)

edu <- res[grepl("^education", res$target), ]
if (nrow(edu)) {
  cat(sprintf("\n=== CLAIM 3.12: earth observation predicts district education ===\n"))
  print(as.data.frame(edu[, c("country","target","n_districts","r_oof")]), row.names = FALSE)
  cat(sprintf("\nrange over %d education cells: %.3f to %.3f (claim states 0.48 to 0.71)\n",
              nrow(edu), min(edu$r_oof), max(edu$r_oof)))
}
