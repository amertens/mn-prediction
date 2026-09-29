# =============================================================================
# explore/scripts/30_within_country_diet.R   [DC-01]
#
# DOES THE DIETARY BLOCK PREDICT WITHIN COUNTRIES?
#
# The LOCO ablation said no - Household diet scored -0.134 alone, MIMI and RTFP
# did not appear at all. DIETARY_COVERAGE_2026-09-29 showed why that test was
# unfair: LOCO pools on columns present in ALL FOUR countries
# (23_domain_ablation_loco.R:81), which deletes 100% of MIMI, 100% of RTFP,
# 87% of Household diet and 67% of Food prices, while deleting 0% of climate,
# soil and satellite. The dietary verdict was computed on 2 of 15 columns, or
# on none.
#
# The in-fill and region estimands need no intersect: each country is fitted on
# its own columns. So the dietary block can be tested fairly today, without
# waiting for Ghana's HCES gap to be filled.
#
# DESIGN. Three predictor sets per country x outcome x target, identical rows
# and identical folds, scored with the estimator of record (zero-tuning domain
# index):
#
#   full     every predictor that survives the leakage filter
#   nodiet   full minus the six dietary domains
#   diet     the six dietary domains alone
#
#   full - nodiet  is what diet ADDS on top of everything else
#   diet           is what diet carries on its own
#
# Read PAIRED, per cell, and blocked by country x target - the two rules this
# folder has learned the hard way. A median over cells hides sign flips.
#
# "Dietary" here is intake and food economics, NOT food production:
# Agricultural production is remotely sensed, is already known load-bearing
# (delta +0.012, 0.242 alone under LOCO) and is left in `nodiet` on purpose, so
# this measures the survey-derived dietary block specifically.
#
# Tiers are forced to open,survey_public,survey_dhs to match the production
# ablation this is compared against (the harness default omits survey_dhs,
# which would silently drop the IYCF and adult-nutrition columns).
#
#   Rscript explore/scripts/30_within_country_diet.R
#   EXP_DC_REPS=10
# -> explore/out/30_within_country_diet.csv        per cell x set x estimand
#    explore/out/30_within_country_diet_solo.csv   each dietary domain alone
# =============================================================================
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public,survey_dhs")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_DC_REPS", "10"))
set.seed(20260929L)
E <- exp_load()

DIET_DOMAINS <- c("Household diet and consumption (HCES)",
                  "Dietary inadequacy (MODELLED SURFACE, HCES)",
                  "Market prices (RTFP)",
                  "Food prices and supply",
                  "Infant and young child feeding",
                  "Adult nutrition")

diet_cols <- names(E$domain_of)[E$domain_of %in% DIET_DOMAINS]
diet_cols <- intersect(diet_cols, E$PREDS)
rest_cols <- setdiff(E$PREDS, diet_cols)
message("predictors ", length(E$PREDS), " | dietary ", length(diet_cols),
        " | rest ", length(rest_cols))
for (d in DIET_DOMAINS)
  message(sprintf("   %-46s %3d col(s)", substr(d, 1, 46),
                  sum(E$domain_of[diet_cols] == d)))

ARM <- exp_baseline_arms("domain_index")
IX  <- exp_cell_index(E)

raw <- list(); solo <- list()
for (target in c("level", "prev")) {
  for (i in seq_len(nrow(IX))) {
    cn <- IX$country[i]; on <- IX$outcome[i]
    sets <- list(full   = list(cols = E$PREDS,   min = 20L),
                 nodiet = list(cols = rest_cols, min = 20L),
                 diet   = list(cols = diet_cols, min = 5L))
    cells <- list()
    for (s in names(sets))
      cells[[s]] <- tryCatch(exp_cell(E, cn, on, target, cols = sets[[s]]$cols,
                                      min_cols = sets[[s]]$min),
                             error = function(e) NULL)
    ok <- !vapply(cells, is.null, TRUE)
    if (!ok[["full"]] || !ok[["nodiet"]]) next
    # identical rows are required for a paired read
    if (cells$full$n != cells$nodiet$n) { message("  n mismatch, skipping ", cn, " ", on); next }

    for (s in names(cells)) {
      cc <- cells[[s]]; if (is.null(cc)) next
      a <- exp_infill(cc, ARM, reps = REPS)
      b <- exp_region(cc, ARM)
      for (z in list(a, b)) if (!is.null(z) && nrow(z)) {
        z$set <- s; z$n_cols <- ncol(cc$X); raw[[length(raw) + 1L]] <- z }
    }
    message(sprintf("  %-12s %-14s %-5s  n=%3d  cols full/nodiet/diet = %d/%d/%s",
                    cn, on, target, cells$full$n, ncol(cells$full$X), ncol(cells$nodiet$X),
                    if (ok[["diet"]]) ncol(cells$diet$X) else "-"))

    # each dietary domain alone, in-fill only
    for (d in DIET_DOMAINS) {
      dc <- intersect(names(E$domain_of)[E$domain_of == d], E$PREDS)
      if (length(dc) < 3) next
      cc <- tryCatch(exp_cell(E, cn, on, target, cols = dc, min_cols = 3L),
                     error = function(e) NULL)
      if (is.null(cc)) next
      z <- exp_infill(cc, ARM, reps = REPS)
      if (!is.null(z) && nrow(z)) { z$domain <- d; z$n_cols <- ncol(cc$X)
                                    solo[[length(solo) + 1L]] <- z }
    }
  }
}

R <- dplyr::bind_rows(raw)
SM <- R |> dplyr::group_by(country, outcome, target, estimand, set) |>
  dplyr::summarise(n_areas = dplyr::first(n_areas), n_cols = dplyr::first(n_cols),
                   spearman = mean(spearman, na.rm = TRUE),
                   wmae = mean(wmae, na.rm = TRUE), .groups = "drop") |> as.data.frame()
exp_write(SM, "30_within_country_diet")

# ── paired read, blocked by country x target ────────────────────────────────
W <- SM |> dplyr::select(country, outcome, target, estimand, set, spearman) |>
  tidyr::pivot_wider(names_from = set, values_from = spearman) |> as.data.frame()
W$adds <- W$full - W$nodiet

cat("\n===== DC-01: does the dietary block add anything WITHIN country? =====\n")
for (es in c("infill", "region")) {
  w <- W[W$estimand == es & is.finite(W$adds), ]
  if (!nrow(w)) next
  cat(sprintf("\n--- %s (%d cells) ---\n", es, nrow(w)))
  cat(sprintf("  full - nodiet : better in %2d of %2d | median %+.4f | mean %+.4f\n",
              sum(w$adds > 0), nrow(w), median(w$adds), mean(w$adds)))
  cat(sprintf("  diet alone    : median %+.3f (full %.3f, nodiet %.3f)\n",
              median(w$diet, na.rm = TRUE), median(w$full), median(w$nodiet)))
  cat("\n  by country x target (the blocking that matters):\n")
  for (cn in sort(unique(w$country))) for (tg in sort(unique(w$target))) {
    z <- w[w$country == cn & w$target == tg, ]
    if (!nrow(z)) next
    cat(sprintf("    %-12s %-5s  adds: %d of %d | median %+.4f | diet alone %+.3f\n",
                cn, tg, sum(z$adds > 0), nrow(z), median(z$adds),
                median(z$diet, na.rm = TRUE)))
  }
}

if (length(solo)) {
  S2 <- dplyr::bind_rows(solo) |>
    dplyr::group_by(country, outcome, target, domain) |>
    dplyr::summarise(n_cols = dplyr::first(n_cols),
                     spearman = mean(spearman, na.rm = TRUE), .groups = "drop") |> as.data.frame()
  exp_write(S2, "30_within_country_diet_solo")
  cat("\n===== each dietary domain ALONE, in-fill, median over cells =====\n")
  d <- S2 |> dplyr::group_by(domain, target) |>
    dplyr::summarise(cells = dplyr::n(), countries = dplyr::n_distinct(country),
                     med_cols = stats::median(n_cols),
                     median_rho = round(stats::median(spearman, na.rm = TRUE), 3),
                     positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
  print(d[order(d$target, -d$median_rho), ], row.names = FALSE)
}
cat("\nDONE\n")
