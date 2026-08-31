# =============================================================================
# scripts/accuracy_impact/ws7a_guards.R
#
# WS7a. Run the three structural guards and write their artefacts.
#
#   PROFILE=smoke  restricts the leakage report to Ghana child_iron.
#
#   Rscript scripts/accuracy_impact/ws7a_guards.R
# -> results/tables/leakage_report.csv          (ranked columns per cell)
# -> results/tables/leakage_report_summary.csv  (one line per cell)
# -> results/tables/admin2_unit_reference.csv   (data-derived unit counts)
# -> results/tables/admin2_unit_check.csv       (every consumed table checked)
# =============================================================================
suppressPackageStartupMessages({library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE   <- here("_targets_full")
SUF     <- if (PROFILE == "smoke") "_SMOKE" else ""
TDIR    <- here("results", "tables")
cat(sprintf("[ws7a] profile=%s\n", PROFILE))

# ── 1. Unit-count reference, derived from the survey targets ────────────────
ref <- admin2_analytic_counts(STORE)
if (is.null(ref)) stop("No svy_admin2_* targets readable from ", STORE)
readr::write_csv(ref, file.path(TDIR, sprintf("admin2_unit_reference%s.csv", SUF)))
cat("\n=== analytic Admin-2 units per country (from the store) ===\n")
print(ref[, c("country", "analytic_units", "n_outcomes")], row.names = FALSE)

# ── 2. Check every table the anchoring, benchmarking and dashboard code reads ─
# `count_col` differs across tables because they were written at different
# times; the check normalises on the country label, not the column name.
CONSUMED <- list(
  list("admin1_arms.csv",                    "n_admin2", "country"),
  list("individual_anchor.csv",              "n_units",  "country"),
  list("area_comparison_all.csv",            "n_admin2", "country"),
  list("national_composition.csv",           "n_admin2", "country"),
  list("resolution_comparison.csv",          "n_units",  "country"),
  list("corrected/protocol_reconciliation.csv", "n_area", "country")
)
chk <- list()
for (s in CONSUMED) {
  p <- file.path(TDIR, s[[1]])
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", s[[1]])); next }
  tb <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
  # individual_anchor mixes district and cluster rows; only the district rows
  # are Admin-2 units, and clusters legitimately outnumber districts.
  if (grepl("individual_anchor", s[[1]]) && "unit" %in% names(tb))
    tb <- tb[tb$unit == "district", , drop = FALSE]
  r <- check_unit_counts(tb, s[[2]], ref, country_col = s[[3]], label = s[[1]])
  if (!is.null(r)) chk[[length(chk) + 1L]] <- r
}
chk <- do.call(rbind, chk)
readr::write_csv(chk, file.path(TDIR, sprintf("admin2_unit_check%s.csv", SUF)))
cat("\n=== unit-count check (OVER is the pair-key fan signature) ===\n")
print(as.data.frame(chk), row.names = FALSE)
cat(sprintf("\nOVER rows: %d of %d\n", sum(chk$status == "OVER"), nrow(chk)))

# ── 3. Leakage report ───────────────────────────────────────────────────────
cs <- if (PROFILE == "smoke") "Ghana" else NULL
os <- if (PROFILE == "smoke") "child_iron" else NULL
rep <- build_leakage_report(STORE, countries = cs, outcomes = os)
if (is.null(rep)) stop("Leakage report produced no rows.")
readr::write_csv(rep, file.path(TDIR, sprintf("leakage_report%s.csv", SUF)))
sm <- leakage_report_summary(rep)
readr::write_csv(sm, file.path(TDIR, sprintf("leakage_report_summary%s.csv", SUF)))

cat("\n=== leakage report: per-cell maxima ===\n")
print(as.data.frame(sm[order(-sm$max_abs_r), ]), row.names = FALSE)
cat(sprintf("\ncells: %d | flagged columns (abs r >= %.2f): %d | overall max abs r: %.4f\n",
            nrow(sm), rep$threshold[1], sum(sm$n_flagged), max(sm$max_abs_r)))
if (any(sm$n_flagged > 0)) {
  cat("\n--- FLAGGED COLUMNS ---\n")
  print(as.data.frame(rep[rep$flagged, c("country","outcome","predictor_set","column","r")]),
        row.names = FALSE)
}
