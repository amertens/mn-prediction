# =============================================================================
# R_corrected/comparison.R — head-to-head corrected vs production assembly
#
# Collects the per-slice corrected targets, pulls the matching PRODUCTION
# metrics for reference, writes CSVs to results/tables/corrected/, and emits a
# single bundle (methods_comparison.rds) consumed by the dashboard
# "Methods comparison (corrected vs production)" tab.
# =============================================================================

# Production binary CV metric for a slice (the optimistic baseline).
prod_cv_for <- function(country_label, outcome_tag) {
  cv <- read_prod("cv_perf")
  if (is.null(cv)) return(NULL)
  sel <- cv[(cv$country == country_label) & (cv$outcome == outcome_tag), ,
            drop = FALSE]
  if ("model_type" %in% colnames(sel))
    sel <- sel[sel$model_type %in% c("binary", "binary_prob"), , drop = FALSE]
  if (nrow(sel) == 0) return(NULL)
  data.frame(country = country_label, outcome = outcome_tag,
             scheme = "PRODUCTION cv_perf (leaky preprocessing, random CV)",
             honest = FALSE,
             auc = round(sel$auc[1], 4),
             brier = if ("brier" %in% colnames(sel)) round(sel$brier[1], 4) else NA,
             auc_ci_lo = NA_real_, auc_ci_hi = NA_real_)  # (Issue 6) no honest CI for the leaky path
}

build_methods_comparison <- function(slices_results) {
  out_dir <- file.path("results", "tables", "corrected")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  bind <- function(field) {
    rows <- lapply(slices_results, function(s) s[[field]])
    rows <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, rows)
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
  }

  cv      <- bind("cv_perf")
  prodcv  <- bind("prod_cv")
  calib   <- bind("calibration")
  err     <- bind("admin2_error")
  dec     <- bind("decision")
  trust   <- bind("trust")
  area    <- bind("area_pp")
  intsum  <- bind("interval_summary")

  # CV honesty table: corrected schemes + production reference, side by side.
  cv_compare <- NULL
  if (!is.null(cv)) {
    cols <- intersect(c("country", "outcome", "scheme", "honest", "auc", "brier",
                        "auc_ci_lo", "auc_ci_hi"), colnames(cv))   # (Issue 6) carry AUC CI
    base <- cv[, cols]
    if (!is.null(prodcv)) base <- rbind(base, prodcv[, intersect(cols, colnames(prodcv))])
    cv_compare <- base[order(base$country, base$outcome, -base$auc), ]
  }

  write_if <- function(x, nm) if (!is.null(x))
    utils::write.csv(x, file.path(out_dir, nm), row.names = FALSE)
  write_if(cv_compare, "cv_honesty_compare.csv")
  write_if(calib,  "calibration_oof_compare.csv")
  write_if(err,    "admin2_error_sampling_adjusted.csv")
  write_if(dec,    "decision_value.csv")
  write_if(trust,  "trust_flags.csv")
  write_if(area,   "area_partial_pooling.csv")
  write_if(intsum, "interval_summary.csv")

  # ── Observability: per-field slice coverage (attempted vs surviving) ───────
  # A degraded run (silent fit failures) previously still wrote complete-looking
  # CSVs; this table makes the drop count explicit. See .log_skip() for the
  # per-fit reasons in the run log.
  n_slices <- length(slices_results)
  cov_fields <- c("cv_perf", "prod_cv", "calibration", "admin2_error",
                  "decision", "trust", "area_pp", "interval_summary")
  coverage <- data.frame(
    field = cov_fields, attempted = n_slices,
    surviving = vapply(cov_fields, function(f)
      sum(vapply(slices_results, function(s) {
        x <- s[[f]]; isTRUE(!is.null(x) && is.data.frame(x) && nrow(x) > 0)
      }, logical(1))), integer(1)),
    row.names = NULL)
  coverage$dropped <- coverage$attempted - coverage$surviving
  utils::write.csv(coverage, file.path(out_dir, "slice_coverage.csv"), row.names = FALSE)
  message(sprintf("[corrected] slice coverage  %s",
                  paste(sprintf("%s=%d/%d", coverage$field,
                                coverage$surviving, coverage$attempted),
                        collapse = "  ")))

  bundle <- list(
    cv_compare = cv_compare, calibration = calib, admin2_error = err,
    decision = dec, trust = trust, area_pp = area, interval_summary = intsum,
    build_time = NA, slices = names(slices_results))
  saveRDS(bundle, file.path(out_dir, "methods_comparison.rds"))
  # also drop a copy where the dashboard prep can find it
  dash <- file.path("dashboard", "data")
  if (dir.exists(dash)) saveRDS(bundle, file.path(dash, "methods_comparison.rds"))
  bundle
}
