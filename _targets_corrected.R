# =============================================================================
# _targets_corrected.R  — PARALLEL "corrected-methods" pipeline (P1–P8)
#
# Runs ALONGSIDE the production pipeline; writes to a SEPARATE store
# (_targets_corrected/) and never touches _targets.R or the production store.
# It reuses already-built PRODUCTION data outputs (outcome datasets, survey
# aggregates, GEE covariates) read from _targets_full/objects via read_prod(),
# so it does NOT recompute the heavy GEE extraction.
#
# Run (validate / smoke / full): see docs/CORRECTED_PIPELINE_RUN.md
#   targets::tar_manifest(script = "_targets_corrected.R", store = "_targets_corrected")
#   targets::tar_make(script = "_targets_corrected.R", store = "_targets_corrected")
#
# Scope is controlled by env var CORRECTED_SCOPE: "smoke" (default; one slice,
# tiny library) or "full" (all available country x outcome slices, fuller
# library). Compute guardrail: keep workers low; this is sequential by default.
# =============================================================================

library(targets)
suppressPackageStartupMessages({ library(here) })

tar_option_set(
  packages = c("dplyr", "glmnet", "rpart", "ranger", "here"),
  memory = "transient",
  garbage_collection = TRUE,
  storage = "main", retrieval = "main",
  # Resilience: a single bad slice returns NULL and the run continues, so the
  # final methods_comparison still builds from the slices that succeeded.
  # tar_errored() then lists any failures for the completion check.
  error = "null"
)

# corrected helpers
tar_source("R_corrected")
# reused production helpers (config + Fay-Herriot + impute/bound utilities)
source("R/config.R")
source("R/benchmark_models.R")

# ── Scope & library ─────────────────────────────────────────────────────────
scope <- Sys.getenv("CORRECTED_SCOPE", "smoke")
if (scope == "full") {
  LIB  <- c("mean", "glmnet", "ranger", "rpart")
  MAXP <- 800L            # bounded honest pool (memory guardrail)
  V    <- 5L
} else {
  LIB  <- c("mean", "glmnet", "ranger")
  MAXP <- 400L
  V    <- 5L
}

cfgs <- get_country_configs()
to_lower <- function(x) tolower(x)

# Available production outcome-data objects define which slices we can build.
prod_objs <- list.files(file.path(Sys.getenv("PROD_STORE", "_targets_full"),
                                  "objects"), pattern = "^outcome_data_")
avail_suffix <- sub("^outcome_data_", "", prod_objs)

# Build the slice table (country_label, country_lower, outcome_tag, suffix).
slices <- list()
for (cl in names(cfgs)) {
  low <- to_lower(cl)
  for (oc in names(cfgs[[cl]]$outcomes)) {
    suf <- paste0(low, "_", oc)
    if (!suf %in% avail_suffix) next
    slices[[suf]] <- list(label = cfgs[[cl]]$country, low = low,
                          outcome = oc, suffix = suf)
  }
}
if (scope != "full") {
  pick <- if ("gambia_women_iron" %in% names(slices)) "gambia_women_iron"
          else names(slices)[1]
  slices <- slices[pick]
}
message(sprintf("[corrected] scope=%s | %d slice(s): %s",
                scope, length(slices), paste(names(slices), collapse = ", ")))

# ── Per-country GEE covariate targets (reused across that country's slices) ──
countries_used <- unique(vapply(slices, function(s) s$low, character(1)))
gee_targets <- lapply(countries_used, function(low) {
  tar_target_raw(paste0("gee_", low),
                 substitute(read_prod(nm), list(nm = paste0("gee_admin2_", low))))
})

# ── Per-slice target factory ────────────────────────────────────────────────
slice_targets <- list()
for (s in slices) {
  suf <- s$suffix; low <- s$low; lab <- s$label; ocn <- s$outcome
  cc_sym <- substitute(get_country_configs()[[L]], list(L = lab))
  oc_sym <- substitute(get_country_configs()[[L]]$outcomes[[O]],
                       list(L = lab, O = ocn))

  slice_targets <- c(slice_targets, list(
    # ingest production outputs
    tar_target_raw(paste0("od_", suf),
      substitute(read_prod(nm), list(nm = paste0("outcome_data_", suf)))),
    tar_target_raw(paste0("svy_", suf),
      substitute(read_prod(nm), list(nm = paste0("svy_admin2_", suf)))),

    # P1: leakage-free SL fit (honest cluster + spatial block + optimistic)
    tar_target_raw(paste0("sl_", suf),
      substitute(fit_corrected_sl(OD, CC, OC, LAB, library = LIBR,
                                  max_pred = MP, V = VV),
        list(OD = as.symbol(paste0("od_", suf)), CC = cc_sym, OC = oc_sym,
             LAB = lab, LIBR = LIB, MP = MAXP, VV = V))),
    tar_target_raw(paste0("cvperf_", suf),
      substitute(extract_cv_perf_corrected(SL),
                 list(SL = as.symbol(paste0("sl_", suf))))),
    tar_target_raw(paste0("prodcv_", suf),
      substitute(prod_cv_for(LAB, OCN), list(LAB = lab, OCN = ocn))),

    # P2 calibration OOF
    tar_target_raw(paste0("calib_", suf),
      substitute(diagnostics_calibrated_oof(SL),
                 list(SL = as.symbol(paste0("sl_", suf))))),
    # P5 sampling-error-aware error
    tar_target_raw(paste0("err_", suf),
      substitute(admin2_error_corrected(SL, SV),
        list(SL = as.symbol(paste0("sl_", suf)),
             SV = as.symbol(paste0("svy_", suf))))),
    # P3 decision value
    tar_target_raw(paste0("dec_", suf),
      substitute(decision_value_corrected(SL, SV),
        list(SL = as.symbol(paste0("sl_", suf)),
             SV = as.symbol(paste0("svy_", suf))))),
    # P4 intervals (split-conformal + design-based)
    tar_target_raw(paste0("int_", suf),
      substitute(intervals_corrected(SL, SV),
        list(SL = as.symbol(paste0("sl_", suf)),
             SV = as.symbol(paste0("svy_", suf))))),
    # P6 trust flags
    tar_target_raw(paste0("trust_", suf),
      substitute(trust_flags_corrected(SL, SV, GEE, INT),
        list(SL = as.symbol(paste0("sl_", suf)),
             SV = as.symbol(paste0("svy_", suf)),
             GEE = as.symbol(paste0("gee_", low)),
             INT = as.symbol(paste0("int_", suf))))),
    # P7 partial-pooling area estimator
    tar_target_raw(paste0("area_", suf),
      substitute(area_partial_pooling_corrected(SL, SV, GEE),
        list(SL = as.symbol(paste0("sl_", suf)),
             SV = as.symbol(paste0("svy_", suf)),
             GEE = as.symbol(paste0("gee_", low)))),),

    # per-slice result bundle (named list element)
    tar_target_raw(paste0("result_", suf),
      substitute(list(cv_perf = CV, prod_cv = PCV, calibration = CAL,
                      admin2_error = ERR, decision = DEC, trust = TR,
                      area_pp = AR, interval_summary = INT$summary),
        list(CV = as.symbol(paste0("cvperf_", suf)),
             PCV = as.symbol(paste0("prodcv_", suf)),
             CAL = as.symbol(paste0("calib_", suf)),
             ERR = as.symbol(paste0("err_", suf)),
             DEC = as.symbol(paste0("dec_", suf)),
             TR = as.symbol(paste0("trust_", suf)),
             AR = as.symbol(paste0("area_", suf)),
             INT = as.symbol(paste0("int_", suf)))))
  ))
}

# ── Final comparison bundle (corrected vs production) ───────────────────────
result_syms <- setNames(lapply(names(slices), function(suf)
  as.symbol(paste0("result_", suf))), names(slices))
result_list_expr <- as.call(c(list(as.symbol("list")), result_syms))
comparison_target <- tar_target_raw(
  "methods_comparison",
  substitute(build_methods_comparison(RL), list(RL = result_list_expr)))

c(gee_targets, slice_targets, list(comparison_target))
