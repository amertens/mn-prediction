# =============================================================================
# scripts/_pipeline_setup.R
#
# Shared launch configuration for every pipeline runner (run_full_mode.R,
# run_tier1.R). Sourced, never run directly.
#
# WHY THIS IS A SEPARATE FILE
# ---------------------------
# It holds the environment defaults AND the four invalidation traps for
# settings {targets} cannot see. Those traps are the only thing standing
# between a run and a silently mixed configuration. When a second launcher was
# added (2026-08-30) the choice was to copy them or to share them; copying
# would have meant that adding a fifth trap to one launcher and not the other
# produces exactly the half-configured run the traps exist to prevent. So they
# live here once.
#
# Defines: say(), invalidate(), all_names, apply_invalidation_traps().
# Sourcing sets the environment; the CALLER must then call
# apply_invalidation_traps(), passing the target set it will actually build
# when it builds only a subset -- see that function's note.
# =============================================================================

# Each launcher sets RUNNER_TAG before sourcing so its log lines are
# attributable; defaults to "full" for callers that do not.
if (!exists("RUNNER_TAG")) RUNNER_TAG <- "full"
say <- function(...) cat(sprintf("[%s] %s  %s\n", RUNNER_TAG,
                                 format(Sys.time(), "%H:%M:%S"), sprintf(...)))

Sys.setenv(PIPELINE_MODE = "full", GEE_COVARIATE_HYGIENE = "true",
           # Read by tar_option_set(error=) in _targets.R. A batch run must
           # not throw away every remaining target because one worker died.
           TARGETS_ERROR_MODE = "continue")

# THIS RUN'S CONFIGURATION LIVES HERE, NOT IN THE CALLER'S SHELL.
#
# launch_full_detached.ps1 starts the pipeline through Win32_Process::Create,
# which parents it to the WMI provider host -- deliberately, so the run survives
# the terminal that launched it. The consequence is that it does NOT inherit the
# launching shell's environment. Anything set with `FOO=bar Rscript ...` or
# `$env:FOO` at the prompt is invisible to the run. Settings that matter must be
# defaults here, in a file under version control, where the log banner below
# also records them.
#
# These reach the future workers because _targets.R calls future::plan() AFTER
# this script's Sys.setenv(), and multisession workers inherit the master's
# environment at spawn time (the same mechanism MN_TARGETS_WORKER relies on).
defaults <- c(
  # Harmonised covariate vocabulary: 208 predictors shared by all five
  # countries, unit defects fixed. Set COVARIATE_VOCAB=legacy to reproduce the
  # published analysis of record from the raster-derived gee_* columns.
  COVARIATE_VOCAB        = "harmonized",
  # Anchor each region to its own design-based total. Beat national
  # benchmarking in 22 of 24 cells; see results/tables/admin1_arms.csv.
  AREA_BENCHMARK_ADMIN1  = "true",
  # Cote d'Ivoire external validation: 2.51 CPU-h and the 41-minute critical
  # path. Off for iteration runs, on for the reported one.
  RUN_CIV_VALIDATION     = "false",
  # Permutation-importance replicates; a diagnostic, not a reported estimate.
  VARIMP_N_PERM          = "2")
for (k in names(defaults))
  if (!nzchar(Sys.getenv(k))) do.call(Sys.setenv, setNames(list(defaults[[k]]), k))

say("PIPELINE_MODE=%s  GEE_COVARIATE_HYGIENE=%s", Sys.getenv("PIPELINE_MODE"),
    Sys.getenv("GEE_COVARIATE_HYGIENE"))
say("COVARIATE_VOCAB=%s  AREA_BENCHMARK_ADMIN1=%s  RUN_CIV_VALIDATION=%s  VARIMP_N_PERM=%s",
    Sys.getenv("COVARIATE_VOCAB"), Sys.getenv("AREA_BENCHMARK_ADMIN1"),
    Sys.getenv("RUN_CIV_VALIDATION"), Sys.getenv("VARIMP_N_PERM"))

if (identical(tolower(Sys.getenv("COVARIATE_VOCAB")), "harmonized")) {
  need <- list.files(file.path("data", "covariates", "country"),
                     pattern = "_predictors_admin2_raw\\.csv$")
  if (!length(need))
    stop("COVARIATE_VOCAB=harmonized but data/covariates/country/ is empty. ",
         "Run scripts/covariates/run_all.R first -- otherwise every area-level ",
         "target silently falls back to the legacy vocabulary mid-run.")
  say("harmonised covariates present for %d country/countries", length(need))
}

n_tif <- length(list.files("data/Tanzania_GEE_rasters", pattern = "\\.tif$"))
say("Tanzania rasters on disk: %d", n_tif)
if (n_tif == 0)
  warning("data/Tanzania_GEE_rasters/ is empty - Tanzania falls back to the ",
          "legacy-parity CSV and the covariate intersection stays small")

# Keep the pre-existing comparison table so the old and new numbers can be
# diffed afterwards rather than silently replaced.
src <- "results/tables/area_comparison_all.csv"
if (file.exists(src)) {
  # Tag the backup per runner so a tier-1 run does not overwrite the snapshot
  # taken before the last full run.
  bak <- sprintf("results/tables/area_comparison_all_PRE%s.csv", toupper(RUNNER_TAG))
  file.copy(src, bak, overwrite = TRUE)
  say("kept the pre-run comparison table at %s", bak)
}

# any_of(), NOT all_of(). tar_invalidate() matches against the METADATA store,
# not the manifest: a target that has never been built is absent there, so
# all_of() aborts the whole call with "Element `x` doesn't exist" and NOTHING
# gets invalidated. That failed silently on the first attempt -- the hygiene
# flag would have applied to some targets and not others, producing a run that
# was half hygienic and looked clean.
invalidate <- function(nms, label) {
  nms <- unique(nms[!is.na(nms) & nzchar(nms)])
  if (!length(nms)) { say("nothing to invalidate for %s", label); return(invisible()) }
  before <- tryCatch(nrow(targets::tar_meta(targets_only = TRUE)),
                     error = function(e) NA_integer_)
  ok <- tryCatch({ targets::tar_invalidate(names = tidyselect::any_of(nms)); TRUE },
                 error = function(e) { say("invalidate(%s) failed: %s", label,
                                           conditionMessage(e)); FALSE })
  after <- tryCatch(nrow(targets::tar_meta(targets_only = TRUE)),
                    error = function(e) NA_integer_)
  say("%s: requested %d, metadata rows %s -> %s (%s)", label, length(nms),
      before, after, if (ok) "ok" else "FAILED")
  if (!ok) stop(sprintf("invalidation for %s failed - refusing to run a ",
                        label), "partially-invalidated pipeline")
}

all_names <- tryCatch(targets::tar_manifest()$name, error = function(e) character(0))

# ── Applying the traps ───────────────────────────────────────────────────────
# SCOPE MATTERS, AND GETTING IT WRONG DELETES RESULTS.
#
# tar_invalidate() removes a target's METADATA; the target is only restored
# when the run rebuilds it. A full run rebuilds everything it invalidates, so
# scope never mattered. A SUBSET run does not: on 2026-08-30 the tier-1 run
# invalidated area_covariates_* (5) and area_model_* (24), neither of which is
# in the area_comparison_all closure, and finished without rebuilding them --
# leaving "target not found" for the dashboard's area layer and for any script
# reading those objects. The run reported success because, from {targets}'
# point of view, nothing it was asked to build had failed.
#
# So a runner that builds a subset must pass that subset as `scope`: traps are
# then applied only to targets the run will actually rebuild. Anything outside
# the scope is reported rather than silently invalidated, because a config
# change that cannot be honoured for those targets is exactly the mixed-state
# hazard these traps exist to prevent.
apply_invalidation_traps <- function(scope = NULL) {
  in_scope <- function(nms) {
    if (is.null(scope)) return(nms)
    keep <- intersect(nms, scope)
    skipped <- setdiff(nms, scope)
    if (length(skipped))
      say("NOT invalidating %d out-of-scope target(s) this run would not rebuild (e.g. %s)",
          length(skipped), paste(utils::head(skipped, 3), collapse = ", "))
    keep
  }
  invalidate(in_scope("pipeline_params"), "PIPELINE_MODE")
  invalidate(in_scope(grep("^(area_loco|area_comparison|transport|area_transport)",
                           all_names, value = TRUE)), "GEE_COVARIATE_HYGIENE")

  # TRAP 3. COVARIATE_VOCAB is read at CALL time by harmonized_vocab_active(),
  # inside extract_gee_admin2() and extract_area_covariates(). {targets} cannot
  # see it. Invalidating the two EXTRACTION targets is enough on its own --
  # everything downstream depends on them and is rebuilt by the dependency graph
  # -- but the covariate set is the input to the entire area-level analysis, so
  # getting this wrong would produce a run mixing two vocabularies.
  invalidate(in_scope(grep("^(gee_admin2_|area_covariates_)", all_names, value = TRUE)),
             "COVARIATE_VOCAB")

  # TRAP 4. AREA_BENCHMARK_ADMIN1 is read at call time inside
  # fit_area_level_model(). Only the area_model_* targets consume it.
  invalidate(in_scope(grep("^area_model_", all_names, value = TRUE)),
             "AREA_BENCHMARK_ADMIN1")
}
