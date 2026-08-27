# =============================================================================
# scripts/run_full_after_fast.R
#
# Wait for the in-flight FAST rebuild to finish, verify it succeeded, then run
# the whole pipeline again in FULL mode with GEE covariate hygiene on.
#
# The fast pass is a correctness check on the 2026-08 fixes (water-body
# snapping, Tanzania's restored rasters, duplicate-key collapsing, the
# Fay-Herriot variance guard, national benchmarking). The full pass produces the
# publication numbers: 13 SuperLearner learners instead of 5 (and 500 area-level
# bootstrap replicates instead of 10). Those are the only two live differences --
# see the table in README.md, which was corrected in 2026-08 after this comment
# was found to disagree with both the code and the pipeline's own banner.
#
# TWO INVALIDATION TRAPS, both handled below.
#
#  1. PIPELINE_MODE is read by get_pipeline_params() when the `pipeline_params`
#     target is BUILT. targets caches that target, so changing the environment
#     variable alone changes nothing -- pipeline_params must be invalidated by
#     hand. (The README says this; it is easy to miss.)
#
#  2. GEE_COVARIATE_HYGIENE is worse: gee_hygiene_enabled() reads it at
#     CALL time inside prune_gee_covariates(), so it never enters the dependency
#     graph at all. targets cannot see that it changed, and any target that
#     consumes it but is not otherwise invalidated would silently keep its
#     un-pruned result. The consumers are build_area_loco_dataset() and
#     assemble_area_transport(), so their targets are invalidated by name.
#
# Launch detached (PowerShell), so it survives the session that starts it:
#
#   Start-Process -FilePath "C:\Program Files\R\R-4.4.2\bin\x64\Rscript.exe" `
#     -ArgumentList "scripts/run_full_after_fast.R" `
#     -RedirectStandardOutput "pipeline_full.log" `
#     -RedirectStandardError  "pipeline_full.err" -WindowStyle Hidden
# =============================================================================

FAST_LOG   <- "pipeline_rerun_2026-08.log"
POLL_SEC   <- 120
MAX_WAIT_H <- 24

say <- function(...) cat(sprintf("[full] %s  %s\n", format(Sys.time(), "%H:%M:%S"),
                                 sprintf(...)))

# ---------------------------------------------------------------------------
# 1. Wait for the fast run
# ---------------------------------------------------------------------------
fast_done <- function() {
  if (!file.exists(FAST_LOG)) return(NA)          # never started
  txt <- tryCatch(readLines(FAST_LOG, warn = FALSE), error = function(e) character(0))
  hit <- grep("^\\[rerun\\] finished", txt, value = TRUE)
  if (!length(hit)) return(NA)                    # still running
  grepl("success = TRUE", utils::tail(hit, 1), fixed = TRUE)
}

say("waiting for the fast rebuild to finish (polling every %d s)", POLL_SEC)
t0 <- Sys.time(); status <- NA
repeat {
  status <- fast_done()
  if (!is.na(status)) break
  if (as.numeric(difftime(Sys.time(), t0, units = "hours")) > MAX_WAIT_H) {
    say("gave up waiting after %d h - fast run never reported completion", MAX_WAIT_H)
    quit(status = 1)
  }
  Sys.sleep(POLL_SEC)
}

if (!isTRUE(status)) {
  say("fast rebuild reported FAILURE - not starting the full run")
  say("inspect %s before retrying", FAST_LOG)
  quit(status = 1)
}
say("fast rebuild completed successfully after %.1f h of waiting",
    as.numeric(difftime(Sys.time(), t0, units = "hours")))

# Record what fast produced, so the two runs can be compared afterwards.
fast_snapshot <- "results/tables/area_comparison_all_FASTMODE.csv"
src <- "results/tables/area_comparison_all.csv"
if (file.exists(src)) {
  dir.create(dirname(fast_snapshot), showWarnings = FALSE, recursive = TRUE)
  file.copy(src, fast_snapshot, overwrite = TRUE)
  say("kept a copy of the fast-mode comparison table at %s", fast_snapshot)
}

# ---------------------------------------------------------------------------
# 2. Hand off to the shared full-mode runner
#
# The mode switch, the two invalidation traps and the post-run verification all
# live in scripts/run_full_mode.R so this script and a direct full-mode launch
# cannot drift apart.
# ---------------------------------------------------------------------------
say("handing off to scripts/run_full_mode.R")
source("scripts/run_full_mode.R")
