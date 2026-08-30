# =============================================================================
# scripts/run_full_mode.R
#
# Rebuild the whole pipeline in FULL mode with GEE covariate hygiene on.
#
# LAUNCH VIA scripts/launch_full_detached.ps1, not Start-Process:
#
#   powershell -ExecutionPolicy Bypass -File scripts/launch_full_detached.ps1
#
# Start-Process leaves the run inside the calling application's process tree,
# so it dies when that application does. On 2026-08-23 the host app logged a
# MoAppHang and took the workers down with it 8 minutes into a rebuild. The
# launcher uses Win32_Process::Create, which parents to the WMI host instead.
#
# WHAT FULL MODE ACTUALLY CHANGES (verified against the code, 2026-08):
#   sl_stack  5 learners -> 12 (adds ridge, ranger_low_mtry, ranger_deep,
#             xgb_deep, bart_small/100/200). BART is the one that holds up on
#             the rare outcomes -- women's vitamin A and B12 sit at 1-3%
#             prevalence -- where the others collapse to the mean.
#             gaussianprocess was the 13th learner until 2026-08-24. It is
#             O(n^3) and made the pooled LOCO targets (n = 10k-13k) effectively
#             non-terminating. The gp_sensitivity target measured what dropping
#             it costs on Ghana x child_iron: AUC 0.7016 with and without, to
#             16 decimal places -- it carried zero ensemble weight.
#   B_area    10 -> 500 area-level bootstrap replicates.
#   Everything else (K, prescreening, conformal intervals, ablation) is
#   IDENTICAL. Class-weighted learners are added for any binary outcome under
#   15% prevalence in BOTH modes, so a rare outcome gets 7 learners in fast and
#   14 in full.
#
# TWO INVALIDATION TRAPS, both handled below.
#
#  1. PIPELINE_MODE is read by get_pipeline_params() when the `pipeline_params`
#     target is BUILT. {targets} caches that target, so setting the environment
#     variable alone changes nothing.
#
#  2. GEE_COVARIATE_HYGIENE is worse: gee_hygiene_enabled() reads it at CALL
#     time inside prune_gee_covariates(), so it never enters the dependency
#     graph. {targets} cannot see that it changed, and any consuming target not
#     otherwise invalidated would silently keep its un-pruned result -- a run
#     that is half hygienic and looks clean. The consumers are
#     build_area_loco_dataset() and assemble_area_transport().
# =============================================================================

RUNNER_TAG <- "full"
# Environment defaults + the four invalidation traps for settings {targets}
# cannot see. Shared with scripts/run_tier1.R so the two launchers can never
# drift apart -- see the header of that file.
source(here::here("scripts", "_pipeline_setup.R"))
# A full run rebuilds everything it invalidates, so no scope restriction.
apply_invalidation_traps()

n_out <- tryCatch(length(targets::tar_outdated(reporter = "silent")),
                  error = function(e) NA_integer_)
say("outdated before the run: %s of %d", n_out, length(all_names))

# USE tar_make_future(), NOT tar_make().
#
# _targets.R sets a future plan and its own banner says "Run with:
# targets::tar_make_future(workers = N)". Plain tar_make() is SEQUENTIAL across
# targets -- it dispatches one at a time, so on this 22-core machine the first
# attempt sat at 12% CPU and projected out to several days. tar_make_future()
# dispatches independent targets concurrently.
#
# Workers are capped against BOTH cores and RAM, and RAM is the binding
# constraint here -- not cores.
#
# MEASURED 2026-08 on this 31.4 GB / 22-core machine, mid-run: individual SL
# workers reached 4.45 and 2.67 GB, so the old "~1.5 GB per worker" figure was
# roughly 3x optimistic. At 6 workers the run died at 22 minutes with
# std::bad_alloc and "could not start R" -- an OOM, not the process explosion
# that killed the previous attempt.
#
# Size from memory that is ACTUALLY FREE, not from total RAM. Sizing off the
# total silently assumes the machine is idle. It is not: with Chrome and the
# editor resident, 4 workers drove commit charge to 93% of the 69.9 GB limit
# with 4.4 GB available -- the same condition that produced std::bad_alloc
# earlier -- even though total RAM (31.4 GB) suggested 4 was fine.
#
# Budget 3.5 GB per worker against free memory, holding 2 GB back for the OS
# and the callr master (itself only ~0.25 GB now that storage/retrieval are
# set to "worker"). Cap at 6: past that the store I/O contends and the memory
# risk outweighs the throughput.
# The detached launcher cannot pass this in -- Win32_Process::Create does not
# inherit the calling shell's environment -- so the value for a run lives here.
# Sized against MEASURED worker commit (6.5-7.6 GB), not the 3.5 GB budget the
# previous autosizer assumed, which was roughly 2x optimistic. Set TARGETS_WORKERS
# to override for a foreground run.
n_workers <- as.integer(Sys.getenv("TARGETS_WORKERS", NA))
if (is.na(n_workers) || n_workers < 1L) n_workers <- 2L
Sys.setenv(TARGETS_WORKERS = n_workers)   # _targets.R reads this for its plan
say("dispatching with tar_make_future(workers = %d)", n_workers)

t1 <- Sys.time()
ok <- tryCatch({ targets::tar_make_future(workers = n_workers, reporter = "summary"); TRUE },
               error = function(e) { say("FAILED: %s", conditionMessage(e)); FALSE })
say("finished after %.1f h, success = %s",
    as.numeric(difftime(Sys.time(), t1, units = "hours")), ok)

# Report failures unconditionally. With TARGETS_ERROR_MODE="continue" the run
# can finish "successfully" while individual targets failed, so gating this on
# !ok would hide exactly the errors that mode exists to let through.
{
  e <- tryCatch(targets::tar_errored(), error = function(e) character(0))
  say("errored targets: %d", length(e))
  if (length(e)) {
    print(utils::head(e, 20))
    m <- tryCatch(targets::tar_meta(fields = "error"), error = function(e) NULL)
    if (!is.null(m)) {
      m <- m[m$name %in% e & !is.na(m$error), c("name", "error")]
      for (i in seq_len(min(nrow(m), 5)))
        cat(sprintf("  %s: %s\n", m$name[i], substr(m$error[i], 1, 300)))
    }
  }
}

# Report outdated EXCLUDING the always-run stamp targets and their downstream.
#
# dhs_stamp_* and gee_parity_stamp_* are declared cue = tar_cue("always") so
# they re-check their input files on every run. tar_outdated() propagates
# downstream, so those 10 targets alone make ~850 of 898 report "outdated"
# FOREVER, even immediately after a clean run that built everything. That made
# every run of this script end with an alarming "targets still outdated: 853"
# next to "errored targets: 0", which is a contradiction only in appearance.
#
# tar_sitrep() gives PER-TARGET staleness without the downstream propagation,
# so counting targets whose own flags are set is the honest number.
left <- tryCatch(targets::tar_outdated(reporter = "silent"),
                 error = function(e) character(0))
genuine <- tryCatch({
  s <- as.data.frame(targets::tar_sitrep())
  fl <- intersect(c("meta", "command", "depend", "format", "repository",
                    "iteration", "file", "seed"), names(s))
  s$name[apply(s[, fl, drop = FALSE], 1, function(r) any(r %in% TRUE))]
}, error = function(e) NA_character_)
if (length(genuine) == 1 && is.na(genuine)) {
  say("targets still outdated: %d (sitrep unavailable)", length(left))
} else {
  say("targets outdated: %d reported, %d genuinely stale", length(left), length(genuine))
  say("  (the difference is downstream of %d always-run stamp targets)",
      length(left) - length(genuine))
  if (length(genuine)) print(utils::head(sort(genuine), 20))
}

# Confirm the switches took effect rather than assuming they did.
pp <- tryCatch(targets::tar_read(pipeline_params), error = function(e) NULL)
if (!is.null(pp))
  say("pipeline_params reports mode=%s hygiene=%s", pp$mode, pp$gee_covariate_hygiene)

say("next: Rscript scripts/refresh_slide_numbers.R, then update docs/mn_prediction_slides.qmd")
