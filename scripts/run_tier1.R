# =============================================================================
# scripts/run_tier1.R
#
# Rebuild only the MODEL-COMPARISON slice of the pipeline: everything
# `area_comparison_all` depends on, and nothing else.
#
# WHY
# ---
# The full graph is ~844 targets. Measured 2026-08-30 against the dependency
# graph, `area_comparison_all` -- the table carrying the five area arms, the
# precision-weighting contrast and the prescreen-optimism contrast -- has a
# closure of ~155 targets. The other ~690 are the corrected_* protocol block
# (265), varimp/shap (50), diagnostics (50), transportability/LOCO (41),
# ablation (26), and several families that feed no manuscript output at all
# (plots 26, conformal_ci 24, area_model 24, gp_sensitivity 3).
#
# So while you are iterating on covariate construction and comparing modelling
# approaches, this runs the part that answers the question and skips the part
# that does not. Use scripts/run_full_mode.R for the reported build.
#
# The closure is COMPUTED from the graph, not hardcoded, so it tracks
# _targets.R automatically -- a hardcoded list would rot the first time a
# dependency is added.
#
# LAUNCH the same way as the full run (see run_full_mode.R's header on why
# Start-Process is the wrong parent):
#   powershell -ExecutionPolicy Bypass -File scripts/launch_full_detached.ps1
# or, for a foreground run:
#   TARGETS_WORKERS=3 Rscript scripts/run_tier1.R
# =============================================================================

RUNNER_TAG <- "tier1"
# Environment defaults + the four invalidation traps for settings {targets}
# cannot see. SHARED with run_full_mode.R on purpose: duplicating the traps
# would mean a trap added to one launcher and not the other yields exactly the
# half-configured run the traps exist to prevent.
source(here::here("scripts", "_pipeline_setup.R"))

# ── Compute the closure of the comparison table ──────────────────────────────
SEEDS <- c("area_comparison_all")

net <- targets::tar_network(targets_only = TRUE)
edges <- net$edges
closure <- intersect(SEEDS, unique(c(edges$from, edges$to)))
if (!length(closure))
  stop("none of the seed targets exist in the graph: ",
       paste(SEEDS, collapse = ", "))
repeat {
  parents <- unique(edges$from[edges$to %in% closure])
  grown <- union(closure, parents)
  if (length(grown) == length(closure)) break
  closure <- grown
}

# Traps AFTER the closure exists, and scoped to it. Applying them to the whole
# graph would invalidate targets this run never rebuilds -- which is exactly
# what deleted area_covariates_* and area_model_* on the first tier-1 run.
apply_invalidation_traps(scope = closure)

say("tier-1 closure: %d of %d targets (%.0f%% of the graph)",
    length(closure), length(all_names),
    100 * length(closure) / max(1L, length(all_names)))

# Report what is deliberately NOT being built, by family, so the omission is a
# visible choice rather than something a reader has to infer from a count.
skipped <- setdiff(all_names, closure)
fam <- sub("_(gambia|ghana|sierraleone|malawi|tanzania|civ)(_.*)?$", "",
           tolower(skipped))
fam <- sub("_(child|women)_(iron|vita|zinc|b12|folate)$", "", fam)
top <- sort(table(fam), decreasing = TRUE)
say("skipping %d targets; largest families: %s", length(skipped),
    paste(sprintf("%s(%d)", names(top)[seq_len(min(8, length(top)))],
                  as.integer(top)[seq_len(min(8, length(top)))]),
          collapse = " "))

# NOTE THE !! -- IT IS LITERAL, NOT DECORATION.
# tar_outdated()/tar_make_future() capture `names` as a QUOSURE and evaluate it
# inside a separate callr process, where this script's locals do not exist.
# Passing any_of(closure) fails there with "object 'closure' not found", and
# because the call is wrapped in tryCatch it degrades to NA / a failed run
# rather than an obvious error. Unquoting splices the character vector itself
# into the quosure, so the callr process receives values, not a name to look up.
n_out <- tryCatch(
  length(targets::tar_outdated(names = tidyselect::any_of(!!closure),
                               reporter = "silent")),
  error = function(e) NA_integer_)
say("outdated within the closure: %s of %d", n_out, length(closure))

# ── Workers ──────────────────────────────────────────────────────────────────
# RAM binds here, not cores. Measured 2026-08-29, a single pipeline worker
# committed 7.6 GB, so run_full_mode.R's 3.5 GB/worker budget is roughly 2x
# optimistic. Tier 1 omits the individual-level SL fits that drove the worst of
# that, so it is lighter than a full run at the same worker count.
#
# The default lives HERE rather than in the launcher because
# launch_full_detached.ps1 starts the run through Win32_Process::Create, which
# does not inherit the calling shell's environment -- its -Workers flag sets
# $env:TARGETS_WORKERS in a process the run never sees, so it silently has no
# effect. Same reason the run's configuration lives in _pipeline_setup.R.
# Override for a foreground run with TARGETS_WORKERS=N Rscript scripts/run_tier1.R
n_workers <- suppressWarnings(as.integer(Sys.getenv("TARGETS_WORKERS", NA)))
if (is.na(n_workers) || n_workers < 1L) n_workers <- 3L
Sys.setenv(TARGETS_WORKERS = n_workers)
say("dispatching tar_make_future(workers = %d) over the tier-1 closure", n_workers)

t1 <- Sys.time()
ok <- tryCatch({
  targets::tar_make_future(names = tidyselect::any_of(!!closure),
                           workers = n_workers, reporter = "summary")
  TRUE
}, error = function(e) { say("FAILED: %s", conditionMessage(e)); FALSE })
say("finished after %.2f h, success = %s",
    as.numeric(difftime(Sys.time(), t1, units = "hours")), ok)

# Report failures unconditionally: TARGETS_ERROR_MODE="continue" lets a run
# finish "successfully" with individual targets errored, so gating this on !ok
# would hide exactly the errors that mode exists to surface.
{
  e <- tryCatch(targets::tar_errored(), error = function(e) character(0))
  e <- intersect(e, closure)
  say("errored targets in the closure: %d", length(e))
  if (length(e)) {
    print(utils::head(e, 20))
    m <- tryCatch(targets::tar_meta(fields = "error"), error = function(e) NULL)
    if (!is.null(m)) {
      m <- m[m$name %in% e & !is.na(m$error), c("name", "error")]
      if (nrow(m)) print(utils::head(m, 10))
    }
  }
}
