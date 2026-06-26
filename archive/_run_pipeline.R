# One-off rebuild runner: brings the _targets_full store in line with current code.
Sys.setenv(PIPELINE_MODE = "full")
# Per-method wall-clock cap (seconds) for the area benchmark runner — bounds
# pathologically slow methods (forward/stacked/two_stage) on rare outcomes
# (e.g. women_vitA) while leaving legit methods (bym2 ~52s) intact. Inherited by
# the per-target callr workers.
Sys.setenv(BENCH_METHOD_TIMEOUT = "180")
suppressWarnings(suppressMessages(library(targets)))
# Clear any stale lock left by a previously-killed run.
tryCatch(tar_unlock(store = "_targets_full"), error = function(e) NULL)
cat("=== tar_make rebuild starting ===\n")
t0 <- proc.time()[3]
ok <- tryCatch({
  # Prefer parallel future backend (pipeline's recommended invocation).
  if (requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.callr", quietly = TRUE)) {
    tar_make_future(workers = 4)
  } else {
    cat("future/future.callr unavailable — falling back to serial tar_make()\n")
    tar_make()
  }
  TRUE
}, error = function(e) {
  cat("tar_make_future failed:", conditionMessage(e), "\n")
  cat("Retrying serially with tar_make() ...\n")
  tryCatch({ tar_make(); TRUE },
           error = function(e2) { cat("tar_make ERROR:", conditionMessage(e2), "\n"); FALSE })
})
cat(sprintf("=== rebuild %s in %.1f min ===\n",
            if (isTRUE(ok)) "COMPLETE" else "FAILED", (proc.time()[3] - t0)/60))
# Post-run staleness check
od <- tryCatch(length(tar_outdated(store = "_targets_full")), error = function(e) NA)
cat(sprintf("Outdated targets remaining after rebuild: %s\n", od))
