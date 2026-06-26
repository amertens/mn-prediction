# Full-scope corrected pipeline run (called by run_corrected_then_hiv.ps1).
# Sequential tar_make; resume-on-restart via the _targets_corrected/ cache.
Sys.setenv(CORRECTED_SCOPE = "full")
library(targets)
cat(sprintf("[run] start %s | CORRECTED_SCOPE=full\n", format(Sys.time())))
tar_make(script = "_targets_corrected.R", store = "_targets_corrected",
         reporter = "summary")
cat(sprintf("[run] tar_make returned %s\n", format(Sys.time())))
