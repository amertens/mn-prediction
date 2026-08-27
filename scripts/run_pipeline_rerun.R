# =============================================================================
# scripts/run_pipeline_rerun.R
#
# Restart the full targets rebuild, detached from any interactive session.
#
# WHY A SCRIPT RATHER THAN `Rscript -e`. This job runs for hours. Started from a
# shell that is tied to an editor/agent session it gets killed when that session
# exits -- which is exactly what happened on 2026-08-20, when the run reached
# only 19 of 895 targets before being terminated and left the store half
# rebuilt. Quoting a `reporter=` argument through PowerShell's argument parser
# also mangles it, so the call lives in a file instead.
#
# Launch (Windows PowerShell), detached and surviving session exit:
#
#   Start-Process -FilePath "C:\Program Files\R\R-4.4.2\bin\x64\Rscript.exe" `
#     -ArgumentList "scripts/run_pipeline_rerun.R" `
#     -RedirectStandardOutput "pipeline_rerun.log" `
#     -RedirectStandardError  "pipeline_rerun.err" `
#     -WindowStyle Hidden
#
# Or from any plain terminal:  Rscript scripts/run_pipeline_rerun.R
#
# WHAT THIS REBUILD INCLUDES (all invalidate large parts of the graph):
#   * Tanzania's 15 newly fetched GEE rasters in data/Tanzania_GEE_rasters/,
#     which take that country from 135 to ~176 gee_ columns and restore six
#     covariate families to the cross-country intersection FOR EVERY COUNTRY
#   * water-body handling: lake polygons dropped from the covariate side,
#     displaced lakeside clusters snapped back to the nearest land district
#     rather than discarded (87 rows in Malawi, 175 in Tanzania)
#   * duplicate Admin-2 key collapsing on every code path
#   * Fay-Herriot / BYM2 degenerate sampling-variance guard
#   * national benchmarking of the area-level predictions (now default on)
#   * r_max / r_share / signal columns, and null + spatial-only benchmark rows
# =============================================================================

t0 <- Sys.time()
cat(sprintf("[rerun] started %s\n", format(t0)))
cat(sprintf("[rerun] PIPELINE_MODE = %s\n",
            Sys.getenv("PIPELINE_MODE", "fast (default)")))

n_tif <- length(list.files("data/Tanzania_GEE_rasters", pattern = "\\.tif$"))
cat(sprintf("[rerun] Tanzania rasters on disk: %d\n", n_tif))
if (n_tif == 0)
  warning("data/Tanzania_GEE_rasters/ is empty - Tanzania will fall back to the ",
          "legacy-parity CSV and the covariate intersection will stay small")

ok <- tryCatch({
  targets::tar_make(reporter = "summary")
  TRUE
}, error = function(e) {
  cat(sprintf("[rerun] FAILED: %s\n", conditionMessage(e)))
  FALSE
})

el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("[rerun] finished %s (%.1f min), success = %s\n",
            format(Sys.time()), el, ok))

# A completed rebuild leaves nothing outdated; anything listed here failed or
# was skipped, and is worth seeing in the log rather than discovering later.
left <- tryCatch(targets::tar_outdated(reporter = "silent"),
                 error = function(e) character(0))
cat(sprintf("[rerun] targets still outdated: %d\n", length(left)))
if (length(left)) print(utils::head(left, 40))
