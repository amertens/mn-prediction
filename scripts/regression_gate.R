# =============================================================================
# scripts/regression_gate.R
#
# WS9. Compare every table that has a frozen baseline against that baseline and
# classify what changed.
#
# This script RECOMPUTES NOTHING. It reads results/tables/frozen_2026-08/
# MANIFEST.csv, finds each baseline's live counterpart at the recorded
# `source_path`, and reports one of:
#
#   unchanged        md5 of the live file equals the frozen md5
#   new_scheme_rows  every frozen row survives unchanged and the only difference
#                    is added rows carrying a `scheme` (or `variant`) value the
#                    baseline did not have. This is the expected shape of an
#                    additive result.
#   changed_baseline a value that existed in the baseline is now different, or a
#                    frozen row has disappeared, or the columns changed. Every
#                    one of these has to be explained line by line before the
#                    branch is finished.
#   missing_live     the baseline's source file is gone
#
# The distinction between the second and third categories is the whole point.
# A run that adds rows is additive; a run that moves an existing number is a
# regression until someone explains it.
#
# Run:
#   Rscript scripts/regression_gate.R
#
# Writes results/tables/frozen_2026-08/REGRESSION_GATE.csv and prints a summary.
# Exits 1 if any table is classified changed_baseline or missing_live, so the
# script can gate a release.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})

FROZEN_DIR <- here::here("results", "tables", "frozen_2026-08")
MANIFEST   <- file.path(FROZEN_DIR, "MANIFEST.csv")
OUT        <- file.path(FROZEN_DIR, "REGRESSION_GATE.csv")

# Columns that identify a result *scheme* rather than a measurement. A row whose
# scheme value is absent from the baseline is a new result, not a changed one.
SCHEME_COLS <- c("scheme", "variant", "approach", "method", "model", "strategy")

if (!file.exists(MANIFEST))
  stop("No baseline manifest. Run scripts/freeze_baselines.R first.", call. = FALSE)

manifest <- readr::read_csv(MANIFEST, show_col_types = FALSE)

# Read a table for comparison. RDS baselines are compared by md5 only, since a
# row-level diff of an arbitrary R object is not well defined.
.rg_read <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE))
    tryCatch(as.data.frame(readr::read_csv(path, show_col_types = FALSE,
                                           progress = FALSE)),
             error = function(e) NULL)
  else NULL
}

#' Which scheme columns does this pair of tables share?
.rg_scheme_cols <- function(a, b) intersect(SCHEME_COLS, intersect(names(a), names(b)))

#' Classify one baseline against its live counterpart.
.rg_classify <- function(frozen_abs, live_abs) {
  if (!file.exists(live_abs))
    return(list(status = "missing_live", detail = "live file does not exist"))

  live_md5   <- unname(tools::md5sum(live_abs))
  frozen_md5 <- unname(tools::md5sum(frozen_abs))
  if (identical(live_md5, frozen_md5))
    return(list(status = "unchanged", detail = ""))

  fz <- .rg_read(frozen_abs); lv <- .rg_read(live_abs)
  if (is.null(fz) || is.null(lv))
    return(list(status = "changed_baseline",
                detail = "content differs and is not row-comparable"))

  if (!identical(sort(names(fz)), sort(names(lv)))) {
    added   <- setdiff(names(lv), names(fz))
    removed <- setdiff(names(fz), names(lv))
    # Added columns alone are additive as long as every frozen column is intact.
    if (length(removed) == 0) {
      common <- names(fz)
      if (isTRUE(all.equal(fz[common], lv[seq_len(nrow(fz)), common, drop = FALSE],
                           check.attributes = FALSE)))
        return(list(status = "new_scheme_rows",
                    detail = paste("columns added:", paste(added, collapse = ", "))))
    }
    return(list(status = "changed_baseline",
                detail = sprintf("columns added: [%s]; removed: [%s]",
                                 paste(added, collapse = ", "),
                                 paste(removed, collapse = ", "))))
  }

  sc <- .rg_scheme_cols(fz, lv)
  if (length(sc)) {
    key_fz <- do.call(paste, c(fz[sc], sep = "\r"))
    key_lv <- do.call(paste, c(lv[sc], sep = "\r"))
    new_schemes <- setdiff(unique(key_lv), unique(key_fz))
    kept <- lv[!(key_lv %in% new_schemes), , drop = FALSE]
    if (length(new_schemes) &&
        nrow(kept) == nrow(fz) &&
        isTRUE(all.equal(fz, kept, check.attributes = FALSE)))
      return(list(status = "new_scheme_rows",
                  detail = sprintf("%d added scheme value(s) via [%s]; all %d baseline rows identical",
                                   length(new_schemes), paste(sc, collapse = "+"), nrow(fz))))
  }

  if (nrow(lv) > nrow(fz) &&
      isTRUE(all.equal(fz, lv[seq_len(nrow(fz)), , drop = FALSE],
                       check.attributes = FALSE)))
    return(list(status = "new_scheme_rows",
                detail = sprintf("%d rows appended, all %d baseline rows identical",
                                 nrow(lv) - nrow(fz), nrow(fz))))

  n_row_delta <- nrow(lv) - nrow(fz)
  detail <- if (n_row_delta != 0)
    sprintf("row count %d -> %d", nrow(fz), nrow(lv))
  else {
    diff_cols <- names(fz)[vapply(names(fz), function(cn)
      !isTRUE(all.equal(fz[[cn]], lv[[cn]], check.attributes = FALSE)), logical(1))]
    sprintf("same row count; %d column(s) differ: %s",
            length(diff_cols), paste(utils::head(diff_cols, 8), collapse = ", "))
  }
  list(status = "changed_baseline", detail = detail)
}

rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  frozen_abs <- here::here(manifest$frozen_path[i])
  live_abs   <- here::here(manifest$source_path[i])
  cl <- .rg_classify(frozen_abs, live_abs)
  rows[[i]] <- data.frame(
    source_path = manifest$source_path[i],
    group       = manifest$group[i],
    status      = cl$status,
    detail      = cl$detail,
    frozen_md5  = manifest$md5[i],
    live_md5    = if (file.exists(live_abs)) unname(tools::md5sum(live_abs)) else NA_character_,
    checked_on  = format(Sys.Date(), "%Y-%m-%d"),
    stringsAsFactors = FALSE
  )
}
gate <- do.call(rbind, rows)
gate <- gate[order(factor(gate$status,
                          levels = c("changed_baseline", "missing_live",
                                     "new_scheme_rows", "unchanged")),
                   gate$source_path), ]
readr::write_csv(gate, OUT)

cat("\n=== regression gate against results/tables/frozen_2026-08/ ===\n")
tab <- table(factor(gate$status, levels = c("unchanged", "new_scheme_rows",
                                            "changed_baseline", "missing_live")))
for (s in names(tab)) cat(sprintf("  %-18s %d\n", s, tab[[s]]))

flagged <- gate[gate$status %in% c("changed_baseline", "missing_live"), ]
if (nrow(flagged)) {
  cat("\nEach of these must be explained before the branch is finished:\n")
  for (i in seq_len(nrow(flagged)))
    cat(sprintf("  [%s] %s\n      %s\n", flagged$status[i],
                flagged$source_path[i], flagged$detail[i]))
}
cat("\nwrote ", OUT, "\n", sep = "")

if (nrow(flagged)) quit(status = 1L)
