# =============================================================================
# scripts/accuracy_impact/ws9_regression_gate.R
#
# WS9. The regression gate.
#
# Guardrail 4 says no existing table under results/tables/ may be overwritten:
# new analyses write new files, or new rows carrying a `scheme` column. This
# script checks that claim rather than asserting it, by comparing every file in
# the WS0 freeze against its live counterpart by md5.
#
# It also lists every table this branch ADDED, so a reviewer can see the full
# surface of new results in one place.
#
# Classification, following results/tables/frozen_2026-08/REGRESSION_GATE.csv:
#   unchanged   live md5 equals frozen md5
#   CHANGED     live md5 differs; every instance must be explained line by line
#   missing     the live file is gone
#
#   Rscript scripts/accuracy_impact/ws9_regression_gate.R
# -> results/tables/frozen_2026-09/REGRESSION_GATE.csv
# =============================================================================
suppressPackageStartupMessages({library(here); library(tools)})
FROZEN <- here("results", "tables", "frozen_2026-09")
man <- utils::read.csv(file.path(FROZEN, "MANIFEST.csv"), stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(man))) {
  src <- here(man$source_path[i])
  live <- if (file.exists(src)) unname(tools::md5sum(src)) else NA_character_
  status <- if (is.na(live)) "missing" else
    if (identical(live, man$md5[i])) "unchanged" else "CHANGED"
  rows[[i]] <- data.frame(
    source_path = man$source_path[i], group = man$group[i], status = status,
    detail = "", frozen_md5 = man$md5[i], live_md5 = live %||% NA_character_,
    checked_on = as.character(Sys.Date()), stringsAsFactors = FALSE)
}
gate <- do.call(rbind, rows)
utils::write.csv(gate, file.path(FROZEN, "REGRESSION_GATE.csv"), row.names = FALSE)

cat("=== WS9 regression gate against the WS0 freeze ===\n")
print(as.data.frame(table(gate$status)), row.names = FALSE)
if (any(gate$status == "CHANGED")) {
  cat("\nCHANGED (each of these must be explained in the WS9 findings):\n")
  print(as.data.frame(gate[gate$status == "CHANGED", c("source_path","group")]),
        row.names = FALSE)
}
if (any(gate$status == "missing")) {
  cat("\nMISSING:\n")
  print(as.data.frame(gate[gate$status == "missing", c("source_path")]), row.names = FALSE)
}

# The other half of the gate: what this branch added.
base <- tryCatch(system2("git", c("show", "--name-only", "--pretty=format:",
                                  "f8b23d9:"), stdout = TRUE), error = function(e) NULL)
live_tables <- list.files(here("results", "tables"), pattern = "[.]csv$",
                          recursive = TRUE, full.names = FALSE)
tracked <- tryCatch(system2("git", c("ls-tree", "-r", "--name-only", "f8b23d9",
                                     "results/tables/"), stdout = TRUE),
                    error = function(e) character(0))
tracked <- sub("^results/tables/", "", tracked)
added <- setdiff(live_tables, tracked)
added <- added[!grepl("^frozen_2026-0", added)]
# Files git IGNORES are not additions by this branch. Seven transportability_*
# tables sit in results/tables/ and in .gitignore, dated before this branch was
# cut; a naive comparison against the git tree counts them as new. Ask git.
if (length(added)) {
  ig <- suppressWarnings(system2(
    "git", c("check-ignore", shQuote(file.path("results/tables", added))),
    stdout = TRUE, stderr = FALSE))
  ig <- sub("^results/tables/", "", ig)
  added <- setdiff(added, ig)
}
cat(sprintf("\n=== tables added by this branch: %d ===\n", length(added)))
for (a in sort(added)) cat("  ", a, "\n")
