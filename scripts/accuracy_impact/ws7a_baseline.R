# Generate / refresh the Admin-2 join lint baseline. Run deliberately, never
# automatically: adding a row here grandfathers a name-only join.
#   Rscript scripts/accuracy_impact/ws7a_baseline.R
suppressPackageStartupMessages(library(here))
source(here("R", "lint_admin2_joins.R"))
s <- scan_admin2_joins()
s$audit <- "not assessed"
readr::write_csv(s, here("tests", "testthat", "admin2_join_baseline.csv"))
cat(sprintf("baseline: %d sites across %d files\n", nrow(s), length(unique(s$file))))
print(as.data.frame(table(s$kind)))
cat("\n--- sites per file ---\n")
tf <- sort(table(s$file), decreasing = TRUE)
print(as.data.frame(tf))
