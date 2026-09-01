# Entry point for the structural guards added in WS7a.
#
# This repository is a targets pipeline rather than an R package, so there is no
# library() to load under test. The helper sources R/ directly.
#
#   Rscript tests/testthat.R
#
# Set R_USER and HOME to the OneDrive Documents directory so ~/.rdhs.json
# resolves, as docs/SESSION_FINDINGS_FOR_REVIEW.md section 12 requires.
library(testthat)
testthat::test_dir(file.path("tests", "testthat"), stop_on_failure = TRUE)
