# =============================================================================
# scripts/protocol_v2/rr14_outcome_datasets.R   [RR-14, 2026-09-27]
# Step 1 of rerun_survey_fixes_2026-09-27.sh: rebuild the pipeline's outcome
# datasets for the four protocol countries in the full-mode store, so the
# 15 September survey-report fixes reach protocol-v2 script 01.
#
# Only the targets that carry those fixes are rebuilt: merged_raw / fsec / ext /
# merged_<country> (the Admin1 + Admin2 join that removed Malawi's repeated
# Traditional Authority names) and outcome_data_<country>_<outcome> (outcome
# variables, the population mask, the vitamin A rule). shortcut = TRUE takes
# every other upstream target from the store as it is, in particular the GEE
# zonal extraction (gee_admin2_*), which protocol v2 does not read (its
# predictors come from data/covariates/harmonized/) and whose Sierra Leone run
# failed on a memory allocation on 27 September.
#
# (A multi-line Rscript -e segfaults under Git Bash on Windows, hence a file;
# tar_make evaluates `names` in its own process, hence the inlined vector.)
#   PIPELINE_MODE=full Rscript scripts/protocol_v2/rr14_outcome_datasets.R
# =============================================================================
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (Sys.getenv("PIPELINE_MODE") != "full") stop("run with PIPELINE_MODE=full (the _targets_full store is full mode)")
suppressPackageStartupMessages(library(targets))
all_names <- tar_manifest(fields = "name")$name
cc <- "(gambia|ghana|malawi|sierraleone)"
merge_steps <- grep(paste0("^(merged_raw_|merged_fsec_|merged_ext_|merged_)", cc, "$"), all_names, value = TRUE)
outcomes <- grep(paste0("^outcome_data_", cc, "_"), all_names, value = TRUE)
nm <- c(merge_steps, outcomes)
cat(length(merge_steps), "merge steps and", length(outcomes), "outcome datasets\n"); print(merge_steps)
stopifnot(length(merge_steps) == 16, length(outcomes) >= 24)
eval(bquote(tar_make(names = tidyselect::all_of(.(nm)), shortcut = TRUE, reporter = "timestamp")))
m <- tar_meta(names = tidyselect::all_of(nm), fields = c("name", "time", "error"))
print(as.data.frame(m))
if (any(!is.na(m$error))) quit(status = 1)
