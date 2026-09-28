# =============================================================================
# scripts/protocol_v2/69_index_person_level.R   [IL-02c]
# Person-level AUC and Brier score of the PCA domain index, with ceilings.
# See 69_index_person_level_lib.R for the methods.
#
#   NW=5 REPS=20 Rscript scripts/protocol_v2/69_index_person_level.R
# -> results/tables/protocol_v2/il02_index_person_raw.csv      one row per cell x scheme x draw x method
# -> results/tables/protocol_v2/il02_index_person_summary.csv  mean over draws
# =============================================================================
suppressPackageStartupMessages({library(parallel); library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
LIB <- normalizePath("scripts/protocol_v2/69_index_person_level_lib.R", winslash = "/"); OUTDIR <- "results/tables/protocol_v2"
NW <- as.integer(Sys.getenv("NW", "5"))
TG0 <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
MAIN <- c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12")
cells <- TG0 |> distinct(country, outcome) |> filter(outcome %in% MAIN)
cl <- makePSOCKcluster(NW)
clusterExport(cl, c("LIB", "cells")); invisible(clusterEvalQ(cl, source(LIB)))
out <- parLapplyLB(cl, seq_len(nrow(cells)), function(i) {
  r <- tryCatch(run_cell(cells$country[i], cells$outcome[i]), error = function(e) list(rows = NULL, chk = NULL, note = paste("ERROR:", conditionMessage(e))))
  r$cell <- paste(cells$country[i], cells$outcome[i]); r })
stopCluster(cl)
for (o in out) cat(sprintf("%-26s %s
", o$cell, o$note))
rows <- lapply(out, `[[`, "rows"); chk <- lapply(out, `[[`, "chk")
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "il02_index_person_raw.csv"), row.names = FALSE)
cat("
respondents reproduce the district targets (r between recomputed and targets_v2 prevalence):
")
print(bind_rows(chk), row.names = FALSE)
SUM <- R |> group_by(country, outcome, scheme, method) |>
  summarise(n = first(n), prev = first(prev), n_districts = first(n_districts), reps = n(),
            auc_pooled = mean(auc_pooled, na.rm = TRUE), auc_within = mean(auc_within, na.rm = TRUE),
            brier = mean(brier, na.rm = TRUE), bss_national = mean(bss_national, na.rm = TRUE),
            bss_vs_null = mean(bss_vs_null, na.rm = TRUE), district_spearman = mean(district_spearman, na.rm = TRUE), .groups = "drop")
write.csv(SUM, file.path(OUTDIR, "il02_index_person_summary.csv"), row.names = FALSE)
options(width = 220)
print(as.data.frame(SUM |> filter(scheme %in% c("kfold5", "none"), method %in% c("index_cal", "index_glm", "index_lvl_glm", "null",
  "oracle_district_insample", "oracle_district_loo_eb", "oracle_cluster_loo_eb")) |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
