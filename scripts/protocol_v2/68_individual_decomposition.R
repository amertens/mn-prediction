# =============================================================================
# scripts/protocol_v2/68_individual_decomposition.R   [IL-02b]
#
# FROM JANUARY'S CROSS-VALIDATION TO IL-01, ONE FACTOR AT A TIME
#
# Survey-variable arms (S) and proxy arms (P), five Ghana outcomes, five fold
# draws each, on a PSOCK cluster. Variants (see VARIANTS in the lib):
#   S1  January survey set, 10-fold cluster-blocked, prescreen on ALL rows
#       (January's own CV; reproduces the saved fits' CV numbers)
#   S2  + prescreen inside the training folds
#   S3  + identifiers, sampling-design / region / team columns, dates and
#       blood-derived columns (RDT, genotypes) removed
#   S4b January set, 10-fold district-blocked
#   S4  clean set, 10-fold district-blocked
#   S5  clean set, 5-fold district-blocked (IL-01's folds)
#   S6  IL-01's own column rule (regex + 80-column coverage cap), 5-fold district
#   S7  clean set, 5-fold district, no prescreen
#   P1-P4 the same ladder for January's proxy set (region and month removed at P3)
# Metrics per fit: pooled out-of-fold AUC (IL-01's convention), within-fold AUC
# (pairs inside one fold only, which removes the fold-to-fold offset of the
# training mean that pushes a no-skill pooled AUC below 0.5), Brier skill
# against the full-sample prevalence (IL-01 and January) and against the
# out-of-fold training mean.
#
#   NW=17 REPS=5 Rscript scripts/protocol_v2/68_individual_decomposition.R
# -> results/tables/protocol_v2/il02_decomposition_raw.csv   (one row per fit x learner)
# -> results/tables/protocol_v2/il02_decomposition_cells.csv (stack, mean over draws)
# -> results/tables/protocol_v2/il02_decomposition_avg.csv   (mean over outcomes)
# =============================================================================
suppressPackageStartupMessages(library(parallel))
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
LIB <- normalizePath("scripts/protocol_v2/68_individual_decomposition_lib.R", winslash = "/"); OUTDIR <- "results/tables/protocol_v2"
NW <- as.integer(Sys.getenv("NW", "16")); REPS <- as.integer(Sys.getenv("REPS", "3"))
ONLY <- Sys.getenv("ONLY", "")
source(LIB)
cat(sprintf("survey_jan %d cols | survey_clean %d | proxy_jan %d | proxy_clean %d\n",
            length(SET$survey_jan), length(SET$survey_clean), length(SET$proxy_jan), length(SET$proxy_clean)))
cat("dropped from survey_jan as ids/design/blood:\n"); print(setdiff(SET$survey_jan, SET$survey_clean))
tasks <- expand.grid(variant = names(VARIANTS), outcome = names(OUTC), rep = seq_len(REPS), stringsAsFactors = FALSE)
if (nzchar(ONLY)) tasks <- tasks[tasks$variant %in% strsplit(ONLY, ",")[[1]], ]
tasks <- tasks[order(!grepl("^P", tasks$variant)), ]           # heavy proxy tasks first
cl <- makePSOCKcluster(NW); on.exit(stopCluster(cl))
clusterExport(cl, c("IL", "tasks")); invisible(clusterEvalQ(cl, source(LIB)))
t0 <- Sys.time()
res <- parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  g <- tasks[i, ]; t1 <- Sys.time()
  out <- tryCatch(run_task(g$variant, g$outcome, g$rep), error = function(e) data.frame(variant = g$variant, outcome = g$outcome, rep = g$rep, error = conditionMessage(e)))
  out$secs <- as.numeric(difftime(Sys.time(), t1, units = "secs")); out })
R <- dplyr::bind_rows(res)
write.csv(R, file.path(OUTDIR, paste0("il02_decomposition_raw", if (nzchar(ONLY)) paste0("_", gsub(",", "", ONLY)) else "", ".csv")), row.names = FALSE)
cat(sprintf("done %d tasks in %.1f min\n", nrow(tasks), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if ("error" %in% names(R)) print(R[!is.na(R$error), c("variant", "outcome", "rep", "error")])

# ---- summary ----
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
if ("error" %in% names(R)) { cat("errors:\n"); print(unique(R[!is.na(R$error), c("variant", "outcome", "error")])); R <- R[is.na(R$error), ] }
LAB <- c(S1 = "S1 January survey set, cluster folds, prescreen on all rows",
         S2 = "S2 + prescreen inside folds",
         S3 = "S3 + ids/design/blood columns removed",
         S4b = "S4b January survey set, district folds",
         S4 = "S4 clean set, district folds (10)",
         S5 = "S5 clean set, district folds (5, as IL-01)",
         S6 = "S6 IL-01 column rule, district folds (5)",
         S7 = "S7 clean set, district folds (5), no prescreen",
         P1 = "P1 January proxy set, cluster folds, prescreen on all rows",
         P2 = "P2 + prescreen inside folds",
         P3 = "P3 region/month removed, district folds (10)",
         P4 = "P4 region/month removed, district folds (5)")
st <- R |> filter(learner == "stack") |> group_by(variant, outcome) |>
  summarise(auc = mean(auc_pooled), auc_w = mean(auc_within), bss = mean(bss_national), bss_null = mean(bss_vs_trainmean),
            bss_sd = sd(bss_national), n_cols = median(n_cols_used_median), .groups = "drop")
options(width = 220)
cat("\n== Stack, mean over fold draws: pooled AUC / within-fold AUC / Brier skill vs national / vs training mean ==\n")
w <- st |> mutate(cell = sprintf("%.3f/%.3f %+.3f", auc, auc_w, bss)) |> select(variant, outcome, cell) |>
  pivot_wider(names_from = outcome, values_from = cell) |> mutate(label = LAB[variant]) |> arrange(match(variant, names(LAB)))
print(as.data.frame(w), row.names = FALSE)
cat("\n== Mean over the five outcomes ==\n")
avg <- st |> group_by(variant) |> summarise(auc = mean(auc), auc_within = mean(auc_w), bss_national = mean(bss), bss_vs_trainmean = mean(bss_null),
                                            fold_draw_sd_bss = mean(bss_sd), cols_used = median(n_cols), .groups = "drop") |>
  mutate(label = LAB[variant]) |> arrange(match(variant, names(LAB)))
print(as.data.frame(avg |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\n== Which learner carries it (mean over outcomes, Brier skill vs national) ==\n")
lw <- R |> group_by(variant, learner) |> summarise(bss = round(mean(bss_national), 3), .groups = "drop") |>
  pivot_wider(names_from = learner, values_from = bss) |> arrange(match(variant, names(LAB)))
print(as.data.frame(lw), row.names = FALSE)
if (!nzchar(ONLY)) { write.csv(st, file.path(OUTDIR, "il02_decomposition_cells.csv"), row.names = FALSE); write.csv(avg, file.path(OUTDIR, "il02_decomposition_avg.csv"), row.names = FALSE) }
