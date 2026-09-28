# =============================================================================
# explore/scripts/00_harness_check.R
#
# THE GATE. No probe result in this folder is trustworthy until the harness
# reproduces the numbers already on the record. This script runs the protocol's
# own arms through explore/R/harness.R and compares, cell by cell, against
# results/tables/protocol_v2/benchmarks_v2_cells.csv.
#
# Reference medians in that table (2026-09-27 23:44 run):
#   infill  level  domain_index  0.394   spatial 0.457   null -0.199   (18 cells)
#   infill  prev   domain_index  0.298   spatial 0.307   null -0.195   (18 cells)
#   country level  domain_index  0.287                                 (22 cells)
#   country prev   domain_index  0.192                                 (22 cells)
#
#   Rscript explore/scripts/00_harness_check.R
# -> explore/out/00_harness_check.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
message("loaded: ", length(E$PREDS), " predictors, ",
        nrow(exp_cell_index(E)), " cells")

ARMS <- ARMS_V2[c("null_train_mean", "spatial", "domain_index")]

# ── estimands A and B ───────────────────────────────────────────────────────
rows <- list()
ix <- exp_cell_index(E)
for (i in seq_len(nrow(ix))) {
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tgt),
                     error = function(e) NULL)
    if (is.null(cell)) next
    rows[[paste(i, tgt, "A")]] <- exp_infill(cell, ARMS, reps = REPS)
    rows[[paste(i, tgt, "B")]] <- exp_region(cell, ARMS)
  }
  message("  ", ix$country[i], " ", ix$outcome[i])
}

# ── estimand C ──────────────────────────────────────────────────────────────
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tgt, outcomes = on)
    if (length(cl) < 3) next
    rows[[paste("C", tgt, on)]] <-
      exp_loco(cl, ARMS_V2[c("null_train_mean", "domain_index")],
               domain_of = E$domain_of)
  }
  message("  LOCO ", tgt, " done")
}

raw <- dplyr::bind_rows(rows)
sm  <- exp_summarise(raw)
exp_write(sm, "00_harness_check")

# ── compare against the record ──────────────────────────────────────────────
REF <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/benchmarks_v2_cells.csv"),
                stringsAsFactors = FALSE)

cmp <- merge(sm[, c("country", "outcome", "target", "estimand", "arm", "spearman")],
             REF[, c("country", "outcome", "target", "estimand", "arm", "spearman")],
             by = c("country", "outcome", "target", "estimand", "arm"),
             suffixes = c("_exp", "_ref"))

cat("\n== harness vs record: median Spearman ==\n")
agg <- aggregate(cbind(spearman_exp, spearman_ref) ~ estimand + target + arm,
                 data = cmp, FUN = function(z) round(median(z, na.rm = TRUE), 3))
agg$n <- aggregate(spearman_exp ~ estimand + target + arm, data = cmp, FUN = length)$spearman_exp
print(agg[order(agg$estimand, agg$target, agg$arm), ], row.names = FALSE)

cmp$diff <- cmp$spearman_exp - cmp$spearman_ref
cat("\n== per-cell absolute difference ==\n")
cat("  max |diff| :", round(max(abs(cmp$diff), na.rm = TRUE), 4), "\n")
cat("  mean|diff| :", round(mean(abs(cmp$diff), na.rm = TRUE), 4), "\n")
cat("  cells      :", nrow(cmp), "\n")

worst <- cmp[order(-abs(cmp$diff)), ][1:10, ]
cat("\n== 10 largest discrepancies ==\n")
print(worst[, c("country", "outcome", "target", "estimand", "arm",
                "spearman_exp", "spearman_ref", "diff")], row.names = FALSE)

exp_write(cmp, "00_harness_check_vs_record")

ok <- mean(abs(cmp$diff) < 0.02, na.rm = TRUE)
cat("\nGATE: share of cells within 0.02 of the record:", round(ok, 3), "\n")
if (ok < 0.95) cat("GATE NOT CLEAN - investigate before trusting any probe.\n") else
  cat("GATE CLEAN.\n")
