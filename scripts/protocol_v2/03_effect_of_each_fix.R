# =============================================================================
# scripts/protocol_v2/03_effect_of_each_fix.R
#
# What did each of the five fixes actually buy? Each row of the output isolates
# one change by comparing two arms/targets/schemes that differ ONLY in that
# change, on identical cells.
#
#   FIX 1  replication   spread of a cell's score across fold draws, i.e. how
#                        much a single-draw report could have been off. This is
#                        the quantity that made the published median r 0.058.
#   FIX 1  precision     unweighted MAE vs n_eff-weighted MAE.
#   FIX 2  target        continuous biomarker vs binary prevalence, same arm.
#   FIX 4  dimension     18 domain scores vs 373 raw columns, same learner.
#   FIX 5  estimand      the same arm scored under in-fill, region
#                        extrapolation and country transport.
#
# FIX 3 (within-country rank-normalisation) cannot be isolated from these
# outputs because every arm here already uses it; 04_fix3_normalisation.R runs
# that comparison directly in the setting where it bites.
#
#   Rscript scripts/protocol_v2/03_effect_of_each_fix.R
# -> results/tables/protocol_v2/effect_of_each_fix.csv
# -> results/tables/protocol_v2/fold_draw_risk.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"

RAW   <- read.csv(file.path(OUTDIR, "benchmarks_v2_raw.csv"), stringsAsFactors = FALSE)
CELLS <- read.csv(file.path(OUTDIR, "benchmarks_v2_cells.csv"), stringsAsFactors = FALSE)

rows <- list()
paircmp <- function(a, b, key, label, metric = "spearman") {
  m <- inner_join(a, b, by = key, suffix = c("_a", "_b"))
  ma <- m[[paste0(metric, "_a")]]; mb <- m[[paste0(metric, "_b")]]
  ok <- is.finite(ma) & is.finite(mb)
  if (!any(ok)) return(NULL)
  data.frame(comparison = label, cells = sum(ok),
             mean_a = round(mean(ma[ok]), 3), mean_b = round(mean(mb[ok]), 3),
             median_a = round(median(ma[ok]), 3), median_b = round(median(mb[ok]), 3),
             mean_gain = round(mean(ma[ok] - mb[ok]), 3),
             a_better_in = sum(ma[ok] > mb[ok]), of = sum(ok))
}

# ── FIX 1a. How much could a single fold draw have moved a cell? ────────────
draw <- RAW |> filter(estimand == "infill", is.finite(spearman)) |>
  group_by(country, outcome, target, arm) |>
  summarise(reps = n(), mean_r = mean(spearman), sd_r = sd(spearman),
            min_r = min(spearman), max_r = max(spearman),
            range_r = max(spearman) - min(spearman), .groups = "drop")
write.csv(draw, file.path(OUTDIR, "fold_draw_risk.csv"), row.names = FALSE)
dr <- draw |> filter(!arm %in% c("null_train_mean"))
rows[["fix1_draw"]] <- data.frame(
  comparison = "FIX 1 replication: spread of one cell's r across fold draws",
  cells = nrow(dr), mean_a = round(mean(dr$sd_r, na.rm = TRUE), 3),
  mean_b = round(mean(dr$range_r, na.rm = TRUE), 3),
  median_a = round(median(dr$sd_r, na.rm = TRUE), 3),
  median_b = round(median(dr$range_r, na.rm = TRUE), 3),
  mean_gain = NA_real_, a_better_in = NA_integer_, of = nrow(dr))

# ── FIX 1b. Precision weighting ────────────────────────────────────────────
pw <- CELLS |> filter(target == "prev", is.finite(mae), is.finite(wmae))
rows[["fix1_prec"]] <- data.frame(
  comparison = "FIX 1 precision: unweighted MAE vs n_eff-weighted MAE (pp)",
  cells = nrow(pw), mean_a = round(mean(pw$mae), 2), mean_b = round(mean(pw$wmae), 2),
  median_a = round(median(pw$mae), 2), median_b = round(median(pw$wmae), 2),
  mean_gain = round(mean(pw$mae - pw$wmae), 3),
  a_better_in = sum(pw$wmae < pw$mae), of = nrow(pw))

# ── FIX 2. Continuous target vs binary prevalence ──────────────────────────
for (es in unique(CELLS$estimand)) {
  a <- CELLS |> filter(estimand == es, target == "level")
  b <- CELLS |> filter(estimand == es, target == "prev")
  r <- paircmp(a, b, c("country", "outcome", "estimand", "arm"),
               paste0("FIX 2 target [", es, "]: continuous biomarker vs binary prevalence"))
  if (!is.null(r)) rows[[paste0("fix2_", es)]] <- r
}

# ── FIX 4. Domain scores vs raw columns (same learner) ─────────────────────
for (es in unique(CELLS$estimand)) {
  for (tg in unique(CELLS$target)) {
    a <- CELLS |> filter(estimand == es, target == tg, arm == "domain_enet")
    b <- CELLS |> filter(estimand == es, target == tg, arm == "raw_enet")
    r <- paircmp(a, b, c("country", "outcome"),
                 paste0("FIX 4 dimension [", es, "/", tg,
                        "]: 18 domain scores vs 373 raw columns, same learner"))
    if (!is.null(r)) rows[[paste0("fix4_", es, tg)]] <- r
    # and the zero-tuning index against the penalised raw fit
    a2 <- CELLS |> filter(estimand == es, target == tg, arm == "domain_index")
    r2 <- paircmp(a2, b, c("country", "outcome"),
                  paste0("FIX 4 dimension [", es, "/", tg,
                         "]: zero-tuning domain index vs 373-column elastic net"))
    if (!is.null(r2)) rows[[paste0("fix4b_", es, tg)]] <- r2
  }
}

# ── FIX 5. The same arm under different estimands ──────────────────────────
for (arm_i in c("domain_index", "spatial", "domain_enet")) {
  a <- CELLS |> filter(estimand == "infill", arm == arm_i)
  b <- CELLS |> filter(estimand == "region", arm == arm_i)
  r <- paircmp(a, b, c("country", "outcome", "target"),
               paste0("FIX 5 estimand [", arm_i,
                      "]: in-fill vs whole-region extrapolation"))
  if (!is.null(r)) rows[[paste0("fix5_", arm_i)]] <- r
}

# ── the comparison the withdrawn baseline used to win ──────────────────────
a <- CELLS |> filter(estimand == "infill", arm == "domain_index")
b <- CELLS |> filter(estimand == "infill", arm == "region_mean_jk")
r <- paircmp(a, b, c("country", "outcome", "target"),
             "FIX 5 baseline: covariates vs JACKKNIFED regional mean (in-fill)")
if (!is.null(r)) rows[["fix5_base"]] <- r
a <- CELLS |> filter(estimand == "infill", arm == "domain_index")
b <- CELLS |> filter(estimand == "infill", arm == "spatial")
r <- paircmp(a, b, c("country", "outcome", "target"),
             "covariates vs covariate-free spatial smoother, FOLD-MATCHED (in-fill)")
if (!is.null(r)) rows[["fix5_spatial"]] <- r

EF <- bind_rows(rows)
write.csv(EF, file.path(OUTDIR, "effect_of_each_fix.csv"), row.names = FALSE)
cat("\n================ EFFECT OF EACH FIX ================\n")
print(as.data.frame(EF), row.names = FALSE)

cat("\n--- worst single-draw risk: cells where one draw could have reported ---\n")
worst <- dr |> arrange(desc(range_r)) |> head(12) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(worst), row.names = FALSE)

cat("\nNOTE the null arm's correlation is not interpretable: a leave-fold-out\n")
cat("training mean is nearly constant, so its r is dominated by trivial\n")
cat("fold-to-fold variation and is mechanically negative. Read its MAE only.\n")
cat("\nDONE\n")
