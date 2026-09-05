# =============================================================================
# dashboard/data-raw/04_build_survey_planning.R
#
# Builds the bundle behind the "How much survey do you need?" tab, which
# replaces the old "Scenarios" tab.
#
# WHAT THE OLD TAB DID WRONG
# --------------------------
# "Scenarios" let a user move sliders over model settings and watch a number
# change. That is not a decision anyone in this project's audience makes. The
# decision they DO make is a budget one: how many regions to visit and how many
# clusters to field in each. This tab answers that question and nothing else.
#
# THE SOURCE, AND WHY IT IS TRUSTWORTHY
# -------------------------------------
# scripts/accuracy_impact/ws5_anchoring_budget.R. Crucially, its regional mean
# is JACKKNIFED - a district's regional anchor is computed from the region's
# OTHER districts' retained clusters. That matters because the headline
# anchoring result elsewhere in this project was WITHDRAWN for exactly the leak
# this script avoids: an anchor built from all of a region's respondents
# includes the scored district's own answer. WS5 was written after that
# withdrawal and does not repeat it.
#
# Whole CLUSTERS are dropped rather than individuals, because a survey planner
# buys clusters, not people.
#
# Scored against the full survey's district estimates - the best available
# stand-in for truth.
#
# THE HONESTY THIS TAB HAS TO CARRY
# ---------------------------------
# The curve is NOT monotone in region coverage. Configurations that sample a
# third or two-thirds of regions carry biases of -10 to -12 pp against -1 to -3
# pp elsewhere. Those settings arise only in the six-region countries, where the
# fallback national mean is estimated from very few regions, so they are most
# likely Monte Carlo noise plus a small-country artefact rather than a real
# property of the design. The tab flags them rather than smoothing them away,
# because a grant document that quotes a point off this curve should know which
# points are soft.
#
#   Rscript dashboard/data-raw/04_build_survey_planning.R
# -> dashboard/data/survey_planning.rds
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
setwd(here::here())
DASH <- here("dashboard", "data")
`%||%` <- function(a, b) if (is.null(a)) b else a

# Prefer the DENSE run: 6 cluster fractions x 9 region shares at 40 replicates,
# against 4 x 4 at 25 for the original. A planning curve read at a budget the
# grid never visited is an interpolation nobody checked.
f_dense  <- here("results", "tables", "anchoring_design_curve_DENSE.csv")
f_sparse <- here("results", "tables", "anchoring_design_curve.csv")
src <- if (file.exists(f_dense)) f_dense else f_sparse
curve <- read.csv(src, stringsAsFactors = FALSE)
cat("source:", basename(src), "-", nrow(curve), "rows\n")

# Round the realised region share back onto the grid the script SAMPLED
# (seq(0.2, 1.0, by = 0.1)). Countries have different region counts - 6 in
# Gambia and Sierra Leone against 27 in Malawi - so "sample 20% of regions"
# realises as 0.1852 in one country and 0.1875 in another. Left unrounded those
# become separate grid points, fragmenting a 9 x 6 design into 126 near-
# duplicates and making the pooled median an average over one or two countries
# rather than all four.
curve <- curve |>
  mutate(region_share = round(n_regions_anchored / n_regions_total, 1),
         country = ifelse(country == "Sierra Leone", "SierraLeone", country))

# ── per-cell grid: median and inter-replicate spread ─────────────────────────
# Median, not mean: a replicate that happens to drop every cluster in a small
# region produces an MAE an order of magnitude out, and the mean chases it.
cell <- curve |>
  group_by(country, outcome, region_share, fraction_clusters) |>
  summarise(
    reps        = dplyr::n(),
    pct_survey  = median(pct_survey_used, na.rm = TRUE),
    mae_a2      = median(mae_admin2_pp,   na.rm = TRUE),
    mae_a2_lo   = quantile(mae_admin2_pp, 0.25, na.rm = TRUE),
    mae_a2_hi   = quantile(mae_admin2_pp, 0.75, na.rm = TRUE),
    bias_a2     = median(bias_admin2_pp,  na.rm = TRUE),
    mae_a1      = median(mae_admin1_pp,   na.rm = TRUE),
    n_regions   = max(n_regions_total),
    .groups = "drop")

# ── pooled curve across cells ───────────────────────────────────────────────
pooled <- cell |>
  group_by(region_share, fraction_clusters) |>
  summarise(cells      = dplyr::n(),
            pct_survey = round(median(pct_survey, na.rm = TRUE), 1),
            mae_a2     = round(median(mae_a2,  na.rm = TRUE), 2),
            mae_a2_lo  = round(median(mae_a2_lo, na.rm = TRUE), 2),
            mae_a2_hi  = round(median(mae_a2_hi, na.rm = TRUE), 2),
            bias_a2    = round(median(bias_a2, na.rm = TRUE), 2),
            mae_a1     = round(median(mae_a1,  na.rm = TRUE), 2),
            .groups = "drop") |>
  # A configuration is flagged soft when its median bias is large relative to
  # the rest of the grid. See the header: these are the six-region artefacts.
  mutate(soft = abs(bias_a2) > 5) |>
  arrange(region_share, fraction_clusters)

# ── the full-survey reference, so every number has something to be read against
full_ref <- pooled |> filter(region_share == 1, fraction_clusters == 1)
baseline_mae <- if (nrow(full_ref)) full_ref$mae_a2[1] else NA_real_

saveRDS(list(
  cell = as.data.frame(cell),
  pooled = as.data.frame(pooled),
  baseline_mae = baseline_mae,
  source = basename(src),
  n_reps = max(cell$reps, na.rm = TRUE),
  note = paste(
    "Regional means are JACKKNIFED: a district's anchor uses the region's",
    "other districts only. The un-jackknifed anchoring result reported",
    "elsewhere in this project was withdrawn for that leak; this analysis was",
    "written afterwards and does not repeat it. Whole clusters are dropped,",
    "not individuals, because a planner buys clusters. Configurations marked",
    "'soft' carry large bias and arise only in six-region countries; treat",
    "them as noise, not as design guidance."),
  build_time = Sys.time()),
  file.path(DASH, "survey_planning.rds"))

cat(sprintf("-- Survey planning --\n  cells %d | grid points %d | reps %d | full-survey MAE %.2f pp\n",
            nrow(cell), nrow(pooled), max(cell$reps), baseline_mae))
print(as.data.frame(pooled |> filter(fraction_clusters %in% c(0.25, 0.6, 1))), row.names = FALSE)
