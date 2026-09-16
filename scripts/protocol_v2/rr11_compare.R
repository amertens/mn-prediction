# =============================================================================
# scripts/protocol_v2/rr11_compare.R   [RR-11, 2026-09-15]
#
# Before / after table for the RR-11 re-run: the RR-10 results copied to
# results/tables/protocol_v2_pre_RR11_<date>/ against the new headline and
# the three other transport arms (_withdhs, _open, _nodhs). Mean Spearman over
# cells (the sandbox log's convention), by estimand x target x arm.
#
#   Rscript -e "source('scripts/protocol_v2/rr11_compare.R')"
# -> results/tables/protocol_v2/rr11_before_after.csv (+ printed table)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
P2 <- "results/tables/protocol_v2"
bk <- sort(list.dirs("results/tables", recursive = FALSE), decreasing = TRUE); bk <- bk[grepl("protocol_v2_pre_RR11_", bk)][1]
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
runs <- list(RR10 = rd(file.path(bk, "benchmarks_v2_summary.csv")),
             headline = rd(file.path(P2, "benchmarks_v2_summary.csv")),
             with_dhs = rd(file.path(P2, "benchmarks_v2_summary_withdhs.csv")),
             open_only = rd(file.path(P2, "benchmarks_v2_summary_open.csv")))
runs <- runs[!vapply(runs, is.null, NA)]
L <- bind_rows(lapply(names(runs), function(n) runs[[n]] |> mutate(run = n)))
# headline and open_only share the no-DHS shards; with_dhs has the every-tier shards - shown for the in-country rows too
W <- L |> filter(!(estimand != "country" & run == "open_only")) |>
  select(run, estimand, target, arm, cells, mean_spearman, cells_positive) |>
  mutate(val = sprintf("%.3f (%d/%d)", mean_spearman, cells_positive, cells)) |>
  select(-mean_spearman, -cells_positive, -cells) |>
  pivot_wider(names_from = run, values_from = val) |>
  arrange(factor(estimand, levels = c("infill", "region", "country")), target, arm)
write.csv(W, file.path(P2, "rr11_before_after.csv"), row.names = FALSE)
cat(sprintf("RR-10 backup: %s | arms present: %s\n\n", basename(bk), paste(names(runs), collapse = ", ")))
print(as.data.frame(W), row.names = FALSE, right = FALSE)
