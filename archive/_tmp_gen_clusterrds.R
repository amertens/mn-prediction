# Generate dashboard/data/cluster_resolution.rds from whatever cluster CSVs
# exist now (mirrors the hook in 01_prepare_dashboard_data.R). Lets the new
# Resolution panel render before the full prep script is re-run.
`%||%` <- function(a, b) if (is.null(a)) b else a
cmp_path  <- "results/cluster_vs_admin2_comparison.csv"
loco_path <- "results/cluster_vs_admin2_LOCO.csv"

cmp <- if (file.exists(cmp_path)) {
  d <- read.csv(cmp_path, stringsAsFactors = FALSE)
  d$delta_r   <- d$r_cluster   - d$r_admin2
  d$delta_mae <- d$mae_cluster - d$mae_admin2
  d
} else NULL
loco <- if (file.exists(loco_path)) read.csv(loco_path, stringsAsFactors = FALSE) else NULL

saveRDS(list(comparison = cmp, loco = loco, build_time = Sys.time()),
        "dashboard/data/cluster_resolution.rds")
cat(sprintf("cluster_resolution.rds: comparison rows=%d, LOCO rows=%d\n",
            nrow(cmp %||% data.frame()), nrow(loco %||% data.frame())))
