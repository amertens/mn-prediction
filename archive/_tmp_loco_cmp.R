suppressWarnings(suppressMessages({ library(targets); targets::tar_source("R/") }))
configs   <- get_country_configs()
countries <- c("Gambia", "Ghana", "SierraLeone", "Malawi")  # all 4 (Malawi GPS now in place)
outcomes  <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
methods   <- c("baseline", "penalized", "spatial_plus_soil", "coral")
out_csv   <- "results/cluster_vs_admin2_LOCO.csv"
dir.create("results", showWarnings = FALSE)

extract_country_cluster_gee <- function(cc, buffer_km = 2) {
  d <- readRDS(cc$data_path)
  if (is.list(d) && !is.data.frame(d) && is.data.frame(d$data)) d <- d$data
  ids <- unique(as.character(d[[cc$cluster_id]]))
  gps <- load_cluster_gps(cc, merged_cluster_ids = ids)
  if (is.null(gps) || nrow(gps) < 5) return(NULL)
  clusters <- data.frame(cluster_id = gps$cluster_id, lon = gps$lon, lat = gps$lat,
                         Admin1 = NA_character_, Admin2 = NA_character_, stringsAsFactors = FALSE)
  extract_gee_cluster_buffers(cc, clusters, buffer_km = 2)
}

cat("== extracting cluster-buffer GEE (once per country) ==\n")
cgee <- list()
for (cn in countries) cgee[[cn]] <- extract_country_cluster_gee(configs[[cn]], 2)
gee_a2 <- setNames(lapply(countries, function(cn)
            targets::tar_read_raw(paste0("gee_admin2_", tolower(cn)))), countries)

run_loco <- function(pooled, gee_vars, cnames, otag) {
  tryCatch(run_area_benchmarks_loco(pooled, gee_vars, cnames, outcome_label = otag,
             methods = methods, model_types = "continuous", augment_features = FALSE),
           error = function(e) { cat("  LOCO err:", conditionMessage(e), "\n"); NULL })
}

rows <- list()
for (otag in outcomes) {
  cat(sprintf("\n== outcome: %s ==\n", otag))
  # ---- cluster-resolution pooled LOCO ----
  cl_frames <- list()
  for (cn in countries) {
    oc <- configs[[cn]]$outcomes[[otag]]; if (is.null(oc)) next
    cl <- tryCatch(aggregate_outcome_to_clusters(configs[[cn]], oc), error = function(e) NULL)
    if (is.null(cl) || nrow(cl) < 8) { cat("  ", cn, "cluster: too few\n"); next }
    gc <- grep("^gee_", names(cgee[[cn]]), value = TRUE)
    m <- merge(cl, cgee[[cn]][, c("cluster_id", gc)], by = "cluster_id")
    m$country <- cn
    cl_frames[[cn]] <- m
  }
  if (length(cl_frames) >= 2) {
    common <- Reduce(intersect, lapply(cl_frames, function(d) grep("^gee_", names(d), value = TRUE)))
    keep <- c("country", "Admin2", "cluster_id", "svy_prev", "n_svy", "lon", "lat", "Admin1", common)
    pooled_cl <- do.call(rbind, lapply(cl_frames, function(d) d[, intersect(keep, names(d))]))
    rc <- run_loco(pooled_cl, common, unique(pooled_cl$country), otag)
    if (!is.null(rc) && nrow(rc) > 0) { rc$resolution <- "cluster"; rc$outcome <- otag; rows[[paste0(otag, "_cl")]] <- rc }
    cat(sprintf("  cluster pooled: %d units, %d common gee, %d countries\n",
                nrow(pooled_cl), length(common), length(cl_frames)))
  }
  # ---- admin-2 pooled LOCO (same 3 countries) ----
  svy_list <- setNames(lapply(countries, function(cn)
                tryCatch(targets::tar_read_raw(paste0("svy_admin2_", tolower(cn), "_", otag)),
                         error = function(e) NULL)), countries)
  pa2 <- tryCatch(build_area_loco_dataset(svy_list, gee_a2), error = function(e) NULL)
  if (!is.null(pa2) && nrow(pa2$pooled_data) > 10) {
    ra <- run_loco(pa2$pooled_data, pa2$common_gee_vars, pa2$country_names, otag)
    if (!is.null(ra) && nrow(ra) > 0) { ra$resolution <- "admin2"; ra$outcome <- otag; rows[[paste0(otag, "_a2")]] <- ra }
  }
  flat <- do.call(rbind, lapply(rows, function(d)
            d[, c("resolution", "outcome", "held_out", "method", "pearson_r", "n_test")]))
  write.csv(flat, out_csv, row.names = FALSE)
}

all <- do.call(rbind, lapply(rows, function(d)
         d[, c("resolution", "outcome", "held_out", "method", "pearson_r")]))
cat("\n=== Mean LOCO pearson_r by method x resolution (3 countries) ===\n")
agg <- aggregate(pearson_r ~ method + resolution, all, function(x) round(mean(x, na.rm = TRUE), 3))
print(reshape(agg, idvar = "method", timevar = "resolution", direction = "wide"), row.names = FALSE)
cat("\n=== CORAL per split: admin2 vs cluster ===\n")
cr <- all[all$method == "coral", ]
print(reshape(cr[, c("outcome", "held_out", "resolution", "pearson_r")],
              idvar = c("outcome", "held_out"), timevar = "resolution", direction = "wide"),
      row.names = FALSE)
cat("\nwritten:", out_csv, "\n")
