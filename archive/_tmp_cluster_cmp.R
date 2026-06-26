# Within-country admin2-vs-cluster resolution comparison, all four countries.
# Rebuilds results/cluster_vs_admin2_comparison.csv. Engine = fit_area_sl()
# (logit-transformed prevalence, LOO for n<=50 else 5-fold), run at admin-2
# resolution (gee_admin2 zonal means) and at survey-cluster resolution
# (2km GEE buffers). Malawi now joins the prior three via its new GPS file.
suppressWarnings(suppressMessages({ library(targets); targets::tar_source("R/") }))

configs      <- get_country_configs()
country_keys <- c("Gambia", "Ghana", "SierraLeone", "Malawi")
slugs        <- c(Gambia = "gambia", Ghana = "ghana",
                  SierraLeone = "sierraleone", Malawi = "malawi")
outcomes     <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
out_csv      <- "results/cluster_vs_admin2_comparison.csv"
OUTCOME_TYPE <- "binomial"   # fit_area_sl default: logit-transformed prevalence
SEED         <- 12345L

# ---- pre-extract cluster-buffer GEE once per country (the slow step) ----
cat("== extracting cluster-buffer GEE (once per country) ==\n")
cgee <- list()
for (ck in country_keys) {
  cc <- configs[[ck]]
  d  <- readRDS(cc$data_path)
  if (is.list(d) && !is.data.frame(d) && is.data.frame(d$data)) d <- d$data
  ids <- unique(as.character(d[[cc$cluster_id]]))
  gps <- load_cluster_gps(cc, merged_cluster_ids = ids)
  if (is.null(gps) || nrow(gps) < 5) { cat("  ", ck, ": no usable GPS\n"); next }
  clusters <- data.frame(cluster_id = gps$cluster_id, lon = gps$lon, lat = gps$lat,
                         Admin1 = NA_character_, Admin2 = NA_character_,
                         stringsAsFactors = FALSE)
  cgee[[ck]] <- extract_gee_cluster_buffers(cc, clusters, buffer_km = 2)
}

safe_fit <- function(svy, gee) {
  m <- tryCatch(fit_area_sl(svy, gee, outcome_type = OUTCOME_TYPE)$metrics,
                error = function(e) { cat("    fit err:", conditionMessage(e), "\n"); NULL })
  if (is.null(m) || is.null(m$n)) return(c(n = NA, r = NA, mae = NA))
  c(n = as.integer(m$n), r = m$pearson_r, mae = m$mae_pp)
}

rows <- list()
for (ck in country_keys) {
  cc <- configs[[ck]]; sl <- slugs[[ck]]
  for (ot in outcomes) {
    oc <- cc$outcomes[[ot]]; if (is.null(oc)) next

    # ---- admin-2 resolution ----
    svy_a2 <- tryCatch(tar_read_raw(paste0("svy_admin2_", sl, "_", ot)), error = function(e) NULL)
    gee_a2 <- tryCatch(tar_read_raw(paste0("gee_admin2_", sl)),           error = function(e) NULL)
    set.seed(SEED)
    a2 <- if (!is.null(svy_a2) && !is.null(gee_a2)) safe_fit(svy_a2, gee_a2) else c(n = NA, r = NA, mae = NA)

    # ---- cluster resolution ----
    cl_res <- c(n = NA, r = NA, mae = NA)
    cl <- tryCatch(aggregate_outcome_to_clusters(cc, oc), error = function(e) NULL)
    if (!is.null(cgee[[ck]]) && !is.null(cl) && nrow(cl) >= 8) {
      gc <- grep("^gee_", names(cgee[[ck]]), value = TRUE)
      m  <- merge(cl, cgee[[ck]][, c("cluster_id", gc)], by = "cluster_id")
      if (nrow(m) >= 8) {
        svy_cl <- data.frame(Admin2 = as.character(m$cluster_id),
                             svy_prev = m$svy_prev, n_svy = m$n_svy,
                             stringsAsFactors = FALSE)
        gee_cl <- m[, c("cluster_id", gc)]; names(gee_cl)[1] <- "Admin2"
        gee_cl$Admin2 <- as.character(gee_cl$Admin2)
        set.seed(SEED)
        cl_res <- safe_fit(svy_cl, gee_cl)
      }
    }

    rows[[paste0(ck, "_", ot)]] <- data.frame(
      country   = ck, outcome = ot,
      n_admin2  = unname(a2["n"]),     r_admin2 = unname(a2["r"]),     mae_admin2 = unname(a2["mae"]),
      n_cluster = unname(cl_res["n"]), r_cluster = unname(cl_res["r"]), mae_cluster = unname(cl_res["mae"]),
      stringsAsFactors = FALSE)
    cat(sprintf("[%s/%s] admin2 n=%s r=%s mae=%s | cluster n=%s r=%s mae=%s\n",
                ck, ot, a2["n"], a2["r"], a2["mae"], cl_res["n"], cl_res["r"], cl_res["mae"]))
    write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
  }
}
cat("\nwritten:", out_csv, "\n")
