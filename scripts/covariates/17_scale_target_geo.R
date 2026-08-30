# =============================================================================
# scripts/covariates/17_scale_target_geo.R
#
# CONFIRM THE THREE SANDBOX FINDINGS ON THE CURRENT DATA, ONE AT A TIME.
#
# archive/sandbox_parsimony (docs/parsimony_findings.md) reported three changes
# worth making. They were measured on an earlier covariate set, so before the
# defaults move they are re-measured here on the 294-predictor harmonised set,
# separately, so each can be attributed:
#
#   1 response scale   logit (clamped at 0.005) vs raw prevalence   (section 7)
#   2 LOCO target      level vs within-country z-score              (section 6)
#   3 geography        with vs without geo_lon / geo_lat            (section 5)
#
# All three are changed one at a time against the same baseline, and the
# combination is reported last. Changing three defaults at once and reporting
# the total would make it impossible to say which one paid.
#
#   Rscript scripts/covariates/17_scale_target_geo.R
# -> results/tables/scale_target_geo.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")
OUTCOMES <- c("child_iron", "child_vitA", "women_iron", "women_vitA")

# Attach centroids on the fly so the comparison does not require rebuilding
# area_covariates_*; extract_area_covariates() now does this in the pipeline.
add_geo <- function(g, gadm) {
  cen <- tryCatch(sf::st_drop_geometry(load_admin2_centroids(gadm)),
                  error = function(e) NULL)
  if (is.null(cen)) return(g)
  by <- admin2_join_by(g, cen)
  out <- dplyr::left_join(g, cen[, c(by, "lon", "lat")], by = by)
  out$geo_lon <- out$lon; out$geo_lat <- out$lat
  out$lon <- NULL; out$lat <- NULL
  out
}

cfgs <- get_country_configs()
svy_list <- list(); gee_list <- list(); gee_geo <- list()
for (cn in names(cfgs)) {
  lc <- tolower(cn)
  g <- tryCatch(targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE),
                error = function(e) NULL)
  if (is.list(g) && !is.data.frame(g) && "gee_admin2" %in% names(g)) g <- g$gee_admin2
  if (is.null(g)) next
  gee_list[[cn]] <- g
  gee_geo[[cn]]  <- add_geo(g, cfgs[[cn]]$gadm_code)
  for (on in OUTCOMES) {
    s <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on), store = STORE),
                  error = function(e) NULL)
    if (!is.null(s) && nrow(s)) svy_list[[paste(cn, on)]] <- s
  }
}

variants <- list(
  baseline        = list(scale = "logit", target = "level", geo = FALSE),
  raw_scale       = list(scale = "raw",   target = "level", geo = FALSE),
  zscore_target   = list(scale = "logit", target = "zscore", geo = FALSE),
  geo_coords      = list(scale = "logit", target = "level", geo = TRUE),
  all_three       = list(scale = "raw",   target = "zscore", geo = TRUE))

rows <- list()
for (on in OUTCOMES) {
  sl <- svy_list[grepl(paste0(" ", on, "$"), names(svy_list))]
  if (length(sl) < 2) next
  names(sl) <- sub(paste0(" ", on, "$"), "", names(sl))
  for (vn in names(variants)) {
    v <- variants[[vn]]
    gl <- if (v$geo) gee_geo[names(sl)] else gee_list[names(sl)]
    # assemble_area_transport(), NOT build_area_loco_dataset(): the two produce
    # different shapes (data/predictors/countries vs
    # pooled_data/common_gee_vars/country_names) and only the former is what
    # run_area_transport_loco() reads. Using the wrong one returns metrics =
    # NULL for every cell without erroring.
    pooled <- tryCatch(assemble_area_transport(sl, gl, outcome = on),
                       error = function(e) NULL)
    if (is.null(pooled)) next
    rec <- utils::modifyList(AREA_TRANSPORT_RECIPE,
                             list(scale = v$scale, target = v$target))
    res <- tryCatch(run_area_transport_loco(pooled, rec), error = function(e) {
      message("  ", on, "/", vn, ": ", conditionMessage(e)); NULL })
    if (is.null(res) || is.null(res$metrics) || !nrow(res$metrics)) {
      message("  ", on, "/", vn, ": no metrics returned"); next }
    m <- res$metrics
    rows[[length(rows) + 1L]] <- data.frame(
      outcome = on, variant = vn, scale = v$scale, target = v$target,
      geo = v$geo, cells = nrow(m),
      spearman = round(mean(m$spearman_r, na.rm = TRUE), 3),
      pearson  = round(mean(m$pearson_r, na.rm = TRUE), 3),
      mae_pp   = round(mean(m$mae_pp, na.rm = TRUE), 2),
      abs_bias = round(mean(abs(m$nat_bias_pp), na.rm = TRUE), 2),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-11s %-14s rho=%+.3f  r=%+.3f  MAE=%5.2f  |bias|=%5.2f\n",
                on, vn, mean(m$spearman_r, na.rm = TRUE),
                mean(m$pearson_r, na.rm = TRUE), mean(m$mae_pp, na.rm = TRUE),
                mean(abs(m$nat_bias_pp), na.rm = TRUE)))
  }
}

res <- dplyr::bind_rows(rows)
readr::write_csv(res, here("results", "tables", "scale_target_geo.csv"))

cat("\n================ ONE CHANGE AT A TIME, POOLED OVER OUTCOMES ==========\n")
s <- res %>% group_by(variant, scale, target, geo) %>%
  summarise(outcomes = n(), spearman = round(mean(spearman), 3),
            pearson = round(mean(pearson), 3), mae_pp = round(mean(mae_pp), 2),
            abs_bias = round(mean(abs_bias), 2), .groups = "drop") %>%
  arrange(desc(spearman))
print(as.data.frame(s), row.names = FALSE)

b <- res %>% filter(variant == "baseline") %>% select(outcome, base = spearman)
cat("\n--- change vs baseline, per outcome (Spearman) ---\n")
d <- res %>% filter(variant != "baseline") %>% inner_join(b, by = "outcome") %>%
  mutate(delta = round(spearman - base, 3))
print(d %>% select(outcome, variant, base, spearman, delta) %>%
        tidyr::pivot_wider(id_cols = outcome, names_from = variant,
                           values_from = delta) %>% as.data.frame(), row.names = FALSE)
cat("\n--- mean delta and win rate ---\n")
print(d %>% group_by(variant) %>%
        summarise(mean_delta = round(mean(delta), 3),
                  better_in = sprintf("%d/%d", sum(delta > 0), n()),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)
cat("\n-> results/tables/scale_target_geo.csv\n")
