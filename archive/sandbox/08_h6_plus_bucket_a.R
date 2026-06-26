# =============================================================================
# sandbox/08_h6_plus_bucket_a.R
#
# Experimental: combine the H6 winner (within-country invariance + outcome-
# shared composite filter) with Bucket A features (spatial smoothing +
# within-domain PCA). The hypothesis is that smoothed copies and PCA
# components have HIGHER within-country variance ratios than the raw
# rasters — so the H6 filter applied to the augmented pool should select
# a more transportable variable set.
#
# Not promoted to the main pipeline: this stacks two experimental
# layers; one-off example only.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "feature_engineering_bucket_a.R"))

# Build adjacency for spatial smoothing.
source(here::here("R", "config.R"))
cc_list <- get_country_configs()

pooled_all <- load_all_pooled()

# Adjacency is per-country and only needs the svy_admin2_lists from one
# outcome's pooled object (any of them — Admin2 keys are the same).
example_svy_list <- pooled_all[[1]]$svy_list
adj <- build_adjacency_list(cc_list, example_svy_list)

augmented_pooled <- list()
for (oc in names(pooled_all)) {
  cat(sprintf("[H6+Bucket-A] augmenting %s ...\n", oc))
  p <- pooled_all[[oc]]
  aug <- apply_bucket_a_features(
    pooled_data         = p$pooled_data,
    gee_vars            = p$common_gee_vars,
    adjacency_list      = adj,
    svy_admin2_list     = example_svy_list,
    merged_country_list = NULL,           # DHS aggregates not safe for LOCO
    spatial_smooth = TRUE, domain_pca = TRUE,
    dhs_aggregates = FALSE, wash_gini = FALSE
  )
  augmented_pooled[[oc]] <- list(
    pooled_data = aug,
    common_gee_vars = attr(aug, "gee_vars"),
    country_names = p$country_names
  )
}

# Build cross-outcome data lookup (pooled_data per outcome) for the
# combined_filter method.
cross_pooled <- setNames(lapply(augmented_pooled, function(x) x$pooled_data),
                          names(augmented_pooled))

# Predict-fn: H6 combined_filter on the augmented variable set.
make_h6_aug_predict <- function(this_outcome) {
  function(train, test, vars) {
    out <- fit_predict_combined_filter(
      train, test, vars,
      cross_outcome_pooled = cross_pooled,
      this_outcome = this_outcome,
      min_within_ratio = 0.70, top_k = 12L)
    if (is.null(out)) NULL else out$pred
  }
}

res <- list()
for (oc in OUTCOMES) {
  cat(sprintf("\n[H6+Bucket-A] LOCO on %s ...\n", oc))
  r <- loco_evaluate(augmented_pooled[[oc]],
                       make_h6_aug_predict(oc),
                       outcome = oc, method_name = "h6_plus_bucketA")
  res[[oc]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "08_h6_plus_bucket_a.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "08_h6_plus_bucket_a.csv")))
