# =============================================================================
# scripts/covariates/06_mrp_comparison.R
#
# Run MRP for every country x outcome cell under BOTH covariate vocabularies and
# tabulate it against the district survey estimate, so the two changes -- the
# harmonised covariate set and the new MRP comparator -- can be judged before
# anyone commits to a full pipeline rebuild.
#
# Reads the built targets store; does not re-extract rasters. The harmonised
# area covariates are attached to the stored admin-2 keys directly, so this runs
# in minutes rather than hours.
#
# Metrics are IN-SAMPLE, matching how the other rows of area_comparison_all.csv
# are computed. They are a like-for-like comparison between methods, NOT an
# estimate of out-of-sample accuracy -- the honest out-of-sample numbers come
# from the block-CV corrected-methods analysis.
#
# Output: results/tables/mrp_comparison.csv
#
#   Rscript scripts/covariates/06_mrp_comparison.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
targets::tar_source("R/")

STORE <- here("_targets_full")
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)

cfgs <- get_country_configs()
rows <- list()

metrics <- function(svy, pred_df, pred_col, country, outcome, method, n_cov) {
  m <- dplyr::inner_join(svy[, c("Admin2", "svy_prev", "n_svy")],
                         pred_df[, c("Admin2", pred_col)], by = "Admin2")
  names(m)[names(m) == pred_col] <- "pred"
  m <- m[is.finite(m$svy_prev) & is.finite(m$pred), , drop = FALSE]
  if (nrow(m) < 5) return(NULL)
  data.frame(
    country = country, outcome = outcome, method = method,
    n_covariates = n_cov, n_admin2 = nrow(m),
    n_predicted = sum(is.finite(pred_df[[pred_col]])),
    pearson_r = round(suppressWarnings(stats::cor(m$svy_prev, m$pred)), 3),
    spearman_r = round(suppressWarnings(stats::cor(m$svy_prev, m$pred, method = "spearman")), 3),
    mae_pp = round(100 * mean(abs(m$svy_prev - m$pred)), 2),
    rmse_pp = round(100 * sqrt(mean((m$svy_prev - m$pred)^2)), 2),
    bias_pp = round(100 * mean(m$pred - m$svy_prev), 2),
    stringsAsFactors = FALSE)
}

for (ctry in names(cfgs)) {
  cc <- cfgs[[ctry]]
  leg <- rd(paste0("gee_admin2_", tolower(ctry)))
  if (is.null(leg)) { message("no stored covariates for ", ctry); next }
  keys <- leg[, intersect(c("Admin1", "Admin2"), names(leg)), drop = FALSE]
  har <- tryCatch(append_harmonized_admin2(keys, ctry, shared_only = FALSE),
                  error = function(e) NULL)

  for (oc_name in names(cc$outcomes)) {
    oc <- cc$outcomes[[oc_name]]
    suffix <- paste0(tolower(ctry), "_", oc_name)
    svy <- rd(paste0("svy_admin2_", suffix))
    od  <- rd(paste0("outcome_data_", suffix))
    if (is.null(svy) || is.null(od) || nrow(svy) < 5) next
    message("\n=== ", ctry, " / ", oc_name, " (", nrow(svy), " surveyed districts) ===")

    # national-mean null, for scale
    nat <- stats::weighted.mean(svy$svy_prev, pmax(svy$n_svy %||% 1, 1), na.rm = TRUE)
    rows[[length(rows) + 1L]] <- metrics(
      svy, data.frame(Admin2 = svy$Admin2, pred = nat), "pred",
      ctry, oc_name, "National mean (null)", 0L)

    for (vocab in c("legacy", "harmonized")) {
      acov <- if (vocab == "legacy") leg else har
      if (is.null(acov)) next
      Sys.setenv(COVARIATE_VOCAB = vocab)
      mrp <- tryCatch(fit_mrp_admin2(od, cc, oc, area_cov = acov, svy_admin2 = svy),
                      error = function(e) { message("  MRP failed: ", conditionMessage(e)); NULL })
      if (is.null(mrp)) next
      rows[[length(rows) + 1L]] <- metrics(
        svy, mrp, "mrp_prev", ctry, oc_name, paste0("MRP (", vocab, ")"),
        ncol(acov) - sum(c("Admin1", "Admin2") %in% names(acov)))
    }
  }
}
Sys.setenv(COVARIATE_VOCAB = "legacy")

out <- dplyr::bind_rows(rows)
dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out, here("results", "tables", "mrp_comparison.csv"))
message("\n-> results/tables/mrp_comparison.csv (", nrow(out), " rows)")

print(as.data.frame(out %>% group_by(method) %>%
  summarise(cells = dplyr::n(), median_r = round(stats::median(pearson_r, na.rm = TRUE), 3),
            median_mae_pp = round(stats::median(mae_pp, na.rm = TRUE), 2),
            .groups = "drop") %>% arrange(desc(median_r))), row.names = FALSE)
