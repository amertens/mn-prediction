# =============================================================================
# scripts/covariates/10_confirm_ablation.R
#
# The ablation in 09 re-implements the area-level ensemble for speed. Before any
# of its conclusions change the pipeline, re-check the shortlisted learner sets
# against the REAL estimator -- fit_area_sl(), through mlr3superlearner, with
# its own inner/outer CV -- on a spread of cells.
#
# A configuration that wins in the fast harness but not here does not get
# adopted: the harness is evidence, the production path is the arbiter.
#
#   Rscript scripts/covariates/10_confirm_ablation.R
# Output: results/tables/estimator_ablation_confirm.csv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
targets::tar_source("R/")
STORE <- here("_targets_full")
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
cfgs <- get_country_configs()

LIBS <- list(
  "current 5-learner stack" = list(
    list("glmnet", alpha = 1, id = "lasso"), list("glmnet", alpha = 0, id = "ridge"),
    list("glmnet", alpha = 0.5, id = "elastic_net"),
    list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
    list("xgboost", max_depth = 2, eta = 0.05, nrounds = 200,
         min_child_weight = 3, subsample = 0.8, id = "xgb")),
  "linear only" = list(
    list("glmnet", alpha = 1, id = "lasso"), list("glmnet", alpha = 0, id = "ridge"),
    list("glmnet", alpha = 0.5, id = "elastic_net")),
  "enet only" = list(list("glmnet", alpha = 0.5, id = "elastic_net"))
)

# A spread of sample sizes and reliabilities, not a convenience sample.
CELLS <- list(c("Gambia", "child_iron"), c("Gambia", "women_iron"),
              c("Ghana", "child_iron"), c("Ghana", "women_folate"),
              c("Malawi", "child_zinc"), c("Malawi", "women_folate"),
              c("SierraLeone", "women_folate"))

rows <- list()
for (cs in CELLS) {
  ctry <- cs[1]; ocn <- cs[2]
  cc <- cfgs[[ctry]]; if (is.null(cc) || !ocn %in% names(cc$outcomes)) next
  svy <- rd(paste0("svy_admin2_", tolower(ctry), "_", ocn)); if (is.null(svy)) next
  hc  <- H[H$country == ctry, ]
  rel <- admin2_reliability(svy, deff = 1.5, boot = 0)
  for (lib_name in names(LIBS)) {
    message("\n### ", ctry, " / ", ocn, " -- ", lib_name)
    t0 <- Sys.time()
    fit <- tryCatch(fit_area_sl(svy, hc[, setdiff(names(hc), "country")],
                                outcome_type = "binomial",
                                sl_library = LIBS[[lib_name]]),
                    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(fit) || is.null(fit$metrics)) next
    m <- fit$metrics
    rows[[length(rows) + 1L]] <- data.frame(
      country = ctry, outcome = ocn, library = lib_name,
      n_learners = length(LIBS[[lib_name]]), seconds = round(el, 1),
      pearson_r = round(m$pearson_r %||% NA_real_, 3),
      r_max = round(rel$r_max %||% NA_real_, 3),
      mae_pp = round(m$mae_pp %||% NA_real_, 2),
      stringsAsFactors = FALSE)
  }
}

out <- dplyr::bind_rows(rows)
if (nrow(out)) {
  out$r_share <- round(out$pearson_r / out$r_max, 2)
  readr::write_csv(out, here("results", "tables", "estimator_ablation_confirm.csv"))
  message("\n-> results/tables/estimator_ablation_confirm.csv")
  print(as.data.frame(out), row.names = FALSE)
  print(as.data.frame(out %>% group_by(library) %>%
    summarise(cells = dplyr::n(), median_r = round(stats::median(pearson_r, na.rm = TRUE), 3),
              median_mae_pp = round(stats::median(mae_pp, na.rm = TRUE), 2),
              median_seconds = round(stats::median(seconds), 1), .groups = "drop")),
    row.names = FALSE)
} else message("no results")
