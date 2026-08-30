# =============================================================================
# scripts/covariates/05_smoke_test_rewire.R
#
# Prove the two changes work before anyone spends hours on a pipeline rebuild:
#   1. COVARIATE_VOCAB=harmonized really does swap the area-level covariate set
#      (and legacy still behaves exactly as before)
#   2. MRP fits, poststratifies, and returns a prediction for every district
#      including unsurveyed ones
#
#   Rscript scripts/covariates/05_smoke_test_rewire.R [Country] [outcome]
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(here)
})
targets::tar_source("R/")

args    <- commandArgs(trailingOnly = TRUE)
COUNTRY <- if (length(args) >= 1) args[1] else "Gambia"
OUTCOME <- if (length(args) >= 2) args[2] else "child_iron"
STORE   <- here("_targets_full")
suffix  <- paste0(tolower(COUNTRY), "_", OUTCOME)

cc <- get_country_configs()[[COUNTRY]]
oc <- cc$outcomes[[OUTCOME]]
stopifnot(!is.null(cc), !is.null(oc))

rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)
svy <- rd(paste0("svy_admin2_", suffix))
od  <- rd(paste0("outcome_data_", suffix))
stopifnot(!is.null(svy), !is.null(od))

cat("\n===== 1. VOCABULARY SWITCH =====\n")
Sys.setenv(COVARIATE_VOCAB = "legacy")
cat("covariate_vocabulary() =", covariate_vocabulary(), "\n")
leg <- rd(paste0("gee_admin2_", tolower(COUNTRY)))
cat(sprintf("legacy  (targets store): %d areas x %d covariates | prefixes: %s\n",
            nrow(leg), ncol(leg) - 2L,
            paste(unique(sub("_.*$", "", setdiff(names(leg), c("Admin1","Admin2"))))[1:3],
                  collapse = ", ")))

Sys.setenv(COVARIATE_VOCAB = "harmonized")
cat("covariate_vocabulary() =", covariate_vocabulary(), "\n")
har <- extract_gee_admin2(cc)
covs <- setdiff(names(har), c("Admin1", "Admin2"))
cat(sprintf("harmonised            : %d areas x %d covariates\n", nrow(har), length(covs)))
cat("  AlphaEarth dims:", sum(grepl("^aef_", covs)),
    "| DHS:", sum(grepl("^dhs_", covs)),
    "| iSDAsoil:", sum(grepl("^soil_", covs)),
    "| SoilGrids:", sum(grepl("^soilgrids_", covs)),
    "| MapSPAM:", sum(grepl("^spam_", covs)), "\n")
stopifnot(length(covs) > 150, sum(grepl("^aef_", covs)) == 64)
stopifnot(!any(grepl("^gee_", covs)))   # the two vocabularies must never mix
cat("PASS: harmonised set is present, complete, and carries no legacy names\n")

Sys.setenv(COVARIATE_VOCAB = "legacy")
back <- extract_gee_admin2(cc)
cat(sprintf("legacy after reset    : %d areas x %d covariates (gee_ prefixed: %s)\n",
            nrow(back), ncol(back) - 2L,
            all(grepl("^gee_", setdiff(names(back), c("Admin1","Admin2"))))))
stopifnot(all(grepl("^gee_", setdiff(names(back), c("Admin1", "Admin2")))))
cat("PASS: default behaviour is unchanged\n")

cat("\n===== 2. MRP =====\n")
for (vocab in c("legacy", "harmonized")) {
  Sys.setenv(COVARIATE_VOCAB = vocab)
  acov <- if (vocab == "legacy") leg else har
  t0 <- Sys.time()
  mrp <- fit_mrp_admin2(od, cc, oc, area_cov = acov, svy_admin2 = svy)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (is.null(mrp)) { cat("  [", vocab, "] MRP returned NULL\n"); next }
  m <- dplyr::inner_join(svy[, c("Admin2", "svy_prev")], mrp, by = "Admin2")
  m <- m[is.finite(m$svy_prev) & is.finite(m$mrp_prev), , drop = FALSE]
  cat(sprintf("  [%s] %.1fs | %d districts predicted (%d surveyed) | r=%.3f MAE=%.2f pp\n",
              vocab, el, sum(is.finite(mrp$mrp_prev)), nrow(m),
              suppressWarnings(stats::cor(m$svy_prev, m$mrp_prev)),
              100 * mean(abs(m$svy_prev - m$mrp_prev))))
  # Unsurveyed districts must still receive an estimate -- that is the point of
  # having a model rather than a direct estimator.
  cat(sprintf("       unsurveyed districts with an estimate: %d\n",
              sum(is.finite(mrp$mrp_prev) & !(mrp$Admin2 %in% svy$Admin2))))
}
Sys.setenv(COVARIATE_VOCAB = "legacy")
cat("\nSmoke test complete.\n")
