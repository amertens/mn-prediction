# =============================================================================
# scripts/prototype_feature_engineering.R
#
# Measure whether engineered features (R/feature_engineering.R) improve
# cross-country transportability, vs the existing harmonized covariates.
# For each outcome: assemble the LOCO pooled dataset with (a) BASE covariates
# and (b) BASE + engineered (fe_*) covariates, run the selected recipe(s), and
# compare the pooled within-country-demeaned spatial r.
#
# Reads from _targets_full. Run AFTER the pipeline rerun completes (avoids
# store contention).
# =============================================================================
suppressMessages({ library(targets); library(dplyr) })
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
try(source(file.path("src", "0-functions.R")), silent = TRUE)
set.seed(20260521)

STORE <- "_targets_full"
COUNTRIES <- c("gambia", "ghana", "malawi", "sierraleone")
OUTCOMES  <- c("child_vitA", "child_iron", "women_vitA", "women_iron")
configs <- get_country_configs()
cc_for  <- function(cn) configs[[ names(configs)[tolower(names(configs)) == cn][1] ]]

cat("Loading inputs + building base & engineered covariates...\n")
merged_cache <- lapply(COUNTRIES, function(cn) tar_read_raw(paste0("merged_", cn), store = STORE))
names(merged_cache) <- COUNTRIES

cov_base <- list(); cov_fe <- list()
for (cn in COUNTRIES) {
  g  <- tar_read_raw(paste0("gee_admin2_", cn), store = STORE)
  cc <- cc_for(cn); sy <- cc$survey_year %||% NA_integer_
  base <- build_admin2_covariates(g, merged_cache[[cn]])
  fe   <- engineer_admin2_features(g, survey_year = sy)
  cov_base[[cn]] <- base
  cov_fe[[cn]]   <- dplyr::full_join(base, fe, by = "Admin2")
  cat(sprintf("  %-12s base=%d cov, +fe=%d cov (added %d engineered)\n",
              cn, ncol(base) - 1, ncol(cov_fe[[cn]]) - 1, ncol(fe) - 1))
}

uniform_svy <- function(cn, oc_tag) {
  cc <- cc_for(cn); oc <- cc$outcomes[[oc_tag]]; if (is.null(oc)) return(NULL)
  od <- tryCatch(build_outcome_dataset(merged_cache[[cn]], cc, oc), error = function(e) NULL)
  if (is.null(od) || nrow(od$data) == 0) return(NULL)
  tryCatch(compute_svy_admin2(od, cc, oc), error = function(e) NULL)
}

pooled_demeaned <- function(preds) {
  if (is.null(preds) || nrow(preds) == 0) return(c(pearson = NA, spearman = NA))
  z <- preds %>% group_by(outcome, country) %>%
    mutate(sc = survey_prev - mean(survey_prev), mc = modeled_prev - mean(modeled_prev)) %>% ungroup()
  c(pearson  = suppressWarnings(cor(z$sc, z$mc, use = "complete.obs")),
    spearman = suppressWarnings(cor(rank(z$sc), rank(z$mc), use = "complete.obs")))
}

recipes <- list(
  enet_corr30 = AREA_TRANSPORT_RECIPE,                 # selected parsimonious
  ridge_all   = AREA_TRANSPORT_RECIPE_RIDGE            # max-performance
)

svy_lists <- lapply(OUTCOMES, function(oc) {
  sl <- lapply(COUNTRIES, function(cn) uniform_svy(cn, oc)); names(sl) <- COUNTRIES; sl
})
names(svy_lists) <- OUTCOMES

rows <- list()
for (rc_name in names(recipes)) {
  for (covset in c("base", "base_plus_fe")) {
    cl <- if (covset == "base") cov_base else cov_fe
    preds <- bind_rows(lapply(OUTCOMES, function(oc) {
      p <- assemble_area_transport(svy_lists[[oc]], cl, oc)
      if (is.null(p)) return(NULL)
      run_area_transport_loco(p, recipes[[rc_name]])$predictions
    }))
    pd <- pooled_demeaned(preds)
    # per-outcome spearman
    z <- preds %>% group_by(outcome, country) %>%
      mutate(sc = survey_prev - mean(survey_prev), mc = modeled_prev - mean(modeled_prev)) %>%
      ungroup() %>% group_by(outcome) %>%
      summarise(sp = suppressWarnings(cor(rank(sc), rank(mc), use = "complete.obs")), .groups = "drop")
    rows[[paste(rc_name, covset)]] <- data.frame(
      recipe = rc_name, covset = covset,
      pooled_pearson = round(pd["pearson"], 3), pooled_spearman = round(pd["spearman"], 3),
      child_iron = round(z$sp[z$outcome == "child_iron"][1], 3),
      child_vitA = round(z$sp[z$outcome == "child_vitA"][1], 3),
      women_iron = round(z$sp[z$outcome == "women_iron"][1], 3),
      women_vitA = round(z$sp[z$outcome == "women_vitA"][1], 3),
      stringsAsFactors = FALSE)
  }
}
res <- bind_rows(rows)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/tables/feature_engineering_prototype.csv", row.names = FALSE)
options(width = 200)
cat("\n############ Feature-engineering prototype: base vs base+FE ############\n")
print(as.data.frame(res), row.names = FALSE)
cat("\nDONE\n")
