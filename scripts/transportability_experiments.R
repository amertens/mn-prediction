# =============================================================================
# scripts/transportability_experiments.R
#
# Model / prescreen bake-off for the AREA-LEVEL cross-country transportability
# model, on the HARMONIZED outcome (uniform WHO cutoffs via compute_svy_admin2)
# and the FULL covariate set (_targets_full, ~1.4k GEE+IHME+MAP+fsec vars).
#
# Re-selects the best universal/parsimonious recipe after outcome harmonization,
# focusing on screen size (correlation screening overfits with p >> n) and
# regularization. Scored by the pooled within-country-demeaned spatial r
# (level removed per country) — the robust single measure of spatial transfer.
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

cat("Loading merged + building covariates + harmonized survey prevalence...\n")
merged_cache <- lapply(COUNTRIES, function(cn) tar_read_raw(paste0("merged_", cn), store = STORE))
names(merged_cache) <- COUNTRIES
cov_list <- lapply(COUNTRIES, function(cn)
  build_admin2_covariates(tar_read_raw(paste0("gee_admin2_", cn), store = STORE), merged_cache[[cn]]))
names(cov_list) <- COUNTRIES

uniform_svy <- function(cn, oc_tag) {
  cc <- cc_for(cn); oc <- cc$outcomes[[oc_tag]]; if (is.null(oc)) return(NULL)
  od <- tryCatch(build_outcome_dataset(merged_cache[[cn]], cc, oc), error = function(e) NULL)
  if (is.null(od) || nrow(od$data) == 0) return(NULL)
  tryCatch(compute_svy_admin2(od, cc, oc), error = function(e) NULL)
}
pooled_list <- lapply(OUTCOMES, function(oc) {
  sl <- lapply(COUNTRIES, function(cn) uniform_svy(cn, oc)); names(sl) <- COUNTRIES
  assemble_area_transport(sl, cov_list, oc)
})
names(pooled_list) <- OUTCOMES
for (oc in OUTCOMES) if (!is.null(pooled_list[[oc]]))
  cat(sprintf("  %-12s %d areas, %d common predictors\n", oc,
              nrow(pooled_list[[oc]]$data), length(pooled_list[[oc]]$predictors)))

# ---- recipe grid -----------------------------------------------------------
R <- function(model, K = NULL, center = FALSE)
  list(model = model, alpha = 0.5, screen_K = K, center = center,
       weight = TRUE, scale = "logit", lambda = "lambda.min", seed = 12345)
recipes <- list(
  # ridge (no hard selection — good for many correlated predictors)
  ridge_all_C     = R("ridge"),                 ridge_all = R("ridge", NULL, FALSE),
  ridge_c120_C    = R("ridge", 120, TRUE),      ridge_c60_C = R("ridge", 60, TRUE),
  ridge_c30_C     = R("ridge", 30, TRUE),
  # elastic net
  enet_all_C      = R("enet", NULL, TRUE),      enet_c120_C = R("enet", 120, TRUE),
  enet_c60_C      = R("enet", 60, TRUE),        enet_c30_C  = R("enet", 30, TRUE),
  enet_c15_C      = R("enet", 15, TRUE),        enet_c30    = R("enet", 30, FALSE),
  # lasso
  lasso_all_C     = R("lasso", NULL, TRUE),     lasso_c60_C = R("lasso", 60, TRUE),
  # random forest
  rf_c120_C       = R("rf", 120, TRUE),         rf_c60_C    = R("rf", 60, TRUE),
  rf_c30_C        = R("rf", 30, TRUE)
)

# pooled within-country-demeaned spatial correlation across all held-out areas
pooled_demeaned <- function(preds) {
  if (is.null(preds) || nrow(preds) == 0) return(c(pearson = NA, spearman = NA))
  z <- preds %>% group_by(outcome, country) %>%
    mutate(sc = survey_prev - mean(survey_prev),
           mc = modeled_prev - mean(modeled_prev)) %>% ungroup()
  c(pearson  = suppressWarnings(cor(z$sc, z$mc, use = "complete.obs")),
    spearman = suppressWarnings(cor(rank(z$sc), rank(z$mc), use = "complete.obs")))
}

cat("\nRunning recipe grid...\n")
rows <- list()
for (nm in names(recipes)) {
  rc <- recipes[[nm]]
  preds <- bind_rows(lapply(pooled_list, function(p)
    if (!is.null(p)) run_area_transport_loco(p, rc)$predictions))
  mets  <- bind_rows(lapply(pooled_list, function(p)
    if (!is.null(p)) run_area_transport_loco(p, rc)$metrics))
  pd <- pooled_demeaned(preds)
  # per-outcome mean spearman (each held-out country equally weighted)
  oc_sp <- if (!is.null(mets) && nrow(mets)) mets %>% group_by(outcome) %>%
    summarise(s = mean(spearman_r, na.rm = TRUE), .groups = "drop") else NULL
  rows[[nm]] <- data.frame(
    recipe = nm,
    pooled_pearson  = round(pd["pearson"], 3),
    pooled_spearman = round(pd["spearman"], 3),
    mean_country_spearman = if (!is.null(mets) && nrow(mets)) round(mean(mets$spearman_r, na.rm = TRUE), 3) else NA,
    iron_child  = if (!is.null(oc_sp)) round(oc_sp$s[oc_sp$outcome == "child_iron"][1], 3) else NA,
    iron_women  = if (!is.null(oc_sp)) round(oc_sp$s[oc_sp$outcome == "women_iron"][1], 3) else NA,
    vitA_child  = if (!is.null(oc_sp)) round(oc_sp$s[oc_sp$outcome == "child_vitA"][1], 3) else NA,
    med_selected = if (!is.null(mets) && nrow(mets)) round(median(mets$n_selected), 0) else NA,
    n_null = if (!is.null(mets) && nrow(mets)) sum(mets$n_selected == 0, na.rm = TRUE) else NA,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-14s pooled r(spear)=%.3f  meanC=%.3f\n", nm,
              pd["spearman"], rows[[nm]]$mean_country_spearman))
}
res <- bind_rows(rows) %>% arrange(desc(pooled_spearman))
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/tables/transportability_recipe_bakeoff.csv", row.names = FALSE)
options(width = 200)
cat("\n############ RECIPE BAKE-OFF (harmonized outcome, full covariates) ############\n")
print(as.data.frame(res), row.names = FALSE)
cat("\nDONE\n")
