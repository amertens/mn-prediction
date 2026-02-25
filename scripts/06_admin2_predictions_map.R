# =============================================================================
# scripts/06_admin2_predictions_map.R
#
# Admin2-level prediction, mapping, bootstrap uncertainty, and error
# evaluation for one outcome from the Gambia or Ghana MN pipeline.
#
# What this script does:
#   1. Loads the merged dataset and (if available) the saved SL model.
#      Falls back to fitting a minimal 2-learner SL inline so the script
#      runs even before the full production pipeline has been executed.
#   2. Aggregates individual-level SL predictions to Admin2 prevalence.
#   3. Downloads GADM level-2 boundaries and plots a choropleth map.
#   4. Runs a cluster bootstrap (low B for speed) to produce 95% CIs at
#      Admin2 level; maps CI width as an uncertainty choropleth.
#   5. Computes design-based (survey-weighted) Admin2 prevalence using srvyr.
#   6. Compares SL Admin2 predictions to survey-weighted prevalence:
#      MAE, RMSE, correlation, scatter plot, and per-Admin2 error table.
#
# Usage: set the three parameters below, then source / Rscript the file.
# All outputs are written to results/admin2/{tables,figures}/.
# =============================================================================

# ── USER-SETTABLE PARAMETERS ──────────────────────────────────────────────────
COUNTRY     <- "Gambia"       # "Gambia" or "Ghana"
OUTCOME_TAG <- "child_vitA"   # see OUTCOME REGISTRY below
B_BOOT      <- 10L            # bootstrap replicates (10 for dev; 200 for pub)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(geodata)
  library(srvyr)
  library(survey)
  library(viridis)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(sl3)
  library(origami)
  library(caret)
  library(data.table)
  library(ck37r)
  library(labelled)
  library(recipes)
  library(future.apply)
})

source(here::here("src", "analysis", "sl_helpers.R"))
source(here::here("src", "0-SL-setup.R"))  # defines slmod, slmod2_bin

options(survey.lonely.psu = "adjust")   # single-PSU strata: use stratum mean

# ── OUTPUT DIRECTORIES ────────────────────────────────────────────────────────
out_tables  <- here::here("results", "admin2", "tables")
out_figures <- here::here("results", "admin2", "figures")
dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_figures, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n[admin2] Country: %s | Outcome: %s | B = %d\n\n",
            COUNTRY, OUTCOME_TAG, B_BOOT))


# =============================================================================
# COUNTRY REGISTRY
# Maps country name → data paths, GADM code, survey design columns.
# =============================================================================
country_cfg <- list(

  Gambia = list(
    data_path   = here::here("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"),
    gadm_code   = "GMB",
    # Survey design (MICS 2018): cluster-only, no strata in merged dataset
    psu_col     = "gw_cnum",
    strata_col  = NULL,          # set NULL → survey design without strata
    weight_col  = "gw_svy_weight",
    admin2_col  = "Admin2",
    admin1_col  = "Admin1",
    child_flag  = "gw_child_flag",
    cluster_id  = "gw_cnum"
  ),

  Ghana = list(
    data_path   = here::here("data", "IPD", "Ghana", "Ghana_merged_dataset.rds"),
    gadm_code   = "GHA",
    psu_col     = "gw_EACode",
    strata_col  = "gw_Strata",
    weight_col  = "gw_sWeight",
    admin2_col  = "Admin2",
    admin1_col  = "Admin1",
    child_flag  = "gw_child_flag",
    cluster_id  = "gw_EACode"
  )
)

if (!COUNTRY %in% names(country_cfg))
  stop("Unknown COUNTRY. Choose 'Gambia' or 'Ghana'.")
cc <- country_cfg[[COUNTRY]]


# =============================================================================
# OUTCOME REGISTRY
# One entry per (country × outcome) with all column names and thresholds.
# =============================================================================
outcome_registry <- list(

  Gambia = list(

    child_vitA = list(
      label          = "Vitamin A deficiency (children)",
      population     = "children",
      child_flag_val = 1L,
      continuous     = "gw_cRBPAdjThurn",
      binary         = "gw_cVAD_Thurn",
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Gambia_SL_child_vitA.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Gambia_SL_child_vitA.rds")
    ),

    women_vitA = list(
      label          = "Vitamin A deficiency (women)",
      population     = "women",
      child_flag_val = 0L,
      continuous     = "gw_wRBPAdjThurn",
      binary         = "gw_wVAD_Thurn",
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Gambia_SL_women_vitA.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Gambia_SL_women_vitA.rds")
    ),

    child_iron = list(
      label          = "Iron deficiency anaemia (children)",
      population     = "children",
      child_flag_val = 1L,
      continuous     = "gw_LogFerAdj",
      binary         = "gw_cIDA_Brinda",
      cutoff         = log(12),
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Gambia_SL_child_iron.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Gambia_SL_child_iron.rds")
    ),

    women_iron = list(
      label          = "Iron deficiency anaemia (women)",
      population     = "women",
      child_flag_val = 0L,
      continuous     = "gw_LogFerAdj",
      binary         = "gw_wIDA_Brinda",
      cutoff         = log(15),
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Gambia_SL_mom_iron.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Gambia_SL_mom_iron.rds")
    )
  ),

  Ghana = list(

    child_vitA = list(
      label          = "Vitamin A deficiency (children)",
      population     = "children",
      child_flag_val = 1L,
      continuous     = "gw_cRBPAdjBrinda",
      binary         = "gw_cVAD",
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Ghana_SL_child_vitA.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Ghana_SL_child_vitA.rds")
    ),

    women_vitA = list(
      label          = "Vitamin A deficiency (women)",
      population     = "women",
      child_flag_val = 0L,
      continuous     = "gw_wRBPAdjThurn",
      binary         = "gw_wVADAdjThurn",
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Ghana_SL_women_vitA.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Ghana_SL_women_vitA.rds")
    ),

    child_iron = list(
      label          = "Iron deficiency anaemia (children)",
      population     = "children",
      child_flag_val = 1L,
      continuous     = "gw_cLogFerAdj",
      binary         = "gw_cIDAdjBrinda",
      cutoff         = log(12),
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Ghana_SL_child_iron.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Ghana_SL_child_iron.rds")
    ),

    women_iron = list(
      label          = "Iron deficiency anaemia (women)",
      population     = "women",
      child_flag_val = 0L,
      continuous     = "gw_wLogFerAdj",
      binary         = "gw_wIDAdjBrinda",
      cutoff         = log(15),
      cutoff_dir     = "less",
      model_bin      = here::here("results", "models",
                                  "res_bin_GW_Ghana_SL_women_iron.rds"),
      model_cont     = here::here("results", "models",
                                  "res_GW_Ghana_SL_women_iron.rds")
    )
  )
)

if (!OUTCOME_TAG %in% names(outcome_registry[[COUNTRY]]))
  stop(sprintf("Unknown OUTCOME_TAG '%s' for %s.", OUTCOME_TAG, COUNTRY))
oc <- outcome_registry[[COUNTRY]][[OUTCOME_TAG]]


# =============================================================================
# HELPERS
# =============================================================================

apply_threshold <- function(x, cutoff, direction = "less") {
  if (direction == "less") as.integer(x < cutoff) else as.integer(x > cutoff)
}

# Predictors: all domain-prefix columns MINUS outcome-leaking gw_ vars
gw_exclude_patterns <- c("RBP", "rbp", "VAD", "LogFer", "logfer",
                         "IDA", "Brinda", "Thurn")

build_Xvars <- function(d, gw_exclude_patterns) {
  domain_prefixes <- c("gw_", "dhs2019_", "dhs_", "mics_", "ihme_",
                       "lsms_", "MAP_", "wfp_", "flunet_", "gee_")
  all_vars <- colnames(d)
  domain_vars <- all_vars[Reduce(`|`, lapply(domain_prefixes,
                                             function(p) startsWith(all_vars, p)))]
  leakage_pat <- paste(gw_exclude_patterns, collapse = "|")
  gw_vars     <- domain_vars[startsWith(domain_vars, "gw_")]
  gw_clean    <- gw_vars[!grepl(leakage_pat, gw_vars)]
  non_gw      <- domain_vars[!startsWith(domain_vars, "gw_")]
  c(gw_clean, non_gw)
}


# =============================================================================
# §01  LOAD DATA
# =============================================================================
cat("[01] Loading data...\n")
d_raw <- readRDS(cc$data_path)

# Drop sf geometry if present
if (inherits(d_raw, "sf")) d_raw <- sf::st_drop_geometry(d_raw)

cat(sprintf("  Dataset: %d rows × %d columns\n", nrow(d_raw), ncol(d_raw)))

if (!cc$admin2_col %in% colnames(d_raw))
  stop(sprintf("Admin2 column '%s' not found. Check data merge.", cc$admin2_col))

n_a2 <- length(unique(d_raw[[cc$admin2_col]]))
cat(sprintf("  Admin2 regions: %d unique values\n", n_a2))

# Per-outcome dataset: right population, non-missing continuous outcome
pop_mask  <- d_raw[[cc$child_flag]] == oc$child_flag_val
cont_mask <- !is.na(d_raw[[oc$continuous]])
d_oc      <- d_raw[pop_mask & cont_mask, ]
cat(sprintf("  %s: n = %d (out of %d total)\n",
            oc$label, nrow(d_oc), nrow(d_raw)))

Xvars_full <- build_Xvars(d_raw, gw_exclude_patterns)
Xvars_oc   <- Xvars_full[Xvars_full %in% colnames(d_oc)]
cat(sprintf("  Predictors available: %d\n", length(Xvars_oc)))


# =============================================================================
# §02  LOAD OR FIT MODEL
# =============================================================================
cat("\n[02] Loading / fitting model...\n")

res_bin  <- NULL
res_cont <- NULL

if (file.exists(oc$model_bin)) {
  res_bin <- readRDS(oc$model_bin)
  cat(sprintf("  Loaded binary model: %s\n", basename(oc$model_bin)))
} else if (file.exists(oc$model_cont)) {
  res_cont <- readRDS(oc$model_cont)
  cat(sprintf("  Loaded continuous model (binary not found): %s\n",
              basename(oc$model_cont)))
} else {
  cat("  No saved model found — fitting minimal SL inline (mean + glmnet).\n")
  cat("  (Run scripts/run_gambia_pipeline.R first for full production models.)\n")

  # Minimal 2-learner stack for fast development runs
  lrnr_mean_  <- Lrnr_mean$new()
  lrnr_glmnet <- Lrnr_glmnet$new()
  stack_min   <- make_learner(Stack, lrnr_mean_, lrnr_glmnet)
  sl_min_bin  <- make_learner(Lrnr_sl,
                              learners      = stack_min,
                              loss_function = loss_loglik_binomial,
                              metalearner   = make_learner(Lrnr_nnls))

  # Binary outcome column must exist and have 2 levels
  d_fit <- d_oc[!is.na(d_oc[[oc$binary]]), ]
  if (length(unique(d_fit[[oc$binary]])) == 2) {
    cat("  Fitting binary model...\n")
    res_bin <- tryCatch(
      DHS_SL_clustered(
        d          = d_fit,
        Xvars      = Xvars_oc[Xvars_oc %in% colnames(d_fit)],
        outcome    = oc$binary,
        population = oc$population,
        id         = cc$cluster_id,
        folds      = 3L,         # 3 folds for speed
        CV         = FALSE,
        prescreen  = TRUE,
        sl         = sl_min_bin
      ),
      error = function(e) { message("  Binary fit failed: ", e$message); NULL }
    )
  }

  if (is.null(res_bin)) {
    sl_min_cont <- make_learner(Lrnr_sl,
                                learners      = stack_min,
                                loss_function = loss_squared_error,
                                metalearner   = make_learner(Lrnr_nnls))
    cat("  Falling back to continuous model...\n")
    res_cont <- tryCatch(
      DHS_SL_clustered(
        d          = d_oc,
        Xvars      = Xvars_oc,
        outcome    = oc$continuous,
        population = oc$population,
        id         = cc$cluster_id,
        folds      = 3L,
        CV         = FALSE,
        prescreen  = TRUE,
        sl         = sl_min_cont
      ),
      error = function(e) { message("  Continuous fit failed: ", e$message); NULL }
    )
  }
}

use_binary  <- !is.null(res_bin)
res_active  <- if (use_binary) res_bin else res_cont
sl_for_boot <- if (use_binary) slmod2_bin else slmod   # from 0-SL-setup.R
outcome_col <- if (use_binary) oc$binary else oc$continuous
model_type  <- if (use_binary) "binary_prob" else "continuous_threshold"

if (is.null(res_active)) stop("Model fitting failed — cannot continue.")

cat(sprintf("  Using model type: %s\n", model_type))


# =============================================================================
# §03  ADMIN2 SL PREDICTIONS
# =============================================================================
cat("\n[03] Aggregating SL predictions to Admin2...\n")

# res$res has: dataid, clusterid, outcome, population, Y, yhat_full
# The d_oc rows correspond 1-to-1 with res$res rows (same population + non-NA filter)
pred_df <- res_active$res %>%
  dplyr::select(dataid, Y, yhat_full)

# Attach Admin2 from d_oc by dataid
geo_df <- d_oc %>%
  dplyr::select(dataid,
                Admin2  = dplyr::all_of(cc$admin2_col),
                Admin1  = dplyr::all_of(cc$admin1_col)) %>%
  dplyr::mutate(across(c(Admin2, Admin1), as.character))

pred_geo <- dplyr::left_join(pred_df, geo_df, by = "dataid")

# Deficiency indicator / probability for each individual
pred_geo <- pred_geo %>%
  dplyr::mutate(
    deficient = if (use_binary) yhat_full
                else apply_threshold(yhat_full, oc$cutoff, oc$cutoff_dir)
  )

# Admin2-level: mean predicted prevalence, n, observed prevalence
sl_admin2 <- pred_geo %>%
  dplyr::filter(!is.na(Admin2)) %>%
  dplyr::group_by(Admin2) %>%
  dplyr::summarise(
    n_sl        = dplyr::n(),
    sl_prev     = mean(deficient, na.rm = TRUE),
    obs_prev_sl = mean(Y,         na.rm = TRUE),  # unweighted observed
    .groups     = "drop"
  )

cat(sprintf("  Admin2 units with SL predictions: %d\n", nrow(sl_admin2)))
print(sl_admin2 %>% arrange(desc(sl_prev)), n = 20)


# =============================================================================
# §04  DOWNLOAD GADM ADMIN2 BOUNDARIES
# =============================================================================
cat("\n[04] Downloading GADM level-2 boundaries...\n")

poly_a2 <- geodata::gadm(country = cc$gadm_code, level = 2, path = tempdir())
poly_a2 <- sf::st_as_sf(poly_a2) %>%
  dplyr::select(NAME_1, NAME_2) %>%
  dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2) %>%
  sf::st_transform(crs = 4326)

cat(sprintf("  %d Admin2 polygons downloaded\n", nrow(poly_a2)))

# Join SL predictions to polygons
map_sl <- poly_a2 %>%
  dplyr::left_join(sl_admin2, by = "Admin2")

n_unmatched <- sum(is.na(map_sl$sl_prev))
if (n_unmatched > 0)
  cat(sprintf("  Note: %d Admin2 polygons have no SL prediction (grey on map)\n",
              n_unmatched))


# =============================================================================
# §05  MAP 1 — Admin2 SL Predicted Prevalence
# =============================================================================
cat("\n[05] Plotting Admin2 SL predicted prevalence map...\n")

p_sl_map <- ggplot(map_sl) +
  geom_sf(aes(fill = sl_prev), colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name     = "Predicted\nprevalence",
    labels   = percent_format(accuracy = 1),
    limits   = c(0, 1),
    na.value = "grey85",
    option   = "plasma"
  ) +
  labs(
    title    = sprintf("SuperLearner Predicted Prevalence — %s", oc$label),
    subtitle = sprintf("%s | model: %s", COUNTRY, model_type),
    caption  = sprintf("n = %d individuals | %d Admin2 units",
                       nrow(pred_geo), nrow(sl_admin2))
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_figures,
    sprintf("%s_%s_admin2_sl_map.png", tolower(COUNTRY), OUTCOME_TAG)),
  plot = p_sl_map, width = 7, height = 6, dpi = 200
)
cat("  Saved: admin2_sl_map.png\n")


# =============================================================================
# §06  BOOTSTRAP UNCERTAINTY AT ADMIN2 LEVEL
# =============================================================================
cat(sprintf("\n[06] Bootstrap (B = %d) for Admin2 CIs...\n", B_BOOT))

# Admin2-level analogue of one_bootstrap from sl_helpers.R
one_bootstrap_admin2 <- function(b, d_boot_orig, Xvars_b, outcome_b,
                                 population_b, id_col, K, sl_obj,
                                 d_predict, geo_predict,
                                 cutoff, cutoff_dir, binary_outcome,
                                 admin2_col, admin1_col,
                                 seed_base = 12345L) {

  set.seed(seed_base + b)

  clusters      <- unique(d_boot_orig[[id_col]])
  boot_clusters <- sample(clusters, size = length(clusters), replace = TRUE)

  d_b <- do.call(rbind, lapply(boot_clusters, function(cl) {
    d_boot_orig[d_boot_orig[[id_col]] == cl, , drop = FALSE]
  }))

  fit_b <- tryCatch(
    DHS_SL_clustered(
      d          = d_b,
      Xvars      = Xvars_b[Xvars_b %in% colnames(d_b)],
      outcome    = outcome_b,
      population = population_b,
      id         = id_col,
      folds      = K,
      CV         = FALSE,
      prescreen  = TRUE,
      sl         = sl_obj
    ),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(NULL)

  # Predict on original data
  final_covars <- fit_b$task$nodes$covariates

  X_pred <- tryCatch({
    X0   <- d_predict %>%
      dplyr::select(dplyr::any_of(Xvars_b)) %>%
      as.data.frame()
    cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
    cov0 <- cov0 %>%
      do(ck37r::impute_missing_values(., type = "standard",
                                      add_indicators = TRUE,
                                      prefix = "missing_")$data) %>%
      as.data.frame()
    for (col in setdiff(final_covars, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, final_covars, drop = FALSE]
    data.table::data.table(cov0)
  }, error = function(e) NULL)

  if (is.null(X_pred)) return(NULL)

  pred_task <- tryCatch(
    sl3::sl3_Task$new(data = X_pred, covariates = final_covars, outcome = NULL),
    error = function(e) NULL
  )
  if (is.null(pred_task)) return(NULL)

  yhat <- tryCatch(
    as.numeric(fit_b$sl_fit$predict(pred_task)),
    error = function(e) NULL
  )
  if (is.null(yhat) || length(yhat) != nrow(d_predict)) return(NULL)

  deficient <- if (binary_outcome) yhat
               else as.numeric(apply_threshold(yhat, cutoff, cutoff_dir))

  # Aggregate to Admin2 and Admin1
  out <- data.frame(
    Admin2 = as.character(geo_predict[[admin2_col]]),
    Admin1 = as.character(geo_predict[[admin1_col]]),
    prev   = deficient
  )

  a2_prev <- out %>%
    dplyr::filter(!is.na(Admin2)) %>%
    dplyr::group_by(Admin2) %>%
    dplyr::summarise(prev = mean(prev, na.rm = TRUE), .groups = "drop")

  a1_prev <- out %>%
    dplyr::filter(!is.na(Admin1)) %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(prev = mean(prev, na.rm = TRUE), .groups = "drop")

  list(
    admin2   = a2_prev,
    admin1   = a1_prev,
    national = mean(deficient, na.rm = TRUE)
  )
}

# Prepare bootstrap inputs
d_fit_boot <- d_oc[!is.na(d_oc[[outcome_col]]), ]
Xvars_boot <- Xvars_oc[Xvars_oc %in% colnames(d_fit_boot)]

# geo_predict provides Admin2/Admin1 for the prediction dataset
geo_predict <- d_oc %>%
  dplyr::select(Admin2 = dplyr::all_of(cc$admin2_col),
                Admin1 = dplyr::all_of(cc$admin1_col)) %>%
  dplyr::mutate(across(everything(), as.character))

plan(sequential)   # use plan(multisession, workers = 4) to parallelise

boot_results <- future.apply::future_lapply(
  seq_len(B_BOOT),
  FUN            = one_bootstrap_admin2,
  d_boot_orig    = d_fit_boot,
  Xvars_b        = Xvars_boot,
  outcome_b      = outcome_col,
  population_b   = oc$population,
  id_col         = cc$cluster_id,
  K              = 3L,            # 3 folds for speed; use 5 for production
  sl_obj         = sl_for_boot,
  d_predict      = d_oc,
  geo_predict    = geo_predict,
  cutoff         = oc$cutoff,
  cutoff_dir     = oc$cutoff_dir,
  binary_outcome = use_binary,
  admin2_col     = cc$admin2_col,
  admin1_col     = cc$admin1_col,
  seed_base      = 12345L,
  future.seed    = TRUE,
  future.globals = TRUE
)

boot_results <- Filter(Negate(is.null), boot_results)
cat(sprintf("  Valid replicates: %d / %d\n", length(boot_results), B_BOOT))

if (length(boot_results) < 3) {
  warning("Too few valid bootstrap replicates. CIs will not be computed.")
} else {

  # Collect Admin2 prevalences across replicates
  all_a2 <- do.call(rbind, lapply(seq_along(boot_results), function(i) {
    df_i <- boot_results[[i]]$admin2
    df_i$rep <- i
    df_i
  }))

  admin2_ci <- all_a2 %>%
    dplyr::group_by(Admin2) %>%
    dplyr::summarise(
      n_boot    = dplyr::n(),
      boot_mean = mean(prev,            na.rm = TRUE),
      ci_lo     = quantile(prev, 0.025, na.rm = TRUE),
      ci_hi     = quantile(prev, 0.975, na.rm = TRUE),
      ci_width  = ci_hi - ci_lo,
      .groups   = "drop"
    )

  # National CI
  nat_vec <- sapply(boot_results, function(x) x$national)
  national_ci <- data.frame(
    outcome    = OUTCOME_TAG,
    model_type = model_type,
    n_reps     = length(nat_vec),
    boot_mean  = mean(nat_vec,            na.rm = TRUE),
    ci_lo      = quantile(nat_vec, 0.025, na.rm = TRUE),
    ci_hi      = quantile(nat_vec, 0.975, na.rm = TRUE)
  )

  cat(sprintf("  National: %.1f%% [%.1f%%, %.1f%%]\n",
              national_ci$boot_mean * 100,
              national_ci$ci_lo     * 100,
              national_ci$ci_hi     * 100))

  # Save bootstrap CIs
  readr::write_csv(admin2_ci,
    file.path(out_tables,
      sprintf("%s_%s_admin2_boot_ci.csv", tolower(COUNTRY), OUTCOME_TAG)))
  readr::write_csv(national_ci,
    file.path(out_tables,
      sprintf("%s_%s_national_boot_ci.csv", tolower(COUNTRY), OUTCOME_TAG)))


  # ── Map 2: CI width (uncertainty) ──────────────────────────────────────────
  cat("\n[06b] Plotting Admin2 uncertainty (CI width) map...\n")

  map_ci <- poly_a2 %>%
    dplyr::left_join(admin2_ci, by = "Admin2") %>%
    dplyr::left_join(sl_admin2  %>% dplyr::select(Admin2, sl_prev), by = "Admin2")

  p_ci_map <- ggplot(map_ci) +
    geom_sf(aes(fill = ci_width), colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(
      name     = "95% CI\nwidth (pp)",
      labels   = function(x) paste0(round(x * 100, 1)),
      na.value = "grey85",
      option   = "inferno",
      direction = -1       # dark = narrow CI = more precise
    ) +
    labs(
      title    = sprintf("Bootstrap CI Width — %s", oc$label),
      subtitle = sprintf("%s | B = %d replicates | darker = more precise",
                         COUNTRY, B_BOOT),
      caption  = "CI width = 97.5th − 2.5th percentile of bootstrap distribution"
    ) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(
    filename = file.path(out_figures,
      sprintf("%s_%s_admin2_ci_width_map.png", tolower(COUNTRY), OUTCOME_TAG)),
    plot = p_ci_map, width = 7, height = 6, dpi = 200
  )
  cat("  Saved: admin2_ci_width_map.png\n")


  # ── Forest plot: Admin2 CIs sorted by predicted prevalence ─────────────────
  cat("\n[06c] Plotting Admin2 forest plot...\n")

  forest_df <- admin2_ci %>%
    dplyr::left_join(sl_admin2 %>% dplyr::select(Admin2, n_sl), by = "Admin2") %>%
    dplyr::arrange(boot_mean) %>%
    dplyr::mutate(Admin2 = factor(Admin2, levels = Admin2))

  p_forest <- ggplot(forest_df,
                     aes(x = boot_mean, y = Admin2,
                         xmin = ci_lo, xmax = ci_hi)) +
    geom_errorbarh(colour = "steelblue", height = 0.4, linewidth = 0.6) +
    geom_point(colour = "steelblue", size = 2) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(
      x       = "Predicted prevalence (bootstrap mean + 95% CI)",
      y       = NULL,
      title   = sprintf("Admin2 Bootstrap CIs — %s", oc$label),
      subtitle = sprintf("%s | B = %d | sorted by predicted prevalence", COUNTRY, B_BOOT)
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 7))

  ggsave(
    filename = file.path(out_figures,
      sprintf("%s_%s_admin2_forest.png", tolower(COUNTRY), OUTCOME_TAG)),
    plot = p_forest,
    width = 7, height = max(5, nrow(forest_df) * 0.25 + 1.5), dpi = 200
  )
  cat("  Saved: admin2_forest.png\n")
}


# =============================================================================
# §07  SURVEY-WEIGHTED ADMIN2 PREVALENCE (design-based)
# =============================================================================
cat("\n[07] Computing survey-weighted Admin2 prevalence...\n")

# Use d_oc (population subset with non-missing continuous outcome)
# Restrict to rows with non-missing binary outcome and Admin2
svy_cols <- c(cc$psu_col, cc$weight_col, cc$admin2_col,
              cc$admin1_col, oc$binary)
if (!is.null(cc$strata_col)) svy_cols <- c(svy_cols, cc$strata_col)
svy_cols <- unique(svy_cols[svy_cols %in% colnames(d_oc)])

svy_df <- d_oc %>%
  dplyr::select(dplyr::all_of(svy_cols)) %>%
  dplyr::filter(!is.na(.data[[oc$binary]]),
                !is.na(.data[[cc$admin2_col]]),
                !is.na(.data[[cc$weight_col]])) %>%
  dplyr::mutate(
    across(dplyr::all_of(cc$psu_col), as.factor),
    !!oc$binary := as.numeric(.data[[oc$binary]]),
    !!cc$weight_col := as.numeric(.data[[cc$weight_col]])
  )

if (!is.null(cc$strata_col))
  svy_df <- svy_df %>%
    dplyr::mutate(across(dplyr::all_of(cc$strata_col), as.factor))

# Build survey design
if (!is.null(cc$strata_col)) {
  svy_des <- srvyr::as_survey_design(
    svy_df,
    ids     = !!rlang::sym(cc$psu_col),
    strata  = !!rlang::sym(cc$strata_col),
    weights = !!rlang::sym(cc$weight_col),
    nest    = TRUE
  )
} else {
  svy_des <- srvyr::as_survey_design(
    svy_df,
    ids     = !!rlang::sym(cc$psu_col),
    weights = !!rlang::sym(cc$weight_col)
  )
}

# Admin2-level design-based prevalence
svy_admin2 <- svy_des %>%
  dplyr::group_by(.data[[cc$admin2_col]]) %>%
  dplyr::summarise(
    svy_prev    = survey_mean(!!rlang::sym(oc$binary),
                              vartype = c("se", "ci"), na.rm = TRUE),
    n_svy       = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::rename(Admin2 = dplyr::all_of(cc$admin2_col)) %>%
  dplyr::mutate(
    svy_cv      = if_else(svy_prev > 0, svy_prev_se / svy_prev, NA_real_),
    svy_ci_width = svy_prev_upp - svy_prev_low
  )

cat(sprintf("  Admin2 units with survey estimate: %d\n", nrow(svy_admin2)))

# National design-based prevalence
svy_national <- svy_des %>%
  dplyr::summarise(
    svy_prev = survey_mean(!!rlang::sym(oc$binary), vartype = "ci", na.rm = TRUE)
  )
cat(sprintf("  Survey national prevalence: %.1f%% [%.1f%%, %.1f%%]\n",
            svy_national$svy_prev * 100,
            svy_national$svy_prev_low * 100,
            svy_national$svy_prev_upp * 100))


# =============================================================================
# §08  ERROR ANALYSIS: SL vs SURVEY-WEIGHTED ADMIN2 PREVALENCE
# =============================================================================
cat("\n[08] Computing Admin2 prediction error...\n")

# Merge SL predictions, bootstrap CIs, and survey estimates
merged_admin2 <- sl_admin2 %>%
  dplyr::left_join(svy_admin2 %>%
                     dplyr::select(Admin2, n_svy, svy_prev, svy_prev_se,
                                   svy_prev_low, svy_prev_upp, svy_cv,
                                   svy_ci_width),
                   by = "Admin2")

if (exists("admin2_ci")) {
  merged_admin2 <- merged_admin2 %>%
    dplyr::left_join(
      admin2_ci %>% dplyr::select(Admin2, boot_mean, ci_lo, ci_hi, ci_width),
      by = "Admin2"
    )
}

# Compute error metrics only where both estimates available
eval_df <- merged_admin2 %>%
  dplyr::filter(!is.na(sl_prev), !is.na(svy_prev))

n_eval <- nrow(eval_df)
cat(sprintf("  Admin2 units with both SL and survey estimate: %d\n", n_eval))

if (n_eval >= 2) {
  mae      <- mean(abs(eval_df$sl_prev - eval_df$svy_prev))
  rmse     <- sqrt(mean((eval_df$sl_prev - eval_df$svy_prev)^2))
  r_cor    <- cor(eval_df$sl_prev, eval_df$svy_prev, use = "complete.obs")
  mean_bias <- mean(eval_df$sl_prev - eval_df$svy_prev)

  cat(sprintf("  MAE:       %.3f (%.1f pp)\n", mae,       mae * 100))
  cat(sprintf("  RMSE:      %.3f (%.1f pp)\n", rmse,      rmse * 100))
  cat(sprintf("  Pearson r: %.3f\n",            r_cor))
  cat(sprintf("  Mean bias: %+.3f (%+.1f pp) [SL − survey]\n",
              mean_bias, mean_bias * 100))

  error_summary <- data.frame(
    outcome    = OUTCOME_TAG,
    country    = COUNTRY,
    model_type = model_type,
    n_admin2   = n_eval,
    mae_pp     = round(mae  * 100, 2),
    rmse_pp    = round(rmse * 100, 2),
    pearson_r  = round(r_cor,      3),
    mean_bias_pp = round(mean_bias * 100, 2)
  )
  readr::write_csv(error_summary,
    file.path(out_tables,
      sprintf("%s_%s_admin2_error_summary.csv", tolower(COUNTRY), OUTCOME_TAG)))
} else {
  cat("  Not enough overlapping Admin2 units for error metrics.\n")
}

# Full per-Admin2 table
merged_admin2 <- merged_admin2 %>%
  dplyr::mutate(
    abs_error = abs(sl_prev - svy_prev),
    bias      = sl_prev - svy_prev
  ) %>%
  dplyr::arrange(desc(abs_error))

readr::write_csv(merged_admin2,
  file.path(out_tables,
    sprintf("%s_%s_admin2_full_table.csv", tolower(COUNTRY), OUTCOME_TAG)))
cat("  Saved: admin2_full_table.csv\n")


# =============================================================================
# §09  SCATTER PLOT: Survey-weighted vs SL-predicted Admin2 prevalence
# =============================================================================
cat("\n[09] Plotting scatter: survey vs SL predicted prevalence...\n")

if (n_eval >= 2) {
  lim_lo <- min(c(eval_df$sl_prev, eval_df$svy_prev), na.rm = TRUE) - 0.03
  lim_hi <- max(c(eval_df$sl_prev, eval_df$svy_prev), na.rm = TRUE) + 0.03
  lim_lo <- max(0, lim_lo); lim_hi <- min(1, lim_hi)

  p_scatter <- ggplot(eval_df,
                      aes(x = svy_prev, y = sl_prev)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = "grey50", linewidth = 0.8) +
    # Survey confidence intervals as horizontal error bars
    {if ("svy_prev_low" %in% colnames(eval_df))
      geom_errorbarh(aes(xmin = svy_prev_low, xmax = svy_prev_upp),
                     colour = "grey70", height = 0, linewidth = 0.5)
    } +
    # Bootstrap CIs as vertical error bars (if available)
    {if (exists("admin2_ci") && "ci_lo" %in% colnames(eval_df))
      geom_linerange(aes(ymin = ci_lo, ymax = ci_hi),
                     colour = "steelblue", alpha = 0.5, linewidth = 0.5)
    } +
    geom_point(aes(size = n_svy), colour = "steelblue", alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = Admin2), size = 2.5,
                             max.overlaps = 15, colour = "grey30") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(lim_lo, lim_hi)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(lim_lo, lim_hi)) +
    scale_size_area(name = "n (survey)", max_size = 6) +
    annotate("text", x = lim_lo + 0.02, y = lim_hi - 0.02,
             label = sprintf("MAE = %.1f pp\nRMSE = %.1f pp\nr = %.2f",
                             mae * 100, rmse * 100, r_cor),
             hjust = 0, vjust = 1, size = 3.5,
             colour = "grey20",
             fontface = "italic") +
    labs(
      x       = "Survey-weighted Admin2 prevalence (design-based)",
      y       = "SL predicted Admin2 prevalence",
      title   = sprintf("Admin2 Prediction Error — %s", oc$label),
      subtitle = sprintf("%s | Horizontal bars = survey 95%% CI | Vertical = bootstrap CI",
                         COUNTRY),
      caption  = "Dashed line = perfect agreement (45°). Point size ∝ survey n."
    ) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(out_figures,
      sprintf("%s_%s_admin2_scatter.png", tolower(COUNTRY), OUTCOME_TAG)),
    plot = p_scatter, width = 7, height = 7, dpi = 200
  )
  cat("  Saved: admin2_scatter.png\n")
}


# =============================================================================
# §10  MAP 3 — Admin2 absolute error (|SL − survey|)
# =============================================================================
cat("\n[10] Plotting Admin2 error magnitude map...\n")

map_err <- poly_a2 %>%
  dplyr::left_join(
    merged_admin2 %>% dplyr::select(Admin2, sl_prev, svy_prev, abs_error, bias),
    by = "Admin2"
  )

p_err_map <- ggplot(map_err) +
  geom_sf(aes(fill = abs_error), colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name     = "|SL − survey|\n(pp)",
    labels   = function(x) paste0(round(x * 100, 1)),
    na.value = "grey85",
    option   = "magma",
    direction = -1
  ) +
  labs(
    title    = sprintf("Admin2 Absolute Prediction Error — %s", oc$label),
    subtitle = sprintf("%s | darker = smaller error", COUNTRY),
    caption  = "Grey = Admin2 units missing survey or SL estimate"
  ) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(
  filename = file.path(out_figures,
    sprintf("%s_%s_admin2_error_map.png", tolower(COUNTRY), OUTCOME_TAG)),
  plot = p_err_map, width = 7, height = 6, dpi = 200
)
cat("  Saved: admin2_error_map.png\n")


# =============================================================================
# SUMMARY
# =============================================================================
cat("\n", strrep("=", 65), "\n")
cat(sprintf("  Admin2 analysis complete: %s — %s\n", COUNTRY, oc$label))
cat(strrep("=", 65), "\n\n")

cat("  Tables:\n")
for (f in list.files(out_tables, full.names = TRUE))
  cat(sprintf("    %s\n", f))

cat("\n  Figures:\n")
for (f in list.files(out_figures, full.names = TRUE))
  cat(sprintf("    %s\n", f))

if (n_eval >= 2) {
  cat(sprintf(
    "\n  Error summary:\n    MAE = %.1f pp | RMSE = %.1f pp | r = %.2f | bias = %+.1f pp\n",
    mae * 100, rmse * 100, r_cor, mean_bias * 100
  ))
}
