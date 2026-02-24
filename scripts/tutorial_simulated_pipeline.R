# =============================================================================
# scripts/tutorial_simulated_pipeline.R
#
# A self-contained walkthrough of the full micronutrient analysis pipeline
# using SIMULATED data.  No external files are needed beyond the helper
# scripts that live in src/analysis/.
#
# Each section mirrors a numbered pipeline script:
#   §01  Simulate data + construct per-outcome datasets
#   §02  Fit SuperLearner models
#   §03  Predict and aggregate to Admin1 prevalence
#   §04  Bootstrap 95% confidence intervals
#   §05  Domain ablation study
#
# Uses a 2-learner stack (mean + GLM) so the whole script runs in ~1-2 min.
# Swap `stack_tut` for `stack` or `stack2` from 0-SL-setup.R for a full run.
# =============================================================================

library(here)
library(dplyr)
library(data.table)
library(sl3)
library(origami)
library(caret)
library(labelled)
library(recipes)
library(ck37r)
library(washb)
library(pROC)
library(readr)
library(ggplot2)

# Load the two helper functions: DHS_SL_clustered and one_bootstrap
source(here::here("src/analysis/sl_helpers.R"))

cat("\n")
cat("=============================================================\n")
cat("  Tutorial: Simulated Micronutrient Prediction Pipeline\n")
cat("=============================================================\n\n")


# =============================================================================
# §00  CONFIGURATION
# =============================================================================
# The real pipeline reads this from src/analysis/config.R.
# Here we define a minimal version with only 2 outcomes (child + women VitA)
# and 3 predictor domains (GW, DHS, MICS) for simplicity.

# --- Fast 2-learner SL stack for the tutorial --------------------------------
# The real pipeline uses stack2 (5 learners inc. screeners) from 0-SL-setup.R.
# Here we use mean + GLM so the tutorial completes in minutes, not hours.
lrnr_mean <- Lrnr_mean$new()
lrnr_glm  <- Lrnr_glm_fast$new()
metalrnr  <- make_learner(Lrnr_nnls)

stack_tut     <- make_learner(Stack, lrnr_mean, lrnr_glm)
slmod_tut     <- make_learner(Lrnr_sl, learners = stack_tut,
                              loss_function = loss_squared_error,
                              metalearner   = metalrnr)
slmod_bin_tut <- make_learner(Lrnr_sl, learners = stack_tut,
                              loss_function = loss_loglik_binomial,
                              metalearner   = metalrnr)

cfg <- list(
  # Column names that the pipeline expects in the dataset
  cluster_id = "gw_cnum",
  admin1     = "Admin1",
  month      = "gw_month",
  child_flag = "gw_child_flag",   # 1 = child, 0 = woman

  # Two outcomes for the tutorial (the real pipeline has four)
  outcomes = list(

    child_vitA = list(
      tag            = "child_vitA",
      label          = "Vitamin A – children",
      population     = "children",
      child_flag_val = 1L,
      continuous     = "gw_cRBP",          # continuous biomarker (RBP µmol/L)
      binary         = "gw_cVAD",          # 1 = deficient (RBP < 0.70)
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_file_cont = "tut_SL_child_vitA.rds",
      model_file_bin  = "tut_bin_SL_child_vitA.rds",
      table_tag      = "children_vitA",
      scatter_tag    = "children_vitA"
    ),

    women_vitA = list(
      tag            = "women_vitA",
      label          = "Vitamin A – women",
      population     = "women",
      child_flag_val = 0L,
      continuous     = "gw_wRBP",
      binary         = "gw_wVAD",
      cutoff         = 0.70,
      cutoff_dir     = "less",
      model_file_cont = "tut_SL_women_vitA.rds",
      model_file_bin  = "tut_bin_SL_women_vitA.rds",
      table_tag      = "women_vitA",
      scatter_tag    = "women_vitA"
    )
  ),

  # Predictor domains: identified by column-name prefix
  # GW   = data collected in the survey (nutrition, anthropometry …)
  # DHS  = externally linked DHS regional estimates
  # MICS = externally linked MICS regional estimates
  domains = list(
    GW   = list(prefix = "gw_",   extra = NULL),
    DHS  = list(prefix = "dhs_",  extra = NULL),
    MICS = list(prefix = "mics_", extra = NULL)
  ),

  # Patterns that mark a gw_ column as an outcome (leak prevention)
  gw_exclude_patterns = c("RBP", "VAD"),

  K    = 2L,    # CV folds  (real pipeline: 5)
  B    = 5L,    # bootstrap replicates  (real pipeline: 200)
  seed = 42L,

  out_models  = here::here("results", "tutorial", "models"),
  out_tables  = here::here("results", "tutorial", "tables"),
  out_figures = here::here("results", "tutorial", "figures")
)

dir.create(cfg$out_models,  showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_figures, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# §01  SIMULATE DATA  (mirrors 01_load_and_construct.R)
# =============================================================================
# In the real pipeline this section loads Gambia_merged_dataset.rds.
# Here we generate a synthetic dataset with the same column structure.

cat("[01] Simulating data...\n")
set.seed(cfg$seed)

N        <- 500     # individuals
n_clust  <- 20      # survey clusters (gw_cnum)
n_admin1 <- 4       # Admin1 regions

# --- Cluster / Admin1 assignment ---------------------------------------------
# Each cluster belongs to exactly one Admin1 region.
clust_ids    <- 1:n_clust
clust_admin1 <- rep(paste0("Region_", LETTERS[1:n_admin1]),
                    length.out = n_clust)

# Assign each individual to a cluster (unequal cluster sizes via Poisson draw)
ind_cluster  <- sample(clust_ids, N, replace = TRUE,
                       prob = rpois(n_clust, 25) + 1)
ind_admin1   <- clust_admin1[ind_cluster]

# --- Child vs woman flag (roughly 60% children) ------------------------------
gw_child_flag <- rbinom(N, 1, 0.6)

# --- Predictor variables -----------------------------------------------------
# gw_ prefix: two survey-level predictors (wealth index, dietary diversity)
gw_wealth <- rnorm(N)
gw_diet   <- rnorm(N)

# dhs_ prefix: two regional DHS predictors (wasting rate, anaemia prevalence)
# These vary by Admin1 region, not by individual (ecological-level data)
region_wasting  <- setNames(rnorm(n_admin1), paste0("Region_", LETTERS[1:n_admin1]))
region_anaemia  <- setNames(rnorm(n_admin1), paste0("Region_", LETTERS[1:n_admin1]))
dhs_wasting  <- region_wasting[ind_admin1]
dhs_anaemia  <- region_anaemia[ind_admin1]

# mics_ prefix: one regional MICS predictor (exclusive breastfeeding rate)
region_ebf  <- setNames(rnorm(n_admin1), paste0("Region_", LETTERS[1:n_admin1]))
mics_ebf    <- region_ebf[ind_admin1]

# gw_month: season of interview (1–12)
gw_month <- sample(1:12, N, replace = TRUE)

# --- Outcomes ----------------------------------------------------------------
# RBP (Retinol Binding Protein, µmol/L) is the continuous VitA biomarker.
# We generate it as a linear function of the predictors so the SL has
# something real to learn.

# True signal (same for children and women in this simulation for simplicity)
latent <- 0.8 + 0.3 * gw_wealth - 0.2 * dhs_wasting + 0.1 * mics_ebf

# Children's RBP (slightly lower mean than women's)
gw_cRBP <- latent - 0.15 + rnorm(N, sd = 0.25)
gw_cVAD <- as.integer(gw_cRBP < cfg$outcomes$child_vitA$cutoff)

# Women's RBP
gw_wRBP <- latent + 0.10 + rnorm(N, sd = 0.20)
gw_wVAD <- as.integer(gw_wRBP < cfg$outcomes$women_vitA$cutoff)

# --- Assemble data frame -----------------------------------------------------
df <- data.frame(
  dataid        = paste0("sim", seq_len(N)),
  gw_cnum       = ind_cluster,
  Admin1        = ind_admin1,
  gw_month      = gw_month,
  gw_child_flag = gw_child_flag,
  # GW domain predictors
  gw_wealth     = gw_wealth,
  gw_diet       = gw_diet,
  # DHS domain predictors
  dhs_wasting   = dhs_wasting,
  dhs_anaemia   = dhs_anaemia,
  # MICS domain predictor
  mics_ebf      = mics_ebf,
  # Continuous outcomes
  gw_cRBP       = gw_cRBP,
  gw_wRBP       = gw_wRBP,
  # Binary outcomes
  gw_cVAD       = gw_cVAD,
  gw_wVAD       = gw_wVAD
)

cat(sprintf("  Simulated %d individuals, %d clusters, %d Admin1 regions\n",
            N, n_clust, n_admin1))
cat(sprintf("  Child VAD prevalence : %.1f%%\n",
            100 * mean(df$gw_cVAD[df$gw_child_flag == 1])))
cat(sprintf("  Women VAD prevalence : %.1f%%\n",
            100 * mean(df$gw_wVAD[df$gw_child_flag == 0])))

# --- Build per-outcome datasets (mirrors build_dataset() in 01_) -------------
# For each outcome we need a dataset that:
#   • contains the relevant population (children or women)
#   • has non-missing values on the continuous outcome
#   • includes all predictor columns + the binary outcome column

# Identify predictor columns by domain prefix
domain_vars <- lapply(cfg$domains, function(dom) {
  vars <- colnames(df)[grepl(dom$prefix, colnames(df), fixed = TRUE)]
  if (!is.null(dom$extra)) vars <- unique(c(vars, dom$extra[dom$extra %in% colnames(df)]))
  vars
})
names(domain_vars) <- names(cfg$domains)

# Remove outcome-leaking gw_ columns from the GW predictor set
leakage_pat <- paste(cfg$gw_exclude_patterns, collapse = "|")
domain_vars[["GW"]] <- domain_vars[["GW"]][
  !grepl(leakage_pat, domain_vars[["GW"]])
]

# Xvars_full: all domain predictors + meta columns
meta_cols  <- c("dataid", cfg$admin1, cfg$cluster_id, cfg$month)
Xvars_full <- unique(c(meta_cols,
                        domain_vars[["GW"]],
                        domain_vars[["DHS"]],
                        domain_vars[["MICS"]]))
Xvars_full <- Xvars_full[Xvars_full %in% colnames(df)]

cat(sprintf("  Xvars_full: %d columns\n", length(Xvars_full)))
cat(sprintf("    GW predictors  : %s\n", paste(domain_vars$GW,   collapse = ", ")))
cat(sprintf("    DHS predictors : %s\n", paste(domain_vars$DHS,  collapse = ", ")))
cat(sprintf("    MICS predictors: %s\n", paste(domain_vars$MICS, collapse = ", ")))

build_dataset <- function(outcome_cfg, df, Xvars_full) {
  pop_flag <- outcome_cfg$child_flag_val
  cont_col <- outcome_cfg$continuous
  bin_col  <- outcome_cfg$binary

  keep_cols <- unique(c(cont_col, cfg$cluster_id, Xvars_full))
  if (!is.null(bin_col) && bin_col %in% colnames(df))
    keep_cols <- unique(c(keep_cols, bin_col))
  keep_cols <- keep_cols[keep_cols %in% colnames(df)]

  df %>%
    dplyr::filter(.data[[cfg$child_flag]] == pop_flag) %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(.data[[cont_col]]))
}

gw_data_list <- list(
  child_vitA = build_dataset(cfg$outcomes$child_vitA, df, Xvars_full),
  women_vitA = build_dataset(cfg$outcomes$women_vitA, df, Xvars_full)
)

cat(sprintf("  child_vitA dataset: %d rows\n", nrow(gw_data_list$child_vitA)))
cat(sprintf("  women_vitA dataset: %d rows\n", nrow(gw_data_list$women_vitA)))
cat("[01] Done.\n\n")


# =============================================================================
# §02  FIT SUPERLEARNER MODELS  (mirrors 02_fit_sl_models.R)
# =============================================================================
# DHS_SL_clustered (from sl_helpers.R) does:
#   1. Selects Xvars from the dataset
#   2. Removes all-NA and near-zero-variance columns
#   3. Imputes missing values and creates missing-indicator columns
#   4. Screens predictors with a p-value filter (washb_prescreen)
#   5. Applies recipes cleaning (zero-var, near-zero-var, correlation, normalise)
#   6. Creates a cluster-blocked sl3 task (no cluster is split across CV folds)
#   7. Trains the SuperLearner and returns cross-validated predictions

cat("[02] Fitting SuperLearner models...\n")

dir.create(cfg$out_models, showWarnings = FALSE, recursive = TRUE)

safe_fit <- function(label, d, outcome, population, Xvars, sl, outfile) {
  cat(sprintf("  Fitting %-40s ... ", label))
  t0  <- proc.time()
  res <- tryCatch(
    DHS_SL_clustered(
      d          = d,
      Xvars      = Xvars[Xvars %in% colnames(d)],
      outcome    = outcome,
      population = population,
      id         = cfg$cluster_id,
      folds      = cfg$K,
      CV         = FALSE,
      prescreen  = TRUE,
      sl         = sl
    ),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  elapsed <- round((proc.time() - t0)["elapsed"], 1)
  if (!is.null(res)) {
    saveRDS(res, file = outfile)
    cat(sprintf("done (%.1fs)\n", elapsed))
  }
  res
}

# ---- A. Continuous models ---------------------------------------------------
# The SL is fitted on the raw continuous biomarker.
# Cross-validated predictions are used to evaluate RMSE and AUC.
sl_results_cont <- list()

sl_results_cont$child_vitA <- safe_fit(
  label      = "children / Vitamin A (continuous)",
  d          = gw_data_list$child_vitA,
  outcome    = cfg$outcomes$child_vitA$continuous,
  population = cfg$outcomes$child_vitA$population,
  Xvars      = Xvars_full,
  sl         = slmod_tut,
  outfile    = file.path(cfg$out_models, cfg$outcomes$child_vitA$model_file_cont)
)

sl_results_cont$women_vitA <- safe_fit(
  label      = "women / Vitamin A (continuous)",
  d          = gw_data_list$women_vitA,
  outcome    = cfg$outcomes$women_vitA$continuous,
  population = cfg$outcomes$women_vitA$population,
  Xvars      = Xvars_full,
  sl         = slmod_tut,
  outfile    = file.path(cfg$out_models, cfg$outcomes$women_vitA$model_file_cont)
)

# ---- B. Binary models -------------------------------------------------------
# Direct modelling of binary deficiency probability.
# Preferred over thresholding continuous predictions for prevalence estimation.
sl_results_bin <- list()

get_bin_data <- function(d, bin_col) {
  if (!bin_col %in% colnames(d)) return(NULL)
  d_bin <- d[!is.na(d[[bin_col]]), ]
  if (length(unique(d_bin[[bin_col]])) < 2) return(NULL)
  d_bin
}

for (tag in names(cfg$outcomes)) {
  oc    <- cfg$outcomes[[tag]]
  d_tmp <- get_bin_data(gw_data_list[[tag]], oc$binary)
  if (!is.null(d_tmp)) {
    sl_results_bin[[tag]] <- safe_fit(
      label      = paste0(tag, " (binary)"),
      d          = d_tmp,
      outcome    = oc$binary,
      population = oc$population,
      Xvars      = Xvars_full,
      sl         = slmod_bin_tut,
      outfile    = file.path(cfg$out_models, oc$model_file_bin)
    )
  }
}

cat("[02] Done.\n\n")


# =============================================================================
# §03  PREDICT AND AGGREGATE TO ADMIN1  (mirrors 03_predict_and_aggregate_admin1.R)
# =============================================================================
# For each outcome we:
#   • Extract the cross-validated predictions stored in res$res$yhat_full
#   • Compute Admin1 prevalence: mean predicted probability per Admin1 unit
#   • Compare to observed Admin1 prevalence (from the binary column)
#   • Compute CV performance: RMSE (continuous), AUC and Brier (binary)

cat("[03] Aggregating to Admin1 and computing CV performance...\n")

# Helper: apply a threshold to turn continuous predictions into 0/1
apply_threshold <- function(x, cutoff, direction = "less") {
  if (direction == "less") as.integer(x < cutoff) else as.integer(x > cutoff)
}

admin1_prev_list <- list()
cv_perf_rows     <- list()

for (tag in names(cfg$outcomes)) {

  oc       <- cfg$outcomes[[tag]]
  d_orig   <- gw_data_list[[tag]]
  res_cont <- sl_results_cont[[tag]]
  res_bin  <- sl_results_bin[[tag]]

  # -- CV performance --
  # Continuous RMSE and AUC (using continuous prediction as a ranking score)
  perf_cont <- NULL
  if (!is.null(res_cont)) {
    r    <- res_cont$res
    rmse <- sqrt(mean((r$Y - r$yhat_full)^2, na.rm = TRUE))
    ybin <- apply_threshold(r$Y, oc$cutoff, oc$cutoff_dir)
    score <- if (oc$cutoff_dir == "less") -r$yhat_full else r$yhat_full
    auc_c <- tryCatch(as.numeric(pROC::auc(
      pROC::roc(ybin[!is.na(ybin)], score[!is.na(ybin)], quiet = TRUE)
    )), error = function(e) NA_real_)
    brier_c <- mean((apply_threshold(r$yhat_full, oc$cutoff, oc$cutoff_dir) - ybin)^2,
                    na.rm = TRUE)
    perf_cont <- data.frame(outcome = tag, model_type = "continuous",
                            rmse = rmse, auc = auc_c, brier = brier_c)
  }

  # Binary AUC and Brier score
  perf_bin <- NULL
  if (!is.null(res_bin)) {
    r     <- res_bin$res
    auc_b <- tryCatch(as.numeric(pROC::auc(
      pROC::roc(r$Y, r$yhat_full, quiet = TRUE)
    )), error = function(e) NA_real_)
    brier_b <- mean((r$yhat_full - r$Y)^2, na.rm = TRUE)
    perf_bin <- data.frame(outcome = tag, model_type = "binary",
                           rmse = NA_real_, auc = auc_b, brier = brier_b)
  }

  cv_perf_rows[[tag]] <- rbind(perf_cont, perf_bin)

  # -- Admin1 prevalence aggregation --
  # Use binary model predicted probabilities when available
  if (!is.null(res_bin)) {
    pred_src <- res_bin$res %>% dplyr::select(dataid, yhat_pred = yhat_full)
  } else if (!is.null(res_cont)) {
    pred_src <- res_cont$res %>%
      dplyr::mutate(yhat_pred = apply_threshold(yhat_full, oc$cutoff, oc$cutoff_dir)) %>%
      dplyr::select(dataid, yhat_pred)
  } else next

  # Observed deficiency (from binary column or thresholded continuous)
  if (!is.null(oc$binary) && oc$binary %in% colnames(d_orig)) {
    obs_src <- d_orig %>%
      dplyr::select(dataid, Admin1, Y_obs = dplyr::all_of(oc$binary))
  } else {
    obs_src <- d_orig %>%
      dplyr::mutate(Y_obs = apply_threshold(.data[[oc$continuous]], oc$cutoff, oc$cutoff_dir)) %>%
      dplyr::select(dataid, Admin1, Y_obs)
  }

  prev_df <- obs_src %>%
    dplyr::left_join(pred_src, by = "dataid") %>%
    dplyr::filter(!is.na(Admin1), !is.na(Y_obs), !is.na(yhat_pred))

  admin1_tbl <- prev_df %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n               = dplyr::n(),
      obs_prevalence  = mean(Y_obs,     na.rm = TRUE),
      pred_prevalence = mean(yhat_pred, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(outcome = tag, label = oc$label)

  admin1_prev_list[[tag]] <- admin1_tbl

  csv_file <- file.path(cfg$out_tables,
                        sprintf("tut_admin1_prevalence_%s.csv", oc$table_tag))
  readr::write_csv(admin1_tbl, csv_file)
  cat(sprintf("  [%s] Admin1 table saved: %s\n", tag, basename(csv_file)))
  print(admin1_tbl[, c("Admin1", "n", "obs_prevalence", "pred_prevalence")],
        digits = 3, row.names = FALSE)
  cat("\n")
}

# Print CV performance summary
cv_perf_table <- do.call(rbind, cv_perf_rows)
cat("  CV performance:\n")
print(cv_perf_table[, c("outcome", "model_type", "rmse", "auc", "brier")],
      digits = 4, row.names = FALSE)

cv_csv <- file.path(cfg$out_tables, "tut_cv_performance.csv")
readr::write_csv(cv_perf_table, cv_csv)
cat(sprintf("\n  Saved: %s\n", basename(cv_csv)))
cat("[03] Done.\n\n")


# =============================================================================
# §04  BOOTSTRAP UNCERTAINTY  (mirrors 04_bootstrap_uncertainty.R)
# =============================================================================
# Cluster bootstrap: resample survey clusters with replacement, refit the SL
# on each resample, predict on the original data, aggregate to Admin1.
# Repeat B times to get a distribution of Admin1 prevalence estimates.
# 95% CI = 2.5th and 97.5th percentiles of that distribution.
#
# one_bootstrap() is defined in sl_helpers.R.  It uses
#   fit_b$task$nodes$covariates   (not fit_b$Xvars)
# to get the correct set of covariate names — the SL task's covariates never
# include the internally-created Y and id columns.

cat("[04] Bootstrap uncertainty (B =", cfg$B, "replicates)...\n")
library(future.apply)
plan(sequential)   # use plan(multisession) to parallelise

for (tag in names(cfg$outcomes)) {

  oc       <- cfg$outcomes[[tag]]
  d_orig   <- gw_data_list[[tag]]
  res_bin  <- sl_results_bin[[tag]]
  res_cont <- sl_results_cont[[tag]]

  use_binary  <- !is.null(res_bin)
  outcome_col <- if (use_binary) oc$binary else oc$continuous
  sl_obj      <- if (use_binary) slmod_bin_tut else slmod_tut
  model_type  <- if (use_binary) "binary_prob" else "continuous_threshold"

  cat(sprintf("\n  Bootstrap for: %s  (model: %s)\n", tag, model_type))

  if (!outcome_col %in% colnames(d_orig)) {
    cat(sprintf("  [skip] column '%s' not in dataset\n", outcome_col)); next
  }

  d_fit   <- d_orig[!is.na(d_orig[[outcome_col]]), ]
  Xvars_b <- Xvars_full[Xvars_full %in% colnames(d_fit)]

  boot_results <- future.apply::future_lapply(
    seq_len(cfg$B),
    FUN            = one_bootstrap,
    d_boot_orig    = d_fit,
    Xvars_b        = Xvars_b,
    outcome_b      = outcome_col,
    population_b   = oc$population,
    id_col         = cfg$cluster_id,
    K              = cfg$K,
    sl_obj         = sl_obj,
    d_predict      = d_orig,
    cutoff         = oc$cutoff,
    cutoff_dir     = oc$cutoff_dir,
    binary_outcome = use_binary,
    seed_base      = cfg$seed,
    future.seed    = TRUE,
    future.globals = TRUE
  )

  boot_results <- Filter(Negate(is.null), boot_results)
  cat(sprintf("  Valid replicates: %d / %d\n", length(boot_results), cfg$B))

  if (length(boot_results) < 2) {
    cat("  [warn] Too few valid replicates; skipping CI output.\n"); next
  }

  # --- Admin1 CIs ---
  all_a1 <- do.call(rbind, lapply(seq_along(boot_results), function(i) {
    df_i <- boot_results[[i]]$admin1; df_i$rep <- i; df_i
  }))

  admin1_ci <- all_a1 %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      boot_mean = mean(prev,            na.rm = TRUE),
      ci_lo     = quantile(prev, 0.025, na.rm = TRUE),
      ci_hi     = quantile(prev, 0.975, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(outcome = tag, model_type = model_type)

  # Merge observed prevalence from §03
  if (!is.null(admin1_prev_list[[tag]])) {
    admin1_ci <- admin1_ci %>%
      dplyr::left_join(
        admin1_prev_list[[tag]] %>% dplyr::select(Admin1, obs_prevalence, pred_prevalence),
        by = "Admin1"
      )
  }

  a1_csv <- file.path(cfg$out_tables,
                      sprintf("tut_admin1_ci_%s.csv", oc$table_tag))
  readr::write_csv(admin1_ci, a1_csv)
  cat(sprintf("  Admin1 CI saved: %s\n", basename(a1_csv)))
  print(admin1_ci[, c("Admin1", "obs_prevalence", "boot_mean", "ci_lo", "ci_hi")],
        digits = 3, row.names = FALSE)

  # --- National CI ---
  nat_vec <- sapply(boot_results, function(x) x$national)
  nat_ci  <- data.frame(outcome    = tag,
                        model_type = model_type,
                        boot_mean  = mean(nat_vec,            na.rm = TRUE),
                        ci_lo      = quantile(nat_vec, 0.025, na.rm = TRUE),
                        ci_hi      = quantile(nat_vec, 0.975, na.rm = TRUE))

  nat_csv <- file.path(cfg$out_tables,
                       sprintf("tut_national_ci_%s.csv", oc$table_tag))
  readr::write_csv(nat_ci, nat_csv)
  cat(sprintf("  National prevalence: %.1f%% [%.1f%%, %.1f%%]\n",
              nat_ci$boot_mean * 100, nat_ci$ci_lo * 100, nat_ci$ci_hi * 100))
}

cat("\n[04] Done.\n\n")


# =============================================================================
# §05  DOMAIN ABLATION  (mirrors 05_domain_ablation.R)
# =============================================================================
# For each domain (GW, DHS, MICS) we:
#   1. Refit the SL with that domain's columns removed
#   2. Compare CV performance to the full model
#   3. Delta AUC = AUC_full - AUC_reduced  (positive → domain helps)
#
# Because the tutorial has only a few predictors per domain the differences
# will be small, but the mechanics are identical to the real pipeline.

cat("[05] Domain ablation...\n")

extract_metrics <- function(res_obj, cutoff, cutoff_dir, model_type) {
  if (is.null(res_obj)) return(NULL)
  r <- res_obj$res; Y <- r$Y; Yhat <- r$yhat_full
  if (model_type == "continuous") {
    Ybin  <- apply_threshold(Y, cutoff, cutoff_dir)
    score <- if (cutoff_dir == "less") -Yhat else Yhat
    auc   <- tryCatch(as.numeric(pROC::auc(
      pROC::roc(Ybin[!is.na(Ybin)], score[!is.na(Ybin)], quiet = TRUE)
    )), error = function(e) NA_real_)
    brier <- mean((apply_threshold(Yhat, cutoff, cutoff_dir) - Ybin)^2, na.rm = TRUE)
    rmse  <- sqrt(mean((Y - Yhat)^2, na.rm = TRUE))
  } else {
    auc   <- tryCatch(as.numeric(pROC::auc(
      pROC::roc(Y, Yhat, quiet = TRUE)
    )), error = function(e) NA_real_)
    brier <- mean((Yhat - Y)^2, na.rm = TRUE)
    rmse  <- NA_real_
  }
  data.frame(rmse = rmse, auc = auc, brier = brier)
}

ablation_rows <- list()

for (tag in names(cfg$outcomes)) {

  oc       <- cfg$outcomes[[tag]]
  d_orig   <- gw_data_list[[tag]]
  use_bin  <- !is.null(sl_results_bin[[tag]])
  res_full <- if (use_bin) sl_results_bin[[tag]] else sl_results_cont[[tag]]
  mtype    <- if (use_bin) "binary" else "continuous"
  out_col  <- if (use_bin) oc$binary else oc$continuous
  sl_obj   <- if (use_bin) slmod_bin_tut else slmod_tut

  if (is.null(res_full)) next

  full_m <- extract_metrics(res_full, oc$cutoff, oc$cutoff_dir, mtype)
  cat(sprintf("\n  Ablation: %s  |  Full AUC = %.3f\n", tag, full_m$auc))

  d_fit   <- d_orig[!is.na(d_orig[[out_col]]), ]
  Xvars_b <- Xvars_full[Xvars_full %in% colnames(d_fit)]

  for (dom_name in names(cfg$domains)) {
    dom_cols    <- domain_vars[[dom_name]]
    remove_cols <- intersect(dom_cols, Xvars_b)

    if (length(remove_cols) == 0) {
      cat(sprintf("    Domain %-6s: [no columns]\n", dom_name)); next
    }

    Xvars_red <- setdiff(Xvars_b, remove_cols)
    Xvars_red <- Xvars_red[Xvars_red %in% colnames(d_fit)]

    if (length(Xvars_red) == 0) {
      cat(sprintf("    Domain %-6s: [no remaining predictors]\n", dom_name)); next
    }

    res_red <- tryCatch(
      DHS_SL_clustered(d = d_fit, Xvars = Xvars_red, outcome = out_col,
                       population = oc$population, id = cfg$cluster_id,
                       folds = cfg$K, CV = FALSE, prescreen = TRUE, sl = sl_obj),
      error = function(e) NULL
    )

    red_m <- extract_metrics(res_red, oc$cutoff, oc$cutoff_dir, mtype)

    if (!is.null(red_m)) {
      delta_auc <- full_m$auc - red_m$auc
      cat(sprintf("    Domain %-6s (remove %d cols): AUC_reduced=%.3f  ΔAUC=%.3f%s\n",
                  dom_name, length(remove_cols), red_m$auc, delta_auc,
                  if (delta_auc > 0.01) "  ← informative" else ""))

      ablation_rows[[paste(tag, dom_name, sep = "_")]] <- data.frame(
        outcome       = tag,
        model_type    = mtype,
        domain        = dom_name,
        n_dom_cols    = length(remove_cols),
        auc_full      = full_m$auc,
        auc_reduced   = red_m$auc,
        delta_auc     = delta_auc,
        brier_full    = full_m$brier,
        brier_reduced = red_m$brier,
        delta_brier   = red_m$brier - full_m$brier
      )
    } else {
      cat(sprintf("    Domain %-6s: [reduced model failed]\n", dom_name))
    }
  }
}

if (length(ablation_rows) > 0) {
  abl_tbl <- do.call(rbind, ablation_rows)
  abl_csv <- file.path(cfg$out_tables, "tut_domain_ablation.csv")
  readr::write_csv(abl_tbl, abl_csv)
  cat(sprintf("\n  Ablation table saved: %s\n", basename(abl_csv)))
}

cat("\n[05] Done.\n\n")


# =============================================================================
# SUMMARY
# =============================================================================
cat("=============================================================\n")
cat("  Tutorial complete.  Output files:\n")
for (f in list.files(cfg$out_tables, full.names = TRUE))
  cat(sprintf("    %s\n", f))
cat("=============================================================\n\n")
