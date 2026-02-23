


# function to calculate model performances
evaluate_sl_performance <- function(res, outcome){
  # Ensure numeric vectors
  Y <- as.numeric(res$Y)
  Yhat <- as.numeric(res$yhat_full)

  # Baseline MSE
  mse_mean <- as.numeric(res$cv_risk_w_sl_revere[res$cv_risk_w_sl_revere$learner == "Lrnr_mean", "MSE"])

  # Model MSE
  mse <- as.numeric(res$cv_risk_w_sl_revere[res$cv_risk_w_sl_revere$learner == "SuperLearner", "MSE"])

  # RMSE and MAE
  rmse <- sqrt(mse)
  mae <- mean(abs(Y - Yhat))

  # Variance of outcome
  varY <- var(Y)

  # Standard R2 (variance benchmark)
  R2 <- 1 - mse / varY

  # R2 relative to mean predictor (cross-validated SL often uses this)
  R2_baseline <- 1 - mse / mse_mean

  # Relative improvement over the baseline MSE
  rel_improve <- 1 - (mse / mse_mean)

  # Correlation between predicted and observed
  corr <- as.numeric(cor.test(Y, Yhat)$estimate)
  corr.pval <- as.numeric(cor.test(Y, Yhat)$p.value)

  # Average prediction bias
  bias <- mean(Yhat - Y)

  # Compile results in a list
  out <- data.frame(
    outcome = outcome,
    mse_model = mse,
    mse_mean = mse_mean,
    rmse = rmse,
    mae = mae,
    outcome_variance = varY,
    R2_variance_based = R2,
    R2_vs_mean = R2_baseline,
    relative_improvement = rel_improve*100,
    correlation = corr,
    correlation_pvalue = corr.pval,
    mean_bias = bias
  )

  return(out)
}




makeVlist <- function(dta) {
  labels <- sapply(dta, function(x) attr(x, "label"))
  tibble(name = names(labels),
         label = as.character(labels))
}


clean_DHS <- function(dfull){
  HR_clean <- clean_HR(dfull$HRdata)
  PR_clean <- clean_PR(dfull$PRdata)

  # Merge household data into the PR data:
  PR_merged <- PR_clean %>%
    left_join(HR_clean, by = c("cluster","hh"))

  IR_clean <- clean_IR(dfull$IRdata)

  # Merge IR (women’s) data into PR_merged for female lines:
  PR_IR_merged <- PR_merged %>%
    left_join(IR_clean, by = c("cluster","hh","line"))
  head(PR_IR_merged)

  table(PR_IR_merged$anemia_cat)

  d <- PR_IR_merged %>%
    mutate(mod_sev_anemia=case_when(anemia_cat==1 | anemia_cat==2 ~ 1,
                                    anemia_cat==3 | anemia_cat==4~ 0,
                                    TRUE ~ NA)) %>%
    filter(!is.na(mod_sev_anemia))

  if(nrow(d)==0){
    cat("No anemia data\n")
  }else{
    return(d)
  }

}



#fix getDHSdata functions to correctly select the country instead of all countries
getDHSdata <- function(country, indicator = NULL, Recode = NULL, year){
  IR_Individual <- c("ancvisit4+", "RH_ANCN_W_N4P", "womananemia",
                     "AN_ANEM_W_ANY", "unmet_family", "FP_NADA_W_UNT", "FP_CUSA_W_MOD",
                     "AN_NUTS_W_THN")
  PR_Household_Member <- c("CN_ANMC_C_ANY", "CN_NUTS_C_WH2",
                           "wasting", "CN_NUTS_C_HA2", "stunting", "ML_PMAL_C_RDT",
                           "WS_TLET_H_IMP", "WS_TLET_P_BAS", "WS_SRCE_P_BAS")
  KR_Children <- c("CH_DIAT_C_ORT", "DPT3", "CH_VACC_C_DP3",
                   "CH_VACC_C_DP1", "CH_VACC_C_BAS", "CH_VACC_C_NON", "CN_BRFS_C_EXB",
                   "CH_VACC_C_MSL", "PCV3", "RotaC1")
  BRdata_Birth <- c("RH_DELA_C_SKP", "CM_ECMR_C_NNR", "nmr")
  HRdata_Household <- c("ML_NETP_H_IT2")
  HIV <- c("HA_HIVP_B_HIV")
  indicator <- indicator
  if (is.null(indicator) & is.null(Recode)) {
    Type <- NULL
  }
  else if (!is.null(Recode)) {
    Type = Recode
  }
  else if (indicator %in% IR_Individual) {
    Type <- c("Individual Recode")
  }
  else if (indicator %in% PR_Household_Member) {
    Type <- c("Household Member Recode")
  }
  else if (indicator %in% KR_Children) {
    Type <- c("Children's Recode")
  }
  else if (indicator %in% BRdata_Birth) {
    Type <- c("Births Recode")
  }
  else if (indicator %in% HRdata_Household) {
    Type <- c("Household Recode")
  }
  else if (indicator %in% HIV) {
    Type <-
      c("Individual Recode", "Men's Recode", "HIV Test Results Recode")
  }
  else {
    Type <- NULL
  }
  if (!is.null(Type)) {
    message(paste(Type, "is used.\n\n"))
  }
  else {
    message("All DHS files are downloaded.\n\n")
  }
  CountryName <- stringr::str_to_title(country)
  country_vec <- rdhs::dhs_countries()$CountryName == CountryName
  dhs_countries_data <- dhs_countries()
  countryId <- dhs_countries_data[country_vec, ]
  potential_surveys <-
    rdhs::dhs_datasets(countryIds = countryId$DHS_CountryCode, surveyYear = year) %>% dplyr::filter(FileFormat == "Stata dataset (.dta)")
  if (length(Type) == 1) {
    surveys <- potential_surveys %>% dplyr::filter(FileType ==
                                                     c(Type))
    data.paths.tmp <- get_datasets(surveys[surveys$SurveyYear ==
                                             year,]$FileName, clear_cache = T)
    Rdata <- readRDS(paste0(data.paths.tmp))
    return(Rdata)
  }
  else if (length(Type) > 1) {
    surveys <- potential_surveys %>% dplyr::filter(FileType %in%
                                                     c(Type))
    data.paths.tmp <- get_datasets(surveys[surveys$SurveyYear ==
                                             year,]$FileName, clear_cache = T)
    all = list()
    listname = surveys$FileType
    for (i in 1:length(Type)) {
      all[[listname[i]]] <- readRDS(paste0(data.paths.tmp[i]))
    }
    return(all)
  }
  else if (is.null(Type)) {
    all <- NULL
    list <- c(
      "Men's Recode",
      "Household Member Recode",
      "Children's Recode",
      "Births Recode",
      "Couples' Recode",
      "Household Recode",
      "Individual Recode"
    )
    listname <- c("MRdata",
                  "PRdata",
                  "KRdata",
                  "BRdata",
                  "CRdata",
                  "HRdata",
                  "IRdata")
    for (i in 1:length(list)) {
      Type <- list[i]
      surveys <- potential_surveys %>% dplyr::filter(FileType ==
                                                       c(Type))
      if (dim(surveys)[1] == 0) {

      }
      else {
        data.paths.tmp <- get_datasets(surveys[surveys$SurveyYear ==
                                                 year,]$FileName, clear_cache = T)
        Rdata <- readRDS(paste0(data.paths.tmp))
        all[[listname[i]]] <- Rdata
      }
    }
    return(all)
  }
}

getDHSgeo <- function(country, year){
  CountryName <- stringr::str_to_title(country)
  country_vec <- rdhs::dhs_countries()$CountryName == CountryName
  dhs_countries_data <- dhs_countries()
  countryId <- dhs_countries_data[country_vec, ]
  surveys <- rdhs::dhs_datasets(countryIds = countryId$DHS_CountryCode,
                                surveyYear = year) %>% dplyr::filter(FileType == "Geographic Data")
  data.paths.tmp <- get_datasets(surveys[surveys$SurveyYear == year, ]$FileName, clear_cache = T)
  geo <- readRDS(paste0(data.paths.tmp))
  return(geo)
}




# auc_from_sl3  <- function(res,
#                           pred_source = c("auto","cv","full","sl_fit"),
#                           weights = NULL,
#                           ci = TRUE,
#                           ci_method = c("delong","bootstrap"),
#                           quiet = TRUE) {
#   stopifnot(all(c("Y","sl_fit") %in% names(res)))
#   pred_source <- match.arg(pred_source)
#   ci_method   <- match.arg(ci_method)
#
#   # --- pull Y and coerce to {0,1}
#   Y <- unclass(res$Y)
#   if (is.factor(Y)) {
#     if (length(levels(Y)) != 2) stop("Outcome factor must have exactly 2 levels.")
#     # assume the *second* level is the event (1) if not numeric
#     Y <- as.integer(Y == levels(Y)[2])
#   } else if (is.logical(Y)) {
#     Y <- as.integer(Y)
#   } else {
#     Y <- as.integer(Y)
#   }
#
#   # --- choose predictions
#   get_preds_from_fit <- function(fit) {
#     # If res$sl_fit is a Lrnr_cv fit object, its predict() are OOF CV preds
#     tryCatch(as.numeric(fit$predict()), error = function(e) NULL)
#   }
#
#   preds <- NULL
#   if (pred_source %in% c("auto","cv","sl_fit")) {
#     preds <- get_preds_from_fit(res$sl_fit)
#   }
#   if (is.null(preds) && pred_source %in% c("auto","full") && !is.null(res$yhat_full)) {
#     preds <- as.numeric(res$yhat_full)
#   }
#   if (is.null(preds))
#     stop("Could not find predictions. Supply CV predictions via sl_fit (Lrnr_cv) or provide res$yhat_full.")
#
#   # --- clean pairs (drop NAs), clamp to [0,1] if needed
#   ok <- is.finite(preds) & !is.na(Y)
#   Y <- Y[ok]; preds <- preds[ok]
#   preds <- pmin(pmax(preds, 0), 1)
#
#   # --- weighted vs unweighted AUC
#   if (!is.null(weights)) {
#     # Use yardstick for weighted ROC AUC
#     if (!requireNamespace("yardstick", quietly = TRUE))
#       stop("Install yardstick for weighted AUC: install.packages('yardstick')")
#     w <- as.numeric(weights)[ok]
#     df_eval <- data.frame(
#       truth   = factor(Y, levels = c(0,1)),
#       .pred_1 = preds,
#       w       = w
#     )
#     est <- yardstick::roc_auc(df_eval, truth = truth, .pred_1, case_weights = w)
#     return(list(
#       auc     = as.numeric(est$.estimate),
#       method  = "ROC AUC (weighted)",
#       preds   = preds,
#       Y       = Y
#     ))
#   } else {
#     # Unweighted AUC + optional DeLong CI via pROC
#     if (!requireNamespace("pROC", quietly = TRUE))
#       stop("Install pROC for ROC AUC: install.packages('pROC')")
#
#     roc_obj <- pROC::roc(response = Y, predictor = preds, quiet = quiet, direction = "<")
#     auc_val <- as.numeric(pROC::auc(roc_obj))
#
#     if (!ci) {
#       return(list(
#         auc    = auc_val,
#         method = "ROC AUC (unweighted)",
#         roc    = roc_obj,
#         preds  = preds,
#         Y      = Y
#       ))
#     }
#
#     if (ci_method == "delong") {
#       ci_vec <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
#     } else {
#       # simple bootstrap CI (faster default B=200; adjust as desired)
#       ci_vec <- as.numeric(pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 200))
#     }
#
#     return(data.frame(
#       auc     = auc_val,
#       ci_low  = ci_vec[1],
#       ci_mid  = ci_vec[2],
#       ci_high = ci_vec[3],
#       method  = sprintf("ROC AUC (unweighted, %s CI)", ci_method)
#     ))
#   }
# }


auc_from_sl3  <- function(res,
                          pred_source = c("auto","cv","full","sl_fit"),
                          ci = TRUE,
                          ci_method = c("delong","bootstrap"),
                          quiet = TRUE) {
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Install pROC: install.packages('pROC')")
  if (!requireNamespace("PRROC", quietly = TRUE))
    stop("Install PRROC: install.packages('PRROC')")

  stopifnot(all(c("Y","sl_fit") %in% names(res)))
  pred_source <- match.arg(pred_source)
  ci_method   <- match.arg(ci_method)

  # outcome -> {0,1}
  Y <- unclass(res$Y)
  if (is.factor(Y)) {
    if (length(levels(Y)) != 2) stop("Outcome factor must have exactly 2 levels.")
    Y <- as.integer(Y == levels(Y)[2])
  } else if (is.logical(Y)) {
    Y <- as.integer(Y)
  } else {
    Y <- as.integer(Y)
  }

  # predictions (auto prefers CV preds from sl_fit, fallback to yhat_full)
  get_preds_from_fit <- function(fit) tryCatch(as.numeric(fit$predict()), error = function(e) NULL)
  preds <- NULL
  if (pred_source %in% c("auto","cv","sl_fit")) preds <- get_preds_from_fit(res$sl_fit)
  if (is.null(preds) && pred_source %in% c("auto","full") && !is.null(res$yhat_full))
    preds <- as.numeric(res$yhat_full)
  if (is.null(preds)) stop("No predictions found: supply Lrnr_cv sl_fit or res$yhat_full.")

  # clean/clamp
  ok <- is.finite(preds) & !is.na(Y)
  Y <- Y[ok]; preds <- pmin(pmax(preds[ok], 0), 1)

  # ROC AUC (+ optional CI)
  roc_obj <- pROC::roc(response = Y, predictor = preds, quiet = quiet, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  if (ci) {
    ci_vec <- if (ci_method == "delong") {
      as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    } else {
      as.numeric(pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 200))
    }
  }

  # PR-AUC (good for rare outcomes)
  pos <- preds[Y == 1]; neg <- preds[Y == 0]
  pr_auc <- if (length(pos) > 0 && length(neg) > 0) {
    PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)$auc.integral
  } else NA_real_

  # Brier + Brier Skill Score
  brier <- mean((preds - Y)^2)
  prevalence <- mean(Y)
  brier_ref <- prevalence * (1 - prevalence)   # null = predict prevalence everywhere
  brier_skill <- 1 - (brier / brier_ref)       # >0 beats null; 0 ~ null; <0 worse

  out <- data.frame(
    roc_auc         = auc_val,
    roc_auc_ci_low  = if (ci) ci_vec[1] else NA_real_,
    roc_auc_ci_mid  = if (ci) ci_vec[2] else NA_real_,
    roc_auc_ci_high = if (ci) ci_vec[3] else NA_real_,
    roc_method      = sprintf("ROC AUC (unweighted, %s CI)", if (ci) ci_method else "no CI"),
    pr_auc          = pr_auc,
    pr_baseline     = prevalence,  # PR-AUC baseline ≈ prevalence
    brier           = brier,
    brier_ref       = brier_ref,
    brier_skill     = brier_skill,
    prevalence      = prevalence,
    n               = length(Y)
  )
  out
}





plot_calibration <- function(res, bins = 10) {
  library(dplyr)
  library(ggplot2)

  Y <- res$Y
  preds <- if (!is.null(res$yhat_full)) as.numeric(res$yhat_full) else res$sl_fit$predict()

  ok <- is.finite(preds) & !is.na(Y)
  Y <- Y[ok]; preds <- preds[ok]
  preds <- pmin(pmax(preds, 0), 1)

  calib_data <- data.frame(
    yhat = preds,
    Y    = Y
  ) %>%
    mutate(bin = ntile(yhat, bins)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(yhat),
      obs_rate  = mean(Y),
      n         = n(),
      .groups   = "drop"
    )

  ggplot(calib_data, aes(x = mean_pred, y = obs_rate, size = n)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    scale_size_continuous(name = "Bin size") +
    labs(
      x = "Mean predicted probability",
      y = "Observed prevalence",
      title = "Calibration plot"
    ) +
    theme_minimal()
}




auc_pr_brier_from_sl3 <- function(
    res,
    pred_source = c("auto","cv","full","sl_fit"),
    ci = TRUE,
    B = 500,
    cluster_id = NULL,     # vector same length as Y; if given, bootstrap at cluster level
    seed = NULL,
    quiet = TRUE
) {
  # deps
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Install pROC: install.packages('pROC')")
  if (!requireNamespace("PRROC", quietly = TRUE))
    stop("Install PRROC: install.packages('PRROC')")

  stopifnot(all(c("res","sl_fit") %in% names(res)))
  pred_source <- match.arg(pred_source)

  # --- outcome to {0,1}
  Y <- unclass(res$res$Y)
  if (is.factor(Y)) {
    if (length(levels(Y)) != 2) stop("Outcome factor must have exactly 2 levels.")
    Y <- as.integer(Y == levels(Y)[2])
  } else if (is.logical(Y)) {
    Y <- as.integer(Y)
  } else {
    Y <- as.integer(Y)
  }

  # --- predictions (auto prefers CV preds from sl_fit; fallback to yhat_full)
  get_preds_from_fit <- function(fit) tryCatch(as.numeric(fit$predict()), error = function(e) NULL)
  preds <- NULL
  # if (pred_source %in% c("auto","cv","sl_fit")) preds <- get_preds_from_fit(res$sl_fit)
  # if (is.null(preds) && pred_source %in% c("auto","full") && !is.null(res$yhat_full))
    preds <- as.numeric(res$res$yhat_full)
  if (is.null(preds)) stop("No predictions found: supply Lrnr_cv sl_fit or res$yhat_full.")

  # --- clean/clamp
  ok <- is.finite(preds) & !is.na(Y)
  Y <- Y[ok]; preds <- pmin(pmax(preds[ok], 0), 1)
  n <- length(Y)
  prevalence <- mean(Y)

  # --- point estimates
  roc_obj <- pROC::roc(response = Y, predictor = preds, quiet = quiet, direction = "<")
  roc_auc <- as.numeric(pROC::auc(roc_obj))

  pr_auc <- {
    pos <- preds[Y == 1]; neg <- preds[Y == 0]
    if (length(pos) > 0 && length(neg) > 0) {
      PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)$auc.integral
    } else NA_real_
  }

  brier <- mean((preds - Y)^2)
  brier_ref <- prevalence * (1 - prevalence)
  brier_skill <- if (brier_ref > 0) 1 - (brier / brier_ref) else NA_real_

  # --- bootstrap CIs (shared replicates for ROC, PR, Brier, BSS)
  roc_ci <- pr_ci <- brier_ci <- brier_skill_ci <- c(NA_real_, NA_real_)
  roc_boot <- pr_boot <- brier_boot <- brier_skill_boot <- NULL

  if (ci) {
    if (!is.null(seed)) set.seed(seed)

    if (is.null(cluster_id)) {
      # i.i.d. bootstrap
      boot_one <- function() {
        idx <- sample.int(n, replace = TRUE)
        yb  <- Y[idx]; pb <- preds[idx]

        # ROC
        roc_b <- as.numeric(pROC::auc(pROC::roc(yb, pb, quiet = TRUE, direction = "<")))

        # PR
        pos_b <- pb[yb == 1]; neg_b <- pb[yb == 0]
        pr_b  <- if (length(pos_b) > 0 && length(neg_b) > 0)
          PRROC::pr.curve(scores.class0 = pos_b, scores.class1 = neg_b, curve = FALSE)$auc.integral
        else NA_real_

        # Brier & BSS
        brier_b <- mean((pb - yb)^2)
        p_b     <- mean(yb)
        bref_b  <- p_b * (1 - p_b)
        bss_b   <- if (bref_b > 0) 1 - (brier_b / bref_b) else NA_real_

        c(roc_b, pr_b, brier_b, bss_b)
      }
      boots <- replicate(B, boot_one())
    } else {
      # cluster bootstrap
      stopifnot(length(cluster_id) == length(res$Y))
      cluster_id <- cluster_id[ok]
      clust <- unique(cluster_id); G <- length(clust)

      boot_one <- function() {
        samp_cl <- sample(clust, size = G, replace = TRUE)
        idx <- unlist(lapply(samp_cl, function(g) which(cluster_id == g)), use.names = FALSE)
        yb  <- Y[idx]; pb <- preds[idx]

        roc_b <- as.numeric(pROC::auc(pROC::roc(yb, pb, quiet = TRUE, direction = "<")))
        pos_b <- pb[yb == 1]; neg_b <- pb[yb == 0]
        pr_b  <- if (length(pos_b) > 0 && length(neg_b) > 0)
          PRROC::pr.curve(scores.class0 = pos_b, scores.class1 = neg_b, curve = FALSE)$auc.integral
        else NA_real_

        brier_b <- mean((pb - yb)^2)
        p_b     <- mean(yb)
        bref_b  <- p_b * (1 - p_b)
        bss_b   <- if (bref_b > 0) 1 - (brier_b / bref_b) else NA_real_

        c(roc_b, pr_b, brier_b, bss_b)
      }
      boots <- replicate(B, boot_one())
    }

    roc_boot        <- boots[1, ]
    pr_boot         <- boots[2, ]
    brier_boot      <- boots[3, ]
    brier_skill_boot<- boots[4, ]

    roc_ci         <- stats::quantile(roc_boot,         probs = c(0.025, 0.975), na.rm = TRUE)
    pr_ci          <- stats::quantile(pr_boot,          probs = c(0.025, 0.975), na.rm = TRUE)
    brier_ci       <- stats::quantile(brier_boot,       probs = c(0.025, 0.975), na.rm = TRUE)
    brier_skill_ci <- stats::quantile(brier_skill_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  }

  data.frame(
    n               = n,
    prevalence      = prevalence,
    # ROC
    roc_auc         = roc_auc,
    roc_ci_low      = if (ci) as.numeric(roc_ci[1]) else NA_real_,
    roc_ci_high     = if (ci) as.numeric(roc_ci[2]) else NA_real_,
    # PR
    pr_auc          = pr_auc,
    pr_baseline     = prevalence,
    pr_gain         = if (!is.na(pr_auc)) pr_auc / prevalence else NA_real_,
    pr_norm_auc     = if (!is.na(pr_auc)) (pr_auc - prevalence) / (1 - prevalence) else NA_real_,
    pr_ci_low       = if (ci) as.numeric(pr_ci[1]) else NA_real_,
    pr_ci_high      = if (ci) as.numeric(pr_ci[2]) else NA_real_,
    # Brier
    brier           = brier,
    brier_ref       = brier_ref,
    brier_skill     = brier_skill,
    brier_ci_low    = if (ci) as.numeric(brier_ci[1]) else NA_real_,
    brier_ci_high   = if (ci) as.numeric(brier_ci[2]) else NA_real_,
    brier_skill_ci_low  = if (ci) as.numeric(brier_skill_ci[1]) else NA_real_,
    brier_skill_ci_high = if (ci) as.numeric(brier_skill_ci[2]) else NA_real_#,
    # # bootstrap draws
    # roc_boot        = if (ci) roc_boot else NULL,
    # pr_boot         = if (ci) pr_boot else NULL,
    # brier_boot      = if (ci) brier_boot else NULL,
    # brier_skill_boot= if (ci) brier_skill_boot else NULL,
    # method          = if (is.null(cluster_id)) sprintf("Bootstrap (%d i.i.d. resamples)", B)
    # else sprintf("Cluster bootstrap (%d resamples over %d clusters)", B, length(unique(cluster_id)))
  )
}



calc_mn_importance <- function(sl_fit, eval_fun = loss_loglik_binomial,
                               importance.metric="ratio", n_vars=20,
                               covariate_groups=NULL){

  set.seed(983)

  fit = sl_fit
  type = "permute"
  fold_number = "validation"
  #importance_metric = c("difference", "ratio"),
  importance_metric = "difference"

  task <- fit$training_task
  d <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  weights <- task$weights
  if (!is.null(covariate_groups)) {
    if (!is.list(covariate_groups)) {
      stop("Covariate groups must be a list.")
    }
    if (!all(unlist(covariate_groups) %in% X)) {
      stop("Groups contain covariates that are not in the task's covariates.")
    }
    # if (!all(X %in% unlist(covariate_groups))) {
    #   missingX <- as.list(X[which(!X %in% unlist(covariate_groups))])
    #   names(missingX) <- X[which(!X %in% unlist(covariate_groups))]
    #   covariate_groups <- c(covariate_groups, missingX)
    # }
    if(any(is.null(names(covariate_groups))) | any(names(covariate_groups) =="") | any(is.na(names(covariate_groups)))){
      no_name <- unique(which(is.null(names(covariate_groups)) |names(covariate_groups) == "" | is.na(names(covariate_groups))))
      if (any(sapply(covariate_groups[unique(no_name)],length)) != 1) {
        stop("Covariate groups with more than one covariate must be named.")
      }else{
        if(all(sapply(covariate_groups[unique(no_name)],length)) == 1){
          names(covariate_groups[no_name]) <- unlist(covariate_groups[no_name])
        }
      }
    }
    X <- covariate_groups
  }else{
    names(X) <- X
  }

  original_pred <- fit$predict_fold(task, fold_number = fold_number)
  original_eval <- eval_fun(original_pred, Y)
  if (!is.null(attr(original_eval, "loss")) && !attr(original_eval,"loss")) {
    original_risk <- original_eval
  }else {
    original_losses <- original_eval
    original_risk <- weighted.mean(original_losses, weights)
  }

  res_list <- lapply(X, function(x){
    if (type == "permute"){
      perm <- sample(1:nrow(d), nrow(d))
      x_perm <- d[perm, x, with = FALSE]
      data.table::setnames(x_perm, x)
      x_perm_name <- task$add_columns(x_perm)
      task_x_perm <- task$next_in_chain(column_names = x_perm_name)
      x_perm_pred <- fit$predict_fold(task_x_perm, fold_number = fold_number)
      x_perm_eval <- eval_fun(x_perm_pred, Y)

      if (!is.null(attr(original_eval, "loss")) && !attr(original_eval,"loss")) {
        no_x_risk <- x_perm_eval
      }else{
        no_x_losses <- x_perm_eval
        no_x_risk <- weighted.mean(no_x_losses, weights)
      }
    }

    if (importance_metric == "ratio") {
      result <- no_x_risk/original_risk
    }
    else if (importance_metric == "difference") {
      result <- no_x_risk - original_risk
    }
    return(result)
  })

  result <- data.table::data.table(covariate = names(X), metric = unlist(res_list))
  if (!is.null(covariate_groups)) {
    colnames(result)[1] <- "covariate_group"
  }
  data.table::setorderv(result, cols = "metric", order = -1L)
  metric_name <- paste0("risk_", importance_metric)
  if (!is.null(attr(original_eval, "name"))) {
    metric_name <- gsub("risk", attr(original_eval, "name"),
                        metric_name)
  }
  colnames(result)[2] <- metric_name


  varimp20=result %>% arrange(.[,2]) %>% tail(n=n_vars)
  p=importance_plot(x = varimp20)
  return(list(varimp=result, varimp20=varimp20, p=p))
}



