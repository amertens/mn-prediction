# =============================================================================
# scripts/protocol_v2/67_january_individual_forensics.R   [IL-02a]
#
# WHAT THE JANUARY 2026 PERSON-LEVEL FIGURE MEASURED
#
# The January Ghana deck's "Model results by outcome and predictor set" figure
# (Brier skill 3-67 percent, survey-only AUC 0.87-1.00) was drawn from
# results/res_full_bin_GW_Ghana_SL_all.rds in the mn-proxies repository, built
# on 22 August 2025 by src/Ghana/model_performance_bin.R from the eighteen sl3
# fits results/models/res_full_bin_GW_Ghana_SL_*.rds. This script re-scores
# those same fits three ways and puts the three side by side:
#   figure     the metrics the figure plotted
#   insample   sl_fit$predict(): the full-data fit predicting its own training
#              rows. auc_pr_brier_from_sl3(pred_source = "auto") called this
#              first; its comment calls it "CV preds", but on an sl3 Lrnr_sl a
#              bare $predict() is resubstitution
#   cv         the fits' own out-of-fold predictions (res$yhat_full), from the
#              10-fold cluster-blocked CV the fits were trained with
# It also records whether any cluster was split across folds, and the
# covariates each fit received (to count identifiers, sampling-design columns
# and blood-derived columns among them).
#
# Reads files outside this repository (set MN_PROXIES to relocate). Needs sl3
# and every learner package the January library used, attached, or
# $predict() fails on the HAL learners.
#
#   Rscript scripts/protocol_v2/67_january_individual_forensics.R
# -> results/tables/protocol_v2/il02_january_figure_vs_cv.csv
# -> results/tables/protocol_v2/il02_january_covariates.csv
# =============================================================================
suppressPackageStartupMessages({
  library(sl3); library(data.table); library(hal9001); library(ranger); library(glmnet)
  library(xgboost); library(randomForest); library(polspline); library(arm); library(gam)
})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
MNP <- Sys.getenv("MN_PROXIES", "C:/Users/andre/OneDrive/Documents/mn-proxies")
OUTDIR <- "results/tables/protocol_v2"
fig <- readRDS(file.path(MNP, "results/res_full_bin_GW_Ghana_SL_all.rds"))

auc <- function(y, p) { ok <- is.finite(p) & is.finite(y); y <- y[ok]; p <- p[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
bss <- function(y, p) { p <- pmin(pmax(p, 0), 1); 1 - mean((y - p)^2) / (mean(y) * (1 - mean(y))) }
# the patterns name identifiers, sampling-design / region / team columns, dates, and blood-derived columns
FLAG <- paste0("cnum|^gw_cn$|bccn|b_cn$|hhln|bchnc|bclnr|bcgln|b_hhn|wlnr|^gw_in1$|^gw_wp1$|b_wp1|indivID|momID|pcn$|",
               "^gw_mcn$|Team|Region|region|Strata|strata|Number_of|probab|PSU|Check_|sWeight|Total.population|Households|",
               "INT_m|a_mon|VAS_m2|VAS_d$|VAS_y$|VAS_date|gw_month|gcmst|MalariaYN|wrmal|malref|sickle|Sickle|thal|Thal|mrdr|",
               "pctopendef|^Admin1$|^dataid$")

cells <- list(c("child_vitA", "Vit-A deficiency", "children"), c("women_vitA", "Vit-A deficiency", "women"),
              c("women_b12", "Vit-B12 deficiency", "women"), c("women_folate", "Folate deficiency", "women"),
              c("child_iron", "Ferratin deficiency", "children"), c("mom_iron", "Ferratin deficiency", "women"))
sets <- data.frame(suf = c("_gwPreds", "_full", ""), set = c("Survey", "All", "Proxy"))
rows <- list(); covs <- list()
for (cc in cells) for (j in seq_len(nrow(sets))) {
  f <- file.path(MNP, "results/models", paste0("res_full_bin_GW_Ghana_SL_", cc[1], sets$suf[j], ".rds"))
  if (!file.exists(f)) next
  r <- readRDS(f)
  y <- as.numeric(unclass(haven::zap_labels(r$Y))); cvp <- as.numeric(unclass(r$yhat_full))
  ins <- tryCatch(as.numeric(r$sl_fit$predict()), error = function(e) rep(NA_real_, length(y)))
  cl <- r$task$id; fo <- r$task$folds; fid <- integer(length(y)); for (k in seq_along(fo)) fid[fo[[k]]$validation_set] <- k
  cv <- r$task$nodes$covariates
  fg <- fig[fig$Y == cc[2] & fig$population == cc[3] & fig$Xvars == sets$set[j], ]
  rows[[length(rows) + 1]] <- data.frame(outcome = cc[1], set = sets$set[j], n = length(y), prevalence = mean(y),
    folds = length(fo), clusters = length(unique(cl)),
    clusters_split_across_folds = sum(tapply(fid, cl, function(v) length(unique(v)) > 1)),
    covariates = length(cv), covariates_flagged = sum(grepl(FLAG, cv)),
    figure_auc = fg$roc_auc, figure_brier_skill = fg$brier_skill,
    insample_auc = auc(y, ins), insample_brier_skill = bss(y, ins),
    cv_auc = auc(y, cvp), cv_brier_skill = bss(y, cvp))
  covs[[length(covs) + 1]] <- data.frame(outcome = cc[1], set = sets$set[j], covariate = cv, flagged = grepl(FLAG, cv))
  rm(r); gc(verbose = FALSE)
}
R <- do.call(rbind, rows)
write.csv(R, file.path(OUTDIR, "il02_january_figure_vs_cv.csv"), row.names = FALSE)
write.csv(do.call(rbind, covs), file.path(OUTDIR, "il02_january_covariates.csv"), row.names = FALSE)
options(width = 200)
print(transform(R, prevalence = round(prevalence, 3), figure_auc = round(figure_auc, 3), figure_brier_skill = round(figure_brier_skill, 3),
                insample_auc = round(insample_auc, 3), insample_brier_skill = round(insample_brier_skill, 3),
                cv_auc = round(cv_auc, 3), cv_brier_skill = round(cv_brier_skill, 3)), row.names = FALSE)
cat(sprintf("\nfigure equals in-sample in %d of %d fits (AUC and Brier skill to 3 decimals)\n",
            sum(abs(R$figure_auc - R$insample_auc) < 5e-4 & abs(R$figure_brier_skill - R$insample_brier_skill) < 5e-4, na.rm = TRUE), nrow(R)))
