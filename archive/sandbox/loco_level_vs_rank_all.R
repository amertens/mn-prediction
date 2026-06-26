# =============================================================================
# sandbox/loco_level_vs_rank_all.R
#
# Extend the Gambia/child-iron re-harmonization finding to ALL held-out
# countries and outcomes. Tests whether the pattern generalises:
#   transport fails on the LEVEL (bias/calibration) offset, not on RANKING.
#
# For every outcome and every held-out country: train an elastic-net area model
# on the OTHER countries' admin-2 (svy_prev ~ common GEE covariates), predict the
# held-out country's admin-2 prevalence, and decompose the error:
#   r            — Pearson corr of predicted vs observed admin-2 prevalence
#                  (scale/shift-invariant -> pure RANKING transfer)
#   bias_pp      — mean(pred - truth) -> the cross-survey LEVEL offset
#   mae_pp       — total absolute error
#   mae_ctr_pp   — MAE after removing the constant offset from BOTH series
#                  (= ranking+spread error only). If mae_ctr << mae, the failure
#                  is a level offset, not a pattern failure.
# =============================================================================
suppressWarnings(suppressMessages({
  library(here); library(dplyr); library(targets); library(glmnet)
}))
STORE <- "_targets_full"
cc_keys  <- c(gambia="Gambia", ghana="Ghana", malawi="Malawi", sierraleone="Sierra Leone")
outcomes <- c("child_vitA","women_vitA","child_iron","women_iron","women_folate","women_b12")
meta_nm  <- tar_meta(store = STORE)$name

# cache GEE covariates per country
gee <- list()
for (k in names(cc_keys)) {
  t <- paste0("gee_admin2_", k)
  if (t %in% meta_nm) gee[[k]] <- tryCatch(tar_read_raw(t, store = STORE), error = function(e) NULL)
}

fit_fold <- function(train, test, vars) {
  Xtr <- as.matrix(train[, vars, drop = FALSE]); Xte <- as.matrix(test[, vars, drop = FALSE])
  med <- apply(Xtr, 2, function(z){ m <- median(z[is.finite(z)], na.rm=TRUE); if (is.finite(m)) m else 0 })
  for (j in seq_along(vars)) { zt<-Xtr[,j]; zt[!is.finite(zt)]<-med[j]; Xtr[,j]<-zt
                               ze<-Xte[,j]; ze[!is.finite(ze)]<-med[j]; Xte[,j]<-ze }
  keep <- which(apply(Xtr, 2, function(z) stats::sd(z) > 0))
  if (length(keep) < 2) return(rep(NA_real_, nrow(test)))
  cv <- tryCatch(cv.glmnet(Xtr[,keep,drop=FALSE], train$svy_prev, alpha = 0.5,
                           weights = pmax(train$n_svy, 1)), error = function(e) NULL)
  if (is.null(cv)) return(rep(NA_real_, nrow(test)))
  pmin(pmax(as.numeric(predict(cv, Xte[,keep,drop=FALSE], s = "lambda.min")), 0), 1)
}
metrics <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred); truth<-truth[ok]; pred<-pred[ok]
  if (length(truth) < 4) return(NULL)
  tc <- truth - mean(truth); pc <- pred - mean(pred)
  data.frame(n = length(truth),
             r          = round(suppressWarnings(cor(truth, pred)), 3),
             truth_pp   = round(mean(truth)*100, 1),
             pred_pp    = round(mean(pred)*100, 1),
             bias_pp    = round(mean(pred - truth)*100, 1),
             mae_pp     = round(mean(abs(pred - truth))*100, 1),
             mae_ctr_pp = round(mean(abs(pc - tc))*100, 1))
}

res <- list()
for (oc in outcomes) {
  avail <- names(cc_keys)[vapply(names(cc_keys),
            function(k) paste0("svy_admin2_",k,"_",oc) %in% meta_nm && !is.null(gee[[k]]), logical(1))]
  if (length(avail) < 2) next
  svy <- setNames(lapply(avail, function(k)
            tryCatch(tar_read_raw(paste0("svy_admin2_",k,"_",oc), store=STORE), error=function(e) NULL)), avail)
  avail <- avail[!vapply(svy, is.null, logical(1))]; svy <- svy[avail]
  common <- Reduce(intersect, lapply(avail, function(k) grep("^gee_", names(gee[[k]]), value=TRUE)))
  if (length(common) < 5) next
  frame <- function(k) dplyr::inner_join(
      svy[[k]][, c("Admin2","svy_prev","n_svy")], gee[[k]][, c("Admin2", common)], "Admin2") |>
      dplyr::filter(is.finite(svy_prev))
  for (h in avail) {
    tr <- dplyr::bind_rows(lapply(setdiff(avail, h), frame))
    te <- frame(h)
    if (nrow(tr) < 10 || nrow(te) < 4) next
    m <- metrics(te$svy_prev, fit_fold(tr, te, common))
    if (!is.null(m)) res[[length(res)+1]] <- cbind(outcome=oc, held_out=cc_keys[[h]], m)
  }
}
out <- do.call(rbind, res)
options(width = 200)
for (oc in unique(out$outcome)) {
  cat("\n=====", oc, "=====\n")
  print(out[out$outcome==oc, c("held_out","n","r","truth_pp","pred_pp","bias_pp","mae_pp","mae_ctr_pp")], row.names = FALSE)
}
cat("\n--- summary across all folds ---\n")
cat(sprintf("median r = %.2f | median |bias| = %.1f pp | median MAE = %.1f pp | median centered MAE = %.1f pp\n",
    median(out$r, na.rm=TRUE), median(abs(out$bias_pp), na.rm=TRUE),
    median(out$mae_pp, na.rm=TRUE), median(out$mae_ctr_pp, na.rm=TRUE)))
cat(sprintf("MAE reduction from removing level offset: %.0f%% (median)\n",
    100*(1 - median(out$mae_ctr_pp,na.rm=TRUE)/median(out$mae_pp,na.rm=TRUE))))
saveRDS(out, here::here("results","prototype_sae_sl","loco_level_vs_rank.rds"))
