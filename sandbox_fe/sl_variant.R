# =============================================================================
# sl_variant.R — drive the PRODUCTION sl3 SuperLearner with swappable recipes.
# Mirrors DHS_SL_clustered() preprocessing exactly through impute+NZV, then
# applies one of several reduction/transform variants, trains the real SL with
# cluster-blocked CV, and returns cross-validated AUC (binary) / Spearman (cont).
#
# Variants:
#   current   : step_corr(0.85) + step_normalize (z-score)         [PRODUCTION]
#   rank      : step_corr(0.85) + rank/percentile transform
#   pca       : per-domain PCA (k PCs/domain) on z-scored features
#   rankpca   : per-domain PCA on rank-transformed features
#   bundle    : restrict Xvars to a supplied set, then current recipe
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(sl3); library(origami)
  library(recipes); library(caret); library(matrixStats)
})

dom_of_vec <- function(cn){
  pref <- c(DHS="^dhs", MICS="^mics_", IHME="^ihme_", LSMS="^lsms_", MAP="^MAP_",
            WFP="^wfp_", FLUNET="^flunet_", GEE="^gee_", FSEC="^fsec_", SOIL="^soil_",
            MISS="^missing_")
  out <- rep("OTHER", length(cn)); for (nm in names(pref)) out[grepl(pref[[nm]], cn)] <- nm
  out
}

.rank_transform <- function(M) {          # column-wise plotting-position ranks
  n <- nrow(M)
  apply(M, 2, function(x) (rank(x, ties.method="average") - 0.5)/n)
}
.per_domain_pca <- function(M, k=3) {     # z-score within domain, PCA, keep k PCs
  doms <- dom_of_vec(colnames(M)); out <- list()
  for (d in unique(doms)) {
    idx <- which(doms==d); Md <- M[, idx, drop=FALSE]
    mu <- colMeans(Md); sd <- matrixStats::colSds(Md); sd[sd<1e-8|!is.finite(sd)] <- 1
    Z <- sweep(sweep(Md,2,mu,"-"),2,sd,"/")
    if (ncol(Z) < 2) { colnames(Z) <- paste0(d,"_",seq_len(ncol(Z))); out[[d]] <- Z; next }
    pc <- tryCatch(prcomp(Z, center=FALSE, scale.=FALSE), error=function(e) NULL)
    if (is.null(pc)) { out[[d]] <- Z; next }
    kk <- min(k, ncol(pc$rotation)); P <- Z %*% pc$rotation[,1:kk,drop=FALSE]
    colnames(P) <- paste0(d,"_PC",1:kk); out[[d]] <- P
  }
  do.call(cbind, out)
}

#' Run one preprocessing variant through the real SL; return CV predictions+AUC.
run_sl_variant <- function(d, Xvars, outcome, id, sl, folds=5L,
                           variant="current", k_pca=3L, bundle=NULL, seed=12345L,
                           prescreen_pval=0.2) {
  if (!is.null(bundle)) Xvars <- intersect(bundle, Xvars)
  X   <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  cov <- labelled::unlabelled(X, user_na_to_na = TRUE)
  Y   <- d[[outcome]]; id_vec <- d[[id]]

  cov <- cov[, !sapply(cov, function(x) all(is.na(x))), drop=FALSE]
  cov <- cov %>% dplyr::select(dplyr::where(~{ nz<-.x[!is.na(.x)]; length(nz)>0 && length(unique(nz))>1 }))
  nzv <- caret::nearZeroVar(cov); if (length(nzv)>0) cov <- cov[,-nzv,drop=FALSE]
  cov <- as.data.frame(ck37r::impute_missing_values(cov, type="standard",
                          add_indicators=TRUE, prefix="missing_")$data)
  nzv <- caret::nearZeroVar(cov); if (length(nzv)>0) cov <- cov[,-nzv,drop=FALSE]

  is_bin <- length(unique(Y[!is.na(Y)]))==2
  fam <- if (is_bin) "binomial" else "gaussian"

  if (variant %in% c("pca","rankpca")) {
    M <- as.matrix(cov)
    if (variant=="rankpca") M <- .rank_transform(M)
    cov <- as.data.frame(.per_domain_pca(M, k=k_pca))
    # light prescreen on PCs
    Wv <- washb::washb_prescreen(Y=Y, Ws=cov, family=fam, pval=prescreen_pval, print=FALSE)
    if (length(Wv)>=2) cov <- cov[, Wv, drop=FALSE]
    rec <- recipes::recipe(~., data=cov) %>% recipes::step_zv(recipes::all_predictors()) %>%
           recipes::step_normalize(recipes::all_numeric()) %>% recipes::prep()
    cov <- data.frame(recipes::bake(rec, new_data=cov))
  } else {
    Wv <- washb::washb_prescreen(Y=Y, Ws=cov, family=fam, pval=prescreen_pval, print=FALSE)
    cov <- cov %>% dplyr::select(dplyr::all_of(Wv))
    rec <- recipes::recipe(~., data=cov) %>%
           recipes::step_zv(recipes::all_predictors()) %>%
           recipes::step_nzv(recipes::all_predictors()) %>%
           recipes::step_corr(recipes::all_numeric(), threshold=0.85) %>%
           recipes::prep()
    cov <- data.frame(recipes::bake(rec, new_data=cov))
    if (variant=="rank") {
      cov <- as.data.frame(.rank_transform(as.matrix(cov)))
    } else { # current = z-score
      mu <- colMeans(cov); sd <- matrixStats::colSds(as.matrix(cov)); sd[sd<1e-8|!is.finite(sd)] <- 1
      cov <- as.data.frame(sweep(sweep(cov,2,mu,"-"),2,sd,"/"))
    }
  }

  covars <- colnames(cov)
  dat <- data.table::data.table(Y=Y, id=id_vec, cov)
  set.seed(seed); fold_obj <- origami::make_folds(cluster_ids=id_vec, V=folds)
  task <- sl3::make_sl3_Task(data=dat, covariates=covars, outcome="Y", id="id", folds=fold_obj)
  suppressMessages(fit <- sl$train(task))
  yhat <- fit$predict_fold(task, "validation")

  auc <- NA_real_; sp <- NA_real_
  if (is_bin) {
    yy <- Y; ok <- !is.na(yy)&!is.na(yhat); n1<-sum(yy[ok]==1); n0<-sum(yy[ok]==0)
    if (n1>0&&n0>0){ r<-rank(yhat[ok]); auc<-(sum(r[yy[ok]==1])-n1*(n1+1)/2)/(n1*n0) }
  } else sp <- suppressWarnings(cor(Y, yhat, use="complete.obs", method="spearman"))
  list(auc=auc, spearman=sp, n=nrow(dat), p=length(covars), variant=variant)
}
