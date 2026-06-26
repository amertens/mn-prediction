# =============================================================================
# harness.R  — shared utilities for feature-engineering experiments (v2, fast)
#
# Speed design:
#   * Each cache object is numified ONCE (prep_obj) -> obj$M numeric matrix.
#   * CV / LOCO subset rows of M; preprocessing operates on numeric matrices.
#   * prep(Mtr, Mte, ytr) -> list(tr, te)  numeric matrices, train-stats only.
# =============================================================================
suppressMessages({
  library(here); library(glmnet); library(ranger); library(matrixStats); library(caret)
})

CACHE <- here::here("sandbox_fe", "cache")
COUNTRIES <- c("Gambia","Ghana","SierraLeone","Malawi")

numify <- function(X) {
  M <- suppressWarnings(matrix(as.numeric(as.matrix(X)), nrow = nrow(X)))
  colnames(M) <- colnames(X); M
}

# Attach numeric matrix once; cache to disk to avoid recomputation.
prep_obj <- function(obj) {
  if (!is.null(obj$M)) return(obj)
  obj$M <- numify(obj$X); obj$X <- NULL; obj
}

load_objs <- function(outcome, countries = COUNTRIES, numeric_cache = TRUE) {
  objs <- list()
  for (cn in countries) {
    fn  <- file.path(CACHE, sprintf("%s__%s.rds", cn, outcome))
    nfn <- file.path(CACHE, sprintf("%s__%s__num.rds", cn, outcome))
    if (numeric_cache && file.exists(nfn)) { objs[[cn]] <- readRDS(nfn); next }
    if (!file.exists(fn)) next
    o <- prep_obj(readRDS(fn))
    if (numeric_cache) saveRDS(o, nfn)
    objs[[cn]] <- o
  }
  objs
}

fast_auc <- function(y, p) {
  y <- as.numeric(y); ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y==1); n0 <- sum(y==0); if (n1==0||n0==0) return(NA_real_)
  r <- rank(p); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)
}

# ---- column hygiene on TRAIN matrix -----------------------------------------
.clean <- function(Mtr, Mte) {
  keep <- colSums(!is.na(Mtr)) >= 5
  Mtr <- Mtr[,keep,drop=FALSE]; Mte <- Mte[,keep,drop=FALSE]
  sds <- matrixStats::colSds(Mtr, na.rm=TRUE)
  keep2 <- is.finite(sds) & sds > 1e-8
  list(tr=Mtr[,keep2,drop=FALSE], te=Mte[,keep2,drop=FALSE])
}

# ---- imputation (train stats only) ------------------------------------------
imp_median <- function(Mtr, Mte, add_ind=FALSE) {
  med <- matrixStats::colMedians(Mtr, na.rm=TRUE); med[!is.finite(med)] <- 0
  fill <- function(M){ idx<-which(is.na(M),arr.ind=TRUE); if(nrow(idx)) M[idx]<-med[idx[,2]]; M }
  Tr<-fill(Mtr); Te<-fill(Mte)
  if (add_ind) {
    frac <- colMeans(is.na(Mtr)); jj <- which(frac>0.02 & frac<0.98)
    if (length(jj)) {
      trI <- (is.na(Mtr)[,jj,drop=FALSE])*1.0; teI <- (is.na(Mte)[,jj,drop=FALSE])*1.0
      colnames(trI)<-paste0("miss_",colnames(Mtr)[jj]); colnames(teI)<-colnames(trI)
      Tr<-cbind(Tr,trI); Te<-cbind(Te,teI)
    }
  }
  list(tr=Tr, te=Te)
}
imp_mean <- function(Mtr, Mte) {
  mu <- colMeans(Mtr, na.rm=TRUE); mu[!is.finite(mu)] <- 0
  fill <- function(M){ idx<-which(is.na(M),arr.ind=TRUE); if(nrow(idx)) M[idx]<-mu[idx[,2]]; M }
  list(tr=fill(Mtr), te=fill(Mte))
}
imp_zero <- function(Mtr, Mte) {  # impute 0 (after zscore=mean); proxy for "no info"
  list(tr=replace(Mtr,is.na(Mtr),0), te=replace(Mte,is.na(Mte),0))
}
# vectorized KNN on standardized, median-filled space (donors = train rows)
imp_knn <- function(Mtr, Mte, k=10) {
  base <- imp_median(Mtr, Mte)
  mu <- colMeans(base$tr); sd <- matrixStats::colSds(base$tr); sd[sd<1e-8|!is.finite(sd)]<-1
  Ztr <- sweep(sweep(base$tr,2,mu,"-"),2,sd,"/")
  fillM <- function(Morig, Bfilled) {
    Z <- sweep(sweep(Bfilled,2,mu,"-"),2,sd,"/")
    res <- Bfilled
    narows <- which(rowSums(is.na(Morig))>0)
    if (length(narows)) {
      D <- proxy_dist(Z[narows,,drop=FALSE], Ztr)  # rows x ntrain
      for (ii in seq_along(narows)) {
        i <- narows[ii]; nc <- which(is.na(Morig[i,]))
        nn <- order(D[ii,])[seq_len(min(k,ncol(D)))]
        res[i,nc] <- colMeans(base$tr[nn,nc,drop=FALSE])
      }
    }
    res
  }
  list(tr=fillM(Mtr, base$tr), te=fillM(Mte, base$te))
}
proxy_dist <- function(A, B) {  # squared euclidean rows(A) x rows(B), fast
  aa <- rowSums(A^2); bb <- rowSums(B^2)
  outer(aa, bb, "+") - 2*(A %*% t(B))
}

# ---- transforms (train stats only) ------------------------------------------
tf_none   <- function(Tr,Te) list(tr=Tr, te=Te)
tf_zscore <- function(Tr,Te){ mu<-colMeans(Tr);sd<-matrixStats::colSds(Tr);sd[sd<1e-8|!is.finite(sd)]<-1
  list(tr=sweep(sweep(Tr,2,mu,"-"),2,sd,"/"), te=sweep(sweep(Te,2,mu,"-"),2,sd,"/")) }
tf_signedlog <- function(Tr,Te){ g<-function(M) sign(M)*log1p(abs(M)); tf_zscore(g(Tr),g(Te)) }
tf_rank <- function(Tr,Te){       # train->plotting position; test via train ECDF (vectorized)
  n<-nrow(Tr); Trr<-Tr; Ter<-Te
  for(j in seq_len(ncol(Tr))){ x<-Tr[,j]
    Trr[,j]<-(rank(x,ties.method="average")-0.5)/n
    Ter[,j]<-findInterval(Te[,j], sort(x))/n }
  list(tr=Trr, te=Ter)
}

# ---- correlation filter (fast, caret C) -------------------------------------
corr_filter <- function(Tr, Te, thresh=0.85) {
  if (ncol(Tr) < 2) return(list(tr=Tr,te=Te))
  C <- suppressWarnings(cor(Tr)); C[!is.finite(C)] <- 0
  drop <- caret::findCorrelation(C, cutoff=thresh)
  if (length(drop)) { Tr<-Tr[,-drop,drop=FALSE]; Te<-Te[,-drop,drop=FALSE] }
  list(tr=Tr, te=Te)
}

# ---- models -----------------------------------------------------------------
fit_predict_glmnet <- function(Xtr, ytr, Xte, alpha=1, family="binomial", fixed_lambda=NULL) {
  ok <- complete.cases(Xtr) & !is.na(ytr)
  if (!is.null(fixed_lambda)) {
    fit <- tryCatch(glmnet(Xtr[ok,,drop=FALSE], ytr[ok], alpha=alpha, family=family), error=function(e) NULL)
    if (is.null(fit)) return(rep(mean(ytr,na.rm=TRUE), nrow(Xte)))
    return(as.numeric(predict(fit, newx=Xte, s=fixed_lambda,
                              type=if(family=="binomial")"response" else "link")))
  }
  cv <- tryCatch(cv.glmnet(Xtr[ok,,drop=FALSE], ytr[ok], alpha=alpha, family=family, nfolds=5),
                 error=function(e) NULL)
  if (is.null(cv)) return(rep(mean(ytr,na.rm=TRUE), nrow(Xte)))
  as.numeric(predict(cv, newx=Xte, s="lambda.min",
                     type=if(family=="binomial")"response" else "link"))
}
fit_predict_ranger <- function(Xtr, ytr, Xte, family="binomial", num.trees=300) {
  ok <- !is.na(ytr); prob <- family=="binomial"
  yy <- if (prob) factor(ytr[ok], levels=c(0,1)) else ytr[ok]
  rf <- ranger::ranger(x=Xtr[ok,,drop=FALSE], y=yy, num.trees=num.trees,
                       probability=prob, min.node.size=ifelse(prob,10,5), num.threads=0)
  pr <- predict(rf, data=Xte)$predictions
  if (prob) pr[,"1"] else pr
}

# ---- within-country cluster-blocked CV --------------------------------------
eval_within <- function(obj, prep, target=c("binary","continuous"),
                        model=c("glmnet","ranger"), alpha=1, K=5, seed=1,
                        cols=NULL, fixed_lambda=NULL) {
  target<-match.arg(target); model<-match.arg(model)
  y <- if(target=="binary") obj$y_bin else obj$y_cont
  ok <- !is.na(y)
  M <- obj$M; if (!is.null(cols)) M <- M[, intersect(cols, colnames(M)), drop=FALSE]
  M <- M[ok,,drop=FALSE]; y <- y[ok]; cl <- obj$cluster[ok]
  fam <- if(target=="binary") "binomial" else "gaussian"
  set.seed(seed); ucl<-unique(cl); fo<-setNames(sample(rep(seq_len(K),length.out=length(ucl))),ucl)
  folds <- fo[cl]; preds <- rep(NA_real_, length(y))
  for (k in seq_len(K)) {
    tr <- folds!=k; te <- !tr
    if (sum(te)==0 || (target=="binary" && length(unique(y[tr]))<2)) next
    cc <- .clean(M[tr,,drop=FALSE], M[te,,drop=FALSE])
    pp <- prep(cc$tr, cc$te, y[tr])
    pr <- if (model=="glmnet") fit_predict_glmnet(pp$tr,y[tr],pp$te,alpha,fam,fixed_lambda)
          else fit_predict_ranger(pp$tr,y[tr],pp$te,fam)
    preds[te] <- pr
  }
  if (target=="binary") fast_auc(y, preds)
  else suppressWarnings(cor(y, preds, use="complete.obs", method="spearman"))
}

# ---- leave-one-country-out transfer -----------------------------------------
eval_loco <- function(objs, prep, shared, target=c("binary","continuous"),
                      model=c("glmnet","ranger"), alpha=1, per_country_norm=FALSE,
                      fixed_lambda=NULL) {
  target<-match.arg(target); model<-match.arg(model)
  fam <- if(target=="binary") "binomial" else "gaussian"
  norm_one <- function(M){ for(j in seq_len(ncol(M))){x<-M[,j]; nn<-sum(!is.na(x)); if(nn>1) M[,j]<-(rank(x,na.last="keep",ties.method="average")-0.5)/nn}; M }
  res <- list()
  for (test_cn in names(objs)) {
    train_cns <- setdiff(names(objs), test_cn)
    getXY <- function(cn){ o<-objs[[cn]]; y<-if(target=="binary")o$y_bin else o$y_cont; ok<-!is.na(y)
      list(X=o$M[ok,shared,drop=FALSE], y=y[ok]) }
    trp <- lapply(train_cns, getXY); tep <- getXY(test_cn)
    if (target=="binary" && length(unique(tep$y))<2) next
    if (per_country_norm){ trp<-lapply(trp,function(p){p$X<-norm_one(p$X);p}); tep$X<-norm_one(tep$X) }
    Xtr <- do.call(rbind, lapply(trp,function(p)p$X)); ytr <- unlist(lapply(trp,function(p)p$y))
    cc <- .clean(Xtr, tep$X)
    pp <- prep(cc$tr, cc$te, ytr)
    pr <- if (model=="glmnet") fit_predict_glmnet(pp$tr,ytr,pp$te,alpha,fam,fixed_lambda)
          else fit_predict_ranger(pp$tr,ytr,pp$te,fam)
    res[[test_cn]] <- if(target=="binary") fast_auc(tep$y,pr)
                      else suppressWarnings(cor(tep$y,pr,use="complete.obs",method="spearman"))
  }
  unlist(res)
}

# Convenience prep builder operating on numeric matrices
build_prep <- function(impute="median", add_ind=FALSE, transform="zscore",
                       reduce="none", corr_thresh=0.85) {
  function(Mtr, Mte, ytr) {
    im <- switch(impute, median=imp_median(Mtr,Mte,add_ind), mean=imp_mean(Mtr,Mte),
                 zero=imp_zero(Mtr,Mte), knn=imp_knn(Mtr,Mte,10))
    Tr<-im$tr; Te<-im$te
    if (reduce=="corr"){ r<-corr_filter(Tr,Te,corr_thresh); Tr<-r$tr; Te<-r$te }
    tt <- switch(transform, none=tf_none(Tr,Te), zscore=tf_zscore(Tr,Te),
                 rank=tf_rank(Tr,Te), signedlog=tf_signedlog(Tr,Te))
    list(tr=as.matrix(tt$tr), te=as.matrix(tt$te))
  }
}
