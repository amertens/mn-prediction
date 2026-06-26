#!/usr/bin/env Rscript
# Experiment: does aggressive, leak-free dimension reduction match/beat the
# full ~1000-feature set, given only 14-87 effective area-points?
#  Strategies: full | per-domain PCA(k) | cluster-resolution-only | global PCA
source(here::here("sandbox_fe","harness.R"))

outcomes <- c("child_vitA","women_iron","child_iron")
dom_of <- function(cn){pref<-c(DHS="^dhs",MICS="^mics_",IHME="^ihme_",LSMS="^lsms_",MAP="^MAP_",WFP="^wfp_",FLUNET="^flunet_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}

# --- preps -------------------------------------------------------------------
prep_full <- build_prep("median",FALSE,"zscore","none")

prep_pca_domain <- function(dommap, k=3) function(Mtr,Mte,ytr){
  im<-imp_median(Mtr,Mte); Tr<-im$tr; Te<-im$te
  doms<-dommap[colnames(Tr)]; outTr<-list(); outTe<-list()
  for(d in unique(doms[!is.na(doms)])){
    idx<-which(doms==d)
    if(length(idx)<2){outTr[[d]]<-Tr[,idx,drop=FALSE];outTe[[d]]<-Te[,idx,drop=FALSE];next}
    z<-tf_zscore(Tr[,idx,drop=FALSE],Te[,idx,drop=FALSE])
    pc<-tryCatch(prcomp(z$tr,center=FALSE,scale.=FALSE),error=function(e)NULL)
    if(is.null(pc)){outTr[[d]]<-z$tr;outTe[[d]]<-z$te;next}
    kk<-min(k,ncol(pc$rotation))
    a<-z$tr%*%pc$rotation[,1:kk,drop=FALSE]; b<-z$te%*%pc$rotation[,1:kk,drop=FALSE]
    colnames(a)<-paste0(d,"_PC",1:kk);colnames(b)<-colnames(a)
    outTr[[d]]<-a; outTe[[d]]<-b
  }
  list(tr=do.call(cbind,outTr), te=do.call(cbind,outTe))
}
prep_pca_global <- function(k=10) function(Mtr,Mte,ytr){
  im<-imp_median(Mtr,Mte); z<-tf_zscore(im$tr,im$te)
  pc<-prcomp(z$tr,center=FALSE,scale.=FALSE); kk<-min(k,ncol(pc$rotation))
  list(tr=z$tr%*%pc$rotation[,1:kk,drop=FALSE], te=z$te%*%pc$rotation[,1:kk,drop=FALSE])
}

rows<-list()
add<-function(...) rows[[length(rows)+1]]<<-data.frame(...)
for (otag in outcomes) {
  objs <- load_objs(otag)
  for (cn in names(objs)) {
    o<-objs[[cn]]; dommap<-setNames(dom_of(colnames(o$M)), colnames(o$M))
    # cluster-resolution-only column set
    ndist<-apply(o$M,2,function(x)length(unique(round(x[!is.na(x)],8))))
    na2<-length(unique(o$admin2)); clcols<-colnames(o$M)[ndist>na2+1]
    strategies <- list(
      full        = list(model="ranger", prep=prep_full, cols=NULL),
      full_lasso  = list(model="glmnet", prep=build_prep("median",FALSE,"rank","none"), cols=NULL, lam=0.02),
      pca_dom3    = list(model="ranger", prep=prep_pca_domain(dommap,3), cols=NULL),
      pca_dom3_ls = list(model="glmnet", prep=prep_pca_domain(dommap,3), cols=NULL, lam=0.02),
      pca_glob10  = list(model="ranger", prep=prep_pca_global(10), cols=NULL),
      clusterres  = list(model="ranger", prep=prep_full, cols=clcols)
    )
    for (sn in names(strategies)) {
      st<-strategies[[sn]]
      a<-tryCatch(eval_within(o, st$prep, "binary", st$model, alpha=1, cols=st$cols,
                              fixed_lambda=if(!is.null(st$lam)) st$lam else NULL),
                  error=function(e){cat("ERR",otag,cn,sn,conditionMessage(e),"\n");NA})
      add(outcome=otag, country=cn, strategy=sn, auc=round(a,3))
    }
  }
  cat("done",otag,"\n")
}
res<-do.call(rbind,rows)
write.csv(res, here::here("sandbox_fe","results_09_reduction.csv"), row.names=FALSE)
cat("\n=== mean AUC across countries (within-country CV) ===\n")
ag<-aggregate(auc~outcome+strategy,res,function(x)round(mean(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag[ag$outcome==o,];s<-s[order(-s$auc),];cat("--",o,"--\n");print(s[,c("strategy","auc")],row.names=FALSE)}
