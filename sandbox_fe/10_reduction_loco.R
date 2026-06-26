#!/usr/bin/env Rscript
# Experiment: does dimension reduction improve LOCO transfer?
# Combine with per-country normalization (the two transfer-stabilizing ideas).
source(here::here("sandbox_fe","harness.R"))

outcomes <- c("child_vitA","women_iron","child_iron")
dom_of <- function(cn){pref<-c(IHME="^ihme_",MAP="^MAP_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}

prep_pca_domain <- function(dommap,k=3) function(Mtr,Mte,ytr){
  im<-imp_median(Mtr,Mte);Tr<-im$tr;Te<-im$te;doms<-dommap[colnames(Tr)];oT<-list();oE<-list()
  for(d in unique(doms[!is.na(doms)])){idx<-which(doms==d)
    if(length(idx)<2){oT[[d]]<-Tr[,idx,drop=FALSE];oE[[d]]<-Te[,idx,drop=FALSE];next}
    z<-tf_zscore(Tr[,idx,drop=FALSE],Te[,idx,drop=FALSE]);pc<-tryCatch(prcomp(z$tr,center=FALSE,scale.=FALSE),error=function(e)NULL)
    if(is.null(pc)){oT[[d]]<-z$tr;oE[[d]]<-z$te;next};kk<-min(k,ncol(pc$rotation))
    oT[[d]]<-z$tr%*%pc$rotation[,1:kk,drop=FALSE];oE[[d]]<-z$te%*%pc$rotation[,1:kk,drop=FALSE]}
  list(tr=do.call(cbind,oT),te=do.call(cbind,oE))
}

rows<-list(); add<-function(...) rows[[length(rows)+1]]<<-data.frame(...)
for (otag in outcomes) {
  objs<-load_objs(otag)
  shared<-Reduce(intersect, lapply(objs,function(o)colnames(o$M)))
  dommap<-setNames(dom_of(shared),shared)
  full_z   <- build_prep("median",FALSE,"zscore","none")
  pca3     <- prep_pca_domain(dommap,3)
  configs<-list(
    list(tag="full lasso pooled",     model="glmnet", prep=full_z, pcn=FALSE, lam=0.02),
    list(tag="full lasso pcNorm",     model="glmnet", prep=build_prep("median",FALSE,"none","none"), pcn=TRUE, lam=0.02),
    list(tag="pcaDom3 lasso pooled",  model="glmnet", prep=pca3, pcn=FALSE, lam=0.02),
    list(tag="pcaDom3 lasso pcNorm",  model="glmnet", prep=pca3, pcn=TRUE, lam=0.02),
    list(tag="full ranger pooled",    model="ranger", prep=build_prep("median",FALSE,"none","none"), pcn=FALSE),
    list(tag="pcaDom3 ranger pooled", model="ranger", prep=pca3, pcn=FALSE),
    list(tag="pcaDom3 ranger pcNorm", model="ranger", prep=pca3, pcn=TRUE)
  )
  for(cfg in configs){
    m<-tryCatch(eval_loco(objs,cfg$prep,shared,"binary",cfg$model,per_country_norm=cfg$pcn,
                          fixed_lambda=if(!is.null(cfg$lam))cfg$lam else NULL),
                error=function(e){cat("ERR",otag,cfg$tag,conditionMessage(e),"\n");NULL})
    if(!is.null(m)) for(cn in names(m)) add(outcome=otag,config=cfg$tag,test_country=cn,auc=round(m[[cn]],3))
  }
  cat("done",otag,"\n")
}
res<-do.call(rbind,rows)
write.csv(res, here::here("sandbox_fe","results_10_reduction_loco.csv"), row.names=FALSE)
cat("\n=== LOCO transfer mean AUC (across held-out countries) ===\n")
ag<-aggregate(auc~outcome+config,res,function(x)round(mean(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag[ag$outcome==o,];s<-s[order(-s$auc),];cat("--",o,"--\n");print(s[,c("config","auc")],row.names=FALSE)}
# stability: min across held-out countries (worst-case transfer)
cat("\n=== WORST-CASE (min) AUC across held-out countries ===\n")
ag2<-aggregate(auc~outcome+config,res,function(x)round(min(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag2[ag2$outcome==o,];s<-s[order(-s$auc),];cat("--",o,"--\n");print(s[,c("config","auc")],row.names=FALSE)}
