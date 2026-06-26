#!/usr/bin/env Rscript
# Which domains carry transferable AREA-LEVEL signal?
# Single-domain lasso models, area-level LOCO, mean Spearman across held-out.
source(here::here("sandbox_fe","harness.R"))
dom_of <- function(cn){pref<-c(IHME="^ihme_",MAP="^MAP_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}
outcomes <- c("child_vitA","women_iron","child_iron")

aggregate_area <- function(o, cols, min_n=5){y<-o$y_bin;ok<-!is.na(y);a2<-o$admin2[ok];y<-y[ok];M<-o$M[ok,cols,drop=FALSE]
  ua<-unique(a2);keep<-ua[sapply(ua,function(u)sum(a2==u)>=min_n)]
  list(prev=sapply(keep,function(u)mean(y[a2==u])), X=t(sapply(keep,function(u){idx<-a2==u;colMeans(M[idx,,drop=FALSE],na.rm=TRUE)})))}

rows<-list();add<-function(...) rows[[length(rows)+1]]<<-data.frame(...)
for(otag in outcomes){
  objs<-load_objs(otag); shared<-Reduce(intersect,lapply(objs,function(o)colnames(o$M)))
  dm<-dom_of(shared); doms<-unique(dm[dm!="OTHER"])
  feature_sets<-c(setNames(lapply(doms,function(d)shared[dm==d]),doms), list(ALL=shared))
  for(fsn in names(feature_sets)){
    cols<-feature_sets[[fsn]]; if(length(cols)<2) next
    areas<-lapply(objs,aggregate_area,cols=cols)
    sps<-c()
    for(test_cn in names(areas)){
      tr<-setdiff(names(areas),test_cn)
      Xtr<-do.call(rbind,lapply(tr,function(c)areas[[c]]$X));ytr<-unlist(lapply(tr,function(c)areas[[c]]$prev))
      Xte<-areas[[test_cn]]$X;yte<-areas[[test_cn]]$prev; if(length(yte)<4) next
      cc<-.clean(Xtr,Xte);im<-imp_median(cc$tr,cc$te);z<-tf_zscore(im$tr,im$te)
      pr<-fit_predict_glmnet(as.matrix(z$tr),ytr,as.matrix(z$te),1,"gaussian",fixed_lambda=0.01)
      sps<-c(sps, suppressWarnings(cor(yte,pr,method="spearman",use="complete.obs")))
    }
    add(outcome=otag,domain=fsn,n_feat=length(cols),mean_spearman=round(mean(sps,na.rm=TRUE),3))
  }
}
res<-do.call(rbind,rows)
write.csv(res, here::here("sandbox_fe","results_12_domain.csv"),row.names=FALSE)
cat("=== Transferable area-level signal by domain (mean Spearman, LOCO) ===\n")
for(o in outcomes){s<-res[res$outcome==o,];s<-s[order(-s$mean_spearman),];cat("--",o,"--\n");print(s[,c("domain","n_feat","mean_spearman")],row.names=FALSE)}
