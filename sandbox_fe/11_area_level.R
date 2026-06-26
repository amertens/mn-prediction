#!/usr/bin/env Rscript
# Experiment: reframe as AREA-LEVEL prediction.
# Aggregate to Admin2 (prevalence + feature means), pool ~203 areas across
# countries. Compare individual-level vs area-level transfer (LOCO).
# Metric: Spearman cor(pred, observed area prevalence) in held-out country
#         = "do we correctly rank high- vs low-burden areas?"
source(here::here("sandbox_fe","harness.R"))

outcomes <- c("child_vitA","women_iron","child_iron")
dom_of <- function(cn){pref<-c(IHME="^ihme_",MAP="^MAP_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}

# aggregate one country's cache to admin2 level on the shared cols
aggregate_area <- function(o, shared, min_n=5) {
  y<-o$y_bin; ok<-!is.na(y); a2<-o$admin2[ok]; y<-y[ok]; M<-o$M[ok,shared,drop=FALSE]
  ua<-unique(a2); keep<-ua[sapply(ua,function(u)sum(a2==u)>=min_n)]
  prev<-sapply(keep,function(u)mean(y[a2==u]))
  nn  <-sapply(keep,function(u)sum(a2==u))
  Xa <- t(sapply(keep,function(u){ idx<-a2==u; colMeans(M[idx,,drop=FALSE],na.rm=TRUE) }))
  list(prev=prev, X=Xa, n=nn, area=keep)
}

rows<-list(); add<-function(...) rows[[length(rows)+1]]<<-data.frame(...)
for (otag in outcomes) {
  objs<-load_objs(otag)
  shared<-Reduce(intersect, lapply(objs,function(o)colnames(o$M)))
  areas<-lapply(objs, aggregate_area, shared=shared)
  na2<-sapply(areas,function(a)length(a$prev))
  cat(sprintf("\n%s areas/country: %s (total=%d)\n", otag, paste(sprintf("%s=%d",names(na2),na2),collapse=" "), sum(na2)))

  dommap<-setNames(dom_of(shared),shared)
  prep_pca<-function(Mtr,Mte){im<-imp_median(Mtr,Mte);Tr<-im$tr;Te<-im$te;doms<-dommap[colnames(Tr)];oT<-list();oE<-list()
    for(d in unique(doms[!is.na(doms)])){idx<-which(doms==d);if(length(idx)<2){oT[[d]]<-Tr[,idx,drop=FALSE];oE[[d]]<-Te[,idx,drop=FALSE];next}
      z<-tf_zscore(Tr[,idx,drop=FALSE],Te[,idx,drop=FALSE]);pc<-tryCatch(prcomp(z$tr,center=FALSE,scale.=FALSE),error=function(e)NULL)
      if(is.null(pc)){oT[[d]]<-z$tr;oE[[d]]<-z$te;next};kk<-min(3,ncol(pc$rotation))
      oT[[d]]<-z$tr%*%pc$rotation[,1:kk,drop=FALSE];oE[[d]]<-z$te%*%pc$rotation[,1:kk,drop=FALSE]}
    list(tr=do.call(cbind,oT),te=do.call(cbind,oE))}

  for (test_cn in names(areas)) {
    tr_cns<-setdiff(names(areas),test_cn)
    Xtr<-do.call(rbind,lapply(tr_cns,function(c)areas[[c]]$X)); ytr<-unlist(lapply(tr_cns,function(c)areas[[c]]$prev))
    Xte<-areas[[test_cn]]$X; yte<-areas[[test_cn]]$prev
    if(length(yte)<4) next
    cc<-.clean(Xtr,Xte)
    for(strat in c("full_z_lasso","pcaDom3_lasso","pcaDom3_ranger")){
      pp<-if(grepl("pca",strat)) prep_pca(cc$tr,cc$te) else { im<-imp_median(cc$tr,cc$te); tf_zscore(im$tr,im$te) }
      Tr<-as.matrix(pp$tr); Te<-as.matrix(pp$te)
      pr<-if(grepl("ranger",strat)) fit_predict_ranger(Tr,ytr,Te,"gaussian")
          else fit_predict_glmnet(Tr,ytr,Te,1,"gaussian",fixed_lambda=0.01)
      sp<-suppressWarnings(cor(yte,pr,method="spearman",use="complete.obs"))
      add(outcome=otag,test_country=test_cn,strategy=strat,n_area=length(yte),spearman=round(sp,3))
    }
  }
}
res<-do.call(rbind,rows)
write.csv(res, here::here("sandbox_fe","results_11_area.csv"), row.names=FALSE)
cat("\n=== AREA-LEVEL LOCO: Spearman cor(pred, obs area prevalence) ===\n")
print(res,row.names=FALSE)
cat("\n=== mean Spearman across held-out countries ===\n")
ag<-aggregate(spearman~outcome+strategy,res,function(x)round(mean(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag[ag$outcome==o,];s<-s[order(-s$spearman),];cat("--",o,"--\n");print(s[,c("strategy","spearman")],row.names=FALSE)}
