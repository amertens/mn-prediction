#!/usr/bin/env Rscript
# ===========================================================================
# 11_distributional_heterosked_transport.R
#
# Extends 10_distributional_vs_binary_prototype.R with two questions:
#   (1a) Does a HETEROSKEDASTIC-sigma distributional model (log-variance modelled
#        on covariates) widen the tail advantage over the homoskedastic one?
#   (1b) Does modelling the continuous distribution transport RANKING better than
#        dichotomising, across countries (leave-one-country-out), even though the
#        biomarker LEVEL offset is known to break transport?
#
# Strategies (identical covariates + honest CV throughout):
#   S1b  individual logistic on binary deficiency
#   S2   individual gaussian on log-biomarker, HOMOskedastic empirical-residual integration
#   S3   individual gaussian on log-biomarker + log-variance model, HETEROskedastic
#        standardized-residual integration
#
# Run from repo root.  Reads _targets_full; writes:
#   results/sensitivity/distributional_within_heterosked.csv
#   results/sensitivity/distributional_transport_childiron.csv
# Prototype only: small fixed covariate set, ridge lambda is stochastic (~+/-0.03).
# ===========================================================================
suppressMessages({library(targets); library(glmnet)})
tar_config_set(store = "_targets_full")
dir.create("results/sensitivity", showWarnings = FALSE, recursive = TRUE)

prep <- function(Xtr, Xte){                 # train-only standardize + median impute
  med<-apply(Xtr,2,median,na.rm=TRUE); med[!is.finite(med)]<-0
  ctr<-colMeans(Xtr,na.rm=TRUE); ctr[!is.finite(ctr)]<-0
  sdv<-apply(Xtr,2,sd,na.rm=TRUE); sdv[!is.finite(sdv)|sdv==0]<-1
  imp<-function(M){for(j in seq_len(ncol(M))){x<-M[,j];x[is.na(x)]<-med[j];M[,j]<-(x-ctr[j])/sdv[j]};M}
  list(tr=imp(Xtr), te=imp(Xte))
}
rfit <- function(x,y,fam,w=NULL){ if(is.null(w)) w<-rep(1,length(y))
  tryCatch(cv.glmnet(x,y,family=fam,alpha=0,weights=w,nfolds=5), error=function(e) NULL) }

# S1b / S2 / S3 predictions for one prepped train/test split
predict_strats <- function(Ptr, Pte, def_tr, logc_tr, w_tr, thr){
  out<-list(S1b=NULL,S2=NULL,S3=NULL)
  f1b<-rfit(Ptr,def_tr,"binomial",w_tr); f2<-rfit(Ptr,logc_tr,"gaussian",w_tr)
  if(is.null(f2)) return(out)
  mu_tr<-as.numeric(predict(f2,Ptr,s="lambda.min")); e<-logc_tr-mu_tr
  mu_te<-as.numeric(predict(f2,Pte,s="lambda.min"))
  if(!is.null(f1b)) out$S1b<-as.numeric(predict(f1b,Pte,s="lambda.min",type="response"))
  out$S2<-sapply(mu_te,function(m) weighted.mean((m+e)<thr,w_tr))
  fv<-rfit(Ptr,log(e^2+1e-6),"gaussian",w_tr)
  if(!is.null(fv)){
    sig_tr<-exp(0.5*as.numeric(predict(fv,Ptr,s="lambda.min"))); sig_tr[!is.finite(sig_tr)|sig_tr<=0]<-sd(e)
    sig_te<-exp(0.5*as.numeric(predict(fv,Pte,s="lambda.min"))); sig_te[!is.finite(sig_te)|sig_te<=0]<-sd(e)
    rstd<-e/sig_tr
    out$S3<-mapply(function(m,s) weighted.mean((m+s*rstd)<thr,w_tr), mu_te, sig_te)
  } else out$S3<-out$S2
  out
}

## ---- PART 1: within-country (Ghana), heteroskedastic S3 ------------------
mg_all<-tar_read_raw("merged_ghana"); acg<-tar_read_raw("area_covariates_ghana")$gee_admin2
pick<-function(rx){h<-grep(rx,colnames(acg),value=TRUE); if(length(h))h[1] else NA}
covg<-na.omit(c(pick("gee_soilzinc_mean_0_20"),pick("gee_soiliron_mean_0_20"),
  pick("gee_soilnitrogen_mean_0_20"),pick("gee_soilphosphorus_mean_0_20"),
  pick("gee_accessibility_2019"),pick("gee_elevation_2000"),pick("gee_popdensity_2010"),
  pick("gee_ghsbuilts_2015_built_surface")))
acg2<-acg[,c("Admin2",covg)]; acg2$Admin2<-as.character(acg2$Admin2)
for(cc in covg) acg2[[cc]]<-suppressWarnings(as.numeric(acg2[[cc]]))

run_within<-function(cont_col,thr,child,wt_col,svy){
  set.seed(20260630)
  mg<-mg_all; mg<-if(child) mg[which(mg$gw_child_flag==1),] else mg[which(mg$gw_child_flag==0),]
  ind<-data.frame(Admin2=as.character(mg$Admin2),cont=suppressWarnings(as.numeric(mg[[cont_col]])),
                  wt=suppressWarnings(as.numeric(mg[[wt_col]])),stringsAsFactors=FALSE)
  ind$wt[is.na(ind$wt)]<-1; ind<-ind[is.finite(ind$cont)&ind$cont>0&!is.na(ind$Admin2)&nzchar(ind$Admin2),]
  ind$logc<-log(ind$cont); ind$def<-as.integer(ind$cont<thr); lthr<-log(thr)
  sv<-as.data.frame(tar_read_raw(svy)); sv$Admin2<-as.character(sv$Admin2)
  areas<-sort(unique(ind$Admin2)); areas<-areas[areas%in%acg2$Admin2 & areas%in%sv$Admin2]
  A<-acg2[match(areas,acg2$Admin2),]; sv<-sv[match(areas,sv$Admin2),]; ind<-ind[ind$Admin2%in%areas,]
  X<-as.matrix(A[,covg,drop=FALSE]); nA<-length(areas)
  pr<-data.frame(Admin2=areas,truth=sv$svy_prev,n_svy=sv$n_svy,S1b=NA,S2=NA,S3=NA)
  for(i in seq_len(nA)){
    itr<-ind[ind$Admin2%in%areas[-i],]; ite<-ind[ind$Admin2==areas[i],]
    P<-prep(X[match(itr$Admin2,areas),,drop=FALSE],X[match(ite$Admin2,areas),,drop=FALSE])
    o<-predict_strats(P$tr,P$te,itr$def,itr$logc,itr$wt,lthr)
    if(!is.null(o$S1b))pr$S1b[i]<-mean(o$S1b); if(!is.null(o$S2))pr$S2[i]<-mean(o$S2); if(!is.null(o$S3))pr$S3[i]<-mean(o$S3)
  }
  ev<-pr[pr$n_svy>=8 & is.finite(pr$truth),]
  m<-function(p)c(pearson=cor(p,ev$truth),spearman=cor(p,ev$truth,method="spearman"),mae=mean(abs(p-ev$truth)))
  list(natprev=weighted.mean(ind$def,ind$wt),n=nrow(ev),tab=rbind(S1b=m(ev$S1b),S2=m(ev$S2),S3=m(ev$S3)))
}
cells<-list(
  child_vitA=list(c="gw_cRBPAdjThurn",t=0.70,ch=TRUE, w="gw_child_weight",s="svy_admin2_ghana_child_vitA"),
  child_iron=list(c="gw_cFerrAdjThurn",t=12, ch=TRUE, w="gw_child_weight",s="svy_admin2_ghana_child_iron"),
  women_iron=list(c="gw_wFerrAdjThurn",t=15, ch=FALSE,w="gw_sWeight",     s="svy_admin2_ghana_women_iron"))
rows1<-list()
cat("######## PART 1: within-country (Ghana), heteroskedastic S3 ########\n")
for(nm in names(cells)){cl<-cells[[nm]]; r<-run_within(cl$c,cl$t,cl$ch,cl$w,cl$s)
  cat(sprintf("\n-- %s (prev=%.1f%%, n=%d) --\n",nm,100*r$natprev,r$n)); print(round(r$tab,3))
  for(s in rownames(r$tab)) rows1[[length(rows1)+1]]<-data.frame(cell=nm,natprev=r$natprev,n=r$n,strategy=s,
    pearson=r$tab[s,"pearson"],spearman=r$tab[s,"spearman"],mae=r$tab[s,"mae"])}
write.csv(do.call(rbind,rows1),"results/sensitivity/distributional_within_heterosked.csv",row.names=FALSE)
rm(mg_all); gc(verbose=FALSE)

## ---- PART 2: LOCO transport (child iron, log-ferritin scale) -------------
cat("\n\n######## PART 2: LOCO transport, child iron ########\n")
ccfg<-list(gambia=list(ch="gw_LogFerAdj",logscale=TRUE,wt="gw_svy_weight"),
           ghana =list(ch="gw_cFerrAdjThurn",logscale=FALSE,wt="gw_sWeight"),
           sierraleone=list(ch="gw_cFerrAdj",logscale=FALSE,wt="gw_svy_weight"))
THR<-log(12)
build<-function(co){cf<-ccfg[[co]]; mg<-tar_read_raw(paste0("merged_",co)); mg<-mg[which(mg$gw_child_flag==1),]
  cont<-suppressWarnings(as.numeric(mg[[cf$ch]])); logc<-if(cf$logscale)cont else log(cont)
  ind<-data.frame(country=co,Admin2=as.character(mg$Admin2),logc=logc,wt=suppressWarnings(as.numeric(mg[[cf$wt]])),stringsAsFactors=FALSE)
  ind$wt[is.na(ind$wt)]<-1; ind<-ind[is.finite(ind$logc)&!is.na(ind$Admin2)&nzchar(ind$Admin2),]; ind$def<-as.integer(ind$logc<THR)
  ac<-tar_read_raw(paste0("area_covariates_",co))$gee_admin2; ac$Admin2<-as.character(ac$Admin2)
  sv<-as.data.frame(tar_read_raw(paste0("svy_admin2_",co,"_child_iron"))); sv$Admin2<-as.character(sv$Admin2)
  list(ind=ind,ac=ac,sv=sv)}
D<-lapply(names(ccfg),build); names(D)<-names(ccfg)
covset<-c("gee_soilzinc_mean_0_20","gee_soiliron_mean_0_20","gee_soilnitrogen_mean_0_20",
          "gee_soilphosphorus_mean_0_20","gee_accessibility_2019","gee_elevation_2000","gee_popdensity_2010")
covset<-covset[sapply(covset,function(cc) all(sapply(D,function(d) cc%in%colnames(d$ac))))]
attach_cov<-function(d){ac<-d$ac[,c("Admin2",covset)]; for(cc in covset)ac[[cc]]<-suppressWarnings(as.numeric(ac[[cc]])); merge(d$ind,ac,by="Admin2",all.x=TRUE)}
IND<-lapply(D,attach_cov)
loco<-function(holdout){set.seed(20260630)
  tr<-do.call(rbind,IND[setdiff(names(IND),holdout)]); te<-IND[[holdout]]
  P<-prep(as.matrix(tr[,covset,drop=FALSE]),as.matrix(te[,covset,drop=FALSE]))
  o<-predict_strats(P$tr,P$te,tr$def,tr$logc,tr$wt,THR); te$S1b<-o$S1b; te$S2<-o$S2; te$S3<-o$S3
  agg<-function(v) tapply(seq_along(v),te$Admin2,function(ix) weighted.mean(v[ix],te$wt[ix]))
  pa<-data.frame(Admin2=names(agg(te$S2)),S1b=agg(te$S1b),S2=agg(te$S2),S3=agg(te$S3))
  pa<-merge(pa,D[[holdout]]$sv[,c("Admin2","svy_prev","n_svy")],by="Admin2"); pa<-pa[pa$n_svy>=8&is.finite(pa$svy_prev),]
  sm<-function(p)c(spearman=cor(p,pa$svy_prev,method="spearman"),pearson=cor(p,pa$svy_prev),mae=mean(abs(p-pa$svy_prev)))
  list(n=nrow(pa),tab=rbind(S1b=sm(pa$S1b),S2=sm(pa$S2),S3=sm(pa$S3)))}
rows2<-list()
for(co in names(IND)){r<-loco(co); cat(sprintf("\n-- held-out: %s (n=%d districts) --\n",co,r$n)); print(round(r$tab,3))
  for(s in rownames(r$tab)) rows2[[length(rows2)+1]]<-data.frame(holdout=co,n=r$n,strategy=s,
    spearman=r$tab[s,"spearman"],pearson=r$tab[s,"pearson"],mae=r$tab[s,"mae"])}
write.csv(do.call(rbind,rows2),"results/sensitivity/distributional_transport_childiron.csv",row.names=FALSE)
cat("\nWrote results/sensitivity/{distributional_within_heterosked,distributional_transport_childiron}.csv\n")
