#!/usr/bin/env Rscript
# ===========================================================================
# 10_distributional_vs_binary_prototype.R
#
# Question: does modelling the CONTINUOUS biomarker distribution and
# integrating to the deficiency cut-point recover statistical power versus
# dichotomising first? And how much of any gain is simply about modelling
# individuals instead of noisy area proportions?
#
# Clean contrast, identical covariates + identical leave-one-AREA-out CV:
#   S1   ridge on the survey-weighted district PROPORTION (logit)  -- the
#        status-quo "model area-level prevalence from proxies" analog.
#   S1b  ridge LOGISTIC on the individual binary deficiency outcome -- fair
#        binary baseline that uses individual data.
#   S2   ridge on the individual CONTINUOUS log-biomarker, prevalence recovered
#        by integrating the (empirical-residual) distribution past the cut-point.
# If S2 > S1b, the only difference is "dichotomise first" vs "model the
# distribution". If S1b >> S1, the lever is individual- vs area-level modelling.
#
# Run from repo root:  Rscript sensitivity/10_distributional_vs_binary_prototype.R
# Reads the _targets_full store; writes results/sensitivity/distributional_prototype_metrics.csv
# ===========================================================================
suppressMessages({library(targets); library(glmnet)})
set.seed(20260630)
tar_config_set(store = "_targets_full")

COUNTRY <- "ghana"   # prototype country (merged_<country>, area_covariates_<country>, svy_admin2_<country>_*)
MIN_NSVY <- 8        # drop areas with tiny survey n from the evaluation
B_BOOT   <- 2000     # area bootstrap reps for the S2-vs-S1b contrast CI

mg_all <- tar_read_raw(paste0("merged_", COUNTRY))
ac     <- tar_read_raw(paste0("area_covariates_", COUNTRY))$gee_admin2

# --- a-priori fixed, modest covariate set (vitA/iron environmental pathway) ---
pick <- function(rx){h<-grep(rx,colnames(ac),value=TRUE); if(length(h)) h[1] else NA}
cov_names <- na.omit(c(pick("gee_soilzinc_mean_0_20"), pick("gee_soiliron_mean_0_20"),
  pick("gee_soilnitrogen_mean_0_20"), pick("gee_soilphosphorus_mean_0_20"),
  pick("gee_accessibility_2019"), pick("gee_elevation_2000"),
  pick("gee_popdensity_2010"), pick("gee_ghsbuilts_2015_built_surface")))
acov <- ac[, c("Admin2", cov_names)]; acov$Admin2 <- as.character(acov$Admin2)
for (cc in cov_names) acov[[cc]] <- suppressWarnings(as.numeric(acov[[cc]]))

# --- in-fold standardize + median impute (train statistics only) ---
prep <- function(Xtr, Xte){
  med<-apply(Xtr,2,median,na.rm=TRUE); med[!is.finite(med)]<-0
  ctr<-colMeans(Xtr,na.rm=TRUE); ctr[!is.finite(ctr)]<-0
  sdv<-apply(Xtr,2,sd,na.rm=TRUE); sdv[!is.finite(sdv)|sdv==0]<-1
  imp<-function(M){for(j in seq_len(ncol(M))){x<-M[,j];x[is.na(x)]<-med[j];M[,j]<-(x-ctr[j])/sdv[j]};M}
  list(tr=imp(Xtr), te=imp(Xte))
}
rfit <- function(x,y,fam,w=NULL){ if(is.null(w)) w<-rep(1,length(y))
  tryCatch(cv.glmnet(x,y,family=fam,alpha=0,weights=w,nfolds=5), error=function(e) NULL) }

run_cell <- function(cont_col, thresh, child, wt_col, svy){
  mg <- mg_all
  if (child && "gw_child_flag" %in% names(mg))  mg <- mg[which(mg$gw_child_flag==1),]
  if (!child && "gw_child_flag" %in% names(mg)) mg <- mg[which(mg$gw_child_flag==0),]
  ind <- data.frame(Admin2=as.character(mg$Admin2),
                    cont=suppressWarnings(as.numeric(mg[[cont_col]])),
                    wt=suppressWarnings(as.numeric(mg[[wt_col]])), stringsAsFactors=FALSE)
  ind$wt[is.na(ind$wt)] <- 1
  ind <- ind[is.finite(ind$cont) & ind$cont>0 & !is.na(ind$Admin2) & nzchar(ind$Admin2),]
  ind$def <- as.integer(ind$cont < thresh); ind$logc <- log(ind$cont)
  sv <- as.data.frame(tar_read_raw(svy)); sv$Admin2 <- as.character(sv$Admin2)
  areas <- sort(unique(ind$Admin2)); areas <- areas[areas %in% acov$Admin2 & areas %in% sv$Admin2]
  A <- acov[match(areas,acov$Admin2),]; sv <- sv[match(areas,sv$Admin2),]; ind <- ind[ind$Admin2 %in% areas,]
  X <- as.matrix(A[,cov_names,drop=FALSE]); nA <- length(areas)
  wprev <- sapply(areas, function(a){d<-ind[ind$Admin2==a,]; sum(d$wt*d$def)/sum(d$wt)})
  pr <- data.frame(Admin2=areas, truth=sv$svy_prev, n_svy=sv$n_svy, S1=NA, S1b=NA, S2=NA)
  for (i in seq_len(nA)){
    Xa <- prep(X[-i,,drop=FALSE], X[i,,drop=FALSE])
    yv <- pmin(pmax(wprev[-i],.01),.99); f1 <- rfit(Xa$tr, qlogis(yv), "gaussian")
    if(!is.null(f1)) pr$S1[i] <- plogis(as.numeric(predict(f1, Xa$te, s="lambda.min")))
    itr <- ind[ind$Admin2 %in% areas[-i],]; ite <- ind[ind$Admin2==areas[i],]
    Pi <- prep(X[match(itr$Admin2,areas),,drop=FALSE], X[match(ite$Admin2,areas),,drop=FALSE])
    f1b <- rfit(Pi$tr, itr$def, "binomial", itr$wt)
    if(!is.null(f1b)) pr$S1b[i] <- mean(as.numeric(predict(f1b, Pi$te, s="lambda.min", type="response")))
    f2 <- rfit(Pi$tr, itr$logc, "gaussian", itr$wt)
    if(!is.null(f2)){ mu_tr<-as.numeric(predict(f2,Pi$tr,s="lambda.min"))
      mu_te<-as.numeric(predict(f2,Pi$te,s="lambda.min")); resid<-itr$logc-mu_tr; thr<-log(thresh)
      pr$S2[i] <- mean(sapply(mu_te, function(m) weighted.mean((m+resid)<thr, itr$wt))) }
  }
  ev <- pr[pr$n_svy>=MIN_NSVY & is.finite(pr$truth),]
  m  <- function(p) c(pearson=cor(p,ev$truth), spearman=cor(p,ev$truth,method="spearman"), mae=mean(abs(p-ev$truth)))
  boot <- replicate(B_BOOT,{ idx<-sample(seq_len(nrow(ev)),replace=TRUE)
    cor(ev$S2[idx],ev$truth[idx]) - cor(ev$S1b[idx],ev$truth[idx]) })
  list(natprev=weighted.mean(ind$def,ind$wt), n=nrow(ev),
       tab=rbind(S1=m(ev$S1), S1b=m(ev$S1b), S2=m(ev$S2)),
       gain=cor(ev$S2,ev$truth)-cor(ev$S1b,ev$truth),
       gain_ci=quantile(boot,c(.025,.975)))
}

cells <- list(
  child_vitA = list(cont="gw_cRBPAdjThurn",  thr=0.70, child=TRUE,  wt="gw_child_weight", svy=paste0("svy_admin2_",COUNTRY,"_child_vitA")),
  child_iron = list(cont="gw_cFerrAdjThurn", thr=12,   child=TRUE,  wt="gw_child_weight", svy=paste0("svy_admin2_",COUNTRY,"_child_iron")),
  women_iron = list(cont="gw_wFerrAdjThurn", thr=15,   child=FALSE, wt="gw_sWeight",      svy=paste0("svy_admin2_",COUNTRY,"_women_iron"))
)

rows <- list()
for (nm in names(cells)){
  cl <- cells[[nm]]
  r <- tryCatch(run_cell(cl$cont,cl$thr,cl$child,cl$wt,cl$svy), error=function(e){cat("ERR",nm,":",conditionMessage(e),"\n");NULL})
  if (is.null(r)) next
  cat(sprintf("\n===== %s  (national prev=%.1f%%, eval n=%d areas) =====\n", nm, 100*r$natprev, r$n))
  print(round(r$tab,3))
  cat(sprintf("  S2(continuous) - S1b(binary) Pearson gain: %+.3f  [%.3f, %.3f]\n",
              r$gain, r$gain_ci[1], r$gain_ci[2]))
  for (s in rownames(r$tab))
    rows[[length(rows)+1]] <- data.frame(cell=nm, natprev=r$natprev, n=r$n, strategy=s,
      pearson=r$tab[s,"pearson"], spearman=r$tab[s,"spearman"], mae=r$tab[s,"mae"],
      s2_minus_s1b_gain=ifelse(s=="S2",r$gain,NA), gain_lo=ifelse(s=="S2",r$gain_ci[1],NA),
      gain_hi=ifelse(s=="S2",r$gain_ci[2],NA))
}
dir.create("results/sensitivity", showWarnings=FALSE, recursive=TRUE)
out <- do.call(rbind, rows)
write.csv(out, "results/sensitivity/distributional_prototype_metrics.csv", row.names=FALSE)
cat("\nWrote results/sensitivity/distributional_prototype_metrics.csv\n")
