#!/usr/bin/env Rscript
# A/B: predictor hygiene (counts->rates + MAP snapshot dedup) ON vs OFF.
# Baseline = all proxy cols; Hygiene = drop prune_predictor_cols().
# Within-country CV (ranger + lasso-rank) and area-level LOCO.
source(here::here("sandbox_fe","harness.R"))
source(here::here("R","config.R"))
source(here::here("R","data_prep.R"))   # prune_predictor_cols()

cfgs <- get_country_configs()
SURVEY_YEAR <- sapply(names(cfgs), function(c) cfgs[[c]]$survey_year)
outcomes <- c("child_vitA","women_iron","child_iron")

rows_w <- list(); add_w <- function(...) rows_w[[length(rows_w)+1]] <<- data.frame(...)
for (otag in outcomes) {
  objs <- load_objs(otag)
  for (cn in names(objs)) {
    o <- objs[[cn]]; allc <- colnames(o$M)
    drop <- prune_predictor_cols(allc, survey_year = SURVEY_YEAR[[cn]])
    keep <- setdiff(allc, drop)
    for (mdl in c("ranger","glmnet")) {
      prep <- if (mdl=="ranger") build_prep("median",FALSE,"none","none") else build_prep("median",FALSE,"rank","none")
      lam  <- if (mdl=="glmnet") 0.02 else NULL
      a_all <- tryCatch(eval_within(o, prep, "binary", mdl, alpha=1, cols=allc, fixed_lambda=lam), error=function(e)NA)
      a_hy  <- tryCatch(eval_within(o, prep, "binary", mdl, alpha=1, cols=keep, fixed_lambda=lam), error=function(e)NA)
      add_w(outcome=otag, country=cn, model=mdl, p_all=length(allc), p_hy=length(keep),
            auc_all=round(a_all,3), auc_hy=round(a_hy,3), delta=round(a_hy-a_all,3))
    }
  }
  cat("within done:", otag, "\n")
}
resw <- do.call(rbind, rows_w)
write.csv(resw, here::here("sandbox_fe","results_23_hygiene_within.csv"), row.names=FALSE)

# ---- area-level LOCO (lasso), baseline vs hygiene on shared set ----
aggregate_area <- function(o,cols,min_n=5){y<-o$y_bin;ok<-!is.na(y);a2<-o$admin2[ok];y<-y[ok];M<-o$M[ok,cols,drop=FALSE]
  ua<-unique(a2);keep<-ua[sapply(ua,function(u)sum(a2==u)>=min_n)]
  list(prev=sapply(keep,function(u)mean(y[a2==u])),X=t(sapply(keep,function(u){idx<-a2==u;colMeans(M[idx,,drop=FALSE],na.rm=TRUE)})))}
rows_l <- list()
for (otag in outcomes) {
  objs <- load_objs(otag); shared <- Reduce(intersect, lapply(objs,function(o)colnames(o$M)))
  # hygiene drop on the shared pooled set (use median survey year for snapshot dedup)
  drop_sh <- prune_predictor_cols(shared, survey_year = as.integer(median(unlist(SURVEY_YEAR))))
  keep_sh <- setdiff(shared, drop_sh)
  loco <- function(cols){
    sp<-c()
    for(test_cn in names(objs)){tr<-setdiff(names(objs),test_cn)
      areas<-lapply(objs,aggregate_area,cols=cols)
      Xtr<-do.call(rbind,lapply(tr,function(c)areas[[c]]$X));ytr<-unlist(lapply(tr,function(c)areas[[c]]$prev))
      Xte<-areas[[test_cn]]$X;yte<-areas[[test_cn]]$prev; if(length(yte)<4) next
      cc<-.clean(Xtr,Xte);im<-imp_median(cc$tr,cc$te);z<-tf_zscore(im$tr,im$te)
      pr<-fit_predict_glmnet(as.matrix(z$tr),ytr,as.matrix(z$te),1,"gaussian",fixed_lambda=0.01)
      sp<-c(sp,suppressWarnings(cor(yte,pr,method="spearman",use="complete.obs")))}
    mean(sp,na.rm=TRUE)
  }
  rows_l[[length(rows_l)+1]] <- data.frame(outcome=otag, p_all=length(shared), p_hy=length(keep_sh),
    sp_all=round(loco(shared),3), sp_hy=round(loco(keep_sh),3))
  cat("loco done:", otag, "\n")
}
resl <- do.call(rbind, rows_l)
resl$delta <- round(resl$sp_hy - resl$sp_all, 3)
write.csv(resl, here::here("sandbox_fe","results_23_hygiene_loco.csv"), row.names=FALSE)

cat("\n=== WITHIN-COUNTRY: hygiene ON vs OFF (CV AUC) ===\n")
print(resw[,c("outcome","country","model","p_all","p_hy","auc_all","auc_hy","delta")], row.names=FALSE)
cat("\n  mean delta by model:\n"); print(aggregate(delta~model, resw, function(x) round(mean(x,na.rm=TRUE),4)), row.names=FALSE)
cat("\n=== AREA-LEVEL LOCO: hygiene ON vs OFF (Spearman) ===\n")
print(resl, row.names=FALSE)
