#!/usr/bin/env Rscript
# Experiment: within-country CV. Two model families.
#  (A) ranger (transform-invariant): test imputation x dimension-reduction
#  (B) glmnet fixed-lambda (linear):  test transforms x dimension-reduction
source(here::here("sandbox_fe","harness.R"))

outcomes <- c("child_vitA","women_iron","child_iron")
all_objs <- lapply(outcomes, load_objs); names(all_objs) <- outcomes

run_grid <- function(grid, model, alpha=1, lam=NULL) {
  rows <- list()
  for (otag in outcomes) {
    objs <- all_objs[[otag]]
    for (g in grid) {
      prep <- build_prep(g$impute, g$add_ind, g$transform, g$reduce)
      for (cn in names(objs)) {
        m <- tryCatch(eval_within(objs[[cn]], prep, "binary", model, alpha, fixed_lambda=lam),
                      error=function(e){cat("ERR",otag,g$tag,cn,conditionMessage(e),"\n");NA})
        rows[[length(rows)+1]] <- data.frame(outcome=otag,variant=g$tag,country=cn,auc=round(m,3))
      }
    }
    cat("  done",model,otag,"\n")
  }
  do.call(rbind, rows)
}

# (A) ranger: imputation x reduction (transform irrelevant for trees -> use none)
gridA <- list(
  list(tag="median",            impute="median",add_ind=FALSE,transform="none",reduce="none"),
  list(tag="median+ind",        impute="median",add_ind=TRUE, transform="none",reduce="none"),
  list(tag="mean",              impute="mean",  add_ind=FALSE,transform="none",reduce="none"),
  list(tag="zero",              impute="zero",  add_ind=FALSE,transform="none",reduce="none"),
  list(tag="knn",               impute="knn",   add_ind=FALSE,transform="none",reduce="none"),
  list(tag="median+corr0.85",   impute="median",add_ind=FALSE,transform="none",reduce="corr"),
  list(tag="median+ind+corr",   impute="median",add_ind=TRUE, transform="none",reduce="corr")
)
resA <- run_grid(gridA, "ranger")
resA$model <- "ranger"

# (B) glmnet fixed-lambda: transforms x reduction
gridB <- list(
  list(tag="zscore+corr",   impute="median",add_ind=FALSE,transform="zscore",   reduce="corr"),
  list(tag="rank+corr",     impute="median",add_ind=FALSE,transform="rank",     reduce="corr"),
  list(tag="signedlog+corr",impute="median",add_ind=FALSE,transform="signedlog",reduce="corr"),
  list(tag="zscore+nored",  impute="median",add_ind=FALSE,transform="zscore",   reduce="none"),
  list(tag="rank+nored",    impute="median",add_ind=FALSE,transform="rank",     reduce="none")
)
resB <- run_grid(gridB, "glmnet", alpha=1, lam=0.02)
resB$model <- "glmnet_lasso"

res <- rbind(resA, resB)
write.csv(res, here::here("sandbox_fe","results_04_within.csv"), row.names=FALSE)

summ <- function(r,label){
  cat("\n===",label,": mean AUC across countries ===\n")
  ag <- aggregate(auc~outcome+variant, r, function(x) round(mean(x,na.rm=TRUE),3))
  for(o in outcomes){ s<-ag[ag$outcome==o,]; s<-s[order(-s$auc),]
    cat("--",o,"--\n"); print(s[,c("variant","auc")], row.names=FALSE) }
}
summ(resA,"RANGER (imputation x reduction)")
summ(resB,"GLMNET LASSO (transform x reduction)")
