#!/usr/bin/env Rscript
# Experiment: leave-one-country-out (LOCO) transfer on the SHARED feature set.
# Central question: do features need to be made country-RELATIVE (per-country
# normalized) to transfer, vs pooled on absolute scale?
source(here::here("sandbox_fe","harness.R"))
`%||%` <- function(a,b) if(is.null(a)) b else a

outcomes <- c("child_vitA","women_iron","child_iron")
rows <- list()

for (otag in outcomes) {
  objs <- load_objs(otag)
  shared <- Reduce(intersect, lapply(objs, function(o) colnames(o$M)))
  cat(sprintf("\n%s: %d countries, %d shared cols\n", otag, length(objs), length(shared)))

  configs <- list(
    list(tag="lasso pooled-zscore",       model="glmnet", pcn=FALSE, prep=build_prep("median",FALSE,"zscore","none"), lam=0.02),
    list(tag="lasso percountry-rank",     model="glmnet", pcn=TRUE,  prep=build_prep("median",FALSE,"none","none"),   lam=0.02),
    list(tag="ranger pooled-raw",         model="ranger", pcn=FALSE, prep=build_prep("median",FALSE,"none","none")),
    list(tag="ranger percountry-rank",    model="ranger", pcn=TRUE,  prep=build_prep("median",FALSE,"none","none"))
  )
  for (cfg in configs) {
    m <- tryCatch(eval_loco(objs, cfg$prep, shared, "binary", cfg$model,
                            per_country_norm=cfg$pcn, fixed_lambda=cfg$lam %||% NULL),
                  error=function(e){cat("ERR",otag,cfg$tag,conditionMessage(e),"\n");NULL})
    if (is.null(m)) next
    for (cn in names(m))
      rows[[length(rows)+1]] <- data.frame(outcome=otag, config=cfg$tag, test_country=cn, auc=round(m[[cn]],3))
  }
}
res <- do.call(rbind, rows)
write.csv(res, here::here("sandbox_fe","results_05_loco.csv"), row.names=FALSE)

cat("\n=== LOCO transfer AUC (mean across held-out countries) ===\n")
ag <- aggregate(auc~outcome+config, res, function(x) round(mean(x,na.rm=TRUE),3))
for(o in outcomes){ s<-ag[ag$outcome==o,]; s<-s[order(-s$auc),]
  cat("--",o,"--\n"); print(s[,c("config","auc")], row.names=FALSE) }
cat("\n=== Full per-country detail ===\n"); print(res, row.names=FALSE)
