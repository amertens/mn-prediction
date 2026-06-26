#!/usr/bin/env Rscript
# Experiment: feature construction & selection (ranger workhorse).
#  - SEASONALITY: replace 12 monthly climatology cols/var with 6 summaries
#  - DROP COUNTS: remove population-scaled count columns (hurt transfer?)
# Evaluated within-country (full features) and LOCO (shared features).
source(here::here("sandbox_fe","harness.R"))
source(here::here("sandbox_fe","fe_features.R"))

outcomes <- c("child_vitA","women_iron","child_iron")

# feature recipes: function(M) -> M (column names preserved/extended)
recipes_fe <- list(
  raw          = function(M) M,
  seasonality  = function(M) {                       # drop raw monthly, add summaries
    mc <- is_monthly_col(colnames(M)); S <- build_seasonality(M)
    cbind(M[,!mc,drop=FALSE], S)
  },
  drop_counts  = function(M) M[, !is_count_col(colnames(M)), drop=FALSE],
  season_nocnt = function(M) {
    mc <- is_monthly_col(colnames(M)); S <- build_seasonality(M)
    M2 <- cbind(M[,!mc,drop=FALSE], S)
    M2[, !is_count_col(colnames(M2)), drop=FALSE]
  }
)

apply_recipe <- function(obj, rf) { obj$M <- rf(obj$M); obj }

rows_w <- list(); rows_l <- list()
prep <- build_prep("median", FALSE, "none", "none")  # ranger: impute only

for (otag in outcomes) {
  objs <- load_objs(otag)
  shared <- Reduce(intersect, lapply(objs, function(o) colnames(o$M)))
  for (rn in names(recipes_fe)) {
    rf <- recipes_fe[[rn]]
    objs_fe <- lapply(objs, apply_recipe, rf=rf)
    # within
    for (cn in names(objs_fe)) {
      a <- tryCatch(eval_within(objs_fe[[cn]], prep, "binary","ranger"),
                    error=function(e){cat("Werr",otag,rn,cn,conditionMessage(e),"\n");NA})
      rows_w[[length(rows_w)+1]] <- data.frame(outcome=otag,recipe=rn,country=cn,auc=round(a,3))
    }
    # LOCO on shared (recompute shared after recipe)
    shared_fe <- Reduce(intersect, lapply(objs_fe, function(o) colnames(o$M)))
    m <- tryCatch(eval_loco(objs_fe, prep, shared_fe, "binary","ranger", per_country_norm=FALSE),
                  error=function(e){cat("Lerr",otag,rn,conditionMessage(e),"\n");NULL})
    if(!is.null(m)) for(cn in names(m))
      rows_l[[length(rows_l)+1]] <- data.frame(outcome=otag,recipe=rn,test_country=cn,auc=round(m[[cn]],3))
  }
  cat("done",otag,"\n")
}
resw <- do.call(rbind, rows_w); resl <- do.call(rbind, rows_l)
write.csv(resw, here::here("sandbox_fe","results_06_within.csv"), row.names=FALSE)
write.csv(resl, here::here("sandbox_fe","results_06_loco.csv"), row.names=FALSE)

cat("\n=== WITHIN-COUNTRY (ranger) mean AUC ===\n")
ag<-aggregate(auc~outcome+recipe,resw,function(x)round(mean(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag[ag$outcome==o,];s<-s[order(-s$auc),];cat("--",o,"--\n");print(s[,c("recipe","auc")],row.names=FALSE)}
cat("\n=== LOCO TRANSFER (ranger) mean AUC ===\n")
ag<-aggregate(auc~outcome+recipe,resl,function(x)round(mean(x,na.rm=TRUE),3))
for(o in outcomes){s<-ag[ag$outcome==o,];s<-s[order(-s$auc),];cat("--",o,"--\n");print(s[,c("recipe","auc")],row.names=FALSE)}
