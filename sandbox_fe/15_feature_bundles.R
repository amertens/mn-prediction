#!/usr/bin/env Rscript
# (b) Outcome-specific feature bundles from the domain-biology finding.
#   vitA_env    = environment + soil + food security  (gee_, soil_, fsec_)
#   iron_health = malaria + modelled health + food security (MAP_, ihme_, fsec_)
# Tests:
#   1. WITHIN-country, real SL (fast): all vs matched bundle vs MISMATCHED bundle
#      (mismatch should be worse -> proves biological specificity, not just size)
#   2. AREA-LEVEL LOCO (lasso): matched bundle vs all (transfer)
source(here::here("sandbox_fe","harness.R"))
source(here::here("sandbox_fe","sl_variant.R"))
source(here::here("R","config.R")); source(here::here("R", "sensitivity", "sl_fitting.R"))

BUNDLES <- list(
  all         = NULL,
  vitA_env    = c("^gee_","^soil_","^fsec_"),
  iron_health = c("^MAP_","^ihme_","^fsec_")
)
bundle_cols <- function(cols, prefixes) if (is.null(prefixes)) cols else
  cols[Reduce(`|`, lapply(prefixes, function(p) grepl(p, cols)))]
matched_bundle <- function(otag) if (grepl("vitA",otag)) "vitA_env" else "iron_health"
mismatched_bundle <- function(otag) if (grepl("vitA",otag)) "iron_health" else "vitA_env"

# ---- 1. within-country real SL --------------------------------------------
params <- get_pipeline_params("fast"); learners <- setup_sl_learners(params)
cells <- list(c("Ghana","child_vitA"), c("Malawi","child_vitA"),
              c("Gambia","women_iron"), c("Ghana","child_iron"))
build_d <- function(o){ d<-as.data.frame(o$M); d$Y<-o$y_bin; d$clusterid<-o$cluster; list(d=d,Xvars=colnames(o$M)) }

rows_w <- list()
for (cell in cells) {
  cn<-cell[1]; otag<-cell[2]
  o <- load_objs(otag, countries=cn)[[cn]]; bd <- build_d(o)
  sets <- c("all", matched_bundle(otag), mismatched_bundle(otag))
  for (bn in sets) {
    cols <- bundle_cols(bd$Xvars, BUNDLES[[bn]])
    r<-tryCatch(run_sl_variant(bd$d, bd$Xvars, "Y","clusterid", sl=learners$slmod2_bin,
                               folds=5L, variant="current", bundle=cols),
                error=function(e){cat("ERR",cn,otag,bn,conditionMessage(e),"\n");NULL})
    if(is.null(r)) next
    role <- if(bn=="all") "all" else if(bn==matched_bundle(otag)) "MATCHED" else "mismatched"
    rows_w[[length(rows_w)+1]] <- data.frame(country=cn,outcome=otag,bundle=bn,role=role,
                                             auc=round(r$auc,3),p=r$p)
    cat(sprintf("%-8s %-11s %-12s(%-10s) AUC=%.3f p=%d\n",cn,otag,bn,role,r$auc,r$p))
  }
}
resw<-do.call(rbind,rows_w)
write.csv(resw, here::here("sandbox_fe","results_15_bundles_within.csv"), row.names=FALSE)

# ---- 2. area-level LOCO (lasso) -------------------------------------------
aggregate_area <- function(o,cols,min_n=5){y<-o$y_bin;ok<-!is.na(y);a2<-o$admin2[ok];y<-y[ok];M<-o$M[ok,cols,drop=FALSE]
  ua<-unique(a2);keep<-ua[sapply(ua,function(u)sum(a2==u)>=min_n)]
  list(prev=sapply(keep,function(u)mean(y[a2==u])),X=t(sapply(keep,function(u){idx<-a2==u;colMeans(M[idx,,drop=FALSE],na.rm=TRUE)})))}
rows_l<-list()
for (otag in c("child_vitA","women_iron","child_iron")) {
  objs<-load_objs(otag); shared<-Reduce(intersect,lapply(objs,function(o)colnames(o$M)))
  for (bn in c("all", matched_bundle(otag))) {
    cols<-bundle_cols(shared, BUNDLES[[bn]]); if(length(cols)<2) next
    areas<-lapply(objs,aggregate_area,cols=cols); sps<-c()
    for(test_cn in names(areas)){tr<-setdiff(names(areas),test_cn)
      Xtr<-do.call(rbind,lapply(tr,function(c)areas[[c]]$X));ytr<-unlist(lapply(tr,function(c)areas[[c]]$prev))
      Xte<-areas[[test_cn]]$X;yte<-areas[[test_cn]]$prev; if(length(yte)<4) next
      cc<-.clean(Xtr,Xte);im<-imp_median(cc$tr,cc$te);z<-tf_zscore(im$tr,im$te)
      pr<-fit_predict_glmnet(as.matrix(z$tr),ytr,as.matrix(z$te),1,"gaussian",fixed_lambda=0.01)
      sps<-c(sps,suppressWarnings(cor(yte,pr,method="spearman",use="complete.obs")))}
    rows_l[[length(rows_l)+1]]<-data.frame(outcome=otag,bundle=bn,n_feat=length(cols),
                                           mean_spearman=round(mean(sps,na.rm=TRUE),3))
  }
}
resl<-do.call(rbind,rows_l)
write.csv(resl, here::here("sandbox_fe","results_15_bundles_loco.csv"), row.names=FALSE)

cat("\n=== WITHIN-COUNTRY (real SL): matched bundle vs all vs mismatched ===\n")
print(resw[,c("country","outcome","role","bundle","auc","p")], row.names=FALSE)
cat("\n=== AREA-LEVEL LOCO transfer (lasso): matched bundle vs all ===\n")
print(resl, row.names=FALSE)
