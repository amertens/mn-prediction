#!/usr/bin/env Rscript
# (a) Full-stack confirmation on two key cells: does rank hold, how does PCA behave?
source(here::here("sandbox_fe","harness.R"))
source(here::here("sandbox_fe","sl_variant.R"))
source(here::here("R","config.R")); source(here::here("R", "sensitivity", "sl_fitting.R"))

params <- get_pipeline_params("full")
learners <- setup_sl_learners(params)
cells <- list(c("Ghana","child_vitA"), c("Gambia","women_iron"))
variants <- c("current","rank","pca")
build_d <- function(o){ d<-as.data.frame(o$M); d$Y<-o$y_bin; d$clusterid<-o$cluster; list(d=d,Xvars=colnames(o$M)) }

rows <- list()
for (cell in cells) {
  o <- load_objs(cell[2], countries=cell[1])[[cell[1]]]; bd <- build_d(o)
  for (v in variants) {
    t0<-proc.time()
    r<-tryCatch(run_sl_variant(bd$d,bd$Xvars,"Y","clusterid",sl=learners$slmod2_bin,folds=5L,variant=v),
                error=function(e){cat("ERR",cell[1],cell[2],v,conditionMessage(e),"\n");NULL})
    el<-(proc.time()-t0)["elapsed"]; if(is.null(r)) next
    rows[[length(rows)+1]]<-data.frame(country=cell[1],outcome=cell[2],variant=v,auc=round(r$auc,3),p=r$p,sec=round(el,0))
    cat(sprintf("[FULL] %-8s %-11s %-8s AUC=%.3f p=%d %.0fs\n",cell[1],cell[2],v,r$auc,r$p,el))
  }
}
res<-do.call(rbind,rows)
write.csv(res, here::here("sandbox_fe","results_14b_fullstack.csv"), row.names=FALSE)
print(res,row.names=FALSE)
