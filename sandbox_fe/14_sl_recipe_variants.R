#!/usr/bin/env Rscript
# (a) Confirm recipe changes with the PRODUCTION sl3 SuperLearner stack.
# Compare current(z-score) vs rank vs per-domain PCA vs rank+PCA on real SL.
source(here::here("sandbox_fe","harness.R"))     # load_objs (numeric caches)
source(here::here("sandbox_fe","sl_variant.R"))
source(here::here("R","config.R"))
source(here::here("R", "sensitivity", "sl_fitting.R"))

stack <- Sys.getenv("SL_STACK","fast")
params <- get_pipeline_params(if(stack=="full")"full" else "fast")
learners <- setup_sl_learners(params)

# representative (country, outcome) cells spanning prevalence + signal regimes
cells <- list(
  c("Gambia","women_iron"), c("Ghana","child_vitA"), c("Ghana","child_iron"),
  c("Malawi","child_vitA"), c("Gambia","child_iron")
)
variants <- c("current","rank","pca","rankpca")

build_d <- function(o){
  d <- as.data.frame(o$M); d$Y <- o$y_bin; d$clusterid <- o$cluster
  list(d=d, Xvars=colnames(o$M))
}

rows <- list()
for (cell in cells) {
  cn <- cell[1]; otag <- cell[2]
  o <- load_objs(otag, countries=cn)[[cn]]
  bd <- build_d(o)
  for (v in variants) {
    t0 <- proc.time()
    r <- tryCatch(run_sl_variant(bd$d, bd$Xvars, "Y", "clusterid",
                                 sl=learners$slmod2_bin, folds=5L, variant=v),
                  error=function(e){cat("ERR",cn,otag,v,conditionMessage(e),"\n");NULL})
    el <- (proc.time()-t0)["elapsed"]
    if (is.null(r)) next
    rows[[length(rows)+1]] <- data.frame(country=cn, outcome=otag, variant=v,
                                         auc=round(r$auc,3), p=r$p, sec=round(el,1))
    cat(sprintf("%-8s %-11s %-8s AUC=%.3f  p=%d  %.0fs\n", cn, otag, v, r$auc, r$p, el))
  }
}
res <- do.call(rbind, rows)
write.csv(res, here::here("sandbox_fe", sprintf("results_14_sl_recipe_%s.csv", stack)), row.names=FALSE)
cat("\n=== mean AUC across cells, by variant (real SL,", stack, "stack) ===\n")
ag <- aggregate(auc~variant, res, function(x) round(mean(x,na.rm=TRUE),3))
ag <- ag[order(-ag$auc),]; print(ag, row.names=FALSE)
cat("\n=== delta vs current, per cell ===\n")
w <- reshape(res[,c("country","outcome","variant","auc")], idvar=c("country","outcome"),
             timevar="variant", direction="wide")
names(w) <- sub("auc.","",names(w))
for (v in setdiff(variants,"current")) w[[paste0("d_",v)]] <- round(w[[v]]-w[["current"]],3)
print(w, row.names=FALSE)
