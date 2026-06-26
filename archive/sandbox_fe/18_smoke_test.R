#!/usr/bin/env Rscript
# Smoke-test the drafted production changes:
#  1. all edited files parse
#  2. bundle_prefixes_for_outcome + build_outcome_dataset$Xvars_bundle
#  3. DHS_SL_clustered runs with normalize_method = zscore AND rank
#  4. rank recipe bakes NEW data with no NAs (prediction-path safety)
suppressMessages({library(here);library(dplyr)})
source(here::here("R","config.R"))
source(here::here("R","data_prep.R"))
source(here::here("R", "sensitivity", "sl_fitting.R"))
source(here::here("src","analysis","sl_helpers.R"))
source(here::here("src","0-functions.R"))
source(here::here("R","feature_engineering_constructs.R"))
cat("[1] all files sourced OK\n")

# bundle mapping
cat("[2] bundle prefixes:\n")
for (t in c("child_vitA","women_iron","women_b12","child_zinc","mystery"))
  cat(sprintf("    %-12s -> %s\n", t, paste(bundle_prefixes_for_outcome(t), collapse=",")))

# build a real outcome dataset (Ghana child_vitA) to check Xvars_bundle
cfgs <- get_country_configs(); cc <- cfgs[["Ghana"]]
merged <- load_merged_data(cc$data_path)
od <- build_outcome_dataset(merged, cc, cc$outcomes[["child_vitA"]])
cat(sprintf("[2] Ghana child_vitA: Xvars=%d  Xvars_bundle=%d (env-forward)\n",
            length(od$Xvars), length(od$Xvars_bundle)))
stopifnot(length(od$Xvars_bundle) >= 2, length(od$Xvars_bundle) <= length(od$Xvars))

# minimal SL stack + cfg global
params <- get_pipeline_params("fast")
cat(sprintf("[cfg] normalize_method=%s use_outcome_bundles=%s\n",
            params$normalize_method, params$use_outcome_bundles))
assign("cfg", list(seed=params$seed, admin1=cc$admin1_col,
                   cluster_id=cc$cluster_id, K=params$K), envir=globalenv())
learners <- setup_sl_learners(params)

d <- od$data
for (nm in c("zscore","rank")) {
  t0 <- proc.time()
  fit <- DHS_SL_clustered(d=d, Xvars=od$Xvars_bundle, outcome=cc$outcomes$child_vitA$binary,
                          population="children", id=cc$cluster_id, folds=5L,
                          sl=learners$slmod2_bin, normalize_method=nm)
  el <- (proc.time()-t0)["elapsed"]
  y <- fit$res$Y; p <- fit$res$yhat_full; ok <- !is.na(y)&!is.na(p)
  n1<-sum(y[ok]==1);n0<-sum(y[ok]==0);r<-rank(p[ok]);auc<-(sum(r[y[ok]==1])-n1*(n1+1)/2)/(n1*n0)
  # prediction-path safety: bake the saved recipe on the SAME cov schema
  nNA <- sum(is.na(recipes::bake(fit$auto_recipe,
              new_data = as.data.frame(setNames(lapply(fit$pre_recipe_cols, function(.) rnorm(20)),
                                                fit$pre_recipe_cols)))))
  cat(sprintf("[3] normalize=%-6s p=%d AUC=%.3f  %.0fs  bake-NA=%d\n",
              nm, length(fit$Xvars), auc, el, nNA))
}
cat("[DONE] smoke test passed\n")
