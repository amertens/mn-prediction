#!/usr/bin/env Rscript
# Diagnostic: are proxy predictors area-level (constant within Admin2/cluster)?
# If so, the EFFECTIVE n for learning proxy->outcome is the number of distinct
# predictor rows (areas), not the number of individuals.
source(here::here("sandbox_fe","harness.R"))
source(here::here("sandbox_fe","fe_features.R"))

outcomes <- c("child_vitA","women_iron","child_iron")
for (otag in outcomes) {
  objs <- load_objs(otag)
  cat(sprintf("\n=== %s ===\n", otag))
  for (cn in names(objs)) {
    o <- objs[[cn]]
    M <- o$M
    # use shared-ish env cols to gauge distinct area signatures
    sig <- apply(round(M[, seq_len(min(50,ncol(M)))], 6), 1, paste, collapse="|")
    n_ind   <- nrow(M)
    n_admin2<- length(unique(o$admin2))
    n_clust <- length(unique(o$cluster))
    n_rows  <- length(unique(sig))
    # fraction of GEE variance that is between-cluster (ICC-like) for a few vars
    geecols <- grep("^gee_", colnames(M), value=TRUE)[1:min(20,sum(grepl("^gee_",colnames(M))))]
    icc <- sapply(geecols, function(g){
      x<-M[,g]; cl<-o$cluster; ok<-!is.na(x)
      if(sum(ok)<10||length(unique(cl[ok]))<2) return(NA)
      vt<-var(x[ok]); if(vt<1e-12) return(NA)
      vb<-var(tapply(x[ok], cl[ok], mean)); vb/vt
    })
    cat(sprintf("  %-12s n_ind=%4d  n_cluster=%3d  n_admin2=%3d  distinct_predrows=%4d  med_GEE_betweenClustVar=%.2f\n",
                cn, n_ind, n_clust, n_admin2, n_rows, median(icc, na.rm=TRUE)))
  }
}
