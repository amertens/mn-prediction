#!/usr/bin/env Rscript
# Profile spatial granularity of each proxy column:
#  distinct non-NA values relative to #clusters and #admin2.
# Classify: ADMIN2-level (<= ~n_admin2 distinct), CLUSTER-level (~n_cluster),
#  INDIVIDUAL/fine (>> n_cluster).
source(here::here("sandbox_fe","harness.R"))

otag <- "child_vitA"
objs <- load_objs(otag)
dom_of <- function(cn){pref<-c(DHS="^dhs",MICS="^mics_",IHME="^ihme_",LSMS="^lsms_",MAP="^MAP_",WFP="^wfp_",FLUNET="^flunet_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}

for (cn in names(objs)) {
  o <- objs[[cn]]; M <- o$M
  nc <- length(unique(o$cluster)); na2 <- length(unique(o$admin2))
  ndist <- apply(M, 2, function(x) length(unique(round(x[!is.na(x)],8))))
  level <- ifelse(ndist <= na2 + 1, "admin2",
            ifelse(ndist <= nc + 1, "cluster", "fine"))
  dom <- dom_of(colnames(M))
  cat(sprintf("\n== %s (clusters=%d admin2=%d) ==\n", cn, nc, na2))
  print(table(level))
  cat("by domain x level:\n")
  print(table(dom, level))
}
