#!/usr/bin/env Rscript
# Critically evaluate the metadata available for per-domain PCA.
#  Q1. Is a "domain" (source prefix) a coherent group? -> variance concentration:
#      fraction of variance in top-3 PCs. Low = heterogeneous grab-bag.
#  Q2. Does a finer CONSTRUCT grouping (parsed from names) concentrate better?
#  Q3. Are PC loadings transportable across countries? -> congruence of PC1
#      loading vectors between country pairs (|cor|). Low = non-transportable.
#  Q4. Effective-n / rank: areas per country vs #cols per group.
source(here::here("sandbox_fe","harness.R"))

otag <- "child_vitA"; objs <- load_objs(otag)
shared <- Reduce(intersect, lapply(objs, function(o) colnames(o$M)))

prefix_domain <- function(cn){pref<-c(IHME="^ihme_",MAP="^MAP_",GEE="^gee_",FSEC="^fsec_",SOIL="^soil_");out<-rep("OTHER",length(cn));for(nm in names(pref))out[grepl(pref[[nm]],cn)]<-nm;out}
construct <- function(cn){
  out<-rep("other",length(cn)); g<-function(p) grepl(p,cn,ignore.case=TRUE)
  out[g("trmm|chirps|precip|rainfall")]<-"precip"
  out[g("ndvi|evi|vegetation|gpp|npp|lai")]<-"vegetation"
  out[g("temp|_lst|tavg|tair|temperature")]<-"temperature"
  out[g("elevation|srtm|slope|dem")]<-"elevation"
  out[g("accessib|traveltime|friction")]<-"accessibility"
  out[g("ccnl|nightlight|_ntl|viirs")]<-"nightlight"
  out[g("worldcereal|cropland|cereal|maize|_crop")]<-"cropland"
  out[g("aerosol|atmosphere|aod|co_|no2")]<-"atmosphere"
  out[g("^MAP_")]<-"malaria"
  out[g("soil")]<-"soil"
  out[g("^fsec")]<-"foodsec"
  out[g("^ihme") & out=="other"]<-"ihme_health"
  out
}

# build area-deduplicated, median-imputed, z-scored matrix per country (shared cols)
prep_area <- function(o){
  y<-o$y_bin; ok<-!is.na(y); a2<-o$admin2[ok]; M<-o$M[ok,shared,drop=FALSE]
  ua<-unique(a2); Xa<-t(sapply(ua,function(u){idx<-a2==u;colMeans(M[idx,,drop=FALSE],na.rm=TRUE)}))
  med<-matrixStats::colMedians(Xa,na.rm=TRUE);med[!is.finite(med)]<-0
  for(j in 1:ncol(Xa)){na<-is.na(Xa[,j]);if(any(na))Xa[na,j]<-med[j]}
  mu<-colMeans(Xa);sd<-matrixStats::colSds(Xa);sd[sd<1e-8|!is.finite(sd)]<-1
  sweep(sweep(Xa,2,mu,"-"),2,sd,"/")
}
Zc <- lapply(objs, prep_area)

var_top3 <- function(Z, grp, g){
  idx<-which(grp==g); if(length(idx)<2) return(c(p1=NA,p3=NA,ncol=length(idx)))
  Zi<-Z[,idx,drop=FALSE]; Zi<-Zi[,matrixStats::colSds(Zi)>1e-8,drop=FALSE]
  if(ncol(Zi)<2) return(c(p1=NA,p3=NA,ncol=ncol(Zi)))
  pc<-prcomp(Zi,center=FALSE,scale.=FALSE); v<-pc$sdev^2/sum(pc$sdev^2)
  c(p1=round(v[1],2), p3=round(sum(v[1:min(3,length(v))]),2), ncol=ncol(Zi))
}

cat("=== Q1/Q4: PREFIX-domain coherence (var in top1/top3 PCs) — pooled areas ===\n")
Zpool <- do.call(rbind, Zc)
dgrp <- prefix_domain(shared)
for(g in sort(unique(dgrp))){ s<-var_top3(Zpool,dgrp,g); cat(sprintf("  %-6s ncol=%3d  PC1=%.2f  PC1-3=%.2f\n",g,s["ncol"],s["p1"],s["p3"])) }
cat(sprintf("  [areas: %s]\n", paste(sprintf("%s=%d",names(objs),sapply(Zc,nrow)),collapse=" ")))

cat("\n=== Q2: CONSTRUCT grouping coherence (parsed from names) ===\n")
cgrp <- construct(shared)
cat("  construct sizes:\n"); print(table(cgrp))
for(g in sort(unique(cgrp))){ s<-var_top3(Zpool,cgrp,g); if(!is.na(s["p1"])) cat(sprintf("  %-12s ncol=%3d  PC1=%.2f  PC1-3=%.2f\n",g,s["ncol"],s["p1"],s["p3"])) }

cat("\n=== Q3: cross-country PC1 loading congruence (|cor| of loadings) ===\n")
pc1_load <- function(Z, idx){ Zi<-Z[,idx,drop=FALSE]; sds<-matrixStats::colSds(Zi); keep<-sds>1e-8
  if(sum(keep)<2) return(NULL); pc<-prcomp(Zi[,keep,drop=FALSE],center=FALSE,scale.=FALSE)
  l<-rep(NA,length(idx)); l[keep]<-pc$rotation[,1]; l }
congruence_for <- function(grp,label){
  cat(sprintf("-- %s --\n",label))
  for(g in sort(unique(grp))){
    idx<-which(grp==g); if(length(idx)<3) next
    L<-lapply(Zc,function(Z)pc1_load(Z,idx)); L<-L[!sapply(L,is.null)]
    if(length(L)<2) next
    cn<-names(L); pairs<-combn(length(L),2); cong<-c()
    for(p in 1:ncol(pairs)){a<-L[[pairs[1,p]]];b<-L[[pairs[2,p]]];ok<-!is.na(a)&!is.na(b)
      if(sum(ok)>2) cong<-c(cong, abs(cor(a[ok],b[ok])))}
    cat(sprintf("  %-12s mean|cor(PC1 loadings)|=%.2f  (range %.2f-%.2f)\n",g,mean(cong,na.rm=TRUE),min(cong,na.rm=TRUE),max(cong,na.rm=TRUE)))
  }
}
congruence_for(prefix_domain(shared),"PREFIX domains")
congruence_for(construct(shared),"CONSTRUCT groups")
