#!/usr/bin/env Rscript
# =============================================================================
# Does a helminth / NTD-ecology covariate improve AREA-LEVEL iron prediction?
#
# Part A (runs now, real data): use the IHME lymphatic-filariasis (LF) prevalence
#   layer already in the predictor set as a crude proxy for NTD/helminth ecology,
#   and test its association with area-level iron-deficiency prevalence (vs a
#   malaria yardstick). LF is NOT the iron-causal helminth (hookworm is) — this
#   is a hypothesis-motivating signal, not the answer.
#
# Part B (runs when ESPEN data present): merge the real STH/schisto Admin2
#   prevalence produced by scripts/build_espen_admin2.R and compare area-level
#   LOO prediction of iron prevalence: baseline (shared features) vs +helminth.
# =============================================================================
source(here::here("sandbox_fe","harness.R"))

iron_outcomes <- c("child_iron","women_iron")

aggregate_area <- function(o, cols, min_n=5){
  y<-o$y_bin; ok<-!is.na(y); a2<-o$admin2[ok]; y<-y[ok]; M<-o$M[ok,cols,drop=FALSE]
  ua<-unique(a2); keep<-ua[sapply(ua,function(u)sum(a2==u)>=min_n)]
  list(prev=sapply(keep,function(u)mean(y[a2==u])),
       X=t(sapply(keep,function(u){idx<-a2==u; colMeans(M[idx,,drop=FALSE],na.rm=TRUE)})),
       area=keep)
}

cat("================ PART A: LF (NTD-ecology) proxy signal ================\n")
cat("Spearman cor(area iron-deficiency prevalence, covariate), by country + pooled.\n")
cat("LF = lymphatic filariasis prevalence (IHME); malaria = MAP/IHME rate (yardstick).\n\n")

for (otag in iron_outcomes) {
  objs <- load_objs(otag)
  pooled_lf <- c(); pooled_mal <- c(); pooled_prev <- c()
  cat(sprintf("---- %s ----\n", otag))
  for (cn in names(objs)) {
    o <- objs[[cn]]; cols <- colnames(o$M)
    lf_col  <- grep("lf_prevalence_rate|_lf_.*rate", cols, value=TRUE, ignore.case=TRUE)[1]
    mal_cands <- grep("Pf_Parasite_Rate|Pf_Incidence_Rate|malaria_prevalence_rate", cols, value=TRUE, ignore.case=TRUE)
    # prefer a column that actually varies across areas (skip broadcast/constant ones)
    mal_col <- NA
    for (mc in mal_cands) { v<-o$M[,mc]; if (length(unique(round(v[!is.na(v)],6)))>5) { mal_col<-mc; break } }
    use <- na.omit(c(lf_col, mal_col)); if (length(use)==0) { cat(sprintf("  %-12s no LF/malaria cols\n",cn)); next }
    ar <- aggregate_area(o, use)
    sp <- function(v) if (v %in% colnames(ar$X)) round(suppressWarnings(cor(ar$prev, ar$X[,v], method="spearman", use="complete.obs")),2) else NA
    cat(sprintf("  %-12s n_area=%2d  LF=%s  malaria=%s\n", cn, length(ar$prev),
                ifelse(is.na(lf_col),"NA",sprintf("%+.2f",sp(lf_col))),
                ifelse(is.na(mal_col),"NA",sprintf("%+.2f",sp(mal_col)))))
    if (!is.na(lf_col) && lf_col %in% colnames(ar$X)) { pooled_lf<-c(pooled_lf,ar$X[,lf_col]); pooled_prev<-c(pooled_prev,ar$prev) }
    if (!is.na(mal_col) && mal_col %in% colnames(ar$X)) pooled_mal<-c(pooled_mal, ar$X[,mal_col])
  }
  if (length(pooled_lf)>3)
    cat(sprintf("  POOLED       LF=%+.2f  malaria=%+.2f\n",
                suppressWarnings(cor(pooled_prev,pooled_lf,method="spearman",use="complete.obs")),
                if(length(pooled_mal)==length(pooled_prev)) suppressWarnings(cor(pooled_prev,pooled_mal,method="spearman",use="complete.obs")) else NA))
  cat("\n")
}

cat("================ PART B: real ESPEN STH/schisto merge test ================\n")
espen_dir <- here::here("data","ESPEN")
have <- file.exists(file.path(espen_dir, sprintf("%s_espen_admin2.csv",
          c("Gambia","Ghana","SierraLeone","Malawi"))))
if (!any(have)) {
  cat("ESPEN Admin2 files not found. To activate this test:\n")
  cat("  1. Get an API key: https://espen.afro.who.int  (set ESPEN_API_KEY)\n")
  cat("  2. Rscript scripts/build_espen_admin2.R   # writes data/ESPEN/<Country>_espen_admin2.csv\n")
  cat("  3. re-run this script — it will compare baseline vs +helminth area-level prediction.\n")
} else {
  dom_helm <- function(cn) grepl("^helminth_|^espen_", cn)
  for (otag in iron_outcomes) {
    objs <- load_objs(otag)
    shared <- Reduce(intersect, lapply(objs, function(o) colnames(o$M)))
    # merge helminth Admin2 columns onto each country's area aggregates
    base_sp <- c(); helm_sp <- c()
    for (test_cn in names(objs)) {
      # build area data per country incl helminth merge
      mk <- function(cn2) {
        o <- objs[[cn2]]; ar <- aggregate_area(o, shared)
        ef <- file.path(espen_dir, sprintf("%s_espen_admin2.csv", cn2))
        if (file.exists(ef)) {
          e <- read.csv(ef); rownames(e) <- as.character(e$Admin2)
          hv <- grep("^helminth_|^espen_", names(e), value=TRUE)
          H <- as.matrix(e[as.character(ar$area), hv, drop=FALSE])
          list(prev=ar$prev, X=cbind(ar$X, H))
        } else list(prev=ar$prev, X=ar$X)
      }
      parts <- lapply(setdiff(names(objs),test_cn), mk); te <- mk(test_cn)
      common <- Reduce(intersect, c(lapply(parts,function(p)colnames(p$X)), list(colnames(te$X))))
      Xtr <- do.call(rbind, lapply(parts,function(p)p$X[,common,drop=FALSE])); ytr<-unlist(lapply(parts,function(p)p$prev))
      Xte <- te$X[,common,drop=FALSE]; yte<-te$prev
      run <- function(cols){ cc<-.clean(Xtr[,cols,drop=FALSE],Xte[,cols,drop=FALSE]); im<-imp_median(cc$tr,cc$te); z<-tf_zscore(im$tr,im$te)
        pr<-fit_predict_glmnet(as.matrix(z$tr),ytr,as.matrix(z$te),1,"gaussian",fixed_lambda=0.01)
        suppressWarnings(cor(yte,pr,method="spearman",use="complete.obs")) }
      helm_cols <- grep("^helminth_|^espen_", common, value=TRUE)
      base_cols <- setdiff(common, helm_cols)
      base_sp <- c(base_sp, run(base_cols)); helm_sp <- c(helm_sp, run(common))
    }
    cat(sprintf("%-12s area-LOCO Spearman: baseline=%.3f  +helminth=%.3f  (delta %+.3f)\n",
                otag, mean(base_sp,na.rm=TRUE), mean(helm_sp,na.rm=TRUE),
                mean(helm_sp,na.rm=TRUE)-mean(base_sp,na.rm=TRUE)))
  }
}
