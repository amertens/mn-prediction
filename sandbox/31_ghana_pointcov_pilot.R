# =============================================================================
# sandbox/31_ghana_pointcov_pilot.R
#
# Follow-on to sandbox/30: re-run the Ghana cluster-level geostatistics pilot
# with TRUE point-level covariates (gee_cl_*, extracted at 2km/5km buffers
# around each cluster's GPS) instead of the admin-2-constant gee_a2_*.
# Head-to-head leave-area-out CV: does point extraction beat area-constant covs?
#
# Run AFTER the regen finishes (terra extraction is IO/CPU heavy).
#   Rscript sandbox/31_ghana_pointcov_pilot.R
# =============================================================================
suppressWarnings(suppressMessages({
  library(here); library(INLA); library(fmesher); library(dplyr); library(sf); library(terra)
  source(here::here("R","config.R")); source(here::here("R","data_prep.R"))
  source(here::here("src","GEE","extract_gee_admin2.R"))
  source(here::here("src","GEE","extract_gee_cluster.R"))
}))
set.seed(1)
cc <- get_country_configs()$Ghana; oc <- cc$outcomes$child_iron
g  <- load_merged_data(cc$data_path)
d  <- g[which(g$gw_child_flag == 1), , drop = FALSE]
d$y <- suppressWarnings(as.numeric(d[[oc$binary]]))
d$lon <- suppressWarnings(as.numeric(d$lon)); d$lat <- suppressWarnings(as.numeric(d$lat))
d <- d[is.finite(d$y) & is.finite(d$lon) & is.finite(d$lat) & !is.na(d[[cc$admin2_col]]), ]

# ---- cluster table -----------------------------------------------------------
clus <- d %>% dplyr::group_by(cl = .data[[cc$psu_col]]) %>%
  dplyr::summarise(n=dplyr::n(), n_def=sum(y==1,na.rm=TRUE), lon=mean(lon), lat=mean(lat),
    Admin2=dplyr::first(.data[[cc$admin2_col]]), Admin1=dplyr::first(.data[[cc$admin1_col]]),
    urban=if ("gw_urban/ rural" %in% names(d)) dplyr::first(.data[["gw_urban/ rural"]]) else NA_character_,
    .groups="drop") %>%
  dplyr::mutate(prev=n_def/n)
cat(sprintf("Ghana child_iron: %d clusters, %d areas, prev=%.1f%%\n", nrow(clus),
            dplyr::n_distinct(clus$Admin2), 100*mean(d$y)))

# ---- POINT-LEVEL covariate extraction (gee_cl_) ------------------------------
clus_sf <- sf::st_as_sf(clus, coords=c("lon","lat"), crs=4326, remove=FALSE)
rasters <- list.files(here::here("data","Ghana_GEE rasters"), pattern="\\.tif$", full.names=TRUE)
rasters <- c(rasters, list.files(here::here("data","GEE"), pattern="\\.tif$", full.names=TRUE))
cat(sprintf("Extracting %d rasters at cluster buffers (2km urban / 5km rural)...\n", length(rasters)))
clcov <- extract_gee_cluster(clus_sf, rasters, id_col="cl",
                             buffer_urban=2000, buffer_rural=5000, buffer_m=4000,
                             urban_col=if (!all(is.na(clus$urban))) "urban" else NULL,
                             country_tokens=c("Ghana"), verbose=FALSE)
gee_cl <- grep("^gee_cl_", names(clcov), value=TRUE)
clus <- dplyr::left_join(clus, clcov[c("cl", gee_cl)], by="cl")
cat(sprintf("  extracted %d gee_cl_ columns\n", length(gee_cl)))

# sanity: points inside their admin-2? + does gee_cl vary WITHIN area (gee_a2 doesn't)?
gee_a2 <- grep("^gee_a2_", colnames(d), value=TRUE)
a2 <- d %>% dplyr::group_by(cl=.data[[cc$psu_col]]) %>% dplyr::summarise(dplyr::across(dplyr::all_of(gee_a2), ~mean(.x,na.rm=TRUE)), .groups="drop")
clus <- dplyr::left_join(clus, a2, by="cl")
within_var <- function(cols){ v<-sapply(cols, function(c) mean(tapply(clus[[c]], clus$Admin2, function(z) sd(z,na.rm=TRUE)),na.rm=TRUE)); mean(v,na.rm=TRUE) }
cat(sprintf("  mean within-area SD: gee_cl_=%.4g vs gee_a2_=%.4g (gee_a2_ should be ~0)\n",
            within_var(head(gee_cl,50)), within_var(head(gee_a2,50))))

# ---- shared CV machinery -----------------------------------------------------
prescreen <- function(tr, pool, k=6){ scr<-sapply(pool,function(v){x<-tr[[v]]; if(is.null(x)||sd(x,na.rm=TRUE)==0||mean(is.finite(x))<.8)return(0); abs(suppressWarnings(cor(x,tr$prev,use="complete.obs")))})
  scr<-sort(scr[is.finite(scr)&scr>0],decreasing=TRUE); pk<-character()
  for(v in names(scr)){ if(length(pk)>=k)break; if(length(pk)==0||max(abs(suppressWarnings(cor(tr[[v]],tr[pk],use="complete.obs"))),na.rm=TRUE)<0.8) pk<-c(pk,v)}; pk }
std_apply <- function(tr,te,cols){ for(v in cols){mu<-mean(tr[[v]],na.rm=TRUE);s<-sd(tr[[v]],na.rm=TRUE);if(!is.finite(s)||s==0)s<-1
  tr[[v]]<-(tr[[v]]-mu)/s;te[[v]]<-(te[[v]]-mu)/s;tr[[v]][!is.finite(tr[[v]])]<-0;te[[v]][!is.finite(te[[v]])]<-0}; list(tr=tr,te=te) }
loc<-cbind(clus$lon,clus$lat); mesh<-fmesher::fm_mesh_2d_inla(loc,max.edge=c(0.35,1.2),cutoff=0.08,offset=c(0.5,1.5))
spde<-INLA::inla.spde2.pcmatern(mesh,prior.range=c(1,0.5),prior.sigma=c(1,0.5)); INLA::inla.setOption(num.threads="1:1")
fit_inla<-function(tr,te,use_space,use_cov,pk){ Atr<-INLA::inla.spde.make.A(mesh,cbind(tr$lon,tr$lat)); Ate<-INLA::inla.spde.make.A(mesh,cbind(te$lon,te$lat)); idx<-INLA::inla.spde.make.index("s",spde$n.spde)
  et<-c(list(b0=rep(1,nrow(tr))),if(use_cov)as.list(tr[pk])); ee<-c(list(b0=rep(1,nrow(te))),if(use_cov)as.list(te[pk]))
  se<-INLA::inla.stack(tag="e",data=list(y=tr$n_def,Ntrials=tr$n),A=list(if(use_space)Atr else 1,1),effects=list(if(use_space)idx else list(idn=rep(1,nrow(tr))),et))
  sp<-INLA::inla.stack(tag="p",data=list(y=rep(NA,nrow(te)),Ntrials=te$n),A=list(if(use_space)Ate else 1,1),effects=list(if(use_space)idx else list(idn=rep(1,nrow(te))),ee))
  st<-INLA::inla.stack(se,sp); tm<-c("0","b0",if(use_cov)pk,if(use_space)"f(s,model=spde)")
  r<-INLA::inla(stats::as.formula(paste("y ~",paste(tm,collapse=" + "))),family="binomial",Ntrials=Ntrials,data=INLA::inla.stack.data(st),
     control.predictor=list(A=INLA::inla.stack.A(st),link=1,compute=TRUE),control.inla=list(int.strategy="eb"),safe=TRUE)
  r$summary.fitted.values$mean[INLA::inla.stack.index(st,"p")$data] }

areas<-unique(clus$Admin2); nf<-10; fold<-setNames(sample(rep_len(1:nf,length(areas))),areas); clus$fold<-fold[as.character(clus$Admin2)]
# specs: covariate pool x model
specs<-list(A2_cov_area=list(pool="a2",sp=FALSE,cov=TRUE), CL_cov_point=list(pool="cl",sp=FALSE,cov=TRUE),
            A2_spacecov_area=list(pool="a2",sp=TRUE,cov=TRUE), CL_spacecov_point=list(pool="cl",sp=TRUE,cov=TRUE))
pools<-list(a2=gee_a2, cl=gee_cl)
oof<-lapply(specs,function(s){clus$pred<-NA_real_;clus}); names(oof)<-names(specs)
for(k in 1:nf){ tr0<-clus[clus$fold!=k,]; te0<-clus[clus$fold==k,]; if(nrow(te0)==0)next
  for(s in names(specs)){ sp<-specs[[s]]; pk<-prescreen(tr0,pools[[sp$pool]]); sa<-std_apply(tr0,te0,pk)
    p<-tryCatch(fit_inla(sa$tr,sa$te,sp$sp,sp$cov,pk),error=function(e){cat("err",s,k,conditionMessage(e),"\n");rep(NA,nrow(te0))})
    oof[[s]]$pred[clus$fold==k]<-p }
  cat(sprintf("  fold %d/%d done\n",k,nf)) }
agg<-function(cl) cl%>%dplyr::group_by(Admin2)%>%dplyr::summarise(obs=stats::weighted.mean(prev,n),pred=stats::weighted.mean(pred,n),.groups="drop")%>%dplyr::filter(is.finite(pred),is.finite(obs))
aggs<-lapply(oof,agg); common<-Reduce(intersect,lapply(aggs,function(a)a$Admin2))
cat(sprintf("\n==== Leave-area-out CV (common %d areas): area-constant vs POINT covariates ====\n",length(common)))
res<-lapply(names(aggs),function(s){a<-aggs[[s]];a<-a[a$Admin2%in%common,]; data.frame(spec=s,MAE_pp=round(100*mean(abs(a$pred-a$obs)),2),r=round(suppressWarnings(cor(a$pred,a$obs)),3))})
print(do.call(rbind,res),row.names=FALSE)
cat("\n(Compare to sandbox/30 gee_a2_ results: A=0.571, C=0.605.)\n")

# ---- buffer sensitivity sweep (point covariates, space+cov spec) ------------
cat("\n==== Buffer sensitivity (CL space+cov, leave-area-out, common areas) ====\n")
sweep <- lapply(c(1000, 2000, 5000, 10000), function(R) {
  cv   <- extract_gee_cluster(clus_sf, rasters, id_col = "cl", buffer_m = R,
                              urban_col = NULL, country_tokens = c("Ghana"), verbose = FALSE)
  cols <- grep("^gee_cl_", names(cv), value = TRUE)
  cl2  <- clus; cl2[intersect(cols, names(cl2))] <- NULL
  cl2  <- dplyr::left_join(cl2, cv[c("cl", cols)], by = "cl"); cl2$pred <- NA_real_
  for (k in 1:nf) {
    tr0 <- cl2[cl2$fold != k, ]; te0 <- cl2[cl2$fold == k, ]; if (nrow(te0) == 0) next
    pk <- prescreen(tr0, cols); sa <- std_apply(tr0, te0, pk)
    cl2$pred[cl2$fold == k] <- tryCatch(fit_inla(sa$tr, sa$te, TRUE, TRUE, pk),
                                        error = function(e) rep(NA, nrow(te0)))
  }
  a <- agg(cl2); a <- a[a$Admin2 %in% common, ]
  data.frame(buffer = paste0(R/1000, "km"), MAE_pp = round(100*mean(abs(a$pred-a$obs)), 2),
             r = round(suppressWarnings(cor(a$pred, a$obs)), 3))
})
print(do.call(rbind, sweep), row.names = FALSE)
cat("(2/5km urban/rural split result is the CL_spacecov_point row above.)\n")
saveRDS(list(clus=clus, res=res, sweep=sweep), here::here("sandbox","_ghana_pointcov_pilot.rds"))
