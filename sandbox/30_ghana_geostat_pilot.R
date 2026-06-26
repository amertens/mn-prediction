# =============================================================================
# sandbox/30_ghana_geostat_pilot.R
#
# PILOT: cluster-level model-based geostatistics for Ghana (vs the current
# admin-2-aggregated approach). Treats each survey cluster (EA) as a binomial
# observation at its GPS location, fits an INLA-SPDE Matern spatial field, and
# aggregates predictions to admin-2. Honest leave-area-out CV compares:
#   (A) covariate-only  (proxy GLM, no space)   — the current style
#   (B) spatial-only    (SPDE field, no covs)   — pure geostatistics
#   (C) spatial + covs  (full geostatistical)
#
# Outcome: child iron deficiency. Run: Rscript sandbox/30_ghana_geostat_pilot.R
# =============================================================================
suppressWarnings(suppressMessages({
  library(here); library(INLA); library(fmesher); library(dplyr)
  source(here::here("R","config.R")); source(here::here("R","data_prep.R"))
}))
set.seed(1)
cc <- get_country_configs()$Ghana
oc <- cc$outcomes$child_iron
cat(sprintf("Outcome: Ghana child_iron  binary=%s continuous=%s cutoff=%s\n",
            oc$binary, oc$continuous %||% "NA", oc$cutoff %||% "NA"))

g <- load_merged_data(cc$data_path)
d <- g[which(g$gw_child_flag == 1), , drop = FALSE]
d$y <- suppressWarnings(as.numeric(d[[oc$binary]]))
d$lon <- suppressWarnings(as.numeric(d$lon)); d$lat <- suppressWarnings(as.numeric(d$lat))
d <- d[is.finite(d$y) & is.finite(d$lon) & is.finite(d$lat) & !is.na(d[[cc$admin2_col]]), ]
cat(sprintf("Individuals: %d | lon[%.2f,%.2f] lat[%.2f,%.2f]\n",
            nrow(d), min(d$lon), max(d$lon), min(d$lat), max(d$lat)))

# ---- cluster aggregation: binomial counts + mean GPS + admin codes ----------
clus <- d %>%
  dplyr::group_by(cl = .data[[cc$psu_col]]) %>%
  dplyr::summarise(
    n = dplyr::n(), n_def = sum(y == 1, na.rm = TRUE),
    lon = mean(lon), lat = mean(lat),
    Admin2 = dplyr::first(.data[[cc$admin2_col]]),
    Admin1 = dplyr::first(.data[[cc$admin1_col]]),
    .groups = "drop") %>%
  dplyr::mutate(prev = n_def / n)
cat(sprintf("Clusters: %d (median n/cluster=%.0f) | admin2 areas=%d | overall prev=%.1f%%\n",
            nrow(clus), median(clus$n), dplyr::n_distinct(clus$Admin2), 100*mean(d$y)))

# ---- attach area-constant covariates (prescreening happens INSIDE each fold) -
gee <- grep("^gee_a2_", colnames(d), value = TRUE)
covdat <- d %>% dplyr::group_by(cl = .data[[cc$psu_col]]) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(gee), ~ mean(.x, na.rm = TRUE)), .groups="drop")
clus <- dplyr::left_join(clus, covdat, by = "cl")

# Train-only prescreen: top-6 |cor| with prev, decorrelated; returns picked names.
prescreen <- function(tr, k = 6) {
  scr <- sapply(gee, function(v){ x<-tr[[v]]; if(sd(x,na.rm=TRUE)==0||mean(is.finite(x))<.8) return(0)
                                  abs(suppressWarnings(cor(x, tr$prev, use="complete.obs"))) })
  scr <- sort(scr[is.finite(scr) & scr>0], decreasing = TRUE)
  picked <- character()
  for (v in names(scr)) { if (length(picked)>=k) break
    if (length(picked)==0 || max(abs(suppressWarnings(cor(tr[[v]], tr[picked], use="complete.obs"))),na.rm=TRUE) < 0.8)
      picked <- c(picked, v) }
  picked
}
# standardize cols of `te` using `tr` mean/sd (train-only); NA -> 0
std_apply <- function(tr, te, cols) {
  for (v in cols) { mu<-mean(tr[[v]],na.rm=TRUE); s<-sd(tr[[v]],na.rm=TRUE); if(!is.finite(s)||s==0) s<-1
    tr[[v]]<-(tr[[v]]-mu)/s; te[[v]]<-(te[[v]]-mu)/s
    tr[[v]][!is.finite(tr[[v]])]<-0; te[[v]][!is.finite(te[[v]])]<-0 }
  list(tr=tr, te=te)
}
INLA::inla.setOption(num.threads = "1:1")  # limit memory under machine load

# ---- INLA-SPDE mesh + model builders ----------------------------------------
loc  <- cbind(clus$lon, clus$lat)
mesh <- fmesher::fm_mesh_2d_inla(loc, max.edge = c(0.35, 1.2), cutoff = 0.08,
                                 offset = c(0.5, 1.5))
spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(1.0, 0.5), prior.sigma = c(1.0, 0.5))
cat(sprintf("Mesh nodes: %d\n", mesh$n))

fit_inla <- function(train, test, use_space, use_cov, picked) {
  Atr <- INLA::inla.spde.make.A(mesh, cbind(train$lon, train$lat))
  Ate <- INLA::inla.spde.make.A(mesh, cbind(test$lon,  test$lat))
  idx <- INLA::inla.spde.make.index("s", spde$n.spde)
  eff_tr <- c(list(b0 = rep(1, nrow(train))), if (use_cov) as.list(train[picked]))
  eff_te <- c(list(b0 = rep(1, nrow(test))),  if (use_cov) as.list(test[picked]))
  stk_e <- INLA::inla.stack(tag="e", data=list(y=train$n_def, Ntrials=train$n),
              A=list(if(use_space) Atr else 1, 1),
              effects=list(if(use_space) idx else list(idn=rep(1,nrow(train))), eff_tr))
  stk_p <- INLA::inla.stack(tag="p", data=list(y=rep(NA,nrow(test)), Ntrials=test$n),
              A=list(if(use_space) Ate else 1, 1),
              effects=list(if(use_space) idx else list(idn=rep(1,nrow(test))), eff_te))
  stk <- INLA::inla.stack(stk_e, stk_p)
  terms <- c("0","b0", if(use_cov) picked, if(use_space) "f(s, model=spde)")
  form <- stats::as.formula(paste("y ~", paste(terms, collapse=" + ")))
  r <- INLA::inla(form, family="binomial", Ntrials=Ntrials,
                  data=INLA::inla.stack.data(stk),
                  control.predictor=list(A=INLA::inla.stack.A(stk), link=1, compute=TRUE),
                  control.inla=list(int.strategy="eb"), safe=TRUE)
  ip <- INLA::inla.stack.index(stk, "p")$data
  r$summary.fitted.values$mean[ip]
}

# ---- leave-area-out CV (10 folds over admin2), prescreen INSIDE each fold ----
areas <- unique(clus$Admin2); nf <- 10
fold <- setNames(sample(rep_len(1:nf, length(areas))), areas)
clus$fold <- fold[as.character(clus$Admin2)]
specs <- list(A_covonly=c(FALSE,TRUE), B_spaceonly=c(TRUE,FALSE), C_space_cov=c(TRUE,TRUE))
oof <- lapply(names(specs), function(s) { clus$pred <- NA_real_; clus }); names(oof) <- names(specs)
for (k in 1:nf) {
  tr0 <- clus[clus$fold != k, ]; te0 <- clus[clus$fold == k, ]
  if (nrow(te0)==0) next
  picked <- prescreen(tr0)                       # train-only covariate selection
  sa <- std_apply(tr0, te0, picked); tr <- sa$tr; te <- sa$te
  for (s in names(specs)) {
    p <- tryCatch(fit_inla(tr, te, specs[[s]][1], specs[[s]][2], picked),
                  error=function(e){cat("  fit err",s,k,conditionMessage(e),"\n");rep(NA,nrow(te))})
    oof[[s]]$pred[clus$fold==k] <- p
  }
  cat(sprintf("  fold %d/%d done (test clusters=%d, picked=%d)\n", k, nf, nrow(te), length(picked)))
}

# ---- aggregate cluster preds -> admin-2, compare to direct survey prev -------
agg <- function(cl) cl %>% dplyr::group_by(Admin2) %>%
  dplyr::summarise(obs = stats::weighted.mean(prev, n), pred = stats::weighted.mean(pred, n), .groups="drop") %>%
  dplyr::filter(is.finite(pred), is.finite(obs))
aggs <- lapply(oof, agg)
common <- Reduce(intersect, lapply(aggs, function(a) a$Admin2))  # fair: common areas
cat(sprintf("\n==== Leave-area-out CV: admin-2 preds (common %d areas across specs) ====\n", length(common)))
res <- lapply(names(aggs), function(s){ a<-aggs[[s]]; a<-a[a$Admin2 %in% common,]
   data.frame(spec=s, n_area=nrow(a), MAE_pp=round(100*mean(abs(a$pred-a$obs)),2),
   RMSE_pp=round(100*sqrt(mean((a$pred-a$obs)^2)),2),
   r=round(suppressWarnings(cor(a$pred,a$obs)),3)) })
print(do.call(rbind,res), row.names=FALSE)
cat("\n(Compare r to the current admin-2 area-SL honest OOF r ~0.12 for iron.)\n")
saveRDS(list(clus=clus, oof=oof, res=res), here::here("sandbox","_ghana_geostat_pilot.rds"))
