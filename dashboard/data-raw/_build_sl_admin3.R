# Build Sierra Leone Admin-3 (chiefdom) prevalence layer for the map.
# Reconstructs chiefdom prevalence from individuals (cluster -> Admin-3 map +
# survey weights + the binary deficiency flag), then fits a Fay-Herriot /
# empirical-Bayes area model on the chiefdom GEE covariates to cover ALL 153
# chiefdoms with a point estimate + 95% interval. Also fetches GADM level-3
# boundaries. Outputs dashboard/data/admin3_predictions.rds + admin3_boundaries.rds.
suppressWarnings(suppressMessages({ library(dplyr); library(sf) }))
root <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
if (!dir.exists(file.path(root, "_targets_full"))) root <- getwd()
S    <- file.path(root, "_targets_full", "objects")
TR   <- file.path(root, "results", "transportability")
DASH <- file.path(root, "dashboard", "data")
read_prod <- function(nm){ p<-file.path(S,nm); if(file.exists(p)) readRDS(p) else NULL }
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path(root, "R", "config.R"))
clip01 <- function(x) pmin(pmax(x,0),1)

who_thresholds <- list(
  child_vitA=c(none=0.02,mild=0.10,moderate=0.20), women_vitA=c(none=0.02,mild=0.10,moderate=0.20),
  child_iron=c(none=0.05,mild=0.20,moderate=0.40), women_iron=c(none=0.05,mild=0.20,moderate=0.40),
  women_folate=c(none=0.05,mild=0.20,moderate=0.40), women_b12=c(none=0.05,mild=0.20,moderate=0.40))
classify_who <- function(prev, oc){ th<-who_thresholds[[oc]]; if(is.null(th)||is.na(prev)) return(NA_character_)
  if(prev<th["none"])return("Low"); if(prev<th["mild"])return("Mild"); if(prev<th["moderate"])return("Moderate"); "Severe" }

greedy_decorr <- function(train, cand, Y, k=4L, cor_max=0.75){
  sc <- vapply(cand, function(v){x<-train[[v]]
    if(sum(is.finite(x))<3||length(unique(x[is.finite(x)]))<2) return(0)
    r<-suppressWarnings(stats::cor(x,Y,use="pairwise.complete.obs")); if(is.na(r))0 else abs(r)}, numeric(1))
  ord <- names(sort(sc, decreasing=TRUE)); ord <- ord[sc[ord]>0]; chosen <- character(0)
  for(v in ord){ if(length(chosen)>=k) break
    ok <- length(chosen)==0 || all(vapply(chosen, function(c){
      rr<-suppressWarnings(abs(stats::cor(train[[v]],train[[c]],use="pairwise.complete.obs"))); is.na(rr)||rr<cor_max}, logical(1)))
    if(ok) chosen<-c(chosen,v) }
  chosen
}

cfg <- get_country_configs()[["SierraLeone"]]
c2a <- readRDS(file.path(TR, "sierraleone_cluster_to_admin3.rds"))   # cluster -> NAME_3
gee <- readRDS(file.path(TR, "sierraleone_admin3_gee.rds"))          # 153 chiefdoms
gee_vars <- grep("^gee_", colnames(gee), value = TRUE)
cluster_id <- cfg$cluster_id %||% "gw_cnum"
weight_col <- cfg$weight_col %||% "gw_svy_weight"

rows <- list()
for (oc_name in names(cfg$outcomes)) {
  oc <- cfg$outcomes[[oc_name]]
  od <- read_prod(paste0("outcome_data_sierraleone_", oc_name))$data
  yb <- oc$binary
  if (is.null(od) || is.null(yb) || !yb %in% colnames(od)) next
  d <- data.frame(y = as.numeric(od[[yb]]),
                  w = if (weight_col %in% colnames(od)) od[[weight_col]] else 1,
                  cl = od[[cluster_id]], stringsAsFactors = FALSE)
  d <- d[is.finite(d$y), , drop = FALSE]
  d$NAME_3 <- c2a$NAME_3[match(d$cl, c2a$cluster)]
  d <- d[!is.na(d$NAME_3), , drop = FALSE]
  if (nrow(d) == 0) next
  agg <- d %>% group_by(NAME_3) %>%
    summarise(svy_prev = sum(w * y, na.rm=TRUE) / sum(w, na.rm=TRUE),
              n_svy = dplyr::n(), .groups = "drop") %>%
    filter(is.finite(svy_prev))
  if (nrow(agg) < 8) { cat(sprintf("  skip %s (only %d surveyed chiefdoms)\n", oc_name, nrow(agg))); next }

  train <- dplyr::inner_join(gee[, c("NAME_3", gee_vars)], agg, by = "NAME_3")
  vars <- greedy_decorr(train, gee_vars, train$svy_prev, k = 4L)
  if (length(vars) < 1) next
  med <- vapply(vars, function(v) stats::median(train[[v]], na.rm=TRUE), numeric(1))
  ctr <- vapply(vars, function(v){x<-train[[v]];x[!is.finite(x)]<-med[[v]];mean(x)}, numeric(1))
  scl <- vapply(vars, function(v){x<-train[[v]];x[!is.finite(x)]<-med[[v]];s<-stats::sd(x);if(!is.finite(s)||s<1e-9)1 else s}, numeric(1))
  zmat <- function(dd){M<-sapply(vars,function(v){x<-dd[[v]];x[!is.finite(x)]<-med[[v]];(x-ctr[[v]])/scl[[v]]})
                       z<-as.data.frame(M); colnames(z)<-paste0("z",seq_along(vars)); z}
  Ztr <- zmat(train); Zall <- zmat(gee)
  p <- pmin(pmax(train$svy_prev,1e-4),1-1e-4); D <- p*(1-p)/pmax(train$n_svy/1.5,1); D <- pmax(D,1e-5)
  fit <- stats::lm(stats::reformulate(colnames(Ztr),"y"), data = cbind(y=train$svy_prev, Ztr))
  n <- nrow(Ztr); pp <- sum(!is.na(stats::coef(fit)))
  A <- max(sum(stats::residuals(fit)^2)/max(n-pp,1) - mean(D), 1e-4)
  prr <- stats::predict(fit, newdata = Zall, se.fit = TRUE)
  est <- clip01(as.numeric(prr$fit)); se <- sqrt(A + as.numeric(prr$se.fit)^2)
  mi <- match(train$NAME_3, gee$NAME_3); gamma <- A/(A+D)
  est[mi] <- clip01(gamma*train$svy_prev + (1-gamma)*as.numeric(stats::predict(fit)))
  se[mi]  <- sqrt(gamma*D)
  obs <- rep(NA_real_, nrow(gee)); obs[mi] <- train$svy_prev
  nsv <- rep(NA_integer_, nrow(gee)); nsv[mi] <- train$n_svy
  rows[[oc_name]] <- data.frame(
    country = "Sierra Leone", outcome = oc_name,
    Admin1 = gee$Admin1, Admin2 = gee$NAME_2, Admin3 = gee$NAME_3,
    pred_prev = est, obs_prev = obs,
    ci_lo = clip01(est-1.96*se), ci_hi = clip01(est+1.96*se), ci_width = clip01(est+1.96*se)-clip01(est-1.96*se),
    n_survey = nsv, who_class = vapply(est, classify_who, character(1), oc=oc_name),
    stringsAsFactors = FALSE)
  cat(sprintf("  %s: %d chiefdoms, %d surveyed\n", oc_name, nrow(gee), nrow(agg)))
}
admin3_pred <- bind_rows(rows)
saveRDS(admin3_pred, file.path(DASH, "admin3_predictions.rds"))
cat(sprintf("admin3_predictions.rds: %d rows, outcomes: %s\n", nrow(admin3_pred),
            paste(unique(admin3_pred$outcome), collapse=", ")))

# ── GADM level-3 boundaries (chiefdoms) ─────────────────────────────────────
if (requireNamespace("geodata", quietly = TRUE)) {
  cache <- file.path(root, "data", "external_cache", "gadm")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  g3 <- tryCatch(geodata::gadm("SLE", level = 3, path = cache), error = function(e) NULL)
  if (!is.null(g3)) {
    s3 <- sf::st_as_sf(g3)
    s3$Admin3 <- s3$NAME_3; s3$Admin2 <- s3$NAME_2; s3$Admin1 <- s3$NAME_1
    s3 <- s3[, c("Admin1","Admin2","Admin3","geometry")]
    s3$country <- "Sierra Leone"
    s3 <- sf::st_make_valid(sf::st_simplify(s3, dTolerance = 0.001, preserveTopology = TRUE))
    saveRDS(list(sierraleone = s3), file.path(DASH, "admin3_boundaries.rds"))
    cat(sprintf("admin3_boundaries.rds: %d chiefdom polygons\n", nrow(s3)))
  } else cat("GADM level-3 download failed — admin3_boundaries not written.\n")
} else cat("geodata not available — admin3_boundaries not written.\n")
