# Build a full-coverage SL -> BYM2 small-area map layer WITH spatially-smoothed
# 95% credible intervals, for every country x outcome.
#
# Motivation (sandbox/sae_sl_hybrid_prototype.R): feeding the area-level
# SuperLearner mean as the covariate into a BYM2 spatial model gives credible
# intervals with near-nominal coverage (~90-93%) at roughly half the width of
# the Fay-Herriot design-based intervals (which over-cover), while leaving the
# point estimates' rank ordering essentially unchanged. This is the "intervals,
# not accuracy" win — BYM2 is the better SAE partner for the dashboard map.
#
# Method:
#   1. Prescreen GEE covariates: greedy top-k by |corr| with the direct
#      estimate, dropping ones collinear with an already-chosen covariate
#      (same routine the FH layer uses) -> well-conditioned, fast SL.
#   2. Area-level SuperLearner (compact glmnet/ranger library, MSE loss) on the
#      surveyed districts. mu_i = SL synthetic mean for EVERY district:
#        - surveyed districts: 5-fold OUT-OF-FOLD SL prediction (no leakage),
#        - unsurveyed districts: SL fit on all surveyed, predicted.
#   3. BYM2 via INLA over ALL districts (graph from the admin-2 polygons):
#        direct_i ~ mu_i + f(area, bym2, graph), Gaussian likelihood with the
#        per-area sampling variance fixed via `scale` (= 1 / D_i). Unsurveyed
#        districts carry NA outcome -> predicted from covariate + spatial field.
#   4. Posterior fitted mean + 2.5%/97.5% quantiles -> pred_prev, ci_lo, ci_hi.
#
# Output schema matches admin2_fh_predictions so the map helpers consume it
# unchanged: country, outcome, Admin2, pred_prev, obs_prev, ci_lo, ci_hi,
# ci_width, n_survey, who_class.
suppressWarnings(suppressMessages({
  library(dplyr); library(sf); library(spdep)
}))
ONE  <- Sys.getenv("BYM2_ONE", "")
root <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
if (!dir.exists(file.path(root, "_targets_full"))) root <- getwd()
S    <- file.path(root, "_targets_full", "objects")
DASH <- file.path(root, "dashboard", "data")
read_prod <- function(nm){ p<-file.path(S,nm); if(file.exists(p)) tryCatch(readRDS(p),error=function(e)NULL) else NULL }

countries <- c(gambia="Gambia", ghana="Ghana", sierraleone="Sierra Leone", malawi="Malawi")
outcome_tags <- c("child_vitA","women_vitA","child_iron","women_iron",
                  "women_folate","women_b12","child_zinc","women_zinc")
who_thresholds <- list(
  child_vitA=c(none=0.02,mild=0.10,moderate=0.20), women_vitA=c(none=0.02,mild=0.10,moderate=0.20),
  child_iron=c(none=0.05,mild=0.20,moderate=0.40), women_iron=c(none=0.05,mild=0.20,moderate=0.40),
  women_folate=c(none=0.05,mild=0.20,moderate=0.40), women_b12=c(none=0.05,mild=0.20,moderate=0.40),
  child_zinc=c(none=0.10,mild=0.20,moderate=0.30), women_zinc=c(none=0.10,mild=0.20,moderate=0.30))
classify_who <- function(prev, oc){ th<-who_thresholds[[oc]]; if(is.null(th)||is.na(prev)) return(NA_character_)
  if(prev<th["none"])return("Low"); if(prev<th["mild"])return("Mild"); if(prev<th["moderate"])return("Moderate"); "Severe" }
clip01 <- function(x) pmin(pmax(x,0),1)

boundaries <- tryCatch(readRDS(file.path(DASH, "admin2_boundaries.rds")), error=function(e) NULL)

# Greedy top-k by |cor with Y|, dropping covariates collinear with chosen ones.
greedy_decorr <- function(train, cand, Y, k = 15L, cor_max = 0.85) {
  sc <- vapply(cand, function(v){ x<-train[[v]]
    if (sum(is.finite(x))<3 || length(unique(x[is.finite(x)]))<2) return(0)
    r <- suppressWarnings(stats::cor(x, Y, use="pairwise.complete.obs")); if (is.na(r)) 0 else abs(r) },
    numeric(1))
  ord <- names(sort(sc, decreasing = TRUE)); ord <- ord[sc[ord] > 0]
  chosen <- character(0)
  for (v in ord) {
    if (length(chosen) >= k) break
    ok <- length(chosen) == 0 || all(vapply(chosen, function(c){
      rr <- suppressWarnings(abs(stats::cor(train[[v]], train[[c]], use="pairwise.complete.obs")))
      is.na(rr) || rr < cor_max }, logical(1)))
    if (ok) chosen <- c(chosen, v)
  }
  chosen
}

SL_LIB <- list(
  list("glmnet", alpha = 1,   id = "lasso"),
  list("glmnet", alpha = 0,   id = "ridge"),
  list("glmnet", alpha = 0.5, id = "elastic_net"),
  list("ranger", num.trees = 400, min.node.size = 5, id = "ranger")
)

# median-impute a frame's `vars` using train medians
impute_with <- function(target, train, vars) {
  for (v in vars) {
    x <- train[[v]]; if (inherits(x,"haven_labelled")) x <- as.double(unclass(x))
    med <- stats::median(x[is.finite(x)], na.rm = TRUE); if (!is.finite(med)) med <- 0
    tx <- target[[v]]; if (inherits(tx,"haven_labelled")) tx <- as.double(unclass(tx))
    tx[!is.finite(tx)] <- med; target[[v]] <- tx
  }
  target
}

# SL synthetic mean for ALL districts: OOF on surveyed, full-fit predict on rest
sl_synthetic <- function(train, all_df, vars, kfold = 5L, seed = 20260615) {
  ok <- requireNamespace("mlr3superlearner", quietly = TRUE) &&
        requireNamespace("mlr3extralearners", quietly = TRUE)
  if (!ok) return(NULL)
  suppressWarnings(suppressMessages(library(mlr3extralearners)))
  fit_sl <- function(tr) tryCatch(mlr3superlearner::mlr3superlearner(
    data = data.frame(Y = tr$svy_prev, tr[, vars, drop = FALSE]),
    target = "Y", library = SL_LIB, outcome_type = "continuous",
    folds = min(5L, max(3L, nrow(tr) - 1L))), error = function(e) NULL)
  pred_sl <- function(fit, nd) tryCatch(as.numeric(predict(
    fit, data.frame(nd[, vars, drop = FALSE]))), error = function(e) NULL)

  # full fit on all surveyed -> predict every district (covariate for unsurveyed)
  tr_full <- impute_with(train, train, vars)
  full    <- fit_sl(tr_full); if (is.null(full)) return(NULL)
  mu_all  <- pred_sl(full, impute_with(all_df, train, vars))
  if (is.null(mu_all)) return(NULL)
  names(mu_all) <- all_df$Admin2

  # OOF for surveyed districts
  set.seed(seed); n <- nrow(train); fold <- sample(rep_len(seq_len(kfold), n))
  oof <- rep(NA_real_, n)
  for (k in seq_len(kfold)) {
    tr <- which(fold != k); te <- which(fold == k)
    if (length(tr) < 4) next
    f <- fit_sl(impute_with(train[tr,,drop=FALSE], train[tr,,drop=FALSE], vars))
    if (is.null(f)) next
    p <- pred_sl(f, impute_with(train[te,,drop=FALSE], train[tr,,drop=FALSE], vars))
    if (!is.null(p)) oof[te] <- p
  }
  mu_all[match(train$Admin2[is.finite(oof)], names(mu_all))] <- oof[is.finite(oof)]
  clip01(as.numeric(mu_all))
}

fit_bym2_country_outcome <- function(low, oc, deff = 1.5) {
  if (!requireNamespace("INLA", quietly = TRUE)) return(NULL)
  svy <- read_prod(paste0("svy_admin2_", low, "_", oc))
  gee <- read_prod(paste0("gee_admin2_", low))
  if (is.null(svy) || is.null(gee) || !"svy_prev" %in% colnames(svy)) return(NULL)
  gee_vars <- grep("^gee_", colnames(gee), value = TRUE)
  if (length(gee_vars) < 2) return(NULL)
  all_df <- gee[, c("Admin2", gee_vars), drop = FALSE]

  bnd <- boundaries[[low]]
  if (is.null(bnd)) return(NULL)
  # align polygons to gee row order; keep only districts present in both
  all_df <- all_df[all_df$Admin2 %in% bnd$Admin2, , drop = FALSE]
  bnd <- bnd[match(all_df$Admin2, bnd$Admin2), ]
  if (nrow(all_df) < 8) return(NULL)

  surv <- svy %>% dplyr::select(Admin2, svy_prev, dplyr::any_of("svy_prev_se"), n_svy) %>%
    dplyr::filter(is.finite(svy_prev))
  surv <- surv[surv$Admin2 %in% all_df$Admin2, , drop = FALSE]
  if (nrow(surv) < 6) return(NULL)
  train <- dplyr::inner_join(all_df, surv, by = "Admin2")

  vars <- greedy_decorr(train, gee_vars, train$svy_prev, k = 15L)
  if (length(vars) < 1) return(NULL)

  mu_all <- sl_synthetic(train, all_df, vars)
  if (is.null(mu_all) || all(!is.finite(mu_all))) return(NULL)

  # spatial graph over ALL districts (in all_df / mu_all order)
  nb <- tryCatch(spdep::poly2nb(bnd, queen = TRUE), error = function(e) NULL)
  if (is.null(nb)) return(NULL)
  W  <- tryCatch(spdep::nb2mat(nb, style = "B", zero.policy = TRUE), error=function(e) NULL)
  g  <- tryCatch(INLA::inla.read.graph(W), error = function(e) NULL)
  if (is.null(g) || g$n != nrow(all_df)) return(NULL)

  # design-based sampling variance D_i for surveyed districts
  mi <- match(train$Admin2, all_df$Admin2)
  p  <- pmin(pmax(train$svy_prev,1e-4),1-1e-4)
  D  <- if ("svy_prev_se" %in% names(train)) train$svy_prev_se^2 else rep(NA_real_, nrow(train))
  bad <- !is.finite(D) | D <= 0
  D[bad] <- p[bad]*(1-p[bad])/pmax(train$n_svy[bad]/deff, 1); D <- pmax(D, 1e-5)

  y     <- rep(NA_real_, nrow(all_df)); y[mi] <- train$svy_prev
  scale <- rep(1, nrow(all_df));        scale[mi] <- 1 / D     # precision multiplier
  dat <- data.frame(y = y, mu = mu_all, area = seq_len(nrow(all_df)), scale = scale)

  fit <- tryCatch(INLA::inla(
    y ~ mu + f(area, model = "bym2", graph = g, scale.model = TRUE,
               hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                            phi  = list(prior = "pc",      param = c(0.5, 0.5)))),
    family = "gaussian", data = dat, scale = dat$scale,
    control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
    control.predictor = list(compute = TRUE),
    control.compute = list(return.marginals.predictor = FALSE)),
    error = function(e) { cat("    [bym2] INLA failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)

  fv  <- fit$summary.fitted.values
  est <- clip01(fv$mean)
  lo  <- clip01(fv$`0.025quant`); hi <- clip01(fv$`0.975quant`)
  obs <- rep(NA_real_, nrow(all_df)); obs[mi] <- train$svy_prev
  nsv <- rep(NA_integer_, nrow(all_df)); nsv[mi] <- train$n_svy
  data.frame(country=unname(countries[low]), outcome=oc, Admin2=all_df$Admin2,
             pred_prev=est, obs_prev=obs, ci_lo=lo, ci_hi=hi, ci_width=hi-lo,
             n_survey=nsv, who_class=vapply(est, classify_who, character(1), oc=oc),
             stringsAsFactors=FALSE)
}

rows <- list(); ok <- 0; fail <- character(0)
combos <- if (nzchar(ONE)) {
  list(c(sub("_.*","",ONE), sub("^[^_]*_","",ONE)))
} else {
  cc<-list(); for(l in names(countries)) for(o in outcome_tags) cc[[length(cc)+1]]<-c(l,o); cc
}
for (z in combos) {
  low <- z[1]; oc <- z[2]
  r <- tryCatch(fit_bym2_country_outcome(low, oc), error=function(e){cat("  ERR",low,oc,e$message,"\n");NULL})
  if (!is.null(r) && nrow(r)>0) { rows[[paste(low,oc)]]<-r; ok<-ok+1
    cat(sprintf("  BYM2 ok %-24s districts=%d surveyed=%d mean CI(pp): surv=%.1f unsurv=%.1f\n",
        paste(low,oc), nrow(r), sum(!is.na(r$obs_prev)),
        mean(r$ci_width[!is.na(r$obs_prev)],na.rm=TRUE)*100,
        mean(r$ci_width[is.na(r$obs_prev)],na.rm=TRUE)*100))
  } else fail <- c(fail, paste(low,oc))
}
out <- if (length(rows)) do.call(rbind, rows) else NULL
if (!nzchar(ONE) && !is.null(out)) {
  saveRDS(out, file.path(DASH, "admin2_bym2_predictions.rds"))
  cat(sprintf("\nadmin2_bym2_predictions.rds: %d rows, %d combos OK\n", nrow(out), ok))
  if (length(fail)) cat("  no BYM2 layer for:", paste(fail, collapse=", "), "\n")
} else if (nzchar(ONE)) {
  cat("\n[ONE-slice test]\n")
  print(utils::head(out[order(-out$ci_width), c("Admin2","pred_prev","ci_lo","ci_hi","obs_prev")], 5))
}
