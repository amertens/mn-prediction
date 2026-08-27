# Build a full-coverage Fay-Herriot / empirical-Bayes small-area map layer WITH
# design-based intervals, for every country x outcome.
#
# Method (area-level random-effects / Fay-Herriot, computed directly so it is
# robust to the near-collinear GEE covariates that make sae::eblupFH singular):
#   1. Pick the top-k GEE covariates by |corr| with the direct estimate, greedily
#      DROPPING ones highly correlated with an already-chosen covariate, then
#      standardize them (training mean/sd) -> well-conditioned design.
#   2. Area-level regression of the survey-weighted ("direct") prevalence on
#      those covariates (the synthetic predictor mu = X beta).
#   3. Fay-Herriot moment estimate of the between-area model variance A
#      (residual variance net of mean sampling variance).
#   4. Surveyed district i: EBLUP = gamma_i * direct_i + (1-gamma_i) * mu_i, with
#      gamma_i = A/(A + D_i) and design sampling variance D_i (from the survey
#      SE, falling back to p(1-p)/n_eff). 95% interval from the EB posterior SD
#      sqrt(gamma_i * D_i).
#   5. Unsurveyed district j: synthetic mu_j, 95% interval from sqrt(A + reg.SE^2)
#      (no local survey -> wider, model-based interval).
# Output schema matches admin2_predictions so the map helpers can consume it.
suppressWarnings(suppressMessages(library(dplyr)))
ONE  <- Sys.getenv("FH_ONE", "")
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

# Greedy top-k by |cor with Y|, dropping covariates collinear with chosen ones.
greedy_decorr <- function(train, cand, Y, k = 3L, cor_max = 0.75) {
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

fit_fh_country_outcome <- function(low, oc, k = 3L, deff = 1.5) {
  svy <- read_prod(paste0("svy_admin2_", low, "_", oc))
  gee <- read_prod(paste0("gee_admin2_", low))
  if (is.null(svy) || is.null(gee) || !"svy_prev" %in% colnames(svy)) return(NULL)
  gee_vars <- grep("^gee_", colnames(gee), value = TRUE)
  if (length(gee_vars) < 2) return(NULL)
  all_df <- gee[, c("Admin2", gee_vars), drop = FALSE]
  surv <- svy %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy) %>%
    dplyr::filter(is.finite(svy_prev))
  if (nrow(surv) < 6) return(NULL)
  train <- dplyr::inner_join(all_df, surv, by = "Admin2")
  vars <- greedy_decorr(train, gee_vars, train$svy_prev, k = k)
  if (length(vars) < 1) return(NULL)

  # standardize covariates on TRAIN stats (impute with train medians first)
  med <- vapply(vars, function(v) stats::median(train[[v]], na.rm=TRUE), numeric(1))
  ctr <- vapply(vars, function(v){ x<-train[[v]]; x[!is.finite(x)]<-med[[v]]; mean(x) }, numeric(1))
  scl <- vapply(vars, function(v){ x<-train[[v]]; x[!is.finite(x)]<-med[[v]]; s<-stats::sd(x)
                                   if(!is.finite(s)||s<1e-9) 1 else s }, numeric(1))
  zmat <- function(d){
    M <- sapply(vars, function(v){ x<-d[[v]]; x[!is.finite(x)]<-med[[v]]; (x-ctr[[v]])/scl[[v]] })
    z <- as.data.frame(M); colnames(z) <- paste0("z", seq_along(vars)); z
  }
  Ztr <- zmat(train); Zall <- zmat(all_df)

  # design-based sampling variance D_i
  #
  # Most Admin-2 areas contain a single PSU, so survey::svymean cannot form a
  # between-cluster variance and returns SE ~ 1e-17 rather than 0. A `D <= 0`
  # guard does not catch that: 1e-17 is finite and positive, survives, and is
  # then floored at 1e-5 -- an implied SE of 0.32 pp for a district measured on
  # six children. The EBLUP weight gamma = A/(A+D) goes to ~1 and Fay-Herriot
  # hands back the raw direct estimate for exactly the districts that most need
  # shrinking. Measured on this project's data: mean gamma 0.71 as coded vs 0.26
  # with a real sampling variance, and 29% of districts got gamma > 0.95.
  #
  # So treat any SE that is implausibly small relative to p(1-p)/n as degenerate
  # too, not just SE <= 0.
  p <- pmin(pmax(train$svy_prev,1e-4),1-1e-4); D <- train$svy_prev_se^2
  D_srs <- p*(1-p)/pmax(train$n_svy/deff, 1)
  bad <- !is.finite(D) | D <= 0 | D < 0.25 * D_srs
  D[bad] <- D_srs[bad]; D <- pmax(D, 1e-5)

  # area-level regression (synthetic predictor)
  fit <- tryCatch(stats::lm(stats::reformulate(colnames(Ztr), "y"),
                            data = cbind(y = train$svy_prev, Ztr)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  n <- nrow(Ztr); pp <- sum(!is.na(stats::coef(fit)))
  A <- max(sum(stats::residuals(fit)^2)/max(n-pp,1) - mean(D), 1e-4)   # FH moment variance
  pr <- stats::predict(fit, newdata = Zall, se.fit = TRUE)
  mu_all <- as.numeric(pr$fit); se_reg <- as.numeric(pr$se.fit)

  est <- clip01(mu_all); se <- sqrt(A + se_reg^2)        # unsurveyed: synthetic + reg. uncertainty
  mi  <- match(train$Admin2, all_df$Admin2)
  gamma <- A/(A + D)
  est[mi] <- clip01(gamma*train$svy_prev + (1-gamma)*as.numeric(stats::predict(fit)))  # EBLUP
  se[mi]  <- sqrt(gamma * D)                              # EB posterior SD
  lo <- clip01(est - 1.96*se); hi <- clip01(est + 1.96*se)
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
  r <- tryCatch(fit_fh_country_outcome(low, oc), error=function(e){cat("  ERR",low,oc,e$message,"\n");NULL})
  if (!is.null(r) && nrow(r)>0) { rows[[paste(low,oc)]]<-r; ok<-ok+1
    cat(sprintf("  FH ok %-26s districts=%d surveyed=%d\n", paste(low,oc), nrow(r), sum(!is.na(r$obs_prev))))
  } else fail <- c(fail, paste(low,oc))
}
out <- if (length(rows)) do.call(rbind, rows) else NULL
if (!nzchar(ONE) && !is.null(out)) {
  saveRDS(out, file.path(DASH, "admin2_fh_predictions.rds"))
  cat(sprintf("\nadmin2_fh_predictions.rds: %d rows, %d combos OK\n", nrow(out), ok))
  if (length(fail)) cat("  no FH layer for:", paste(fail, collapse=", "), "\n")
  cat("  mean CI width (pp): surveyed=",
      round(mean(out$ci_width[!is.na(out$obs_prev)],na.rm=TRUE)*100,1),
      " unsurveyed=", round(mean(out$ci_width[is.na(out$obs_prev)],na.rm=TRUE)*100,1), "\n")
} else if (nzchar(ONE)) {
  cat("\n[ONE-slice test]\n")
  print(utils::head(out[order(-out$ci_width), c("Admin2","pred_prev","ci_lo","ci_hi","obs_prev")], 5))
}
