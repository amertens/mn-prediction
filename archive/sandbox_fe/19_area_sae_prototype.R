#!/usr/bin/env Rscript
# =============================================================================
# Prototype: Fay-Herriot area-level small-area estimation (SAE) with honest
# uncertainty, vs the current "treat area prevalence as a fixed number" approach.
#
# For each Admin2 we form:
#   * a DESIGN-BASED direct estimate p_dir and its sampling variance psi
#     (survey weights + cluster design; smoothed where a single PSU makes the
#      direct variance unstable)
#   * area-level proxy covariates (bundle-reduced, supervised top-k)
# Then fit FH:  p_dir_i = x_i'b + u_i + e_i,  e_i ~ N(0, psi_i),  u_i ~ N(0, s2u)
# EBLUP shrinks unreliable (few-cluster) areas toward the covariate model.
#
# Outputs: shrinkage behaviour, leave-area-out CV (FH vs synthetic vs direct),
# and the honest-uncertainty gap vs a naive individual-level binomial interval.
# =============================================================================
suppressMessages({library(here);library(dplyr);library(survey);library(sae)
  library(glmnet);library(matrixStats)})
source(here::here("R","config.R")); source(here::here("R","data_prep.R"))
options(survey.lonely.psu = "adjust")

# ---- design-based direct estimates + sampling variance per Admin2 ----------
area_direct <- function(merged, cc, oc) {
  d <- merged
  pc <- cc$child_flag
  if (!is.null(pc) && pc %in% colnames(d)) d <- d[d[[pc]] == oc$child_flag_val, ]
  if (is.null(oc$binary) || !oc$binary %in% colnames(d)) return(NULL)
  d <- d[!is.na(d[[oc$binary]]), ]
  y <- as.numeric(d[[oc$binary]]); yv <- sort(unique(y[!is.na(y)]))
  if (length(yv) == 2 && all(yv == c(1, 2))) y <- ifelse(y == 1, 1, 0)
  d$.y <- y
  d$.a2 <- as.character(d[[cc$admin2_col]]); d <- d[!is.na(d$.a2), ]
  d$.id <- d[[cc$cluster_id]]
  d$.w  <- if (!is.null(cc$weight_col) && cc$weight_col %in% colnames(d)) as.numeric(d[[cc$weight_col]]) else 1
  d$.w[is.na(d$.w) | d$.w <= 0] <- median(d$.w[d$.w > 0], na.rm = TRUE)
  des <- tryCatch(svydesign(ids = ~.id, weights = ~.w, data = d, nest = TRUE),
                  error = function(e) svydesign(ids = ~1, weights = ~.w, data = d))
  est <- svyby(~.y, ~.a2, des, svymean, na.rm = TRUE, vartype = "se")
  out <- data.frame(Admin2 = as.character(est$.a2), p_dir = as.numeric(est$.y),
                    se_dir = as.numeric(est$se), stringsAsFactors = FALSE)
  ntab <- as.data.frame(table(d$.a2)); names(ntab) <- c("Admin2", "n_ind")
  ncl  <- tapply(d$.id, d$.a2, function(x) length(unique(x)))
  out <- merge(out, ntab, by = "Admin2"); out$nclust <- as.integer(ncl[out$Admin2])
  # smooth sampling variance: psi = se^2; where unstable (NA/0/1-PSU) use
  # binomial p(1-p)/n inflated by the median design effect.
  out$psi_raw <- out$se_dir^2
  binom <- pmax(out$p_dir * (1 - out$p_dir), 1e-4) / out$n_ind
  ok <- is.finite(out$psi_raw) & out$psi_raw > 0 & out$nclust >= 2
  deff <- if (any(ok)) median(out$psi_raw[ok] / binom[ok], na.rm = TRUE) else 1.5
  deff <- max(deff, 1)
  out$psi <- ifelse(ok, out$psi_raw, binom * deff)
  out$deff_used <- deff
  out
}

# ---- area-level proxy covariate means (bundle, proxy-only) -----------------
area_covars <- function(merged, cc, oc, bundle = TRUE) {
  od <- build_outcome_dataset(merged, cc, oc)
  Xv <- if (bundle) od$Xvars_bundle else od$Xvars
  d <- od$data; a2 <- as.character(d[[cc$admin2_col]])
  M <- suppressWarnings(matrix(as.numeric(as.matrix(d[, Xv, drop = FALSE])), nrow = nrow(d)))
  colnames(M) <- Xv
  ua <- unique(a2[!is.na(a2)])
  Xa <- t(sapply(ua, function(u) colMeans(M[a2 == u, , drop = FALSE], na.rm = TRUE)))
  data.frame(Admin2 = ua, Xa, check.names = FALSE)
}

# supervised top-k covariate selection (correlation with p_dir), then DECORRELATE
# (drop near-duplicates |r|>0.9) so the FH design matrix is non-singular.
select_k <- function(df, ycol, xcols, k = 5) {
  cors <- sapply(xcols, function(v) {
    x <- df[[v]]; if (all(is.na(x)) || sd(x, na.rm = TRUE) < 1e-8) return(0)
    abs(suppressWarnings(cor(x, df[[ycol]], use = "complete.obs")))
  })
  ranked <- names(sort(cors[cors > 0], decreasing = TRUE))
  kept <- character(0)
  for (v in ranked) {
    if (length(kept) >= k) break
    if (length(kept) == 0) { kept <- v; next }
    rmax <- max(abs(suppressWarnings(cor(df[, kept, drop = FALSE], df[[v]],
                                         use = "complete.obs"))), na.rm = TRUE)
    if (is.finite(rmax) && rmax < 0.9) kept <- c(kept, v)
  }
  kept
}

run_country <- function(country, otag, k = 6) {
  cc <- get_country_configs()[[country]]; oc <- cc$outcomes[[otag]]
  merged <- load_merged_data(cc$data_path)
  dir <- area_direct(merged, cc, oc); cov <- area_covars(merged, cc, oc, bundle = TRUE)
  dat <- merge(dir, cov, by = "Admin2"); dat <- dat[is.finite(dat$p_dir) & dat$psi > 0, ]
  # median-impute covariates
  xcols <- setdiff(names(cov), "Admin2")
  for (v in xcols) { x <- dat[[v]]; if (any(is.na(x))) dat[[v]][is.na(x)] <- median(x, na.rm = TRUE) }
  xcols <- xcols[sapply(xcols, function(v) sd(dat[[v]], na.rm = TRUE) > 1e-8)]
  m <- nrow(dat)
  cat(sprintf("\n=== %s %s: %d areas (nclust: %d-%d), deff=%.2f ===\n",
              country, otag, m, min(dat$nclust), max(dat$nclust), dat$deff_used[1]))

  sel <- select_k(dat, "p_dir", xcols, k)
  # standardize selected covariates
  Z <- scale(as.matrix(dat[, sel, drop = FALSE]))
  fdat <- data.frame(p_dir = dat$p_dir, psi = dat$psi, Z)
  form <- as.formula(paste("p_dir ~", paste(colnames(Z), collapse = " + ")))
  fh <- tryCatch(sae::eblupFH(form, vardir = psi, data = fdat, method = "REML"),
                 error = function(e) { cat("  FH error:", conditionMessage(e), "\n"); NULL })
  if (is.null(fh)) return(NULL)
  mse <- tryCatch(sae::mseFH(form, vardir = psi, data = fdat, method = "REML"),
                  error = function(e) NULL)
  s2u <- fh$fit$refvar
  gamma <- s2u / (s2u + dat$psi)              # shrinkage weight on the direct est
  dat$eblup <- fh$eblup; dat$gamma <- gamma
  dat$mse_fh <- if (!is.null(mse)) mse$mse else NA
  fitted_is <- as.numeric(cbind(1, Z) %*% fh$fit$estcoef$beta)  # in-sample synthetic
  cor_is <- suppressWarnings(cor(dat$p_dir, fitted_is))
  cat(sprintf("  model var s2u=%.5f | mean shrinkage gamma=%.2f (range %.2f-%.2f)\n",
              s2u, mean(gamma), min(gamma), max(gamma)))
  cat(sprintf("  shrinkage: few-cluster areas pulled MORE to model (gamma nclust<=med=%.2f vs >med=%.2f)\n",
              mean(gamma[dat$nclust <= median(dat$nclust)]),
              mean(gamma[dat$nclust >  median(dat$nclust)])))

  # leave-area-out CV: predict held-out direct estimate (honest, out-of-sample)
  cv <- data.frame()
  for (i in seq_len(m)) {
    tr <- setdiff(seq_len(m), i)
    seli <- select_k(dat[tr, ], "p_dir", xcols, k)
    if (length(seli) < 1) next
    Zi <- scale(as.matrix(dat[, seli, drop = FALSE]))
    fdi <- data.frame(p_dir = dat$p_dir, psi = dat$psi, Zi)
    fo <- as.formula(paste("p_dir ~", paste(colnames(Zi), collapse = " + ")))
    fhi <- tryCatch(sae::eblupFH(fo, vardir = psi, data = fdi[tr, ], method = "REML"),
                    error = function(e) NULL)
    if (is.null(fhi)) next
    b <- fhi$fit$estcoef$beta
    synth <- sum(b * c(1, as.numeric(Zi[i, ])))
    cv <- rbind(cv, data.frame(p_dir = dat$p_dir[i], synth = synth,
                               mean_base = mean(dat$p_dir[tr])))
  }
  rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
  cor_cv <- suppressWarnings(cor(cv$p_dir, cv$synth))
  cat(sprintf("  OPTIMISM GAP  cor(p_dir, synthetic): in-sample=%.2f  ->  LOAO-CV=%.2f\n",
              cor_is, cor_cv))
  cat(sprintf("  LOAO-CV RMSE: covariate-model=%.3f   grand-mean=%.3f   (skill = %+0.1f%%)\n",
              rmse(cv$p_dir, cv$synth), rmse(cv$p_dir, cv$mean_base),
              100*(1 - rmse(cv$p_dir, cv$synth)/rmse(cv$p_dir, cv$mean_base))))

  # uncertainty: per-area direct estimate is too noisy to map; FH borrows strength
  dir_hw <- 1.96 * sqrt(dat$psi)   # smoothed direct half-width (design ~non-estimable: ~1 PSU/area)
  fh_hw  <- 1.96 * sqrt(dat$mse_fh)
  cat(sprintf("  median 95%% CI half-width: direct(smoothed)=%.3f  ->  FH-EBLUP=%.3f  (%.0f%% narrower)\n",
              median(dir_hw, na.rm = TRUE), median(fh_hw, na.rm = TRUE),
              100*(1 - median(fh_hw,na.rm=TRUE)/median(dir_hw,na.rm=TRUE))))
  dat$country <- country; dat$outcome <- otag; dat$fh_hw <- fh_hw; dat$dir_hw <- dir_hw
  list(dat = dat,
       summary = data.frame(country = country, outcome = otag, m_areas = m,
                            s2u = s2u, gamma_mean = mean(gamma),
                            cor_insample = round(cor_is,2), cor_loaocv = round(cor_cv,2),
                            rmse_synth = round(rmse(cv$p_dir, cv$synth),3),
                            rmse_mean = round(rmse(cv$p_dir, cv$mean_base),3),
                            hw_direct = round(median(dir_hw, na.rm = TRUE),3),
                            hw_fh = round(median(fh_hw, na.rm = TRUE),3)))
}

cells <- list(c("Ghana","child_vitA"), c("Malawi","child_vitA"),
              c("Ghana","child_iron"), c("Gambia","women_iron"))
res <- list(); summ <- list()
for (cl in cells) {
  r <- tryCatch(run_country(cl[1], cl[2]), error = function(e){cat("ERR",cl[1],cl[2],conditionMessage(e),"\n");NULL})
  if (is.null(r)) next
  res[[paste(cl,collapse="_")]] <- r$dat; summ[[length(summ)+1]] <- r$summary
}
S <- do.call(rbind, summ)
saveRDS(res, here::here("sandbox_fe","results_19_sae.rds"))
write.csv(S, here::here("sandbox_fe","results_19_sae_summary.csv"), row.names = FALSE)
cat("\n=== SUMMARY ===\n"); print(S, row.names = FALSE)
