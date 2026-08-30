# =============================================================================
# sandbox_parsimony/R/35_hierarchical_biomarker.R
#
# STEP 5: the one untested route to the extra signal in the continuous target.
#
# FINDINGS.md sections 13 and 16: the district mean log-biomarker carries about
# twice the signal of the district prevalence (r_max 0.63 vs 0.31), but the
# two-step route -- predict the mean, convert with a plug-in normal CDF -- gives
# the gain straight back, and modelling the district SD made it worse because
# that SD is itself estimated from ~13 people.
#
# The remaining idea was a HIERARCHICAL model fitted to the INDIVIDUAL
# biomarkers, with a district random effect:
#
#   log(biomarker)_ij = X_j beta + u_j + e_ij ,  u_j ~ N(0, tau^2), e_ij ~ N(0, sigma^2)
#   p_j = Phi( (cut - (X_j beta + u_j)) / sigma )
#
# Why this could beat the two-step even though the moment version did not:
#   * sigma is estimated from EVERY individual in the country, not from the ~13
#     in one district, so the conversion no longer divides by a noisy number;
#   * districts are shrunk toward the covariate model by the ratio
#     tau^2/(tau^2 + sigma^2/n_j) automatically, with no plug-in step;
#   * the individual likelihood uses how far each person sits from the cutoff,
#     which thresholding throws away.
#
# For an UNSURVEYED district u_j is unknown and set to 0, which is what the
# product requires. Evaluation is leave-district-out, so the held-out district
# contributes no individuals to the fit and its u_j is genuinely unavailable.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr); library(lme4); library(ranger)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi")
CFG_KEY <- c(gambia = "Gambia", ghana = "Ghana",
             sierraleone = "SierraLeone", malawi = "Malawi")
WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"
`%|%` <- function(a, b) if (is.null(a)) b else a

sink(stdout(), type = "message")
source("R/config.R")
cfg <- get_country_configs()
cov_ext <- readRDS("sandbox_parsimony/out/cov_ext.rds")

.num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  suppressWarnings(as.numeric(v))
}
should_log <- function(x) {
  xf <- x[is.finite(x)]
  if (!length(xf) || any(xf <= 0)) return(FALSE)
  stats::quantile(xf, .99) > 20
}

#' Individual biomarker rows joined to district covariates
build_individual <- function(ctry, oc, cc) {
  f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
  if (!file.exists(f) || is.null(oc$continuous) || is.null(oc$cutoff)) return(NULL)
  od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) return(NULL)
  d <- od$data
  if (!all(c(cc$admin2_col, cc$weight_col, oc$binary) %in% names(d)) ||
      !oc$continuous %in% names(d)) return(NULL)
  x <- .num(d[[oc$continuous]]); lg <- should_log(x); if (lg) x <- log(x)
  out <- data.frame(
    Admin2 = as.character(d[[cc$admin2_col]]),
    x = x, y = .num(d[[oc$binary]]),
    # NORMALISE the survey weight to mean 1 within country.
    #
    # The weight columns are NOT on a common scale: Ghana's median is 0.98 and
    # Gambia's 0.65, but Malawi's is 790,000 (the raw DHS 1e6 convention).
    # A weighted MEAN is scale-invariant so compute_svy_admin2() is unaffected,
    # but lme4 treats `weights` as PRECISION -- each Malawian would count as
    # ~790,000 observations, collapsing the residual SD and destroying the fit.
    # (Left unnormalised this route scored MAE 25.5 pp, worse than doing
    # nothing.) Anywhere a survey weight is used as anything other than a
    # weighted average, it has to be rescaled first.
    w = { ww <- .num(d[[cc$weight_col]]); ww[!is.finite(ww) | ww <= 0] <- NA
          ww <- ww / mean(ww, na.rm = TRUE); ww[!is.finite(ww)] <- 1; ww },
    stringsAsFactors = FALSE)
  out <- out[is.finite(out$x) & !is.na(out$Admin2) &
               !grepl(WATER, out$Admin2, ignore.case = TRUE), , drop = FALSE]
  attr(out, "cut") <- if (lg) log(as.numeric(oc$cutoff)) else as.numeric(oc$cutoff)
  attr(out, "dir_less") <- (oc$cutoff_dir %|% "less") == "less"
  out
}

run_cell <- function(ctry, oc, cc, label) {
  ind <- build_individual(ctry, oc, cc); if (is.null(ind) || nrow(ind) < 200) return(NULL)
  cut <- attr(ind, "cut"); dir_less <- attr(ind, "dir_less")

  # observed district prevalence (survey weighted) -- the evaluation target
  obs <- ind |> group_by(Admin2) |>
    summarise(prev = stats::weighted.mean(y, w, na.rm = TRUE),
              n_svy = sum(is.finite(y)), .groups = "drop") |>
    filter(is.finite(prev)) |> as.data.frame()

  cv <- cov_ext[[ctry]]; if (is.null(cv)) return(NULL)
  dd <- dplyr::inner_join(obs, cv, by = "Admin2")
  dd <- dd[is.finite(dd$lon), , drop = FALSE]
  if (nrow(dd) < 25) return(NULL)
  cur <- intersect(curated_vars(names(dd)), names(dd))
  if (length(cur) < 6) return(NULL)

  ind <- ind[ind$Admin2 %in% dd$Admin2, , drop = FALSE]
  n <- nrow(dd); set.seed(11)
  fold <- sample(rep_len(1:5, n))
  p_h <- p_2s <- p_ind <- rep(NA_real_, n)

  for (f in 1:5) {
    te <- which(fold == f); tr <- setdiff(seq_len(n), te)
    pp <- prep_X(dd[tr, , drop = FALSE], dd[te, , drop = FALSE], cur)
    Xtr <- cbind(pp$Xtr, lon = dd$lon[tr], lat = dd$lat[tr])
    Xte <- cbind(pp$Xte, lon = dd$lon[te], lat = dd$lat[te])
    nm <- paste0("z", seq_len(ncol(Xtr)))
    colnames(Xtr) <- colnames(Xte) <- nm

    # ---- (1) hierarchical individual model ------------------------------
    # District covariates are attached to every individual in that district;
    # the random intercept then absorbs whatever the covariates miss.
    tr_names <- dd$Admin2[tr]
    it <- ind[ind$Admin2 %in% tr_names, , drop = FALSE]
    Xmap <- as.data.frame(Xtr); Xmap$Admin2 <- tr_names
    it <- dplyr::left_join(it, Xmap, by = "Admin2")
    form <- stats::as.formula(paste("x ~", paste(nm, collapse = " + "), "+ (1 | Admin2)"))
    m <- tryCatch(lme4::lmer(form, data = it, weights = it$w,
                             control = lme4::lmerControl(calc.derivs = FALSE)),
                  error = function(e) NULL)
    if (!is.null(m)) {
      nd <- as.data.frame(Xte)
      mu <- tryCatch(as.numeric(stats::predict(m, newdata = nd, re.form = ~0)),
                     error = function(e) NULL)   # u_j = 0: unsurveyed district
      if (!is.null(mu)) {
        # The predictive spread for an UNSURVEYED district is NOT the residual
        # SD alone. With u_j unknown,
        #   P(x < cut) = E_u[ Phi((cut - Xb - u)/sigma) ] = Phi((cut - Xb)/sqrt(sigma^2 + tau^2))
        # so the between-district variance has to be added back. Using sigma
        # alone integrates the random effect out of the MEAN but not out of the
        # VARIANCE, which makes every district look alike and inflates MAE
        # roughly two-fold. (Measured: 25.5 pp vs 12.4 pp once corrected.)
        sig <- stats::sigma(m)
        tau2 <- tryCatch({
          vc <- as.data.frame(lme4::VarCorr(m))
          v <- vc$vcov[vc$grp == "Admin2"]
          if (length(v)) v[1] else 0
        }, error = function(e) 0)
        sig_marg <- sqrt(sig^2 + max(tau2, 0))
        p_h[te] <- if (dir_less) stats::pnorm(cut, mu, sig_marg) else
          1 - stats::pnorm(cut, mu, sig_marg)
      }
    }

    # ---- (2) two-step control: district means, pooled SD ----------------
    mom <- it |> group_by(Admin2) |>
      summarise(mu = stats::weighted.mean(x, w), sd = stats::sd(x), .groups = "drop")
    mi <- match(dd$Admin2[tr], mom$Admin2)
    rf <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = mom$mu[mi], Xtr),
                                  num.trees = 800, min.node.size = 5, seed = 1),
                   error = function(e) NULL)
    if (!is.null(rf)) {
      mu2 <- stats::predict(rf, data = as.data.frame(Xte))$predictions
      s2 <- stats::median(mom$sd[mi], na.rm = TRUE)
      p_2s[te] <- if (dir_less) stats::pnorm(cut, mu2, s2) else
        1 - stats::pnorm(cut, mu2, s2)
    }

    # ---- (3) indicator control: predict the prevalence directly ---------
    rf3 <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = .logit(dd$prev[tr]), Xtr),
                                   num.trees = 800, min.node.size = 5,
                                   case.weights = pmax(dd$n_svy[tr], 1), seed = 1),
                    error = function(e) NULL)
    if (!is.null(rf3))
      p_ind[te] <- .ilogit(stats::predict(rf3, data = as.data.frame(Xte))$predictions)
  }

  sc <- function(p, tag) {
    ok <- is.finite(p)
    if (sum(ok) < 10) return(NULL)
    data.frame(cell = label, route = tag, n = sum(ok),
               mae_pp = round(100 * mean(abs(p[ok] - dd$prev[ok])), 2),
               pearson = round(suppressWarnings(cor(p[ok], dd$prev[ok])), 3),
               spearman = round(suppressWarnings(
                 cor(p[ok], dd$prev[ok], method = "spearman")), 3),
               stringsAsFactors = FALSE)
  }
  bind_rows(sc(p_ind, "indicator (production)"),
            sc(p_2s, "two-step mean + pooled SD"),
            sc(p_h, "hierarchical individual"))
}

res <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  for (oc in cc$outcomes) {
    lab <- paste(ctry, oc$tag)
    r <- tryCatch(run_cell(ctry, oc, cc, lab), error = function(e) {
      message("  ", lab, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) { res[[lab]] <- r; message("  ", lab, " done") }
  }
}

out <- bind_rows(res)
write.csv(out, "sandbox_parsimony/out/hierarchical_biomarker.csv", row.names = FALSE)

cat("\n=== Three routes to an Admin-2 prevalence, same folds and covariates ===\n")
print(as.data.frame(out |> group_by(route) |>
  summarise(cells = n(),
            mae_pp = round(mean(mae_pp, na.rm = TRUE), 2),
            pearson = round(mean(pearson, na.rm = TRUE), 3),
            spearman = round(mean(spearman, na.rm = TRUE), 3),
            .groups = "drop") |> arrange(mae_pp)), row.names = FALSE)

w <- out |> select(cell, route, spearman, mae_pp) |>
  tidyr::pivot_wider(names_from = route, values_from = c(spearman, mae_pp))
w <- w[stats::complete.cases(w), ]
if (nrow(w) > 2) {
  for (v in c("two-step mean + pooled SD", "hierarchical individual")) {
    ds <- w[[paste0("spearman_", v)]] - w[["spearman_indicator (production)"]]
    dm <- w[["mae_pp_indicator (production)"]] - w[[paste0("mae_pp_", v)]]
    cat(sprintf("\n%-28s vs indicator: rho %+0.3f (t=%.2f, %d/%d) | MAE %+0.2f pp (%d/%d)",
                v, mean(ds), mean(ds)/(sd(ds)/sqrt(nrow(w))), sum(ds > 0), nrow(w),
                mean(dm), sum(dm > 0), nrow(w)))
  }
  cat("\n")
}
message("\nSaved -> sandbox_parsimony/out/hierarchical_biomarker.csv")
