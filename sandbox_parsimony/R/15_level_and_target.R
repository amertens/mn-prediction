# =============================================================================
# sandbox_parsimony/R/15_level_and_target.R
#
# Two structural questions that matter more than the choice of learner.
#
# A. WHAT SPATIAL LEVEL should the model be fitted and reported at?
#    Reliability of the survey estimate at cluster / Admin-2 / Admin-1, so the
#    ceiling at each level is visible. Also counts the design points available
#    at each level, since that -- not the number of individuals -- is what a
#    proxy-only model actually has to learn from.
#
# B. WHAT TARGET: the deficiency INDICATOR or the CONTINUOUS biomarker?
#    Thresholding a biomarker at a WHO cutoff throws away everything about how
#    far from the cutoff each person sits. If the district mean of the (log)
#    biomarker is a more reliable quantity than the district prevalence, then
#    modelling it and converting to prevalence afterwards should beat modelling
#    the prevalence directly.
# =============================================================================
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi")
`%|%` <- function(a, b) if (is.null(a)) b else a

sink(stdout(), type = "message")
source("R/config.R")
cfg <- get_country_configs()
cfg_key <- c(gambia = "Gambia", ghana = "Ghana",
             sierraleone = "SierraLeone", malawi = "Malawi")

.num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  suppressWarnings(as.numeric(v))
}

#' Should this biomarker be modelled on the log scale?
#' Ferritin, RBP, folate, B12 and zinc are all right-skewed concentrations, so
#' the district MEAN is better behaved in logs. Some configs already supply a
#' logged column (Gambia gw_LogFerAdj), which must not be logged twice.
#' Decided once, here, and applied to BOTH the biomarker and its WHO cutoff --
#' getting that pairing wrong silently produces nonsense conversions.
should_log <- function(x) {
  xf <- x[is.finite(x)]
  if (!length(xf) || any(xf <= 0)) return(FALSE)
  stats::quantile(xf, .99) > 20
}

#' Weighted prevalence and weighted mean (log) biomarker by an arbitrary grouping
#' @return data.frame with attribute "logged" recording the scale decision
agg_by <- function(d, grp, bin, cont, wt) {
  g <- as.character(d[[grp]])
  y <- .num(d[[bin]])
  x <- if (!is.null(cont) && cont %in% names(d)) .num(d[[cont]]) else rep(NA_real_, nrow(d))
  logged <- should_log(x)
  if (logged) x <- log(x)
  w <- .num(d[[wt]]); w[!is.finite(w) | w <= 0] <- 1
  ok <- !is.na(g)
  out <- data.frame(g = g[ok], y = y[ok], x = x[ok], w = w[ok]) |>
    group_by(g) |>
    summarise(
      n_bin  = sum(is.finite(y)),
      prev   = if (sum(is.finite(y))) stats::weighted.mean(y[is.finite(y)], w[is.finite(y)]) else NA_real_,
      n_cont = sum(is.finite(x)),
      mu     = if (sum(is.finite(x))) stats::weighted.mean(x[is.finite(x)], w[is.finite(x)]) else NA_real_,
      sd_w   = if (sum(is.finite(x)) > 2) stats::sd(x[is.finite(x)]) else NA_real_,
      .groups = "drop")
  attr(out, "logged") <- logged
  out
}

#' Reliability of a MEAN: lambda = 1 - mean(sigma^2/n) / var(mu)
rel_mean <- function(mu, sd_w, n, deff = 1.5) {
  ok <- is.finite(mu) & is.finite(sd_w) & is.finite(n) & n > 1
  if (sum(ok) < 5) return(NA_real_)
  s2 <- stats::median(sd_w[ok]^2, na.rm = TRUE)          # pooled within-unit variance
  v_obs <- stats::var(mu[ok])
  max(0, v_obs - mean(deff * s2 / n[ok])) / v_obs
}

rows <- list(); conv <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[cfg_key[[ctry]]]]; if (is.null(cc)) next
  for (oc in cc$outcomes) {
    f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
    if (!file.exists(f)) next
    od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) next
    d <- od$data
    need <- c(cc$admin1_col, cc$admin2_col, cc$weight_col, oc$binary)
    if (!all(need %in% names(d))) next
    clus <- if (!is.null(cc$cluster_id) && cc$cluster_id %in% names(d)) cc$cluster_id else NULL

    lv <- list(Admin1 = cc$admin1_col, Admin2 = cc$admin2_col)
    if (!is.null(clus)) lv$cluster <- clus

    for (lvl in names(lv)) {
      a <- agg_by(d, lv[[lvl]], oc$binary, oc$continuous, cc$weight_col)
      a <- a[is.finite(a$prev), ]
      if (nrow(a) < 5) next
      rows[[paste(ctry, oc$tag, lvl)]] <- data.frame(
        country = ctry, outcome = oc$tag, level = lvl,
        n_units = nrow(a), median_n = stats::median(a$n_bin),
        # ceiling for the PREVALENCE target
        lambda_prev = round(reliability(a$prev, a$n_bin, 1.5)$lambda, 3),
        r_max_prev  = round(reliability(a$prev, a$n_bin, 1.5)$r_max, 3),
        # ceiling for the CONTINUOUS target
        lambda_cont = round(rel_mean(a$mu, a$sd_w, a$n_cont, 1.5), 3),
        r_max_cont  = round(sqrt(pmax(rel_mean(a$mu, a$sd_w, a$n_cont, 1.5), 0)), 3),
        stringsAsFactors = FALSE)
    }

    # --- can the continuous mean reproduce the prevalence at all? ------------
    # If prevalence is a deterministic function of the district mean and a
    # common within-district SD, then predicting the mean is strictly more
    # informative. Check how well that mapping holds empirically at Admin-2.
    a2_all <- agg_by(d, cc$admin2_col, oc$binary, oc$continuous, cc$weight_col)
    logged <- isTRUE(attr(a2_all, "logged"))
    a2 <- a2_all[is.finite(a2_all$prev) & is.finite(a2_all$mu) & is.finite(a2_all$sd_w), ]
    if (nrow(a2) >= 10 && !is.null(oc$cutoff)) {
      cut_raw <- as.numeric(oc$cutoff)
      # the cutoff MUST take the same transform the biomarker took
      cut_use <- if (logged) log(cut_raw) else cut_raw
      s <- stats::median(a2$sd_w, na.rm = TRUE)
      dir_less <- (oc$cutoff_dir %|% "less") == "less"
      p_hat <- if (dir_less) stats::pnorm(cut_use, a2$mu, s) else
        1 - stats::pnorm(cut_use, a2$mu, s)
      conv[[paste(ctry, oc$tag)]] <- data.frame(
        country = ctry, outcome = oc$tag, n_areas = nrow(a2),
        r_mu_vs_prev = round(suppressWarnings(cor(a2$mu, a2$prev)), 3),
        rho_mu_vs_prev = round(suppressWarnings(cor(a2$mu, a2$prev, method = "spearman")), 3),
        r_convert = round(suppressWarnings(cor(p_hat, a2$prev)), 3),
        mae_convert_pp = round(mean(abs(p_hat - a2$prev)) * 100, 2),
        bias_convert_pp = round(mean(p_hat - a2$prev) * 100, 2),
        stringsAsFactors = FALSE)
    }
  }
}

lev <- bind_rows(rows)
write.csv(lev, "sandbox_parsimony/out/level_reliability.csv", row.names = FALSE)
cv <- bind_rows(conv)
write.csv(cv, "sandbox_parsimony/out/continuous_vs_binary.csv", row.names = FALSE)

cat("\n=== A. Ceiling by spatial level ===\n")
cat("n_units is the number of DESIGN POINTS a proxy-only model gets to learn from.\n\n")
print(as.data.frame(lev |> arrange(country, outcome, level)), row.names = FALSE)

cat("\n=== A2. Averages by level ===\n")
print(as.data.frame(lev |> group_by(level) |>
  summarise(cells = n(), median_units = round(median(n_units)),
            median_n_per_unit = round(median(median_n)),
            mean_r_max_prev = round(mean(r_max_prev, na.rm = TRUE), 3),
            mean_r_max_cont = round(mean(r_max_cont, na.rm = TRUE), 3),
            .groups = "drop")), row.names = FALSE)

cat("\n=== B. Continuous biomarker mean vs the deficiency indicator (Admin-2) ===\n")
cat("r_convert: prevalence implied by the OBSERVED district mean + a common within-\n")
cat("district SD, correlated with observed prevalence. That is the ceiling on the\n")
cat("two-step route, before any covariate model is involved.\n\n")
print(as.data.frame(cv |> arrange(desc(abs(r_mu_vs_prev)))), row.names = FALSE)

# ===========================================================================
# C. Head-to-head: predict the indicator, or predict the biomarker and convert?
#    Same learner (curated 16 + lon/lat random forest), same folds, same
#    evaluation target (observed Admin-2 prevalence). The only difference is
#    what the model is asked to learn.
# ===========================================================================
source("sandbox_parsimony/R/02_features.R")
suppressPackageStartupMessages(library(ranger))
pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")

h2h <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[cfg_key[[ctry]]]]; if (is.null(cc)) next
  for (oc in cc$outcomes) {
    P <- pooled_all[[oc$tag]]; if (is.null(P)) next
    f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
    if (!file.exists(f) || is.null(oc$cutoff)) next
    od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) next
    d <- od$data
    if (!all(c(cc$admin2_col, cc$weight_col, oc$binary) %in% names(d))) next

    a2 <- agg_by(d, cc$admin2_col, oc$binary, oc$continuous, cc$weight_col)
    logged <- isTRUE(attr(a2, "logged"))
    cut_use <- if (logged) log(as.numeric(oc$cutoff)) else as.numeric(oc$cutoff)
    dir_less <- (oc$cutoff_dir %|% "less") == "less"

    cov <- P$data[P$data$country == ctry, , drop = FALSE]
    m <- dplyr::inner_join(
      a2 |> rename(Admin2 = g) |> filter(is.finite(prev), is.finite(mu), is.finite(sd_w)),
      cov |> select(-any_of(c("svy_prev", "n_svy", "svy_prev_se"))), by = "Admin2")
    if (nrow(m) < 25) next

    cur <- intersect(curated_vars(P$predictors), names(m))
    n <- nrow(m); s_pool <- stats::median(m$sd_w, na.rm = TRUE)
    reps <- list()
    for (rep_i in 1:10) {
      set.seed(800L + rep_i)
      fold <- sample(rep_len(1:5, n))
      p_bin <- rep(NA_real_, n); p_cont <- rep(NA_real_, n); p_cal <- rep(NA_real_, n); p_cal2 <- rep(NA_real_, n)
      for (fi in 1:5) {
        te <- which(fold == fi); tr <- setdiff(seq_len(n), te)
        pp <- prep_X(m[tr, , drop = FALSE], m[te, , drop = FALSE], cur)
        Xtr <- cbind(pp$Xtr, lon = m$lon[tr], lat = m$lat[tr])
        Xte <- cbind(pp$Xte, lon = m$lon[te], lat = m$lat[te])
        w <- pmax(m$n_bin[tr], 1)
        # (a) learn the prevalence directly, on the logit scale
        f1 <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = .logit(m$prev[tr]), Xtr,
                                                               check.names = TRUE),
                                      num.trees = 800, min.node.size = 5,
                                      case.weights = w, seed = rep_i), error = function(e) NULL)
        p_bin[te] <- if (is.null(f1)) mean(m$prev[tr]) else
          .ilogit(stats::predict(f1, data = data.frame(Xte, check.names = TRUE))$predictions)
        # (b) learn the district mean biomarker, then convert with a common SD
        f2 <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = m$mu[tr], Xtr,
                                                               check.names = TRUE),
                                      num.trees = 800, min.node.size = 5,
                                      case.weights = w, seed = rep_i), error = function(e) NULL)
        mu_hat <- if (is.null(f2)) mean(m$mu[tr]) else
          stats::predict(f2, data = data.frame(Xte, check.names = TRUE))$predictions
        s_tr <- stats::median(m$sd_w[tr], na.rm = TRUE)
        conv_fn <- function(mu) if (dir_less) stats::pnorm(cut_use, mu, s_tr) else
          1 - stats::pnorm(cut_use, mu, s_tr)
        p_cont[te] <- conv_fn(mu_hat)

        # (c) same, then recalibrate the mean -> prevalence conversion on the
        # TRAINING fold. Assuming one within-district SD and normality is what
        # puts the level off by 5-20 pp in some cells.
        #
        # Two flavours, neither of which is free:
        #   slope+intercept (p_cal2): at ~24 training districts the fitted slope
        #     is noisy and occasionally negative, which flips the ranking;
        #   intercept-only (p_cal): a monotone shift WITHIN a fold, so it looks
        #     like it should leave Spearman untouched -- but the shift is
        #     re-estimated per fold, so the assembled out-of-fold vector is only
        #     piecewise monotone and the ranking still moves (and .logit clamps
        #     at 0.005/0.995, which can tie values together).
        # Empirically both recover the level and give back most of the ranking
        # gain. The conclusion is that the fix is a better generative model of
        # the district biomarker distribution, not a better post-hoc patch.
        p_tr_raw <- conv_fn(m$mu[tr])
        shift <- mean(.logit(m$prev[tr])) - mean(.logit(p_tr_raw))
        p_cal[te] <- .ilogit(.logit(p_cont[te]) + shift)
        cal2 <- tryCatch(stats::lm(.logit(m$prev[tr]) ~ .logit(p_tr_raw)),
                         error = function(e) NULL)
        p_cal2[te] <- if (is.null(cal2)) p_cont[te] else
          .ilogit(stats::predict(cal2, newdata = data.frame(p_tr_raw = p_cont[te])))
      }
      reps[[rep_i]] <- data.frame(
        r_bin  = suppressWarnings(cor(p_bin, m$prev)),
        rho_bin = suppressWarnings(cor(p_bin, m$prev, method = "spearman")),
        mae_bin = mean(abs(p_bin - m$prev)) * 100,
        r_cont  = suppressWarnings(cor(p_cont, m$prev)),
        rho_cont = suppressWarnings(cor(p_cont, m$prev, method = "spearman")),
        mae_cont = mean(abs(p_cont - m$prev)) * 100,
        bias_cont = mean(p_cont - m$prev) * 100,
        rho_cal = suppressWarnings(cor(p_cal, m$prev, method = "spearman")),
        mae_cal = mean(abs(p_cal - m$prev)) * 100,
        bias_cal = mean(p_cal - m$prev) * 100,
        rho_cal2 = suppressWarnings(cor(p_cal2, m$prev, method = "spearman")),
        mae_cal2 = mean(abs(p_cal2 - m$prev)) * 100)
    }
    r <- bind_rows(reps)
    h2h[[paste(ctry, oc$tag)]] <- data.frame(
      country = ctry, outcome = oc$tag, n_areas = n, logged = logged,
      rho_indicator = round(mean(r$rho_bin), 3),
      rho_biomarker = round(mean(r$rho_cont), 3),
      rho_gain = round(mean(r$rho_cont - r$rho_bin), 3),
      mae_indicator_pp = round(mean(r$mae_bin), 2),
      mae_biomarker_pp = round(mean(r$mae_cont), 2),
      bias_biomarker_pp = round(mean(r$bias_cont), 2),
      rho_biomarker_shift = round(mean(r$rho_cal), 3),
      mae_biomarker_shift_pp = round(mean(r$mae_cal), 2),
      bias_biomarker_shift_pp = round(mean(r$bias_cal), 2),
      rho_biomarker_lmcal = round(mean(r$rho_cal2), 3),
      mae_biomarker_lmcal_pp = round(mean(r$mae_cal2), 2),
      stringsAsFactors = FALSE)
  }
}
h <- bind_rows(h2h)
write.csv(h, "sandbox_parsimony/out/head_to_head_target.csv", row.names = FALSE)

cat("\n=== C. Head-to-head, same learner and folds, scored on observed prevalence ===\n")
print(as.data.frame(h |> arrange(desc(rho_gain))), row.names = FALSE)
cat("\npaired over cells:\n")
if (nrow(h) > 2) {
  di <- h$rho_gain
  cat(sprintf("  mean rho gain from modelling the biomarker = %+.3f (SE %.3f, t = %.2f), better in %d of %d cells\n",
              mean(di), sd(di)/sqrt(length(di)), mean(di)/(sd(di)/sqrt(length(di))),
              sum(di > 0), length(di)))
  dm <- h$mae_indicator_pp - h$mae_biomarker_pp
  cat(sprintf("  mean MAE change = %+.2f pp (positive = biomarker route better), better in %d of %d cells\n",
              mean(dm), sum(dm > 0), length(dm)))
  for (v in c("shift", "lmcal")) {
    rc <- paste0("rho_biomarker_", v); mc <- paste0("mae_biomarker_", v, "_pp")
    if (!all(c(rc, mc) %in% names(h))) next
    dc <- h[[rc]] - h$rho_indicator
    dmc <- h$mae_indicator_pp - h[[mc]]
    cat(sprintf("  %-6s recalibration: rho gain %+.3f (SE %.3f, t = %.2f, %d/%d cells) | MAE gain %+.2f pp (%d/%d)\n",
                v, mean(dc), sd(dc) / sqrt(length(dc)),
                mean(dc) / (sd(dc) / sqrt(length(dc))), sum(dc > 0), length(dc),
                mean(dmc), sum(dmc > 0), length(dmc)))
  }
}
message("\nSaved -> out/level_reliability.csv, out/continuous_vs_binary.csv, out/head_to_head_target.csv")
