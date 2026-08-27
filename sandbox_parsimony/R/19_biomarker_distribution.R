# =============================================================================
# sandbox_parsimony/R/19_biomarker_distribution.R
#
# NEXT STEP 2: model the district biomarker DISTRIBUTION, not just its mean.
#
# FINDINGS.md section 13 found that the district mean log-biomarker carries
# roughly twice the signal of the district prevalence (r_max 0.63 vs 0.31), but
# that converting a predicted mean back to a prevalence with ONE pooled
# within-district SD gives the gain straight back: the level goes wrong by
# 5-20 pp, and every post-hoc recalibration that fixes the level also destroys
# the ranking.
#
# The diagnosis there was that the two-step route is the wrong shape. So model
# location AND scale:
#
#   mu_a    = district mean of log(biomarker)      <- covariate model
#   sigma_a = district SD of log(biomarker)        <- covariate model
#   p_a     = Phi((cut - mu_a) / sigma_a)          <- no free parameters
#
# Both mu and sigma are learned from the proxies, so a district whose biomarker
# distribution is unusually wide is allowed to have a higher deficiency rate at
# the same mean -- which the single-sigma version cannot express.
#
# Four contenders, all scored on the SAME folds against the SAME observed
# Admin-2 prevalence:
#   indicator          predict the prevalence directly (production behaviour)
#   mean_pooled_sd     predict mu, convert with one pooled sigma (section 13)
#   loc_scale          predict mu AND sigma, convert
#   loc_scale_emp      predict mu and sigma, but convert with the EMPIRICAL
#                      within-district distribution shape instead of a normal
#                      (skewness of log-biomarkers is not always zero)
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr); library(ranger)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi")
CFG_KEY <- c(gambia = "Gambia", ghana = "Ghana",
             sierraleone = "SierraLeone", malawi = "Malawi")
`%|%` <- function(a, b) if (is.null(a)) b else a
WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"

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

#' District-level location, scale, skew and prevalence of one biomarker
district_moments <- function(ctry, oc, cc) {
  f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
  if (!file.exists(f) || is.null(oc$continuous) || is.null(oc$cutoff)) return(NULL)
  od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) return(NULL)
  d <- od$data
  if (!all(c(cc$admin2_col, cc$weight_col, oc$binary) %in% names(d)) ||
      !oc$continuous %in% names(d)) return(NULL)
  x <- .num(d[[oc$continuous]]); lg <- should_log(x); if (lg) x <- log(x)
  y <- .num(d[[oc$binary]])
  w <- .num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- 1
  g <- as.character(d[[cc$admin2_col]])
  ok <- !is.na(g)
  out <- data.frame(Admin2 = g[ok], x = x[ok], y = y[ok], w = w[ok]) |>
    group_by(Admin2) |>
    summarise(
      n_svy    = sum(is.finite(y)),
      svy_prev = if (sum(is.finite(y))) stats::weighted.mean(y[is.finite(y)], w[is.finite(y)]) else NA_real_,
      n_cont   = sum(is.finite(x)),
      mu       = if (sum(is.finite(x))) stats::weighted.mean(x[is.finite(x)], w[is.finite(x)]) else NA_real_,
      sigma    = if (sum(is.finite(x)) > 3) stats::sd(x[is.finite(x)]) else NA_real_,
      # standardised residuals pooled later to build an empirical shape
      .groups = "drop") |> as.data.frame()
  out <- out[!grepl(WATER, out$Admin2, ignore.case = TRUE), , drop = FALSE]
  attr(out, "logged") <- lg
  attr(out, "cut") <- if (lg) log(as.numeric(oc$cutoff)) else as.numeric(oc$cutoff)
  attr(out, "dir_less") <- (oc$cutoff_dir %|% "less") == "less"
  # pooled empirical distribution of within-district standardised biomarker,
  # used by the loc_scale_emp variant instead of assuming normality
  z <- data.frame(Admin2 = g[ok], x = x[ok]) |> group_by(Admin2) |>
    mutate(z = (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)) |>
    ungroup() |> pull(z)
  attr(out, "z_pool") <- z[is.finite(z)]
  out
}

conv_normal <- function(mu, sigma, cut, dir_less) {
  sigma <- pmax(sigma, 1e-6)
  if (dir_less) stats::pnorm(cut, mu, sigma) else 1 - stats::pnorm(cut, mu, sigma)
}
conv_empirical <- function(mu, sigma, cut, dir_less, z_pool) {
  sigma <- pmax(sigma, 1e-6)
  zc <- (cut - mu) / sigma
  # share of the pooled standardised distribution below each district's cutoff
  f <- stats::ecdf(z_pool)
  p <- f(zc)
  if (dir_less) p else 1 - p
}

rows <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  cv <- cov_ext[[ctry]]; if (is.null(cv)) next
  for (oc in cc$outcomes) {
    dm <- district_moments(ctry, oc, cc); if (is.null(dm)) next
    cut <- attr(dm, "cut"); dir_less <- attr(dm, "dir_less"); zp <- attr(dm, "z_pool")
    m <- dplyr::inner_join(dm, cv, by = "Admin2")
    m <- m[is.finite(m$svy_prev) & is.finite(m$mu) & is.finite(m$sigma) &
             is.finite(m$lon), , drop = FALSE]
    if (nrow(m) < 25) next

    vars <- setdiff(names(m), c("Admin2", "Admin1", "svy_prev", "n_svy", "n_cont",
                                "mu", "sigma", "lon", "lat", ".w"))
    vars <- vars[vapply(vars, function(v) is.numeric(m[[v]]) &&
                          sum(is.finite(m[[v]])) > 2 &&
                          stats::sd(m[[v]], na.rm = TRUE) > 1e-8, logical(1))]
    cur <- intersect(curated_vars(vars), vars)
    n <- nrow(m)

    reps <- list()
    for (rep_i in 1:10) {
      set.seed(950L + rep_i)
      fold <- sample(rep_len(1:5, n))
      P <- matrix(NA_real_, n, 4,
                  dimnames = list(NULL, c("indicator", "mean_pooled_sd",
                                          "loc_scale", "loc_scale_emp")))
      for (f in 1:5) {
        te <- which(fold == f); tr <- setdiff(seq_len(n), te)
        pp <- prep_X(m[tr, , drop = FALSE], m[te, , drop = FALSE], cur)
        Xtr <- cbind(pp$Xtr, lon = m$lon[tr], lat = m$lat[tr])
        Xte <- cbind(pp$Xte, lon = m$lon[te], lat = m$lat[te])
        w <- pmax(m$n_svy[tr], 1)
        rf <- function(target) {
          g <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = target, Xtr,
                                                                check.names = TRUE),
                                       num.trees = 800, min.node.size = 5,
                                       case.weights = w, seed = rep_i),
                        error = function(e) NULL)
          if (is.null(g)) rep(mean(target), length(te)) else
            stats::predict(g, data = data.frame(Xte, check.names = TRUE))$predictions
        }
        P[te, "indicator"] <- .ilogit(rf(.logit(m$svy_prev[tr])))
        mu_hat <- rf(m$mu[tr])
        # sigma must stay positive; learn it in logs
        sg_hat <- exp(rf(log(pmax(m$sigma[tr], 1e-6))))
        s_pool <- stats::median(m$sigma[tr], na.rm = TRUE)
        P[te, "mean_pooled_sd"] <- conv_normal(mu_hat, s_pool, cut, dir_less)
        P[te, "loc_scale"]      <- conv_normal(mu_hat, sg_hat, cut, dir_less)
        P[te, "loc_scale_emp"]  <- conv_empirical(mu_hat, sg_hat, cut, dir_less, zp)
      }
      reps[[rep_i]] <- as.data.frame(lapply(colnames(P), function(k) {
        p <- P[, k]
        data.frame(rho = suppressWarnings(cor(p, m$svy_prev, method = "spearman")),
                   r = suppressWarnings(cor(p, m$svy_prev)),
                   mae = mean(abs(p - m$svy_prev)) * 100,
                   bias = mean(p - m$svy_prev) * 100)
      }) |> setNames(colnames(P)) |> bind_rows(.id = "route"))
    }
    a <- bind_rows(reps) |> group_by(route) |>
      summarise(rho = round(mean(rho, na.rm = TRUE), 3),
                r = round(mean(r, na.rm = TRUE), 3),
                mae_pp = round(mean(mae, na.rm = TRUE), 2),
                bias_pp = round(mean(bias, na.rm = TRUE), 2), .groups = "drop")
    a$country <- ctry; a$outcome <- oc$tag; a$n_areas <- n
    rows[[paste(ctry, oc$tag)]] <- a
    message(sprintf("  %-11s %-13s n=%3d  %s", ctry, oc$tag, n,
                    paste(sprintf("%s=%.3f", a$route, a$rho), collapse = "  ")))
  }
}

res <- bind_rows(rows)
write.csv(res, "sandbox_parsimony/out/biomarker_distribution.csv", row.names = FALSE)

cat("\n=== Routes to an Admin-2 prevalence, same folds and learner ===\n")
print(as.data.frame(res |> group_by(route) |>
  summarise(cells = n(), rho = round(mean(rho), 3), r = round(mean(r), 3),
            mae_pp = round(mean(mae_pp), 2),
            abs_bias_pp = round(mean(abs(bias_pp)), 2), .groups = "drop") |>
  arrange(desc(rho))), row.names = FALSE)

cat("\n=== paired against the indicator route ===\n")
w <- res |> select(country, outcome, route, rho, mae_pp) |>
  tidyr::pivot_wider(names_from = route, values_from = c(rho, mae_pp))
for (v in c("mean_pooled_sd", "loc_scale", "loc_scale_emp")) {
  d <- w[[paste0("rho_", v)]] - w$rho_indicator
  dm <- w$mae_pp_indicator - w[[paste0("mae_pp_", v)]]
  cat(sprintf("  %-15s rho %+0.3f (SE %.3f, t=%.2f, %d/%d cells) | MAE %+0.2f pp (%d/%d)\n",
              v, mean(d), sd(d)/sqrt(length(d)), mean(d)/(sd(d)/sqrt(length(d))),
              sum(d > 0), length(d), mean(dm), sum(dm > 0), length(dm)))
}
message("\nSaved -> sandbox_parsimony/out/biomarker_distribution.csv")
