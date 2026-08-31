# =============================================================================
# scripts/accuracy_impact/ws4b_skill_curve.R
#
# WS4b. The reliability-skill curve: do the micronutrient cells sit where their
# reliability says they should?
#
# THE CLAIM THIS QUANTIFIES
# -------------------------
# Section 11 claim 2 says the constraint is the TARGET, not the predictors, and
# supports it with a positive control: earth observation predicts district
# education at r 0.48 to 0.71. Stage 0 could not locate a committed table
# containing that number anywhere under results/, sensitivity/,
# sandbox_parsimony/ or docs/, so it is recorded in the claims register as
# unsourced. This script replaces it with something stronger than a single
# control: the same pipeline, unchanged, run on a range of DHS indicators whose
# reliabilities differ, so the relationship between reliability and achieved
# skill can be seen rather than asserted.
#
# If the micronutrient cells fall on the same curve as education, stunting,
# vaccination and the rest, then their low skill is explained by their low
# reliability and the predictors are not at fault. If they fall BELOW the curve,
# something specific to micronutrient targets is also wrong.
#
# THE RELIABILITY USED HERE, AND ITS LIMIT
# ----------------------------------------
# WS1 established that the analytic ceiling with a design effect fixed at 1.5 is
# biased low by a mean of 0.161. That criticism does not apply here: the DHS
# direct estimates carry their OWN design-based variance, so
#
#   lambda = 1 - mean(direct.var) / Var(direct.est)
#
# uses no assumed design effect at all.
#
# It has a different limit, and it is the same single-PSU problem. A district
# holding one cluster cannot form a between-cluster variance, and DHS returns a
# variance near zero for it (measured: 2.3e-34 for Ghana's Asunafo South). Those
# districts make the mean sampling variance too small and the reliability too
# high. Districts whose reported variance is implausibly small are therefore
# excluded and the count is reported, so the direction of any residual bias is
# on the record.
#
#   Rscript scripts/accuracy_impact/ws4b_skill_curve.R
# -> results/tables/reliability_skill_curve.csv
# -> results/figures/reliability_skill_curve.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

SEED <- 20260907L; K_SCREEN <- 20L; MIN_DISTRICTS <- 12L
set.seed(SEED)
TDIR <- here("results", "tables"); FDIR <- here("results", "figures")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS_ALL <- setdiff(names(H), c("country", "Admin1", "Admin2"))

# THE DHS COVARIATES MUST GO, AND NOT ONLY FOR THE DHS TARGETS.
#
# Measured: 97 of the harmonised Admin-2 covariates are DHS-derived aggregates,
# including dhs_AN_ANEM_W_ANY. Using them to predict the DHS indicator
# AN_ANEM_W_ANY is the outcome predicting itself, and the first run of this
# script duly returned r 0.946 for Gambia anaemia with r_share 1.05. That is the
# same defect class Section 7 records three times over, arriving here.
#
# They are dropped for BOTH families, not only for the DHS targets. Dropping
# them only where they leak would leave the micronutrient cells fitted on 294
# predictors and the DHS cells on 197, so the two families would no longer be
# on the same curve and the comparison the workstream exists to make would be
# between different models rather than different targets.
COVS <- grep("^dhs", COVS_ALL, value = TRUE, invert = TRUE)
cat(sprintf("[ws4b] covariates: %d total, %d DHS-derived removed, %d used\n",
            length(COVS_ALL), length(COVS_ALL) - length(COVS), length(COVS)))

DHS_FILES <- list(
  Gambia        = "Gambia_2019",
  Ghana         = "Ghana_2014",
  Malawi        = "Malawi_2015",
  `Sierra Leone` = "Sierra Leone_2013"
)

rows <- list()
for (cn in names(DHS_FILES)) {
  stem <- DHS_FILES[[cn]]
  p <- here("data", "DHS", "clean", paste0(stem, "_dhs_admin2_direct.rds"))
  if (!file.exists(p)) { cat("[skip]", cn, "no direct file\n"); next }
  dd <- readRDS(p)
  hc <- H[H$country == gsub(" ", "", cn) | H$country == cn, , drop = FALSE]
  if (!nrow(hc)) { cat("[skip]", cn, "no covariates\n"); next }

  for (ind in sort(unique(dd$indicator))) {
    z <- dd[dd$indicator == ind, , drop = FALSE]
    z <- z[is.finite(z$direct.est) & is.finite(z$direct.var), , drop = FALSE]
    if (nrow(z) < MIN_DISTRICTS) next

    # Exclude districts whose reported design variance is degenerate: a single
    # PSU cannot produce a between-cluster variance and DHS returns ~0.
    v_floor <- 1e-8
    ok <- z$direct.var > v_floor
    n_degen <- sum(!ok)
    zz <- z[ok, , drop = FALSE]
    if (nrow(zz) < MIN_DISTRICTS) next

    v_obs <- stats::var(zz$direct.est)
    if (!is.finite(v_obs) || v_obs <= 0) next
    lambda <- 1 - mean(zz$direct.var) / v_obs
    lambda_t <- min(max(lambda, 0), 1)
    r_max <- sqrt(lambda_t)

    # Model skill: the same area-level learner and the same leave-one-region-out
    # blocking the micronutrient cells get, with the prescreen inside the fold.
    m <- dplyr::inner_join(
      data.frame(Admin1 = trimws(as.character(zz$admin1.name)),
                 Admin2 = trimws(as.character(zz$admin2.name)),
                 y = zz$direct.est, stringsAsFactors = FALSE),
      hc, by = admin2_join_by(
        data.frame(Admin1 = zz$admin1.name, Admin2 = zz$admin2.name), hc))
    m <- m[is.finite(m$y), , drop = FALSE]
    if (nrow(m) < MIN_DISTRICTS || dplyr::n_distinct(m$Admin1) < 3) next
    X <- as.matrix(m[, COVS, drop = FALSE])
    oof <- rep(NA_real_, nrow(m))
    for (rg in unique(m$Admin1)) {
      i <- which(m$Admin1 == rg)
      fit <- .ds_fit(X[-i, , drop = FALSE], m$y[-i], k_screen = K_SCREEN)
      pr <- .ds_predict(fit, X[i, , drop = FALSE])
      if (!is.null(pr)) oof[i] <- pr
    }
    fin <- is.finite(oof) & is.finite(m$y)
    if (sum(fin) < MIN_DISTRICTS || stats::sd(m$y[fin]) == 0) next
    r <- suppressWarnings(stats::cor(m$y[fin], oof[fin]))

    rows[[length(rows) + 1L]] <- data.frame(
      country = cn, target = ind, family = "DHS indicator",
      n_districts = sum(fin), n_degenerate_var = n_degen,
      mean_est = round(mean(zz$direct.est), 4),
      lambda = round(lambda, 4), r_max = round(r_max, 4),
      r_oof = round(r, 4),
      r_share = if (r_max > 0.05) round(r / r_max, 3) else NA_real_,
      stringsAsFactors = FALSE)
    cat(sprintf("  %-13s %-16s districts=%3d  r_max=%.3f  r=%+.3f\n",
                cn, ind, sum(fin), r_max, r))
  }
}
dhs <- dplyr::bind_rows(rows)

# The micronutrient cells, overlaid. Skill is the unanchored leave-one-region-out
# arm, which is the like-for-like comparison; reliability is the empirical
# split-half ceiling WS1a measured, NOT the analytic one.
# The micronutrient cells, fitted HERE on the same reduced covariate set and the
# same learner and folds, rather than read from anchor_controls.csv, which used
# all 294 covariates including the 97 DHS ones. Reading that table would put the
# two families on different models and make the curve uninterpretable.
mn <- NULL
re <- file.path(TDIR, "reliability_empirical.csv")
if (file.exists(re)) {
  rel <- read.csv(re, stringsAsFactors = FALSE); rel <- rel[rel$scheme == "within", ]
  kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
  cfgs <- get_country_configs(); mrows <- list()
  for (cn2 in names(cfgs)) {
    cc <- cfgs[[cn2]]
    hc <- H[H$country == cn2, , drop = FALSE]
    if (!nrow(hc)) next
    for (ocn in names(cc$outcomes)) {
      sv <- tryCatch(targets::tar_read_raw(
        paste0("svy_admin2_", tolower(cn2), "_", ocn),
        store = here("_targets_full")), error = function(e) NULL)
      if (is.null(sv) || nrow(sv) < MIN_DISTRICTS) next
      m <- dplyr::inner_join(
        sv[, intersect(c("Admin1","Admin2","svy_prev"), names(sv)), drop = FALSE],
        hc, by = admin2_join_by(sv, hc))
      m <- m[is.finite(m$svy_prev), , drop = FALSE]
      if (nrow(m) < MIN_DISTRICTS || dplyr::n_distinct(m$Admin1) < 3) next
      X <- as.matrix(m[, COVS, drop = FALSE])
      oof <- rep(NA_real_, nrow(m))
      for (rg in unique(m$Admin1)) {
        i <- which(m$Admin1 == rg)
        fit <- .ds_fit(X[-i, , drop = FALSE], m$svy_prev[-i], k_screen = K_SCREEN)
        pr <- .ds_predict(fit, X[i, , drop = FALSE])
        if (!is.null(pr)) oof[i] <- pr
      }
      fin <- is.finite(oof)
      if (sum(fin) < MIN_DISTRICTS || stats::sd(m$svy_prev[fin]) == 0) next
      rm_e <- rel$r_max_emp[match(paste(kk(cc$country), ocn),
                                  paste(kk(rel$country), rel$outcome))]
      r <- suppressWarnings(stats::cor(m$svy_prev[fin], oof[fin]))
      mrows[[length(mrows) + 1L]] <- data.frame(
        country = cc$country, target = ocn, family = "micronutrient",
        n_districts = sum(fin), n_degenerate_var = NA_integer_,
        mean_est = round(mean(m$svy_prev[fin]), 4),
        lambda = round(rm_e^2, 4), r_max = round(rm_e, 4), r_oof = round(r, 4),
        r_share = if (is.finite(rm_e) && rm_e > 0.05) round(r / rm_e, 3) else NA_real_,
        stringsAsFactors = FALSE)
      cat(sprintf("  %-13s %-16s districts=%3d  r_max=%.3f  r=%+.3f\n",
                  cc$country, ocn, sum(fin), rm_e, r))
    }
  }
  if (length(mrows)) mn <- dplyr::bind_rows(mrows)
}
out <- dplyr::bind_rows(dhs, mn)
readr::write_csv(out, file.path(TDIR, "reliability_skill_curve.csv"))

cat("\n=== WS4b: skill against reliability ===\n")
print(as.data.frame(out |> group_by(family) |>
  summarise(cells = dplyr::n(),
            med_r_max = round(stats::median(r_max, na.rm = TRUE), 3),
            med_r = round(stats::median(r_oof, na.rm = TRUE), 3),
            med_r_share = round(stats::median(r_share, na.rm = TRUE), 3),
            .groups = "drop")), row.names = FALSE)

fin <- is.finite(out$r_max) & is.finite(out$r_oof)
if (sum(fin) > 5) {
  cat(sprintf("\ncorrelation between r_max and achieved r across all %d targets: %.3f\n",
              sum(fin), stats::cor(out$r_max[fin], out$r_oof[fin])))
  d <- out[fin, ]
  fitl <- stats::lm(r_oof ~ r_max, data = d[d$family == "DHS indicator", ])
  pred_mn <- stats::predict(fitl, newdata = d[d$family == "micronutrient", ])
  resid_mn <- d$r_oof[d$family == "micronutrient"] - pred_mn
  cat(sprintf("micronutrient cells sit %+.3f from the DHS-indicator line on average\n",
              mean(resid_mn, na.rm = TRUE)))
  cat(sprintf("  (below the line in %d of %d cells)\n",
              sum(resid_mn < 0, na.rm = TRUE), sum(is.finite(resid_mn))))
}

# The figure. Base graphics so the script has no ggplot dependency.
png(file.path(FDIR, "reliability_skill_curve.png"), width = 1100, height = 850, res = 130)
op <- par(mar = c(4.5, 4.5, 3, 1))
d <- out[is.finite(out$r_max) & is.finite(out$r_oof), ]
isr <- d$family == "micronutrient"
plot(d$r_max, d$r_oof, type = "n", xlim = c(0, 1), ylim = range(c(-0.6, 1)),
     xlab = "Reliability ceiling  r_max",
     ylab = "Achieved out-of-fold correlation",
     main = "Skill against reliability, DHS indicators and micronutrient cells")
abline(0, 1, col = "grey70", lty = 2)
abline(h = 0, col = "grey85")
points(d$r_max[!isr], d$r_oof[!isr], pch = 1, col = "grey35")
points(d$r_max[isr],  d$r_oof[isr],  pch = 19, col = "firebrick")
if (sum(!isr) > 3) abline(stats::lm(r_oof ~ r_max, data = d[!isr, ]), col = "grey35")
legend("topleft", bty = "n",
       legend = c("DHS indicator", "micronutrient cell", "y = x (ceiling)"),
       pch = c(1, 19, NA), lty = c(NA, NA, 2),
       col = c("grey35", "firebrick", "grey70"))
par(op); dev.off()
cat(sprintf("\n-> %s\n", file.path("results", "figures", "reliability_skill_curve.png")))
