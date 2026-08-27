# =============================================================================
# sandbox_parsimony/R/09_loco_spatial.R
#
# LOCO, round 2: combine the three things that each helped separately.
#
# Context. results/tables/benchmarks_all.csv already shows that the spatial
# methods top the LOCO table (spatial_coords rho = 0.267) while the penalised
# regressions closest to the production recipe sit near rho = 0.09. Round 1
# (05_loco.R) showed that within-country centering lifts rank transfer from
# rho ~ 0.17 to ~ 0.22. Nobody has combined them.
#
# NB on spatial_plus_soil (rho = 0.305, the current table leader): its eight
# locked SoilGrids features were chosen as the top-8 by LOCO r in
# archive/sandbox/13_univariate_loco_ranking.R -- the SAME four outcomes and the
# SAME four held-out countries it is then scored on. That number is selection on
# the evaluation metric, so spatial_coords is the honest spatial comparator.
#
# Variants here
#   spatial_tps          thin-plate spline on lon/lat, weighted, logit scale
#                        (a reproduction of fit_predict_spatial_coords)
#   spatial + covars     TPS plus within-country-centered parsimonious covariates
#   zscore_target        train on the WITHIN-COUNTRY z-score of logit prevalence,
#                        so the fit cannot chase between-country level or spread
#                        at all; predictions are re-scaled to the held-out
#                        country using an anchor
#   country_RE           country random intercept (lme4), predicted at RE = 0
#
# Anchors: train_mean (production) vs oracle_national (held-out national
# prevalence supplied). Level and pattern are reported separately.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
suppressPackageStartupMessages(library(mgcv))

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_folate")

source("sandbox_parsimony/R/09_loco_fns.R")

f_decorr20 <- function(X, y, v) decorr_reps(X, v, k = 20L)
f_screen30 <- function(X, y, v) screen_topK(X, y, v, 30L)
f_curated  <- function(X, y, v) intersect(curated_vars(v), v)

VARIANTS <- list(
  list(name = "spatial_tps",            fn = "sp",  feat = NULL),
  list(name = "spatial + decorr20",     fn = "sp",  feat = f_decorr20),
  list(name = "spatial + curated16",    fn = "sp",  feat = f_curated),
  list(name = "zscore + screen30",      fn = "z",   feat = f_screen30),
  list(name = "zscore + decorr20",      fn = "z",   feat = f_decorr20),
  list(name = "zscore + ridge_all",     fn = "z",   feat = NULL, alpha = 0),
  list(name = "zscore + curated16",     fn = "z",   feat = f_curated)
)

res <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P) || length(P$countries) < 3) next
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  vars <- P$predictors
  for (ho in unique(d$country)) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    message(sprintf("\n== %s | held out %s (train %d, test %d) ==", oc, ho, nrow(tr), nrow(te)))
    for (v in VARIANTS) for (anch in c("train_mean", "oracle_national")) {
      pr <- tryCatch(
        if (v$fn == "sp") loco_spatial(tr, te, vars, v$feat, anchor = anch)
        else loco_zscore(tr, te, vars, v$feat, alpha = v$alpha %||% 0.5, anchor = anch),
        error = function(e) NULL)
      if (is.null(pr)) next
      m <- loco_metrics(pr, te); if (is.null(m)) next
      m$outcome <- oc; m$held_out <- ho; m$variant <- v$name; m$anchor <- anch
      res[[paste(oc, ho, v$name, anch)]] <- m
      if (anch == "train_mean")
        message(sprintf("   %-22s r=%+.3f rho=%+.3f rmse=%5.1f pat_rmse=%5.1f bias=%+6.1f",
                        v$name, m$pearson, m$spearman, m$rmse_pp,
                        m$pattern_rmse_pp, m$level_bias_pp))
    }
  }
}

out <- dplyr::bind_rows(res)
write.csv(out, "sandbox_parsimony/out/loco_spatial.csv", row.names = FALSE)
cat("\n=== mean over 16 outcome x held-out cells ===\n")
print(as.data.frame(
  out |> dplyr::group_by(variant, anchor) |>
    dplyr::summarise(cells = dplyr::n(),
                     rho = round(mean(spearman, na.rm = TRUE), 3),
                     r = round(mean(pearson, na.rm = TRUE), 3),
                     rmse = round(mean(rmse_pp, na.rm = TRUE), 1),
                     pat_rmse = round(mean(pattern_rmse_pp, na.rm = TRUE), 1),
                     abs_bias = round(mean(abs(level_bias_pp), na.rm = TRUE), 1),
                     .groups = "drop") |>
    dplyr::arrange(anchor, dplyr::desc(rho))), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/loco_spatial.csv")
