# =============================================================================
# sandbox_parsimony/R/05_loco.R
#
# Leave-one-country-out transportability bake-off.
#
# The production LOCO (R/transportability_area.R) reports one number that mixes
# two very different questions:
#   (1) LEVEL   -- what is the held-out country national prevalence?
#   (2) PATTERN -- which districts inside it are worst?
# (1) is dominated by cross-survey biomarker/assay offsets and is essentially
# untransportable; (2) is the product a program actually uses. Pooling them into
# one Pearson r on the raw prevalence scale lets a level miss of 40 pp swamp a
# perfectly good spatial ranking (Gambia child_iron: r = 0.52, bias = -44 pp).
#
# Variants tested
#   prod_enet30      reproduction of AREA_TRANSPORT_RECIPE
#   prod_ridge_all   the documented max-performance alternative
#   centered_train   the centering as currently coded: train centered per
#                    country, TEST centered on pooled-train means
#   centered_own     centering done properly: the held-out country covariates
#                    are centered on ITS OWN means. This uses no held-out
#                    OUTCOME, only its covariates, which are observable
#                    everywhere without a survey -- that is the whole premise
#                    of the method, so it is not leakage.
#   centered_own_*   the same, on parsimonious predictor sets
#
# For every variant we report BOTH:
#   pearson/spearman on the raw scale (level + pattern), and
#   pattern_rmse_pp after removing the level from both sides, plus the level
#   bias separately. (Pearson r is level-invariant, so there is no distinct
#   "pattern r" -- the level shows up in RMSE and bias, not in r.)
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_folate")

source("sandbox_parsimony/R/05a_loco_fns.R")

VARIANTS <- list(
  list(name = "PROD enet30",            alpha = .5, centering = "none",
       feat = function(X, y, v) screen_topK(X, y, v, 30L)),
  list(name = "PROD ridge_all",         alpha = 0,  centering = "none", feat = NULL),
  list(name = "centered_train enet30",  alpha = .5, centering = "train",
       feat = function(X, y, v) screen_topK(X, y, v, 30L)),
  list(name = "centered_own enet30",    alpha = .5, centering = "own",
       feat = function(X, y, v) screen_topK(X, y, v, 30L)),
  list(name = "centered_own ridge_all", alpha = 0,  centering = "own", feat = NULL),
  list(name = "centered_own decorr20",  alpha = .5, centering = "own",
       feat = function(X, y, v) decorr_reps(X, v, k = 20L)),
  list(name = "centered_own d40s12",    alpha = .5, centering = "own",
       feat = function(X, y, v) decorr_then_screen(X, y, v, 40L, 12L)),
  list(name = "centered_own curated16", alpha = .5, centering = "own",
       feat = function(X, y, v) intersect(curated_vars(v), v)),
  list(name = "curated16 uncentered",   alpha = .5, centering = "none",
       feat = function(X, y, v) intersect(curated_vars(v), v))
)

res <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P) || length(P$countries) < 3) next
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  vars <- P$predictors
  for (ho in unique(d$country)) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    message(sprintf("\n== %s | held out %s (train %d, test %d, p=%d) ==",
                    oc, ho, nrow(tr), nrow(te), length(vars)))
    for (v in VARIANTS) {
      for (anch in c("train_mean", "oracle_national")) {
        pr <- tryCatch(loco_one(tr, te, vars, v$feat, v$alpha, v$centering, anch),
                       error = function(e) NULL)
        if (is.null(pr)) next
        m <- loco_metrics(pr, te)
        if (is.null(m)) next
        m$outcome <- oc; m$held_out <- ho; m$variant <- v$name; m$anchor <- anch
        res[[paste(oc, ho, v$name, anch)]] <- m
        if (anch == "train_mean")
          message(sprintf("   %-24s r=%+.3f rho=%+.3f  rmse=%5.1f  pat_rmse=%5.1f  level_bias=%+6.1f",
                          v$name, m$pearson, m$spearman, m$rmse_pp,
                          m$pattern_rmse_pp, m$level_bias_pp))
      }
    }
  }
}

out <- dplyr::bind_rows(res) |>
  dplyr::select(outcome, held_out, variant, anchor, n_test, r_max, pearson,
                spearman, rmse_pp, pattern_rmse_pp, level_bias_pp, calib)
write.csv(out, "sandbox_parsimony/out/loco_bakeoff.csv", row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/loco_bakeoff.csv")
