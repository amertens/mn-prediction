# =============================================================================
# sandbox_parsimony/R/11_loco_headline.R
#
# One comparable LOCO table on exactly the configuration
# results/tables/benchmarks_all.csv uses: 4 outcomes x 4 held-out countries
# (Gambia, Ghana, Sierra Leone, Malawi), 4-country covariate intersection.
#
# Tanzania is excluded here on purpose. It is 4,000 km from the West African
# three, so a thin-plate spline trained on the others can only extrapolate into
# it, and its incomplete GEE extraction shrinks the shared covariate pool from
# ~366 to ~237 columns (see 06_pool_composition.R). Including it is a separate
# question from "which model transports best", so it is held out of this table.
#
# Contenders
#   PROD enet30 / ridge_all   the production AREA_TRANSPORT_RECIPE and its
#                             documented max-performance alternative
#   spatial_tps               reproduction of fit_predict_spatial_coords, the
#                             honest leader of the existing benchmark table
#   centered_own *            within-country centering, done on the held-out
#                             country own covariate means
#   zscore *                  trained on the within-country z-score of logit
#                             prevalence, so between-country level and spread
#                             are structurally unlearnable
#   spatial + zscore          the spline and the centered covariate block together
#
# Both level anchors are reported: train_mean (what production does) and
# oracle_national (the held-out national prevalence supplied from outside).
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
source("sandbox_parsimony/R/09_loco_fns.R")
suppressPackageStartupMessages(library(mgcv))

cov_list <- readRDS("sandbox_parsimony/out/cov_list.rds")
FOUR <- c("gambia", "ghana", "sierraleone", "malawi")
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_vitA")

f_decorr20 <- function(X, y, v) decorr_reps(X, v, k = 20L)
f_screen30 <- function(X, y, v) screen_topK(X, y, v, 30L)
f_curated  <- function(X, y, v) intersect(curated_vars(v), v)

VARIANTS <- list(
  list(name = "PROD enet30",             kind = "lin", cen = "none", alpha = .5, feat = f_screen30),
  list(name = "PROD ridge_all",          kind = "lin", cen = "none", alpha = 0,  feat = NULL),
  list(name = "spatial_tps",             kind = "sp",  feat = NULL),
  list(name = "spatial + decorr20",      kind = "sp",  feat = f_decorr20),
  list(name = "centered_own enet30",     kind = "lin", cen = "own",  alpha = .5, feat = f_screen30),
  list(name = "centered_own ridge_all",  kind = "lin", cen = "own",  alpha = 0,  feat = NULL),
  list(name = "centered_own decorr20",   kind = "lin", cen = "own",  alpha = .5, feat = f_decorr20),
  list(name = "zscore ridge_all",        kind = "z",   alpha = 0,  feat = NULL),
  list(name = "zscore screen30",         kind = "z",   alpha = .5, feat = f_screen30),
  list(name = "zscore decorr20",         kind = "z",   alpha = .5, feat = f_decorr20),
  list(name = "zscore curated16",        kind = "z",   alpha = .5, feat = f_curated)
)

res <- list()
for (oc in OUTCOMES) {
  P <- assemble_outcome(oc, cov_list[FOUR]); if (is.null(P)) next
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  vars <- P$predictors
  message(sprintf("\n### %s  (%d areas, %d predictors)", oc, nrow(d), length(vars)))
  for (ho in FOUR) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    for (v in VARIANTS) for (anch in c("train_mean", "oracle_national")) {
      pr <- tryCatch(switch(v$kind,
        lin = loco_one(tr, te, vars, v$feat, v$alpha, v$cen, anch),
        sp  = loco_spatial(tr, te, vars, v$feat, anchor = anch),
        z   = loco_zscore(tr, te, vars, v$feat, alpha = v$alpha, anchor = anch)),
        error = function(e) NULL)
      if (is.null(pr)) next
      m <- loco_metrics(pr, te); if (is.null(m)) next
      m$outcome <- oc; m$held_out <- ho; m$variant <- v$name; m$anchor <- anch
      res[[paste(oc, ho, v$name, anch)]] <- m
    }
  }
}

out <- dplyr::bind_rows(res)
write.csv(out, "sandbox_parsimony/out/loco_headline.csv", row.names = FALSE)

cat("\n=== LOCO, 4 outcomes x 4 held-out countries (same cells as benchmarks_all.csv) ===\n")
s <- out |> dplyr::group_by(variant, anchor) |>
  dplyr::summarise(cells = dplyr::n(),
                   rho = round(mean(spearman, na.rm = TRUE), 3),
                   rho_sd = round(stats::sd(spearman, na.rm = TRUE), 3),
                   r = round(mean(pearson, na.rm = TRUE), 3),
                   rmse = round(mean(rmse_pp, na.rm = TRUE), 1),
                   pat_rmse = round(mean(pattern_rmse_pp, na.rm = TRUE), 1),
                   abs_bias = round(mean(abs(level_bias_pp), na.rm = TRUE), 1),
                   .groups = "drop") |>
  dplyr::arrange(anchor, dplyr::desc(rho))
print(as.data.frame(s), row.names = FALSE)

cat("\n=== per outcome, Spearman (train_mean anchor) ===\n")
t2 <- out |> dplyr::filter(anchor == "train_mean") |>
  dplyr::group_by(outcome, variant) |>
  dplyr::summarise(rho = round(mean(spearman, na.rm = TRUE), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = outcome, values_from = rho)
print(as.data.frame(t2), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/loco_headline.csv")
