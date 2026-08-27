# =============================================================================
# sandbox_parsimony/R/06_pool_composition.R
#
# What does adding a country with an incomplete GEE extraction cost everyone?
#
# The pooled / LOCO models take the INTERSECTION of covariate names across
# countries. Tanzania has only 20 of the 35 GEE families the other four share.
# Adding it to the child_vitA / women_vitA pool therefore deletes, for every
# country:
#
#   accessibility (travel time to cities), lst (land surface temperature),
#   terraclimate (precip / PET / PDSI / soil moisture / VPD), productivity
#   (GPP / NPP), dailyevi, lai8days, landcovertype + landcoverlayers (cropland,
#   tree, urban cover fractions), esa worldcereal, wapor, ghsbuilts, ghspop,
#   fldas, aerosoloptical, atmosphere
#
# -- essentially every climate, land-cover, vegetation-dynamics and market-
# access predictor, leaving soil chemistry plus a few static urbanicity layers.
#
# This script measures the cost directly: LOCO over the SAME four countries
# (Gambia, Ghana, Sierra Leone, Malawi) and the SAME held-out sets, run once on
# the 5-country (Tanzania-constrained) covariate pool and once on the richer
# 4-country pool.
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")

cov_list   <- readRDS("sandbox_parsimony/out/cov_list.rds")
pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
FOUR <- c("gambia", "ghana", "sierraleone", "malawi")

VARIANTS <- list(
  list(name = "PROD enet30",           alpha = .5, cen = "none",
       feat = function(X, y, v) screen_topK(X, y, v, 30L)),
  list(name = "PROD ridge_all",        alpha = 0,  cen = "none", feat = NULL),
  list(name = "centered_own enet30",   alpha = .5, cen = "own",
       feat = function(X, y, v) screen_topK(X, y, v, 30L)),
  list(name = "centered_own ridge_all",alpha = 0,  cen = "own", feat = NULL),
  list(name = "centered_own decorr20", alpha = .5, cen = "own",
       feat = function(X, y, v) decorr_reps(X, v, k = 20L))
)

res <- list()
for (oc in c("child_vitA", "women_vitA")) {
  # pool A: 5-country intersection, then restrict ROWS to the four countries
  P5 <- pooled_all[[oc]]; if (is.null(P5)) next
  dA <- P5$data[P5$data$country %in% FOUR, ]
  varsA <- P5$predictors

  # pool B: intersection over the four countries only
  P4 <- assemble_outcome(oc, cov_list[FOUR])
  if (is.null(P4)) next
  dB <- P4$data; varsB <- P4$predictors

  message(sprintf("\n### %s  |  5-country pool: %d predictors   4-country pool: %d predictors",
                  oc, length(varsA), length(varsB)))

  for (pool in c("with_tanzania_pool", "four_country_pool")) {
    d    <- if (pool == "with_tanzania_pool") dA else dB
    vars <- if (pool == "with_tanzania_pool") varsA else varsB
    d <- d[is.finite(d$svy_prev) & is.finite(d$n_svy), ]
    for (ho in FOUR) {
      tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
      if (nrow(tr) < 30 || nrow(te) < 8) next
      for (v in VARIANTS) {
        pr <- tryCatch(loco_one(tr, te, vars, v$feat, v$alpha, v$cen, "train_mean"),
                       error = function(e) NULL)
        if (is.null(pr)) next
        m <- loco_metrics(pr, te); if (is.null(m)) next
        m$outcome <- oc; m$pool <- pool; m$held_out <- ho
        m$variant <- v$name; m$n_pred_pool <- length(vars)
        res[[paste(oc, pool, ho, v$name)]] <- m
      }
    }
  }
}

out <- dplyr::bind_rows(res)
write.csv(out, "sandbox_parsimony/out/pool_composition.csv", row.names = FALSE)

cat("\n=== LOCO on the SAME four countries, two covariate pools ===\n")
s <- out |> dplyr::group_by(outcome, pool, n_pred_pool) |>
  dplyr::summarise(cells = dplyr::n(),
                   mean_r   = round(mean(pearson, na.rm = TRUE), 3),
                   mean_rho = round(mean(spearman, na.rm = TRUE), 3),
                   mean_rmse = round(mean(rmse_pp, na.rm = TRUE), 1),
                   .groups = "drop")
print(as.data.frame(s), row.names = FALSE)

cat("\n=== by variant ===\n")
s2 <- out |> dplyr::group_by(outcome, variant, pool) |>
  dplyr::summarise(rho = round(mean(spearman, na.rm = TRUE), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = pool, values_from = rho) |>
  dplyr::mutate(gain_from_dropping_TZ_pool =
                  round(four_country_pool - with_tanzania_pool, 3))
print(as.data.frame(s2), row.names = FALSE)
