# =============================================================================
# sandbox_parsimony/R/34_rescore_spatial_plus_soil.R
#
# STEP 4: re-score the one benchmark method whose features were chosen on the
# metric it is reported against.
#
# FINDINGS.md section 3: `spatial_plus_soil` leads results/tables/benchmarks_all.csv
# at LOCO rho = 0.305. Its eight locked SoilGrids variables come from
# archive/sandbox/13_univariate_loco_ranking.R, which ranked every candidate by
# LOCO Pearson r across the SAME four outcomes and the SAME four held-out
# countries the method is then scored on. That is selection on the evaluation
# metric, so 0.305 is not an out-of-sample number.
#
# The fix is to move the selection inside the loop: for each held-out country,
# rank the candidates on the TRAINING countries only, take the top 8, and fit
# the same spatial + soil model with those. Nothing else changes.
#
# Reported alongside:
#   locked_published   the eight locked variables, as published (optimistic)
#   selected_in_fold   top 8 chosen on training countries only (honest)
#   spatial_only       no soil variables at all, the untainted comparator
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
suppressPackageStartupMessages({library(dplyr); library(mgcv)})

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_vitA")
FOUR <- c("gambia", "ghana", "sierraleone", "malawi")

# The eight published, locked features (R/benchmark_models.R)
LOCKED <- c("gee_soilaluminium_annual_min", "gee_soilaluminium_stdev_20_50",
            "gee_soilcalcium_stdev_0_20", "gee_soilaluminium_stdev_0_20",
            "gee_soilmagnesium_stdev_0_20", "gee_soilcalcium_stdev_20_50",
            "gee_soiltotalcarbon_mean_0_20", "gee_soilzinc_mean_20_50")

#' The published method: thin-plate spline on lon/lat + a few linear extras,
#' fitted on the raw prevalence scale (fit_predict_spatial_coords defaults to
#' "continuous", which FINDINGS section 7 showed matters a lot).
sp_soil <- function(tr, te, extras, k = 30L) {
  ok <- is.finite(tr$lon) & is.finite(tr$lat)
  if (sum(ok) < 20) return(rep(mean(tr$svy_prev), nrow(te)))
  extras <- intersect(extras, names(tr))
  d_tr <- data.frame(Y = tr$svy_prev, lon = tr$lon, lat = tr$lat)
  d_te <- data.frame(lon = te$lon, lat = te$lat)
  if (length(extras)) {
    pp <- prep_X(tr, te, extras)
    if (ncol(pp$Xtr)) {
      d_tr <- cbind(d_tr, as.data.frame(pp$Xtr))
      d_te <- cbind(d_te, as.data.frame(pp$Xte))
      extras <- pp$vars
    } else extras <- character(0)
  }
  rhs <- sprintf("s(lon, lat, k = %d, bs = 'tp')", min(k, nrow(tr) - 5))
  if (length(extras)) rhs <- paste(c(rhs, sprintf("`%s`", extras)), collapse = " + ")
  g <- tryCatch(mgcv::gam(stats::as.formula(paste("Y ~", rhs)), data = d_tr,
                          weights = pmax(tr$n_svy, 1), method = "REML"),
                error = function(e) NULL)
  if (is.null(g)) return(rep(mean(tr$svy_prev), nrow(te)))
  p <- tryCatch(as.numeric(stats::predict(g, newdata = d_te)),
                error = function(e) rep(mean(tr$svy_prev), nrow(te)))
  pmin(pmax(p, 0), 1)
}

#' Rank candidate soil variables by univariate LOCO r WITHIN the training set.
#' Mirrors 13_univariate_loco_ranking.R, but the inner leave-one-out is over
#' TRAINING countries only, so the held-out country never informs the choice.
select_top8 <- function(tr, cands) {
  cands <- intersect(cands, names(tr))
  if (length(cands) < 8) return(cands)
  ctys <- unique(tr$country)
  if (length(ctys) < 3) return(utils::head(cands, 8))
  sc <- vapply(cands, function(v) {
    pr <- rep(NA_real_, nrow(tr))
    for (ho in ctys) {
      i_te <- which(tr$country == ho); i_tr <- setdiff(seq_len(nrow(tr)), i_te)
      x_tr <- suppressWarnings(as.numeric(tr[[v]][i_tr]))
      x_te <- suppressWarnings(as.numeric(tr[[v]][i_te]))
      m <- stats::median(x_tr[is.finite(x_tr)])
      x_tr[!is.finite(x_tr)] <- m; x_te[!is.finite(x_te)] <- m
      if (!is.finite(stats::sd(x_tr)) || stats::sd(x_tr) < 1e-9) next
      fit <- tryCatch(stats::lm(tr$svy_prev[i_tr] ~ x_tr,
                                weights = pmax(tr$n_svy[i_tr], 1)),
                      error = function(e) NULL)
      if (!is.null(fit))
        pr[i_te] <- stats::coef(fit)[1] + stats::coef(fit)[2] * x_te
    }
    o <- is.finite(pr)
    if (sum(o) < 10) return(0)
    abs(suppressWarnings(stats::cor(pr[o], tr$svy_prev[o])))
  }, numeric(1))
  sc[!is.finite(sc)] <- 0
  names(sort(sc, decreasing = TRUE))[1:8]
}

rows <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  d <- P$data[P$data$country %in% FOUR, ]
  d <- d[is.finite(d$svy_prev) & is.finite(d$n_svy) & is.finite(d$lon), ]
  if (!nrow(d)) next
  soil_cands <- grep("^gee_soil", P$predictors, value = TRUE)
  for (ho in FOUR) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 5) next
    variants <- list(
      locked_published = LOCKED,
      selected_in_fold = select_top8(tr, soil_cands),
      spatial_only     = character(0))
    for (nm in names(variants)) {
      pr <- tryCatch(sp_soil(tr, te, variants[[nm]]), error = function(e) NULL)
      if (is.null(pr)) next
      m <- loco_metrics(pr, te); if (is.null(m)) next
      m$outcome <- oc; m$held_out <- ho; m$variant <- nm
      m$n_features <- length(variants[[nm]])
      rows[[paste(oc, ho, nm)]] <- m
    }
  }
}

res <- bind_rows(rows)
write.csv(res, "sandbox_parsimony/out/spatial_plus_soil_rescored.csv", row.names = FALSE)

cat("\n=== spatial_plus_soil, features locked on the test metric vs chosen in-fold ===\n")
print(as.data.frame(res |> group_by(variant) |>
  summarise(cells = n(),
            spearman = round(mean(spearman, na.rm = TRUE), 3),
            pearson = round(mean(pearson, na.rm = TRUE), 3),
            rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 1),
            .groups = "drop") |> arrange(desc(spearman))), row.names = FALSE)

w <- res |> select(outcome, held_out, variant, spearman) |>
  tidyr::pivot_wider(names_from = variant, values_from = spearman)
w <- w[stats::complete.cases(w), ]
if (nrow(w) > 2) {
  d1 <- w$locked_published - w$selected_in_fold
  cat(sprintf("\noptimism from locking the features on the test metric: %+0.3f (SE %.3f, t = %.2f, %d/%d cells)\n",
              mean(d1), sd(d1) / sqrt(nrow(w)),
              mean(d1) / (sd(d1) / sqrt(nrow(w))), sum(d1 > 0), nrow(w)))
  d2 <- w$selected_in_fold - w$spatial_only
  cat(sprintf("honest value of the soil block over spatial alone:      %+0.3f (SE %.3f, t = %.2f, %d/%d cells)\n",
              mean(d2), sd(d2) / sqrt(nrow(w)),
              mean(d2) / (sd(d2) / sqrt(nrow(w))), sum(d2 > 0), nrow(w)))
}
message("\nSaved -> sandbox_parsimony/out/spatial_plus_soil_rescored.csv")
