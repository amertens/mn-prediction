# =============================================================================
# sandbox_parsimony/R/12_scale_sweep.R
#
# Does AREA_TRANSPORT_RECIPE$scale = "logit" help or hurt?
#
# Reproducing fit_predict_spatial_coords exactly (see the reproduction check in
# 11/scratch) turned up a large, unadvertised effect: the same thin-plate spline
# scores LOCO rho = 0.274 fitted on the raw prevalence scale and rho = 0.106
# fitted on logit(prevalence). The production area-transport recipe fixes
# scale = "logit"; the benchmark suite runs some methods on "continuous" and
# some on "logit" but never isolates the choice.
#
# Mechanism to suspect: on the logit scale the clamp at 0.005 / 0.995 turns the
# zero-prevalence districts -- and at a median n_svy of 6-14 there are many --
# into extreme outliers at logit = -5.3 that the fit then chases. On the raw
# scale those districts are just 0.
#
# This sweeps scale x model x eval so the choice can be made on evidence.
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
suppressPackageStartupMessages(library(mgcv))

cov_list <- readRDS("sandbox_parsimony/out/cov_list.rds")
FOUR <- c("gambia", "ghana", "sierraleone", "malawi")
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_vitA")

#' glmnet LOCO fit with the response scale as an explicit argument
loco_scaled <- function(tr, te, vars, feature_fn, alpha, scale = "logit",
                        centering = "none", seed = 12345L) {
  pp <- prep_X(tr, te, vars); Xtr <- pp$Xtr; Xte <- pp$Xte; vv <- pp$vars
  y <- if (scale == "logit") .logit(tr$svy_prev) else tr$svy_prev
  yf <- y; level <- mean(y)
  if (centering == "own") {
    Xtr <- center_by(Xtr, tr$country)
    Xte <- sweep(Xte, 2, colMeans(Xte), "-")
    for (g in unique(tr$country)) { i <- tr$country == g; yf[i] <- y[i] - mean(y[i]) }
  }
  if (!is.null(feature_fn) && length(vv) > 1) {
    sel <- tryCatch(feature_fn(Xtr, yf, vv), error = function(e) vv)
    if (length(sel) >= 1) { Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE] }
  }
  if (ncol(Xtr) < 2) return(rep(mean(tr$svy_prev), nrow(te)))
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, yf, alpha = alpha,
                                   weights = pmax(tr$n_svy, 1), nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(tr$svy_prev), nrow(te)))
  p <- as.numeric(stats::predict(cv, newx = Xte, s = "lambda.min"))
  if (centering == "own") p <- p + level
  if (scale == "logit") p <- .ilogit(p)
  pmin(pmax(p, 0), 1)
}

loco_tps <- function(tr, te, scale = "logit", k = 30L) {
  Y <- if (scale == "logit") .logit(tr$svy_prev) else tr$svy_prev
  g <- tryCatch(mgcv::gam(Y ~ s(lon, lat, k = k, bs = "tp"),
                          data = data.frame(Y = Y, lon = tr$lon, lat = tr$lat),
                          weights = pmax(tr$n_svy, 1), method = "REML"),
                error = function(e) NULL)
  if (is.null(g)) return(rep(mean(tr$svy_prev), nrow(te)))
  p <- as.numeric(stats::predict(g, newdata = data.frame(lon = te$lon, lat = te$lat)))
  if (scale == "logit") p <- .ilogit(p)
  pmin(pmax(p, 0), 1)
}

f_screen30 <- function(X, y, v) screen_topK(X, y, v, 30L)
f_decorr20 <- function(X, y, v) decorr_reps(X, v, k = 20L)
f_curated  <- function(X, y, v) intersect(curated_vars(v), v)

MODELS <- list(
  list(name = "PROD enet30",            alpha = .5, cen = "none", feat = f_screen30),
  list(name = "ridge_all",              alpha = 0,  cen = "none", feat = NULL),
  list(name = "curated16_ridge",        alpha = 0,  cen = "none", feat = f_curated),
  list(name = "centered_own ridge_all", alpha = 0,  cen = "own",  feat = NULL),
  list(name = "centered_own decorr20",  alpha = .5, cen = "own",  feat = f_decorr20)
)

res <- list()
for (oc in OUTCOMES) {
  P <- assemble_outcome(oc, cov_list[FOUR]); if (is.null(P)) next
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  vars <- P$predictors
  message(sprintf("### %s (%d areas, %d predictors)", oc, nrow(d), length(vars)))
  for (ho in FOUR) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    for (sc in c("continuous", "logit")) {
      pr <- loco_tps(tr, te, sc)
      m <- loco_metrics(pr, te)
      if (!is.null(m)) {
        m$outcome <- oc; m$held_out <- ho; m$model <- "spatial_tps"; m$scale <- sc
        res[[paste(oc, ho, "tps", sc)]] <- m
      }
      for (mo in MODELS) {
        pr <- tryCatch(loco_scaled(tr, te, vars, mo$feat, mo$alpha, sc, mo$cen),
                       error = function(e) NULL)
        if (is.null(pr)) next
        m <- loco_metrics(pr, te); if (is.null(m)) next
        m$outcome <- oc; m$held_out <- ho; m$model <- mo$name; m$scale <- sc
        res[[paste(oc, ho, mo$name, sc)]] <- m
      }
    }
  }
}

out <- dplyr::bind_rows(res)
write.csv(out, "sandbox_parsimony/out/scale_sweep.csv", row.names = FALSE)

cat("\n=== LOCO: response scale, paired on the same 16 cells ===\n")
s <- out |> dplyr::group_by(model, scale) |>
  dplyr::summarise(cells = dplyr::n(),
                   rho = round(mean(spearman, na.rm = TRUE), 3),
                   r = round(mean(pearson, na.rm = TRUE), 3),
                   rmse = round(mean(rmse_pp, na.rm = TRUE), 1),
                   abs_bias = round(mean(abs(level_bias_pp), na.rm = TRUE), 1),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = scale, values_from = c(rho, r, rmse, abs_bias)) |>
  dplyr::mutate(rho_gain_from_continuous = round(rho_continuous - rho_logit, 3)) |>
  dplyr::arrange(dplyr::desc(rho_continuous))
print(as.data.frame(s), row.names = FALSE)

cat("\n=== how many districts get clamped on the logit scale? ===\n")
cl <- list()
for (oc in OUTCOMES) {
  P <- assemble_outcome(oc, cov_list[FOUR]); if (is.null(P)) next
  d <- P$data[is.finite(P$data$svy_prev), ]
  for (ctry in FOUR) {
    s2 <- d[d$country == ctry, ]; if (!nrow(s2)) next
    cl[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, n_areas = nrow(s2),
      pct_prev_zero = round(100 * mean(s2$svy_prev <= 0.005), 0),
      pct_prev_one  = round(100 * mean(s2$svy_prev >= 0.995), 0),
      stringsAsFactors = FALSE)
  }
}
print(as.data.frame(dplyr::bind_rows(cl)), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/scale_sweep.csv")
