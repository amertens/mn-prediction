# =============================================================================
# sandbox_parsimony/R/10_hygiene_effect.R
#
# Does turning GEE_COVARIATE_HYGIENE on actually help?
#
# R/gee_band_semantics.R already documents the problem precisely: ~46% of the
# gee_ Admin-2 columns are cross-band `_annual_*` summaries, many of them
# averaging non-commensurable bands (TerraClimate aet + pdsi + surface pressure
# in one number), and ~243 columns on Sierra Leone are exact duplicates of
# another column. The pruning is written and gated OFF by default.
#
# This scores the flag on the same repeated-CV protocol as 04_within_country.R
# and on LOCO, so the decision to flip the default has evidence behind it.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
source("R/gee_band_semantics.R")

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES  <- c("child_vitA", "child_iron", "women_iron", "women_folate")
MIN_AREAS <- 25L

# Names in the sandbox pool have had year tokens stripped (harmonize_gee_names),
# so the static-year-duplicate rule has nothing to match; the non-commensurable
# cross-band summary rule is the part that bites, and it is name-based, so it
# applies unchanged.
hygiene_drop <- function(vars) {
  Sys.setenv(GEE_COVARIATE_HYGIENE = "true")
  d <- prune_gee_covariates(vars, verbose = FALSE)
  Sys.unsetenv("GEE_COVARIATE_HYGIENE")
  d
}

SPECS <- list(
  list(name = "enet_screen30",
       features = function(X, y, v) screen_topK(X, y, v, 30L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "ridge_all", features = NULL,
       fit = function(...) fit_glmnet_logit(..., alpha = 0)),
  list(name = "decorr20_enet",
       features = function(X, y, v) decorr_reps(X, v, k = 20L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5))
)

res <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  drop <- hygiene_drop(P$predictors)
  vars_on  <- setdiff(P$predictors, drop)
  vars_off <- P$predictors
  message(sprintf("\n### %s: %d predictors -> %d after hygiene (dropped %d)",
                  oc, length(vars_off), length(vars_on), length(drop)))

  # ---- within-country -----------------------------------------------------
  for (ctry in P$countries) {
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < MIN_AREAS) next
    for (sp in SPECS) for (hy in c("off", "on")) {
      v <- if (hy == "on") vars_on else vars_off
      r <- run_cv(dat, v, sp, n_rep = 5L, k = 5L)
      if (is.null(r)) next
      s <- summarise_cv(r)
      s$outcome <- oc; s$country <- ctry; s$model <- sp$name
      s$hygiene <- hy; s$eval <- "within"
      res[[paste("w", oc, ctry, sp$name, hy)]] <- s
    }
    message(sprintf("   [within] %s done", ctry))
  }

  # ---- LOCO ---------------------------------------------------------------
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  for (ho in unique(d$country)) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    for (spn in c("enet_screen30", "ridge_all", "decorr20_enet")) {
      feat <- switch(spn,
        enet_screen30 = function(X, y, v) screen_topK(X, y, v, 30L),
        ridge_all     = NULL,
        decorr20_enet = function(X, y, v) decorr_reps(X, v, k = 20L))
      alpha <- if (spn == "ridge_all") 0 else 0.5
      for (hy in c("off", "on")) {
        v <- if (hy == "on") vars_on else vars_off
        pr <- tryCatch(loco_one(tr, te, v, feat, alpha, "own", "train_mean"),
                       error = function(e) NULL)
        if (is.null(pr)) next
        m <- loco_metrics(pr, te); if (is.null(m)) next
        m$outcome <- oc; m$country <- ho; m$model <- spn
        m$hygiene <- hy; m$eval <- "loco"
        res[[paste("l", oc, ho, spn, hy)]] <- m
      }
    }
  }
  message("   [loco] done")
}

out <- dplyr::bind_rows(res)
write.csv(out, "sandbox_parsimony/out/hygiene_effect.csv", row.names = FALSE)

cat("\n=== GEE covariate hygiene: mean Spearman, paired on the same cells ===\n")
s <- out |> dplyr::group_by(eval, model, hygiene) |>
  dplyr::summarise(cells = dplyr::n(),
                   rho = round(mean(spearman, na.rm = TRUE), 3),
                   r = round(mean(pearson, na.rm = TRUE), 3),
                   rmse = round(mean(rmse_pp, na.rm = TRUE), 2),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = hygiene, values_from = c(rho, r, rmse, cells))
print(as.data.frame(s), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/hygiene_effect.csv")
