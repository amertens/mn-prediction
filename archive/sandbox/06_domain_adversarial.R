# =============================================================================
# sandbox/06_domain_adversarial.R
#
# Hypothesis H5: train predictions that are predictive of the OUTCOME but
# NOT predictive of the country. Domain-invariant features should transport
# better than country-specific features by construction.
#
# Approximation (no full DANN setup): use a two-step weighting scheme.
#   step 1: for each variable, compute residual variance after regressing it
#           on the country dummy. Variables with high *between-country*
#           variance contribute little to within-country admin-2 ranking
#           AND will pick up country fixed effects in LOCO.
#   step 2: keep variables whose within-country variance dominates
#           (variance ratio within / total > 0.5), then fit a penalised
#           regression on those.
#
# Rationale: this is the variance-decomposition argument from §5.5 of
# docs/data_and_approach_ideas.qmd. A variable whose effect is mostly a
# country-level shift won't help predict admin-2 ranking in a held-out
# country; a variable whose effect mostly varies WITHIN country captures
# transferable within-country structure.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

# Compute within/between variance ratios for each variable on the training set.
country_variance_decomp <- function(train, vars) {
  out <- vapply(vars, function(v) {
    x <- train[[v]]; ctry <- train$country
    grand_var <- stats::var(x, na.rm = TRUE)
    if (!is.finite(grand_var) || grand_var == 0) return(0)
    # Within-country variance: pooled variance after subtracting country mean.
    ctry_means <- tapply(x, ctry, mean, na.rm = TRUE)
    x_centered <- x - ctry_means[ctry]
    within_var <- stats::var(x_centered, na.rm = TRUE)
    within_var / grand_var
  }, numeric(1))
  data.frame(variable = vars, within_ratio = out,
              stringsAsFactors = FALSE)
}

fit_predict_domain_invariant <- function(train, test, vars,
                                            min_within_ratio = 0.5,
                                            top_k = 8) {
  vars <- screen_vars(train, vars)
  if (length(vars) == 0) return(NULL)
  scores <- country_variance_decomp(train, vars)
  # Keep variables with high within-country variance ratio; fallback to top by
  # within_ratio if filter is too aggressive.
  keep <- scores[scores$within_ratio >= min_within_ratio, , drop = FALSE]
  if (nrow(keep) < 2) keep <- scores
  # Among the kept, rank by univariate correlation with outcome on training.
  imp_pre <- impute_median(train, test, keep$variable)
  keep$abs_r <- vapply(keep$variable, function(v)
    suppressWarnings(abs(cor(imp_pre$train[[v]], train$svy_prev,
                                use = "pairwise"))), numeric(1))
  selected <- head(keep[order(-keep$abs_r), ]$variable, top_k)
  if (length(selected) < 2) return(NULL)
  imp <- impute_median(train, test, selected)
  Y   <- train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  fit <- tryCatch(glmnet::cv.glmnet(Xtr, Y, alpha = 0.5,
                                        nfolds = min(nrow(Xtr), 10L),
                                        weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  p <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.min"))
  pmin(pmax(p, 0), 1)
}

pooled_all <- load_all_pooled()

res <- list()
configs <- list(
  list(min_within_ratio = 0.30, top_k = 8,  name = "domain_inv_w30_k8"),
  list(min_within_ratio = 0.50, top_k = 8,  name = "domain_inv_w50_k8"),
  list(min_within_ratio = 0.70, top_k = 8,  name = "domain_inv_w70_k8"),
  list(min_within_ratio = 0.50, top_k = 12, name = "domain_inv_w50_k12")
)
for (cfg in configs) {
  cat(sprintf("\n[H5] %s ...\n", cfg$name))
  r <- loco_evaluate_all(pooled_all,
    predict_fn  = function(train, test, vars)
      fit_predict_domain_invariant(train, test, vars,
                                      min_within_ratio = cfg$min_within_ratio,
                                      top_k = cfg$top_k),
    method_name = cfg$name)
  res[[cfg$name]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "06_domain_adversarial.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "06_domain_adversarial.csv")))
