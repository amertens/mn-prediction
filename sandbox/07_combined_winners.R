# =============================================================================
# sandbox/07_combined_winners.R
#
# Hypothesis H6: combine H5 (within-country variance filter) + H4 (outcome-
# shared features). Pre-screen variables on:
#   (i)  high within-country variance ratio (don't encode country fixed effects)
#   (ii) high mean absolute correlation across all 4 outcomes (index shared
#        upstream causal factors)
# Then fit a small penalised regression.
#
# Rationale: both filters target transportability from different angles —
# H5 targets feature INVARIANCE (variation that exists within every country)
# and H4 targets feature SHARED MEANING (variables that predict multiple
# outcomes). Their intersection should be the most-transportable subset of
# all 152 GEE candidates.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

country_variance_decomp <- function(train, vars) {
  out <- vapply(vars, function(v) {
    x <- train[[v]]; ctry <- train$country
    grand_var <- stats::var(x, na.rm = TRUE)
    if (!is.finite(grand_var) || grand_var == 0) return(0)
    ctry_means <- tapply(x, ctry, mean, na.rm = TRUE)
    x_centered <- x - ctry_means[ctry]
    within_var <- stats::var(x_centered, na.rm = TRUE)
    within_var / grand_var
  }, numeric(1))
  data.frame(variable = vars, within_ratio = out)
}

shared_feature_scores <- function(pooled_all, training_countries, vars) {
  cor_by_outcome <- vapply(names(pooled_all), function(oc) {
    pd <- pooled_all[[oc]]$pooled_data
    tr <- pd[pd$country %in% training_countries, , drop = FALSE]
    vapply(vars, function(v)
      suppressWarnings(abs(cor(tr[[v]], tr$svy_prev, use = "pairwise"))),
      numeric(1))
  }, numeric(length(vars)))
  data.frame(variable = vars,
              mean_abs_r = rowMeans(cor_by_outcome, na.rm = TRUE),
              min_abs_r  = apply(cor_by_outcome, 1, min, na.rm = TRUE))
}

make_combined_predict <- function(pooled_all,
                                    min_within_ratio = 0.5,
                                    top_k = 8) {
  function(train, test, vars) {
    vars <- screen_vars(train, vars)
    if (length(vars) < 3) return(NULL)
    held_out <- unique(test$country)
    training_countries <- setdiff(unique(train$country), held_out)
    # Filter 1: within-country variance ratio
    inv_scores <- country_variance_decomp(train, vars)
    pass_invariance <- inv_scores[inv_scores$within_ratio >= min_within_ratio,
                                    , drop = FALSE]
    if (nrow(pass_invariance) < 2) pass_invariance <- inv_scores
    # Filter 2: outcome-shared score (using only vars that passed filter 1)
    shr_scores <- shared_feature_scores(pooled_all, training_countries,
                                          pass_invariance$variable)
    combined <- merge(pass_invariance, shr_scores, by = "variable")
    # Composite score: product of within-ratio and mean-abs-r (both want high).
    combined$composite <- combined$within_ratio * combined$mean_abs_r
    selected <- head(combined[order(-combined$composite), ]$variable, top_k)
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
}

pooled_all <- load_all_pooled()

res <- list()
configs <- list(
  list(min_within_ratio = 0.50, top_k = 8,  name = "combined_w50_k8"),
  list(min_within_ratio = 0.70, top_k = 8,  name = "combined_w70_k8"),
  list(min_within_ratio = 0.70, top_k = 12, name = "combined_w70_k12"),
  list(min_within_ratio = 0.50, top_k = 15, name = "combined_w50_k15")
)
for (cfg in configs) {
  cat(sprintf("\n[H6] %s ...\n", cfg$name))
  per_outcome <- list()
  for (oc in OUTCOMES) {
    fn <- make_combined_predict(pooled_all,
                                  min_within_ratio = cfg$min_within_ratio,
                                  top_k = cfg$top_k)
    r <- loco_evaluate(pooled_all[[oc]], fn, outcome = oc,
                        method_name = cfg$name)
    per_outcome[[oc]] <- r
  }
  res[[cfg$name]] <- dplyr::bind_rows(per_outcome)
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "07_combined.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "07_combined.csv")))
