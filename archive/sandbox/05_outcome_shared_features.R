# =============================================================================
# sandbox/05_outcome_shared_features.R
#
# Hypothesis H4: features that are jointly predictive of ALL four LOCO
# outcomes (iron child, iron women, vitA child, vitA women) are more likely
# to transport than features that only predict one. Compute per-feature
# average absolute correlation with all 4 outcomes, keep the top-K
# multi-outcome predictors, fit a small penalised model per outcome on
# those shared features.
#
# Rationale: micronutrient deficiencies share underlying risk pathways
# (diet quality, parasite burden, sanitation, livelihood) — features
# that index those shared upstream causes are by construction more
# likely to be biologically meaningful and to transport. Features that
# correlate with only one outcome are likely noise or one-off coincidences.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

# For each variable, compute mean absolute correlation across all 4 outcomes
# on training data (which has all 4 outcomes in their respective pooled
# datasets keyed by country × Admin2).
shared_feature_scores <- function(pooled_all, training_countries, vars) {
  cor_by_outcome <- vapply(names(pooled_all), function(oc) {
    pd <- pooled_all[[oc]]$pooled_data
    tr <- pd[pd$country %in% training_countries, , drop = FALSE]
    vapply(vars, function(v)
      suppressWarnings(abs(cor(tr[[v]], tr$svy_prev, use = "pairwise"))),
      numeric(1))
  }, numeric(length(vars)))
  # cor_by_outcome is variables × outcomes
  data.frame(
    variable           = vars,
    mean_abs_r         = rowMeans(cor_by_outcome, na.rm = TRUE),
    min_abs_r          = apply(cor_by_outcome, 1, min, na.rm = TRUE),
    n_outcomes_above_0_1 = rowSums(cor_by_outcome > 0.1, na.rm = TRUE),
    stringsAsFactors   = FALSE
  )
}

# Predict-fn closure: for the held-out country in this outcome, build the
# shared-feature ranking on the TRAINING countries across ALL 4 outcomes,
# select top-K, fit a penalised regression on the current outcome.
make_shared_predict <- function(pooled_all, this_outcome, top_k = 10) {
  function(train, test, vars) {
    vars <- screen_vars(train, vars)
    if (length(vars) < 3) return(NULL)
    held_out <- unique(test$country)
    training_countries <- setdiff(unique(train$country), held_out)
    scores <- shared_feature_scores(pooled_all, training_countries, vars)
    # Prefer features predictive of at least 3 of 4 outcomes; else fall back
    # to top mean_abs_r.
    keep <- scores[scores$n_outcomes_above_0_1 >= 3, , drop = FALSE]
    if (nrow(keep) < 2) keep <- scores[order(-scores$mean_abs_r), , drop = FALSE]
    selected <- head(keep[order(-keep$mean_abs_r), ]$variable, top_k)
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
  list(top_k = 6,  name = "shared_top6"),
  list(top_k = 10, name = "shared_top10"),
  list(top_k = 15, name = "shared_top15")
)
for (cfg in configs) {
  cat(sprintf("\n[H4] %s ...\n", cfg$name))
  per_outcome <- list()
  for (oc in OUTCOMES) {
    fn <- make_shared_predict(pooled_all, this_outcome = oc, top_k = cfg$top_k)
    r <- loco_evaluate(pooled_all[[oc]], fn, outcome = oc,
                        method_name = cfg$name)
    per_outcome[[oc]] <- r
  }
  res[[cfg$name]] <- dplyr::bind_rows(per_outcome)
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "05_shared_features.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "05_shared_features.csv")))
