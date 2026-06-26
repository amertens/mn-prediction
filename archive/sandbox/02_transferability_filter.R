# =============================================================================
# sandbox/02_transferability_filter.R
#
# Hypothesis H1: features whose univariate correlation with the outcome is
# STABLE across training countries transport better than features that flip
# sign or magnitude. Pre-screen variables by cross-country correlation
# stability, then fit a small penalised regression.
#
# Rationale: the SuperLearner's LOCO failures (e.g., child_iron Sierra Leone
# at r = -0.19) come from variables whose effect direction differs between
# training and held-out country. If we drop those *before* fitting, the
# remaining model is by construction transportable.
#
# Algorithm:
#   For each candidate variable v:
#     For each training country c:
#       r_c = cor(v, svy_prev) on country c's data
#     stability(v) = 1 - sd(r_c) / (|mean(r_c)| + eps)
#     OR signed-agreement(v) = fraction of training countries where sign(r_c)
#                              agrees with sign(mean(r_c))
#   Keep top-K by stability AND signed-agreement >= threshold.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

# Per-variable transferability scores on a training set.
transferability_scores <- function(train, vars, eps = 0.05) {
  countries <- unique(train$country)
  cor_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(countries),
                     dimnames = list(vars, countries))
  for (c in countries) {
    sub <- train[train$country == c, , drop = FALSE]
    if (nrow(sub) < 5) next
    for (v in vars) {
      cor_mat[v, c] <- suppressWarnings(cor(sub[[v]], sub$svy_prev,
                                              use = "pairwise"))
    }
  }
  mean_r  <- rowMeans(cor_mat, na.rm = TRUE)
  sd_r    <- apply(cor_mat, 1, stats::sd,  na.rm = TRUE)
  sign_agree <- rowMeans(sign(cor_mat) == sign(mean_r), na.rm = TRUE)
  stab    <- abs(mean_r) / (sd_r + eps)
  data.frame(variable = vars,
              mean_r = mean_r, sd_r = sd_r, abs_mean_r = abs(mean_r),
              stability = stab, signed_agreement = sign_agree,
              stringsAsFactors = FALSE)
}

# Predict-fn: filter then penalised regression.
fit_predict_transferable <- function(train, test, vars,
                                       top_k = 8, min_agreement = 0.75) {
  vars <- screen_vars(train, vars)
  if (length(vars) == 0) return(NULL)
  scores <- transferability_scores(train, vars)
  # Keep variables with high signed agreement AND high stability score.
  keep <- scores[scores$signed_agreement >= min_agreement, , drop = FALSE]
  if (nrow(keep) < 2) {
    # Fall back: just use top-K stable variables regardless of sign agreement.
    keep <- scores[order(-scores$stability), , drop = FALSE]
  }
  keep <- keep[order(-keep$stability), , drop = FALSE]
  selected <- head(keep$variable, top_k)
  if (length(selected) < 2) return(NULL)
  imp <- impute_median(train, test, selected)
  Y   <- train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  fit <- tryCatch(
    glmnet::cv.glmnet(Xtr, Y, alpha = 0.5, nfolds = min(nrow(Xtr), 10L),
                       weights = wt),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  p <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.min"))
  pmin(pmax(p, 0), 1)
}

pooled_all <- load_all_pooled()

res <- list()
# Try a couple of (top_k, agreement) settings to see sensitivity.
configs <- list(
  list(top_k = 8,  min_agreement = 0.75, name = "transferable_k8_a75"),
  list(top_k = 5,  min_agreement = 1.00, name = "transferable_k5_a100"),
  list(top_k = 12, min_agreement = 0.50, name = "transferable_k12_a50")
)
for (cfg in configs) {
  cat(sprintf("\n[H1] %s ...\n", cfg$name))
  r <- loco_evaluate_all(pooled_all,
    predict_fn  = function(train, test, vars)
      fit_predict_transferable(train, test, vars,
                                top_k = cfg$top_k,
                                min_agreement = cfg$min_agreement),
    method_name = cfg$name)
  res[[cfg$name]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "02_transferability.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "02_transferability.csv")))
