# =============================================================================
# _helpers.R  — shared utilities for the methods scripts in this folder
# -----------------------------------------------------------------------------
# Self-contained: depends only on the simplified-subset CSVs and `glmnet`.
# Source this from each method script. Run the scripts from the repo root, e.g.
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/national_anchor_calibration.R"
# =============================================================================

suppressWarnings(suppressMessages({ library(glmnet) }))

# Load one aggregation level (cluster / admin2 / admin1). Tries a few paths so it
# works whether you run from the repo root or from inside the folder.
load_level <- function(level = "admin2") {
  cands <- c(
    file.path("simplified subset", "data", paste0("mn_", level, ".csv")),
    file.path("..", "data", paste0("mn_", level, ".csv")),
    file.path("data", paste0("mn_", level, ".csv"))
  )
  p <- cands[file.exists(cands)][1]
  if (is.na(p)) stop("Could not find mn_", level, ".csv from ", getwd())
  read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
}

# The 16 proxy predictors = every column that is not an id or a per-outcome field.
predictor_cols <- function(df) {
  meta <- c("country", "admin1", "admin2", "cluster_id", "n_clusters")
  oc   <- grep("^(prev_|n_|ndef_|var_)", names(df), value = TRUE)
  setdiff(names(df), c(meta, oc))
}

# Modeling frame for one outcome: area rows with a finite prevalence + sample.
# Predictors are median-imputed (the only preprocessing needed here).
prep_outcome <- function(df, outcome) {
  pcols <- predictor_cols(df)
  y    <- df[[paste0("prev_",  outcome)]]
  n    <- df[[paste0("n_",     outcome)]]
  ndef <- df[[paste0("ndef_",  outcome)]]
  keep <- is.finite(y) & is.finite(n) & n > 0
  X <- as.matrix(df[keep, pcols, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    v <- X[, j]; m <- suppressWarnings(median(v[is.finite(v)], na.rm = TRUE))
    if (!is.finite(m)) m <- 0
    v[!is.finite(v)] <- m; X[, j] <- v
  }
  list(X = X, y = y[keep], n = n[keep], ndef = ndef[keep],
       country = df$country[keep], area = df$admin2[keep], predictors = pcols)
}

# ── Base area-level learner: weighted elastic-net (regularised for p >> n) ────
# family = "gaussian" on prevalence (weighted by n), or "binomial" on the
# (deficient, not-deficient) counts. Returns a predict() closure giving
# prevalence probabilities.
fit_area <- function(X, y, n, ndef, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  nf <- max(3, min(5, floor(nrow(X) / 3)))
  if (family == "gaussian") {
    cv <- cv.glmnet(X, y, weights = n, family = "gaussian", alpha = 0.5, nfolds = nf)
    function(newX) as.numeric(predict(cv, newX, s = "lambda.min", type = "response"))
  } else {
    # glmnet 2-column binomial response is (failures, successes); the SECOND
    # column is the modelled "success" class, so deficiency counts go second.
    yb <- cbind(pmax(n - ndef, 0), ndef)
    cv <- cv.glmnet(X, yb, family = "binomial", alpha = 0.5, nfolds = nf)
    function(newX) as.numeric(predict(cv, newX, s = "lambda.min", type = "response"))
  }
}

# ── Metrics ──────────────────────────────────────────────────────────────────
rmse_pp <- function(a, b) round(100 * sqrt(mean((a - b)^2)), 2)
mae_pp  <- function(a, b) round(100 * mean(abs(a - b)), 2)
bias_pp <- function(a, b) round(100 * mean(a - b), 2)            # pred - truth
spear   <- function(a, b) round(suppressWarnings(stats::cor(a, b, method = "spearman")), 3)
clip01  <- function(p, e = 1e-4) pmin(pmax(p, e), 1 - e)
