# =============================================================================
# run_superlearner_example.R
# -----------------------------------------------------------------------------
# A self-contained, runnable example of a SuperLearner analysis on the
# simplified subset. It predicts an area's micronutrient-deficiency PREVALENCE
# from geospatial proxy predictors and reports cross-validated performance.
# It also shows a minimal MULTI-OUTCOME (sequential) example — predicting one
# deficiency using another (predicted) deficiency as an extra covariate.
#
# You do NOT need anything from the wider mn-prediction project to run this --
# just this folder and a few CRAN packages.
#
# HOW TO RUN (note the quotes — the folder name contains a space):
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/run_superlearner_example.R"
#
# REQUIRED PACKAGES (all on CRAN):
#   install.packages(c("SuperLearner", "glmnet", "ranger"))
# =============================================================================

set.seed(20240624)

# ── 0. Settings you may want to change ───────────────────────────────────────
LEVEL      <- "cluster"     # "cluster", "admin2", or "admin1"
OUTCOME    <- "women_iron"  # any outcome present in the data (see column names)
N_FOLDS    <- 5L
QUICK_TEST <- FALSE         # TRUE = tiny library, fewer folds (smoke test)

# ── 1. Locate and load the dataset ───────────────────────────────────────────
this_dir <- tryCatch(dirname(sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "."
data_dir <- file.path(this_dir, "data")
if (!dir.exists(data_dir)) data_dir <- file.path("simplified subset", "data")

file_map <- c(cluster = "mn_cluster.csv",
              admin2  = "mn_admin2.csv",
              admin1  = "mn_admin1.csv")
csv_path <- file.path(data_dir, file_map[[LEVEL]])
stopifnot(file.exists(csv_path))
dat <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Loaded %s: %d rows x %d columns\n", basename(csv_path),
            nrow(dat), ncol(dat)))

# ── 2. Identify predictor columns (shared across all outcomes) ───────────────
# The 16 proxy predictors are every column that is neither an identifier nor a
# per-outcome field (prev_* / n_* / ndef_* / var_*).
meta_cols      <- c("country", "admin1", "admin2", "cluster_id", "n_clusters")
outcome_cols   <- grep("^(prev_|n_|ndef_|var_)", names(dat), value = TRUE)
predictor_cols <- setdiff(names(dat), c(meta_cols, outcome_cols))
cat(sprintf("%d predictors: %s\n", length(predictor_cols),
            paste(predictor_cols, collapse = ", ")))

# Per-outcome columns follow a fixed naming convention:
prev_col <- function(oc) paste0("prev_", oc)   # survey-weighted prevalence (Y)
n_col    <- function(oc) paste0("n_",    oc)   # sample size  (use as weight)
var_col  <- function(oc) paste0("var_",  oc)   # design variance (precision wt)

# Simple median imputation for the few missing predictor cells.
X_all <- dat[, predictor_cols, drop = FALSE]
for (j in names(X_all))
  if (anyNA(X_all[[j]])) X_all[[j]][is.na(X_all[[j]])] <- median(X_all[[j]], na.rm = TRUE)

if (!requireNamespace("SuperLearner", quietly = TRUE))
  stop("Please install.packages('SuperLearner')")
library(SuperLearner)
have <- function(p) requireNamespace(p, quietly = TRUE)
sl_lib <- c("SL.mean", "SL.glm")
if (have("glmnet")) sl_lib <- c(sl_lib, "SL.glmnet")
if (have("ranger")) sl_lib <- c(sl_lib, "SL.ranger")
if (QUICK_TEST) { sl_lib <- c("SL.mean", "SL.glm"); N_FOLDS <- 3L }

# ── Helper: cross-validated SuperLearner on a given (Y, X, weights) ──────────
cv_report <- function(Y, X, w, label) {
  keep <- is.finite(Y)            # drop areas where this outcome was not measured
  Y <- Y[keep]; X <- X[keep, , drop = FALSE]; w <- w[keep]
  cv <- CV.SuperLearner(
    Y = Y, X = X, family = gaussian(), SL.library = sl_lib, obsWeights = w,
    V = N_FOLDS, innerCvControl = list(list(V = N_FOLDS)), verbose = FALSE)
  p <- cv$SL.predict
  rmse <- sqrt(mean((Y - p)^2)); r2 <- 1 - sum((Y - p)^2) / sum((Y - mean(Y))^2)
  cat(sprintf("\n== %s ==\n  areas=%d  CV RMSE=%.4f (baseline %.4f)  CV R2=%.3f\n",
              label, length(Y), rmse, sqrt(mean((Y - mean(Y))^2)), r2))
  w_avg <- colMeans(do.call(rbind, lapply(cv$AllSL, function(f) f$coef)))
  cat("  ensemble weights: ",
      paste(sprintf("%s=%.2f", names(w_avg), w_avg), collapse = "  "), "\n")
  invisible(p)
}

cat(sprintf("\nLearner library: %s | CV folds: %d\n",
            paste(sl_lib, collapse = ", "), N_FOLDS))

# ── 3. Single-outcome example: predict OUTCOME from the 16 proxies ───────────
cv_report(dat[[prev_col(OUTCOME)]], X_all, dat[[n_col(OUTCOME)]],
          sprintf("Single-outcome: %s prevalence (%s level)", OUTCOME, LEVEL))

# ── 4. Multi-outcome (sequential) example ────────────────────────────────────
# Borrow strength across deficiencies: first predict an "upstream" outcome
# (women_iron), then predict a "downstream" one (women_vitA) with the PREDICTED
# upstream prevalence added as an extra covariate. Compare to proxies-only.
# (This is the within-country form of the sequential idea; for transport you
#  would feed the upstream *prediction*, not the truth — see ANDY_KIM plan.)
up <- "women_iron"; down <- "women_vitA"
if (all(c(prev_col(up), prev_col(down)) %in% names(dat))) {
  ok <- is.finite(dat[[prev_col(up)]]) & is.finite(dat[[prev_col(down)]])
  Xd <- X_all[ok, , drop = FALSE]
  Yd <- dat[[prev_col(down)]][ok]; wd <- dat[[n_col(down)]][ok]
  cat(sprintf("\n--- sequential demo: %s -> %s (%d areas) ---\n", up, down, sum(ok)))
  cv_report(Yd, Xd, wd, sprintf("%s from proxies only", down))
  # Add the upstream prevalence as a covariate. In a real transport workflow this
  # column would be the cross-validated PREDICTION of the upstream outcome, not
  # its observed value; we use the observed value here purely to illustrate the
  # structure and the potential covariance gain.
  Xd_seq <- cbind(Xd, upstream_iron = dat[[prev_col(up)]][ok])
  cv_report(Yd, Xd_seq, wd, sprintf("%s from proxies + %s", down, up))
}

cat("\nDone. Edit LEVEL / OUTCOME at the top, or loop cv_report() over every\n")
cat("outcome. See data_dictionary.csv for the full column reference.\n")
