# =============================================================================
# explore/scripts/08_np_estimators.R   [probe NP-01]
#
# QUESTION. The project has tried elastic net, the zero-tuning domain index,
# HAL/PCHAL and a SuperLearner. The standard chemometrics and genomics answers
# to n << p are untried here. Do any of them beat the index?
#
# ARMS
#   pls2         partial least squares, 2 components, no tuning. PLS was built
#                for n << p with collinear predictors (spectra); it finds the
#                directions of X that covary with y, unlike PCA which ignores y
#   pcr2         principal components regression, 2 components - the
#                unsupervised twin of pls2, to separate "low rank helps" from
#                "supervision helps"
#   spca         supervised PCA (Bair-Tibshirani): screen columns by marginal
#                correlation ON TRAINING ROWS, then PC1 of the survivors
#   mcp          MCP penalty (ncvreg) - less biased than lasso at large signals
#   stabsel      stability selection: elastic net over 50 subsamples, keep
#                columns chosen in >= 60%, refit ridge on those
#   ridge_cv     plain CV-tuned ridge on all predictors. Pairs with probe 01's
#                REML-ridge: same estimator, shrinkage chosen by CV instead of
#                estimated, which isolates how much the TUNING costs
#   index_std    the domain index with each axis standardised before summing.
#                XO-01 found the index sums UN-standardised axes, so any axis
#                with small spread is silently under-weighted; this is the fix
#   + domain_index, spatial, null on identical folds
#
#   Rscript explore/scripts/08_np_estimators.R
# -> explore/out/08_np_cells.csv, 08_np_loco.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
suppressPackageStartupMessages({library(glmnet)})

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()

has <- function(p) requireNamespace(p, quietly = TRUE)

# ── arms ────────────────────────────────────────────────────────────────────

arm_pls <- function(ncomp = 2L) function(tr, te, y, X, D, aux) {
  if (!has("pls")) return(rep(mean(y[tr]), length(te)))
  d <- data.frame(y = y[tr]); d$X <- X[tr, , drop = FALSE]
  nc <- min(ncomp, length(tr) - 2L, ncol(X))
  if (nc < 1) return(rep(mean(y[tr]), length(te)))
  fit <- try(pls::plsr(y ~ X, ncomp = nc, data = d, scale = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(mean(y[tr]), length(te)))
  nd <- list(X = X[te, , drop = FALSE])
  as.numeric(stats::predict(fit, newdata = nd, ncomp = nc))
}

arm_pcr <- function(ncomp = 2L) function(tr, te, y, X, D, aux) {
  nc <- min(ncomp, length(tr) - 2L, ncol(X))
  if (nc < 1) return(rep(mean(y[tr]), length(te)))
  pc <- try(stats::prcomp(X[tr, , drop = FALSE], center = TRUE, scale. = FALSE),
            silent = TRUE)
  if (inherits(pc, "try-error")) return(rep(mean(y[tr]), length(te)))
  Ztr <- pc$x[, seq_len(nc), drop = FALSE]
  Zte <- scale(X[te, , drop = FALSE], center = pc$center, scale = FALSE) %*%
    pc$rotation[, seq_len(nc), drop = FALSE]
  cf <- stats::lm.fit(cbind(1, Ztr), y[tr])$coefficients
  cf[!is.finite(cf)] <- 0
  as.numeric(cbind(1, Zte) %*% cf)
}

#' Supervised PCA: screen on training rows only, then PC1 of the survivors
arm_spca <- function(topk = 20L) function(tr, te, y, X, D, aux) {
  r <- suppressWarnings(apply(X[tr, , drop = FALSE], 2, function(z)
    if (stats::sd(z) == 0) 0 else stats::cor(z, y[tr], method = "spearman")))
  r[!is.finite(r)] <- 0
  keep <- order(abs(r), decreasing = TRUE)[seq_len(min(topk, ncol(X)))]
  Xs <- X[, keep, drop = FALSE]
  # orient each survivor by the sign of its training correlation, so members
  # of the screened set cannot cancel (the project's own domain-score logic)
  Xs <- sweep(Xs, 2, sign(r[keep] + 1e-12), "*")
  pc <- try(stats::prcomp(Xs[tr, , drop = FALSE], center = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) return(rep(mean(y[tr]), length(te)))
  ztr <- pc$x[, 1]
  zte <- as.numeric(scale(Xs[te, , drop = FALSE], center = pc$center, scale = FALSE) %*%
                      pc$rotation[, 1])
  cf <- stats::lm.fit(cbind(1, ztr), y[tr])$coefficients
  cf[!is.finite(cf)] <- 0
  as.numeric(cbind(1, zte) %*% cf)
}

# cv.ncvreg on 383 columns is by far the most expensive arm here (it was on
# pace for ~6 hours over the full grid). nlambda is cut to 25 and the inner CV
# to 3 folds: this is a screen for whether a non-convex penalty is even in the
# running, not a tuned production fit.
arm_mcp <- function(tr, te, y, X, D, aux) {
  if (!has("ncvreg")) return(rep(mean(y[tr]), length(te)))
  fit <- try(ncvreg::cv.ncvreg(X[tr, , drop = FALSE], y[tr], penalty = "MCP",
                               nfolds = 3, nlambda = 25), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(mean(y[tr]), length(te)))
  as.numeric(stats::predict(fit, X = X[te, , drop = FALSE], lambda = fit$lambda.min))
}

#' Stability selection: keep what is chosen often across subsamples
arm_stabsel <- function(B = 25L, thresh = 0.6, alpha = 0.5) function(tr, te, y, X, D, aux) {
  n <- length(tr)
  cnt <- rep(0L, ncol(X)); names(cnt) <- colnames(X)
  for (b in seq_len(B)) {
    s <- sample(tr, max(5L, floor(n / 2)))
    if (length(unique(y[s])) < 3) next
    f <- try(glmnet::glmnet(X[s, , drop = FALSE], y[s], alpha = alpha,
                            nlambda = 25), silent = TRUE)
    if (inherits(f, "try-error")) next
    # the solution roughly half-way down the path: a fixed sparsity, not a tuned one
    j <- max(1L, floor(ncol(f$beta) / 2))
    nz <- which(as.numeric(f$beta[, j]) != 0)
    cnt[nz] <- cnt[nz] + 1L
  }
  keep <- which(cnt / B >= thresh)
  if (length(keep) < 2) keep <- order(cnt, decreasing = TRUE)[1:min(5, ncol(X))]
  .v2_enet(X[tr, keep, drop = FALSE], y[tr], X[te, keep, drop = FALSE], alpha = 0)
}

arm_ridge_cv <- function(tr, te, y, X, D, aux) {
  f <- try(glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], alpha = 0, nfolds = 5),
           silent = TRUE)
  if (inherits(f, "try-error")) return(rep(mean(y[tr]), length(te)))
  as.numeric(stats::predict(f, newx = X[te, , drop = FALSE], s = "lambda.min"))
}

#' The domain index with each axis standardised (on TRAINING rows) before the
#' weighted sum. XO-01: the index sums un-standardised axes, so an axis with
#' small spread contributes little however strong its association.
arm_index_std <- function(tr, te, y, X, D, aux) {
  mu <- colMeans(D[tr, , drop = FALSE])
  sdv <- apply(D[tr, , drop = FALSE], 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Ds <- sweep(sweep(D, 2, mu, "-"), 2, sdv, "/")
  arm_domain_index_v2(tr, te, y, X, Ds, aux)
}

ARMS <- c(exp_baseline_arms(),
          list(pls2 = arm_pls(2L), pcr2 = arm_pcr(2L), spca = arm_spca(20L),
               mcp = arm_mcp, stabsel = arm_stabsel(), ridge_cv = arm_ridge_cv,
               index_std = arm_index_std))

for (p in c("pls", "ncvreg"))
  if (!has(p)) message("NOTE: package '", p, "' missing - its arm returns the train mean")

rows <- list()
ix <- exp_cell_index(E)
for (i in seq_len(nrow(ix))) {
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tgt),
                     error = function(e) NULL)
    if (is.null(cell)) next
    rows[[paste(i, tgt, "A")]] <- exp_infill(cell, ARMS, reps = REPS)
    rows[[paste(i, tgt, "B")]] <- exp_region(cell, ARMS)
  }
  message("  ", ix$country[i], " ", ix$outcome[i])
}

loco <- list()
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tgt, outcomes = on)
    if (length(cl) < 3) next
    loco[[paste(tgt, on)]] <- exp_loco(cl, ARMS[setdiff(names(ARMS), "spatial")],
                                       domain_of = E$domain_of)
    message("  LOCO ", tgt, " ", on)
  }
}

raw <- dplyr::bind_rows(rows); sm <- exp_summarise(raw)
exp_write(sm, "08_np_cells")
lc <- dplyr::bind_rows(loco); exp_write(lc, "08_np_loco")

cat("\n== in-country: median Spearman over cells ==\n")
a <- aggregate(spearman ~ estimand + target + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$estimand, a$target, -a$spearman), ], row.names = FALSE)

cat("\n== transport (LOCO) ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
b$positive <- aggregate(spearman ~ target + arm, data = lc,
                        FUN = function(z) sum(z > 0, na.rm = TRUE))$spearman
b$cells <- aggregate(spearman ~ target + arm, data = lc,
                     FUN = function(z) sum(is.finite(z)))$spearman
print(b[order(b$target, -b$spearman), ], row.names = FALSE)

cat("\n== index vs index_std, paired by cell (in-fill, level) ==\n")
p <- sm[sm$estimand == "infill" & sm$target == "level" &
          sm$arm %in% c("domain_index", "index_std"), ]
w <- reshape(p[, c("country", "outcome", "arm", "spearman")],
             idvar = c("country", "outcome"), timevar = "arm", direction = "wide")
names(w) <- sub("^spearman\\.", "", names(w))
w$gain <- round(w$index_std - w$domain_index, 3)
print(w[order(-w$gain), ], row.names = FALSE)
cat("mean gain from standardising the axes:", round(mean(w$gain, na.rm = TRUE), 4),
    " cells improved:", sum(w$gain > 0, na.rm = TRUE), "of", sum(is.finite(w$gain)), "\n")
