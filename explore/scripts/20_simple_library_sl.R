# =============================================================================
# explore/scripts/20_simple_library_sl.R   [probe SL-07]
#
# A SUPERLEARNER OF MANY *SIMPLE* MODELS
#
# The project has tested a SuperLearner over its protocol ARMS (SL-06: every
# candidate inside one cross-validated SuperLearner) and over tuned learners
# (SL-01/02/03), and the zero-tuning index beat all of them. What it has NOT
# tested is the other library design: very many very simple candidates - one
# variable each, pairs, single principal components - which is a different
# bias/variance trade from a handful of flexible learners.
#
# THE THEORY, AND WHY IT IS NOT OBVIOUSLY GOOD HERE. The SuperLearner oracle
# inequality (van der Laan & Dudoit 2003; van der Vaart, Dudoit & van der Laan
# 2006) bounds the CV-selector's risk by the oracle's risk plus a term of order
# (1 + log K)/n, and permits the library K to grow polynomially in n. At
# n = 14-87 that penalty is not small: with K = 120 and n = 30,
# log(120)/30 = 0.16, before the constant. The asymptotic promise that "adding
# candidates is nearly free" is exactly the promise that fails at this n.
#
# THE COUNTER-ARGUMENT, which is why it is worth running anyway. Simple
# candidates have far lower variance than tuned ones, so the CV selection is
# more stable even though the oracle risk is higher. That is the
# componentwise-boosting intuition.
#
# THE PREDICTION, stated before running. A non-negative convex combination of
# univariate least-squares fits is algebraically close to a shrunken ridge on
# the same variables - and WS-01 already found the index IS a max-shrinkage
# ridge that nothing beats in-country. So this should land near the index
# rather than above it. Recorded so the result can contradict it.
#
# ARMS
#   sl_uni      NNLS-weighted SuperLearner over one-variable OLS on every
#               domain axis, plus an intercept-only learner
#   sl_uni_pair + all pairs among the top-10 screened axes
#   sl_pc       + single principal components of the full predictor matrix
#   sl_rank     the same library, weights fitted to rank loss rather than MSE
#               (SL-02 found a rank meta-learner recovers the index; this asks
#               whether it does so on a simple library too)
#   + domain_index and null_train_mean on identical folds
#
#   Rscript explore/scripts/20_simple_library_sl.R
# -> explore/out/20_simple_library_sl.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
suppressPackageStartupMessages({library(nnls)})

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
ix <- exp_cell_index(E)
HAVE_NNLS <- requireNamespace("nnls", quietly = TRUE)
if (!HAVE_NNLS) message("NOTE: nnls missing; falling back to a projected least-squares solve")

#' non-negative least squares meta-weights, normalised to sum to 1
nnls_w <- function(Z, y) {
  w <- if (HAVE_NNLS) {
    tryCatch(nnls::nnls(Z, y)$x, error = function(e) NULL)
  } else {
    b <- tryCatch(stats::lm.fit(Z, y)$coefficients, error = function(e) NULL)
    if (!is.null(b)) pmax(b, 0) else NULL
  }
  if (is.null(w) || !any(is.finite(w)) || sum(w, na.rm = TRUE) <= 0)
    return(rep(1 / ncol(Z), ncol(Z)))
  w[!is.finite(w)] <- 0
  w / sum(w)
}

#' Build the library: each entry is a list(fit = function(tr), predict = ...)
#' represented compactly as a column-index set; every learner is OLS on those
#' columns, which keeps the whole library uniform and fast.
build_library <- function(D, X, y, tr, kind) {
  p <- ncol(D)
  sets <- lapply(seq_len(p), function(j) j)               # univariate on axes
  if (kind %in% c("pair", "pc")) {
    r <- apply(D[tr, , drop = FALSE], 2, function(z)
      suppressWarnings(abs(stats::cor(z, y[tr], method = "spearman"))))
    r[!is.finite(r)] <- 0
    top <- order(r, decreasing = TRUE)[seq_len(min(10L, p))]
    prs <- utils::combn(top, 2, simplify = FALSE)
    sets <- c(sets, prs)
  }
  extra <- NULL
  if (kind == "pc") {
    pc <- tryCatch(stats::prcomp(X[tr, , drop = FALSE], center = TRUE),
                   error = function(e) NULL)
    if (!is.null(pc)) {
      k <- min(8L, ncol(pc$rotation))
      Z <- scale(X, center = pc$center, scale = FALSE) %*% pc$rotation[, seq_len(k), drop = FALSE]
      extra <- Z
    }
  }
  list(sets = sets, extra = extra)
}

#' out-of-fold predictions from one learner set on rows `rows`
.ols_pred <- function(M, y, itr, ite) {
  Xt <- cbind(1, M[itr, , drop = FALSE])
  b <- tryCatch(stats::lm.fit(Xt, y[itr])$coefficients, error = function(e) NULL)
  if (is.null(b)) return(rep(mean(y[itr]), length(ite)))
  b[!is.finite(b)] <- 0
  as.numeric(cbind(1, M[ite, , drop = FALSE]) %*% b)
}

make_simple_sl <- function(kind = c("uni", "pair", "pc"), meta = c("mse", "rank"),
                           inner = 5L) {
  kind <- match.arg(kind); meta <- match.arg(meta)
  function(tr, te, y, X, D, aux) {
    lib <- build_library(D, X, y, tr, kind)
    mats <- lapply(lib$sets, function(s) D[, s, drop = FALSE])
    if (!is.null(lib$extra))
      mats <- c(mats, lapply(seq_len(ncol(lib$extra)),
                             function(j) lib$extra[, j, drop = FALSE]))
    mats <- c(mats, list(matrix(0, nrow(D), 1)))   # intercept-only candidate
    K <- length(mats)
    if (K < 2 || length(tr) < 12) return(rep(mean(y[tr]), length(te)))

    # inner CV over the training rows to get honest candidate predictions
    f <- rep_len(seq_len(min(inner, length(tr))), length(tr))
    Z <- matrix(NA_real_, length(tr), K)
    for (j in unique(f)) {
      itr <- tr[f != j]; ite <- tr[f == j]
      if (length(itr) < 6) next
      for (k in seq_len(K)) Z[f == j, k] <- .ols_pred(mats[[k]], y, itr, ite)
    }
    ok <- stats::complete.cases(Z)
    if (sum(ok) < 8) return(rep(mean(y[tr]), length(te)))
    ytr <- y[tr]

    w <- if (meta == "mse") {
      nnls_w(Z[ok, , drop = FALSE], ytr[ok])
    } else {
      # rank meta-learner: weight candidates by their inner-CV rank correlation,
      # non-negative and normalised (SL-02's recipe, applied to this library)
      rr <- apply(Z[ok, , drop = FALSE], 2, function(z)
        suppressWarnings(stats::cor(z, ytr[ok], method = "spearman")))
      rr[!is.finite(rr)] <- 0
      rr <- pmax(rr, 0)
      if (sum(rr) <= 0) rep(1 / K, K) else rr / sum(rr)
    }

    # refit every candidate on the full training set, combine with those weights
    P <- vapply(seq_len(K), function(k) .ols_pred(mats[[k]], y, tr, te),
                numeric(length(te)))
    if (is.null(dim(P))) P <- matrix(P, nrow = length(te))
    as.numeric(P %*% w)
  }
}

ARMS <- c(exp_baseline_arms()[c("null_train_mean", "domain_index")],
          list(sl_uni      = make_simple_sl("uni",  "mse"),
               sl_uni_pair = make_simple_sl("pair", "mse"),
               sl_pc       = make_simple_sl("pc",   "mse"),
               sl_rank     = make_simple_sl("pair", "rank")))

rows <- list()
for (i in seq_len(nrow(ix))) {
  for (tg in c("level", "prev")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tg),
                     error = function(e) NULL)
    if (is.null(cell)) next
    rows[[paste(i, tg)]] <- exp_infill(cell, ARMS, reps = REPS)
  }
  message("  ", ix$country[i], " ", ix$outcome[i])
}

loco <- list()
for (tg in c("level", "prev")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tg, outcomes = on)
    if (length(cl) < 3) next
    loco[[paste(tg, on)]] <- exp_loco(cl, ARMS, domain_of = E$domain_of)
    message("  LOCO ", tg, " ", on)
  }
}

SM <- exp_summarise(dplyr::bind_rows(rows)); exp_write(SM, "20_simple_library_sl")
LC <- dplyr::bind_rows(loco); exp_write(LC, "20_simple_library_sl_loco")

cat("\n== in-fill: median Spearman over cells ==\n")
a <- aggregate(spearman ~ target + arm, data = SM[SM$estimand == "infill", ],
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$target, -a$spearman), ], row.names = FALSE)

cat("\n== transport (LOCO): mean over held-out cells ==\n")
b <- aggregate(spearman ~ target + arm, data = LC,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
print(b[order(b$target, -b$spearman), ], row.names = FALSE)

cat("\n== paired against the index, blocks = country x target ==\n")
pair <- function(df, arm, ref, lbl) {
  d <- df[df$arm %in% c(arm, ref), c("country", "outcome", "target", "arm", "spearman")]
  w <- reshape(d, idvar = c("country", "outcome", "target"), timevar = "arm",
               direction = "wide")
  names(w) <- sub("^spearman[.]", "", names(w))
  g <- w[[arm]] - w[[ref]]; ok <- is.finite(g)
  bk <- paste(w$country[ok], w$target[ok])
  bl <- stats::aggregate(list(g = g[ok]), by = list(b = bk), FUN = mean)
  cat(sprintf("  %-8s %-12s mean %+.4f  cells %2d/%2d  blocks %d/%d  block p=%.4f\n",
      lbl, arm, mean(g[ok]), sum(g[ok] > 0), sum(ok),
      sum(bl$g > 0), nrow(bl),
      stats::binom.test(sum(bl$g > 0), nrow(bl), 0.5)$p.value))
}
for (a2 in c("sl_uni", "sl_uni_pair", "sl_pc", "sl_rank"))
  pair(SM[SM$estimand == "infill", ], a2, "domain_index", "infill")
for (a2 in c("sl_uni", "sl_uni_pair", "sl_pc", "sl_rank"))
  pair(LC, a2, "domain_index", "transport")
