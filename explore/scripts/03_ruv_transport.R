# =============================================================================
# explore/scripts/03_ruv_transport.R   [probe RV-01]
#
# QUESTION. The project's transport estimand is a RANK claim only, because
# biomarker levels carry large cross-survey offsets (raw ferritin 6x across
# countries; AS-01 shows it is not the assay, and it is not the adjustment
# method). Standardising the outcome within country removes the offset by
# throwing the level away.
#
# Genomics has a name for this and a different answer: it is a BATCH EFFECT,
# and RUV / SVA remove it by estimating the unwanted-variation SUBSPACE and
# projecting it out, rather than discarding the signal that shares a scale
# with it.
#
# THE PREDICTOR-SIDE VERSION (what this probe tests). Rank-normalising each
# column within country (the project's fix 3) removes each column's country
# MEAN and SCALE. It does not remove latent MULTIVARIATE country structure -
# directions in predictor space along which countries differ systematically
# for reasons of raster vintage, survey year or agro-ecological regime rather
# than nutrition. Those directions are precisely what a model fitted on three
# countries extrapolates along when it meets a fourth.
#
#   ruv_bc_k   project out the top-k BETWEEN-COUNTRY DISCRIMINANT directions,
#              estimated from TRAINING countries only, using country labels and
#              never the outcome. Keeps within-country district variation,
#              removes what makes countries look different.
#   ruv_pc_k   project out the top-k principal components of the pooled
#              training matrix (the unsupervised comparator, to separate
#              "remove country structure" from "remove the biggest directions")
#
# THE OUTCOME-SIDE MEASUREMENT. Also reports how much of the level variance is
# between-country, per outcome - the size of the thing being removed.
#
#   Rscript explore/scripts/03_ruv_transport.R
# -> explore/out/03_ruv_loco.csv, 03_ruv_level_variance.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

E <- exp_load()
DOM <- E$domain_of
CLIM <- "Climate and weather"; SOIL <- "Soil characteristics"

# ── projections ─────────────────────────────────────────────────────────────

#' Basis of the top-k between-country discriminant directions
#'
#' Between-country scatter B = sum_c n_c (mu_c - mu)(mu_c - mu)', whitened by
#' the pooled within-country scatter. With 4 countries, B has rank <= 3, so
#' k is capped at (number of training countries - 1).
#'
#' Estimated from TRAINING rows and their country labels only. The outcome is
#' never used, so this is a predictor-side transform, not a fitted model.
bc_basis <- function(X, ctry, k = 1L, ridge = 1e-3) {
  cs <- unique(ctry)
  k <- min(k, length(cs) - 1L)
  if (k < 1) return(NULL)
  p <- ncol(X)
  mu <- colMeans(X)
  Bm <- matrix(0, p, p)
  W <- matrix(0, p, p)
  for (c in cs) {
    idx <- which(ctry == c)
    if (length(idx) < 2) next
    mc <- colMeans(X[idx, , drop = FALSE])
    d <- mc - mu
    Bm <- Bm + length(idx) * tcrossprod(d)
    Xc <- sweep(X[idx, , drop = FALSE], 2, mc, "-")
    W <- W + crossprod(Xc)
  }
  W <- W / max(nrow(X) - length(cs), 1)
  W <- W + diag(ridge * mean(diag(W)) + 1e-8, p)
  Wi <- tryCatch(chol2inv(chol(W)), error = function(e) NULL)
  if (is.null(Wi)) return(NULL)
  eg <- eigen(Wi %*% Bm, symmetric = FALSE)
  V <- Re(eg$vectors[, seq_len(k), drop = FALSE])
  qr.Q(qr(V))                       # orthonormal basis of the subspace
}

#' Project a matrix onto the orthogonal complement of a basis
project_out <- function(X, Bs) {
  if (is.null(Bs) || !ncol(Bs)) return(X)
  X - (X %*% Bs) %*% t(Bs)
}

#' Arm factory: strip a subspace from X, rebuild the domain axes, then run the
#' given inner arm. The subspace is learned on the TRAINING rows only.
make_stripped_arm <- function(how = c("bc", "pc"), k = 1L,
                              inner = arm_domain_index_v2, doms = NULL) {
  how <- match.arg(how)
  function(tr, te, y, X, D, aux) {
    cc <- if (is.null(doms)) colnames(X) else
      colnames(X)[which(DOM[colnames(X)] %in% doms)]
    if (length(cc) < 5) return(rep(mean(y[tr]), length(te)))
    Xs <- X[, cc, drop = FALSE]
    Bs <- if (how == "bc") {
      bc_basis(Xs[tr, , drop = FALSE], aux$country[tr], k = k)
    } else {
      pc <- try(stats::prcomp(Xs[tr, , drop = FALSE], center = TRUE), silent = TRUE)
      if (inherits(pc, "try-error")) NULL else
        pc$rotation[, seq_len(min(k, ncol(pc$rotation))), drop = FALSE]
    }
    Xp <- project_out(Xs, Bs)
    D2 <- domain_representation_v2(Xp, DOM, sign_rows = tr)
    if (!ncol(D2)) return(rep(mean(y[tr]), length(te)))
    inner(tr, te, y, Xp, D2, aux)
  }
}

CS <- c(CLIM, SOIL)
ARMS <- list(
  null_train_mean = ARMS_V2$null_train_mean,
  index           = ARMS_V2$domain_index,
  index_cs        = function(tr, te, y, X, D, aux) {
    cc <- colnames(X)[which(DOM[colnames(X)] %in% CS)]
    if (length(cc) < 5) return(rep(mean(y[tr]), length(te)))
    D2 <- domain_representation_v2(X[, cc, drop = FALSE], DOM, sign_rows = tr)
    arm_domain_index_v2(tr, te, y, X, D2, aux)
  },
  ruv_bc_k1 = make_stripped_arm("bc", 1L),
  ruv_bc_k2 = make_stripped_arm("bc", 2L),
  ruv_pc_k1 = make_stripped_arm("pc", 1L),
  ruv_pc_k3 = make_stripped_arm("pc", 3L),
  ruv_bc_k1_cs = make_stripped_arm("bc", 1L, doms = CS),
  ruv_bc_k2_cs = make_stripped_arm("bc", 2L, doms = CS)
)

# ── transport, the estimand this probe is about ─────────────────────────────
ix <- exp_cell_index(E)
loco <- list()
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tgt, outcomes = on)
    if (length(cl) < 3) next
    loco[[paste(tgt, on)]] <- exp_loco(cl, ARMS, domain_of = DOM)
    message("  LOCO ", tgt, " ", on)
  }
}
lc <- dplyr::bind_rows(loco); exp_write(lc, "03_ruv_loco")

# ── the within-country control ──────────────────────────────────────────────
# Removing between-country structure should be roughly NEUTRAL within a
# country. If it improves both, it is absorbing signal, not batch.
rows <- list()
for (i in seq_len(nrow(ix))) {
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tgt),
                     error = function(e) NULL)
    if (is.null(cell)) next
    cell$aux$country <- rep(cell$country, cell$n)   # single country: bc is a no-op
    rows[[paste(i, tgt)]] <- exp_infill(
      cell, ARMS[c("null_train_mean", "index", "index_cs", "ruv_pc_k1", "ruv_pc_k3")],
      reps = as.integer(Sys.getenv("EXP_REPS", "10")))
  }
  message("  in-fill ", ix$country[i], " ", ix$outcome[i])
}
sm <- exp_summarise(dplyr::bind_rows(rows)); exp_write(sm, "03_ruv_infill")

# ── how big is the thing being removed? ─────────────────────────────────────
TG <- E$TG
vr <- list()
for (on in unique(TG$outcome)) {
  t <- TG[TG$outcome == on & is.finite(TG$y_level), ]
  if (dplyr::n_distinct(t$country) < 2) next
  gm <- tapply(t$y_level, t$country, mean)
  vb <- stats::var(gm)
  vw <- mean(tapply(t$y_level, t$country, stats::var), na.rm = TRUE)
  vr[[on]] <- data.frame(outcome = on, countries = dplyr::n_distinct(t$country),
                         between_country_var = vb, mean_within_var = vw,
                         between_share = vb / (vb + vw),
                         range_of_country_means = diff(range(gm)),
                         mean_within_sd = sqrt(vw), stringsAsFactors = FALSE)
}
VR <- dplyr::bind_rows(vr); exp_write(VR, "03_ruv_level_variance")

cat("\n== transport (LOCO): mean Spearman over held-out cells ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
b$positive <- aggregate(spearman ~ target + arm, data = lc,
                        FUN = function(z) sum(z > 0, na.rm = TRUE))$spearman
b$cells <- aggregate(spearman ~ target + arm, data = lc,
                     FUN = function(z) sum(is.finite(z)))$spearman
print(b[order(b$target, -b$spearman), ], row.names = FALSE)

cat("\n== within-country control (in-fill, level): should be roughly neutral ==\n")
a <- sm[sm$estimand == "infill" & sm$target == "level", ]
print(aggregate(spearman ~ arm, data = a,
                FUN = function(z) round(median(z, na.rm = TRUE), 3)), row.names = FALSE)

cat("\n== between-country share of level variance, per outcome ==\n")
print(VR[order(-VR$between_share),
         c("outcome", "countries", "between_share", "range_of_country_means",
           "mean_within_sd")], row.names = FALSE, digits = 3)
