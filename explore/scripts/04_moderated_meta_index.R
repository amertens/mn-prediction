# =============================================================================
# explore/scripts/04_moderated_meta_index.R   [probe EB-01]
#
# QUESTION. The index sets its weights from ONE cell's 14-87 districts. The
# memory note fe_effective_n names exactly this as the root cause of weak and
# unstable prediction: the effective n is the number of AREAS, not individuals.
# But there are 24 cells sharing the same predictor vocabulary - the
# many-features x many-contrasts layout that empirical-Bayes moderation (limma,
# and random-effects meta-analysis) was built for.
#
# Does shrinking each cell's weights toward the cross-cell consensus beat
# letting each cell estimate its own?
#
# THE METHOD. For axis j in cell c, the index's own statistic is
# z_cj = fisher-z(spearman) * sqrt(n-3), which has unit variance by
# construction. So the hierarchy is textbook:
#       z_cj ~ N(theta_cj, 1),   theta_cj ~ N(mu_j, tau_j^2)
# and the posterior mean is  mu_j + (tau^2/(tau^2+1)) (z_cj - mu_j).
# mu_j and tau_j^2 are estimated from the OTHER cells, never from cell c, so
# the shrinkage target cannot contain the cell it is applied to.
#
# THE NESTING THAT MATTERS. Cells of the same country share districts: Ghana
# child_vitA and Ghana women_iron are the same 75 units. Using a companion
# cell's full data to weight a prediction for a held-out Ghana district would
# read that district's own survey through another biomarker. So companion rows
# are always filtered to exclude the held-out district keys. Other countries'
# cells contribute in full. That matches what a deployment actually has: the
# other biomarkers ARE measured in the surveyed districts, and are not measured
# in the unsurveyed one.
#
#   Rscript explore/scripts/04_moderated_meta_index.R
# -> explore/out/04_moderated_cells.csv, 04_moderated_loco.csv,
#    04_moderated_axis_stats.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
ix <- exp_cell_index(E)

#' Fisher-z statistic per axis, unit variance by construction
axis_z <- function(D, y) {
  n <- length(y)
  out <- apply(D, 2, function(x) {
    if (stats::sd(x) == 0) return(0)
    r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
    if (!is.finite(r)) return(0)
    r <- max(min(r, 0.999), -0.999)
    0.5 * log((1 + r) / (1 - r)) * sqrt(max(n - 3, 1))
  })
  out[!is.finite(out)] <- 0
  out
}

#' Build every cell once per target, keeping district keys for the nesting
build_all <- function(tgt) {
  cs <- list()
  for (i in seq_len(nrow(ix))) {
    cc <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tgt),
                   error = function(e) NULL)
    if (is.null(cc)) next
    cc$key <- paste(cc$country, cc$Admin1, cc$Admin2, sep = "|")
    cs[[paste(ix$country[i], ix$outcome[i], sep = "|")]] <- cc
  }
  cs
}

#' Empirical-Bayes moderated weights for one cell, given companions
#'
#' @param own_z    the cell's own z on its training rows
#' @param comp     list of companion cells
#' @param drop_key district keys that must not contribute (the held-out fold)
moderated_weights <- function(own_z, comp, drop_key, axes) {
  Z <- list()
  for (cc in comp) {
    keep <- which(!(cc$key %in% drop_key))
    if (length(keep) < 12) next
    a <- intersect(axes, colnames(cc$D))
    if (length(a) < 1) next
    z <- axis_z(cc$D[keep, a, drop = FALSE], cc$y_mod[keep])
    v <- stats::setNames(rep(NA_real_, length(axes)), axes)
    v[a] <- z
    Z[[length(Z) + 1L]] <- v
  }
  if (!length(Z)) return(own_z)
  M <- do.call(rbind, Z)
  mu <- colMeans(M, na.rm = TRUE)
  # tau^2 = observed between-cell variance minus the sampling variance (1),
  # floored at 0: the standard method-of-moments estimator
  vv <- apply(M, 2, stats::var, na.rm = TRUE)
  tau2 <- pmax(vv - 1, 0)
  mu[!is.finite(mu)] <- 0; tau2[!is.finite(tau2)] <- 0
  shrink <- tau2 / (tau2 + 1)
  out <- mu[axes] + shrink[axes] * (own_z[axes] - mu[axes])
  out[!is.finite(out)] <- 0
  out
}

#' The moderated index as an arm, closing over the companion cells
make_meta_index <- function(cells, self_name) {
  comp <- cells[setdiff(names(cells), self_name)]
  self <- cells[[self_name]]
  function(tr, te, y, X, D, aux) {
    axes <- colnames(D)
    own <- axis_z(D[tr, , drop = FALSE], y[tr])
    w <- moderated_weights(own, comp, drop_key = self$key[te], axes = axes)
    itr <- as.numeric(D[tr, , drop = FALSE] %*% w)
    ite <- as.numeric(D[te, , drop = FALSE] %*% w)
    if (stats::sd(itr) == 0) return(rep(mean(y[tr]), length(te)))
    ((ite - mean(itr)) / stats::sd(itr)) * stats::sd(y[tr]) + mean(y[tr])
  }
}

# ── in-country ──────────────────────────────────────────────────────────────
rows <- list()
for (tgt in c("prev", "level")) {
  cells <- build_all(tgt)
  message("built ", length(cells), " cells for target ", tgt)
  for (nm in names(cells)) {
    cell <- cells[[nm]]
    arms <- c(exp_baseline_arms(),
              list(meta_index = make_meta_index(cells, nm)))
    rows[[paste(nm, tgt, "A")]] <- exp_infill(cell, arms, reps = REPS)
    rows[[paste(nm, tgt, "B")]] <- exp_region(cell, arms)
    message("  ", nm, " ", tgt)
  }
}
sm <- exp_summarise(dplyr::bind_rows(rows)); exp_write(sm, "04_moderated_cells")

# ── transport ───────────────────────────────────────────────────────────────
# Under LOCO the companions are the cells of the TRAINING countries only; the
# held-out country contributes nothing, which is enforced by dropping every key
# of the held-out country.
loco <- list()
for (tgt in c("prev", "level")) {
  cells <- build_all(tgt)
  for (on in unique(ix$outcome)) {
    cl <- cells[grepl(paste0("\\|", on, "$"), names(cells))]
    if (length(cl) < 3) next
    comp <- cells[setdiff(names(cells), names(cl))]
    keys <- unlist(lapply(cl, function(z) z$key))
    ctry_of_key <- sub("\\|.*$", "", keys)
    meta_arm <- function(tr, te, y, X, D, aux) {
      axes <- colnames(D)
      own <- axis_z(D[tr, , drop = FALSE], y[tr])
      held <- unique(aux$country[te])
      drop <- unlist(lapply(comp, function(z) z$key[sub("\\|.*$", "", z$key) %in% held]))
      w <- moderated_weights(own, comp, drop_key = drop, axes = axes)
      itr <- as.numeric(D[tr, , drop = FALSE] %*% w)
      ite <- as.numeric(D[te, , drop = FALSE] %*% w)
      if (stats::sd(itr) == 0) return(rep(mean(y[tr]), length(te)))
      ((ite - mean(itr)) / stats::sd(itr)) * stats::sd(y[tr]) + mean(y[tr])
    }
    loco[[paste(tgt, on)]] <- exp_loco(
      cl, list(null_train_mean = ARMS_V2$null_train_mean,
               domain_index = ARMS_V2$domain_index, meta_index = meta_arm),
      domain_of = E$domain_of)
    message("  LOCO ", tgt, " ", on)
  }
}
lc <- dplyr::bind_rows(loco); exp_write(lc, "04_moderated_loco")

# ── descriptive: which RAW predictors survive cross-cell moderation ─────────
# Not a prediction, a screen. Reported so the annotation work has a ranked list.
cells <- build_all("level")
PRED <- Reduce(intersect, lapply(cells, function(z) colnames(z$X)))
ZM <- do.call(rbind, lapply(cells, function(z) axis_z(z$X[, PRED, drop = FALSE], z$y_mod)))
mu <- colMeans(ZM); vv <- apply(ZM, 2, stats::var); tau2 <- pmax(vv - 1, 0)
se_mu <- sqrt((tau2 + 1) / nrow(ZM))
zbar <- mu / se_mu
p <- 2 * stats::pnorm(-abs(zbar))
AX <- data.frame(predictor = PRED, domain = E$domain_of[PRED],
                 cells = nrow(ZM), mean_z = mu, tau2 = tau2,
                 pooled_z = zbar, p = p, fdr = stats::p.adjust(p, "BH"),
                 sign_consistency = colMeans(sign(ZM) == matrix(sign(mu),
                   nrow(ZM), length(mu), byrow = TRUE)),
                 stringsAsFactors = FALSE)
AX <- AX[order(AX$fdr, -abs(AX$pooled_z)), ]
exp_write(AX, "04_moderated_axis_stats")

cat("\n== in-country: median Spearman over cells ==\n")
a <- aggregate(spearman ~ estimand + target + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$estimand, a$target, -a$spearman), ], row.names = FALSE)

cat("\n== transport (LOCO) ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
b$positive <- aggregate(spearman ~ target + arm, data = lc,
                        FUN = function(z) sum(z > 0, na.rm = TRUE))$spearman
print(b[order(b$target, -b$spearman), ], row.names = FALSE)

cat("\n== predictors surviving cross-cell moderation at FDR < 0.05 (top 25) ==\n")
print(head(AX[AX$fdr < 0.05, c("predictor", "domain", "pooled_z", "fdr",
                               "sign_consistency")], 25), row.names = FALSE, digits = 3)
cat("\ntotal at FDR < 0.05:", sum(AX$fdr < 0.05), "of", nrow(AX), "\n")
