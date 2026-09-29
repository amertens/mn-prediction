# =============================================================================
# explore/scripts/23_plasmode_countries.R   [probe PL-01]
#
# WOULD 20 COUNTRIES FIX THIS, OR IS THE DISTRICT GROUND TRUTH THE PROBLEM?
#
# A plasmode: real covariates, imposed truth. The real pooled Admin-2 domain
# axes are kept exactly as they are - their collinearity and spatial structure
# are the hard part to simulate and are the reason for using a plasmode at all -
# and only the OUTCOME is generated, from a known model.
#
# WHAT THIS CAN AND CANNOT ANSWER, stated before any result.
#
#   CAN:  the SHAPE of the country curve. TC-02 measured +0.05 per country over
#         1 -> 2 -> 3, which extrapolated naively gives 1.10 at twenty countries,
#         so the curve must bend somewhere and the data cannot say where.
#   CAN:  which of the two constraints binds, because a simulation can break
#         them independently and the real data never can - hold measurement
#         error fixed and vary countries, then hold countries fixed and vary
#         measurement error.
#   CANNOT: what twenty REAL countries would give. The quantity that governs
#         transport is how the predictor-outcome mapping varies BETWEEN
#         countries, and with four countries there are three degrees of freedom
#         to estimate it. Every synthetic country here is drawn from a
#         between-country distribution ESTIMATED ON THOSE FOUR. If a twenty-
#         country pool spanned wider agro-ecologies or food systems, its
#         heterogeneity would be larger than anything simulated here and the
#         curve would sit lower. This is an if-then, not a forecast.
#
# THE ASSUMPTION DOING THE WORK, named so it can be attacked: Sigma_beta, the
# between-country covariance of the index weights, estimated from four fitted
# weight vectors. Everything about the country dimension follows from it.
#
# DESIGN
#   truth      y_c = D_c (beta_bar + delta_c) + country_offset_c,  delta_c ~ N(0, Sigma_beta)
#   observed   y_obs = y_true + e,  e ~ N(0, noise_mult * sigma2_samp), sigma2_samp
#              taken from each district's own measured n_eff (DE-01)
#   districts  sampled from the real 206 without replacement within a country,
#              with a small multivariate jitter so twenty countries are not
#              twenty copies of the same rows; overlap across countries is
#              reported because it cannot be removed
#   calibrate  the signal scale is set so that 4 countries at real noise
#              reproduces the observed transport Spearman of about 0.30
#   factorial  countries in {3,4,6,10,20} x noise multiplier in {0, 0.5, 1, 2}
#
#   Rscript explore/scripts/23_plasmode_countries.R
#   EXP_PL_REPS=40
# -> explore/out/23_plasmode.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_PL_REPS", "40"))
set.seed(20260929L)
E <- exp_load()

# ── real pooled axes and the real district design ───────────────────────────
cells <- exp_all_cells(E, "level", outcomes = "child_vitA")
common <- Reduce(intersect, lapply(cells, function(z) colnames(z$X)))
Xr <- do.call(rbind, lapply(cells, function(z) z$X[, common, drop = FALSE]))
ctry0 <- rep(vapply(cells, function(z) z$country, ""), vapply(cells, function(z) z$n, 0L))
D0 <- domain_representation_v2(Xr, E$domain_of, sign_rows = seq_len(nrow(Xr)))
message("real pool: ", nrow(D0), " districts x ", ncol(D0), " axes, ",
        length(unique(ctry0)), " countries")

# district sampling variance on the level scale, from the real design
TG <- E$TG
sv_pool <- local({
  s <- c()
  for (cn in unique(ctry0)) {
    t <- TG[TG$country == cn & TG$outcome == "child_vitA", ]
    v <- (t$sd_level^2) / pmax(t$n_eff_cont, 1)
    s <- c(s, v[is.finite(v)])
  }
  s[s > 0]
})
message("real district sampling variance: median ", signif(median(sv_pool), 3))

# ── between-country heterogeneity of the index weights ──────────────────────
# Sigma_beta from the four real countries' own fitted weights (diagonal: with
# four countries a full covariance is not estimable).
W <- do.call(rbind, lapply(unique(ctry0), function(cn) {
  k <- which(ctry0 == cn)
  y <- unlist(lapply(cells, function(z) z$y_mod))[k]
  .index_weights_v2(D0[k, , drop = FALSE], y)
}))
beta_bar <- colMeans(W)
# The raw across-country sd of the fitted weights is NOT the between-country
# heterogeneity: it is that plus the sampling noise of each country's own
# weight estimate. The index weights are fisher-z * sqrt(n-3), which have unit
# sampling variance by construction, so the moment correction is
# tau^2 = max(var - 1, 0) - the same estimator EB-01 used. Both are run,
# because with four countries the corrected value is itself barely identified
# and the two bracket the answer.
sd_raw  <- apply(W, 2, stats::sd)
sd_corr <- sqrt(pmax(apply(W, 2, stats::var) - 1, 0))
HET <- Sys.getenv("EXP_PL_HET", "corrected")
sd_beta <- if (HET == "raw") sd_raw else sd_corr
message("heterogeneity estimator: ", HET,
        " | mean sd raw ", signif(mean(sd_raw), 3),
        " vs corrected ", signif(mean(sd_corr), 3),
        " | axes with tau2 floored to 0: ", sum(sd_corr == 0), " of ", length(sd_corr))
message("index weights: ", ncol(W), " axes; between-country sd/mean = ",
        signif(mean(sd_beta) / max(mean(abs(beta_bar)), 1e-9), 3))

N_DIST <- c(30, 75, 87, 14)   # the real district counts

#' one synthetic world
simulate <- function(n_ctry, noise_mult, signal, jitter = 0.25) {
  Dl <- list(); yl <- list(); cl <- c()
  for (c0 in seq_len(n_ctry)) {
    nd <- sample(N_DIST, 1)
    idx <- sample(nrow(D0), min(nd, nrow(D0)))
    Dc <- D0[idx, , drop = FALSE]
    # jitter so many countries are not many copies of the same rows
    Dc <- Dc + matrix(stats::rnorm(length(Dc), 0, jitter), nrow(Dc))
    delta <- stats::rnorm(ncol(D0), 0, sd_beta)
    lin <- as.numeric(Dc %*% (beta_bar + delta))
    lin <- (lin - mean(lin)) / max(stats::sd(lin), 1e-9)
    # `signal` is the CORRELATION between the linear predictor and the true
    # district value: the rest is irreducible district heterogeneity the
    # covariates cannot reach. Without this term the truth is a deterministic
    # function of D and every signal level scores ~0.87, which is what the
    # first version of this script did.
    y_true <- signal * lin +
      sqrt(max(1 - signal^2, 0)) * stats::rnorm(nrow(Dc)) +
      stats::rnorm(1, 0, 0.5)                            # + country offset
    e <- stats::rnorm(nrow(Dc), 0, sqrt(noise_mult * sample(sv_pool, nrow(Dc), TRUE)))
    Dl[[c0]] <- Dc; yl[[c0]] <- y_true + e; cl <- c(cl, rep(c0, nrow(Dc)))
    attr(Dl[[c0]], "truth") <- y_true
  }
  list(D = do.call(rbind, Dl), y = unlist(yl), ctry = cl,
       truth = unlist(lapply(Dl, attr, "truth")))
}

#' leave-one-country-out transport, scored against the TRUE district value
score_world <- function(w) {
  out <- c()
  for (h in unique(w$ctry)) {
    te <- which(w$ctry == h); tr <- which(w$ctry != h)
    if (length(tr) < 20 || length(te) < 8) next
    # outcomes standardised within country before pooling, as the protocol does
    ys <- w$y
    for (c0 in unique(w$ctry[tr])) {
      k <- tr[w$ctry[tr] == c0]; ys[k] <- as.numeric(scale(ys[k]))
    }
    p <- tryCatch(arm_domain_index_v2(tr, te, ys, w$D, w$D,
                                      list(target = "level")),
                  error = function(e) rep(NA_real_, length(te)))
    if (length(p) != length(te)) next
    out <- c(out, suppressWarnings(stats::cor(w$truth[te], p, method = "spearman")))
  }
  mean(out, na.rm = TRUE)
}

# ── calibrate the signal so 4 countries at real noise matches reality ───────
target_real <- 0.298
SIG_GRID <- c(0.15, 0.25, 0.35, 0.45, 0.55, 0.70)
cal <- sapply(SIG_GRID, function(s)
  mean(replicate(12, score_world(simulate(4, 1, s))), na.rm = TRUE))
sig <- SIG_GRID[which.min(abs(cal - target_real))]
message("calibration: signal grid ", paste(round(cal, 3), collapse = " "),
        " -> chose signal = ", sig, " (target ", target_real, ")")

# ── the factorial ───────────────────────────────────────────────────────────
grid <- expand.grid(n_ctry = c(3, 4, 6, 10, 20),
                    noise = c(0, 0.5, 1, 2))
res <- list()
for (i in seq_len(nrow(grid))) {
  v <- replicate(REPS, score_world(simulate(grid$n_ctry[i], grid$noise[i], sig)))
  res[[i]] <- data.frame(n_countries = grid$n_ctry[i], noise_mult = grid$noise[i],
                         spearman = mean(v, na.rm = TRUE),
                         se = stats::sd(v, na.rm = TRUE) / sqrt(sum(is.finite(v))),
                         reps = sum(is.finite(v)))
  message(sprintf("  countries %2d  noise %.1f -> %.3f",
                  grid$n_ctry[i], grid$noise[i], res[[i]]$spearman))
}
R <- dplyr::bind_rows(res); R$signal <- sig; R$het <- HET
exp_write(R, paste0("23_plasmode_", HET))

cat("\n== transport Spearman against the TRUE district value ==\n")
w <- reshape(R[, c("n_countries", "noise_mult", "spearman")],
             idvar = "n_countries", timevar = "noise_mult", direction = "wide")
names(w) <- sub("^spearman\\.", "noise_", names(w))
print(w, row.names = FALSE, digits = 3)

cat("\n== what each lever buys, from the simulation ==\n")
base <- R$spearman[R$n_countries == 4 & R$noise_mult == 1]
c20   <- R$spearman[R$n_countries == 20 & R$noise_mult == 1]
n0    <- R$spearman[R$n_countries == 4 & R$noise_mult == 0]
both  <- R$spearman[R$n_countries == 20 & R$noise_mult == 0]
cat(sprintf("  baseline, 4 countries at real noise      : %.3f\n", base))
cat(sprintf("  20 countries, real noise                 : %.3f  (%+.3f)\n", c20, c20 - base))
cat(sprintf("  4 countries, NO district noise           : %.3f  (%+.3f)\n", n0, n0 - base))
cat(sprintf("  20 countries AND no district noise       : %.3f  (%+.3f)\n", both, both - base))
cat("\n  (scored against the TRUE district value, so the noise arms are not\n")
cat("   merely scoring against a cleaner yardstick - the truth is fixed)\n")
