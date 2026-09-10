# =============================================================================
# scripts/policy_deck/06_civ_rank_uncertainty.R
#
# HOW CERTAIN IS THE COTE D'IVOIRE RANKING?
#
# There is no stored method for uncertainty on a transported district RANKING.
# The conformal machinery in R/conformal.R produces intervals on a LEVEL, and
# for the production CIV output those intervals are degenerate (constant width,
# clamped at zero). So this script builds the uncertainty the estimand actually
# admits: resample the training set, refit, re-rank, and record how far each
# Cote d'Ivoire district moves.
#
#   A. BOOTSTRAP (B = 400). Resample training districts with replacement,
#      stratified within country so the country balance is preserved. Refit the
#      domain PCs (orientation from the resampled rows only) and the index, then
#      re-rank CIV's 33 districts.
#   B. LEAVE-ONE-TRAINING-COUNTRY-OUT (4 refits). Coarser, but answers a
#      different question: does the ranking depend on WHICH countries we learned
#      from, rather than on how many districts they happen to contain.
#
# WHAT THIS UNCERTAINTY COVERS: sampling variability in the training set.
# WHAT IT DOES NOT COVER: whether the model transports to Cote d'Ivoire at all.
# There is no CIV ground truth, so nothing internal to CIV can test that. The
# external bound on it is the held-out transport accuracy (0.37), not this.
#
#   Rscript scripts/policy_deck/06_civ_rank_uncertainty.R
# -> results/tables/policy_deck/civ_rank_uncertainty.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
set.seed(20260909L)

P2   <- "results/tables/protocol_v2"
OUTT <- "results/tables/policy_deck"; dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)
B    <- 400L
ON   <- "child_iron"

TG  <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
CIV <- readRDS("results/transportability/civ_canonical_admin2_full.rds")

PREDS <- intersect(MD$column[MD$domain %in% c("Climate and weather", "Soil characteristics")],
                   names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")

build_cell <- function(cn) {
  t <- TG[TG$country == cn & TG$outcome == ON, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)],
                  by = c("Admin1", "Admin2"))
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  list(country = cn, n = nrow(m), y = m$y_level, X = Xr)
}
cl <- Filter(Negate(is.null), stats::setNames(lapply(COUNTRIES, build_cell), COUNTRIES))

civX <- prep_predictors_v2(as.matrix(CIV[, intersect(PREDS, names(CIV)), drop = FALSE]))
pool <- c(lapply(cl, function(z) z$X), list(CoteDIvoire = civX))
common <- Reduce(intersect, lapply(pool, colnames))
cat(sprintf("common climate + soil columns: %d\n", length(common)))

Xm   <- do.call(rbind, c(lapply(cl, function(z) z$X[, common, drop = FALSE]),
                         list(civX[, common, drop = FALSE])))
ctry <- c(rep(names(cl), vapply(cl, function(z) z$n, 0L)), rep("CoteDIvoire", nrow(civX)))
Y    <- c(unlist(lapply(cl, function(z) as.numeric(scale(z$y)))), rep(NA_real_, nrow(civX)))
tr_all <- which(ctry != "CoteDIvoire")
te     <- which(ctry == "CoteDIvoire")
nD <- length(te)

rank_from <- function(tr) {
  Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  p  <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, Dm, NULL),
                 error = function(e) rep(NA_real_, length(te)))
  if (all(is.na(p))) return(rep(NA_real_, length(te)))
  rank(-p, ties.method = "average")
}

# ── A. bootstrap ─────────────────────────────────────────────────────────────
cat("bootstrapping", B, "resamples ")
R <- matrix(NA_real_, nrow = nD, ncol = B)
for (b in seq_len(B)) {
  tr <- unlist(lapply(names(cl), function(cn) {
    idx <- which(ctry == cn); sample(idx, length(idx), replace = TRUE) }))
  R[, b] <- rank_from(tr)
  if (b %% 50 == 0) cat(".")
}
cat(" done\n")

ok <- colSums(is.finite(R)) == nD
R  <- R[, ok, drop = FALSE]
cat(sprintf("usable resamples: %d of %d\n", ncol(R), B))

# ── B. leave-one-training-country-out ────────────────────────────────────────
L <- sapply(names(cl), function(drop_cn) rank_from(which(ctry != "CoteDIvoire" & ctry != drop_cn)))
cat("leave-one-training-country-out refits:", ncol(L), "\n")

worst_k <- ceiling(nD / 3)   # "worst third" = 11 of 33 districts
U <- data.frame(
  Admin1     = CIV$Admin1,
  Admin2     = CIV$Admin2,
  rank_med   = apply(R, 1, stats::median),
  rank_lo    = apply(R, 1, stats::quantile, 0.05),
  rank_hi    = apply(R, 1, stats::quantile, 0.95),
  p_worst3rd = rowMeans(R <= worst_k),
  loco_range = apply(L, 1, function(x) diff(range(x))),
  stringsAsFactors = FALSE)
U$rank_width <- U$rank_hi - U$rank_lo
U$width_pct  <- 100 * U$rank_width / nD
U <- U[order(U$rank_med), ]
write.csv(U, file.path(OUTT, "civ_rank_uncertainty.csv"), row.names = FALSE)

cat(sprintf("\n90%% rank interval width: median %.0f of %d districts (%.0f%% of the list)\n",
            stats::median(U$rank_width), nD, stats::median(U$width_pct)))
cat(sprintf("districts whose 90%% interval stays inside the worst third: %d of %d\n",
            sum(U$rank_hi <= worst_k), nD))
cat(sprintf("districts with P(worst third) >= 0.80: %d;  <= 0.20: %d\n",
            sum(U$p_worst3rd >= 0.8), sum(U$p_worst3rd <= 0.2)))
cat("\nten worst-ranked districts, with their 90% rank interval:\n")
print(head(U[, c("Admin1", "Admin2", "rank_med", "rank_lo", "rank_hi", "p_worst3rd")], 10),
      row.names = FALSE, digits = 2)
