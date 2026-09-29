# =============================================================================
# scripts/policy_deck/25_mnf15_v2_ghana_heldout.R
#
# "Without a single Ghanaian blood sample": per-district predictions for Ghana
# from the index trained on The Gambia, Sierra Leone and Malawi only, for the
# v2 MNF15 talk (docs/slides/MNF15-talk-2026-09-v2.qmd, slide 6).
#
# This is estimand C of scripts/protocol_v2/02b_merge_and_loco.R run for ONE
# held-out country and ONE cell, with nothing changed: the same predictor
# policy (tiers open + survey_public, the transport headline), the same
# within-country rank normalisation, the same pooled domain components oriented
# on the training countries only, the same domain_index arm and score_v2. The
# only difference is that it keeps the per-district predictions, with district
# names, which 02b discards. The check below requires it to reproduce the
# committed benchmarks_v2_cells.csv value (0.5537) before anything is written.
#
# Two outputs:
#   ghana_heldout_<outcome>_<target>.csv      the 75 surveyed districts: survey
#                                             value, prediction, ranks (scored)
#   ghana_heldout_<outcome>_<target>_all.csv  every Ghana district: the same
#                                             trained model applied to all 260
#                                             (rank-normalised over all 260, as
#                                             deployment would), for the map only
#
#   Rscript scripts/policy_deck/25_mnf15_v2_ghana_heldout.R
#   V2H_OUTCOME=women_iron V2H_TARGET=level Rscript scripts/policy_deck/25_mnf15_v2_ghana_heldout.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")   # as 02b (TP-01)
source("R/protocol_v2.R")

OUTT <- "results/tables/policy_deck"; dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)
P2   <- "results/tables/protocol_v2"
ON   <- Sys.getenv("V2H_OUTCOME", "child_iron")
TGT  <- Sys.getenv("V2H_TARGET", "level")
HELD <- "Ghana"

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")

# as 02b's build_cell (no centroids needed: domain_index does not read them), plus the names
build_cell <- function(cn) {
  t <- TG[TG$country == cn & TG$outcome == ON, ]
  ycol <- if (TGT == "prev") "y_prev" else "y_level"
  wcol <- if (TGT == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y_nat <- m[[ycol]]
  list(country = cn, n = nrow(m), y_nat = y_nat,
       y_mod = if (TGT == "prev") .v2_logit(y_nat) else y_nat,
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1, Admin2 = m$Admin2, n_psu = m$n_psu)
}
cl <- Filter(Negate(is.null), setNames(lapply(COUNTRIES, build_cell), COUNTRIES))
stopifnot(HELD %in% names(cl), length(cl) >= 3)
common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
ynat <- unlist(lapply(cl, function(z) z$y_nat))
wv   <- unlist(lapply(cl, function(z) z$w))
aux  <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
te <- which(ctry == HELD); tr <- which(ctry != HELD)
Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
pred <- ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, aux)
s <- score_v2(ynat[te], pred, wv[te], scale = if (TGT == "prev") "prev" else "level")

# the committed number this must reproduce
CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"))
ref <- CELL$spearman[CELL$country == HELD & CELL$outcome == ON & CELL$target == TGT &
                     CELL$estimand == "country" & CELL$arm == "domain_index"]
cat(sprintf("%s %s %s held out: spearman %.4f (committed %.4f), %d districts, %d common columns, trained on %s\n",
            HELD, ON, TGT, s$spearman, ref, length(te), length(common), paste(setdiff(names(cl), HELD), collapse = ", ")))
if (length(ref) != 1 || abs(s$spearman - ref) > 0.005) stop("does not reproduce benchmarks_v2_cells.csv; nothing written")

g <- cl[[HELD]]
d <- data.frame(Admin1 = g$Admin1, Admin2 = g$Admin2, n_psu = g$n_psu, y_obs = g$y_nat, pred = pred)
d$r_survey <- rank(-d$y_obs, ties.method = "average")   # 1 = worst (higher y = worse status, as in fig7)
d$r_model  <- rank(-d$pred,  ties.method = "average")
tag <- sprintf("%s_%s", ON, TGT)
write.csv(d, file.path(OUTT, sprintf("ghana_heldout_%s.csv", tag)), row.names = FALSE)
writeLines(sprintf("%.4f", s$spearman), file.path(OUTT, sprintf("ghana_heldout_%s_rho.txt", tag)))

# ---- deployment view: the same training rows, Ghana rank-normalised over ALL its districts
all_g <- S[S$country == HELD, c("Admin1", "Admin2", PREDS)]
Xg <- prep_predictors_v2(as.matrix(all_g[, PREDS, drop = FALSE]))
cc <- intersect(common, colnames(Xg))
Xtr <- do.call(rbind, lapply(cl[setdiff(names(cl), HELD)], function(z) z$X[, cc, drop = FALSE]))
Ytr <- unlist(lapply(cl[setdiff(names(cl), HELD)], function(z) as.numeric(scale(z$y_mod))))
Xa  <- rbind(Xtr, Xg[, cc, drop = FALSE])
tr2 <- seq_len(nrow(Xtr)); te2 <- nrow(Xtr) + seq_len(nrow(Xg))
Ya  <- c(Ytr, rep(NA_real_, nrow(Xg)))
Da  <- domain_representation_v2(Xa, domain_of, sign_rows = tr2)
pa  <- ARMS_V2[["domain_index"]](tr2, te2, Ya, Xa, Da, list(y_nat = Ya))
A <- data.frame(Admin1 = all_g$Admin1, Admin2 = all_g$Admin2, pred_all = pa)
A$q_all <- 100 * (rank(-A$pred_all, ties.method = "average") - 0.5) / nrow(A)
chk <- inner_join(d, A, by = c("Admin1", "Admin2"))
cat(sprintf("deployment view: %d districts, %d common columns; on the %d surveyed districts it agrees with the scored fit at %.3f and with the survey at %.3f\n",
            nrow(A), length(cc), nrow(chk), cor(chk$pred, chk$pred_all, method = "spearman"),
            cor(chk$y_obs, chk$pred_all, method = "spearman")))
write.csv(A, file.path(OUTT, sprintf("ghana_heldout_%s_all.csv", tag)), row.names = FALSE)
cat("wrote", sprintf("ghana_heldout_%s.csv and _all.csv", tag), "\n")
