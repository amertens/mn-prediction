# =============================================================================
# scripts/protocol_v2/50_hapc_smoke_data.R   [HP-01, data step]
#
# Exports the design matrices the hapc smoke test needs, on the protocol's
# own representation, so the Python fit is scored on exactly the rows,
# targets and folds the domain index is scored on:
#   (a) in-fill, one cell: Ghana child iron, 75 districts, domain PCs built
#       on all rows (sign orientation is irrelevant for a kernel method), and
#       the protocol's ten 5-fold district draws as fold ids;
#   (b) leave-one-country-out, child iron and child vitamin A, all four
#       countries pooled: within-country rank-normalised columns, domain PCs
#       from the common columns, country labels for the folds.
# Y is the modelling target (biomarker level, or logit prevalence) and y_nat
# the natural-scale target for scoring, plus the effective n weights.
#
#   Rscript scripts/protocol_v2/50_hapc_smoke_data.R
# -> results/tables/protocol_v2/hapc_smoke/*.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUT <- "results/tables/protocol_v2/hapc_smoke"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
TG <- read.csv("results/tables/protocol_v2/targets_v2.csv", stringsAsFactors = FALSE)
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)

cell <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  list(country = cn, Admin1 = m$Admin1, Admin2 = m$Admin2, n = nrow(m), y_nat = m[[ycol]], y_mod = if (target == "prev") .v2_logit(m[[ycol]]) else m[[ycol]], w = m[[wcol]], X = Xr)
}
# (a) in-fill: Ghana child iron, both targets, domain PCs on all rows, fold ids for ten draws
for (target in c("level", "prev")) {
  z <- cell("Ghana", "child_iron", target)
  D <- domain_representation_v2(z$X, domain_of, sign_rows = seq_len(z$n))
  folds <- sapply(1:10, function(r) make_folds_v2("kfold_district", z$n, k = 5, rep_id = r))
  out <- data.frame(country = z$country, Admin1 = z$Admin1, Admin2 = z$Admin2, y_mod = z$y_mod, y_nat = z$y_nat, w = z$w, folds, check.names = FALSE)
  names(out)[7:16] <- paste0("fold_rep", 1:10)
  write.csv(cbind(out, as.data.frame(D)), file.path(OUT, sprintf("infill_ghana_child_iron_%s.csv", target)), row.names = FALSE)
  cat(sprintf("in-fill Ghana child_iron %s: %d districts, %d domain PCs\n", target, z$n, ncol(D)))
}
# (b) LOCO, pooled: child iron and child vitamin A, level target
for (on in c("child_iron", "child_vitA")) for (target in c("level", "prev")) {
  cl <- list(); for (cn in c("Gambia", "Ghana", "Malawi", "SierraLeone")) { z <- tryCatch(cell(cn, on, target), error = function(e) NULL); if (!is.null(z) && z$n >= 12) cl[[cn]] <- z }
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  D <- domain_representation_v2(Xm, domain_of, sign_rows = seq_len(nrow(Xm)))
  out <- data.frame(country = ctry, Admin1 = unlist(lapply(cl, `[[`, "Admin1")), Admin2 = unlist(lapply(cl, `[[`, "Admin2")),
                    y_mod = unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))), y_nat = unlist(lapply(cl, `[[`, "y_nat")), w = unlist(lapply(cl, `[[`, "w")), check.names = FALSE)
  write.csv(cbind(out, as.data.frame(D)), file.path(OUT, sprintf("loco_%s_%s.csv", on, target)), row.names = FALSE)
  cat(sprintf("LOCO %s %s: %d districts in %d countries, %d common columns, %d domain PCs\n", on, target, nrow(out), length(cl), length(common), ncol(D)))
}
cat("DONE\n")
