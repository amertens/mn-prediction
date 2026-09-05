# =============================================================================
# scripts/protocol_v2/20_sl_rank_loss.R
#
# DOES A RANK-ALIGNED META-LEARNER RESCUE THE SUPERLEARNER?
#
# THE PROBLEM (from 19_sl_with_domain_index.R)
# --------------------------------------------
# With domain_index in the library, the SuperLearner selected it in 1% of
# fits and scored BELOW it (Spearman 0.173 discrete / 0.231 NNLS vs 0.285 for
# the index alone). The constant "mean" learner took the largest NNLS weight.
# The meta-learner minimises cross-validated squared error; at 14-87 noisy
# districts that rewards shrinkage to the mean, and the mean ranks nothing.
# The project's decision is a ranking, so the loss and the decision disagree.
#
# THE TEST
# --------
# The meta-learner only ever sees Z (the cross-validated base-learner
# predictions) and Y. So ONE SuperLearner fit per fold yields Z, and four
# meta-learners are derived from the identical Z:
#
#   mse_discrete   argmin CV squared error            (production default)
#   mse_nnls       NNLS on Z                          (SuperLearner default)
#   rank_discrete  argmax CV Spearman                 (method.asl_rank pick)
#   rank_nnls      NNLS on ranks of Z vs ranks of Y   (method.asl_rank ensemble)
#
# Same base fits, same folds, same seeds: any difference is the loss. The
# production method object (method.asl_rank in R/area_superlearner.R) is then
# refitted on one fold to confirm it reproduces the derived coefficients.
#
# PROTOCOL: identical to script 19. Admin-2, prevalence on the logit scale,
# 5-fold district folds replicated SL_REPS times, survey weights on the base
# learners, meta-CV nested inside each outer fold.
#
#   Rscript scripts/protocol_v2/20_sl_rank_loss.R
# -> results/tables/protocol_v2/sl_rank_loss_scores.csv
# -> results/tables/protocol_v2/sl_rank_loss_selection.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(SuperLearner)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
source("R/area_superlearner.R")      # SL.asl_* wrappers, .asl_rank_coef, method.asl_rank

OUTDIR <- "results/tables/protocol_v2"
REPS   <- as.integer(Sys.getenv("SL_REPS", "5"))
set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)],
                  by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1","Admin2","lon","lat")],
               by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  Xr[!is.finite(Xr)] <- 0
  list(n = nrow(m), y = m$y_prev, w = m$n_eff, X = Xr)
}

DOMAIN_OF_DF <- NULL
SL.domain_index <- function(Y, X, newX, family, obsWeights, id, ...) {
  Xall <- rbind(as.matrix(X), as.matrix(newX)); ntr <- nrow(X)
  D <- domain_representation_v2(Xall, DOMAIN_OF_DF, sign_rows = seq_len(ntr))
  p <- arm_domain_index_v2(seq_len(ntr), ntr + seq_len(nrow(newX)),
                           c(Y, rep(NA_real_, nrow(newX))), Xall, D, NULL)
  fit <- list(); class(fit) <- "SL.domain_index"
  list(pred = p, fit = fit)
}
predict.SL.domain_index <- function(object, newdata, ...)
  stop("SL.domain_index predicts via newX at fit time")
LIB <- c("SL.mean", "SL.asl_enet", "SL.asl_ranger", "SL.domain_index")
SHORT <- c(SL.mean = "mean", SL.asl_enet = "enet", SL.asl_ranger = "rf",
           SL.domain_index = "domain_index")

scores <- list(); picks <- list(); validated <- NULL
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)

for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  ymod <- .v2_logit(cl$y)
  Xdf <- as.data.frame(cl$X)
  orig <- colnames(cl$X); names(Xdf) <- make.names(orig, unique = TRUE)
  DOMAIN_OF_DF <<- stats::setNames(domain_of[orig], names(Xdf))

  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    arms <- c("mse_discrete", "mse_nnls", "rank_discrete", "rank_nnls",
              "domain_index", "enet", "rf", "mean")
    pred <- lapply(arms, function(a) rep(NA_real_, cl$n)); names(pred) <- arms

    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12) next
      V <- max(3L, min(5L, floor(length(tr) / 5)))
      seed <- 1000 * i + 10 * r + f
      set.seed(seed)
      fit <- tryCatch(suppressWarnings(SuperLearner(
        Y = ymod[tr], X = Xdf[tr, , drop = FALSE], newX = Xdf[te, , drop = FALSE],
        family = gaussian(), SL.library = LIB, obsWeights = cl$w[tr],
        cvControl = list(V = V), method = "method.NNLS")),
        error = function(e) { cat("  SL error:", conditionMessage(e), "\n"); NULL })
      if (is.null(fit)) next

      lp <- fit$library.predict; colnames(lp) <- SHORT[sub("_All$", "", fit$libraryNames)]
      Z  <- fit$Z; colnames(Z) <- colnames(lp)
      risk <- fit$cvRisk; names(risk) <- colnames(lp)
      coef_mse <- fit$coef; names(coef_mse) <- colnames(lp)
      rk <- .asl_rank_coef(Z, ymod[tr]); names(rk$rho) <- colnames(lp)
      coef_rank <- rk$coef_nnls; names(coef_rank) <- colnames(lp)

      pick_mse  <- names(which.min(risk))
      pick_rank <- names(which.max(rk$rho))
      pred$mse_discrete[te]  <- lp[, pick_mse]
      pred$mse_nnls[te]      <- as.numeric(lp %*% coef_mse)
      pred$rank_discrete[te] <- lp[, pick_rank]
      pred$rank_nnls[te]     <- as.numeric(lp %*% coef_rank)
      for (nm in c("domain_index", "enet", "rf", "mean")) pred[[nm]][te] <- lp[, nm]

      picks[[length(picks) + 1L]] <- data.frame(
        country = cn, outcome = on, rep = r, fold = f,
        pick_mse = pick_mse, pick_rank = pick_rank,
        w_mse_index = coef_mse[["domain_index"]], w_rank_index = coef_rank[["domain_index"]],
        w_mse_mean = coef_mse[["mean"]], w_rank_mean = coef_rank[["mean"]],
        stringsAsFactors = FALSE)

      # one-time validation: the production method object must reproduce the
      # derived coefficients on the identical folds (same seed -> same validRows)
      if (is.null(validated)) {
        set.seed(seed)
        fit2 <- tryCatch(suppressWarnings(SuperLearner(
          Y = ymod[tr], X = Xdf[tr, , drop = FALSE], newX = Xdf[te, , drop = FALSE],
          family = gaussian(), SL.library = LIB, obsWeights = cl$w[tr],
          cvControl = list(V = V), method = method.asl_rank)),
          error = function(e) { cat("  method.asl_rank error:", conditionMessage(e), "\n"); NULL })
        validated <- if (is.null(fit2)) "FAILED (fit error)" else {
          c2 <- fit2$coef; names(c2) <- colnames(lp)
          sprintf("max |coef diff| = %.2e; pick agrees = %s",
                  max(abs(c2 - coef_rank)),
                  identical(names(which.min(fit2$cvRisk)) |> sub("_All$", "", x = _) |> (\(z) SHORT[[z]])(),
                            pick_rank))
        }
      }
    }
    for (nm in arms) {
      s <- score_v2(cl$y, pred[[nm]], cl$w, scale = "prev")
      scores[[length(scores) + 1L]] <- data.frame(
        country = cn, outcome = on, rep = r, arm = nm, n_areas = cl$n,
        spearman = s$spearman, topk = s$topk, stringsAsFactors = FALSE)
    }
  }
  cat("done", cn, on, "\n")
}

SC <- bind_rows(scores); PK <- bind_rows(picks)
write.csv(SC, file.path(OUTDIR, "sl_rank_loss_scores.csv"), row.names = FALSE)
write.csv(PK, file.path(OUTDIR, "sl_rank_loss_selection.csv"), row.names = FALSE)

cat("\n===== production method.asl_rank reproduces derived coefficients? =====\n  ", validated, "\n")

cat("\n===== DISCRETE SELECTION FREQUENCY (%), n =", nrow(PK), "fits =====\n")
cat("-- MSE loss --\n");  print(round(100 * prop.table(table(PK$pick_mse)), 1))
cat("-- RANK loss --\n"); print(round(100 * prop.table(table(PK$pick_rank)), 1))
cat("\n===== mean ensemble weight on domain_index / on the constant =====\n")
cat(sprintf("  MSE-NNLS : index %.3f | mean %.3f\n  RANK-NNLS: index %.3f | mean %.3f\n",
            mean(PK$w_mse_index), mean(PK$w_mse_mean), mean(PK$w_rank_index), mean(PK$w_rank_mean)))

cat("\n===== OUT-OF-FOLD SPEARMAN, prevalence (cell medians over reps) =====\n")
CELL <- SC |> group_by(country, outcome, arm) |>
  summarise(spearman = median(spearman, na.rm = TRUE), .groups = "drop")
print(as.data.frame(CELL |> group_by(arm) |>
  summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
            median_rho = round(median(spearman, na.rm = TRUE), 3),
            cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_rho))), row.names = FALSE)

W <- tidyr::pivot_wider(CELL, names_from = arm, values_from = spearman)
hh <- function(a, b) sprintf("%-14s vs %-14s: better in %2d of %2d cells (median diff %+.3f)",
  a, b, sum(W[[a]] > W[[b]], na.rm = TRUE), sum(is.finite(W[[a]] - W[[b]])),
  median(W[[a]] - W[[b]], na.rm = TRUE))
cat("\n===== HEAD-TO-HEAD =====\n")
cat(hh("rank_discrete", "mse_discrete"), "\n"); cat(hh("rank_nnls", "mse_nnls"), "\n")
cat(hh("rank_discrete", "domain_index"), "\n"); cat(hh("rank_nnls", "domain_index"), "\n")
cat("\nDONE\n")
