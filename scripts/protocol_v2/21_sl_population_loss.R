# =============================================================================
# scripts/protocol_v2/21_sl_population_loss.R
#
# CAN A POPULATION-AWARE META-LEARNER BEAT THE INDEX ALONE?
#
# Script 20 showed that aligning the meta-learner's loss with ranking (Spearman
# / rank-NNLS) lets the SuperLearner recover the domain index's performance but
# not exceed it. Neither loss knows that districts differ in population by two
# orders of magnitude, while the decision - reach the fifth of districts with
# the most deficiency BURDEN - does. Two population-aware meta-learners:
#
#   wrank    weighted Spearman per learner (discrete); weighted NNLS on ranks
#            (ensemble), weights = population
#   burden   share of prevalence x population captured in the top fifth of the
#            training rows, per learner (discrete) and by Nelder-Mead over
#            softmax weights (ensemble)
#
# Same design as scripts 19/20: ONE SuperLearner fit per fold, every meta-
# learner derived from the identical Z, so differences are the loss and nothing
# else. The production factories (make_method_asl_wrank / _burden) are refitted
# on one fold to confirm they reproduce the derived coefficients.
#
# Evaluated out of fold on: Spearman, top-quintile overlap, BURDEN CAPTURED
# (the NCE metric, using test-row population) and MAE on the prevalence scale,
# so the decision to change the production default can weigh level accuracy
# as well as ranking.
#
#   Rscript scripts/protocol_v2/21_sl_population_loss.R
# -> results/tables/protocol_v2/sl_population_loss_scores.csv
# -> results/tables/protocol_v2/sl_population_loss_selection.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(SuperLearner)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
source("R/area_superlearner.R")

OUTDIR <- "results/tables/protocol_v2"
REPS   <- as.integer(Sys.getenv("SL_REPS", "5"))
TOPFRAC <- 0.20
set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
POP <- readRDS("dashboard/data/admin2_population.rds")
POP$country <- gsub(" ", "", POP$country)  # FIX 2026-09-04: the file spells "Sierra Leone" with a space; without this the join silently dropped the country
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"

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
  pp <- POP[POP$country == cn, c("Admin2", pop_for(on))]; names(pp)[2] <- "pop"
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)],
                  by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1","Admin2","lon","lat")],
               by = c("Admin1", "Admin2")) |>
    left_join(pp, by = "Admin2")
  m <- m[is.finite(m$pop) & m$pop > 0, ]
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  Xr[!is.finite(Xr)] <- 0
  list(n = nrow(m), y = m$y_prev, w = m$n_eff, pop = m$pop, X = Xr)
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

capture <- function(y, pop, score, frac = TOPFRAC) {
  ok <- is.finite(y) & is.finite(pop) & is.finite(score)
  if (sum(ok) < 5) return(NA_real_)
  y <- y[ok]; pop <- pop[ok]; score <- score[ok]; b <- y * pop
  k <- max(1L, round(frac * length(y))); sel <- order(score, decreasing = TRUE)[seq_len(k)]
  sum(b[sel]) / sum(b)
}

META <- c("mse_discrete", "mse_nnls", "rank_discrete", "rank_nnls",
          "wrank_discrete", "wrank_nnls", "burden_discrete", "burden_ens")
BASE <- c("domain_index", "enet", "rf", "mean")
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
    pred <- lapply(c(META, BASE), function(a) rep(NA_real_, cl$n)); names(pred) <- c(META, BASE)

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
      Ytr <- ymod[tr]; pop_tr <- cl$pop[tr]; burden_tr <- cl$y[tr] * pop_tr
      risk <- fit$cvRisk; names(risk) <- colnames(lp)
      cm <- fit$coef; names(cm) <- colnames(lp)
      rk <- .asl_rank_coef(Z, Ytr)
      wr <- .asl_wrank_coef(Z, Ytr, pop_tr)
      bd <- .asl_burden_coef(Z, Ytr, burden_tr, TOPFRAC)
      pk <- list(mse = names(which.min(risk)), rank = colnames(lp)[which.max(rk$rho)],
                 wrank = colnames(lp)[which.max(wr$rho)], burden = colnames(lp)[which.max(bd$cap)])
      pred$mse_discrete[te]    <- lp[, pk$mse];    pred$mse_nnls[te]   <- as.numeric(lp %*% cm)
      pred$rank_discrete[te]   <- lp[, pk$rank];   pred$rank_nnls[te]  <- as.numeric(lp %*% rk$coef_nnls)
      pred$wrank_discrete[te]  <- lp[, pk$wrank];  pred$wrank_nnls[te] <- as.numeric(lp %*% wr$coef_nnls)
      pred$burden_discrete[te] <- lp[, pk$burden]; pred$burden_ens[te] <- as.numeric(lp %*% bd$coef_ens)
      for (nm in BASE) pred[[nm]][te] <- lp[, nm]

      picks[[length(picks) + 1L]] <- data.frame(
        country = cn, outcome = on, rep = r, fold = f,
        pick_mse = pk$mse, pick_rank = pk$rank, pick_wrank = pk$wrank, pick_burden = pk$burden,
        w_mse_idx = cm[["domain_index"]], w_rank_idx = rk$coef_nnls[colnames(lp) == "domain_index"],
        w_wrank_idx = wr$coef_nnls[colnames(lp) == "domain_index"],
        w_burden_idx = bd$coef_ens[colnames(lp) == "domain_index"], stringsAsFactors = FALSE)

      if (is.null(validated)) {
        set.seed(seed)
        f_w <- tryCatch(suppressWarnings(SuperLearner(Y = Ytr, X = Xdf[tr, , drop = FALSE],
                 newX = Xdf[te, , drop = FALSE], family = gaussian(), SL.library = LIB,
                 obsWeights = cl$w[tr], cvControl = list(V = V),
                 method = make_method_asl_wrank(pop_tr))), error = function(e) NULL)
        set.seed(seed)
        f_b <- tryCatch(suppressWarnings(SuperLearner(Y = Ytr, X = Xdf[tr, , drop = FALSE],
                 newX = Xdf[te, , drop = FALSE], family = gaussian(), SL.library = LIB,
                 obsWeights = cl$w[tr], cvControl = list(V = V),
                 method = make_method_asl_burden(burden_tr, TOPFRAC))), error = function(e) NULL)
        validated <- sprintf("wrank factory: %s | burden factory: %s",
          if (is.null(f_w)) "FIT ERROR" else sprintf("max|coef diff| = %.1e", max(abs(unname(f_w$coef) - wr$coef_nnls))),
          if (is.null(f_b)) "FIT ERROR" else sprintf("max|coef diff| = %.1e", max(abs(unname(f_b$coef) - bd$coef_ens))))
      }
    }
    for (nm in c(META, BASE)) {
      s <- score_v2(cl$y, pred[[nm]], cl$w, scale = "prev")
      scores[[length(scores) + 1L]] <- data.frame(
        country = cn, outcome = on, rep = r, arm = nm, n_areas = cl$n,
        spearman = s$spearman, topk = s$topk,
        capture = capture(cl$y, cl$pop, pred[[nm]]),
        mae_pp = 100 * mean(abs(cl$y - .v2_expit(pred[[nm]])), na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
  cat("done", cn, on, "\n")
}

SC <- bind_rows(scores); PK <- bind_rows(picks)
write.csv(SC, file.path(OUTDIR, "sl_population_loss_scores.csv"), row.names = FALSE)
write.csv(PK, file.path(OUTDIR, "sl_population_loss_selection.csv"), row.names = FALSE)

cat("\n===== factories reproduce derived coefficients? =====\n  ", validated, "\n")
cat("\n===== DISCRETE SELECTION OF domain_index (%), n =", nrow(PK), "=====\n")
for (k in c("pick_mse", "pick_rank", "pick_wrank", "pick_burden"))
  cat(sprintf("  %-12s %5.1f%%\n", sub("pick_", "", k), 100 * mean(PK[[k]] == "domain_index")))
cat("\n===== mean ensemble weight on domain_index =====\n")
for (k in c("w_mse_idx", "w_rank_idx", "w_wrank_idx", "w_burden_idx"))
  cat(sprintf("  %-13s %.3f\n", sub("w_|_idx", "", k), mean(PK[[k]], na.rm = TRUE)))

CELL <- SC |> group_by(country, outcome, arm) |>
  summarise(across(c(spearman, topk, capture, mae_pp), ~ median(.x, na.rm = TRUE)), .groups = "drop")
cat("\n===== OUT-OF-FOLD, cell medians over reps (random capture = 0.20) =====\n")
print(as.data.frame(CELL |> group_by(arm) |>
  summarise(cells = dplyr::n(), spearman = round(mean(spearman, na.rm = TRUE), 3),
            topk = round(mean(topk, na.rm = TRUE), 3),
            capture = round(mean(capture, na.rm = TRUE), 3),
            mae_pp = round(mean(mae_pp, na.rm = TRUE), 2), .groups = "drop") |>
  arrange(desc(capture))), row.names = FALSE)

hh <- function(metric, a, b) {
  W <- tidyr::pivot_wider(CELL[, c("country","outcome","arm", metric)], names_from = arm, values_from = all_of(metric))
  d <- W[[a]] - W[[b]]
  sprintf("%-8s %-16s vs %-13s better in %2d of %2d | median %+.3f", metric, a, b,
          sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE))
}
cat("\n===== HEAD-TO-HEAD =====\n")
for (m in c("capture", "spearman")) {
  for (a in c("rank_nnls", "wrank_nnls", "burden_ens", "burden_discrete")) cat(hh(m, a, "mse_discrete"), "\n")
  for (a in c("rank_nnls", "wrank_nnls", "burden_ens", "burden_discrete")) cat(hh(m, a, "domain_index"), "\n")
}
cat("\nMAE (lower is better):\n")
for (a in c("rank_nnls", "wrank_nnls", "burden_ens")) {
  W <- tidyr::pivot_wider(CELL[, c("country","outcome","arm","mae_pp")], names_from = arm, values_from = mae_pp)
  d <- W[[a]] - W[["mse_discrete"]]
  cat(sprintf("  %-12s vs mse_discrete: lower MAE in %2d of %2d | median %+.2f pp\n", a,
              sum(d < 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE)))
}
cat("\nDONE\n")
