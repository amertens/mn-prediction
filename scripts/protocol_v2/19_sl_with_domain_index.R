# =============================================================================
# scripts/protocol_v2/19_sl_with_domain_index.R
#
# CAN THE DOMAIN INDEX BE A SUPERLEARNER LIBRARY MEMBER, AND WHICH PACKAGE?
#
# THE QUESTION
# ------------
# domain_index beat every tuned learner under the corrected protocol. If it is
# added to a SuperLearner library, does the SuperLearner SELECT it? That turns
# "the simple index beats the ensemble" into a result the ensemble produces
# itself, which is a much stronger sentence. The number to report is the
# SELECTION FREQUENCY, not the ensemble's score.
#
# WHY TWO PACKAGES
# ----------------
# The project's area-level SL (R/benchmark_models.R, fit_predict_sl_prescreened)
# uses mlr3superlearner. Reading its source (0.1.2) shows three things that
# matter for THIS project and are not documented prominently:
#
#   1. The library is a hard whitelist (available_learners_regr: 17 keys).
#      lookup() filters by name, so a custom learner cannot be added without
#      modifying the package. domain_index CANNOT go in.
#   2. There is no observation-weight argument. Area-level rows carry very
#      different effective n (6 to 500+); the production SL fits them equally.
#   3. make_mlr3_resampling() returns plain rsmp("cv") for regression BEFORE it
#      checks the group column, so `group=` is silently ignored for continuous
#      outcomes. District/cluster blocking never happens in the regression case.
#      And set_folds() picks leave-one-out for n < 30.
#
# The classic SuperLearner package has obsWeights, cvControl(validRows=) for
# custom folds, id= for clusters, and a custom wrapper is a plain function.
# sl3 + origami can do all of this too, but is heavier and the project already
# retired its sl3 path. So this script runs the SAME outer folds through both
# mlr3superlearner (native library only, no weights) and SuperLearner (native
# library + domain_index, n_eff weights), and reports side by side.
#
# PROTOCOL. Admin-2, prevalence target on the logit scale, 5-fold district
# folds replicated SL_REPS times (same make_folds_v2 as the benchmarks). The
# SuperLearner's own CV for the meta-learner is nested INSIDE each outer fold,
# so nothing the held-out districts contain reaches learner selection.
#
#   Rscript scripts/protocol_v2/19_sl_with_domain_index.R
# -> results/tables/protocol_v2/sl_domain_index_scores.csv
# -> results/tables/protocol_v2/sl_domain_index_selection.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(SuperLearner)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
REPS   <- as.integer(Sys.getenv("SL_REPS", "5"))
set.seed(20260903L)
HAS_MLR3 <- requireNamespace("mlr3superlearner", quietly = TRUE)

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

# ── SuperLearner wrappers ────────────────────────────────────────────────────
# DOMAIN_OF_DF is set per cell: the wrapper must know each column's domain, and
# SuperLearner passes X as a data.frame whose names must be syntactic.
DOMAIN_OF_DF <- NULL
SL.domain_index <- function(Y, X, newX, family, obsWeights, ...) {
  Xall <- rbind(as.matrix(X), as.matrix(newX)); ntr <- nrow(X)
  D <- domain_representation_v2(Xall, DOMAIN_OF_DF, sign_rows = seq_len(ntr))
  p <- arm_domain_index_v2(seq_len(ntr), ntr + seq_len(nrow(newX)),
                           c(Y, rep(NA_real_, nrow(newX))), Xall, D, NULL)
  fit <- list(); class(fit) <- "SL.domain_index"
  list(pred = p, fit = fit)
}
predict.SL.domain_index <- function(object, newdata, ...)
  stop("SL.domain_index predicts via newX at fit time")
SL.enet <- function(...) SL.glmnet(..., alpha = 0.5, nfolds = 3, useMin = TRUE)
SL.rf   <- function(...) SL.ranger(..., num.trees = 250, min.node.size = 5)
LIB_SL <- c("SL.mean", "SL.enet", "SL.rf", "SL.domain_index")

# mlr3superlearner: the same three native learners. domain_index cannot be
# added (whitelist), and there is no weights argument.
LIB_MLR3 <- list(list("mean", id = "mean"),
                 list("glmnet", alpha = 0.5, id = "enet"),
                 list("ranger", num.trees = 250, min.node.size = 5, id = "rf"))

scores <- list(); picks <- list(); mlr3_custom_msg <- NA_character_
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)

for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  ymod <- .v2_logit(cl$y)
  Xdf <- as.data.frame(cl$X)
  orig <- colnames(cl$X); names(Xdf) <- make.names(orig, unique = TRUE)
  DOMAIN_OF_DF <<- stats::setNames(domain_of[orig], names(Xdf))

  # one-time demonstration that mlr3superlearner rejects a custom learner
  if (HAS_MLR3 && is.na(mlr3_custom_msg)) {
    mlr3_custom_msg <- tryCatch({
      suppressMessages(mlr3superlearner::mlr3superlearner(
        data = data.frame(Y = ymod, Xdf), target = "Y",
        library = c(LIB_MLR3, list(list("domain_index", id = "domain_index"))),
        outcome_type = "continuous", folds = 3))
      "ACCEPTED (unexpected)"
    }, error = function(e) paste("REJECTED:", conditionMessage(e)))
  }

  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    pred <- list(sl_discrete = rep(NA_real_, cl$n), sl_nnls = rep(NA_real_, cl$n),
                 mlr3_discrete = rep(NA_real_, cl$n),
                 domain_index = rep(NA_real_, cl$n), enet = rep(NA_real_, cl$n),
                 rf = rep(NA_real_, cl$n), mean = rep(NA_real_, cl$n))
    t_sl <- 0; t_m3 <- 0
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12) next
      V <- max(3L, min(5L, floor(length(tr) / 5)))
      set.seed(1000 * i + 10 * r + f)

      # --- classic SuperLearner: weights, custom learner, nested CV ----------
      t0 <- Sys.time()
      fit <- tryCatch(suppressWarnings(SuperLearner(
        Y = ymod[tr], X = Xdf[tr, , drop = FALSE], newX = Xdf[te, , drop = FALSE],
        family = gaussian(), SL.library = LIB_SL, obsWeights = cl$w[tr],
        cvControl = list(V = V), method = "method.NNLS")),
        error = function(e) NULL)
      t_sl <- t_sl + as.numeric(Sys.time() - t0, units = "secs")
      if (!is.null(fit)) {
        lp <- fit$library.predict
        colnames(lp) <- sub("_All$", "", colnames(lp))
        pick <- names(which.min(fit$cvRisk)); pick <- sub("_All$", "", pick)
        pred$sl_discrete[te] <- lp[, pick]
        pred$sl_nnls[te]     <- as.numeric(fit$SL.predict)
        for (nm in c("domain_index", "enet", "rf", "mean"))
          pred[[nm]][te] <- lp[, paste0("SL.", nm)]
        co <- fit$coef; names(co) <- sub("_All$", "", names(co))
        picks[[length(picks) + 1L]] <- data.frame(
          country = cn, outcome = on, rep = r, fold = f, package = "SuperLearner",
          discrete_pick = sub("^SL\\.", "", pick),
          w_domain_index = unname(co["SL.domain_index"]), w_enet = unname(co["SL.enet"]),
          w_rf = unname(co["SL.rf"]), w_mean = unname(co["SL.mean"]),
          stringsAsFactors = FALSE)
      }

      # --- mlr3superlearner: native library only, no weights -----------------
      if (HAS_MLR3) {
        t0 <- Sys.time()
        m3 <- tryCatch(suppressMessages(suppressWarnings(
          mlr3superlearner::mlr3superlearner(
            data = data.frame(Y = ymod[tr], Xdf[tr, , drop = FALSE]), target = "Y",
            library = LIB_MLR3, outcome_type = "continuous", folds = V,
            discrete = TRUE))), error = function(e) NULL)
        t_m3 <- t_m3 + as.numeric(Sys.time() - t0, units = "secs")
        if (!is.null(m3)) {
          pred$mlr3_discrete[te] <- tryCatch(
            as.numeric(stats::predict(m3, Xdf[te, , drop = FALSE])),
            error = function(e) rep(NA_real_, length(te)))
          picks[[length(picks) + 1L]] <- data.frame(
            country = cn, outcome = on, rep = r, fold = f, package = "mlr3superlearner",
            discrete_pick = tryCatch(m3$learners[[1]]$id, error = function(e) NA_character_),
            w_domain_index = NA_real_, w_enet = NA_real_, w_rf = NA_real_, w_mean = NA_real_,
            stringsAsFactors = FALSE)
        }
      }
    }
    for (nm in names(pred)) {
      s <- score_v2(cl$y, pred[[nm]], cl$w, scale = "prev")
      scores[[length(scores) + 1L]] <- data.frame(
        country = cn, outcome = on, rep = r, arm = nm, n_areas = cl$n,
        spearman = s$spearman, topk = s$topk, secs_sl = t_sl, secs_mlr3 = t_m3,
        stringsAsFactors = FALSE)
    }
  }
  cat("done", cn, on, "\n")
}

SC <- bind_rows(scores); PK <- bind_rows(picks)
write.csv(SC, file.path(OUTDIR, "sl_domain_index_scores.csv"), row.names = FALSE)
write.csv(PK, file.path(OUTDIR, "sl_domain_index_selection.csv"), row.names = FALSE)

cat("\n===== mlr3superlearner: custom learner attempt =====\n ", mlr3_custom_msg, "\n")

cat("\n===== DISCRETE SELECTION FREQUENCY (cell x rep x fold) =====\n")
for (pk in unique(PK$package)) {
  d <- PK[PK$package == pk, ]
  cat("--", pk, ": n =", nrow(d), "--\n")
  print(round(100 * prop.table(table(d$discrete_pick)), 1))
}
if (any(PK$package == "SuperLearner")) {
  cat("\n===== NNLS ensemble weights, SuperLearner (mean over fits) =====\n")
  print(round(colMeans(PK[PK$package == "SuperLearner",
                          c("w_domain_index","w_enet","w_rf","w_mean")], na.rm = TRUE), 3))
}

cat("\n===== OUT-OF-FOLD SPEARMAN, prevalence (cell medians over reps) =====\n")
CELL <- SC |> group_by(country, outcome, arm) |>
  summarise(spearman = median(spearman, na.rm = TRUE), .groups = "drop")
print(as.data.frame(CELL |> group_by(arm) |>
  summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
            median_rho = round(median(spearman, na.rm = TRUE), 3),
            cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_rho))), row.names = FALSE)

W <- tidyr::pivot_wider(CELL, names_from = arm, values_from = spearman)
hh <- function(a, b) sprintf("%-14s vs %-13s: better in %d of %d cells (median diff %+.3f)",
  a, b, sum(W[[a]] > W[[b]], na.rm = TRUE), sum(is.finite(W[[a]] - W[[b]])),
  median(W[[a]] - W[[b]], na.rm = TRUE))
cat("\n===== HEAD-TO-HEAD =====\n")
cat(hh("sl_discrete", "domain_index"), "\n")
cat(hh("sl_nnls", "domain_index"), "\n")
cat(hh("domain_index", "enet"), "\n")
if (HAS_MLR3) cat(hh("sl_discrete", "mlr3_discrete"), "\n")

cat(sprintf("\nruntime per cell-rep: SuperLearner %.1fs | mlr3superlearner %.1fs\n",
            mean(SC$secs_sl[SC$arm == "mean"]), mean(SC$secs_mlr3[SC$arm == "mean"])))
cat("\nDONE\n")
