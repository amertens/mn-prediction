# =============================================================================
# explore/R/harness.R
#
# A thin, generalised re-implementation of the protocol-v2 scoring loop, so
# that any arm and any predictor matrix can be scored on the SAME folds and the
# SAME metrics as the numbers already on the record.
#
# It deliberately reuses R/protocol_v2.R rather than reimplementing it: the cell
# construction, rank-normalisation (fix 3), domain representation (fix 4),
# fold seeds (fix 1) and scoring are the project's, not this sandbox's. What is
# new here is only that the arm and the feature matrix are arguments.
#
# Arm signature is protocol-v2's, unchanged:
#     function(tr, te, y, X, D, aux) -> numeric(length(te))
# so every existing arm drops in without a wrapper.
#
# READ-ONLY on the main project. Writes nothing outside explore/.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

EXP_ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
EXP_OUT  <- file.path(EXP_ROOT, "explore", "out")

source(file.path(EXP_ROOT, "R", "protocol_v2.R"))

# ── data ────────────────────────────────────────────────────────────────────

#' Load the 24-cell target table, the shared predictor store and centroids
#'
#' Applies the project's own leakage rule (drop_near_outcome_v2) so this
#' sandbox cannot accidentally be more permissive than the pipeline.
exp_load <- function() {
  # TP-01: the headline benchmark run (benchmarks_v2_cells.csv) uses the
  # open + public-survey tiers throughout, so this sandbox defaults to the same
  # set. Without it, every number here is quietly computed on a larger
  # predictor set than the record it is being compared against.
  if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS")))
    Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")

  TG <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/targets_v2.csv"),
                 stringsAsFactors = FALSE)
  S  <- read.csv(file.path(EXP_ROOT, "data/covariates/harmonized/predictors_admin2_shared.csv"),
                 check.names = FALSE)
  MD <- read.csv(file.path(EXP_ROOT, "data/covariates/harmonized/predictors_admin2_shared_metadata.csv"),
                 stringsAsFactors = FALSE)
  PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
  domain_of <- stats::setNames(MD$domain, MD$column)

  BND <- readRDS(file.path(EXP_ROOT, "dashboard/data/admin2_boundaries.rds"))
  countries <- c(gambia = "Gambia", ghana = "Ghana",
                 malawi = "Malawi", sierraleone = "SierraLeone")
  CENT <- do.call(rbind, lapply(names(countries), function(lc) {
    b  <- BND[[lc]]
    xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
    data.frame(country = countries[[lc]],
               Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
               Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
               lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
  }))

  list(TG = TG, S = S, MD = MD, PREDS = PREDS, domain_of = domain_of,
       CENT = CENT, countries = unname(countries))
}

#' All 24 country x outcome cells present in the target table
exp_cell_index <- function(E) {
  E$TG |> dplyr::distinct(country, outcome) |>
    dplyr::filter(country %in% E$countries) |> as.data.frame()
}

#' Build one cell, exactly as scripts/protocol_v2/02 does
#'
#' @param cols optional subset of predictor columns (default: the full set).
#'   A probe that wants its own feature block passes `extra` instead.
#' @param extra optional data.frame with country/Admin1/Admin2 plus new
#'   columns, joined on the pair key and appended to the predictor matrix.
#' @param extra_domain optional domain label for `extra`'s columns, so they
#'   take part in the domain representation.
#' @param prep TRUE applies the project's rank-normalise-within-country and
#'   median-impute (fix 3); FALSE returns the raw columns.
exp_cell <- function(E, cn, on, target, cols = NULL, extra = NULL,
                     extra_domain = NULL, prep = TRUE, min_cols = 20L) {
  preds <- if (is.null(cols)) E$PREDS else intersect(cols, E$PREDS)
  domain_of <- E$domain_of

  t <- E$TG[E$TG$country == cn & E$TG$outcome == on, ]
  ycol     <- if (target == "prev") "y_prev" else "y_level"
  ncol_eff <- if (target == "prev") "n_eff"  else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[ncol_eff]]), ]
  if (nrow(t) < 12) return(NULL)

  sc <- E$S[E$S$country == cn, c("Admin1", "Admin2", preds)]
  m  <- dplyr::inner_join(t, sc, by = c("Admin1", "Admin2"))

  if (!is.null(extra)) {
    ex <- extra[extra$country == cn, setdiff(names(extra), "country"), drop = FALSE]
    newc <- setdiff(names(ex), c("Admin1", "Admin2"))
    m <- dplyr::left_join(m, ex, by = c("Admin1", "Admin2"))
    preds <- c(preds, newc)
    if (!is.null(extra_domain))
      domain_of <- c(domain_of, stats::setNames(rep(extra_domain, length(newc)), newc))
  }

  m <- dplyr::inner_join(m, E$CENT[E$CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")],
                         by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)

  X0 <- as.matrix(m[, preds, drop = FALSE])
  Xr <- if (prep) prep_predictors_v2(X0) else X0
  if (ncol(Xr) < min_cols) return(NULL)
  D <- domain_representation_v2(Xr, domain_of)

  y_nat <- m[[ycol]]
  y_mod <- if (target == "prev") .v2_logit(y_nat) else y_nat

  list(country = cn, outcome = on, target = target,
       y_nat = y_nat, y_mod = y_mod, X = Xr, D = D,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1,
                  y_nat = y_nat, target = target, w = m[[ncol_eff]]),
       w = m[[ncol_eff]], Admin1 = m$Admin1, Admin2 = m$Admin2, n = nrow(m))
}

#' Build every cell for one target, as a named list
exp_all_cells <- function(E, target, outcomes = NULL, countries = NULL, ...) {
  ix <- exp_cell_index(E)
  if (!is.null(outcomes))  ix <- ix[ix$outcome %in% outcomes, ]
  if (!is.null(countries)) ix <- ix[ix$country %in% countries, ]
  out <- list()
  for (i in seq_len(nrow(ix))) {
    cc <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], target, ...),
                   error = function(e) { message("  cell failed: ", ix$country[i],
                                                 " ", ix$outcome[i], ": ", conditionMessage(e)); NULL })
    if (!is.null(cc)) out[[paste(ix$country[i], ix$outcome[i], sep = "|")]] <- cc
  }
  out
}

# ── estimand A: in-fill within country ──────────────────────────────────────

#' Replicated 5-fold over districts; out-of-fold predictions; weighted scoring
#'
#' Fold seeds come from make_folds_v2(), so a probe's folds are identical to
#' the ones behind the numbers on the record.
exp_infill <- function(cell, arms, reps = 10L, k = 5L, quiet = TRUE) {
  stopifnot(is.list(arms), length(arms) > 0)
  rows <- list()
  for (r in seq_len(reps)) {
    folds <- make_folds_v2("kfold_district", cell$n, k = k, rep_id = r)
    for (a in names(arms)) {
      fn <- arms[[a]]
      pred_mod <- rep(NA_real_, cell$n)
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12 || !length(te)) next
        p <- tryCatch(fn(tr, te, cell$y_mod, cell$X, cell$D, cell$aux),
                      error = function(e) { if (!quiet) message("   arm ", a, ": ",
                                                               conditionMessage(e)); rep(NA_real_, length(te)) })
        if (length(p) == length(te)) pred_mod[te] <- p
      }
      pred <- if (cell$target == "prev") .v2_expit(pred_mod) else pred_mod
      if (a == "null_train_mean") {        # null defined on the natural scale
        pred <- rep(NA_real_, cell$n)
        for (f in unique(folds)) {
          te <- which(folds == f); tr <- which(folds != f)
          if (length(tr) < 12 || !length(te)) next
          pred[te] <- mean(cell$y_nat[tr])
        }
      }
      s <- score_v2(cell$y_nat, pred, cell$w,
                    scale = if (cell$target == "prev") "prev" else "level")
      rows[[paste(a, r)]] <- cbind(
        data.frame(country = cell$country, outcome = cell$outcome,
                   target = cell$target, estimand = "infill", arm = a,
                   rep = r, n_areas = cell$n, stringsAsFactors = FALSE), s)
    }
  }
  dplyr::bind_rows(rows)
}

#' Region extrapolation (estimand B): exhaustive leave-one-region-out
exp_region <- function(cell, arms, quiet = TRUE) {
  if (dplyr::n_distinct(cell$Admin1) < 3) return(NULL)
  folds <- make_folds_v2("loro", cell$n, blocks = cell$Admin1)
  rows <- list()
  for (a in names(arms)) {
    fn <- arms[[a]]
    pred_mod <- rep(NA_real_, cell$n)
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12 || !length(te)) next
      p <- tryCatch(fn(tr, te, cell$y_mod, cell$X, cell$D, cell$aux),
                    error = function(e) { if (!quiet) message("   arm ", a, ": ",
                                                             conditionMessage(e)); rep(NA_real_, length(te)) })
      if (length(p) == length(te)) pred_mod[te] <- p
    }
    pred <- if (cell$target == "prev") .v2_expit(pred_mod) else pred_mod
    s <- score_v2(cell$y_nat, pred, cell$w,
                  scale = if (cell$target == "prev") "prev" else "level")
    rows[[a]] <- cbind(data.frame(country = cell$country, outcome = cell$outcome,
                                  target = cell$target, estimand = "region", arm = a,
                                  rep = 1L, n_areas = cell$n, stringsAsFactors = FALSE), s)
  }
  dplyr::bind_rows(rows)
}

# ── estimand C: transport to an unseen country ──────────────────────────────

#' Leave-one-country-out, scored WITHIN each held-out country
#'
#' Outcomes are within-country standardised before pooling (the outcome-side
#' twin of fix 3), so only rank metrics are meaningful and only those are kept.
#'
#' @param cells named list from exp_all_cells(), for ONE outcome
#' @param domain_of column -> domain map; when supplied, the domain axes are
#'   rebuilt INSIDE each fold on the pooled matrix with the orientation learned
#'   from the training countries only. This matters: building domain scores per
#'   country lets each country learn its own PC1 sign, so a domain score means
#'   the opposite thing in two countries and transport is destroyed. 02b does
#'   it this way and so must anything compared against 02b's numbers.
#'   Pass NULL to hand the arm the pooled feature matrix as both X and D
#'   (for probes whose feature block has no domain structure).
exp_loco <- function(cells, arms, domain_of = NULL, min_common = 20L, quiet = TRUE) {
  if (length(cells) < 3) return(NULL)
  common <- Reduce(intersect, lapply(cells, function(z) colnames(z$X)))
  if (length(common) < min_common) return(NULL)

  Y    <- unlist(lapply(cells, function(z) as.numeric(scale(z$y_mod))))
  Xm   <- do.call(rbind, lapply(cells, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(vapply(cells, function(z) z$country, ""), vapply(cells, function(z) z$n, 0L))
  ynat <- unlist(lapply(cells, function(z) z$y_nat))
  wv   <- unlist(lapply(cells, function(z) z$w))
  target <- cells[[1]]$target
  aux <- list(lon = unlist(lapply(cells, function(z) z$aux$lon)),
              lat = unlist(lapply(cells, function(z) z$aux$lat)),
              Admin1 = paste(ctry, unlist(lapply(cells, function(z) z$Admin1))),
              y_nat = Y, target = "level", w = wv, country = ctry)

  folds <- as.integer(factor(ctry))
  rows <- list()
  for (a in names(arms)) {
    fn <- arms[[a]]
    pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 20) next
      Dm <- if (is.null(domain_of)) Xm else
        domain_representation_v2(Xm, domain_of, sign_rows = tr)
      p <- tryCatch(fn(tr, te, Y, Xm, Dm, aux),
                    error = function(e) { if (!quiet) message("   arm ", a, ": ",
                                                              conditionMessage(e)); rep(NA_real_, length(te)) })
      if (length(p) == length(te)) pred[te] <- p
    }
    for (cn in unique(ctry)) {
      k <- which(ctry == cn)
      s <- score_v2(ynat[k], pred[k], wv[k],
                    scale = if (target == "prev") "prev" else "level")
      s$mae <- NA_real_; s$wmae <- NA_real_; s$bias <- NA_real_; s$rmse_sd <- NA_real_
      rows[[paste(a, cn)]] <- cbind(
        data.frame(country = cn, outcome = cells[[1]]$outcome, target = target,
                   estimand = "country", arm = a, rep = 1L, n_areas = length(k),
                   n_common = length(common), stringsAsFactors = FALSE), s)
    }
  }
  dplyr::bind_rows(rows)
}

# ── comparators every probe must report ─────────────────────────────────────

#' The baseline arms, so a probe's number is never read without its comparator
exp_baseline_arms <- function(which = c("null_train_mean", "spatial", "domain_index")) {
  ARMS_V2[which]
}

#' Summarise a raw score table to one row per cell x arm (mean over reps)
exp_summarise <- function(raw) {
  raw |>
    dplyr::group_by(country, outcome, target, estimand, arm) |>
    dplyr::summarise(n_areas = dplyr::first(n_areas),
                     reps = dplyr::n(),
                     spearman = mean(spearman, na.rm = TRUE),
                     pearson  = mean(pearson,  na.rm = TRUE),
                     mae      = mean(mae,      na.rm = TRUE),
                     wmae     = mean(wmae,     na.rm = TRUE),
                     rmse_sd  = mean(rmse_sd,  na.rm = TRUE),
                     topk     = mean(topk,     na.rm = TRUE),
                     .groups = "drop") |>
    as.data.frame()
}

#' Write a probe's tables to explore/out/ and echo the headline
exp_write <- function(df, name) {
  dir.create(EXP_OUT, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(EXP_OUT, paste0(name, ".csv"))
  utils::write.csv(df, f, row.names = FALSE)
  message("-> ", f, "  (", nrow(df), " rows)")
  invisible(f)
}
