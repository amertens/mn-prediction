# =============================================================================
# scripts/protocol_v2/17_model_augmented_survey.R   [GAP 4]
#
# DOES THE MODEL LET A COUNTRY FIELD A SMALLER SURVEY?
#
# THE CLAIM WITHOUT EVIDENCE
# --------------------------
# The NCE says the survey sample needed for subnational estimates "may be
# substantially smaller when supplemented by our modeling approach". Our design
# curve (ws5_anchoring_budget.R) measures only the SURVEY-ONLY baseline: give
# every district its region's jackknifed survey mean. The model has never been
# put on top of it, so the sentence currently asserts something we have not
# tested. This script tests it.
#
# THE TWO ARMS, SCORED ON IDENTICAL SUB-SAMPLES
#   survey_only  every district gets its region's mean, computed from the
#                sub-sampled survey. Flat within a region by construction.
#   model_aug    the same regional mean supplies the LEVEL, and the model
#                supplies the WITHIN-REGION PATTERN:
#                   pred_d = region_mean_sub + (mhat_d - mean(mhat in region))
#                This is the "survey supplies the level, model supplies the
#                pattern" design, and it is the only way the model can help a
#                survey-anchored estimate: a monotone shift cannot change a
#                ranking, so the model has to move districts RELATIVE to their
#                regional mean or it adds nothing.
#
# mhat comes from an in-fill 5-fold protocol-v2 fit in which a district never
# sees its own outcome, so the model arm gets no information the survey arm is
# denied beyond the covariates themselves.
#
# HOW THE SMALLER SURVEY IS SIMULATED, AND WHY IT IS A SIMULATION
# ---------------------------------------------------------------
# Fielding a fraction f of clusters reduces a district's EFFECTIVE sample size
# to about f x n_eff. We therefore redraw each district's observed prevalence as
#     y_sub ~ Binomial(round(f x n_eff), y_full) / round(f x n_eff)
# n_eff is the Kish effective n already divided by the MEASURED design effect,
# so the cluster structure is carried rather than assumed. This is a parametric
# resample, not a re-subsampling of clusters from raw data as ws5 does, and it
# is labelled as such. VALIDATION: the survey_only arm here should reproduce the
# empirical ws5 curve; the script prints both so the reader can check, and if
# they disagree the simulation is wrong and the model comparison is void.
#
# Scored against the FULL survey's district estimates, the same stand-in for
# truth ws5 uses.
#
#   Rscript scripts/protocol_v2/17_model_augmented_survey.R
# -> results/tables/protocol_v2/model_augmented_survey.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR    <- "results/tables/protocol_v2"
FRACTIONS <- c(0.15, 0.25, 0.40, 0.60, 0.80, 1.00)
REPS      <- as.integer(Sys.getenv("G4_REPS", "40"))
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
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1","Admin2","lon","lat")],
               by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 2) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  list(n = nrow(m), y = m$y_prev, n_eff = m$n_eff, Admin1 = m$Admin1,
       X = Xr, D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1,
                  y_nat = m$y_prev))
}

#' honest in-fill model prediction: 5-fold over districts, averaged over
#' replicates so a single fold draw does not decide the answer.
#'
#' Parameterised by arm because the first run of this script used
#' spatial_plus_domain while the in-fill benchmark's winner is domain_index,
#' and the two disagreed about whether the model helps. Running both here
#' removes the arm as an explanation for that disagreement.
model_pred <- function(cl, arm, nrep = 5L) {
  ym <- .v2_logit(cl$y)
  acc <- matrix(NA_real_, nrow = cl$n, ncol = nrep)
  for (r in seq_len(nrep)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12) next
      p <- tryCatch(ARMS_V2[[arm]](tr, te, ym, cl$X, cl$D, cl$aux),
                    error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) acc[te, r] <- p
    }
  }
  rowMeans(acc, na.rm = TRUE)          # on the logit scale
}

rows <- list()
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)

for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  mh <- list(sp  = model_pred(cl, "spatial_plus_domain"),
             idx = model_pred(cl, "domain_index"))
  if (all(!is.finite(mh$sp)) && all(!is.finite(mh$idx))) next
  reg <- factor(cl$Admin1)

  for (f in FRACTIONS) {
    for (r in seq_len(REPS)) {
      nn <- pmax(round(f * cl$n_eff), 1)
      # INCREMENTAL noise only. y is ALREADY a noisy survey estimate with
      # effective n; resampling from it as though it were truth double-counts
      # sampling error and made the f = 1 case disagree with the full survey by
      # 2.6 pp. The extra variance from fielding only a fraction f is
      #     var_extra = (1/f - 1) * p(1-p) / n_eff
      # which is zero at f = 1, so the curve now anchors correctly.
      pv <- pmin(pmax(cl$y, 1e-4), 1 - 1e-4)
      var_extra <- pmax((1 / f - 1) * pv * (1 - pv) / cl$n_eff, 0)
      y_sub <- pmin(pmax(cl$y + stats::rnorm(cl$n, 0, sqrt(var_extra)), 0), 1)
      # jackknifed regional mean: a district never contributes to its own anchor
      rs <- tapply(y_sub * nn, reg, sum)[as.character(reg)]
      rn <- tapply(nn, reg, sum)[as.character(reg)]
      denom <- rn - nn
      rmean <- ifelse(denom > 0, (rs - y_sub * nn) / denom,
                      stats::weighted.mean(y_sub, nn))
      p_survey <- rmean

      # DECOMPOSITION. The first run found regional-mean + model-pattern WORSE
      # than the regional mean alone, while the in-fill benchmark has the model
      # BEATING the regional mean. Those are only compatible if the model's
      # value is between-region and its within-region pattern is noise. So
      # score three things, not two:
      #   survey    regional mean only            (flat within region)
      #   aug_*     regional mean + model pattern (survey level, model pattern)
      #   alone_*   the model's own prediction    (model level AND pattern)
      anchor <- .v2_logit(pmin(pmax(rmean, 1e-4), 1 - 1e-4))
      p <- list(survey = p_survey)
      for (nm in names(mh)) {
        mhat <- mh[[nm]]
        mreg <- tapply(mhat, reg, mean)[as.character(reg)]
        p[[paste0("aug_", nm)]]   <- .v2_expit(anchor + (mhat - mreg))
        p[[paste0("alone_", nm)]] <- .v2_expit(mhat)
      }
      base <- is.finite(cl$y)
      for (nm in names(p)) {
        ok <- base & is.finite(p[[nm]])
        if (!any(ok)) next
        rows[[length(rows) + 1L]] <- data.frame(
          country = cn, outcome = on, fraction = f, rep = r, n_areas = cl$n,
          arm = nm,
          mae = 100 * mean(abs(cl$y[ok] - p[[nm]][ok])),
          rho = suppressWarnings(stats::cor(cl$y[ok], p[[nm]][ok], method = "spearman")),
          stringsAsFactors = FALSE)
      }
    }
  }
  cat("done", cn, on, "\n")
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "model_augmented_survey.csv"), row.names = FALSE)

CELL <- R |> group_by(country, outcome, fraction, arm) |>
  summarise(mae = mean(mae, na.rm = TRUE), rho = mean(rho, na.rm = TRUE),
            .groups = "drop")
SUMM <- CELL |> group_by(fraction, arm) |>
  summarise(cells = dplyr::n(),
            mae = round(median(mae, na.rm = TRUE), 2),
            rho = round(median(rho, na.rm = TRUE), 3), .groups = "drop")
write.csv(SUMM, file.path(OUTDIR, "model_augmented_survey_summary.csv"),
          row.names = FALSE)

cat("
===== GAP 4: DOES THE MODEL LET YOU FIELD A SMALLER SURVEY? =====
")
cat("median district MAE (pp) and rank correlation, by share of clusters fielded

")
print(as.data.frame(tidyr::pivot_wider(SUMM[, c("fraction","arm","mae")],
      names_from = arm, values_from = mae)), row.names = FALSE)
cat("
rank correlation:
")
print(as.data.frame(tidyr::pivot_wider(SUMM[, c("fraction","arm","rho")],
      names_from = arm, values_from = rho)), row.names = FALSE)

cat("
--- head-to-head vs survey-only, per cell (paired on identical draws) ---
")
W <- tidyr::pivot_wider(CELL[, c("country","outcome","fraction","arm","mae")],
                        names_from = arm, values_from = mae)
for (a in setdiff(unique(CELL$arm), "survey")) {
  d <- W[[a]] - W$survey                     # negative = model better
  cat(sprintf("%-12s better in %2d of %2d cell-fractions | median %+0.2f pp
",
              a, sum(d < 0, na.rm = TRUE), sum(is.finite(d)),
              median(d, na.rm = TRUE)))
}

cat("
--- VALIDATION: survey-only here vs the empirical ws5 curve ---
")
cat("MATCHED on the (country, outcome) cells present in BOTH, because a median
")
cat("over different cell sets is not a validation.
")
w <- tryCatch(read.csv("results/tables/anchoring_design_curve_DENSE.csv"),
              error = function(e) NULL)
if (!is.null(w)) {
  w$country <- ifelse(w$country == "Sierra Leone", "SierraLeone", w$country)
  w <- w[round(w$n_regions_anchored / w$n_regions_total, 1) == 1, ]
  ws <- w |> group_by(country, outcome, fraction_clusters) |>
    summarise(ws5 = median(mae_admin2_pp, na.rm = TRUE), .groups = "drop")
  mine <- CELL[CELL$arm == "survey", c("country","outcome","fraction","mae")]
  cmp <- merge(mine, ws, by.x = c("country","outcome","fraction"),
               by.y = c("country","outcome","fraction_clusters"))
  if (nrow(cmp)) {
    cmp$diff <- cmp$mae - cmp$ws5
    s <- cmp |> group_by(fraction) |>
      summarise(cells = dplyr::n(), this = round(median(mae), 2),
                ws5 = round(median(ws5), 2), diff = round(median(diff), 2),
                .groups = "drop")
    print(as.data.frame(s), row.names = FALSE)
    cat(sprintf("
matched cells: %d | median |diff| = %.2f pp
",
                nrow(cmp), median(abs(cmp$diff), na.rm = TRUE)))
    cat("PASS if the f = 1.00 difference is near zero: the two are then the same
")
    cat("quantity and the arm comparison above is calibrated.
")
  } else cat("no overlapping cells -- cannot validate
")
}
cat("
DONE
")
