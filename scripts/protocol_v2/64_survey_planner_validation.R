# =============================================================================
# scripts/protocol_v2/64_survey_planner_validation.R   [SP-01, 2026-09-27]
#
# IF A SURVEY CAN ONLY REACH k DISTRICTS, WHICH k SHOULD THEY BE?
#
# Retrospective validation of model-guided district selection, pre-registered
# in docs/findings/SP-01_SURVEY_PLANNER_DESIGN_2026-09-27.md (design and
# readings written before this ran). For each country x outcome cell: pretend
# only k of the surveyed districts get biomarkers, choose them by one of four
# rules, fit the zero-tuning index on those k, and score the resulting map on
# the districts held back. Selection never sees the target country's outcomes:
# the model-guided arms use the transported climate + soil index fitted on the
# other three countries (the pre-registered candidate), which is what a real
# planner would hold before fieldwork.
#
# Arms (40 draws each):
#   random          simple random sample of k                      (reference)
#   pps             k drawn with probability proportional to population
#   spread_model    transported score cut into k strata, one per stratum
#   extremes_model  k/2 from the top and k/2 from the bottom of the score
#   transport_only  no in-country fit at all: the transported score itself,
#                   scored on the same held-out districts
#
# Held-out sets differ by arm by construction (they are the unselected
# districts, which is the set a planner actually leaves unsurveyed); the
# extremes arm is therefore scored on a range-restricted middle, which is the
# decision-relevant comparison, not an artefact.
#
#   Rscript scripts/protocol_v2/64_survey_planner_validation.R
# -> results/tables/protocol_v2/survey_planner_validation.csv          (cell x arm x fraction)
#    results/tables/protocol_v2/survey_planner_validation_summary.csv  (arm x fraction over cells)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
set.seed(20260927L)

P2   <- "results/tables/protocol_v2"
REPS <- as.integer(Sys.getenv("SP_REPS", "40"))
FRACS <- c(0.20, 0.35, 0.50, 0.75)
CS_DOMS <- c("Climate and weather", "Soil characteristics")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
POP <- readRDS("dashboard/data/admin2_population.rds")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
CS    <- drop_near_outcome_v2(intersect(MD$column[MD$domain %in% CS_DOMS], names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
LBL <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
.key <- function(a, b) paste(trimws(a), trimws(b), sep = "|")

# ── the transported selection score, one per cell ────────────────────────────
# Pooled other-country fit at district level on climate + soil, outcome
# z-scored within country, basis oriented on training rows only (block C of
# scripts/policy_deck/10_viz_tables.R).
transported_score <- function(target_cn, on) {
  cl <- list()
  for (cn in COUNTRIES) {
    t <- TG[TG$country == cn & TG$outcome == on & is.finite(TG$y_level), ]
    m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", CS)], by = c("Admin1", "Admin2"))
    if (nrow(m) < 12) next
    cl[[cn]] <- list(n = nrow(m), y = m$y_level, X = prep_predictors_v2(as.matrix(m[, CS])),
                     Admin1 = m$Admin1, Admin2 = m$Admin2)
  }
  if (is.null(cl[[target_cn]]) || length(cl) < 3) return(NULL)   # need >= 2 training countries
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y))))
  tr <- which(ctry != target_cn); te <- which(ctry == target_cn)
  D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  p <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, D, NULL), error = function(e) rep(NA_real_, length(te)))
  if (all(!is.finite(p))) return(NULL)
  data.frame(Admin1 = cl[[target_cn]]$Admin1, Admin2 = cl[[target_cn]]$Admin2, tscore = p)
}

pick <- function(arm, n, k, tscore, popw, rep_seed) {
  set.seed(rep_seed)
  switch(arm,
    random = sample.int(n, k),
    pps = { w <- popw
            if (!any(is.finite(w) & w > 0)) w <- rep(1, n)
            w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0])
            sample.int(n, k, prob = w) },
    spread_model = { o <- order(tscore, sample.int(n))            # random tie-break inside the sort
                     grp <- cut(seq_len(n), breaks = k, labels = FALSE)
                     vapply(seq_len(k), function(g) { i <- o[grp == g]; i[sample.int(length(i), 1L)] }, 0L) },
    extremes_model = { o <- order(tscore, sample.int(n))
                       c(head(o, floor(k / 2)), tail(o, ceiling(k / 2))) })
}

rows <- list()
for (cn in COUNTRIES) for (on in setdiff(unique(TG$outcome[TG$country == cn]), c("child_zinc", "women_zinc"))) {
  t <- TG[TG$country == cn & TG$outcome == on & is.finite(TG$y_level) & is.finite(TG$y_prev), ]
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  n <- nrow(m); if (n < 14) next
  ts <- transported_score(cn, on)
  if (is.null(ts)) { cat(sprintf("skip %s %s: no transported score\n", cn, on)); next }
  tscore <- ts$tscore[match(.key(m$Admin1, m$Admin2), .key(ts$Admin1, ts$Admin2))]
  if (mean(is.finite(tscore)) < 0.9) { cat(sprintf("skip %s %s: score joins %d of %d\n", cn, on, sum(is.finite(tscore)), n)); next }
  tscore[!is.finite(tscore)] <- stats::median(tscore, na.rm = TRUE)
  pop <- POP[POP$country == LBL[[cn]], ]
  popw <- (if (startsWith(on, "child_")) pop$pop_child else pop$pop_women)[match(.key(m$Admin1, m$Admin2), .key(pop$Admin1, pop$Admin2))]
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS]))
  D  <- domain_representation_v2(Xr, domain_of)                   # unsupervised basis: X only, every district has X
  Yl <- as.numeric(scale(m$y_level)); Yp <- .v2_logit(m$y_prev)
  k5 <- function(x) max(1L, ceiling(0.2 * length(x)))             # "worst fifth" of the held-out set
  # The national prevalence is the survey's primary product (SP-01 amendment):
  # the reference is the population-weighted mean over ALL surveyed districts,
  # and each draw's estimate is the same estimator over the k visited ones.
  pw <- popw; if (!any(is.finite(pw) & pw > 0)) pw <- rep(1, n); pw[!is.finite(pw) | pw <= 0] <- min(pw[is.finite(pw) & pw > 0])
  nat_full <- sum(pw * m$y_prev) / sum(pw)
  nat_of <- function(sel) sum(pw[sel] * m$y_prev[sel]) / sum(pw[sel])
  cat(sprintf("%-12s %-13s n %3d  national %.3f\n", cn, on, n, nat_full))
  for (f in FRACS) {
    k <- max(5L, ceiling(f * n)); if (n - k < 5) next
    # fixed model-score strata for the stratified national estimator (spread
    # selection is one-per-stratum stratified sampling, a probability design;
    # inclusion probability 1/n_g inside stratum g, so the Hajek weight is
    # n_g x district population)
    strat_id <- as.integer(cut(rank(tscore, ties.method = "first"), breaks = k))
    n_g <- table(strat_id)
    for (arm in c("random", "pps", "spread_model", "extremes_model")) {
      sp <- cap <- mae <- rep(NA_real_, REPS)
      sp0 <- rep(NA_real_, REPS)                                   # transport_only on the same held-out set
      nat_naive <- nat_strat <- rep(NA_real_, REPS)
      for (r in seq_len(REPS)) {
        tr <- unique(pick(arm, n, k, tscore, popw, rep_seed = 90000L + r * 37L + match(arm, c("random", "pps", "spread_model", "extremes_model"))))
        te <- setdiff(seq_len(n), tr); if (length(te) < 5) next
        pl <- tryCatch(arm_domain_index_v2(tr, te, Yl, NULL, D, NULL), error = function(e) rep(NA_real_, length(te)))
        if (all(is.finite(pl))) sp[r] <- suppressWarnings(cor(pl, Yl[te], method = "spearman"))
        aux <- list(Admin1 = m$Admin1, y_nat = Yp, target = "prev")
        pp <- tryCatch(ARMS_V2[["domain_index_cal"]](tr, te, Yp, NULL, D, aux), error = function(e) rep(NA_real_, length(te)))
        if (all(is.finite(pp))) {
          mae[r] <- mean(abs(.v2_expit(pp) - m$y_prev[te])) * 100
          # worst fifth on the DEFICIENCY-PREVALENCE scale (orientation is
          # unambiguous there; the level column's sign follows the biomarker)
          kw <- k5(te)
          worst_true <- order(-m$y_prev[te])[seq_len(kw)]; worst_pred <- order(-pp)[seq_len(kw)]
          cap[r] <- length(intersect(worst_true, worst_pred)) / kw
        }
        sp0[r] <- suppressWarnings(cor(tscore[te], Yl[te], method = "spearman"))
        nat_naive[r] <- nat_of(tr)
        if (arm == "spread_model") {
          wht <- as.numeric(n_g[as.character(strat_id[tr])]) * pw[tr]
          nat_strat[r] <- sum(wht * m$y_prev[tr]) / sum(wht)
        }
      }
      rows[[length(rows) + 1]] <- data.frame(country = cn, outcome = on, n_surveyed = n, fraction = f, k = k, arm = arm,
        reps = sum(is.finite(sp)), spearman = mean(sp, na.rm = TRUE), spearman_sd = sd(sp, na.rm = TRUE),
        capture = mean(cap, na.rm = TRUE), mae_pp = mean(mae, na.rm = TRUE),
        spearman_transport_only = mean(sp0, na.rm = TRUE),
        nat_full_pp = 100 * nat_full,
        nat_bias_pp = 100 * (mean(nat_naive, na.rm = TRUE) - nat_full),
        nat_rmse_pp = 100 * sqrt(mean((nat_naive - nat_full)^2, na.rm = TRUE)),
        nat_ci_pp = 100 * 1.96 * sd(nat_naive, na.rm = TRUE),
        nat_strat_bias_pp = if (arm == "spread_model") 100 * (mean(nat_strat, na.rm = TRUE) - nat_full) else NA_real_,
        nat_strat_ci_pp = if (arm == "spread_model") 100 * 1.96 * sd(nat_strat, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
}
V <- bind_rows(rows)
write.csv(V, file.path(P2, "survey_planner_validation.csv"), row.names = FALSE)

SU <- V |> group_by(fraction, arm) |>
  summarise(cells = n(), spearman = mean(spearman, na.rm = TRUE), capture = mean(capture, na.rm = TRUE),
            mae_pp = mean(mae_pp, na.rm = TRUE),
            nat_abs_bias_pp = mean(abs(nat_bias_pp), na.rm = TRUE),
            nat_ci_pp = mean(nat_ci_pp, na.rm = TRUE),
            nat_strat_abs_bias_pp = mean(abs(nat_strat_bias_pp), na.rm = TRUE),
            nat_strat_ci_pp = mean(nat_strat_ci_pp, na.rm = TRUE), .groups = "drop")
ref <- V |> filter(arm == "random") |> select(country, outcome, fraction, sp_ref = spearman)
DL <- V |> filter(arm != "random") |> left_join(ref, by = c("country", "outcome", "fraction")) |>
  mutate(delta = spearman - sp_ref) |> group_by(fraction, arm) |>
  summarise(mean_delta = mean(delta, na.rm = TRUE), median_delta = median(delta, na.rm = TRUE),
            cells_better = sum(delta > 0, na.rm = TRUE), .groups = "drop")
SU <- left_join(SU, DL, by = c("fraction", "arm"))
write.csv(SU, file.path(P2, "survey_planner_validation_summary.csv"), row.names = FALSE)
cat("\nsummary (mean over cells; delta = against random on the same cells):\n")
print(as.data.frame(SU), row.names = FALSE, digits = 3)
cat(sprintf("\nwritten %s and _summary: %d cell rows\n", file.path(P2, "survey_planner_validation.csv"), nrow(V)))
