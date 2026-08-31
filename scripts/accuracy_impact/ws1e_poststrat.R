# =============================================================================
# scripts/accuracy_impact/ws1e_poststrat.R
#
# WS1e. Does the individual-level model predict the DISTRICT, or the SAMPLE the
# survey happened to draw from it?
#
# WHY THIS MATTERS FOR THE CEILING
# --------------------------------
# The analytic ceiling bounds a predictor that is independent of the sampling
# realisation. An estimator that sees respondent-level or cluster-level
# information and is then averaged over the REALISED sample is not independent
# of it: its district aggregate tracks who happened to be interviewed, and the
# observed district prevalence tracks the same people. Two quantities sharing a
# sampling realisation can correlate above the ceiling without any of that
# correlation being transportable skill.
#
# WHAT STAGE 0 ESTABLISHED, AND WHAT IT LEAVES
# --------------------------------------------
# The proxy arm's predictors are Admin-2 and cluster resolved: Xvars is
# Xvars_all minus the entire GW domain, so it holds no respondent-level
# covariate. Individual predictions are therefore CONSTANT WITHIN A CLUSTER, and
# the only way the realised sample can enter the district aggregate is through
# the weights the clusters receive. That is a narrower channel than the general
# case, and this script measures it rather than assuming it is negligible.
#
# THE 2 x 2
# ---------
# Aggregation of the out-of-fold predictions:
#   realised   unweighted mean over the respondents actually sampled, which is
#              what scripts/covariates/18_individual_anchor.R does
#   poststrat  survey-weighted mean, which estimates the district POPULATION
#              rather than the achieved sample
# Target scored against:
#   sample_mean   the unweighted district mean of the outcome
#   svy_weighted  the survey-weighted district mean of the outcome
#   svy_prev      the design-based district prevalence, which is the estimand
#
# WorldPop supplies total population only for these countries, with no age-sex
# structure, so the population composition used is the survey's own weighted
# composition. That is labelled here rather than presented as a census frame.
#
#   PROFILE=smoke   Ghana child_iron only
#   Rscript scripts/accuracy_impact/ws1e_poststrat.R
# -> results/tables/poststratified_aggregation.csv
# =============================================================================
suppressPackageStartupMessages({library(here); library(dplyr)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full")
SUF <- switch(PROFILE, smoke = "_SMOKE", sierraleone = "_SLE", "")
K_SCREEN <- 40L; SEED <- 20260903L; MIN_N <- 5L
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# The same light NNLS stack the individual anchor uses, so this measures the
# aggregation and not a change of learner.
fit_pred <- function(Xtr, ytr, Xte) {
  sel <- .awsl_screen(Xtr, ytr, K_SCREEN)
  s <- tryCatch(.awsl_stack(Xtr[, sel, drop = FALSE], ytr, rep(1, length(ytr))),
                error = function(e) NULL)
  if (is.null(s)) return(rep(NA_real_, nrow(Xte)))
  .awsl_predict(s, Xte[, sel, drop = FALSE])
}
prep <- function(d, vars) {
  vars <- intersect(vars, names(d))
  vars <- vars[vapply(vars, function(v) is.numeric(d[[v]]) || is.logical(d[[v]]) ||
                        inherits(d[[v]], "haven_labelled"), logical(1))]
  X <- vapply(vars, function(v) num(d[[v]]), numeric(nrow(d)))
  if (is.null(dim(X))) X <- matrix(X, nrow = nrow(d))
  colnames(X) <- vars
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  X[, apply(X, 2, function(z) stats::sd(z) > 0), drop = FALSE]
}

cfgs <- get_country_configs()
# PROFILE=sierraleone runs only the country where post-stratification CAN
# change anything: Gambia, Ghana and Malawi have a within-district survey-weight
# coefficient of variation of exactly 0, so the weighted and unweighted district
# means are identical there by construction and the full sweep spends hours
# confirming an identity. See docs/TARGET_ESTIMAND.md section 4.
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]
if (PROFILE == "sierraleone") cfgs <- cfgs["SierraLeone"]
rows <- list()
set.seed(SEED)
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  ocs <- names(cc$outcomes)
  if (PROFILE == "smoke") ocs <- intersect(ocs, "child_iron")
  for (on in ocs) {
    oc <- cc$outcomes[[on]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE),
                   error = function(e) NULL)
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on), store = STORE),
                   error = function(e) NULL)
    if (is.null(od) || is.null(oc$binary)) next
    d <- od$data
    if (!all(c("Admin2", oc$binary) %in% names(d))) next
    y <- num(d[[oc$binary]]); keep <- is.finite(y)
    d <- d[keep, , drop = FALSE]; y <- y[keep]
    if (length(unique(y)) < 2 || nrow(d) < 100) next
    blk <- if ("Admin1" %in% names(d)) as.character(d$Admin1) else rep("a", nrow(d))
    if (dplyr::n_distinct(blk) < 3) next
    w <- if (cc$weight_col %in% names(d)) num(d[[cc$weight_col]]) else rep(1, nrow(d))
    w[!is.finite(w) | w <= 0] <- NA_real_
    w[is.na(w)] <- stats::median(w, na.rm = TRUE)

    X <- prep(d, od$Xvars); if (ncol(X) < 5) next
    oof <- rep(NA_real_, nrow(X))
    for (f in unique(blk)) {
      i <- which(blk == f); tr <- setdiff(seq_len(nrow(X)), i)
      if (!length(tr) || length(i) >= nrow(X) || length(unique(y[tr])) < 2) next
      oof[i] <- fit_pred(X[tr, , drop = FALSE], y[tr], X[i, , drop = FALSE])
    }
    ok <- is.finite(oof)
    if (sum(ok) < 50) next
    agg <- data.frame(a1 = blk[ok], a2 = as.character(d$Admin2)[ok],
                      obs = y[ok], pred = oof[ok], w = w[ok],
                      stringsAsFactors = FALSE) |>
      group_by(a1, a2) |>
      summarise(n = dplyr::n(),
                obs_sample = mean(obs),
                obs_wt     = stats::weighted.mean(obs, w),
                pred_real  = mean(pred),
                pred_ps    = stats::weighted.mean(pred, w),
                .groups = "drop") |>
      filter(n >= MIN_N)
    if (nrow(agg) < 8) next

    # svy_prev is the design-based estimand. Join on the pair key where the
    # survey table carries Admin1, per R/admin2_key_hygiene.R.
    agg$svy_prev <- NA_real_
    if (!is.null(sv)) {
      aggk <- data.frame(Admin1 = agg$a1, Admin2 = agg$a2, stringsAsFactors = FALSE)
      by <- admin2_join_by(sv, aggk)
      j <- if (identical(by, c("Admin1", "Admin2")))
        match(paste(agg$a1, agg$a2), paste(sv$Admin1, sv$Admin2)) else
        match(agg$a2, sv$Admin2)
      agg$svy_prev <- sv$svy_prev[j]
    }

    for (aggn in c("realised", "poststrat")) {
      p <- if (aggn == "realised") agg$pred_real else agg$pred_ps
      for (tg in c("sample_mean", "svy_weighted", "svy_prev")) {
        o <- switch(tg, sample_mean = agg$obs_sample,
                    svy_weighted = agg$obs_wt, svy_prev = agg$svy_prev)
        fin <- is.finite(p) & is.finite(o)
        if (sum(fin) < 8 || stats::sd(o[fin]) == 0) next
        rows[[length(rows) + 1L]] <- data.frame(
          country = cc$country, outcome = on, aggregation = aggn, target = tg,
          n_units = sum(fin),
          r = round(suppressWarnings(stats::cor(p[fin], o[fin])), 4),
          mae_pp = round(100 * mean(abs(p[fin] - o[fin])), 2),
          stringsAsFactors = FALSE)
      }
    }
    cat(sprintf("  [ok] %-13s %-13s units=%d\n", cn, on, nrow(agg)))
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
readr::write_csv(res, here("results", "tables",
                           sprintf("poststratified_aggregation%s.csv", SUF)))

cat("\n=== WS1e: aggregation against target ===\n")
s <- res |> group_by(aggregation, target) |>
  summarise(cells = dplyr::n(), mean_r = round(mean(r, na.rm = TRUE), 4),
            med_r = round(stats::median(r, na.rm = TRUE), 4),
            mae = round(mean(mae_pp, na.rm = TRUE), 2), .groups = "drop")
print(as.data.frame(s), row.names = FALSE)

cat("\n--- paired: post-stratified minus realised, same cells and target ---\n")
w <- res |> select(country, outcome, target, aggregation, r) |>
  tidyr::pivot_wider(names_from = aggregation, values_from = r)
if (all(c("realised", "poststrat") %in% names(w))) {
  w$delta <- round(w$poststrat - w$realised, 4)
  print(as.data.frame(w |> group_by(target) |>
    summarise(cells = dplyr::n(), mean_delta = round(mean(delta, na.rm = TRUE), 4),
              med_delta = round(stats::median(delta, na.rm = TRUE), 4),
              n_lower = sum(delta < 0, na.rm = TRUE), .groups = "drop")),
    row.names = FALSE)
}
