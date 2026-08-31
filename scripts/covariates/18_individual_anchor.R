# =============================================================================
# scripts/covariates/18_individual_anchor.R
#
# THE INDIVIDUAL-LEVEL ANCHOR, AND THE CLUSTER-LINKAGE QUESTION
#
# The area-level comparison scores ten approaches that all fit at the district.
# It has no upper bound: nothing says what would be achievable if the geospatial
# proxies were replaced by a household survey. This script adds that anchor, and
# uses the same run to test whether linking at the survey CLUSTER rather than
# the district buys anything.
#
# WHY AN ANCHOR IS NEEDED
# -----------------------
# A district-level r of 0.19 is uninterpretable on its own. Against a reliability
# ceiling of 0.13 it is most of what is attainable; against a questionnaire model
# reaching 0.60 it is a poor showing. Fitting the same learner on the same rows
# with the household questionnaire added settles which reading is right.
#
# The anchor must exclude the ASSAYS, not merely the outcome column. A first
# version used Xvars_full, which drops the outcome's own binary and continuous
# columns but keeps gw_cFER, gw_cSTFR, gw_cCRP and gw_cAGP -- and since iron
# deficiency is a ferritin cutoff, it scored r = 0.973, the outcome predicting
# itself. The question worth answering is what a household survey WITHOUT a
# blood draw buys, so every assay column is removed.
#
# The SECOND version scored r = 0.986, which is worse, because its local regex
# anchored the analyte tokens with `(^|_)` and the tokens sit after a
# population prefix -- `(^|_)FER` does not match `gw_wFER`. This now calls
# is_biomarker_column() from R/data_prep.R, the same predicate production uses,
# so the two cannot drift apart again.
#
# WHY THE CLUSTER QUESTION BELONGS HERE
# -------------------------------------
# Measured 2026-08-31: the proxy predictors are cluster-resolved (222 buffer
# covariates at 10/25/50 km), but they vary within a district in only 7% of
# Ghana's districts, because 62 of its 75 districts contain exactly ONE survey
# cluster. Individuals in a district therefore share one predictor vector, and
# the effective sample size for learning the predictor-outcome map is the number
# of CLUSTERS, not individuals: 70 / 90 / 60 / 103 against 30 / 75 / 14 / 87
# districts. Cluster linkage should help most where clusters outnumber districts
# by the widest margin, which is Sierra Leone (60 vs 14) -- the country whose
# admin-2 results are worst. This tests that directly.
#
# THE ARMS
#   1 proxy, district      proxy-only predictors, scored at admin-2
#   2 questionnaire, district  the SAME learner with the survey QUESTIONNAIRE
#                              items added -- wealth, education, diet, WASH,
#                              illness, supplementation. The anchor: what a
#                              household survey WITHOUT a blood draw achieves.
#   3 proxy, cluster       proxy-only, scored at the survey cluster
#   4 questionnaire, cluster
#
# All four use leave-one-region-out folds, so no district or cluster is
# predicted by a model that has seen its own region.
#
#   Rscript scripts/covariates/18_individual_anchor.R
# -> results/tables/individual_anchor.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")
OUTCOMES <- c("child_iron", "child_vitA", "women_iron", "women_vitA")
K_SCREEN <- 40L

num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# A light, honest learner. The production 12-learner stack would take hours over
# 16 cells x 4 arms, and the question here is the GAP between arms, which a
# consistent learner answers as well as an expensive one.
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

rows <- list()
for (cn in names(get_country_configs())) {
  cc <- get_country_configs()[[cn]]; lc <- tolower(cn)
  for (on in OUTCOMES) {
    oc <- cc$outcomes[[on]]; if (is.null(oc)) next
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE),
                   error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    if (!all(c("Admin2", oc$binary) %in% names(d))) next
    y <- num(d[[oc$binary]])
    keep <- is.finite(y)
    d <- d[keep, , drop = FALSE]; y <- y[keep]
    if (length(unique(y)) < 2 || nrow(d) < 100) next
    clu <- if (cc$cluster_id %in% names(d)) as.character(d[[cc$cluster_id]]) else NA
    blk <- if ("Admin1" %in% names(d)) as.character(d$Admin1) else rep("a", nrow(d))
    if (dplyr::n_distinct(blk) < 3) next

    for (arm in c("proxy", "quest")) {
      # THE ANCHOR MUST EXCLUDE THE BIOMARKERS, NOT JUST THE OUTCOME COLUMN.
      # Xvars_full drops the outcome's own binary and continuous columns but
      # KEEPS the raw assay panel -- gw_cFER, gw_cSTFR, gw_cCRP, gw_cAGP and
      # haemoglobin. Since iron deficiency IS a ferritin cutoff, the first run
      # of this script scored r = 0.973, which is the outcome predicting itself.
      # An anchor is only meaningful if it answers "what does a household survey
      # WITHOUT a blood draw buy", so every assay column goes.
      # Use the PRODUCTION guard, not a local regex. The first version here
      # anchored its tokens with `(^|_)` and scored r = 0.986 on Gambia women's
      # iron, because the analyte tokens sit AFTER a population prefix
      # (gw_wFER, gw_cSTFR) so the anchor never matched. is_biomarker_column()
      # in R/data_prep.R is the single definition of "the blood draw" and is
      # verified against all 1,297 gw_ columns; this must not drift from it.
      #
      # The filter is still applied HERE as well as in data_prep.R because the
      # stored outcome_data_* targets were built before that fix, so their
      # Xvars_full still carries the assay panel.
      vars <- if (arm == "proxy") od$Xvars else {
        full <- od$Xvars_full %||% od$Xvars
        full[!is_biomarker_column(full)]
      }
      X <- prep(d, vars); if (ncol(X) < 5) next
      oof <- rep(NA_real_, nrow(X))
      for (f in unique(blk)) {
        i <- which(blk == f); tr <- setdiff(seq_len(nrow(X)), i)
        if (!length(tr) || length(i) >= nrow(X)) next
        if (length(unique(y[tr])) < 2) next
        oof[i] <- fit_pred(X[tr, , drop = FALSE], y[tr], X[i, , drop = FALSE])
      }
      for (unit in c("district", "cluster")) {
        key <- if (unit == "district") paste(blk, d$Admin2) else clu
        if (all(is.na(key))) next
        ok <- is.finite(oof)
        agg <- data.frame(key = key[ok], obs = y[ok], pred = oof[ok]) %>%
          group_by(key) %>% summarise(n = dplyr::n(), obs = mean(obs),
                                      pred = mean(pred), .groups = "drop") %>%
          filter(n >= 5)
        if (nrow(agg) < 8) next
        rows[[length(rows) + 1L]] <- data.frame(
          country = cn, outcome = on, arm = arm, unit = unit,
          n_units = nrow(agg), n_pred = ncol(X),
          r = round(suppressWarnings(stats::cor(agg$obs, agg$pred)), 3),
          mae_pp = round(100 * mean(abs(agg$obs - agg$pred)), 2),
          stringsAsFactors = FALSE)
        cat(sprintf("  %-12s %-11s %-5s %-8s units=%3d p=%4d  r=%+.3f MAE=%5.2f\n",
                    cn, on, arm, unit, nrow(agg), ncol(X),
                    stats::cor(agg$obs, agg$pred), 100 * mean(abs(agg$obs - agg$pred))))
      }
    }
  }
}

res <- dplyr::bind_rows(rows)
readr::write_csv(res, here("results", "tables", "individual_anchor.csv"))

cat("\n========= THE ANCHOR: proxy vs household questionnaire (no assay) =====\n")
print(res %>% group_by(unit, arm) %>%
        summarise(cells = n(), mean_r = round(mean(r, na.rm = TRUE), 3),
                  med_r = round(stats::median(r, na.rm = TRUE), 3),
                  mae = round(mean(mae_pp, na.rm = TRUE), 2), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

w <- res %>% select(country, outcome, unit, arm, r) %>%
  tidyr::pivot_wider(names_from = arm, values_from = r)
if (all(c("proxy", "full") %in% names(w))) {
  cat("\n--- how much does having the survey buy, per unit? ---\n")
  print(w %>% mutate(gap = round(full - proxy, 3)) %>% group_by(unit) %>%
          summarise(cells = n(), mean_gap = round(mean(gap, na.rm = TRUE), 3),
                    quest_better = sprintf("%d/%d", sum(gap > 0, na.rm = TRUE), n()),
                    .groups = "drop") %>% as.data.frame(), row.names = FALSE)
}

cat("\n================ DOES CLUSTER LINKAGE HELP? ==========================\n")
u <- res %>% select(country, outcome, arm, unit, r, n_units) %>%
  tidyr::pivot_wider(names_from = unit, values_from = c(r, n_units))
if (all(c("r_district", "r_cluster") %in% names(u))) {
  print(u %>% mutate(gain = round(r_cluster - r_district, 3)) %>%
          group_by(country) %>%
          summarise(cells = n(), districts = max(n_units_district, na.rm = TRUE),
                    clusters = max(n_units_cluster, na.rm = TRUE),
                    mean_gain = round(mean(gain, na.rm = TRUE), 3),
                    better = sprintf("%d/%d", sum(gain > 0, na.rm = TRUE), n()),
                    .groups = "drop") %>% as.data.frame(), row.names = FALSE)
}
cat("\n-> results/tables/individual_anchor.csv\n")
