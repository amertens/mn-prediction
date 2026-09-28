# =============================================================================
# explore/scripts/14_timing_artefacts.R   [probes MA-01 and PT-01]
#
# MA-01  Do collection-timing artefacts inject NOISE into district estimates,
#        and does removing them raise the district split-half reliability?
# PT-01  Can collection timing improve INDIVIDUAL-level prediction?
#
# WHAT THE EXPLORATION ALREADY ESTABLISHED (2026-09-28), which sets the
# expectation for both:
#   - a team visits a cluster in ~1 day (median within-cluster date spread 1,
#     max 6), so date has almost no within-cluster variance and is 90-100%
#     explained by place
#   - the artefacts that DO vary within cluster are small: adjusting Malawi
#     district means for time of blood draw moves them by at most 0.013 log
#     units and leaves the ranking at Spearman 0.999
#   - Gambia days-since-VAS is mechanistically right (+3.9% RBP if VAS within
#     60 days) but n = 276 gives t = 1.0
#   - fasting is unusable: 65 of 3,099 respondents fasted
# Both probes are therefore expected to be null. They are run to CLOSE the
# family explicitly rather than leave it inferred.
#
# MA-01 DESIGN. Split respondents within each district at random into halves,
# take each half's district mean, and correlate the two halves across
# districts: that is the reliability of a district estimate at half sample.
# Do it on the raw biomarker and on the biomarker with the artefact block
# (time of draw, fasting, days since VAS) regressed out within cluster. If the
# artefacts are noise, removing them should RAISE the split-half correlation.
# 200 random splits.
#
# PT-01 DESIGN. Individual-level ridge predicting the deficiency indicator,
# cluster-blocked 5-fold CV, with and without the timing block. Scored by AUC
# and Brier skill against the training-fold prevalence.
#
#   Rscript explore/scripts/14_timing_artefacts.R
# -> explore/out/14_split_half_reliability.csv, 14_person_timing.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
suppressPackageStartupMessages({library(glmnet)})

NSPLIT <- as.integer(Sys.getenv("EXP_SPLITS", "200"))
set.seed(20260928L)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# ── assemble individual-level records with timing, per country ──────────────
recs <- list()

## Malawi: date, time of blood draw, fasting
m <- readRDS(file.path(EXP_ROOT, "data/IPD/Malawi/clean_malawi_mn_data.RDS"))
m$day <- as.numeric(m$date_interview - min(m$date_interview, na.rm = TRUE))
for (bm in c("rbp", "vitb12", "zn_gdl")) {
  if (!bm %in% names(m)) next
  y <- num(m[[bm]])
  keep <- is.finite(y) & !is.na(m$cluster) & !is.na(m$Admin2)
  if (sum(keep) < 200) next
  recs[[paste("Malawi", bm)]] <- data.frame(
    country = "Malawi", marker = bm,
    y = log(pmax(y[keep], 1e-6)),
    cluster = as.character(m$cluster[keep]), district = as.character(m$Admin2[keep]),
    pm = as.integer(num(m$time_blood_draw[keep]) == 2),
    fast = as.integer(num(m$fast[keep]) == 1),
    vasgap = NA_real_, day = m$day[keep], stringsAsFactors = FALSE)
}

## Gambia: interview date and days since vitamin A supplementation
g <- readRDS(file.path(EXP_ROOT, "data/IPD/Gambia/Gambia_merged_dataset.rds"))
ci <- as.Date(num(g$gw_cIntDate), origin = "1960-01-01")
vd <- as.Date(num(g$gw_cVASDate), origin = "1960-01-01")
gap <- as.numeric(ci - vd); gap[!is.finite(gap) | gap < 0 | gap > 365] <- NA
rbp <- num(g$gw_cRBP)
dcol <- if ("Admin2" %in% names(g)) "Admin2" else grep("Admin2|district", names(g), value = TRUE)[1]
keep <- is.finite(rbp) & !is.na(g$gw_MICS_Cluster_Number) & !is.na(g[[dcol]])
if (sum(keep) > 100) recs[["Gambia rbp"]] <- data.frame(
  country = "Gambia", marker = "gw_cRBP",
  y = log(pmax(rbp[keep], 1e-6)),
  cluster = as.character(g$gw_MICS_Cluster_Number[keep]),
  district = as.character(g[[dcol]][keep]),
  pm = NA_integer_, fast = NA_integer_, vasgap = gap[keep],
  day = as.numeric(ci[keep] - min(ci, na.rm = TRUE)), stringsAsFactors = FALSE)

R <- dplyr::bind_rows(recs)
message("assembled ", nrow(R), " records over ", dplyr::n_distinct(R$marker), " markers")

# ── MA-01: split-half reliability, raw vs artefact-adjusted ─────────────────
#' remove the artefact block within cluster (cluster fixed effects), so the
#' adjustment cannot absorb between-district signal
adjust_within_cluster <- function(d) {
  vars <- c("pm", "fast", "vasgap")
  use <- vars[vapply(vars, function(v) sum(is.finite(d[[v]])) > 50 &&
                       length(unique(stats::na.omit(d[[v]]))) > 1, TRUE)]
  if (!length(use)) return(list(y = d$y, used = character(0)))
  f <- stats::as.formula(paste("y ~ factor(cluster) +", paste(use, collapse = " + ")))
  dd <- d; for (v in use) dd[[v]][!is.finite(dd[[v]])] <- stats::median(dd[[v]], na.rm = TRUE)
  fit <- tryCatch(stats::lm(f, dd), error = function(e) NULL)
  if (is.null(fit)) return(list(y = d$y, used = character(0)))
  cf <- stats::coef(fit); cf <- cf[names(cf) %in% use]
  cf[!is.finite(cf)] <- 0
  adj <- d$y
  for (v in names(cf)) adj <- adj - cf[[v]] * (dd[[v]] - mean(dd[[v]], na.rm = TRUE))
  list(y = adj, used = names(cf))
}

split_half <- function(y, district, nsplit = NSPLIT, min_n = 6L) {
  d <- data.frame(y = y, district = district)
  ok <- table(d$district) >= min_n
  d <- d[d$district %in% names(ok)[ok], ]
  if (dplyr::n_distinct(d$district) < 8) return(NA_real_)
  out <- numeric(nsplit)
  for (b in seq_len(nsplit)) {
    h <- stats::ave(seq_len(nrow(d)), d$district,
                    FUN = function(i) sample(rep(c(1, 2), length.out = length(i))))
    a <- tapply(d$y[h == 1], d$district[h == 1], mean)
    bb <- tapply(d$y[h == 2], d$district[h == 2], mean)
    k <- intersect(names(a), names(bb))
    out[b] <- suppressWarnings(stats::cor(a[k], bb[k], method = "spearman"))
  }
  mean(out, na.rm = TRUE)
}

mrows <- list()
for (key in unique(paste(R$country, R$marker))) {
  d <- R[paste(R$country, R$marker) == key, ]
  aj <- adjust_within_cluster(d)
  r_raw <- split_half(d$y, d$district)
  r_adj <- split_half(aj$y, d$district)
  mrows[[key]] <- data.frame(
    country = d$country[1], marker = d$marker[1], n = nrow(d),
    districts = dplyr::n_distinct(d$district),
    artefacts_used = paste(aj$used, collapse = "+"),
    split_half_raw = r_raw, split_half_adjusted = r_adj,
    gain = r_adj - r_raw, stringsAsFactors = FALSE)
  message("  split-half ", key)
}
M <- dplyr::bind_rows(mrows); exp_write(M, "14_split_half_reliability")

# ── PT-01: does timing help INDIVIDUAL prediction? ──────────────────────────
auc <- function(y, p) {
  o <- order(p); y <- y[o]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p[o], ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

prows <- list()
for (key in unique(paste(R$country, R$marker))) {
  d <- R[paste(R$country, R$marker) == key, ]
  d$def <- as.integer(d$y < stats::quantile(d$y, 0.25, na.rm = TRUE))  # bottom quartile
  cls <- unique(d$cluster)
  if (length(cls) < 10) next
  for (blk in c("none", "timing")) {
    Xl <- if (blk == "none") NULL else {
      v <- c("pm", "fast", "vasgap", "day")
      v <- v[vapply(v, function(z) sum(is.finite(d[[z]])) > 50 &&
                      length(unique(stats::na.omit(d[[z]]))) > 1, TRUE)]
      if (!length(v)) NULL else {
        mm <- as.matrix(d[, v, drop = FALSE])
        for (j in seq_len(ncol(mm))) mm[!is.finite(mm[, j]), j] <-
          stats::median(mm[, j], na.rm = TRUE)
        mm
      }
    }
    set.seed(11)
    fold_of <- stats::setNames(sample(rep(1:5, length.out = length(cls))), cls)
    folds <- fold_of[d$cluster]
    pred <- rep(NA_real_, nrow(d))
    for (f in 1:5) {
      te <- which(folds == f); tr <- which(folds != f)
      if (!length(te) || length(tr) < 50) next
      if (is.null(Xl)) { pred[te] <- mean(d$def[tr]); next }
      fit <- tryCatch(glmnet::cv.glmnet(Xl[tr, , drop = FALSE], d$def[tr],
                                        family = "binomial", alpha = 0, nfolds = 5),
                      error = function(e) NULL)
      pred[te] <- if (is.null(fit)) mean(d$def[tr]) else
        as.numeric(stats::predict(fit, newx = Xl[te, , drop = FALSE],
                                  s = "lambda.min", type = "response"))
    }
    ok <- is.finite(pred)
    base <- mean(d$def[ok])
    prows[[paste(key, blk)]] <- data.frame(
      country = d$country[1], marker = d$marker[1], block = blk, n = sum(ok),
      auc = auc(d$def[ok], pred[ok]),
      brier = mean((d$def[ok] - pred[ok])^2),
      brier_skill = 1 - mean((d$def[ok] - pred[ok])^2) / mean((d$def[ok] - base)^2),
      cols = if (is.null(Xl)) 0L else ncol(Xl), stringsAsFactors = FALSE)
  }
  message("  person-level ", key)
}
P <- dplyr::bind_rows(prows); exp_write(P, "14_person_timing")

cat("\n== MA-01: district split-half reliability, raw vs timing-adjusted ==\n")
print(M, row.names = FALSE, digits = 3)
cat(sprintf("\nmean gain from removing the timing artefacts: %+.4f (%d of %d markers improved)\n",
            mean(M$gain, na.rm = TRUE), sum(M$gain > 0, na.rm = TRUE), sum(is.finite(M$gain))))

cat("\n== PT-01: individual-level prediction from timing alone ==\n")
print(P, row.names = FALSE, digits = 3)
cat("\n(AUC 0.5 and Brier skill 0 are the no-information values.)\n")
