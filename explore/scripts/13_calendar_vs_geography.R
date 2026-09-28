# =============================================================================
# explore/scripts/13_calendar_vs_geography.R   [probe CG-01]
#
# IS THE SPATIAL SMOOTHER FITTING GEOGRAPHY, OR THE FIELDWORK CALENDAR?
#
# The individual-level exploration (2026-09-28) established that a survey team
# visits a cluster in about one day - median within-cluster date spread 1 day,
# max 6 - so collection date is almost the same variable as place:
#   R2(date ~ cluster)  = 0.999 Malawi, 0.989 Gambia
#   R2(date ~ district) = 0.908 Malawi
# and the raw date slopes are not zero: Malawi log RBP -0.0009/day (t = -2.3),
# AGP +0.0013/day (t = 2.0), i.e. -6% and +9% across the 70-day window.
#
# Fieldwork teams move through space contiguously, so the fieldwork calendar is
# itself a SMOOTH SPATIAL SURFACE. That raises a specific threat to the
# project's standing conclusion 2 - "within a surveyed country, geography alone
# captures nearly all district-level accuracy; covariates add nothing on top of
# a spatial smoother". The smoother may be fitting the survey schedule.
#
# THE TEST
#   1. How much between-district outcome variance does the district's own
#      fieldwork date explain, per cell?
#   2. How spatial is the calendar? R2 of date on a thin-plate spline in
#      lon/lat - if high, a spatial smoother can represent it.
#   3. THE DECISIVE ONE. Re-run the in-fill comparison on the outcome
#      RESIDUALISED ON FIELDWORK DATE (residualisation inside the training fold
#      only, so no leakage). If `spatial` beats `domain_index` on the raw
#      outcome but not on the date-residualised outcome, its advantage was
#      calendar, not geography.
#
#   Rscript explore/scripts/13_calendar_vs_geography.R
# -> explore/out/13_calendar_variance.csv, 13_calendar_arms.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
suppressPackageStartupMessages({library(mgcv)})

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
ix <- exp_cell_index(E)

FW <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/fieldwork_windows_admin2.csv"),
               stringsAsFactors = FALSE)
FW$date_med <- as.Date(FW$date_med)
FW$day <- NA_real_
for (cn in unique(FW$country)) {
  k <- FW$country == cn
  FW$day[k] <- as.numeric(FW$date_med[k] - min(FW$date_med[k], na.rm = TRUE))
}

#' attach the district's fieldwork day to a cell, on the pair key
with_day <- function(cell) {
  f <- FW[FW$country == cell$country, c("Admin1", "Admin2", "day")]
  i <- match(paste(cell$Admin1, cell$Admin2), paste(f$Admin1, f$Admin2))
  cell$day <- f$day[i]
  cell
}

# ── 1 & 2. how much variance is calendar, and how spatial is the calendar ────
vrows <- list()
for (i in seq_len(nrow(ix))) {
  cn <- ix$country[i]; on <- ix$outcome[i]
  for (tgt in c("level", "prev")) {
    cell <- tryCatch(with_day(exp_cell(E, cn, on, tgt)), error = function(e) NULL)
    if (is.null(cell)) next
    ok <- is.finite(cell$day) & is.finite(cell$y_mod)
    if (sum(ok) < 12) next
    y <- cell$y_mod[ok]; d <- cell$day[ok]
    lo <- cell$aux$lon[ok]; la <- cell$aux$lat[ok]

    r2_date <- summary(stats::lm(y ~ d))$r.squared
    # how spatial is the calendar itself?
    k <- max(5, min(30, floor(sum(ok) / 3)))
    r2_date_geo <- tryCatch(
      summary(mgcv::gam(d ~ s(lo, la, k = k)))$r.sq, error = function(e) NA_real_)
    # how much of the outcome does geography explain, and does date add to it?
    r2_geo <- tryCatch(summary(mgcv::gam(y ~ s(lo, la, k = k)))$r.sq,
                       error = function(e) NA_real_)
    r2_geo_date <- tryCatch(summary(mgcv::gam(y ~ s(lo, la, k = k) + d))$r.sq,
                            error = function(e) NA_real_)
    vrows[[paste(cn, on, tgt)]] <- data.frame(
      country = cn, outcome = on, target = tgt, n = sum(ok),
      r2_outcome_on_date = r2_date,
      r2_date_on_geography = r2_date_geo,
      r2_outcome_on_geography = r2_geo,
      r2_outcome_on_geo_plus_date = r2_geo_date,
      date_adds_over_geo = r2_geo_date - r2_geo,
      stringsAsFactors = FALSE)
  }
  message("  variance ", cn, " ", on)
}
V <- dplyr::bind_rows(vrows); exp_write(V, "13_calendar_variance")

# ── 3. the decisive test: arms on the date-residualised outcome ─────────────
# Residualising uses the TRAINING rows only: the slope of y on day is fitted on
# tr and applied to te, so the held-out district's own outcome never informs it.
resid_arm <- function(inner) function(tr, te, y, X, D, aux) {
  d <- aux$day
  ok <- is.finite(d[tr]) & is.finite(y[tr])
  if (sum(ok) < 10) return(inner(tr, te, y, X, D, aux))
  b <- tryCatch(stats::lm(y[tr][ok] ~ d[tr][ok])$coefficients, error = function(e) NULL)
  if (is.null(b) || !all(is.finite(b))) return(inner(tr, te, y, X, D, aux))
  yr <- y - (b[1] + b[2] * ifelse(is.finite(d), d, mean(d[tr], na.rm = TRUE)))
  inner(tr, te, yr, X, D, aux) + (b[1] + b[2] * ifelse(is.finite(d[te]), d[te],
                                                       mean(d[tr], na.rm = TRUE)))
}
#' an arm that sees ONLY the fieldwork date
arm_date_only <- function(tr, te, y, X, D, aux) {
  d <- aux$day
  ok <- is.finite(d[tr]) & is.finite(y[tr])
  if (sum(ok) < 10) return(rep(mean(y[tr]), length(te)))
  b <- stats::lm(y[tr][ok] ~ d[tr][ok])$coefficients
  if (!all(is.finite(b))) return(rep(mean(y[tr]), length(te)))
  as.numeric(b[1] + b[2] * ifelse(is.finite(d[te]), d[te], mean(d[tr][ok])))
}

ARMS <- list(
  null_train_mean = ARMS_V2$null_train_mean,
  spatial         = ARMS_V2$spatial,
  domain_index    = ARMS_V2$domain_index,
  date_only       = arm_date_only,
  # the same two arms, scored on the outcome with the calendar trend removed
  spatial_dtr     = resid_arm(ARMS_V2$spatial),
  index_dtr       = resid_arm(ARMS_V2$domain_index)
)

rows <- list()
for (i in seq_len(nrow(ix))) {
  cn <- ix$country[i]; on <- ix$outcome[i]
  for (tgt in c("level", "prev")) {
    cell <- tryCatch(with_day(exp_cell(E, cn, on, tgt)), error = function(e) NULL)
    if (is.null(cell) || !sum(is.finite(cell$day))) next
    cell$aux$day <- cell$day
    rows[[paste(i, tgt)]] <- exp_infill(cell, ARMS, reps = REPS)
  }
  message("  arms ", cn, " ", on)
}
A <- exp_summarise(dplyr::bind_rows(rows)); exp_write(A, "13_calendar_arms")

# ── report ──────────────────────────────────────────────────────────────────
cat("\n== 1. how much between-district outcome variance is the FIELDWORK CALENDAR? ==\n")
p <- aggregate(cbind(r2_outcome_on_date, r2_date_on_geography, r2_outcome_on_geography,
                     date_adds_over_geo) ~ country + target, data = V,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
print(p[order(p$target, -p$r2_outcome_on_date), ], row.names = FALSE)

cat("\n== 2. is the calendar itself a smooth spatial surface? (per country) ==\n")
q <- aggregate(r2_date_on_geography ~ country, data = V,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
print(q[order(-q$r2_date_on_geography), ], row.names = FALSE)

cat("\n== 3. THE TEST: does the smoother's edge survive removing the calendar? ==\n")
for (tgt in c("level", "prev")) {
  s <- A[A$estimand == "infill" & A$target == tgt, ]
  w <- reshape(s[, c("country", "outcome", "arm", "spearman")],
               idvar = c("country", "outcome"), timevar = "arm", direction = "wide")
  names(w) <- sub("^spearman[.]", "", names(w))
  cat(sprintf("\n-- %s --\n", tgt))
  cat(sprintf("  raw outcome        : spatial %.3f  index %.3f  gap %+.3f  (spatial wins %d/%d)\n",
      mean(w$spatial, na.rm = TRUE), mean(w$domain_index, na.rm = TRUE),
      mean(w$spatial - w$domain_index, na.rm = TRUE),
      sum(w$spatial > w$domain_index, na.rm = TRUE), sum(is.finite(w$spatial))))
  cat(sprintf("  calendar removed   : spatial %.3f  index %.3f  gap %+.3f  (spatial wins %d/%d)\n",
      mean(w$spatial_dtr, na.rm = TRUE), mean(w$index_dtr, na.rm = TRUE),
      mean(w$spatial_dtr - w$index_dtr, na.rm = TRUE),
      sum(w$spatial_dtr > w$index_dtr, na.rm = TRUE), sum(is.finite(w$spatial_dtr))))
  cat(sprintf("  date ALONE         : %.3f\n", mean(w$date_only, na.rm = TRUE)))
  cat(sprintf("  smoother's loss from removing the calendar: %+.3f\n",
      mean(w$spatial_dtr - w$spatial, na.rm = TRUE)))
}
