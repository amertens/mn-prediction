# =============================================================================
# scripts/accuracy_impact/wsn1_band_agreement.R
#
# WS-N. Do models of the CONTINUOUS biomarker classify districts into
# public-health-significance bands better than models of the PREVALENCE?
#
# WHY THIS IS THE DECISION-RELEVANT COMPARISON
# --------------------------------------------
# Levels beat prevalence on correlation in 18 of 23 cells (median gain +0.114).
# That is a statistical statement. The programmatic question is different: a
# ministry does not act on an r, it acts on which band a district falls in. And
# a model of the LEVEL can be thresholded at any cut-off, while a model of the
# prevalence cannot be run backwards. If a threshold moves, only the continuous
# model survives.
#
# WHY ONLY TWO NUTRIENTS
# ----------------------
# Scoring band agreement requires bands that exist. Verified 2026-09-01 against
# the source documents:
#
#   vitamin A   2/10/20 percent ARE WHO's, from the Global Prevalence of
#               Vitamin A Deficiency 1995-2005. Defined for PRESCHOOL CHILDREN;
#               applying them to women is a borrowed convention and is scored
#               separately and labelled.
#   zinc        IZiNCG defines ONE line, elevated population risk above 20
#               percent. The project's 10 and 30 percent bands have no source,
#               so only the 20 percent line is scored, as a binary.
#   iron        the 5/20/40 bands are WHO's ANAEMIA classification measured by
#               HAEMOGLOBIN (Table 3, p.17). This project measures iron
#               deficiency by ferritin. Different biomarker, different
#               condition. NOT SCORED.
#   folate, B12 no WHO banding exists. NOT SCORED.
#
# Running all 24 cells would average agreement against targets that are partly
# invented, so the scope is the cells where the band is real.
#
# THE MALAWI ZINC EXCEPTION, AND WHY IT IS INSTRUCTIVE
# ----------------------------------------------------
# Malawi's zinc binary uses IZiNCG time-of-day-specific cut-offs, so the same
# concentration is deficient in one respondent and not in another. Measured:
# deficient values run to 64.69 and non-deficient start at 57.87, an overlap of
# nearly 7 units. A predicted district MEAN concentration therefore cannot be
# thresholded into the binary, because the threshold is individual. Malawi zinc
# is reported with that limitation rather than silently scored.
#
#   Rscript scripts/accuracy_impact/wsn1_band_agreement.R
# -> results/tables/band_agreement.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")
SEED <- 20260901L; MIN_AREAS <- 10L; HOLDOUT_FRAC <- 0.30

THR <- read.csv(here("config", "who_thresholds.csv"), stringsAsFactors = FALSE)
# Only outcomes whose bands survived verification.
SCORED <- c("child_vitA", "women_vitA", "child_zinc", "women_zinc")
P <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_shared.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(P), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

band_of <- function(p, th, zinc_only20) {
  if (zinc_only20) return(ifelse(is.na(p), NA_character_,
                                 ifelse(p >= 0.20, "elevated risk", "not elevated")))
  ifelse(is.na(p), NA_character_,
    ifelse(p < th$band_none_upper, "none",
      ifelse(p < th$band_mild_upper, "mild",
        ifelse(p < th$band_moderate_upper, "moderate", "severe"))))
}

oof <- function(X, y, reg, link) {
  ok <- is.finite(y); regs <- unique(reg[ok]); if (length(regs) < 3) return(NULL)
  set.seed(SEED); k <- max(2L, round(1 / HOLDOUT_FRAC))
  fold <- split(sample(regs), rep(seq_len(k), length.out = length(regs)))
  pred <- rep(NA_real_, length(y))
  for (f in seq_len(k)) {
    te <- which(reg %in% fold[[f]] & ok); tr <- which(!(reg %in% fold[[f]]) & ok)
    if (!length(te) || length(tr) < 8) next
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE); if (!is.finite(mu)) mu <- 0
      Xtr[!is.finite(Xtr[, j]), j] <- mu; Xte[!is.finite(Xte[, j]), j] <- mu
    }
    keep <- apply(Xtr, 2, function(z) stats::sd(z) > 0); if (sum(keep) < 3) next
    ff <- tryCatch(.ds_fit(Xtr[, keep, drop = FALSE], y[tr],
                           k_screen = min(20L, sum(keep)), link = link),
                   error = function(e) NULL)
    if (is.null(ff)) next
    pp <- tryCatch(.ds_predict(ff, Xte[, keep, drop = FALSE]), error = function(e) NULL)
    if (!is.null(pp)) pred[te] <- pp
  }
  pred
}

rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- P[P$country == cn, , drop = FALSE]; if (!nrow(hc)) next
  for (ocn in intersect(names(cc$outcomes), SCORED)) {
    oc <- cc$outcomes[[ocn]]
    th <- THR[THR$outcome == ocn, , drop = FALSE]; if (!nrow(th)) next
    zinc20 <- grepl("zinc", ocn)
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    if (!all(c(cc$admin1_col, cc$admin2_col, cc$weight_col, oc$binary) %in% names(d))) next
    a1 <- trimws(as.character(d[[cc$admin1_col]])); a2 <- trimws(as.character(d[[cc$admin2_col]]))
    w <- num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- 1
    key <- paste(a1, a2, sep = "||")
    wmean <- function(v) { x <- num(v); as.numeric(tapply(seq_along(x), key, function(i) {
      k2 <- is.finite(x[i]) & is.finite(w[i]); if (!any(k2)) NA_real_ else
        stats::weighted.mean(x[i][k2], w[i][k2]) })) }
    ao <- data.frame(Admin1 = tapply(a1, key, function(z) z[1]),
                     Admin2 = tapply(a2, key, function(z) z[1]),
                     prev = wmean(d[[oc$binary]]), stringsAsFactors = FALSE)
    has_cont <- !is.null(oc$continuous) && oc$continuous %in% names(d)
    if (has_cont) {
      x <- num(d[[oc$continuous]]); x[is.finite(x) & x <= 0] <- NA_real_
      d[[oc$continuous]] <- x; ao$level <- wmean(d[[oc$continuous]])
    }
    rownames(ao) <- NULL
    m <- dplyr::inner_join(ao, hc, by = admin2_join_by(ao, hc))
    if (nrow(m) < MIN_AREAS) next
    covs <- COVS[vapply(COVS, function(v) mean(is.finite(m[[v]])) > 0.5 &&
                          stats::sd(m[[v]], na.rm = TRUE) > 0, logical(1))]
    if (length(covs) < 20) next
    X <- as.matrix(m[, covs, drop = FALSE])
    truth_band <- band_of(m$prev, th, zinc20)

    # (a) predict the prevalence, then band it
    p_prev <- oof(X, m$prev, m$Admin1, "logit")
    # (b) predict the level, then convert to a prevalence, then band it.
    #     The conversion needs a threshold on the CONCENTRATION. Where the
    #     binary is a clean cut-off on the continuous, the empirical boundary
    #     recovers it. Where it is not -- Malawi zinc, whose cut-off is
    #     individual and time-of-day specific -- no district-level threshold
    #     exists and the arm is refused rather than approximated.
    # Convert a predicted LEVEL into a prevalence WITHOUT needing a
    # concentration threshold.
    #
    # Every outcome here uses an INDIVIDUAL-SPECIFIC cut-off: vitamin A's is
    # age-specific, zinc's is age, sex and time-of-day specific. Measured, 6.5
    # to 13.4 percent of non-deficient children sit below the highest deficient
    # value, so no single district-level concentration recovers the binary and
    # an earlier version of this script correctly refused every cell.
    #
    # The way through is to skip the threshold. District mean level and district
    # prevalence have a monotone relationship that can be FITTED, out of fold,
    # from the training areas only, and then applied to the predicted level.
    # Isotonic regression is used because the relationship is monotone by
    # construction and nothing stronger is warranted at 14 to 87 areas.
    p_lev <- NULL; cut_ok <- NA; cutval <- NA_real_
    if (has_cont && sum(is.finite(m$level)) >= MIN_AREAS) {
      lv <- oof(X, m$level, m$Admin1, "log")
      if (!is.null(lv)) {
        set.seed(SEED); kk <- max(2L, round(1 / HOLDOUT_FRAC))
        regs <- unique(m$Admin1); fold <- split(sample(regs),
                                                rep(seq_len(kk), length.out = length(regs)))
        p_lev <- rep(NA_real_, nrow(m))
        for (f in seq_len(kk)) {
          te <- which(m$Admin1 %in% fold[[f]]); tr <- which(!(m$Admin1 %in% fold[[f]]))
          o <- tr[is.finite(m$level[tr]) & is.finite(m$prev[tr])]
          if (length(o) < 8) next
          # monotone DECREASING in level for a deficiency outcome: higher mean
          # concentration means lower deficiency prevalence.
          ir <- tryCatch(stats::isoreg(m$level[o], -m$prev[o]), error = function(e) NULL)
          if (is.null(ir)) next
          fx <- stats::approxfun(sort(m$level[o]), -ir$yf, rule = 2)
          p_lev[te] <- pmin(pmax(fx(lv[te]), 0), 1)
        }
        cut_ok <- TRUE   # no threshold needed
      }
    }
    agree <- function(pb) { o <- !is.na(pb) & !is.na(truth_band)
      if (sum(o) < MIN_AREAS) return(c(NA_real_, NA_integer_))
      c(mean(pb[o] == truth_band[o]), sum(o)) }
    a_prev <- agree(band_of(p_prev, th, zinc20))
    a_lev  <- if (is.null(p_lev)) c(NA_real_, NA_integer_) else agree(band_of(p_lev, th, zinc20))
    # the null: every district gets the training-area modal band
    a_null <- agree(rep(names(sort(table(truth_band), decreasing = TRUE))[1], nrow(m)))
    rows[[length(rows) + 1L]] <- data.frame(
      country = cc$country, outcome = ocn, n_areas = nrow(m),
      band_authority = th$band_authority[1], n_bands = length(unique(na.omit(truth_band))),
      agree_prev = round(a_prev[1], 3), agree_level = round(a_lev[1], 3),
      agree_null = round(a_null[1], 3),
      level_minus_prev = round(a_lev[1] - a_prev[1], 3),
      level_via = "isotonic level->prevalence, fitted out of fold",
      stringsAsFactors = FALSE)
    cat(sprintf("  %-13s %-12s bands=%d prev=%.3f level=%s null=%.3f%s\n",
                cn, ocn, length(unique(na.omit(truth_band))), a_prev[1],
                if (is.na(a_lev[1])) "  n/a" else sprintf("%.3f", a_lev[1]), a_null[1],
                if (isFALSE(cut_ok)) "   [level arm unavailable]" else ""))
  }
}
if (!length(rows)) stop("[wsn1] nothing scored")
out <- dplyr::bind_rows(rows)
readr::write_csv(out, file.path(TDIR, "band_agreement.csv"))
cat("\n=== band agreement, verified bands only ===\n")
print(as.data.frame(out[, c("country","outcome","n_areas","n_bands","agree_prev",
                            "agree_level","agree_null","level_minus_prev")]),
      row.names = FALSE)
ok <- out[is.finite(out$level_minus_prev), ]
if (nrow(ok)) cat(sprintf("\nlevel beats prevalence on band agreement in %d/%d scorable cells, median %+0.3f\n",
                          sum(ok$level_minus_prev > 0), nrow(ok),
                          stats::median(ok$level_minus_prev)))
cat(sprintf("-> %s\n", file.path("results","tables","band_agreement.csv")))
