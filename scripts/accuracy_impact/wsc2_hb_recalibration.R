# =============================================================================
# scripts/accuracy_impact/wsc2_hb_recalibration.R
#
# WS-C2. The field-haemoglobin arm buys ranking and costs level. Recalibrate it.
#
# Measured at sixteen cells: the haemoglobin arm reaches mean r 0.327 against
# 0.154 for proxy, and district MAE 11.07 pp against 8.33. It is the strongest
# ranking instrument in the project and the worst of the three as a prevalence.
#
# Recalibration is done OUT OF FOLD under the same region-blocked protocol: for
# each held-out region, a linear map from prediction to observed is fitted on
# the OTHER regions only and applied to the held-out one. Fitting the map on all
# regions would be the same all-data-preprocessing defect Section 2.3 measures
# at +0.182.
#
#   Rscript scripts/accuracy_impact/wsc2_hb_recalibration.R
# -> results/tables/hb_arm_recalibration.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results","tables")
K_SCREEN <- 40L; SEED <- 20260925L; MIN_N <- 5L
OUTCOMES <- c("child_iron","child_vitA","women_iron","women_vitA")
set.seed(SEED)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

fit_pred <- function(Xtr, ytr, Xte) {
  sel <- .awsl_screen(Xtr, ytr, K_SCREEN)
  s <- tryCatch(.awsl_stack(Xtr[, sel, drop=FALSE], ytr, rep(1, length(ytr))),
                error=function(e) NULL)
  if (is.null(s)) return(rep(NA_real_, nrow(Xte)))
  .awsl_predict(s, Xte[, sel, drop=FALSE])
}
prep <- function(d, vars) {
  vars <- intersect(vars, names(d))
  vars <- vars[vapply(vars, function(v) is.numeric(d[[v]]) || is.logical(d[[v]]) ||
                        inherits(d[[v]], "haven_labelled"), logical(1))]
  if (!length(vars)) return(matrix(numeric(0), nrow=nrow(d), ncol=0))
  X <- vapply(vars, function(v) num(d[[v]]), numeric(nrow(d)))
  if (is.null(dim(X))) X <- matrix(X, nrow=nrow(d))
  colnames(X) <- vars
  for (j in seq_len(ncol(X))) { m <- stats::median(X[,j], na.rm=TRUE)
    X[!is.finite(X[,j]), j] <- if (is.finite(m)) m else 0 }
  X[, apply(X, 2, function(z) stats::sd(z) > 0), drop=FALSE]
}

cfgs <- get_country_configs(); rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  for (on in OUTCOMES) {
    oc <- cc$outcomes[[on]]; if (is.null(oc)) next
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_",lc,"_",on), store=STORE),
                   error=function(e) NULL)
    if (is.null(od)) next
    full <- od$Xvars_full %||% od$Xvars
    if (length(grep("^gw_", full)) == 0L) next          # Malawi: no survey block
    d <- od$data
    if (!all(c("Admin2", oc$binary) %in% names(d))) next
    y <- num(d[[oc$binary]]); keep <- is.finite(y)
    d <- d[keep,,drop=FALSE]; y <- y[keep]
    if (length(unique(y)) < 2 || nrow(d) < 100) next
    blk <- if ("Admin1" %in% names(d)) as.character(d$Admin1) else rep("a", nrow(d))
    if (dplyr::n_distinct(blk) < 3) next
    X <- prep(d, full[allowed_under_arm(full, "questionnaire_hb")])
    if (ncol(X) < 5) next

    oof <- rep(NA_real_, nrow(X)); oof_cal <- rep(NA_real_, nrow(X))
    for (f in unique(blk)) {
      i <- which(blk == f); tr <- setdiff(seq_len(nrow(X)), i)
      if (!length(tr) || length(i) >= nrow(X) || length(unique(y[tr])) < 2) next
      p_te <- fit_pred(X[tr,,drop=FALSE], y[tr], X[i,,drop=FALSE])
      oof[i] <- p_te
      # The recalibration map is fitted on the TRAINING regions only, by an
      # inner leave-one-region-out pass so the map sees out-of-fold predictions
      # rather than in-sample ones.
      inner <- rep(NA_real_, length(tr)); btr <- blk[tr]
      for (g in unique(btr)) {
        ii <- which(btr == g); tt <- setdiff(seq_along(tr), ii)
        if (length(tt) < 30 || length(unique(y[tr][tt])) < 2) next
        inner[ii] <- fit_pred(X[tr[tt],,drop=FALSE], y[tr][tt], X[tr[ii],,drop=FALSE])
      }
      okc <- is.finite(inner)
      if (sum(okc) > 30 && stats::sd(inner[okc]) > 1e-8) {
        cal <- stats::lm(y[tr][okc] ~ inner[okc])
        oof_cal[i] <- pmin(pmax(stats::coef(cal)[1] +
                                  stats::coef(cal)[2] * p_te, 0), 1)
      } else oof_cal[i] <- p_te
    }
    ok <- is.finite(oof) & is.finite(oof_cal)
    if (sum(ok) < 50) next
    key <- paste(blk, d$Admin2)
    agg <- data.frame(key=key[ok], obs=y[ok], raw=oof[ok], cal=oof_cal[ok]) |>
      group_by(key) |> summarise(n=dplyr::n(), obs=mean(obs), raw=mean(raw),
                                 cal=mean(cal), .groups="drop") |> filter(n >= MIN_N)
    if (nrow(agg) < 8) next
    rows[[length(rows)+1L]] <- data.frame(
      country=cc$country, outcome=on, n_units=nrow(agg),
      r_raw=round(stats::cor(agg$obs, agg$raw),4),
      r_cal=round(stats::cor(agg$obs, agg$cal),4),
      mae_raw=round(100*mean(abs(agg$raw-agg$obs)),2),
      mae_cal=round(100*mean(abs(agg$cal-agg$obs)),2),
      bias_raw=round(100*mean(agg$raw-agg$obs),2),
      bias_cal=round(100*mean(agg$cal-agg$obs),2),
      stringsAsFactors=FALSE)
    cat(sprintf("  %-13s %-11s r %.3f->%.3f  MAE %.2f->%.2f\n", cn, on,
                stats::cor(agg$obs,agg$raw), stats::cor(agg$obs,agg$cal),
                100*mean(abs(agg$raw-agg$obs)), 100*mean(abs(agg$cal-agg$obs))))
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows.")
readr::write_csv(res, file.path(TDIR,"hb_arm_recalibration.csv"))
cat("\n=== WS-C2: field-haemoglobin arm, before and after out-of-fold recalibration ===\n")
cat(sprintf("cells: %d\n", nrow(res)))
cat(sprintf("mean r   %.3f -> %.3f\n", mean(res$r_raw), mean(res$r_cal)))
cat(sprintf("mean MAE %.2f -> %.2f pp\n", mean(res$mae_raw), mean(res$mae_cal)))
cat(sprintf("mean |bias| %.2f -> %.2f pp\n", mean(abs(res$bias_raw)), mean(abs(res$bias_cal))))
cat(sprintf("MAE improved in %d of %d cells\n", sum(res$mae_cal < res$mae_raw), nrow(res)))
