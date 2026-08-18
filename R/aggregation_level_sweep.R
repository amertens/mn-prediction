# =============================================================================
# R/aggregation_level_sweep.R
#
# Does the AGGREGATION LEVEL at which we fit change NATIONAL-PREVALENCE accuracy?
# For each country x outcome and each LEVEL in {cluster, admin1, admin2}: fit a
# spatial-block out-of-fold ridge on level-aggregated proxies, roll the OOF unit
# predictions up to a national estimate, and compare to the observed national
# prevalence. Shared by:
#   - the `aggregation_level_national` target in _targets.R (refreshes w/ pipeline)
#   - sensitivity/12_aggregation_level_national.R (ad-hoc runs)
#
# See the header of sensitivity/12 for the full rationale and caveats (cluster
# rows share ~admin-2-resolution predictors; ridge not the full SuperLearner —
# consistent across levels, which is what makes the comparison fair).
# =============================================================================

#' Aggregate individual data to a level: prevalence, weight, block, mean predictors.
.agg_to_level <- function(d, level_col, block_col, ycol, wcol, Xvars) {
  d <- d[!is.na(d[[level_col]]) & !is.na(d[[ycol]]) & !is.na(d[[wcol]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  w <- d[[wcol]]; lev <- as.character(d[[level_col]]); blk <- as.character(d[[block_col]])
  wmean <- function(x, wt) { ok <- !is.na(x); if (!any(ok)) NA_real_ else sum(x[ok]*wt[ok])/sum(wt[ok]) }
  units <- unique(lev)
  prev  <- vapply(units, function(u) wmean(d[[ycol]][lev==u], w[lev==u]), numeric(1))
  wsum  <- vapply(units, function(u) sum(w[lev==u], na.rm = TRUE), numeric(1))
  block <- vapply(units, function(u) { b <- blk[lev==u]; b <- b[!is.na(b)]
    if (length(b)) names(sort(table(b), decreasing = TRUE))[1] else NA_character_ }, character(1))
  X <- vapply(Xvars, function(v) vapply(units, function(u) wmean(d[[v]][lev==u], w[lev==u]),
                                        numeric(1)), numeric(length(units)))
  if (length(units) == 1L) X <- matrix(X, nrow = 1, dimnames = list(NULL, Xvars))
  list(unit = units, prev = prev, w = wsum, block = block, X = X)
}

#' Spatial-block out-of-fold ridge predictions (folds = blocks).
.oof_ridge_national <- function(X, y, w, block) {
  keep <- apply(X, 2, function(c) { v <- stats::var(c, na.rm = TRUE); !is.na(v) && v > 0 })
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 1) return(rep(stats::weighted.mean(y, w, na.rm = TRUE), length(y)))
  for (j in seq_len(ncol(X))) { m <- stats::median(X[, j], na.rm = TRUE); X[is.na(X[, j]), j] <- m }

  blocks <- unique(block)
  if (length(blocks) < 3) {  # too few spatial blocks -> random 5-fold
    set.seed(12345); block <- sample(rep_len(1:min(5, length(y)), length(y)))
    blocks <- sort(unique(block))
  }
  oof <- rep(NA_real_, length(y))
  for (b in blocks) {
    tr <- which(block != b); te <- which(block == b)
    if (length(tr) < 3 || length(te) < 1) next
    fit <- tryCatch(
      glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], alpha = 0,
                        weights = w[tr], nfolds = min(5, length(tr))),
      error = function(e) NULL)
    if (is.null(fit)) next
    oof[te] <- as.numeric(stats::predict(fit, X[te, , drop = FALSE], s = "lambda.min"))
  }
  oof[is.na(oof)] <- stats::weighted.mean(y, w, na.rm = TRUE)
  pmin(pmax(oof, 0), 1)
}

#' Aggregation-level national-prevalence sweep for ONE country x outcome.
#' @param od Output of build_outcome_dataset() (list with $data, $Xvars).
#' @return data.frame (one row per level) or NULL.
aggregation_sweep_one <- function(od, cc, oc) {
  if (is.null(od) || is.null(od$data)) return(NULL)
  d <- od$data; Xvars <- intersect(od$Xvars, colnames(d))
  # Only survey-weight-average NUMERIC predictors (some proxies are character/
  # factor and would error in the weighted mean).
  Xvars <- Xvars[vapply(Xvars, function(v) is.numeric(d[[v]]), logical(1))]
  ycol <- oc$binary; wcol <- cc$weight_col
  if (is.null(ycol) || !ycol %in% colnames(d) || !wcol %in% colnames(d) || length(Xvars) < 2)
    return(NULL)

  levels <- list(cluster = cc$cluster_id, admin1 = cc$admin1_col, admin2 = cc$admin2_col)
  rows <- lapply(names(levels), function(lv) {
    level_col <- levels[[lv]]
    if (is.null(level_col) || !level_col %in% colnames(d)) return(NULL)
    a <- .agg_to_level(d, level_col, cc$admin1_col, ycol, wcol, Xvars)
    if (is.null(a) || length(a$unit) < 2) return(NULL)
    oof <- .oof_ridge_national(a$X, a$prev, a$w, a$block)
    obs_nat  <- stats::weighted.mean(a$prev, a$w, na.rm = TRUE)
    pred_nat <- stats::weighted.mean(oof,   a$w, na.rm = TRUE)
    data.frame(
      country = cc$country, outcome = oc$tag, level = lv, n_units = length(a$unit),
      obs_national_pct  = 100 * obs_nat,
      pred_national_pct = 100 * pred_nat,
      abs_error_pp      = 100 * abs(pred_nat - obs_nat),
      signed_error_pp   = 100 * (pred_nat - obs_nat),
      rank_spearman     = suppressWarnings(stats::cor(oof, a$prev, method = "spearman",
                                                      use = "complete.obs")),
      stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(rows)
}

#' Run the sweep across a named list of outcome_data objects.
#' @param od_list named list; names are "<lower(country_key)>_<outcome_tag>".
#' @param configs get_country_configs().
#' @return combined data.frame across all country x outcome x level.
run_aggregation_level_sweep <- function(od_list, configs) {
  rows <- list()
  for (ck in names(configs)) {
    cc <- configs[[ck]]
    for (tag in names(cc$outcomes)) {
      sfx <- paste0(tolower(ck), "_", tag)
      od  <- od_list[[sfx]]
      if (is.null(od)) next
      rows[[sfx]] <- tryCatch(aggregation_sweep_one(od, cc, cc$outcomes[[tag]]),
                              error = function(e) { message("[agg-sweep] ", sfx, ": ",
                                                            conditionMessage(e)); NULL })
    }
  }
  dplyr::bind_rows(rows)
}
