# =============================================================================
# scripts/accuracy_impact/wsm1_variable_importance.R
#
# WS-M. Variable importance, reported BY DOMAIN rather than by variable.
#
# WHY NOT A RANKED LIST OF VARIABLES
# ----------------------------------
# A single ranked importance list over ~370 predictors on 14 to 87 areas would
# look authoritative and would not replicate. The effective sample size is the
# number of areas, the predictors are heavily collinear within domain, and this
# project has already measured what happens when many predictors meet few
# areas: 2 survivors from 294, and 0 from 31 chosen in advance on mechanism.
#
# So importance is measured three ways, none of which is a per-variable ranking:
#
#   DOMAIN DROP    refit with one whole domain removed and measure the loss.
#                  This answers the question a programme actually asks -- is it
#                  worth buying soil data, or DHS data -- and is stable because
#                  a domain has many correlated members.
#   DOMAIN ALONE   refit using ONLY that domain. Drop and alone disagree when a
#                  domain is redundant with others: informative alone, free to
#                  remove.
#   SELECTION RATE how often a domain's members survive the in-fold screen,
#                  across folds and cells. A frequency, not a coefficient.
#
# Everything is scored under the consolidated protocol: region-blocked folds,
# training-area-mean null, the shared predictor set with its domain metadata.
#
#   Rscript scripts/accuracy_impact/wsm1_variable_importance.R
# -> results/tables/variable_importance_by_domain.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")
SEED <- 20260930L
MIN_AREAS <- 10L
HOLDOUT_FRAC <- 0.30
MIN_DOMAIN_N <- 3L      # a domain with fewer members cannot be dropped meaningfully

SHARED <- here("data", "covariates", "harmonized", "predictors_admin2_shared.csv")
META   <- here("data", "covariates", "harmonized", "predictors_admin2_shared_metadata.csv")
if (!file.exists(SHARED) || !file.exists(META))
  stop("[wsm1] shared predictor set or metadata not found; run build_shared_predictor_set.R")
P <- suppressMessages(readr::read_csv(SHARED, show_col_types = FALSE)) |> as.data.frame()
M <- suppressMessages(readr::read_csv(META, show_col_types = FALSE)) |> as.data.frame()
cfgs <- get_country_configs()
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# Only predictors that vary WITHIN a country can move a district ranking, so a
# national constant is excluded here rather than credited to its domain.
M <- M[M$subnational %in% c(TRUE, "TRUE"), , drop = FALSE]
DOMAINS <- sort(unique(M$domain))
cat(sprintf("[wsm1] %d subnational predictors across %d domains\n",
            nrow(M), length(DOMAINS)))

fit_predict <- function(X, y, n_y, reg, link) {
  ok <- is.finite(y)
  regs <- unique(reg[ok]); if (length(regs) < 3) return(NULL)
  set.seed(SEED)
  k <- max(2L, round(1 / HOLDOUT_FRAC))
  fold <- split(sample(regs), rep(seq_len(k), length.out = length(regs)))
  pred <- rep(NA_real_, length(y)); sel_all <- character(0)
  for (f in seq_len(k)) {
    te <- which(reg %in% fold[[f]] & ok); tr <- which(!(reg %in% fold[[f]]) & ok)
    if (!length(te) || length(tr) < 8) next
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE); if (!is.finite(mu)) mu <- 0
      Xtr[!is.finite(Xtr[, j]), j] <- mu; Xte[!is.finite(Xte[, j]), j] <- mu
    }
    keep <- apply(Xtr, 2, function(z) stats::sd(z) > 0)
    if (sum(keep) < 3) next
    ff <- tryCatch(.ds_fit(Xtr[, keep, drop = FALSE], y[tr],
                           k_screen = min(20L, sum(keep)), link = link),
                   error = function(e) NULL)
    if (is.null(ff)) next
    sel_all <- c(sel_all, ff$sel)
    pp <- tryCatch(.ds_predict(ff, Xte[, keep, drop = FALSE]), error = function(e) NULL)
    if (!is.null(pp)) pred[te] <- pp
  }
  o <- is.finite(y) & is.finite(pred)
  if (sum(o) < MIN_AREAS) return(NULL)
  list(mae = mean(abs(pred[o] - y[o])), n = sum(o), selected = sel_all)
}

rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- P[P$country == cn, , drop = FALSE]; if (!nrow(hc)) next
  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od) || is.null(oc$binary)) next
    d <- od$data
    if (!all(c(cc$admin1_col, cc$admin2_col, cc$weight_col, oc$binary) %in% names(d))) next
    a1 <- trimws(as.character(d[[cc$admin1_col]])); a2 <- trimws(as.character(d[[cc$admin2_col]]))
    w <- num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- 1
    key <- paste(a1, a2, sep = "||"); yb <- num(d[[oc$binary]])
    ao <- data.frame(Admin1 = tapply(a1, key, function(z) z[1]),
                     Admin2 = tapply(a2, key, function(z) z[1]),
                     prev = as.numeric(tapply(seq_along(yb), key, function(i) {
                       k2 <- is.finite(yb[i]) & is.finite(w[i])
                       if (!any(k2)) NA_real_ else stats::weighted.mean(yb[i][k2], w[i][k2]) })),
                     n_prev = as.integer(tapply(seq_along(yb), key,
                                                function(i) sum(is.finite(yb[i])))),
                     stringsAsFactors = FALSE)
    rownames(ao) <- NULL
    m <- dplyr::inner_join(ao, hc, by = admin2_join_by(ao, hc))
    if (nrow(m) < MIN_AREAS) next
    usable <- intersect(M$column, names(m))
    usable <- usable[vapply(usable, function(v)
      mean(is.finite(m[[v]])) > 0.5 && stats::sd(m[[v]], na.rm = TRUE) > 0, logical(1))]
    if (length(usable) < 20) next
    X <- as.matrix(m[, usable, drop = FALSE])
    base <- fit_predict(X, m$prev, m$n_prev, m$Admin1, "logit")
    if (is.null(base)) next
    cat(sprintf("  %-13s %-13s %d areas, %d predictors, base MAE %.4f\n",
                cn, ocn, nrow(m), length(usable), base$mae))
    dm <- M$domain[match(usable, M$column)]
    for (dom in unique(dm)) {
      idx <- which(dm == dom)
      if (length(idx) < MIN_DOMAIN_N) next
      drop_r  <- fit_predict(X[, -idx, drop = FALSE], m$prev, m$n_prev, m$Admin1, "logit")
      alone_r <- fit_predict(X[,  idx, drop = FALSE], m$prev, m$n_prev, m$Admin1, "logit")
      rows[[length(rows) + 1L]] <- data.frame(
        country = cc$country, outcome = ocn, domain = dom,
        n_in_domain = length(idx), n_areas = base$n,
        base_mae = round(base$mae, 5),
        drop_mae = if (is.null(drop_r)) NA_real_ else round(drop_r$mae, 5),
        # Positive means removing the domain HURT, i.e. it was carrying signal.
        drop_cost = if (is.null(drop_r)) NA_real_ else round(drop_r$mae - base$mae, 5),
        alone_mae = if (is.null(alone_r)) NA_real_ else round(alone_r$mae, 5),
        sel_rate = round(mean(base$selected %in% usable[idx]), 4),
        stringsAsFactors = FALSE)
    }
  }
}

if (!length(rows)) stop("[wsm1] nothing scored")
out <- dplyr::bind_rows(rows)
readr::write_csv(out, file.path(TDIR, "variable_importance_by_domain.csv"))

cat("\n=== domain importance, pooled over cells ===\n")
s <- out |> group_by(domain) |>
  summarise(cells = dplyr::n(), n_vars = max(n_in_domain),
            mean_drop_cost = round(mean(drop_cost, na.rm = TRUE), 5),
            hurts_to_drop = sprintf("%d/%d", sum(drop_cost > 0, na.rm = TRUE), dplyr::n()),
            mean_alone_mae = round(mean(alone_mae, na.rm = TRUE), 5),
            mean_sel_rate = round(mean(sel_rate, na.rm = TRUE), 3),
            .groups = "drop") |> arrange(desc(mean_drop_cost))
print(as.data.frame(s), row.names = FALSE)
cat(sprintf("\n-> %s\n", file.path("results","tables","variable_importance_by_domain.csv")))
