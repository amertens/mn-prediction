# =============================================================================
# shap_kernel/shap_kernel_lib.R  — Option B: model-agnostic SHAP for the SL
#
# Computes per-district SHAP "top factors" for the ACTUAL SuperLearner
# predictions (including the BART winner) — not a surrogate. Self-contained:
# a standard Monte-Carlo / sampling SHAP estimator (Strumbelj & Kononenko 2014),
# vectorised over instances, using only base R + the SL's own predict path.
# No fastshap/kernelshap dependency.
#
# Key design choices that make this correct and tractable:
#   - Predict via predict(fit$sl_fit$mlr3_fit, newdata) — the working path used
#     by the ablation modules (the wrapper's $predict has a known cluster_id bug).
#   - Explain the SAVED, model-ready train_data (already imputed + recipe-baked),
#     so we never have to reconstruct the ck37r imputation that the recipe omits.
#   - Recover each row's district by POSITION alignment with the non-missing-
#     outcome rows of outcome_data$data (with a hard length check).
#   - Explain only the top-K most important features (from single_var_importance)
#     and a per-district sample of individuals, to bound cost. These caps are
#     logged, not silent.
#
# NOTE: heavy. Intended to run unattended/checkpointed via run_shap_kernel.R.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
PROD_STORE <- Sys.getenv("PROD_STORE", "_targets_full")
.read_obj <- function(nm) {
  p <- file.path(PROD_STORE, "objects", nm)
  if (file.exists(p)) tryCatch(readRDS(p), error = function(e) NULL) else NULL
}

# Top-K features by permutation importance (falls back to first K covars).
shap_top_features <- function(country_label, oc_tag, covars, top_k) {
  vp <- file.path("results", "tables", "single_var_importance.csv")
  if (file.exists(vp)) {
    vi <- utils::read.csv(vp, stringsAsFactors = FALSE)
    vi <- vi[vi$country == country_label & vi$outcome == oc_tag &
               vi$variable %in% covars, , drop = FALSE]
    if (nrow(vi) > 0) {
      vi <- vi[order(-vi$drop_mean), ]
      return(utils::head(vi$variable, top_k))
    }
  }
  utils::head(covars, top_k)
}

# Monte-Carlo sampling SHAP for features `feats` over instances X (n x covars),
# background bg (m x covars), prediction fn pfun(covars_df) -> numeric length n.
mc_shap <- function(X, bg, covars, feats, pfun, nsim = 20L, seed = 1L) {
  set.seed(seed)
  n <- nrow(X)
  phi <- matrix(0, n, length(feats), dimnames = list(NULL, feats))
  for (j in feats) {
    acc <- numeric(n)
    for (s in seq_len(nsim)) {
      z <- bg[sample(nrow(bg), n, replace = TRUE), covars, drop = FALSE]
      perm <- sample(covars)
      jpos <- match(j, perm)
      S <- if (jpos > 1) perm[seq_len(jpos - 1)] else character(0)
      xw <- z                              # rest of features from reference
      keep <- c(S, j)
      xw[, keep] <- X[, keep, drop = FALSE]  # S and j from the instance
      xo <- xw
      xo[, j] <- z[, j]                      # remove j (take from reference)
      acc <- acc + (pfun(xw) - pfun(xo))
    }
    phi[, j] <- acc / nsim
  }
  phi
}

shap_aggregate_district <- function(phi, a2, a1, country_label, oc_tag, n_top) {
  rows <- list()
  for (dd in unique(a2)) {
    ii <- which(a2 == dd)
    m  <- colMeans(phi[ii, , drop = FALSE])
    o  <- order(abs(m), decreasing = TRUE)
    tk <- o[seq_len(min(n_top, length(o)))]
    rows[[dd]] <- data.frame(
      country = country_label, outcome = oc_tag,
      Admin1 = a1[ii][1], Admin2 = dd, rank = seq_along(tk),
      variable = names(m)[tk], shap = as.numeric(m[tk]),
      direction = ifelse(m[tk] > 0, "increases risk", "decreases risk"),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# One country x outcome slice. Returns list(district_factors, global_importance,
# method, ...) or list(error=...) so the driver can log and continue.
compute_shap_slice <- function(low, cc, oc,
                               n_top = 5L, top_k = 20L,
                               n_per_district = 8L, n_explain_max = 400L,
                               bg_size = 100L, nsim = 20L, seed = 1L,
                               smoke = FALSE) {
  if (smoke) { top_k <- 5L; n_per_district <- 1L; n_explain_max <- 30L
               bg_size <- 20L; nsim <- 2L }
  oc_tag <- oc$tag %||% oc$label
  country_label <- cc$country
  admin2_col <- cc$admin2_col %||% "Admin2"
  admin1_col <- cc$admin1_col %||% "Admin1"

  sf <- .read_obj(paste0("sl_fit_", low, "_", oc_tag))
  od <- .read_obj(paste0("outcome_data_", low, "_", oc_tag))
  if (is.null(sf) || is.null(od)) return(list(error = "missing sl_fit/outcome_data"))

  use_binary <- !is.null(sf$bin_fit)
  fit <- if (use_binary) sf$bin_fit else sf$cont_fit
  mlr3_fit <- tryCatch(fit$sl_fit$mlr3_fit, error = function(e) NULL)
  if (is.null(mlr3_fit) || is.null(fit$train_data)) return(list(error = "no mlr3_fit/train_data"))

  tmpl <- as.data.frame(fit$train_data)
  covars <- intersect(fit$Xvars %||% fit$sl_fit$covars, colnames(tmpl))
  if (length(covars) < 2) return(list(error = "too few covars"))

  ycol <- if (use_binary) oc$binary else oc$continuous
  d <- od$data
  if (is.null(ycol) || !ycol %in% colnames(d)) return(list(error = "outcome col missing"))
  valid <- !is.na(d[[ycol]])
  if (sum(valid) != nrow(tmpl))
    return(list(error = sprintf("row alignment mismatch (valid=%d, train=%d)",
                                sum(valid), nrow(tmpl))))
  a2 <- as.character(d[[admin2_col]][valid])
  a1 <- if (admin1_col %in% colnames(d)) as.character(d[[admin1_col]][valid])
        else rep(NA_character_, sum(valid))

  Xfull <- tmpl[, covars, drop = FALSE]
  keep_row <- stats::complete.cases(Xfull) & !is.na(a2)
  Xfull <- Xfull[keep_row, , drop = FALSE]; a2 <- a2[keep_row]; a1 <- a1[keep_row]
  if (nrow(Xfull) < 5) return(list(error = "too few usable rows"))

  # explain rows: up to n_per_district per district, capped at n_explain_max
  set.seed(seed)
  idx <- unlist(lapply(split(seq_along(a2), a2), function(ii)
    if (length(ii) <= n_per_district) ii else sample(ii, n_per_district)))
  if (length(idx) > n_explain_max) idx <- sort(sample(idx, n_explain_max))
  Xe <- Xfull[idx, , drop = FALSE]; a2e <- a2[idx]; a1e <- a1[idx]
  bg <- Xfull[sample(nrow(Xfull), min(bg_size, nrow(Xfull))), , drop = FALSE]

  feats <- intersect(shap_top_features(country_label, oc_tag, covars, top_k), covars)
  if (length(feats) == 0) return(list(error = "no features selected"))

  # predict closure: a cached n-row template carrying cluster_id/Y; overwrite covars
  n <- nrow(Xe)
  tmpl_n <- tmpl[rep(1L, n), , drop = FALSE]
  pfun <- function(Cmat) {
    tmpl_n[, covars] <- Cmat[, covars]
    as.numeric(stats::predict(mlr3_fit, tmpl_n))
  }

  phi <- mc_shap(Xe, bg, covars, feats, pfun, nsim = nsim, seed = seed)

  district_factors <- shap_aggregate_district(phi, a2e, a1e, country_label, oc_tag, n_top)
  global <- data.frame(country = country_label, outcome = oc_tag,
                       variable = colnames(phi),
                       mean_abs_shap = colMeans(abs(phi)),
                       mean_shap = colMeans(phi), stringsAsFactors = FALSE)
  global <- global[order(-global$mean_abs_shap), ]; global$rank <- seq_len(nrow(global))

  list(district_factors = district_factors, global_importance = global,
       method = "mc_sampling_shap", outcome_type = if (use_binary) "binary" else "continuous",
       n_explained = nrow(Xe), n_features = length(feats), nsim = nsim,
       n_districts = length(unique(a2e)))
}
