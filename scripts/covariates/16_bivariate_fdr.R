# =============================================================================
# scripts/covariates/16_bivariate_fdr.R
#
# IS THERE *ANY* SIGNAL? BIVARIATE ASSOCIATION WITH FDR CONTROL.
#
# The area-level models barely beat a national-mean null and are matched by a
# covariate-free spatial smoother. That invites a blunt question: are any of
# these predictors associated with the biomarker outcomes at all, or is the
# whole covariate set noise?
#
# This script answers it the simplest defensible way -- one bivariate test per
# (predictor x country x outcome), then Benjamini-Hochberg FDR control -- and
# reports what survives. It deliberately does NOT cross-validate: the question
# here is association, not prediction. A predictor can be genuinely associated
# and still useless out of sample, and separating those two things is the point.
#
# TWO CORRECTIONS THAT CHANGE THE ANSWER
# --------------------------------------
# 1. SPATIAL AUTOCORRELATION. Districts are not independent. A naive Pearson
#    p-value on 87 spatially autocorrelated districts is anticonservative --
#    two smooth maps correlate "significantly" by construction. So alongside
#    the naive test this reports a RESTRICTED PERMUTATION null that shuffles the
#    outcome only WITHIN Admin-1 regions, holding the regional pattern fixed and
#    asking what the covariate adds on top of it. The gap between the two tests
#    is itself the finding.
# 2. THE FAMILY MATTERS. FDR over one outcome's 294 predictors answers a
#    different question from FDR over all 24 outcomes x 294. Both are reported.
#
# Also fits a penalised regression (elastic net, CV-selected lambda) per cell
# and records which predictors keep non-zero coefficients -- the multivariable
# counterpart to the bivariate screen.
#
#   Rscript scripts/covariates/16_bivariate_fdr.R
# -> results/tables/bivariate_fdr.csv        one row per predictor x cell
# -> results/tables/bivariate_fdr_summary.csv
# -> results/tables/penalized_retained.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")
N_PERM <- 999L
set.seed(20260830L)

cfgs <- get_country_configs()
rows <- list(); pen_rows <- list()

for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  cov <- tryCatch(targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE),
                  error = function(e) NULL)
  if (is.list(cov) && !is.data.frame(cov) && "gee_admin2" %in% names(cov))
    cov <- cov$gee_admin2
  if (is.null(cov)) next

  for (on in names(cc$outcomes)) {
    svy <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on), store = STORE),
                    error = function(e) NULL)
    if (is.null(svy) || !nrow(svy)) next
    by <- admin2_join_by(svy, cov)
    d <- dplyr::inner_join(svy[, c(by, "svy_prev", "n_svy")], cov, by = by)
    d <- d[is.finite(d$svy_prev), , drop = FALSE]
    if (nrow(d) < 12) next

    covs <- setdiff(names(d), c("Admin1", "Admin2", "svy_prev", "n_svy"))
    covs <- covs[vapply(covs, function(v) is.numeric(d[[v]]), logical(1))]
    X <- as.matrix(d[, covs, drop = FALSE]); X[!is.finite(X)] <- NA
    for (j in seq_len(ncol(X))) { m <- stats::median(X[, j], na.rm = TRUE)
      X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0 }
    X <- X[, apply(X, 2, function(z) stats::sd(z) > 0), drop = FALSE]
    if (ncol(X) < 2) next
    y <- stats::qlogis(pmin(pmax(d$svy_prev, .002), .998))

    r_obs <- suppressWarnings(as.numeric(stats::cor(X, y)))
    n <- nrow(X)
    # naive two-sided t test on the correlation
    tstat <- r_obs * sqrt((n - 2) / pmax(1e-12, 1 - r_obs^2))
    p_naive <- 2 * stats::pt(-abs(tstat), df = n - 2)

    # RESTRICTED PERMUTATION, SHUFFLING WITHIN ADMIN-1 BLOCKS.
    #
    # Districts are not independent: two smooth maps correlate "significantly"
    # on a naive t-test by construction, so that p-value is anticonservative.
    # Shuffling the outcome only WITHIN each Admin-1 region holds the regional
    # mean structure fixed -- the coarse spatial pattern a covariate-free
    # smoother already captures -- and asks whether the covariate tracks
    # district-level deviations on top of it. That is the signal the covariates
    # would have to supply to be worth their place, and it is the conservative
    # test. One line, so there is nothing subtle to get wrong.
    blk <- if ("Admin1" %in% names(d)) as.character(d$Admin1) else rep("all", n)
    if (dplyr::n_distinct(blk) >= 2 && max(table(blk)) >= 2) {
      cnt <- integer(length(r_obs))
      for (b in seq_len(N_PERM)) {
        yb <- stats::ave(y, blk, FUN = function(z) if (length(z) > 1) sample(z) else z)
        rp <- abs(suppressWarnings(as.numeric(stats::cor(X, yb))))
        cnt <- cnt + (rp >= abs(r_obs))
      }
      p_block <- (cnt + 1) / (N_PERM + 1)
    } else p_block <- rep(NA_real_, length(r_obs))

    rows[[length(rows) + 1L]] <- data.frame(
      country = cn, outcome = on, n_areas = n, variable = colnames(X),
      r = round(r_obs, 4), p_naive = p_naive, p_block = p_block,
      stringsAsFactors = FALSE)

    # ── penalised regression: what survives a multivariable fit ─────────────
    if (requireNamespace("glmnet", quietly = TRUE) && n >= 15) {
      Xs <- scale(X)
      Xs <- Xs[, apply(Xs, 2, function(z) all(is.finite(z))), drop = FALSE]
      fit <- tryCatch(glmnet::cv.glmnet(Xs, y, alpha = 0.5,
                        nfolds = min(10, n), family = "gaussian"),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        for (lam in c("lambda.min", "lambda.1se")) {
          co <- as.matrix(stats::coef(fit, s = fit[[lam]]))
          nz <- rownames(co)[co[, 1] != 0]; nz <- setdiff(nz, "(Intercept)")
          pen_rows[[length(pen_rows) + 1L]] <- data.frame(
            country = cn, outcome = on, n_areas = n, lambda = lam,
            n_retained = length(nz),
            dev_ratio = round(fit$glmnet.fit$dev.ratio[which(fit$lambda == fit[[lam]])][1], 4),
            retained = paste(utils::head(nz[order(-abs(co[nz, 1]))], 12), collapse = "|"),
            stringsAsFactors = FALSE)
        }
      }
    }
    cat(sprintf("  %-12s %-13s n=%3d p=%3d | max|r|=%.2f | naive p<.05: %3d | block p<.05: %3d\n",
                cn, on, n, ncol(X), max(abs(r_obs)), sum(p_naive < .05),
                sum(p_block < .05, na.rm = TRUE)))
  }
}

res <- dplyr::bind_rows(rows)

# ── FDR, two families ───────────────────────────────────────────────────────
res <- res %>% group_by(country, outcome) %>%
  mutate(q_within_naive = stats::p.adjust(p_naive, "BH"),
         q_within_block = stats::p.adjust(p_block, "BH")) %>% ungroup() %>%
  mutate(q_global_naive = stats::p.adjust(p_naive, "BH"),
         q_global_block = stats::p.adjust(p_block, "BH"))

readr::write_csv(res, here("results", "tables", "bivariate_fdr.csv"))
pen <- dplyr::bind_rows(pen_rows)
readr::write_csv(pen, here("results", "tables", "penalized_retained.csv"))

cat("\n================ BIVARIATE ASSOCIATION, WITH FDR =====================\n")
cat(sprintf("tests: %d predictor x cell across %d cells\n\n",
            nrow(res), dplyr::n_distinct(res$country, res$outcome)))
tab <- data.frame(
  test = c("naive p < 0.05", "naive, BH within cell q<0.05", "naive, BH across all q<0.05",
           "block-permutation p < 0.05", "block, BH within cell q<0.05", "block, BH across all q<0.05"),
  n = c(sum(res$p_naive < .05), sum(res$q_within_naive < .05), sum(res$q_global_naive < .05),
        sum(res$p_block < .05, na.rm = TRUE), sum(res$q_within_block < .05, na.rm = TRUE),
        sum(res$q_global_block < .05, na.rm = TRUE)))
tab$pct <- sprintf("%.2f%%", 100 * tab$n / nrow(res))
print(tab, row.names = FALSE)
cat(sprintf("\nexpected by chance at p<0.05: %.0f\n", 0.05 * nrow(res)))

cat("\n--- cells with ANY survivor after within-cell BH (block null) ---\n")
sv <- res %>% filter(q_within_block < .05) %>% count(country, outcome, name = "n_sig")
print(if (nrow(sv)) as.data.frame(sv) else "none", row.names = FALSE)

cat("\n--- strongest 15 associations (block permutation) ---\n")
print(res %>% arrange(p_block, -abs(r)) %>%
        select(country, outcome, variable, r, p_block, q_within_block) %>%
        utils::head(15) %>% as.data.frame(), row.names = FALSE)

cat("\n================ PENALISED REGRESSION: WHAT IS RETAINED ==============\n")
print(pen %>% group_by(lambda) %>%
        summarise(cells = n(), median_retained = stats::median(n_retained),
                  mean_retained = round(mean(n_retained), 1),
                  cells_with_zero = sum(n_retained == 0),
                  mean_dev_ratio = round(mean(dev_ratio, na.rm = TRUE), 3),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)
cat("\n--- lambda.1se, per cell ---\n")
print(pen %>% filter(lambda == "lambda.1se") %>%
        select(country, outcome, n_retained, dev_ratio, retained) %>%
        mutate(retained = substr(retained, 1, 70)) %>% as.data.frame(), row.names = FALSE)
cat("\n-> results/tables/bivariate_fdr.csv\n-> results/tables/penalized_retained.csv\n")
