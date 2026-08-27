# =============================================================================
# sandbox_parsimony/R/33_vmnis_standardised.R
#
# STEP 1: residualise the VMNIS outcome out-of-fold, then predict the
# STANDARDISED quantity.
#
# FINDINGS.md section 24 showed that removing survey-method effects raises the
# vitamin A ceiling from r_max 0.505 to 0.633, but that adding method terms as
# PREDICTORS barely improves the fit. The reason is that the two things are
# different: you cannot know what method a future survey will use, so method
# terms do not help you predict the number that survey will report. What they
# do is define a comparable target.
#
# So here the correction is applied to the OUTCOME, not the design matrix:
#
#   1. hold out a country
#   2. estimate method + biomarker effects on the TRAINING countries only,
#      with country fixed effects so the contrast is identified within country
#   3. subtract those effects from the training outcomes AND from the held-out
#      country's outcome. This uses the held-out survey's METHOD LABEL, which is
#      observable metadata, never its prevalence -- so it is not leakage
#   4. fit on the standardised training outcome, predict, and score against the
#      standardised held-out outcome
#
# Everything is scored on the same countries as 31/32 so the numbers are
# comparable, and the raw-outcome model is refitted here as the control.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ranger); library(glmnet)})
source("sandbox_parsimony/R/37_national_covariates.R")

# Same covariate switch as 31: NAT_COV picks wdi / panel / both, NAT_COMMON
# restricts to countries both sources cover, NAT_SCREEN sets the per-fold
# prescreen width and NAT_KNN the number of imputation donors.
COV_SOURCE <- Sys.getenv("NAT_COV", "wdi")
SCREEN_K   <- as.integer(Sys.getenv("NAT_SCREEN", "150"))
KNN_K      <- as.integer(Sys.getenv("NAT_KNN", "5"))

.logit  <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.ilogit <- stats::plogis

nat <- readRDS("sandbox_parsimony/out/vmnis_national_rep.rds")

cov <- load_national_covariates(COV_SOURCE)
V        <- cov$vars
LOG_VARS <- cov$log_vars
message(sprintf("covariate source: %s | %d covariates | prescreen to %d per fold",
                cov$source, length(V), SCREEN_K))

adj_level <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("inflam", x, ignore.case = TRUE) ~ "inflammation_adjusted",
    grepl("^none$", x, ignore.case = TRUE) ~ "unadjusted",
    x == "" | grepl("not specified", x, ignore.case = TRUE) ~ "unspecified",
    TRUE ~ "other")
}

dat <- nat |>
  mutate(iso3c = suppressWarnings(countrycode::countrycode(
           country, "country.name", "iso3c", warn = FALSE)),
         pop = as.character(Population),
         adj = adj_level(Dataadjustedfor),
         biomarker = ifelse(grepl("Retinol", as.character(Indicator)),
                            "retinol", "other")) |>
  filter(!is.na(iso3c)) |>
  inner_join(cov$df |> select(iso3c, year, all_of(V)), by = c("iso3c", "year"))

if (nzchar(Sys.getenv("NAT_COMMON"))) {
  keep_cty <- common_countries()
  dat <- dat |> filter(iso3c %in% keep_cty)
  message(sprintf("common support: restricted to %d countries", length(keep_cty)))
}
message(sprintf("VMNIS x %s joined: %d rows, %d countries",
                cov$source, nrow(dat), n_distinct(dat$iso3c)))

#' Method + biomarker offsets, identified WITHIN country.
#'
#' Returns a named vector of offsets on the logit scale relative to the
#' reference level (inflammation-adjusted retinol). Estimated from `train`
#' only. Levels the training data cannot identify get an offset of 0, which is
#' the conservative choice -- it leaves those surveys uncorrected rather than
#' correcting them with a number the data did not support.
method_offsets <- function(train) {
  train <- train |> mutate(lg = .logit(prev),
                           adj = factor(adj, levels = c("inflammation_adjusted",
                                                        "unadjusted", "unspecified", "other")),
                           biomarker = factor(biomarker, levels = c("retinol", "other")))
  # need at least one country observed under >1 method for the contrast to be
  # identified at all
  ident <- train |> group_by(iso3c) |> summarise(k = n_distinct(adj), .groups = "drop")
  if (max(ident$k, na.rm = TRUE) < 2) return(NULL)
  fit <- tryCatch(stats::lm(lg ~ factor(iso3c) + adj + biomarker, data = train),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- stats::coef(fit)
  keep <- cf[grepl("^adj|^biomarker", names(cf))]
  keep[!is.finite(keep)] <- 0
  keep
}

#' Apply the offsets to a data frame, returning standardised logit prevalence.
apply_offsets <- function(d, off) {
  lg <- .logit(d$prev)
  if (is.null(off) || !length(off)) return(lg)
  d2 <- d |> mutate(adj = factor(adj, levels = c("inflammation_adjusted",
                                                 "unadjusted", "unspecified", "other")),
                    biomarker = factor(biomarker, levels = c("retinol", "other")))
  mm <- stats::model.matrix(~ adj + biomarker, data = d2)
  common <- intersect(colnames(mm), names(off))
  if (!length(common)) return(lg)
  lg - as.numeric(mm[, common, drop = FALSE] %*% off[common])
}

build_X <- function(d) {
  X <- as.matrix(d[, V, drop = FALSE])
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- NA_real_
  keep <- colSums(!is.na(X)) > 0L &
    apply(X, 2, function(x) length(unique(x[!is.na(x)])) > 1L)
  X <- X[, keep, drop = FALSE]
  for (v in intersect(LOG_VARS, colnames(X)))
    X[, v] <- log1p(pmax(X[, v], 0))
  cbind(X, year = d$year)
}

run_panel <- function(d, label, screen_k = SCREEN_K) {
  .mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  d <- d |> group_by(iso3c, year, adj, biomarker) |>
    summarise(prev = mean(prev), across(all_of(V), .mean), .groups = "drop") |>
    as.data.frame()
  # Complete-case on the OUTCOME only; missing covariates are handled by the
  # in-fold KNN imputation below rather than by discarding the survey.
  d <- d[is.finite(d$prev), , drop = FALSE]
  if (n_distinct(d$iso3c) < 12) return(NULL)
  X <- build_X(d); ctys <- unique(d$iso3c)

  pr_raw <- pr_std <- rep(NA_real_, nrow(d))
  y_std_all <- rep(NA_real_, nrow(d))
  n_ident <- 0L
  n_used <- integer(0)

  fit_predict <- function(ytr, Xtr, Xte) {
    rf <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = ytr, Xtr),
                                  num.trees = 800, min.node.size = 5, seed = 1),
                   error = function(e) NULL)
    if (is.null(rf)) return(rep(NA_real_, nrow(Xte)))
    stats::predict(rf, data = data.frame(Xte))$predictions
  }

  for (ho in ctys) {
    te <- which(d$iso3c == ho); tr <- setdiff(seq_len(nrow(d)), te)
    if (length(tr) < 10) next

    # Indicators, KNN imputation, screening -- all fitted on training rows.
    MI <- missing_indicators(X, tr)
    Xi <- knn_impute_fold(X, tr, k = KNN_K)
    Xf <- if (is.null(MI)) Xi else cbind(Xi, MI)

    # --- control: raw outcome -------------------------------------------
    yr  <- .logit(d$prev)
    s1  <- prescreen_cols(Xf[tr, , drop = FALSE], yr[tr], screen_k)
    n_used <- c(n_used, length(s1))
    pr_raw[te] <- fit_predict(yr[tr], Xf[tr, s1, drop = FALSE],
                              Xf[te, s1, drop = FALSE])

    # --- standardised: offsets from TRAINING countries only --------------
    off <- method_offsets(d[tr, , drop = FALSE])
    if (!is.null(off)) n_ident <- n_ident + 1L
    ys_tr <- apply_offsets(d[tr, , drop = FALSE], off)
    ys_te <- apply_offsets(d[te, , drop = FALSE], off)   # metadata only
    y_std_all[te] <- ys_te
    # screened against the standardised target, so each model gets the
    # covariates that matter for the quantity it is actually predicting
    s2 <- prescreen_cols(Xf[tr, , drop = FALSE], ys_tr, screen_k)
    pr_std[te] <- fit_predict(ys_tr, Xf[tr, s2, drop = FALSE],
                              Xf[te, s2, drop = FALSE])
  }

  score <- function(pred_logit, truth_logit, tag) {
    ok <- is.finite(pred_logit) & is.finite(truth_logit)
    if (sum(ok) < 10) return(NULL)
    p <- .ilogit(pred_logit[ok]); o <- .ilogit(truth_logit[ok])
    data.frame(panel = label, source = cov$source, target = tag, n = sum(ok),
               n_countries = n_distinct(d$iso3c[ok]),
               n_cov_pool = ncol(X),
               n_cov_used = if (length(n_used)) round(mean(n_used)) else NA_integer_,
               mae_pp = round(100 * mean(abs(p - o)), 2),
               pearson = round(suppressWarnings(cor(p, o)), 3),
               spearman = round(suppressWarnings(cor(p, o, method = "spearman")), 3),
               folds_identified = n_ident, stringsAsFactors = FALSE)
  }
  bind_rows(
    score(pr_raw, .logit(d$prev), "raw outcome (control)"),
    score(pr_std, y_std_all,      "standardised outcome"))
}

PANELS <- list(c("Vitamin A", "Preschool-age children"),
               c("Zinc", "Preschool-age children"),
               c("Folate", "Non-pregnant women (NPW)"))

res <- list()
for (pp in PANELS) {
  lab <- paste(pp, collapse = " / ")
  d <- dat |> filter(mn_group == pp[1], pop == pp[2])
  r <- tryCatch(run_panel(d, lab), error = function(e) {message("  ", lab, ": ", conditionMessage(e)); NULL})
  if (!is.null(r)) res[[lab]] <- r
}
out <- bind_rows(res)
OUT_CSV <- sprintf("sandbox_parsimony/out/vmnis_standardised_%s%s.csv", cov$source,
                   if (nzchar(Sys.getenv("NAT_COMMON"))) "_common" else "")
write.csv(out, OUT_CSV, row.names = FALSE)

cat("\n=== Predicting the raw vs the method-standardised national prevalence ===\n")
cat("Offsets estimated on TRAINING countries only, applied to both sides.\n")
cat("folds_identified = held-out folds where the training data could identify\n")
cat("the method contrast at all (needs a country surveyed under >1 method).\n\n")
print(as.data.frame(out), row.names = FALSE)
message("Saved -> ", OUT_CSV)
