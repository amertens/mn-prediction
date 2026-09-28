# =============================================================================
# scripts/protocol_v2/68_individual_decomposition_lib.R   [IL-02b]
# Sourced by 68_individual_decomposition.R and by each of its workers. Not in
# R/ because it reads the January-era data at top level.
#
# One factor at a time, on the January-era Ghana data, from January's own
# cross-validation to the IL-01 protocol. ONE fixed learner stack throughout
# (mean, lasso, ridge, ranger; NNLS weights on the out-of-fold predictions, as
# the January sl3 Lrnr_sl did), so a change in skill is a protocol change, not
# a learner change. Reads the mn-proxies copy of Ghana_merged_dataset.rds
# (13 January 2026), which reproduces the n and prevalence of five of the six
# saved January fits (women's vitamin A has no column in it).
# =============================================================================
suppressPackageStartupMessages({library(glmnet); library(ranger); library(nnls)})
DATA <- file.path(Sys.getenv("MN_PROXIES", "C:/Users/andre/OneDrive/Documents/mn-proxies"), "data/IPD/Ghana/Ghana_merged_dataset.rds")
df <- readRDS(DATA)
num <- function(v) suppressWarnings(as.numeric(unclass(haven::zap_labels(v))))

# ---- January predictor sets, rules copied from src/Ghana/GW Ghana SL_bin_full.R ----
gw_vars <- names(df)[grepl("gw_", names(df))]
for (p in c("wID","cID","VAD","RBP","rbp","Ferr","TFR","crp","agp","cRBP","cVAD","VAI","Folate","B12",
            "AGP","CRP","RDT","fer","bis","nflam","nemia","globin","cHb","wHb","gchb","hbc","Anemia","wm_wmst"))
  gw_vars <- gw_vars[!grepl(p, gw_vars)]
gw_vars <- gw_vars[!gw_vars %in% c("gw_cn", "gw_hhid", "gw_childid")]
grab <- function(p) names(df)[grepl(p, names(df))]
dhs <- names(df)[grepl("dhs2014_", names(df)) | grepl("dhs2016_", names(df)) | grepl("dhs2017_", names(df))]
proxy_vars <- intersect(unique(c(dhs, grab("mics_"), grab("ihme_"), grab("lsms_"), grab("MAP_"),
                                 "nearest_market_id", grab("wfp_"), grab("flunet_"), grab("gee_"))), names(df))
# dataid is a unique string per row; it never reached the learners in the saved fits, so it is left out
SET <- list(survey_jan = c("Admin1", "gw_month", gw_vars),
            proxy_jan  = c("Admin1", "gw_month", proxy_vars),
            proxy_clean = proxy_vars)
# identifiers, sampling design, region, field-team / phlebotomist codes, dates, and blood-derived
# columns (malaria RDT and referral, sickle / thalassaemia genotype, MRDR selection), plus the
# cluster-level open-defecation share
CLEAN_DROP <- paste0("cnum|^gw_cn$|bccn|b_cn$|hhln|bchnc|bclnr|bcgln|b_hhn|wlnr|^gw_in1$|^gw_wp1$|b_wp1|indivID|momID|",
  "pcn$|^gw_mcn$|Team|Region|region|Strata|strata|Number_of|probab|PSU|Check_|sWeight|Total.population|Households|",
  "cttp|datasource|INT_m|a_mon|VAS_m2|VAS_d$|VAS_y$|VAS_date|gw_month|gcmst|MalariaYN|wrmal|malref|",
  "sickle|Sickle|thal|Thal|mrdr|pctopendef|^Admin1$|^dataid$")
SET$survey_clean <- SET$survey_jan[!grepl(CLEAN_DROP, SET$survey_jan)]
# IL-01's own column rule (scripts/protocol_v2/46_individual_level_models.R), applied to these columns
IL_LEAK <- paste0("RBP|rbp|VAD|vad|LogFer|logfer|Ferr|ferr|IDA|ida$|Brinda|BRINDA|Thurn|Folate|folate|Fol|B12|b12|Zinc|zinc|zn_|Hb|hb$|Hgb|hgb|Haem|Hem|anaem|anem|Anem|CRP|crp|AGP|agp|infl|Infl|UIC|uic|Iod|iod|Salt|salt|MUAC|muac|weight|Weight|Wt$|wt$|_id$|ID$|Id$|cnum|clust|Clust|Date|date|month|Month|year|Year|psu|PSU|strat|Strat|Admin|admin|lat|Lat|lon|Lon|GPS|gps|Team|team|line|Line|hhid|caseid|Retinol|retinol|Vit|vit|Anemia|Malaria|malaria|RDT|rdt|Plasmod|Sickle|G6PD|Transferrin|sTfR|stfr|ZPP|zpp|Ret$|Def$|def$|Adj|adj",
                  "|^sf_|ferritin|_nmol|^fol|^zn|^map2_|^rbp|^vit|_def$")

OUTC <- list(
  child_vitA   = list(y = "gw_cVAD",         rows = function(d) is.finite(num(d$gw_cVAD))),
  women_b12    = list(y = "gw_wB12Def",      rows = function(d) is.finite(num(d$gw_wB12Def))),
  women_folate = list(y = "gw_wFolateDef",   rows = function(d) is.finite(num(d$gw_wFolateDef))),
  child_iron   = list(y = "gw_cIDAdjBrinda", rows = function(d) !is.na(d$gw_childid) & is.finite(num(d$gw_cIDAdjBrinda))),
  women_iron   = list(y = "gw_wIDAdjBrinda", rows = function(d) is.na(d$gw_childid) & is.finite(num(d$gw_wIDAdjBrinda))))

il01_cols <- function(d) {
  il <- names(d)[grepl("^gw_", names(d)) & !grepl(IL_LEAK, names(d))]
  il <- il[vapply(il, function(k) { v <- num(d[[k]]); mean(is.finite(v)) >= 0.7 && isTRUE(stats::sd(v, na.rm = TRUE) > 0) && length(unique(v[is.finite(v)])) >= 2 }, TRUE)]
  cov <- vapply(il, function(k) mean(is.finite(num(d[[k]]))), 0)
  il[order(-cov)][seq_len(min(80, length(il)))]
}

# ---- unsupervised preprocessing, January order (one-hot, drop constant / near-zero, median impute + indicators) ----
prep_X <- function(d, cols) {
  cols <- intersect(cols, names(d)); out <- list()
  for (k in cols) {
    v <- d[[k]]
    if (is.character(v) || is.factor(v)) {
      v <- as.character(v); lv <- sort(unique(v[!is.na(v)]))
      if (length(lv) < 2 || length(lv) > 30) next
      for (l in lv[-1]) out[[paste0(k, "__", make.names(l))]] <- ifelse(is.na(v), NA_real_, as.numeric(v == l))
    } else out[[k]] <- num(v)
  }
  X <- as.data.frame(out, check.names = FALSE)
  X <- X[, vapply(X, function(v) length(unique(v[is.finite(v)])) > 1, TRUE), drop = FALSE]
  nz <- caret::nearZeroVar(X); if (length(nz)) X <- X[, -nz, drop = FALSE]
  for (k in names(X)) { v <- X[[k]]; m <- !is.finite(v)
    if (any(m)) { X[[paste0("missing_", k)]] <- as.numeric(m); v[m] <- stats::median(v, na.rm = TRUE); X[[k]] <- v } }
  nz <- caret::nearZeroVar(X); if (length(nz)) X <- X[, -nz, drop = FALSE]
  X <- as.matrix(X); colnames(X) <- make.names(colnames(X), unique = TRUE); X
}
# washb_prescreen(family = "gaussian", pval = 0.2): per-column LRT of Y ~ W against Y ~ 1,
# which for a numeric column is LR = -n log(1 - r^2) on 1 df
prescreen <- function(y, X, pval = 0.2) {
  r <- suppressWarnings(stats::cor(X, y)); r[!is.finite(r)] <- 0
  p <- stats::pchisq(-length(y) * log(pmax(1 - r^2, 1e-12)), 1, lower.tail = FALSE)
  keep <- which(p < pval); if (length(keep) < 2) keep <- order(p)[1:min(2, ncol(X))]
  colnames(X)[keep]
}
corr_filter <- function(X, cutoff = 0.9) {
  if (ncol(X) < 3) return(colnames(X))
  C <- suppressWarnings(stats::cor(X)); C[!is.finite(C)] <- 0
  drop <- caret::findCorrelation(C, cutoff = cutoff); if (length(drop)) colnames(X)[-drop] else colnames(X)
}
block_folds <- function(blocks, K, seed) {
  set.seed(seed); u <- unique(blocks); f <- sample(rep(seq_len(min(K, length(u))), length.out = length(u))); f[match(blocks, u)]
}
auc <- function(y, p) { ok <- is.finite(p) & is.finite(y); y <- y[ok]; p <- p[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
# AUC over pairs that share a fold only: removes the fold-to-fold offset in the training mean
auc_within <- function(y, p, fold) { num_ <- 0; den <- 0
  for (f in unique(fold)) { i <- fold == f; n1 <- sum(y[i] == 1); n0 <- sum(y[i] == 0); if (n1 == 0 || n0 == 0) next
    a <- auc(y[i], p[i]); if (!is.finite(a)) next; num_ <- num_ + a * n1 * n0; den <- den + n1 * n0 }
  if (den == 0) NA_real_ else num_ / den }
metrics <- function(y, p, p_null, fold) {
  b <- mean((y - p)^2); pbar <- mean(y)
  c(auc_pooled = auc(y, p), auc_within = auc_within(y, p, fold), brier = b,
    bss_national = 1 - b / (pbar * (1 - pbar)), bss_vs_trainmean = 1 - b / mean((y - p_null)^2))
}

fit_learners <- function(Xtr, ytr, Xte, inner) {
  P <- matrix(NA_real_, nrow(Xte), 4, dimnames = list(NULL, c("mean", "lasso", "ridge", "ranger")))
  P[, "mean"] <- mean(ytr)
  for (a in c(lasso = 1, ridge = 0)) {
    nm <- if (a == 1) "lasso" else "ridge"
    fit <- tryCatch(glmnet::cv.glmnet(Xtr, ytr, family = "binomial", alpha = a, foldid = inner, nlambda = 50),
                    error = function(e) NULL)
    P[, nm] <- if (is.null(fit)) mean(ytr) else as.numeric(stats::predict(fit, Xte, s = "lambda.min", type = "response"))
  }
  rf <- tryCatch(ranger::ranger(x = as.data.frame(Xtr), y = factor(ytr, levels = c(0, 1)), probability = TRUE,
                                num.trees = 500, num.threads = 1, seed = 1L), error = function(e) NULL)
  P[, "ranger"] <- if (is.null(rf)) mean(ytr) else stats::predict(rf, as.data.frame(Xte))$predictions[, "1"]
  P
}

run_task <- function(variant, outcome, rep) {
  V <- VARIANTS[[variant]]; oc <- OUTC[[outcome]]
  d <- df[oc$rows(df), , drop = FALSE]; y <- num(d[[oc$y]])
  cols <- if (V$set == "survey_il01") il01_cols(d) else SET[[V$set]]
  X <- prep_X(d, cols)
  blocks <- if (V$scheme == "cluster") as.character(d$gw_cnum) else paste(d$Admin1, d$Admin2)
  seed <- 20260927L + 100L * rep + V$K + if (V$scheme == "cluster") 0L else 50L
  fold <- block_folds(blocks, V$K, seed)
  keep_global <- if (V$screen == "global") prescreen(y, X) else colnames(X)
  if (V$screen == "global") keep_global <- corr_filter(X[, keep_global, drop = FALSE])
  Z <- matrix(NA_real_, length(y), 4, dimnames = list(NULL, c("mean", "lasso", "ridge", "ranger"))); p_null <- rep(NA_real_, length(y)); nsel <- c()
  for (f in sort(unique(fold))) {
    tr <- which(fold != f); te <- which(fold == f)
    sel <- keep_global
    if (V$screen == "infold") { sel <- prescreen(y[tr], X[tr, , drop = FALSE]); sel <- corr_filter(X[tr, sel, drop = FALSE]) }
    if (V$screen == "none") sel <- corr_filter(X[tr, , drop = FALSE])
    nsel <- c(nsel, length(sel))
    Xtr <- X[tr, sel, drop = FALSE]; Xte <- X[te, sel, drop = FALSE]
    mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    Xtr <- sweep(sweep(Xtr, 2, mu), 2, sdv, "/"); Xte <- sweep(sweep(Xte, 2, mu), 2, sdv, "/")
    inner <- block_folds(blocks[tr], 5, seed + f)
    Z[te, ] <- fit_learners(Xtr, y[tr], Xte, inner)
    p_null[te] <- mean(y[tr])
  }
  w <- tryCatch(nnls::nnls(Z, y)$x, error = function(e) c(1, 0, 0, 0)); if (sum(w) <= 0) w <- c(1, 0, 0, 0); w <- w / sum(w)
  p_sl <- as.numeric(Z %*% w)
  res <- rbind(stack = metrics(y, p_sl, p_null, fold), t(sapply(colnames(Z), function(k) metrics(y, Z[, k], p_null, fold))))
  data.frame(variant = variant, set = V$set, scheme = paste0(V$scheme, V$K), screen = V$screen, outcome = outcome, rep = rep,
             learner = rownames(res), n = length(y), prev = mean(y), n_cols_in = ncol(X), n_cols_used_median = stats::median(nsel),
             w_lasso = w[2], w_ridge = w[3], w_ranger = w[4], res, row.names = NULL, check.names = FALSE)
}

VARIANTS <- list(
  S1  = list(set = "survey_jan",   scheme = "cluster",  K = 10, screen = "global"),   # January as run (CV numbers, not the figure)
  S2  = list(set = "survey_jan",   scheme = "cluster",  K = 10, screen = "infold"),   # + prescreen inside folds
  S3  = list(set = "survey_clean", scheme = "cluster",  K = 10, screen = "infold"),   # + ids / design / blood columns removed
  S4b = list(set = "survey_jan",   scheme = "district", K = 10, screen = "infold"),   # January columns, district-blocked
  S4  = list(set = "survey_clean", scheme = "district", K = 10, screen = "infold"),   # clean, district-blocked
  S5  = list(set = "survey_clean", scheme = "district", K = 5,  screen = "infold"),   # IL-01's 5 folds
  S6  = list(set = "survey_il01",  scheme = "district", K = 5,  screen = "infold"),   # IL-01's column rule
  S7  = list(set = "survey_clean", scheme = "district", K = 5,  screen = "none"),     # no prescreen at all
  P1  = list(set = "proxy_jan",    scheme = "cluster",  K = 10, screen = "global"),
  P2  = list(set = "proxy_jan",    scheme = "cluster",  K = 10, screen = "infold"),
  P3  = list(set = "proxy_clean",  scheme = "district", K = 10, screen = "infold"),
  P4  = list(set = "proxy_clean",  scheme = "district", K = 5,  screen = "infold"))
