# =============================================================================
# scripts/protocol_v2/57_index_importance.R   [WS-02]
#
# WHICH PREDICTORS DOES THE INDEX ACTUALLY WEIGHT, AND DO THEY REPLICATE?
#
# Projects the zero-tuning index's axis weights back onto the original
# rank-normalised predictors (R/protocol_v2_importance.R) for
#   (a) each country's own fit on all of its surveyed districts, with a
#       stability check over the five training folds of one in-fill draw;
#   (b) the pooled four-country fit that would be deployed to a fifth country,
#       with the four leave-one-country-out fits as its replication check.
# For every predictor and outcome the tables carry: beta (weight per unit of
# rank-normal score), beta_std (times the training SD), share (exact share of
# the index's variance), the rank of |beta_std| within the fit, and how many
# of the four LOCO fits and of the four in-country fits agree with the pooled
# sign. Domain shares are exact and sum to one.
#
#   Rscript scripts/protocol_v2/57_index_importance.R
# -> results/tables/protocol_v2/index_importance_columns.csv   long table, every fit
#    results/tables/protocol_v2/index_importance_domains.csv   domain shares, every fit
#    results/tables/protocol_v2/index_importance_top.csv       pooled top 20 per outcome x target with replication
#    results/tables/protocol_v2/index_importance_patterns.csv  predictors in the pooled top 20 of >= 2 outcomes
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R")

OUTDIR <- "results/tables/protocol_v2"
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")

build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y_nat <- m[[ycol]]
  list(country = cn, outcome = on, target = target, n = nrow(m), y_nat = y_nat,
       y_mod = if (target == "prev") .v2_logit(y_nat) else y_nat, X = Xr, Admin1 = m$Admin1)
}

one_fit <- function(tr, y, X, D, tag) {
  im <- index_importance_v2(tr, y, X, D)
  cols <- im$columns |> mutate(rank = rank(-abs(beta_std), ties.method = "first"))
  list(columns = cbind(tag, cols, stringsAsFactors = FALSE), domains = cbind(tag, im$domains, stringsAsFactors = FALSE))
}

COLS <- list(); DOMS <- list()
# ── (a) in-country fits, with fold stability ─────────────────────────────────
INC_SIGN <- list()   # per country x outcome x target: sign of beta in the full fit
for (cn in COUNTRIES) for (on in unique(TG$outcome)) for (target in c("level", "prev")) {
  cell <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (is.null(cell)) next
  D <- domain_representation_v2(cell$X, domain_of)
  tag <- data.frame(scope = "country", fit = cn, country = cn, outcome = on, target = target, n_train = cell$n)
  f <- one_fit(seq_len(cell$n), cell$y_mod, cell$X, D, tag)
  folds <- make_folds_v2("kfold_district", cell$n, k = 5, rep_id = 1L)
  B <- sapply(sort(unique(folds)), function(k) { tr <- which(folds != k); if (length(tr) < 12) return(rep(NA_real_, ncol(cell$X)))
    index_backproject_v2(.ws_z_pooled(tr, cell$y_mod, D), attr(D, "basis"), colnames(cell$X)) })
  agree <- rowMeans(sign(B) == sign(f$columns$beta), na.rm = TRUE)
  f$columns$fold_sign_agree <- as.numeric(agree)
  COLS[[paste(cn, on, target)]] <- f$columns; DOMS[[paste(cn, on, target)]] <- f$domains
  INC_SIGN[[paste(on, target, cn)]] <- stats::setNames(sign(f$columns$beta), f$columns$column)
  cat("in-country", cn, on, target, "\n"); flush.console()
}

# ── (b) pooled and leave-one-country-out fits ────────────────────────────────
TOP <- list()
for (on in unique(TG$outcome)) for (target in c("level", "prev")) {
  cl <- list()
  for (cn in COUNTRIES) { cc <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(cc)) cl[[cn]] <- cc }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  Dall <- domain_representation_v2(Xm, domain_of)
  tag <- data.frame(scope = "pooled", fit = "all", country = "all", outcome = on, target = target, n_train = nrow(Xm))
  f <- one_fit(seq_len(nrow(Xm)), Y, Xm, Dall, tag)
  # LOCO fits: orientation and weights from the training countries only
  L <- sapply(names(cl), function(cn) { tr <- which(ctry != cn); Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    b <- index_backproject_v2(.ws_z_pooled(tr, Y, Dm), attr(Dm, "basis"), colnames(Xm))
    fl <- one_fit(tr, Y, Xm, Dm, data.frame(scope = "loco", fit = paste0("holdout_", cn), country = cn, outcome = on, target = target, n_train = length(tr)))
    COLS[[paste("loco", cn, on, target)]] <<- fl$columns; DOMS[[paste("loco", cn, on, target)]] <<- fl$domains
    b })
  f$columns$loco_sign_agree <- rowSums(sign(L) == sign(f$columns$beta))
  f$columns$loco_fits <- ncol(L)
  inc <- sapply(names(cl), function(cn) { s <- INC_SIGN[[paste(on, target, cn)]]; if (is.null(s)) rep(NA, nrow(f$columns)) else s[f$columns$column] })
  f$columns$incountry_sign_agree <- rowSums(inc == sign(f$columns$beta), na.rm = TRUE)
  f$columns$incountry_fits <- rowSums(!is.na(inc))
  f$columns$domain <- unname(domain_of[f$columns$column])
  COLS[[paste("pooled", on, target)]] <- f$columns; DOMS[[paste("pooled", on, target)]] <- f$domains
  TOP[[paste(on, target)]] <- f$columns |> arrange(rank) |> filter(rank <= 20)
  cat("pooled", on, target, "\n"); flush.console()
}

COLS <- bind_rows(COLS); DOMS <- bind_rows(DOMS)
COLS$domain <- unname(domain_of[COLS$column])
write.csv(COLS, file.path(OUTDIR, "index_importance_columns.csv"), row.names = FALSE)
write.csv(DOMS, file.path(OUTDIR, "index_importance_domains.csv"), row.names = FALSE)
TOP <- bind_rows(TOP) |> select(outcome, target, rank, column, domain, beta, beta_std, share, loco_sign_agree, loco_fits, incountry_sign_agree, incountry_fits)
write.csv(TOP, file.path(OUTDIR, "index_importance_top.csv"), row.names = FALSE)
PAT <- TOP |> group_by(target, column, domain) |>
  summarise(outcomes_in_top20 = n(), outcomes = paste(outcome, collapse = ", "), signs = paste(ifelse(beta > 0, "+", "-"), collapse = ""),
            mean_rank = round(mean(rank), 1), min_loco_agree = min(loco_sign_agree), .groups = "drop") |>
  filter(outcomes_in_top20 >= 2) |> arrange(target, desc(outcomes_in_top20), mean_rank)
write.csv(PAT, file.path(OUTDIR, "index_importance_patterns.csv"), row.names = FALSE)

# ── console report ───────────────────────────────────────────────────────────
chk <- COLS |> group_by(scope, fit, outcome, target) |> summarise(sum_share = sum(share), .groups = "drop")
cat(sprintf("\nshare sums to one: min %.3f max %.3f over %d fits\n", min(chk$sum_share), max(chk$sum_share), nrow(chk)))
for (target in c("level", "prev")) {
  cat("\n==================== target:", target, "(pooled four-country fit; + = more deficiency) ====================\n")
  for (on in unique(TOP$outcome)) {
    d <- TOP |> filter(outcome == on, target == !!target, rank <= 10) |>
      transmute(rank, predictor = column, domain = substr(domain, 1, 22), sign = ifelse(beta > 0, "+", "-"), beta_std = round(beta_std, 3),
                share = round(share, 3), loco = paste0(loco_sign_agree, "/", loco_fits), in_country = paste0(incountry_sign_agree, "/", incountry_fits))
    ds <- DOMS |> filter(scope == "pooled", outcome == on, target == !!target) |> arrange(desc(share)) |> head(5)
    cat("\n--", on, "| domain shares:", paste0(substr(ds$domain, 1, 18), " ", round(ds$share, 2), collapse = "; "), "\n")
    print(as.data.frame(d), row.names = FALSE)
  }
  cat("\n-- predictors in the pooled top 20 of two or more outcomes --\n")
  print(as.data.frame(PAT |> filter(target == !!target) |> head(25) |> select(-target)), row.names = FALSE)
}
