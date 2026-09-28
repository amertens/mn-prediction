# =============================================================================
# scripts/protocol_v2/69_index_person_level_lib.R   [IL-02c]
# Sourced by 69_index_person_level.R and each of its workers.
#
# PERSON-LEVEL AUC AND BRIER SCORE OF THE PROTOCOL'S PCA DOMAIN INDEX
#
# The index is fitted exactly as in 02_run_benchmarks_v2.R (build_cell, domain
# PCs to 80 percent of each domain's variance, make_folds_v2("kfold_district"),
# arm_domain_index_v2 / arm_domain_index_cal_v2) on the district targets in
# targets_v2.csv. The only addition: each held-out district's out-of-fold
# prediction is given to every respondent in that district, and the
# respondents' own deficiency flags are scored.
#
# Methods scored at person level
#   index_cal      the protocol's calibrated-level index (IS-01), expit -> predicted district prevalence
#   index_rank     the ranking index (rho = 1), expit; same ranking inside a fit, over-dispersed levels
#   index_glm      person-level logistic recalibration of the index score, fitted on training respondents
#   index_lvl_glm  the index trained on the LEVEL target (mean -log biomarker), logistic recalibration
#   null           the training respondents' prevalence
# Ceilings (not predictions: each reads other respondents' held-out outcomes)
#   oracle_district_insample  the district's survey-weighted prevalence, self included
#   oracle_district_loo_eb    the district's other respondents, shrunk to the national rate
#   oracle_cluster_loo_eb     the cluster's other respondents, shrunk
# Schemes: kfold5 = the protocol's in-fill folds, REPS draws; lodo = leave one district out.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
invisible(capture.output(suppressMessages(targets::tar_source("R"))))
STORE <- "_targets_full"; OUTDIR <- "results/tables/protocol_v2"
REPS <- as.integer(Sys.getenv("REPS", "20"))
MAIN <- c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12")

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))
# verbatim from 02_run_benchmarks_v2.R, plus the district key in the return value
build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  ncol_eff <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[ncol_eff]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  D  <- domain_representation_v2(Xr, domain_of)
  y_nat <- m[[ycol]]
  y_mod <- if (target == "prev") .v2_logit(y_nat) else y_nat
  list(country = cn, outcome = on, target = target, y_nat = y_nat, y_mod = y_mod, X = Xr, D = D,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = y_nat, target = target),
       w = m[[ncol_eff]], Admin1 = m$Admin1, n = nrow(m), key = paste(m$Admin1, m$Admin2, sep = "||"))
}

cfgs <- get_country_configs()
load_persons <- function(cn, on) {
  lc <- tolower(cn); cc <- cfgs[[cn]]; oc <- cc$outcomes[[on]]
  od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL)
  if (is.null(od)) return(NULL)
  d <- od$data
  derived <- tryCatch(suppressMessages(resolve_uniform_outcome(d, cc, oc, label = "[il-index]")), error = function(e) NULL)
  y <- if (!is.null(derived)) .v2_num(derived) else .v2_num(d[[oc$binary]])
  w <- .v2_num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- NA
  psu <- if (!is.null(cc$psu_col) && cc$psu_col %in% names(d)) as.character(d[[cc$psu_col]]) else NA_character_
  if (length(y) != nrow(d)) { message("  ", cn, " ", on, ": outcome column unavailable (length ", length(y), ")"); return(NULL) }
  if (length(w) != nrow(d)) { message("  ", cn, " ", on, ": weight column unavailable, unit weights used"); w <- rep(1, nrow(d)) }
  if (length(psu) != nrow(d)) psu <- rep(NA_character_, nrow(d))
  P <- data.frame(key = paste(as.character(d$Admin1), as.character(d$Admin2), sep = "||"), y = y, w = w, psu = psu, stringsAsFactors = FALSE)
  P[is.finite(P$y) & is.finite(P$w), ]
}

auc <- function(y, p) { ok <- is.finite(p) & is.finite(y); y <- y[ok]; p <- p[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auc_within <- function(y, p, fold) { a_ <- 0; den <- 0
  for (f in unique(fold)) { i <- which(fold == f); n1 <- sum(y[i] == 1); n0 <- sum(y[i] == 0); if (n1 == 0 || n0 == 0) next
    a <- auc(y[i], p[i]); if (!is.finite(a)) next; a_ <- a_ + a * n1 * n0; den <- den + n1 * n0 }
  if (den == 0) NA_real_ else a_ / den }
score <- function(y, p, p_null = NULL, fold = NULL) {
  ok <- is.finite(p); y <- y[ok]; p <- pmin(pmax(p[ok], 1e-4), 1 - 1e-4); pbar <- mean(y); b <- mean((y - p)^2)
  c(auc_pooled = auc(y, p), auc_within = if (is.null(fold)) NA_real_ else auc_within(y, p, fold[ok]), brier = b,
    bss_national = 1 - b / (pbar * (1 - pbar)),
    bss_vs_null = if (is.null(p_null)) NA_real_ else 1 - b / mean((y - p_null[ok])^2), coverage = mean(ok))
}
eb_loo <- function(y, g) {   # leave-self-out group rate shrunk to the overall rate (beta-binomial moments)
  s <- ave(y, g, FUN = sum); n <- ave(y, g, FUN = length); pbar <- mean(y)
  gm <- tapply(y, g, mean); gn <- tapply(y, g, length)
  tau2 <- max(stats::var(gm) - mean(pbar * (1 - pbar) / gn), 1e-6)
  a <- max(pbar * (1 - pbar) / tau2 - 1, 1)
  (s - y + a * pbar) / (n - 1 + a)
}

index_scores <- function(cell, tr) {   # standardised index score for every district, weights from training districts
  z <- .index_weights_v2(cell$D[tr, , drop = FALSE], cell$y_mod[tr]); s <- as.numeric(cell$D %*% z)
  sdv <- stats::sd(s[tr]); if (!is.finite(sdv) || sdv == 0) return(rep(0, length(s))); (s - mean(s[tr])) / sdv
}
glm_recal <- function(yp, sp, trp, tep) {
  fit <- tryCatch(suppressWarnings(stats::glm(y ~ s, family = stats::binomial(), data = data.frame(y = yp[trp], s = sp[trp]))), error = function(e) NULL)
  if (is.null(fit)) return(rep(mean(yp[trp]), length(tep)))
  as.numeric(stats::predict(fit, newdata = data.frame(s = sp[tep]), type = "response"))
}


run_cell <- function(cn, on) {
  rows <- list()
  cp <- build_cell(cn, on, "prev"); if (is.null(cp)) return(list(rows = NULL, chk = NULL, note = "skip: no prevalence cell"))
  cl <- build_cell(cn, on, "level")
  P <- load_persons(cn, on); if (is.null(P)) return(list(rows = NULL, chk = NULL, note = "skip: no usable respondent outcome"))
  P <- P[P$key %in% cp$key, ]; if (nrow(P) < 50 || sum(P$y) < 5) return(list(rows = NULL, chk = NULL, note = "skip: too few respondents or cases"))
  pd <- match(P$key, cp$key); pl <- if (!is.null(cl)) match(P$key, cl$key) else rep(NA_integer_, nrow(P))
  # sanity: persons' weighted district prevalence reproduces the target
  wp <- tapply(seq_len(nrow(P)), pd, function(ii) stats::weighted.mean(P$y[ii], P$w[ii]))
  chk <- data.frame(country = cn, outcome = on, n_persons = nrow(P), n_districts = cp$n,
    districts_with_persons = length(wp), r_target = round(stats::cor(wp, cp$y_nat[as.integer(names(wp))]), 4))
  y <- P$y
  # ---- ceilings -----------------------------------------------------------
  orc <- list(oracle_district_insample = cp$y_nat[pd],
              oracle_district_loo_eb = eb_loo(y, P$key),
              oracle_cluster_loo_eb = if (all(is.na(P$psu))) rep(NA_real_, nrow(P)) else eb_loo(y, paste(P$key, P$psu)))
  for (m in names(orc)) rows[[length(rows) + 1]] <- data.frame(country = cn, outcome = on, scheme = "none", rep = 0, method = m,
    n = nrow(P), prev = mean(y), n_districts = cp$n, t(score(y, orc[[m]])))
  # ---- honest out-of-fold index -------------------------------------------
  schemes <- list(kfold5 = seq_len(REPS), lodo = 1)
  for (sch in names(schemes)) for (r in schemes[[sch]]) {
    folds <- if (sch == "kfold5") make_folds_v2("kfold_district", cp$n, rep_id = r) else seq_len(cp$n)
    p_cal <- p_rank <- rep(NA_real_, cp$n)
    pp <- matrix(NA_real_, nrow(P), 3, dimnames = list(NULL, c("index_glm", "index_lvl_glm", "null"))); fold_p <- folds[pd]
    for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
      p_rank[te] <- arm_domain_index_v2(tr, te, cp$y_mod, cp$X, cp$D, cp$aux)
      p_cal[te]  <- arm_domain_index_cal_v2(tr, te, cp$y_mod, cp$X, cp$D, cp$aux)
      trp <- which(pd %in% tr); tep <- which(pd %in% te); if (!length(tep)) next
      s <- index_scores(cp, tr); pp[tep, "index_glm"] <- glm_recal(y, s[pd], trp, tep)
      pp[tep, "null"] <- mean(y[trp])
      if (!is.null(cl)) {   # level-trained index: same held-out districts, trained on the level target
        trl <- which(cl$key %in% cp$key[tr]); if (length(trl) >= 12) {
          sl <- index_scores(cl, trl); ok_l <- is.finite(pl)
          trp_l <- trp[ok_l[trp]]; tep_l <- tep[ok_l[tep]]
          if (length(tep_l)) pp[tep_l, "index_lvl_glm"] <- glm_recal(y, sl[pl], trp_l, tep_l) } }
    }
    preds <- list(index_cal = .v2_expit(p_cal)[pd], index_rank = .v2_expit(p_rank)[pd], index_glm = pp[, "index_glm"],
                  index_lvl_glm = pp[, "index_lvl_glm"], null = pp[, "null"])
    for (m in names(preds)) rows[[length(rows) + 1]] <- data.frame(country = cn, outcome = on, scheme = sch, rep = r, method = m,
      n = nrow(P), prev = mean(y), n_districts = cp$n, t(score(y, preds[[m]], pp[, "null"], if (sch == "lodo") NULL else fold_p)),
      district_spearman = suppressWarnings(stats::cor(if (m %in% c("index_cal", "index_rank")) p_rank else tapply(preds[[m]], pd, mean)[as.character(seq_len(cp$n))], cp$y_nat, method = "spearman", use = "complete.obs")))
  }
  cat(sprintf("%-12s %-12s persons %5d districts %3d prev %.3f done\n", cn, on, nrow(P), cp$n, mean(y)))

  list(rows = bind_rows(rows), chk = chk, note = "ok")
}
