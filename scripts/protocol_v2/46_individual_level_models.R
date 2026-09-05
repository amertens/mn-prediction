# =============================================================================
# scripts/protocol_v2/46_individual_level_models.R   [IL-01]
#
# INDIVIDUAL-LEVEL MODELS: SURVEY-ONLY vs PROXIES vs BOTH, UNDER HONEST FOLDS
#
# The January 2026 deck compared person-level SuperLearners fitted on (a) the
# respondent's own survey variables, (b) district-level proxies, (c) both,
# reporting Brier skill over the null, ROC AUC and precision-recall gain
# (women's B12: AUC 0.689, PR gain 3.09x, Brier skill 0.095). Those numbers
# came from the old fitting path (no survey weights, folds not blocked by
# district). This is the same comparison under the corrected conventions:
#   learner    classic SuperLearner through fit_area_superlearner()
#              (R/area_superlearner.R): survey weights as obsWeights, inner
#              folds blocked by district, library mean / elastic net / ranger,
#              binomial family, discrete selection
#   folds      OUTER 5-fold by district, so a district's respondents are never
#              split across training and test (estimand A at person level)
#   survey     numeric gw_ columns with >= 70% coverage, minus every column
#              that names a biomarker, an adjustment, an id, a weight, a date
#              or a location (the config's gw_exclude_patterns plus the list
#              below); printed per country so the set can be audited
#   proxies    the district's domain principal components (built on the
#              country's surveyed districts, outcome-independent), joined to
#              each respondent by district
#   metrics    on pooled out-of-fold predictions: ROC AUC (Mann-Whitney),
#              Brier skill vs the national prevalence, PR-AUC and PR gain
#              (PR-AUC / prevalence), normalised PR-AUC
#
#   Rscript scripts/protocol_v2/46_individual_level_models.R      IL_CELLS=all|main
# -> results/tables/protocol_v2/individual_level_models.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; STORE <- "_targets_full"; set.seed(20260904L)
CELLS <- Sys.getenv("IL_CELLS", "main")
# IL_COUNTRY (comma list) shards the run by country; each shard writes its own
# file (individual_level_models_<tag>.csv) and scratchpad/merge_il.R joins them.
SHARD <- Sys.getenv("IL_COUNTRY", ""); MAXCOL <- as.integer(Sys.getenv("IL_MAXCOL", "80"))
cfg <- get_country_configs()
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)
LEAK <- "RBP|rbp|VAD|vad|LogFer|logfer|Ferr|ferr|IDA|ida$|Brinda|BRINDA|Thurn|Folate|folate|Fol|B12|b12|Zinc|zinc|zn_|Hb|hb$|Hgb|hgb|Haem|Hem|anaem|anem|Anem|CRP|crp|AGP|agp|infl|Infl|UIC|uic|Iod|iod|Salt|salt|MUAC|muac|weight|Weight|Wt$|wt$|_id$|ID$|Id$|cnum|clust|Clust|Date|date|month|Month|year|Year|psu|PSU|strat|Strat|Admin|admin|lat|Lat|lon|Lon|GPS|gps|Team|team|line|Line|hhid|caseid|Retinol|retinol|Vit|vit|Anemia|Malaria|malaria|RDT|rdt|Plasmod|Sickle|G6PD|Transferrin|sTfR|stfr|ZPP|zpp|Ret$|Def$|def$|Adj|adj"
auc <- function(y, p) { ok <- is.finite(p); y <- y[ok]; p <- p[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0); if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
prauc <- function(y, p) { ok <- is.finite(p); y <- y[ok]; p <- p[ok]; o <- order(-p); y <- y[o]; tp <- cumsum(y); prec <- tp / seq_along(y); rec <- tp / sum(y); if (sum(y) == 0) return(NA_real_)
  rec0 <- c(0, rec); prec0 <- c(1, prec); sum(diff(rec0) * (prec0[-1] + prec0[-length(prec0)]) / 2) }
metrics <- function(y, p, prev) { b <- mean((y - p)^2, na.rm = TRUE); b0 <- mean((y - prev)^2); pa <- prauc(y, p)
  data.frame(auc = auc(y, p), brier = b, brier_skill = 1 - b / b0, pr_auc = pa, pr_gain = pa / prev, pr_norm = (pa - prev) / (1 - prev)) }
group_folds <- function(groups, k, rep_id) { set.seed(20260951L + rep_id); g <- unique(groups); f <- sample(rep(seq_len(min(k, length(g))), length.out = length(g))); f[match(groups, g)] }

rows <- list()
for (cn in (if (nzchar(SHARD)) trimws(strsplit(SHARD, ",")[[1]]) else names(cfg))) { cc <- cfg[[cn]]; lc <- tolower(cn)
  outs <- names(cc$outcomes); if (CELLS == "main") outs <- intersect(outs, c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12"))
  # district proxies: domain PCs on the country's surveyed districts
  Sx <- S[S$country == cn, ]; Xr <- prep_predictors_v2(as.matrix(Sx[, PREDS, drop = FALSE])); Dm <- domain_representation_v2(Xr, domain_of)
  colnames(Dm) <- paste0("px_", colnames(Dm)); key_s <- paste(Sx$Admin1, Sx$Admin2)
  for (on in outs) { oc <- cc$outcomes[[on]]
    if (cn == "Malawi") {
      # Malawi's store data carries only the outcome and district aggregates; the
      # respondent questionnaire items (m01-m126, hunger score, fortification
      # exposure) live in the clean file. Outcome definitions there are the
      # survey's own, not the uniform resolver's -- noted in the log.
      MW <- readRDS("data/IPD/Malawi/clean_malawi_mn_data.RDS")
      defcol <- c(child_vitA = "vitA_def", women_vitA = "vitA_def", child_iron = "iron_def", women_iron = "iron_def", child_zinc = "zinc_def", women_zinc = "zinc_def")[on]
      if (is.na(defcol) || !defcol %in% names(MW)) { cat("  Malawi", on, "skip: no outcome column in the clean file\n"); next }
      grp <- if (grepl("^child", on)) is.finite(MW$psc_agecat) else is.finite(MW$women_agecat)
      d <- MW[grp, ]; y <- .v2_num(d[[defcol]]); w <- .v2_num(d$svy_weight); w[!is.finite(w) | w <= 0] <- NA
      dist <- paste(as.character(d$Admin1), as.character(d$Admin2)); ok <- is.finite(y) & is.finite(w) & !is.na(dist)
      gw <- names(d)[grepl("^m[0-9]+[a-g]?$|^mvisits$|^mlang|^mtype$|^hhs_|^oil_vita$|^sugar_vita$|^salt$|^fast$", names(d))]
    } else {
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL); if (is.null(od)) next
    d <- od$data; yb <- tryCatch(resolve_uniform_outcome(d, cc, oc), error = function(e) NULL)
    y <- if (!is.null(yb)) .v2_num(yb) else .v2_num(d[[oc$binary]]); w <- .v2_num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- NA
    dist <- paste(as.character(d$Admin1), as.character(d$Admin2)); ok <- is.finite(y) & is.finite(w) & !is.na(dist)
    outcome_vars <- unique(unlist(lapply(cc$outcomes, function(o) c(o$binary, o$continuous))))
    LEAK2 <- paste0(LEAK, "|^sf_|ferritin|_nmol|^fol|^zn|^map2_|^rbp|^vit|_def$")
    gw <- names(d)[grepl("^gw_", names(d)) & !grepl(LEAK2, names(d)) & !names(d) %in% outcome_vars]
    }
    gw <- gw[vapply(gw, function(k) { v <- .v2_num(d[[k]]); mean(is.finite(v)) >= 0.7 && stats::sd(v, na.rm = TRUE) > 0 && length(unique(v[is.finite(v)])) >= 2 }, TRUE)]
    # cap the survey set at IL_MAXCOL columns by coverage (outcome-independent), for runtime
    cov <- vapply(gw, function(k) mean(is.finite(.v2_num(d[[k]]))), 0); gw <- gw[order(-cov)][seq_len(min(MAXCOL, length(gw)))]
    if (length(gw) < 2) { cat(sprintf("  %-12s %-12s skip: %d usable survey columns\n", cn, on, length(gw))); next }
    Xs <- as.matrix(as.data.frame(lapply(d[gw], .v2_num))); Xs <- apply(Xs, 2, function(v) { v[!is.finite(v)] <- stats::median(v, na.rm = TRUE); v }); if (is.null(dim(Xs))) Xs <- matrix(Xs, ncol = 1, dimnames = list(NULL, gw))
    Xp <- Dm[match(dist, key_s), , drop = FALSE]; ok <- ok & is.finite(Xp[, 1])
    y <- y[ok]; w <- w[ok]; dist <- dist[ok]; Xs <- Xs[ok, , drop = FALSE]; Xp <- Xp[ok, , drop = FALSE]
    if (length(y) < 100 || sum(y) < 10 || dplyr::n_distinct(dist) < 6) { cat(sprintf("  %-12s %-12s skip (n %d, cases %d, districts %d)\n", cn, on, length(y), sum(y), dplyr::n_distinct(dist))); next }
    cat(sprintf("== %s %s: n %d, prevalence %.2f, districts %d | survey columns %d (e.g. %s) | proxy PCs %d\n", cn, on, length(y), mean(y), dplyr::n_distinct(dist), ncol(Xs), paste(head(colnames(Xs), 6), collapse = ", "), ncol(Xp)))
    SETS <- list(survey_only = Xs, proxies_only = Xp, both = cbind(Xs, Xp))
    folds <- group_folds(dist, 5, 1L)
    for (sname in names(SETS)) { X <- SETS[[sname]]; p <- rep(NA_real_, length(y)); picks <- c()
      for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f)
        fit <- tryCatch(fit_area_superlearner(y[tr], X[tr, , drop = FALSE], newX = X[te, , drop = FALSE], weights = w[tr], block = dist[tr], library = c("mean", "enet", "ranger"), V = 3L, discrete = TRUE, family = stats::binomial(), meta = "mse"), error = function(e) NULL)
        if (!is.null(fit) && !is.null(fit$pred_new) && length(fit$pred_new) == length(te)) { p[te] <- pmin(pmax(fit$pred_new, 1e-4), 1 - 1e-4); picks <- c(picks, fit$pick) } }
      m <- metrics(y, p, mean(y))
      rows[[length(rows) + 1L]] <- cbind(data.frame(country = cn, outcome = on, set = sname, n = length(y), n_districts = dplyr::n_distinct(dist), prevalence = mean(y), n_pred = ncol(X), picks = paste(names(table(picks)), table(picks), collapse = ";"), stringsAsFactors = FALSE), m)
      cat(sprintf("   %-13s AUC %.3f | Brier skill %+.3f | PR gain %.2fx | picks %s\n", sname, m$auc, m$brier_skill, m$pr_gain, paste(names(table(picks)), table(picks), collapse = ";"))) }
  } }
R <- bind_rows(rows); if (!nrow(R)) stop("no cells produced"); write.csv(R, file.path(OUTDIR, paste0("individual_level_models", if (nzchar(SHARD)) paste0("_", gsub("[^A-Za-z]", "", SHARD)) else "", ".csv")), row.names = FALSE)
cat("\n===== IL-01: individual-level models, 5-fold by district, survey-weighted SuperLearner =====\n")
print(as.data.frame(R |> group_by(set) |> summarise(cells = dplyr::n(), auc = round(mean(auc, na.rm = TRUE), 3), brier_skill = round(mean(brier_skill, na.rm = TRUE), 3), pr_gain = round(mean(pr_gain, na.rm = TRUE), 2), pr_norm = round(mean(pr_norm, na.rm = TRUE), 3), .groups = "drop")), row.names = FALSE)
W <- R |> select(country, outcome, set, brier_skill) |> tidyr::pivot_wider(names_from = set, values_from = brier_skill)
if (all(c("survey_only", "proxies_only", "both") %in% names(W))) cat(sprintf("Brier skill: survey-only beats proxies in %d of %d cells; both beats survey-only in %d of %d\n", sum(W$survey_only > W$proxies_only, na.rm = TRUE), sum(is.finite(W$survey_only - W$proxies_only)), sum(W$both > W$survey_only, na.rm = TRUE), sum(is.finite(W$both - W$survey_only))))
cat("\nper cell:\n"); print(as.data.frame(R |> select(country, outcome, set, n, prevalence, auc, brier_skill, pr_gain) |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nDONE\n")
