# =============================================================================
# scripts/protocol_v2/71_person_level_honest.R   [IL-02e]
#
# The January person-level Brier figure, updated to the current data and models
# and scored out of fold, with the index added; plus its continuous twin (MSE on
# log concentrations). See the lib for the arms.
#
# Skill = 1 - loss / loss of the null, where the null predicts the training
# folds' mean for everyone (the survey prevalence, or mean log concentration) on
# the same folds. Predictions are averaged over REPS fold draws; intervals are
# 95 percent cluster (PSU) bootstrap intervals of that averaged skill. Also
# reported per outcome: the between-district share of person-level variance
# (tau^2 / total, method of moments) and, for flags, the AUC a predictor that
# knew every district's TRUE prevalence would reach.
#
#   NW=17 REPS=5 BOOT=1000 IL_CACHE=<local dir> Rscript scripts/protocol_v2/71_person_level_honest.R
#   IL_HONEST_COUNTRY=Ghana (default)
# -> results/tables/protocol_v2/il02_honest_person_level.csv          aggregate metrics only
# -> results/tables/protocol_v2/il02_honest_survey_columns.csv        the survey columns used
# -> $IL_CACHE/il02_honest_person_predictions.rds                     respondent-level; never commit
# Then scripts/protocol_v2/72_plot_person_level_honest.R draws the figures.
# =============================================================================
suppressPackageStartupMessages({library(parallel); library(dplyr)})
ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"; setwd(ROOT)
LIB <- normalizePath("scripts/protocol_v2/71_person_level_honest_lib.R", winslash = "/")
OUTDIR <- "results/tables/protocol_v2"; CACHE <- Sys.getenv("IL_CACHE", tempdir())
NW <- as.integer(Sys.getenv("NW", "17")); REPS <- as.integer(Sys.getenv("REPS", "5")); B <- as.integer(Sys.getenv("BOOT", "1000"))
source(LIB)
TYPES <- c("bin", "cont")
cat(sprintf("%s: outcomes %s\n", COUNTRY71, paste(OUTS71, collapse = ", ")))
cols <- bind_rows(lapply(OUTS71, function(o) { cl <- load_cell71(o)
  cat(sprintf("  %-12s survey columns %3d | proxy PCs %3d | index PCs %3d | districts %d | concentration: %s\n", o, length(cl$survey_cols),
              ncol(cl$Xp), ncol(cl$Dix), length(cl$keys), cl$cont_src))
  data.frame(country = COUNTRY71, outcome = o, column = cl$survey_cols, concentration_source = cl$cont_src) }))
write.csv(cols, file.path(OUTDIR, "il02_honest_survey_columns.csv"), row.names = FALSE)
# workers read these compact per-outcome matrices instead of each loading the store's
# full outcome datasets (ten workers doing that ran the machine out of memory)
CELLS <- file.path(CACHE, paste0("il02_honest_cells_", COUNTRY71, ".rds"))
saveRDS(mget(OUTS71, envir = .cache71), CELLS)

cache_file <- file.path(CACHE, paste0("il02_honest_person_predictions_", COUNTRY71, ".rds"))
if (file.exists(cache_file) && Sys.getenv("IL_REUSE", "0") == "1") {
  K <- readRDS(cache_file); P <- K$P; IDX <- K$IDX; CEIL <- K$CEIL
} else {
  tasks <- expand.grid(type = TYPES, outcome = OUTS71, set = c("survey", "proxies", "both"), rep = seq_len(REPS), stringsAsFactors = FALSE)
  tasks <- tasks[order(tasks$set == "proxies"), ]
  cl <- makePSOCKcluster(NW)
  clusterExport(cl, c("LIB", "ROOT", "tasks", "CELLS"))
  invisible(clusterEvalQ(cl, { setwd(ROOT); source(LIB); K <- readRDS(CELLS); for (nm in names(K)) assign(nm, K[[nm]], envir = .cache71); rm(K); gc() }))
  t0 <- Sys.time()
  res <- parLapplyLB(cl, seq_len(nrow(tasks)), function(i) { g <- tasks[i, ]
    tryCatch(run_sl71(g$type, g$outcome, g$set, g$rep), error = function(e)
      data.frame(type = g$type, outcome = g$outcome, set = g$set, rep = g$rep, error = conditionMessage(e))) })
  stopCluster(cl)
  cat(sprintf("SL arms: %d fits in %.1f min\n", nrow(tasks), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  P <- bind_rows(res)
  if ("error" %in% names(P) && any(!is.na(P$error))) { print(unique(P[!is.na(P$error), c("type", "outcome", "set", "error")])); P <- P[is.na(P$error), ] }
  IDX  <- bind_rows(lapply(TYPES, function(tp) bind_rows(lapply(OUTS71, function(o) bind_rows(lapply(seq_len(REPS), function(r) run_index71(tp, o, r)))))))
  CEIL <- bind_rows(lapply(TYPES, function(tp) bind_rows(lapply(OUTS71, function(o) run_ceiling71(tp, o)))))
  saveRDS(list(P = P, IDX = IDX, CEIL = CEIL), cache_file)
}
cat("\nlearner picks (survey / proxies / both, all folds and draws):\n")
print(as.data.frame(P |> distinct(type, outcome, set, rep, picks) |> group_by(type, outcome, set) |> summarise(picks = paste(picks, collapse = " "), .groups = "drop")), row.names = FALSE)

# ---- metrics -----------------------------------------------------------------
auc1 <- function(y, p) { n1 <- sum(y == 1); n0 <- sum(y == 0); if (!n1 || !n0) return(NA_real_); r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auc_w <- function(y, p, fold) { a_ <- 0; den <- 0
  for (f in unique(fold)) { i <- which(fold == f); n1 <- sum(y[i] == 1); n0 <- sum(y[i] == 0); if (!n1 || !n0) next
    a_ <- a_ + auc1(y[i], p[i]) * n1 * n0; den <- den + n1 * n0 }; if (den == 0) NA_real_ else a_ / den }
avg <- bind_rows(P, IDX) |> group_by(type, outcome, set, row) |> summarise(y = first(y), pred = mean(pred), .groups = "drop")
nulls <- P |> filter(set == "survey") |> group_by(type, outcome, row) |> summarise(null = mean(null), .groups = "drop")
within_auc <- bind_rows(P, IDX) |> filter(type == "bin") |> group_by(outcome, set, rep) |>
  summarise(a = auc_w(y, pred, fold), .groups = "drop") |> group_by(outcome, set) |>
  summarise(auc_within = mean(a, na.rm = TRUE), .groups = "drop") |> mutate(type = "bin")
A <- bind_rows(avg, CEIL |> select(type, outcome, set, row, y, pred)) |> left_join(nulls, by = c("type", "outcome", "row"))
set.seed(20260927L); out <- list()
for (tp in TYPES) for (o in OUTS71) {
  ce <- CEIL[CEIL$type == tp & CEIL$outcome == o, ]; if (!nrow(ce)) next
  ce <- ce[order(ce$row), ]; psu <- psu71(tp, o); stopifnot(length(psu) == nrow(ce))
  idx_by_cl <- split(seq_along(psu), psu); cl_ids <- names(idx_by_cl)
  draws <- replicate(B, unlist(idx_by_cl[sample(cl_ids, length(cl_ids), replace = TRUE)], use.names = FALSE), simplify = FALSE)
  mu <- mean(ce$y); true_auc <- NA_real_
  g <- dist71(tp, o); share <- mom_share71(ce$y, g, tp)
  if (tp == "bin") {   # AUC of a predictor that knew every district's TRUE prevalence (Beta with the moment variance)
    t2 <- max(share * mu * (1 - mu), 1e-6); k <- max(mu * (1 - mu) / t2 - 1, 1e-3); gr <- seq(0.00025, 0.99975, by = 0.0005)
    f <- stats::dbeta(gr, mu * k, (1 - mu) * k); f <- f / sum(f); pc <- gr * f; pn <- (1 - gr) * f
    true_auc <- sum(pc * (cumsum(pn) - pn / 2)) / (sum(pc) * sum(pn)) }
  for (s in c("survey", "proxies", "both", "index", "ceiling", "district_others")) {
    if (s == "ceiling") {          # tau^2 / total: the most any district-constant predictor can reach
      # bootstrapped over DISTRICTS: resampling clusters duplicates them inside a district and biases a
      # between-district variance upward; a resampled district keeps its own label so duplicates stay separate
      idx_by_d <- split(seq_along(g), g)
      bs <- vapply(seq_len(B), function(b) { pk <- sample(length(idx_by_d), replace = TRUE); ii <- unlist(idx_by_d[pk], use.names = FALSE)
        mom_share71(ce$y[ii], rep(seq_along(pk), lengths(idx_by_d[pk])), tp) }, 0)
      est <- share; a <- ce; auc_s <- true_auc
    } else {
      a <- A[A$type == tp & A$outcome == o & A$set == (if (s == "district_others") "ceiling" else s), ]; a <- a[order(a$row), ]
      if (nrow(a) != nrow(ce) || any(!is.finite(a$pred))) next
      sk <- function(i) 1 - mean((a$y[i] - a$pred[i])^2) / mean((a$y[i] - a$null[i])^2)
      bs <- vapply(draws, sk, 0); est <- sk(seq_len(nrow(a))); auc_s <- if (tp == "bin") auc1(a$y, a$pred) else NA_real_
    }
    out[[length(out) + 1]] <- data.frame(country = COUNTRY71, type = tp, outcome = o, set = s, n = nrow(a), n_clusters = length(cl_ids),
      mean_outcome = mu, skill = est, lo = unname(stats::quantile(bs, 0.025)), hi = unname(stats::quantile(bs, 0.975)),
      auc = auc_s, between_district_share = share, true_prevalence_ceiling_auc = true_auc)
  }
}
M <- bind_rows(out) |> left_join(within_auc, by = c("type", "outcome", "set"))
write.csv(M, file.path(OUTDIR, "il02_honest_person_level.csv"), row.names = FALSE)
options(width = 200)
print(as.data.frame(M |> mutate(across(c(mean_outcome, skill, lo, hi, auc, auc_within, between_district_share, true_prevalence_ceiling_auc), ~ round(.x, 3)))), row.names = FALSE)
