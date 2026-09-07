# =============================================================================
# scripts/cluster_level/06_matched_vocabulary.R   [CL-06]
#
# MATCHED VOCABULARY: is the cluster track's deficit the fitting unit, or the
# layers it sees?
#
# The cluster track fits on ~150 buffer covariates; the district track on the
# 473-column shared set. Comparing them confounds the unit with the
# vocabulary. Here the SAME district-level protocol arms are fitted on three
# vocabularies, with the same district targets (targets_v2), the same
# districts (those holding at least one GPS cluster) and the same folds:
#   full         the shared Admin-2 set (473 columns, 24 domains)
#   cluster_agg  the cluster buffer covariates averaged over the clusters in
#                each district (transportable role only), with the cluster
#                metadata's domain labels: the cluster vocabulary, fitted at
#                the district
#   both         the union
# Read next to the cluster-fitted models aggregated to districts (rho_agg in
# benchmarks_cluster_cells / _loco): if cluster_agg at the district matches
# the full set, the layers are not the problem and the fitting unit is; if
# cluster_agg trails the full set by the same margin the cluster fit does, the
# vocabulary explains the gap.
#
# Estimands: in-fill (5-fold by district x 10 draws), leave-one-region-out,
# leave-one-country-out scored at Admin-2 (pooled, within-country rank-
# normalised, the script-43 recipe). Arms: the protocol's zero-tuning domain
# index (estimator of record) plus null, jackknifed regional mean, spatial
# smoother and spatial + domain for the in-country estimands.
#
#   Rscript scripts/cluster_level/06_matched_vocabulary.R
# -> results/tables/cluster_level/matched_vocabulary_cells.csv, _summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/cluster_level"; CDIR <- "data/covariates/cluster"; HDIR <- "data/covariates/harmonized"
REPS <- as.integer(Sys.getenv("CL_REPS", "10")); MIN_TRAIN <- 20L; set.seed(20260907L)
PT <- Sys.getenv("CL_PRED_TAG", "")

TG <- read.csv("results/tables/protocol_v2/targets_v2.csv", stringsAsFactors = FALSE)
TC <- read.csv(file.path(OUTDIR, "targets_cluster.csv"), stringsAsFactors = FALSE); TC$cluster <- as.character(TC$cluster)
S  <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE); SM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
P  <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, ".csv")), check.names = FALSE, stringsAsFactors = FALSE); P$cluster <- as.character(P$cluster)
PM <- read.csv(file.path(CDIR, paste0("predictors_cluster", PT, "_metadata.csv")), stringsAsFactors = FALSE); PM <- PM[PM$column %in% names(P), ]
FULL_COLS <- drop_near_outcome_v2(intersect(SM$column, names(S)), SM)
CL_COLS <- PM$column[PM$role != "fieldwork"]
domain_full <- stats::setNames(SM$domain, SM$column); domain_cl <- stats::setNames(PM$domain, PM$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds"); LC <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
cent <- function(cn) { b <- BND[[LC[[cn]]]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }

# cluster covariates averaged to the district (unweighted over clusters; the district "as its clusters")
CLAGG <- P |> distinct(country, cluster, .keep_all = TRUE) |> inner_join(TC |> distinct(country, cluster, Admin1, Admin2), by = c("country", "cluster")) |>
  group_by(country, Admin1, Admin2) |> summarise(n_clusters = dplyr::n(), across(all_of(CL_COLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
CLAGG[CL_COLS] <- lapply(CLAGG[CL_COLS], function(v) { v[!is.finite(v)] <- NA; v })
cat(sprintf("cluster vocabulary: %d columns aggregated to %d districts (%s)\n", length(CL_COLS), nrow(CLAGG), paste(sprintf("%s %d", names(table(CLAGG$country)), table(CLAGG$country)), collapse = ", ")))
cat(sprintf("full vocabulary: %d columns\n", length(FULL_COLS)))
VOCABS <- list(full = list(cols = FULL_COLS, dom = domain_full), cluster_agg = list(cols = CL_COLS, dom = domain_cl),
               both = list(cols = c(FULL_COLS, CL_COLS), dom = c(domain_full, domain_cl)))
# the cluster columns carry their own names; a name shared with the district set would collide in "both", so prefix on the way in
names_cl <- stats::setNames(paste0("cl__", CL_COLS), CL_COLS)
CLAGG2 <- CLAGG; names(CLAGG2)[match(CL_COLS, names(CLAGG2))] <- names_cl[CL_COLS]
VOCABS$cluster_agg$cols <- unname(names_cl); VOCABS$cluster_agg$dom <- stats::setNames(domain_cl[CL_COLS], names_cl[CL_COLS])
VOCABS$both$cols <- c(FULL_COLS, unname(names_cl)); VOCABS$both$dom <- c(domain_full, VOCABS$cluster_agg$dom)

build_cell <- function(cn, on, target, vocab) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  m <- t |> inner_join(CLAGG2[CLAGG2$country == cn, c("Admin1", "Admin2", unname(names_cl))], by = c("Admin1", "Admin2")) |>
    inner_join(S[S$country == cn, c("Admin1", "Admin2", FULL_COLS)], by = c("Admin1", "Admin2")) |> inner_join(cent(cn), by = c("Admin1", "Admin2"))
  m <- m[is.finite(m$lon), ]; if (nrow(m) < 12) return(NULL)
  cols <- VOCABS[[vocab]]$cols; Xr <- prep_predictors_v2(as.matrix(m[, cols, drop = FALSE])); if (ncol(Xr) < 10) return(NULL)
  yn <- m[[ycol]]; list(country = cn, n = nrow(m), y_nat = yn, y_mod = if (target == "prev") .v2_logit(yn) else yn, X = Xr, w = m[[wcol]], Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone"); OUTCOMES <- unique(TG$outcome)
IN_ARMS <- c("null_train_mean", "region_mean_jk", "spatial", "domain_index", "spatial_plus_domain")
rows <- list()
for (target in c("level", "prev")) for (vocab in names(VOCABS)) { dom <- VOCABS[[vocab]]$dom
  # ── in-fill and region, per country ──
  for (cn in COUNTRIES) for (on in OUTCOMES) { z <- tryCatch(build_cell(cn, on, target, vocab), error = function(e) NULL); if (is.null(z)) next
    Y <- z$y_mod; n <- z$n; aux <- list(lon = z$lon, lat = z$lat, Admin1 = z$Admin1, y_nat = Y); sc <- if (target == "prev") "prev" else "level"
    for (est in c("infill", "region")) { if (est == "infill" && n < 15) next; if (est == "region" && dplyr::n_distinct(z$Admin1) < 3) next
      for (arm in IN_ARMS) { if (arm == "region_mean_jk" && est == "region") next
        sp <- c(); for (r in seq_len(if (est == "infill") REPS else 1L)) { folds <- if (est == "infill") make_folds_v2("kfold_district", n, k = 5, rep_id = r) else as.integer(factor(z$Admin1)); pred <- rep(NA_real_, n)
          for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 8) next
            D <- domain_representation_v2(z$X, dom, sign_rows = tr); p <- tryCatch(ARMS_V2[[arm]](tr, te, Y, z$X, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) pred[te] <- p }
          sp <- c(sp, score_v2(z$y_nat, pred, z$w, scale = sc)$spearman) }
        rows[[length(rows) + 1L]] <- data.frame(vocab = vocab, target = target, estimand = est, country = cn, outcome = on, arm = arm, n_units = n, spearman = median(sp, na.rm = TRUE), stringsAsFactors = FALSE) } }
    cat(sprintf("  %-11s %-5s %-12s %-12s n=%3d done\n", vocab, target, cn, on, n)) }
  # ── leave-one-country-out at Admin-2, pooled (script-43 recipe) ──
  for (on in OUTCOMES) { cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target, vocab), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
    if (length(cl) < 3) next; common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 10) next
    Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
    aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y); folds <- as.integer(factor(ctry)); pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < MIN_TRAIN) next
      Dm <- domain_representation_v2(Xm, dom, sign_rows = tr); if (!ncol(Dm)) next
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) pred[te] <- p }
    for (cn in names(cl)) { k <- which(ctry == cn); rows[[length(rows) + 1L]] <- data.frame(vocab = vocab, target = target, estimand = "country", country = cn, outcome = on, arm = "domain_index", n_units = length(k),
      spearman = score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")$spearman, stringsAsFactors = FALSE) }
    cat(sprintf("  %-11s %-5s LOCO %-12s done (%d cols)\n", vocab, target, on, length(common))) } }
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, paste0("matched_vocabulary_cells", PT, ".csv")), row.names = FALSE)

# ── the comparison: district fits on three vocabularies, next to the cluster fits aggregated to districts ──
CC <- tryCatch(read.csv(file.path(OUTDIR, paste0(if (nzchar(PT)) paste0(sub("^_", "", PT), "/") else "", "benchmarks_cluster_cells.csv")), stringsAsFactors = FALSE), error = function(e) NULL)
CLO <- tryCatch(read.csv(file.path(OUTDIR, paste0(if (nzchar(PT)) paste0(sub("^_", "", PT), "/") else "", "benchmarks_cluster_loco.csv")), stringsAsFactors = FALSE), error = function(e) NULL)
S1 <- R |> filter(arm == "domain_index") |> group_by(vocab, target, estimand) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3), median_rho = round(median(spearman, na.rm = TRUE), 3), positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop")
if (!is.null(CC) && "rho_agg" %in% names(CC)) { cc <- CC |> filter(set == "transportable", arm == "domain_index", estimand %in% c("infill", "region")) |> group_by(target, estimand) |>
    summarise(cells = dplyr::n(), mean_rho = round(mean(rho_agg, na.rm = TRUE), 3), median_rho = round(median(rho_agg, na.rm = TRUE), 3), positive = sum(rho_agg > 0, na.rm = TRUE), .groups = "drop") |> mutate(vocab = "cluster_fitted_agg")
  S1 <- bind_rows(S1, cc) }
if (!is.null(CLO)) { rc <- intersect(c("rho_admin2", "rho_agg", "rho_a2"), names(CLO)); if (length(rc)) { lo <- CLO |> filter(set == "transportable", arm == "domain_index") |> group_by(target) |>
    summarise(cells = dplyr::n(), mean_rho = round(mean(.data[[rc[1]]], na.rm = TRUE), 3), median_rho = round(median(.data[[rc[1]]], na.rm = TRUE), 3), positive = sum(.data[[rc[1]]] > 0, na.rm = TRUE), .groups = "drop") |> mutate(vocab = "cluster_fitted_agg", estimand = "country")
  S1 <- bind_rows(S1, lo) } else cat("cluster LOCO file has columns:", paste(names(CLO), collapse = ", "), "\n") }
S1 <- S1 |> arrange(target, estimand, vocab); write.csv(S1, file.path(OUTDIR, paste0("matched_vocabulary_summary", PT, ".csv")), row.names = FALSE)
cat("\n===== CL-06: matched vocabulary (zero-tuning domain index; district-fitted on three vocabularies vs cluster-fitted aggregated) =====\n")
print(as.data.frame(S1), row.names = FALSE)
W <- S1 |> select(target, estimand, vocab, mean_rho) |> pivot_wider(names_from = vocab, values_from = mean_rho)
cat("\n-- gap decomposition (mean rho): full - cluster_agg = vocabulary; cluster_agg - cluster_fitted_agg = fitting unit --\n"); print(as.data.frame(W), row.names = FALSE)
cat("\nDONE\n")
