# =============================================================================
# scripts/external_validation/04_transport_test.R              [XV-01, XV-02]
#
# DOES THE CLIMATE + SOIL INDEX RANK THE SUB-NATIONAL UNITS OF A COUNTRY IT
# HAS NEVER SEEN, AGAINST LABELS THIS PROJECT DID NOT MAKE?
#
# The transport claim (P1) has only ever been scored leave-one-country-out
# INSIDE the four-country panel -- same harmonisation, same BRINDA adjustment,
# our own weighting. This scores it against WHO VMNIS sub-national deposits for
# six countries that are not in the panel and never enter training:
#
#   XV-01, Africa       Zambia 2023 (9 provinces), Ethiopia 2015 (11 regions),
#                       Sudan 2018 (15 states), Nigeria 2021 (6 NFCMS zones)
#   XV-02, off-continent Pakistan 2018 (8 provinces), India CNNS 2016-18
#                       (29 states; CNNS sampled no women, so child outcomes
#                       only). Both deposit prevalence but no biomarker means,
#                       so the off-continent arm scores PREVALENCE only.
#
# TWO SOIL BLOCKS, AND WHY BOTH RUN ON AFRICA. The pre-registered index uses
# iSDAsoil, which is Africa-only, so the off-continent arm must substitute
# global SoilGrids. SoilGrids carries no plant-available micronutrients (Zn,
# Fe, Ca, Mg, P, K, S) -- mechanistically the interesting half of iSDA for a
# micronutrient outcome -- so the substitution is not neutral. Running BOTH
# blocks on the four African countries prices it on the same cells, so a weak
# off-continent result can be read against that price instead of being
# confounded with it. `XV_SOIL` restricts to one block if wanted.
#
# WHY A RANKING TEST AND NOT A LEVEL TEST. Cross-survey biomarker LEVELS are
# not comparable (raw ferritin spans 6x across our own four surveys;
# fe_transport_level_offset), and VMNIS cut-offs and adjustments vary by
# deposit. A WITHIN-country Spearman is invariant to any country-constant
# offset and to any assay or cut-off choice held fixed inside one survey.
#
# SIGN. targets_v2.csv's `y_level` rises with DEFICIENCY (its correlation with
# y_prev is +0.6 to +0.95 over all 24 panel cells), i.e. the biomarker is
# negated. VMNIS reports the biomarker itself, so its level is negated here to
# match. Getting this backwards silently inverts every result
# (signal_scan_sign_convention).
#
# PROTOCOL, following 16_admin1_transport.R: predictors area-averaged as a
# SIMPLE MEAN OF DISTRICTS, rank-normalised WITHIN country, domain PCs to 80%
# variance oriented from the TRAINING rows only; outcome z-scored within
# country. Training is all four panel countries at their admin-1 rung.
#
#   Rscript scripts/external_validation/04_transport_test.R
# -> results/tables/external_validation/xv_transport.csv
#    results/tables/external_validation/xv_transport_summary.csv
#    results/tables/external_validation/xv_transport_pooled.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

XV <- "data/external_validation"
OUTDIR <- "results/tables/external_validation"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20260922L)
NPERM <- 2000L

PANEL <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
KEY <- c("country", "Admin1", "Admin2")

TG <- read.csv("results/tables/protocol_v2/targets_v2.csv", stringsAsFactors = FALSE)
VM <- read.csv(file.path(XV, "vmnis_admin1_targets.csv"), stringsAsFactors = FALSE)
XW <- read.csv(file.path(XV, "gadm_to_vmnis_crosswalk.csv"), stringsAsFactors = FALSE)

CL <- read.csv(file.path(XV, "gee_clim_admin2.csv"), check.names = FALSE)
clim_cols <- setdiff(names(CL), KEY)

SOILS <- list(
  isda  = list(file = "gee_isda_admin2.csv",  label = "iSDAsoil (Africa only)"),
  sgrid = list(file = "gee_sgrid_admin2.csv", label = "SoilGrids v2.0 (global)"))
want <- Sys.getenv("XV_SOIL", "")
if (nzchar(want)) SOILS <- SOILS[strsplit(want, ",")[[1]]]

rows <- list(); nulls <- list(); cells <- list()

for (soil in names(SOILS)) {
  sf_path <- file.path(XV, SOILS[[soil]]$file)
  if (!file.exists(sf_path)) { cat("skip", soil, "- not extracted\n"); next }
  IS <- read.csv(sf_path, check.names = FALSE)
  soil_cols <- setdiff(names(IS), KEY)
  PREDS <- c(clim_cols, soil_cols)
  COV <- dplyr::inner_join(CL, IS, by = KEY)
  # GADM carries some districts as several polygons under one name (19 rows in
  # India: Arunachal Pradesh and Himachal Pradesh appear 2-3 times over disputed
  # borders). Left as is, the crosswalk join fans out and those districts count
  # 2-3x in their state's mean. Collapse to one row per named district first, so
  # "simple mean of districts" means what it says.
  COV <- COV |> group_by(across(all_of(KEY))) |>
    summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    as.data.frame()
  XWu <- unique(XW[, c("country", "Admin1", "Admin2", "vmnis_unit")])
  domain_of <- stats::setNames(
    c(rep("Climate and weather", length(clim_cols)),
      rep("Soil characteristics", length(soil_cols))), PREDS)
  # iSDA has no off-continent coverage, so those countries simply are not in
  # COV and drop out of NEW by the inner join below
  NEW <- intersect(unique(VM$country), unique(COV$country))
  cat(sprintf("\n[%s] %s: %d climate + %d soil cols, %d districts, test: %s\n",
              soil, SOILS[[soil]]$label, length(clim_cols), length(soil_cols),
              nrow(COV), paste(NEW, collapse = ", ")))

  build_panel <- function(cn, on, target) {
    t <- TG[TG$country == cn & TG$outcome == on, ]
    ycol <- if (target == "prev") "y_prev" else "y_level"
    wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
    t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]
    if (!nrow(t)) return(NULL)
    a1 <- t |> group_by(unit = Admin1) |>
      summarise(y = stats::weighted.mean(.data[[ycol]], .data[[wcol]]),
                w = sum(.data[[wcol]]), .groups = "drop")
    x1 <- COV[COV$country == cn, ] |> group_by(unit = Admin1) |>
      summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    m <- dplyr::inner_join(a1, x1, by = "unit")
    if (nrow(m) < 3) return(NULL)
    list(country = cn, unit = m$unit, y = m$y, w = m$w,
         X = as.matrix(m[, PREDS, drop = FALSE]))
  }

  build_new <- function(cn, on, target) {
    v <- VM[VM$country == cn & VM$outcome == on, ]
    if (!nrow(v)) return(NULL)
    v$y <- if (target == "prev") v$prev / 100 else -v$level
    v <- v[is.finite(v$y), c("unit", "y")]
    if (nrow(v) < 5) return(NULL)
    cv <- COV[COV$country == cn, ] |>
      dplyr::inner_join(XWu, by = KEY, relationship = "many-to-one")
    x1 <- cv |> group_by(unit = vmnis_unit) |>
      summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    m <- dplyr::inner_join(v, x1, by = "unit")
    if (nrow(m) < 5) return(NULL)
    list(country = cn, unit = m$unit, y = m$y, w = rep(1, nrow(m)),
         X = as.matrix(m[, PREDS, drop = FALSE]))
  }

  for (target in c("level", "prev")) {
    for (on in sort(unique(VM$outcome))) {
      tr_list <- list()
      for (cn in PANEL) {
        z <- tryCatch(build_panel(cn, on, target), error = function(e) NULL)
        if (!is.null(z)) tr_list[[cn]] <- z
      }
      if (length(tr_list) < 2) next
      for (cn in NEW) {
        z <- tryCatch(build_new(cn, on, target), error = function(e) NULL)
        if (is.null(z)) next
        cl <- c(tr_list, stats::setNames(list(z), cn))
        Xs <- lapply(cl, function(q) prep_predictors_v2(q$X))
        common <- Reduce(intersect, lapply(Xs, colnames))
        if (length(common) < 20) next
        Xm <- do.call(rbind, lapply(Xs, function(q) q[, common, drop = FALSE]))
        Y <- unlist(lapply(cl, function(q) as.numeric(scale(q$y))))
        ctry <- rep(names(cl), vapply(cl, function(q) length(q$y), 0L))
        wv <- unlist(lapply(cl, function(q) q$w))
        te <- which(ctry == cn); tr <- which(ctry != cn)
        Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        aux <- list(lon = rep(0, length(Y)), lat = rep(0, length(Y)), y_nat = Y)
        for (a in c("null_train_mean", "domain_index", "domain_enet")) {
          p <- tryCatch(ARMS_V2[[a]](tr, te, Y, Xm, Dm, aux),
                        error = function(e) rep(NA_real_, length(te)))
          if (length(p) != length(te)) p <- rep(NA_real_, length(te))
          s <- score_v2(Y[te], p, wv[te],
                        scale = if (target == "prev") "prev" else "level")
          pnull <- NA_real_; pval <- NA_real_
          k <- paste(soil, target, a, cn, on, sep = "|")
          if (is.finite(s$spearman)) {
            nullv <- replicate(NPERM, suppressWarnings(
              stats::cor(sample(Y[te]), p, method = "spearman")))
            pnull <- stats::quantile(nullv, 0.95, na.rm = TRUE)
            pval <- (1 + sum(nullv >= s$spearman, na.rm = TRUE)) / (NPERM + 1)
            nulls[[k]] <- nullv
            cells[[k]] <- list(country = cn, unit = z$unit, obs = Y[te], pred = p)
          }
          rows[[length(rows) + 1L]] <- data.frame(
            soil = soil, country = cn, outcome = on, target = target, arm = a,
            arm_group = VM$arm[match(cn, VM$country)],
            n_units = length(te), n_train = length(tr), n_cols = length(common),
            spearman = s$spearman, topk = s$topk,
            null_p95 = unname(pnull), perm_p = pval,
            thin = length(te) < 8, stringsAsFactors = FALSE)
        }
      }
    }
    cat("  ", soil, target, "done\n")
  }
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "xv_transport.csv"), row.names = FALSE)

#' Pooled permutation test over a set of cells.
#'
#' A single cell of 6-15 units cannot clear its own null (at n = 9 the 95th
#' percentile sits near 0.55), so the claim is made on the cells together.
#' TWO NULLS, because they assume different things. The INDEPENDENT null draws
#' one permutation per cell; it is too generous, because a country's outcomes
#' are measured on the same units and are correlated (NC-01). The COUNTRY-BLOCK
#' null draws ONE relabelling of a country's units per replicate and applies it
#' to every outcome of that country, preserving that correlation. The
#' country-block p is the one to quote.
pooled <- function(d, keys, label) {
  keys <- keys[keys %in% names(nulls)]
  if (length(keys) < 2 || !nrow(d)) return(NULL)
  obs <- mean(d$spearman); pos <- sum(d$spearman > 0); n <- nrow(d)
  p_sign <- stats::binom.test(pos, n, 0.5, alternative = "greater")$p.value
  nd_ind <- rowMeans(do.call(cbind, nulls[keys]), na.rm = TRUE)
  p_ind <- (1 + sum(nd_ind >= obs, na.rm = TRUE)) / (length(nd_ind) + 1)
  cc <- cells[keys]
  by_country <- split(seq_along(cc), vapply(cc, function(q) q$country, ""))
  nd_blk <- numeric(NPERM)
  for (r in seq_len(NPERM)) {
    rho <- numeric(0)
    for (cn in names(by_country)) {
      u <- unique(unlist(lapply(cc[by_country[[cn]]], function(q) q$unit)))
      # one random priority over the country's units, shared by all of its
      # outcomes; order() restricted to a cell keeps this valid when a cell
      # covers only some units (Zambia's women_vitA has 6 of 9)
      pri <- stats::setNames(sample(length(u)), u)
      for (i in by_country[[cn]]) {
        q <- cc[[i]]
        rho <- c(rho, suppressWarnings(stats::cor(
          q$obs[order(pri[q$unit])], q$pred, method = "spearman")))
      }
    }
    nd_blk[r] <- mean(rho, na.rm = TRUE)
  }
  p_blk <- (1 + sum(nd_blk >= obs, na.rm = TRUE)) / (NPERM + 1)
  cat(sprintf("%-34s cells=%2d  mean rho=%+.3f  positive %2d/%2d\n",
              label, n, obs, pos, n))
  cat(sprintf("%-34s   country-block null95=%.3f p=%.4f | indep p=%.4f | sign p=%.5f\n",
              "", stats::quantile(nd_blk, 0.95), p_blk, p_ind, p_sign))
  data.frame(set = label, cells = n, mean_rho = obs, positive = pos,
             block_null_p95 = unname(stats::quantile(nd_blk, 0.95)),
             block_p = p_blk, indep_p = p_ind, sign_p = p_sign,
             stringsAsFactors = FALSE)
}

cat("\n===== PER-CELL SUMMARY =====\n")
print(as.data.frame(R |> filter(arm != "null_train_mean", is.finite(spearman)) |>
  group_by(soil, arm_group, target, arm) |>
  summarise(cells = dplyr::n(), mean_rho = round(mean(spearman), 3),
            positive = sum(spearman > 0), .groups = "drop")), row.names = FALSE)

cat("\n===== POOLED =====\n")
pool <- list()
idx <- R |> filter(arm == "domain_index", is.finite(spearman))
# (a) the pre-registered soil block on the Africa arm
for (tg in c("level", "prev")) {
  d <- idx |> filter(soil == "isda", arm_group == "africa", target == tg)
  pool[[length(pool) + 1L]] <- pooled(d, paste("isda", tg, "domain_index",
    d$country, d$outcome, sep = "|"), sprintf("Africa / iSDA / %s", tg))
}
# (b) the same cells with the substituted block -- the price of substituting
for (tg in c("level", "prev")) {
  d <- idx |> filter(soil == "sgrid", arm_group == "africa", target == tg)
  pool[[length(pool) + 1L]] <- pooled(d, paste("sgrid", tg, "domain_index",
    d$country, d$outcome, sep = "|"), sprintf("Africa / SoilGrids / %s", tg))
}
# (c) the off-continent test, which can only use the substituted block
d <- idx |> filter(soil == "sgrid", arm_group == "offcontinent", target == "prev")
pool[[length(pool) + 1L]] <- pooled(d, paste("sgrid", "prev", "domain_index",
  d$country, d$outcome, sep = "|"), "Off-continent / SoilGrids / prev")
for (cn in sort(unique(d$country))) {
  dd <- d |> filter(country == cn)
  pool[[length(pool) + 1L]] <- pooled(dd, paste("sgrid", "prev", "domain_index",
    dd$country, dd$outcome, sep = "|"), paste0("   ", cn, " / prev"))
}

pool <- bind_rows(pool)
if (nrow(pool)) write.csv(pool, file.path(OUTDIR, "xv_transport_pooled.csv"),
                          row.names = FALSE)
write.csv(R |> filter(arm != "null_train_mean", is.finite(spearman)) |>
            group_by(soil, target, arm, country) |>
            summarise(cells = dplyr::n(), mean_rho = mean(spearman),
                      positive = sum(spearman > 0), .groups = "drop"),
          file.path(OUTDIR, "xv_transport_summary.csv"), row.names = FALSE)
cat("\nwrote xv_transport.csv, xv_transport_summary.csv, xv_transport_pooled.csv\n")
