# =============================================================================
# scripts/covariates/review_predictor_set.R   [RV-01, 2026-09-15]
#
# CRITICAL REVIEW OF THE SHARED PREDICTOR SET BEFORE A MODEL RE-RUN
#
# The audit (audit_predictor_set.R) reports provenance and coverage per column.
# This review asks the modelling questions:
#   1. what actually reaches the design matrix in each arm - per-country
#      prep_predictors_v2() (coverage >= 0.7, sd > 0) and, for transport, the
#      four-country intersection - by tier and domain
#   2. value sanity: shares outside [0, 1], negative rates, |skew| > 5, near-
#      constant columns, per-country missingness
#   3. redundancy: pooled rank-normalised |r| >= 0.98 pairs across sources
#   4. cross-source agreement: the same construct from two programmes (MICS vs
#      DHS WASH, wealth, stunting; MODIS NDVI vs EVI; HCES food share vs wealth)
#      as within-country Spearman - a low value flags a defect in one of them
#   5. broadcast structure: columns whose values are shared by many districts
#      of a country (parent-level estimates broadcast down)
#
#   Rscript -e "source('scripts/covariates/review_predictor_set.R')"
# -> results/tables/predictor_review_<date>.md (+ _pairs.csv, _design.csv)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
HDIR <- "data/covariates/harmonized"; ODIR <- "results/tables"; STAMP <- format(Sys.Date(), "%Y-%m-%d")
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
M <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
cols <- setdiff(names(S), c("country", "Admin1", "Admin2")); stopifnot(setequal(cols, M$column)); M <- M[match(cols, M$column), ]
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
md <- c(sprintf("# Predictor set review, %s", STAMP), "", sprintf("%d units x %d predictors. Tiers: %s.", nrow(S), length(cols), paste(sprintf("%s %d", names(table(M$tier)), table(M$tier)), collapse = ", ")), "")

# ── 1. what reaches the design matrix ────────────────────────────────────────
design <- list()
for (tiers in list(c("open"), c("open", "survey_public"), c("open", "survey_public", "survey_dhs"))) {
  preds <- withr::with_envvar(c(V2_PREDICTOR_TIERS = paste(tiers, collapse = ","), V2_KEEP_NATIONAL = ""), suppressMessages(drop_near_outcome_v2(cols, M)))
  kept <- lapply(COUNTRIES, function(cn) { X <- as.matrix(S[S$country == cn, preds, drop = FALSE]); colnames(prep_predictors_v2(X)) })
  names(kept) <- COUNTRIES
  common <- Reduce(intersect, kept)
  design[[paste(tiers, collapse = "+")]] <- data.frame(tiers = paste(tiers, collapse = "+"), declared = length(preds),
    Gambia = length(kept$Gambia), Ghana = length(kept$Ghana), Malawi = length(kept$Malawi), SierraLeone = length(kept$SierraLeone), transport_common = length(common),
    common_domains = n_distinct(M$domain[M$column %in% common]), stringsAsFactors = FALSE)
  if (identical(tiers, c("open", "survey_public"))) {
    dom_common <- M |> filter(column %in% common) |> count(domain, name = "columns") |> arrange(desc(columns))
    dropped_transport <- setdiff(preds, common)
    dd <- M |> filter(column %in% dropped_transport) |> mutate(why = vapply(column, function(v) { k <- vapply(COUNTRIES, function(cn) v %in% kept[[cn]], NA); paste(COUNTRIES[!k], collapse = ";") }, "")) |> count(source, why, name = "columns") |> arrange(desc(columns))
  }
}
DS <- bind_rows(design); write.csv(DS, file.path(ODIR, sprintf("predictor_review_%s_design.csv", STAMP)), row.names = FALSE)
md <- c(md, "## 1. Design matrix by arm (per-country prep: coverage >= 0.7 and sd > 0; transport = 4-country intersection)", "",
        "| tiers | declared (after policy) | Gambia | Ghana | Malawi | Sierra Leone | transport common | domains in common |", "|---|---:|---:|---:|---:|---:|---:|---:|",
        sprintf("| %s | %d | %d | %d | %d | %d | %d | %d |", DS$tiers, DS$declared, DS$Gambia, DS$Ghana, DS$Malawi, DS$SierraLeone, DS$transport_common, DS$common_domains), "",
        "Transport headline arm (open + survey_public): columns in the common matrix by domain", "", "| domain | columns |", "|---|---:|", sprintf("| %s | %d |", dom_common$domain, dom_common$columns), "",
        "Columns of the headline arm that do NOT reach the transport matrix (missing or constant in the country named)", "", "| source | absent in | columns |", "|---|---|---:|",
        sprintf("| %s | %s | %d |", dd$source, dd$why, dd$columns), "")

# ── 2. value sanity ──────────────────────────────────────────────────────────
share_like <- grepl("share|_pct|_prev|_frac|^dhs_c_|^dhs_w_|^dhs_hh_|^mics_|^hces_any|hces_food_share|hces_own_prod|_coverage|_rate$|^who_|^ihme_", cols) & !grepl("_km|_mean_haz|log_|score|_n$|_days|years|_kg|_g$|_mg|hhsize|_pae|volatility|rel_|inflation|range|_int|density|ratio|_ha$|_total$|_count|parity|age_|delivery|per10k|per_week|attainment|trend|_sd$|_min$|_max$", cols)
stats <- lapply(cols, function(v) { z <- S[[v]]; z <- z[is.finite(z)]; if (length(z) < 4) return(c(NA, NA, NA, NA, NA)); m <- mean(z); s <- sd(z); sk <- if (s > 0) mean((z - m)^3) / s^3 else 0; c(min(z), max(z), s / abs(m + 1e-9), sk, mean(is.finite(S[[v]]))) })
ST <- as.data.frame(do.call(rbind, stats)); names(ST) <- c("min", "max", "cv", "skew", "finite"); ST$column <- cols
ST$out_of_unit <- share_like & (ST$min < -1e-6 | ST$max > 1 + 1e-6)
ST$near_constant <- is.finite(ST$cv) & ST$cv < 0.01 & M$subnational
ST$heavy_skew <- is.finite(ST$skew) & abs(ST$skew) > 5
miss_c <- sapply(COUNTRIES, function(cn) vapply(cols, function(v) mean(!is.finite(S[[v]][S$country == cn])), 0))
ST$partial_missing <- apply(miss_c, 1, function(m) any(m > 0.05 & m < 0.95))   # a country with some but not all units missing
md <- c(md, "## 2. Value sanity", "",
        sprintf("- share-like columns outside [0, 1]: %d %s", sum(ST$out_of_unit), if (any(ST$out_of_unit)) paste0("(", paste(head(ST$column[ST$out_of_unit], 12), collapse = ", "), ")") else ""),
        sprintf("- near-constant subnational columns (CV < 1%%): %d %s", sum(ST$near_constant), if (any(ST$near_constant)) paste0("(", paste(head(ST$column[ST$near_constant], 12), collapse = ", "), ")") else ""),
        sprintf("- |skew| > 5 (rank-normalised before any fit, so cosmetic): %d", sum(ST$heavy_skew)),
        sprintf("- columns with partial missingness inside a country (5-95%% of units): %d %s", sum(ST$partial_missing), if (any(ST$partial_missing)) paste0("(", paste(head(ST$column[ST$partial_missing], 15), collapse = ", "), ")") else ""), "")

# ── 3. redundancy: pooled rank-normalised |r| >= 0.98 ───────────────────────
Xr <- do.call(rbind, lapply(COUNTRIES, function(cn) { X <- as.matrix(S[S$country == cn, cols, drop = FALSE]); Xn <- apply(X, 2, rank_normalize_v2); Xn }))
ok <- colSums(is.finite(Xr)) >= 100
C <- suppressWarnings(stats::cor(Xr[, ok], use = "pairwise.complete.obs"))
pairs <- which(abs(C) >= 0.98 & upper.tri(C), arr.ind = TRUE)
PR <- data.frame(a = rownames(C)[pairs[, 1]], b = colnames(C)[pairs[, 2]], r = round(C[pairs], 3), stringsAsFactors = FALSE)
PR$source_a <- M$source[match(PR$a, M$column)]; PR$source_b <- M$source[match(PR$b, M$column)]; PR$cross_source <- PR$source_a != PR$source_b
PR <- PR[order(-abs(PR$r)), ]; write.csv(PR, file.path(ODIR, sprintf("predictor_review_%s_pairs.csv", STAMP)), row.names = FALSE)
md <- c(md, "## 3. Redundancy (pooled within-country rank-normalised |r| >= 0.98)", "",
        sprintf("- %d pairs, %d across sources; by source pair:", nrow(PR), sum(PR$cross_source)), "",
        if (nrow(PR)) { t <- PR |> count(source_a, source_b, name = "pairs") |> arrange(desc(pairs)); sprintf("- %s x %s: %d", t$source_a, t$source_b, t$pairs) } else "- none", "",
        if (any(PR$cross_source)) c("Cross-source pairs:", "", sprintf("- `%s` ~ `%s` (r = %.3f)", PR$a[PR$cross_source], PR$b[PR$cross_source], PR$r[PR$cross_source])) else NULL, "")

# ── 4. cross-source agreement (within-country Spearman) ─────────────────────
agree <- list(
  c("mics_water_improved", "dhs_hh_improved_water"), c("mics_sanitation_improved", "dhs_hh_improved_sanitation"), c("mics_open_defecation", "dhs_hh_open_defecation"),
  c("mics_c_stunted", "dhs_c_stunted"), c("mics_c_stunted", "ihme_stuntingprevalence"), c("mics_wealth_score_mean", "dhs_hh_wealth_mean"), c("mics_wealth_score_mean", "rwi_mean"),
  c("mics_w_secondary_plus", "dhs_w_secondary_plus"), c("mics_w_literate", "dhs_w_literate"), c("mics_c_diarrhoea_2wk", "dhs_c_diarrhea_2wk"), c("mics_c_fever_2wk", "dhs_c_fever_2wk"),
  c("mics_c_mdd", "dhs_c_mdd_4plus"), c("mics_electricity", "dhs_hh_electricity"), c("mics_heat_vmsl_sy", "dhs_c_measles1"),
  c("ndvi_modis_t0", "evi_t0"), c("ndvi_modis_win", "ndvi_modis_t0"), c("hces_food_share", "dhs_hh_wealth_mean"), c("hces_log_cons_pae_rel", "rwi_mean"), c("hces_hdds", "dhs_c_diet_diversity_score"),
  c("map_sy_pf_parasite_rate", "ihme_malaria_pfpr"), c("map_sy_pf_incidence_rate", "ihme_malaria_incidence_rate"), c("wpop_log_density_survey_year", "ghs_pop"), c("rtfp_fpi_rel_national", "fprice_staple_rel"))
AG <- bind_rows(lapply(agree, function(p) { if (!all(p %in% cols)) return(NULL)
  r <- vapply(COUNTRIES, function(cn) { d <- S[S$country == cn, p]; ok <- is.finite(d[[1]]) & is.finite(d[[2]]); if (sum(ok) < 8) NA_real_ else suppressWarnings(stats::cor(d[[1]][ok], d[[2]][ok], method = "spearman")) }, 0)
  data.frame(a = p[1], b = p[2], Gambia = round(r[1], 2), Ghana = round(r[2], 2), Malawi = round(r[3], 2), SierraLeone = round(r[4], 2), stringsAsFactors = FALSE) }))
md <- c(md, "## 4. Cross-source agreement (within-country Spearman; NA = fewer than 8 units with both)", "", "| a | b | Gambia | Ghana | Malawi | Sierra Leone |", "|---|---|---:|---:|---:|---:|",
        sprintf("| `%s` | `%s` | %s | %s | %s | %s |", AG$a, AG$b, AG$Gambia, AG$Ghana, AG$Malawi, AG$SierraLeone), "")

# ── 5. broadcast structure ───────────────────────────────────────────────────
bc <- sapply(COUNTRIES, function(cn) vapply(cols, function(v) { z <- S[[v]][S$country == cn]; z <- z[is.finite(z)]; if (length(z) < 8) NA_real_ else length(unique(round(z, 8))) / length(z) }, 0))
BC <- data.frame(column = cols, source = M$source, tier = M$tier, distinct_share_min = apply(bc, 1, function(x) suppressWarnings(min(x, na.rm = TRUE))), stringsAsFactors = FALSE)
BC$distinct_share_min[!is.finite(BC$distinct_share_min)] <- NA
bsum <- BC |> filter(!is.na(distinct_share_min), distinct_share_min < 0.5) |> count(source, name = "columns") |> arrange(desc(columns))
md <- c(md, "## 5. Broadcast structure (columns where, in some country, fewer than half the units carry distinct values: parent-level estimates broadcast to districts)", "",
        "| source | columns |", "|---|---:|", sprintf("| %s | %d |", bsum$source, bsum$columns), "")

writeLines(md, file.path(ODIR, sprintf("predictor_review_%s.md", STAMP)))
cat(paste(md, collapse = "\n")); cat("\n-> ", file.path(ODIR, sprintf("predictor_review_%s.md", STAMP)), "\nDONE\n")
