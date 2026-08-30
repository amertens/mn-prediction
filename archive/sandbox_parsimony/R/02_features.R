# =============================================================================
# sandbox_parsimony/R/02_features.R
#
# Parsimonious predictor sets for the Admin-2 proxy models.
#
# Motivation: the cross-country common predictor pool is NOT 237 independent
# measurements. For child_vitA it is ~99 SoilGrids columns (10 nutrients x
# ~9 near-identical summaries of the SAME static raster), ~108 IHME columns
# (mostly HIV counts, which are population-size proxies), and ~27 MAP columns.
# Correlation screening "top-K" on that pool spends its whole budget on
# duplicates of one or two underlying constructs.
# =============================================================================

# ── Domain map: regex -> conceptual family ──────────────────────────────────
DOMAIN_MAP <- list(
  soil_chem   = "^gee_soil(aluminium|calcium|cec|iron|magnesium|nitrogen|phosphorus|potassium|sulfur|totalcarbon|zinc)",
  vegetation  = "^gee_(ndvi|evi|productivity|wapor|lai|fpar)",
  climate     = "^gee_(terraclimate|trmm|chirps|era5|fldas|prec|temp|lst)",
  topography  = "^gee_(elevation|slope|aspect|srtm|topo)",
  urban_pop   = "^gee_(popdensity|gpw|ghsl|ghsbuilts|wsf|ccnl|ntl|viirs|globalhumanmodification|accessibility|traveltime)",
  cropland    = "^gee_(esa_worldcereal|cropland|spam|glc)",
  malaria     = "^MAP_[0-9]+_(Pf|Pv)_(Parasite_Rate|Incidence_Rate|Mortality_Rate|Reproductive_Number)",
  mal_interv  = "^MAP_[0-9]+_(Insecticide|IRS|Antimalarial)",
  blood_gen   = "^MAP_[0-9]+_(Sickle|HbC|G6PD)",
  ihme_health = "^ihme_",
  food_sec    = "^fsec_"
)

domain_of <- function(v) {
  out <- rep("other", length(v))
  for (d in names(DOMAIN_MAP)) out[out == "other" & grepl(DOMAIN_MAP[[d]], v)] <- d
  out
}

# ── Set 1: CURATED — one variable per distinct mechanism ────────────────────
# Chosen a priori from the causal story, not from the data:
#   diet quality / agro-ecology -> vegetation, climate, cropland, soil fertility
#   market access & poverty     -> nightlights, built-up, population, remoteness
#   inflammation / blood loss   -> malaria burden + nets
#   genetic anaemia             -> HbS / HbC / G6PD (mechanistic for iron/anaemia)
CURATED_PATTERNS <- c(
  ndvi           = "^gee_ndvi$",
  elevation      = "^gee_elevation$",
  pop_density    = "^gee_popdensity$",
  nightlights    = "^gee_ccnl$",
  human_mod      = "^gee_globalhumanmodification$",
  settlement_deg = "^gee_ghsl_smod$",
  built_surface  = "^gee_wsf$",
  precip         = "^gee_trmm$",
  soil_carbon    = "^gee_soiltotalcarbon_mean_0_20$",
  soil_zinc      = "^gee_soilzinc_mean_0_20$",
  soil_iron      = "^gee_soiliron_mean_0_20$",
  soil_nitrogen  = "^gee_soilnitrogen_mean_0_20$",
  malaria_pfpr   = "MAP_[0-9]+_Pf_Parasite_Rate$",
  itn_use        = "MAP_[0-9]+_Insecticide_Treated_Net_Use$",
  hbs_allele     = "MAP_[0-9]+_Sickle_Haemoglobin_HbS_Allele_Frequency$",
  g6pd_allele    = "MAP_[0-9]+_G6PDd_Allele_Frequency$"
)

curated_vars <- function(all_vars) {
  hits <- lapply(CURATED_PATTERNS, function(p) grep(p, all_vars, value = TRUE))
  # take the first match per construct (patterns are anchored enough to be ~unique)
  unlist(lapply(hits, function(h) if (length(h)) h[1] else NULL), use.names = FALSE)
}

# ── Set 2: DOMAIN PCs — first k PCs within each conceptual family ───────────
#' @return list(X = matrix of PC scores, rotate = function(newdata) -> scores)
domain_pcs <- function(Xtr, Xte, vars, npc = 1L) {
  dom <- domain_of(vars)
  tr_out <- list(); te_out <- list()
  for (d in unique(dom)) {
    cols <- vars[dom == d]
    if (length(cols) == 0) next
    A <- Xtr[, cols, drop = FALSE]; B <- Xte[, cols, drop = FALSE]
    sdv <- apply(A, 2, stats::sd); keep <- is.finite(sdv) & sdv > 1e-8
    if (!any(keep)) next
    A <- A[, keep, drop = FALSE]; B <- B[, keep, drop = FALSE]
    if (ncol(A) == 1) {
      tr_out[[d]] <- scale(A); te_out[[d]] <- scale(B, attr(tr_out[[d]], "scaled:center"),
                                                    attr(tr_out[[d]], "scaled:scale"))
      colnames(tr_out[[d]]) <- colnames(te_out[[d]]) <- d
      next
    }
    pr <- stats::prcomp(A, center = TRUE, scale. = TRUE)
    k <- min(npc, ncol(pr$rotation))
    tr_out[[d]] <- stats::predict(pr, A)[, seq_len(k), drop = FALSE]
    te_out[[d]] <- stats::predict(pr, B)[, seq_len(k), drop = FALSE]
    colnames(tr_out[[d]]) <- colnames(te_out[[d]]) <- paste0(d, "_pc", seq_len(k))
  }
  if (!length(tr_out)) return(NULL)
  list(Xtr = do.call(cbind, tr_out), Xte = do.call(cbind, te_out))
}

# ── Set 3: DECORRELATED REPRESENTATIVES — hclust on 1-|r|, cut at k ─────────
# Data-driven parsimony that, unlike top-K screening, cannot spend its whole
# budget on K near-duplicates of one construct. Cluster membership uses ONLY
# the training X (no outcome), so it is safe inside CV.
decorr_reps <- function(Xtr, vars, k = 20L, y = NULL) {
  A <- Xtr[, vars, drop = FALSE]
  sdv <- apply(A, 2, stats::sd); keep <- is.finite(sdv) & sdv > 1e-8
  A <- A[, keep, drop = FALSE]; vv <- vars[keep]
  if (ncol(A) <= k) return(vv)
  R <- suppressWarnings(stats::cor(A)); R[!is.finite(R)] <- 0
  hc <- stats::hclust(stats::as.dist(1 - abs(R)), method = "average")
  cl <- stats::cutree(hc, k = min(k, ncol(A)))
  # representative = member most correlated with the cluster's own first PC
  # (outcome-free), so the pick is stable across folds
  vapply(sort(unique(cl)), function(g) {
    idx <- which(cl == g)
    if (length(idx) == 1) return(vv[idx])
    S <- A[, idx, drop = FALSE]
    pc1 <- stats::prcomp(S, center = TRUE, scale. = TRUE)$x[, 1]
    vv[idx[which.max(abs(suppressWarnings(stats::cor(S, pc1))))]]
  }, character(1))
}

# ── Set 4: top-K correlation screen (the CURRENT production recipe) ─────────
screen_topK <- function(Xtr, ytr, vars, K = 30L) {
  r <- abs(suppressWarnings(apply(Xtr[, vars, drop = FALSE], 2,
                                  function(z) stats::cor(z, ytr))))
  r[!is.finite(r)] <- 0
  vars[order(r, decreasing = TRUE)[seq_len(min(K, length(vars)))]]
}

# ── Set 5: decorrelated-then-screened (parsimony + relevance) ──────────────
# Cluster first (outcome-free), then keep the K clusters whose representative
# correlates best with the outcome. Caps redundancy AND targets signal.
decorr_then_screen <- function(Xtr, ytr, vars, k_clust = 40L, K = 12L) {
  reps <- decorr_reps(Xtr, vars, k = k_clust)
  screen_topK(Xtr, ytr, reps, K = K)
}
