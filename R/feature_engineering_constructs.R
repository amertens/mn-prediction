# =============================================================================
# R/feature_engineering_constructs.R   (EXPERIMENTAL — not wired into _targets)
#
# Construct-level feature engineering for proxy predictors. Provided as a
# reviewed draft; deliberately NOT called by the primary fit path.
#
# WHY EXPERIMENTAL — read before using (see sandbox_fe/FINDINGS.md, §metadata):
#   * The only per-column metadata we have is the NAME STRING. There is no
#     curated variable dictionary mapping the ~600 gee_ / ~120 ihme_ columns to
#     constructs, so construct labels here are heuristically parsed and an
#     "other" bucket remains.
#   * Source-prefix "domains" are incoherent for PCA (GEE top-3 PCs hold only
#     ~55% of variance). Construct grouping concentrates variance for some
#     (cropland/elevation ~rank-1) but climate constructs (precip/temp) are
#     genuinely multi-dimensional.
#   * Unsupervised PCA loadings are NOT transportable across countries for the
#     climate/malaria constructs (cross-country |cor(PC1 loadings)| ~0.22–0.46),
#     which is the mechanism behind observed negative transfer. Only 4 countries
#     exist, too few to learn which components transfer.
#
# RECOMMENDED USE:
#   1. `construct_of()` + `build_construct_summaries()` produce INTERPRETABLE,
#      fixed-definition seasonal features (mean / amplitude / phase). These are
#      transportable BY CONSTRUCTION (formula does not depend on the data) and
#      are the safe way to compress the 12-month climatologies.
#   2. `reduce_pca_pooled()` is offered only for the coherent + transportable
#      constructs, and MUST be fit on POOLED multi-country training data so the
#      loadings are shared. Treat outputs as experimental.
# =============================================================================

MONTHS_ABBR <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

#' Map proxy column names to measurement constructs (heuristic, name-based).
#' @return character vector parallel to `cols`.
construct_of <- function(cols) {
  out <- rep("other", length(cols)); g <- function(p) grepl(p, cols, ignore.case = TRUE)
  out[g("trmm|chirps|precip|rainfall")]            <- "precip"
  out[g("ndvi|evi|vegetation|gpp|npp|\\blai\\b")]  <- "vegetation"
  out[g("temp|_lst|tavg|tair|temperature")]        <- "temperature"
  out[g("elevation|srtm|slope|\\bdem\\b")]         <- "elevation"
  out[g("accessib|traveltime|friction")]           <- "accessibility"
  out[g("ccnl|nightlight|_ntl|viirs")]             <- "nightlight"
  out[g("worldcereal|cropland|cereal|maize|_crop")] <- "cropland"
  out[g("aerosol|atmosphere|\\baod\\b|no2")]       <- "atmosphere"
  out[g("^MAP_")]                                  <- "malaria"
  out[g("soil")]                                   <- "soil"
  out[g("^fsec")]                                  <- "foodsec"
  out[g("^ihme") & out == "other"]                 <- "ihme_health"
  out
}

# Constructs whose PCA loadings were transportable enough across countries to
# be candidates for unsupervised reduction (vs the non-transportable rest).
TRANSPORTABLE_CONSTRUCTS <- c("cropland", "elevation", "vegetation", "ihme_health")

#' Build interpretable, transportable seasonal summaries from monthly columns
#' named like  gee_<var>_<Mon>_<suffix>  (e.g. gee_ndvi_Jan_10km).
#' For each variable with >=8 months present: mean, amplitude (max-min), CV,
#' and the first annual Fourier harmonic (amplitude + sin/cos of phase).
#' Returns a data.frame keyed row-for-row with `df` (new fe_ columns only).
build_construct_summaries <- function(df) {
  cols <- colnames(df)
  pat  <- paste0("_(", paste(MONTHS_ABBR, collapse = "|"), ")_")
  monthly <- cols[grepl(pat, cols)]
  if (!length(monthly)) return(data.frame(row.names = seq_len(nrow(df))))
  mo   <- regmatches(monthly, regexpr(paste0("(", paste(MONTHS_ABBR, collapse="|"), ")"), monthly))
  stem <- sub(pat, "_MO_", monthly)
  groups <- split(data.frame(col = monthly, m = match(mo, MONTHS_ABBR),
                             stringsAsFactors = FALSE), stem)
  ang <- function(m) 2 * pi * m / 12
  out <- list()
  for (st in names(groups)) {
    g <- groups[[st]]; if (nrow(g) < 8) next
    g <- g[order(g$m), ]
    M <- as.matrix(sapply(g$col, function(c) suppressWarnings(as.numeric(df[[c]]))))
    if (is.null(dim(M))) next
    mu  <- rowMeans(M, na.rm = TRUE)
    amp <- apply(M, 1, function(r) { r <- r[is.finite(r)]; if (length(r) < 2) NA else max(r) - min(r) })
    cv  <- apply(M, 1, function(r) { r <- r[is.finite(r)]; if (length(r) < 2 || mean(r) == 0) NA else sd(r)/abs(mean(r)) })
    a <- as.numeric(M %*% cos(ang(g$m))) * 2 / nrow(g)
    b <- as.numeric(M %*% sin(ang(g$m))) * 2 / nrow(g)
    base <- gsub("[^A-Za-z0-9]+", "_", sub("_MO_.*$", "", sub("^gee_", "", st)))
    out[[paste0("fe_", base, "_mean")]]  <- mu
    out[[paste0("fe_", base, "_amp")]]   <- amp
    out[[paste0("fe_", base, "_cv")]]    <- cv
    out[[paste0("fe_", base, "_h1")]]    <- sqrt(a^2 + b^2)
    out[[paste0("fe_", base, "_phsin")]] <- sin(atan2(b, a))
    out[[paste0("fe_", base, "_phcos")]] <- cos(atan2(b, a))
  }
  as.data.frame(out)
}

#' EXPERIMENTAL: pooled, construct-restricted PCA reducer.
#' Fit ONLY on pooled multi-country training data; apply to new data with the
#' returned `project()` closure so loadings are shared (transportable).
#'
#' @param train_df  pooled training predictors (numeric, imputed)
#' @param k         PCs per construct
#' @param constructs which constructs to reduce (default: transportable set)
#' @param min_keep  floor: if a construct has < min_keep cols, keep raw columns
#' @return list(features = matrix, project = function(new_df) -> matrix)
reduce_pca_pooled <- function(train_df, k = 3L,
                              constructs = TRANSPORTABLE_CONSTRUCTS,
                              min_keep = 8L) {
  stopifnot(is.data.frame(train_df) || is.matrix(train_df))
  Xtr <- as.matrix(train_df); cn <- colnames(Xtr)
  con <- construct_of(cn)
  blocks <- list()  # per-construct fitted transforms
  for (cc in constructs) {
    idx <- which(con == cc); if (length(idx) < 2) next
    Xb <- Xtr[, idx, drop = FALSE]
    mu <- colMeans(Xb, na.rm = TRUE); sdv <- apply(Xb, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv < 1e-8] <- 1
    Z <- sweep(sweep(Xb, 2, mu, "-"), 2, sdv, "/")
    if (ncol(Z) < min_keep) {            # too few cols -> keep raw (z-scored)
      blocks[[cc]] <- list(cols = cn[idx], mu = mu, sd = sdv, rot = NULL)
    } else {
      pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
      kk <- min(k, ncol(pc$rotation))
      blocks[[cc]] <- list(cols = cn[idx], mu = mu, sd = sdv,
                           rot = pc$rotation[, seq_len(kk), drop = FALSE])
    }
  }
  apply_blocks <- function(new_df) {
    Xn <- as.matrix(new_df); out <- list()
    for (cc in names(blocks)) {
      bl <- blocks[[cc]]; have <- intersect(bl$cols, colnames(Xn))
      if (length(have) < length(bl$cols)) next  # require identical schema
      Xb <- Xn[, bl$cols, drop = FALSE]
      Z  <- sweep(sweep(Xb, 2, bl$mu, "-"), 2, bl$sd, "/")
      if (is.null(bl$rot)) { colnames(Z) <- paste0(cc, "__", bl$cols); out[[cc]] <- Z }
      else { P <- Z %*% bl$rot; colnames(P) <- paste0(cc, "_PC", seq_len(ncol(P))); out[[cc]] <- P }
    }
    if (!length(out)) return(matrix(numeric(0), nrow = nrow(Xn)))
    do.call(cbind, out)
  }
  list(features = apply_blocks(train_df), project = apply_blocks, blocks = blocks)
}
