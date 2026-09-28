# =============================================================================
# explore/scripts/02_alphaearth_kernel.R   [probe AE-01]
#
# QUESTION. AlphaEarth gives 64 learned embedding dimensions per Admin-2
# (`aef_A00..aef_A63`). Scored as "a domain" the record puts it at +0.012 /
# -0.012, i.e. nothing. But the domain representation collapses a domain to its
# leading principal components - and for a LEARNED EMBEDDING, whose whole
# design is that the INNER PRODUCT encodes similarity, taking PC1 is arguably
# the one representation guaranteed to throw the information away.
#
# Does the embedding carry signal when used as a KERNEL instead of as columns?
#
# ARMS
#   aef_linear     BLUP on K = XX'/p over the 64 embedding dims
#   aef_cosine     BLUP on the cosine kernel (the embedding's native metric)
#   aef_index      the CURRENT representation: the index on the embedding's
#                  domain PCs alone - the like-for-like comparator
#   aef_plus_cs    embedding and climate+soil as separate kernels
#   cs_linear      climate+soil kernel alone
#   + the record's comparators on identical folds
#
#   Rscript explore/scripts/02_alphaearth_kernel.R
# -> explore/out/02_alphaearth_cells.csv, 02_alphaearth_loco.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
DOM <- E$domain_of
EMB  <- "Satellite embedding"
CLIM <- "Climate and weather"
SOIL <- "Soil characteristics"

cols_in <- function(X, doms) colnames(X)[which(DOM[colnames(X)] %in% doms)]

#' The index restricted to one domain's axes: the representation on the record
make_index_on <- function(doms) {
  function(tr, te, y, X, D, aux) {
    keep <- grep(paste(make.names(substr(doms, 1, 12)), collapse = "|"),
                 colnames(D), value = TRUE)
    if (length(keep) < 1) return(rep(mean(y[tr]), length(te)))
    arm_domain_index_v2(tr, te, y, X, D[, keep, drop = FALSE], aux)
  }
}

ARMS <- c(
  exp_baseline_arms(),
  list(
    aef_index  = make_index_on(EMB),
    aef_linear = make_blup_arm(function(X, D, aux) {
      cc <- cols_in(X, EMB); if (length(cc) < 5) return(list())
      list(aef = k_linear(X[, cc, drop = FALSE]))
    }),
    aef_cosine = make_blup_arm(function(X, D, aux) {
      cc <- cols_in(X, EMB); if (length(cc) < 5) return(list())
      list(aef = k_cosine(X[, cc, drop = FALSE]))
    }),
    cs_linear = make_blup_arm(function(X, D, aux) {
      cc <- cols_in(X, c(CLIM, SOIL)); if (length(cc) < 5) return(list())
      list(cs = k_linear(X[, cc, drop = FALSE]))
    }),
    aef_plus_cs = make_blup_arm(function(X, D, aux) {
      a <- cols_in(X, EMB); b <- cols_in(X, c(CLIM, SOIL))
      k <- list()
      if (length(a) >= 5) k$aef <- k_linear(X[, a, drop = FALSE])
      if (length(b) >= 5) k$cs  <- k_linear(X[, b, drop = FALSE])
      k
    })
  ))

rows <- list()
ix <- exp_cell_index(E)
for (i in seq_len(nrow(ix))) {
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tgt),
                     error = function(e) NULL)
    if (is.null(cell)) next
    rows[[paste(i, tgt, "A")]] <- exp_infill(cell, ARMS, reps = REPS)
    rows[[paste(i, tgt, "B")]] <- exp_region(cell, ARMS)
  }
  message("  ", ix$country[i], " ", ix$outcome[i])
}

loco <- list()
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tgt, outcomes = on)
    if (length(cl) < 3) next
    loco[[paste(tgt, on)]] <- exp_loco(cl, ARMS[setdiff(names(ARMS), "spatial")],
                                       domain_of = DOM)
    message("  LOCO ", tgt, " ", on)
  }
}

raw <- dplyr::bind_rows(rows); sm <- exp_summarise(raw)
exp_write(sm, "02_alphaearth_cells")
lc <- dplyr::bind_rows(loco); exp_write(lc, "02_alphaearth_loco")

cat("\n== in-country: median Spearman over cells ==\n")
a <- aggregate(spearman ~ estimand + target + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$estimand, a$target, -a$spearman), ], row.names = FALSE)

cat("\n== transport (LOCO) ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
b$positive <- aggregate(spearman ~ target + arm, data = lc,
                        FUN = function(z) sum(z > 0, na.rm = TRUE))$spearman
print(b[order(b$target, -b$spearman), ], row.names = FALSE)
