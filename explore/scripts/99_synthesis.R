# =============================================================================
# explore/scripts/99_synthesis.R
#
# One table of every probe: the arm, its honest number on each estimand, the
# comparator it must beat, and whether the comparison is paired-positive.
#
# READ THE PAIRED COLUMN, NOT THE MEDIAN. AE-01 found the median across cells
# and the paired per-cell comparison disagreeing in DIRECTION, because the
# cells differ enormously in difficulty. `wins` is the number of cells where
# the arm beat the reference on the same folds; that is the number to trust.
#
#   Rscript explore/scripts/99_synthesis.R
# -> explore/out/leads.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

rd <- function(f) {
  p <- file.path(EXP_OUT, f)
  if (!file.exists(p)) { message("missing: ", f); return(NULL) }
  read.csv(p, stringsAsFactors = FALSE)
}

#' Paired comparison of one arm against a reference, on the cells both have
#'
#' TWO UNITS, AND THE BLOCK ONE IS THE VERDICT. Cells of the same country share
#' districts and predictors, so 18-22 cells are not 18-22 independent
#' comparisons and the cell-level sign test is anti-conservative about
#' correlation in one direction and badly under-powered in the other. Both
#' MT-01 and KB-01 read "no effect" on cells and are unambiguous on blocks
#' (6/6 and 8/8). The country is the honest unit: for in-country estimands it
#' is the country the cells belong to, for transport the held-out country.
paired <- function(df, arm, ref, key = c("country", "outcome", "target"),
                   filt = function(d) d) {
  d <- filt(df)
  key <- intersect(key, names(df))
  d <- d[d$arm %in% c(arm, ref), c(key, "arm", "spearman")]
  if (!nrow(d)) return(NULL)
  w <- stats::reshape(d, idvar = key, timevar = "arm", direction = "wide")
  names(w) <- sub("^spearman[.]", "", names(w))
  if (!all(c(arm, ref) %in% names(w))) return(NULL)
  g <- w[[arm]] - w[[ref]]
  ok <- is.finite(g)
  if (!any(ok)) return(NULL)
  w <- w[ok, ]; g <- g[ok]

  # Blocks are country x target. Country alone gives at most 4 blocks, and
  # binom.test(4, 4) is p = 0.125 - unanimity could never clear 0.05, so the
  # verdict column would read "none" however strong the effect. Pooling the two
  # targets gives 6-8 blocks, which is the unit the MT-01 and KB-01 write-ups
  # used and the smallest one that can express significance here.
  bkey <- if ("target" %in% names(w)) paste(w$country, w$target) else w$country
  b <- stats::aggregate(list(gain = g), by = list(block = bkey), FUN = mean)
  bp <- stats::binom.test(sum(b$gain > 0), nrow(b), 0.5)$p.value

  data.frame(arm = arm, reference = ref,
             arm_mean = mean(w[[arm]]), ref_mean = mean(w[[ref]]),
             mean_gain = mean(g), wins = sum(g > 0), cells = length(g),
             cell_p = stats::binom.test(sum(g > 0), length(g), 0.5)$p.value,
             blocks_pos = sum(b$gain > 0), blocks = nrow(b),
             block_gain = mean(b$gain), block_p = bp,
             stringsAsFactors = FALSE)
}

rows <- list()
add <- function(probe, estimand, target, x) {
  if (is.null(x)) return(invisible())
  rows[[length(rows) + 1L]] <<- cbind(
    data.frame(probe = probe, estimand = estimand, target = target,
               stringsAsFactors = FALSE), x)
}

# Filters select the ESTIMAND only. Both targets go into one comparison, and
# the blocks are country x target, so a probe gets 6-8 blocks instead of 4.
IF <- function(d) d[d$estimand == "infill", ]
RG <- function(d) d[d$estimand == "region", ]
LO <- function(d) d                      # the loco files hold only estimand C

# ── KB-01 kernel BLUP ───────────────────────────────────────────────────────
k <- rd("01_kernel_blup_cells.csv"); kl <- rd("01_kernel_blup_loco.csv")
KARMS <- c("blup_all", "blup_cs", "blup_cs_spatial", "blup_5k", "blup_spatial")
if (!is.null(k)) for (a in KARMS) {
  add("KB-01", "infill", "both", paired(k, a, "domain_index", filt = IF))
  add("KB-01", "infill", "both", paired(k, a, "spatial",      filt = IF))
}
if (!is.null(kl)) for (a in KARMS)
  add("KB-01", "transport", "both", paired(kl, a, "domain_index", filt = LO))

# THE NESTED TEST. The project's standing conclusion is that within a country
# covariates add nothing on top of a spatial smoother. Comparing a covariate
# BLUP to the GAM smoother confounds "covariates help" with "the kernel is a
# better smoother". blup_cs_spatial vs blup_spatial is the clean contrast:
# identical machinery, identical folds, covariates added as a second kernel.
if (!is.null(k)) for (a in c("blup_cs_spatial", "blup_5k", "blup_cs"))
  add("KB-01 nested", "infill", "both", paired(k, a, "blup_spatial", filt = IF))

# ── AE-01 AlphaEarth ────────────────────────────────────────────────────────
a2 <- rd("02_alphaearth_cells.csv"); a2l <- rd("02_alphaearth_loco.csv")
if (!is.null(a2)) {
  add("AE-01", "infill", "both", paired(a2, "aef_linear", "aef_index", filt = IF))
  add("AE-01", "infill", "both", paired(a2, "aef_cosine", "aef_index", filt = IF))
  add("AE-01", "infill", "both", paired(a2, "cs_linear",  "domain_index", filt = IF))
}
if (!is.null(a2l)) {
  add("AE-01", "transport", "both", paired(a2l, "aef_linear",  "aef_index", filt = LO))
  add("AE-01", "transport", "both", paired(a2l, "aef_plus_cs", "cs_linear", filt = LO))
  add("AE-01", "transport", "both", paired(a2l, "cs_linear",   "domain_index", filt = LO))
}

# ── RV-01 RUV ───────────────────────────────────────────────────────────────
r3 <- rd("03_ruv_loco.csv")
if (!is.null(r3)) for (a in c("index_cs", "ruv_bc_k1", "ruv_bc_k2", "ruv_pc_k1", "ruv_pc_k3"))
  add("RV-01", "transport", "both", paired(r3, a, "index", filt = LO))

# ── EB-01 moderated meta-index ──────────────────────────────────────────────
m4 <- rd("04_moderated_cells.csv"); m4l <- rd("04_moderated_loco.csv")
if (!is.null(m4))  add("EB-01", "infill", "both",
                       paired(m4, "meta_index", "domain_index", filt = IF))
if (!is.null(m4l)) add("EB-01", "transport", "both",
                       paired(m4l, "meta_index", "domain_index", filt = LO))

# ── MT-01 multi-trait ───────────────────────────────────────────────────────
m5 <- rd("05_multitrait_cells.csv")
if (!is.null(m5)) add("MT-01", "infill", "both",
                      paired(m5, "mt_blup", "st_blup", filt = IF))

# ── MX-01 mechanistic ───────────────────────────────────────────────────────
m7 <- rd("07_mechanistic_cells.csv")
if (!is.null(m7)) {
  add("MX-01", "infill", "both", paired(m7, "mech_matched", "mech_mismatched", filt = IF))
  add("MX-01", "infill", "both", paired(m7, "mech_matched", "domain_index", filt = IF))
  add("MX-01", "infill", "both", paired(m7, "index_plus_mech", "domain_index", filt = IF))
}

# ── NP-01 n<<p estimators ───────────────────────────────────────────────────
NPARMS <- c("pls2", "pcr2", "spca", "mcp", "stabsel", "ridge_cv", "index_std")
n8 <- rd("08_np_cells.csv"); n8l <- rd("08_np_loco.csv")
if (!is.null(n8))  for (a in NPARMS)
  add("NP-01", "infill", "both", paired(n8, a, "domain_index", filt = IF))
if (!is.null(n8l)) for (a in NPARMS)
  add("NP-01", "transport", "both", paired(n8l, a, "domain_index", filt = LO))

# ── TM-01 temporal (cluster track; blocks are country only, one target) ─────
t9 <- rd("09_temporal_cluster.csv")
if (!is.null(t9)) {
  f <- function(d) d[d$unit == "district", ]
  add("TM-01", "cluster->district", "level", paired(t9, "clim_fw", "clim", filt = f))
  add("TM-01", "cluster->district", "level", paired(t9, "clim_fw_phase", "clim_fw", filt = f))
}
t11 <- rd("11_lag_stack.csv")
if (!is.null(t11)) {
  f <- function(d) d[d$unit == "district", ]
  for (a in c("clim_lag_raw", "clim_lag_spline", "clim_lag_window", "clim_prev_seas"))
    add("TM-01b", "cluster->district", "level", paired(t11, a, "clim", filt = f))
}

L <- dplyr::bind_rows(rows)
for (v in c("mean_gain", "block_gain")) L[[v]] <- round(L[[v]], 4)
L$arm_mean <- round(L$arm_mean, 3); L$ref_mean <- round(L$ref_mean, 3)
L$cell_p <- round(L$cell_p, 4); L$block_p <- round(L$block_p, 4)

# Verdict on the BLOCK test, with the whole-block unanimity requirement that
# the small number of countries makes necessary: with 4-6 blocks the sign test
# cannot reach 0.05 unless every block agrees.
L$verdict <- ifelse(L$block_p < 0.05 & L$block_gain > 0, "candidate",
             ifelse(L$block_p < 0.05 & L$block_gain < 0, "harmful", "no effect"))
L <- L[order(L$probe, L$estimand, L$target, -L$block_gain), ]
exp_write(L, "leads")

cat("\n== every arm against its comparator ==\n")
cat("   cells = country x outcome (correlated within country) | blocks = countries\n\n")
print(L[, c("probe", "estimand", "target", "arm", "reference", "arm_mean",
            "ref_mean", "mean_gain", "wins", "cells", "cell_p",
            "blocks_pos", "blocks", "block_p", "verdict")],
      row.names = FALSE)

cat("\n== CANDIDATES (every country block agrees, p < 0.05) ==\n")
cand <- L[L$verdict == "candidate", ]
if (nrow(cand)) print(cand[order(-cand$block_gain),
  c("probe", "estimand", "target", "arm", "reference", "block_gain",
    "blocks_pos", "blocks", "block_p")], row.names = FALSE) else cat("  none\n")

cat("\n== HARMFUL (every country block agrees the arm is worse) ==\n")
harm <- L[L$verdict == "harmful", ]
if (nrow(harm)) print(harm[order(harm$block_gain),
  c("probe", "estimand", "target", "arm", "reference", "block_gain",
    "blocks_pos", "blocks", "block_p")], row.names = FALSE) else cat("  none\n")
