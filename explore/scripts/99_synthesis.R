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
paired <- function(df, arm, ref, key = c("country", "outcome"),
                   filt = function(d) d) {
  d <- filt(df)
  d <- d[d$arm %in% c(arm, ref), c(key, "arm", "spearman")]
  if (!nrow(d)) return(NULL)
  w <- stats::reshape(d, idvar = key, timevar = "arm", direction = "wide")
  names(w) <- sub("^spearman[.]", "", names(w))
  if (!all(c(arm, ref) %in% names(w))) return(NULL)
  g <- w[[arm]] - w[[ref]]
  ok <- is.finite(g)
  if (!any(ok)) return(NULL)
  data.frame(arm = arm, reference = ref,
             arm_mean = mean(w[[arm]][ok]), ref_mean = mean(w[[ref]][ok]),
             mean_gain = mean(g[ok]), wins = sum(g[ok] > 0), cells = sum(ok),
             sign_p = stats::binom.test(sum(g[ok] > 0), sum(ok), 0.5)$p.value,
             stringsAsFactors = FALSE)
}

rows <- list()
add <- function(probe, estimand, target, x) {
  if (is.null(x)) return(invisible())
  rows[[length(rows) + 1L]] <<- cbind(
    data.frame(probe = probe, estimand = estimand, target = target,
               stringsAsFactors = FALSE), x)
}

infill_lvl <- function(d) d[d$estimand == "infill" & d$target == "level", ]
infill_prv <- function(d) d[d$estimand == "infill" & d$target == "prev", ]
loco_lvl   <- function(d) d[d$target == "level", ]
loco_prv   <- function(d) d[d$target == "prev", ]

# ── KB-01 kernel BLUP ───────────────────────────────────────────────────────
k <- rd("01_kernel_blup_cells.csv"); kl <- rd("01_kernel_blup_loco.csv")
if (!is.null(k)) for (a in c("blup_all", "blup_cs", "blup_cs_spatial", "blup_5k", "blup_spatial")) {
  add("KB-01", "infill", "level", paired(k, a, "domain_index", filt = infill_lvl))
  add("KB-01", "infill", "level", paired(k, a, "spatial",      filt = infill_lvl))
  add("KB-01", "infill", "prev",  paired(k, a, "domain_index", filt = infill_prv))
}
if (!is.null(kl)) for (a in c("blup_all", "blup_cs", "blup_cs_spatial", "blup_5k", "blup_spatial")) {
  add("KB-01", "transport", "level", paired(kl, a, "domain_index", filt = loco_lvl))
  add("KB-01", "transport", "prev",  paired(kl, a, "domain_index", filt = loco_prv))
}
# THE NESTED TEST. The project's standing conclusion is that within a country
# covariates add nothing on top of a spatial smoother. Comparing a covariate
# BLUP to the GAM smoother confounds "covariates help" with "the kernel is a
# better smoother". blup_cs_spatial vs blup_spatial is the clean contrast:
# identical machinery, identical folds, covariates added as a second kernel.
if (!is.null(k)) {
  add("KB-01 nested", "infill", "level", paired(k, "blup_cs_spatial", "blup_spatial", filt = infill_lvl))
  add("KB-01 nested", "infill", "prev",  paired(k, "blup_cs_spatial", "blup_spatial", filt = infill_prv))
  add("KB-01 nested", "infill", "level", paired(k, "blup_5k", "blup_spatial", filt = infill_lvl))
}

# ── AE-01 AlphaEarth ────────────────────────────────────────────────────────
a2 <- rd("02_alphaearth_cells.csv"); a2l <- rd("02_alphaearth_loco.csv")
if (!is.null(a2)) {
  add("AE-01", "infill", "level", paired(a2, "aef_linear", "aef_index", filt = infill_lvl))
  add("AE-01", "infill", "level", paired(a2, "cs_linear",  "domain_index", filt = infill_lvl))
  add("AE-01", "infill", "level", paired(a2, "cs_linear",  "spatial", filt = infill_lvl))
}
if (!is.null(a2l)) {
  add("AE-01", "transport", "level", paired(a2l, "aef_linear",  "aef_index", filt = loco_lvl))
  add("AE-01", "transport", "level", paired(a2l, "aef_plus_cs", "cs_linear", filt = loco_lvl))
  add("AE-01", "transport", "level", paired(a2l, "cs_linear",   "domain_index", filt = loco_lvl))
}

# ── RV-01 RUV ───────────────────────────────────────────────────────────────
r3 <- rd("03_ruv_loco.csv")
if (!is.null(r3)) for (a in c("ruv_bc_k1", "ruv_bc_k2", "ruv_pc_k1", "ruv_pc_k3", "index_cs")) {
  add("RV-01", "transport", "level", paired(r3, a, "index", filt = loco_lvl))
}

# ── EB-01 moderated meta-index ──────────────────────────────────────────────
m4 <- rd("04_moderated_cells.csv"); m4l <- rd("04_moderated_loco.csv")
if (!is.null(m4)) {
  add("EB-01", "infill", "level", paired(m4, "meta_index", "domain_index", filt = infill_lvl))
  add("EB-01", "infill", "prev",  paired(m4, "meta_index", "domain_index", filt = infill_prv))
}
if (!is.null(m4l)) {
  add("EB-01", "transport", "level", paired(m4l, "meta_index", "domain_index", filt = loco_lvl))
  add("EB-01", "transport", "prev",  paired(m4l, "meta_index", "domain_index", filt = loco_prv))
}

# ── MT-01 multi-trait ───────────────────────────────────────────────────────
m5 <- rd("05_multitrait_cells.csv")
if (!is.null(m5)) {
  add("MT-01", "infill", "level", paired(m5, "mt_blup", "st_blup",
                                         filt = function(d) d[d$target == "level", ]))
  add("MT-01", "infill", "prev",  paired(m5, "mt_blup", "st_blup",
                                         filt = function(d) d[d$target == "prev", ]))
}

# ── MX-01 mechanistic ───────────────────────────────────────────────────────
m7 <- rd("07_mechanistic_cells.csv")
if (!is.null(m7)) {
  add("MX-01", "infill", "level", paired(m7, "mech_matched", "mech_mismatched", filt = infill_lvl))
  add("MX-01", "infill", "level", paired(m7, "mech_matched", "domain_index", filt = infill_lvl))
  add("MX-01", "infill", "level", paired(m7, "index_plus_mech", "domain_index", filt = infill_lvl))
}

# ── NP-01 n<<p estimators ───────────────────────────────────────────────────
n8 <- rd("08_np_cells.csv"); n8l <- rd("08_np_loco.csv")
if (!is.null(n8)) for (a in c("pls2", "pcr2", "spca", "mcp", "stabsel", "ridge_cv", "index_std")) {
  add("NP-01", "infill", "level", paired(n8, a, "domain_index", filt = infill_lvl))
}
if (!is.null(n8l)) for (a in c("pls2", "pcr2", "spca", "mcp", "stabsel", "ridge_cv", "index_std")) {
  add("NP-01", "transport", "level", paired(n8l, a, "domain_index", filt = loco_lvl))
}

# ── TM-01 temporal ──────────────────────────────────────────────────────────
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
L$mean_gain <- round(L$mean_gain, 4)
L$arm_mean <- round(L$arm_mean, 3); L$ref_mean <- round(L$ref_mean, 3)
L$sign_p <- round(L$sign_p, 4)
L$verdict <- ifelse(L$sign_p < 0.05 & L$mean_gain > 0, "candidate",
             ifelse(L$sign_p < 0.05 & L$mean_gain < 0, "harmful", "no effect"))
L <- L[order(L$probe, L$estimand, L$target, -L$mean_gain), ]
exp_write(L, "leads")

cat("\n== every arm against its comparator, paired ==\n")
print(L[, c("probe", "estimand", "target", "arm", "reference", "arm_mean",
            "ref_mean", "mean_gain", "wins", "cells", "sign_p", "verdict")],
      row.names = FALSE)

cat("\n== candidates (paired-positive at p < 0.05) ==\n")
cand <- L[L$verdict == "candidate", ]
if (nrow(cand)) print(cand[order(-cand$mean_gain),
  c("probe", "estimand", "target", "arm", "reference", "mean_gain", "wins", "cells", "sign_p")],
  row.names = FALSE) else cat("  none\n")
