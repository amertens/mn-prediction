# =============================================================================
# scripts/protocol_v2/34_variance_components_ceiling.R   [VC-01]
#
# THE RELIABILITY CEILING FROM VARIANCE COMPONENTS, NOT SPLIT HALVES
#
# The ceiling of record (headroom_by_cell.csv, r_max_emp) is a split-half
# estimate: two random halves of each district's respondents, correlated
# across districts, Spearman-Brown corrected. CE-01 showed that where a
# district is a single survey cluster -- 57-85% of units in three countries --
# the two halves share the cluster (village, field team, day, assay batch), so
# the ceiling counts cluster effects as geography; on multi-cluster units
# 16-27% of it was cluster effect. That test could only run where a unit has
# two or more clusters, and says nothing about the single-cluster units the
# ceiling is mostly made of.
#
# A variance-components model uses every unit. Respondent i in cluster c in
# district d in region r:
#     y = mu + u_r + u_d + u_c + e
# fitted by REML (lme4). The cluster variance is identified by the units that
# have several clusters and, under the model, applies to the single-cluster
# units too. From the four components and each district's cluster count m_d
# and respondent count n_d:
#     var(observed district mean)   = s2_r + s2_d + s2_c/m_d + s2_e/n_d
#     geography a covariate can see = s2_r + s2_d
#     ceiling_vc        = sqrt((s2_r + s2_d) / var_obs)            honest
#     ceiling_within_vc = sqrt((s2_r + s2_d + s2_c/m_d) / var_obs)
#                         what a respondent split-half measures (cluster
#                         counted as geography); validated against r_max_emp
# The gap between the two is the cluster-effect share of the published
# ceiling, measured on ALL units. Prevalence uses a linear-probability model
# (the scale the split-half uses); the continuous target is the negated log
# biomarker of targets_v2. Fits are unweighted: survey weights are not
# precision weights and would mis-state the components.
#
# Run at the Admin-2 rung for all four countries, and additionally at Malawi's
# district rung (27 districts in 3 real regions), the consistent rung of R6-02.
#
#   Rscript scripts/protocol_v2/34_variance_components_ceiling.R
# -> results/tables/protocol_v2/variance_components_ceiling.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets); library(lme4)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; STORE <- "_targets_full"
cfg <- get_country_configs()
MALAWI_REGION <- c(
  Chitipa = "Northern", Karonga = "Northern", Likoma = "Northern", Mzimba = "Northern",
  `Nkhata Bay` = "Northern", Rumphi = "Northern",
  Dedza = "Central", Dowa = "Central", Kasungu = "Central", Lilongwe = "Central", Mchinji = "Central",
  Nkhotakota = "Central", Ntcheu = "Central", Ntchisi = "Central", Salima = "Central",
  Balaka = "Southern", Blantyre = "Southern", Chikwawa = "Southern", Chiradzulu = "Southern",
  Machinga = "Southern", Mangochi = "Southern", Mulanje = "Southern", Mwanza = "Southern",
  Neno = "Southern", Nsanje = "Southern", Phalombe = "Southern", Thyolo = "Southern", Zomba = "Southern")

vc_fit <- function(y, region, unit, cluster) {
  d <- data.frame(y = y, region = factor(region), unit = factor(paste(region, unit, sep = "|")),
                  cluster = factor(paste(region, unit, cluster, sep = "|")))
  d <- d[is.finite(d$y) & !is.na(d$region) & !is.na(d$unit), ]
  d <- droplevels(d)
  if (nrow(d) < 40 || nlevels(d$unit) < 6) return(NULL)
  fit <- tryCatch(suppressMessages(suppressWarnings(lme4::lmer(
    y ~ 1 + (1 | region) + (1 | unit) + (1 | cluster), data = d, REML = TRUE,
    control = lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)))), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  vc <- as.data.frame(lme4::VarCorr(fit)); v <- stats::setNames(vc$vcov, vc$grp)
  m_u <- tapply(as.character(d$cluster), d$unit, function(z) length(unique(z))); n_u <- tapply(d$y, d$unit, length)
  s2r <- unname(v["region"]); s2d <- unname(v["unit"]); s2c <- unname(v["cluster"]); s2e <- unname(v["Residual"])
  tot <- s2r + s2d + s2c + s2e
  geo <- s2r + s2d; nc <- mean(s2c / m_u); ne <- mean(s2e / n_u); vo <- geo + nc + ne
  c_vc <- sqrt(geo / vo); c_within <- sqrt((geo + nc) / vo)
  data.frame(n_resp = nrow(d), n_units = length(m_u), n_regions = nlevels(d$region), n_clusters = nlevels(d$cluster),
             mean_clusters_per_unit = round(mean(m_u), 2), share_single_cluster = round(mean(m_u == 1), 3),
             median_n_per_unit = stats::median(n_u),
             s2_region = s2r, s2_district = s2d, s2_cluster = s2c, s2_resid = s2e, singular = lme4::isSingular(fit),
             icc_region = s2r / tot, icc_district = s2d / tot, icc_cluster = s2c / tot,
             ceiling_vc = c_vc, ceiling_within_vc = c_within,
             cluster_share_of_ceiling = if (c_within > 0) 1 - c_vc / c_within else NA_real_)
}

rows <- list()
for (cn in names(cfg)) { cc <- cfg[[cn]]; lc <- tolower(cn); clu <- cc$cluster_id
  for (on in names(cc$outcomes)) { oc <- cc$outcomes[[on]]
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL)
    if (is.null(od)) { cat(sprintf("  %-12s %-12s skip: not in store\n", cn, on)); next }
    d <- od$data; if (!all(c(clu, "Admin1", "Admin2") %in% names(d))) { cat(sprintf("  %-12s %-12s skip: missing columns\n", cn, on)); next }
    yb <- tryCatch(resolve_uniform_outcome(d, cc, oc), error = function(e) NULL)
    ybin <- if (!is.null(yb)) as.numeric(yb) else .v2_num(d[[oc$binary]])
    ycont <- rep(NA_real_, nrow(d))
    if (!is.null(oc$continuous) && oc$continuous %in% names(d)) {
      v <- .v2_num(d[[oc$continuous]]); t <- if (identical(oc$cutoff_scale, "log")) v else { v[!is.finite(v) | v <= 0] <- NA; log(v) }; ycont <- -t }
    a1 <- as.character(d$Admin1); a2 <- as.character(d$Admin2); cl <- as.character(d[[clu]])
    rungs <- list(admin2 = list(region = a1, unit = a2))
    if (cn == "Malawi") rungs$district <- list(region = unname(MALAWI_REGION[a1]), unit = a1)
    for (rg in names(rungs)) for (tg in c("prev", "level")) {
      y <- if (tg == "prev") ybin else ycont
      if (sum(is.finite(y)) < 40) { cat(sprintf("  %-12s %-12s %-8s %-5s skip: too few finite outcomes\n", cn, on, rg, tg)); next }
      r <- tryCatch(vc_fit(y, rungs[[rg]]$region, rungs[[rg]]$unit, cl), error = function(e) { cat("  fit error", cn, on, rg, tg, conditionMessage(e), "\n"); NULL })
      if (is.null(r)) { cat(sprintf("  %-12s %-12s %-8s %-5s skip: fit failed\n", cn, on, rg, tg)); next }
      rows[[length(rows) + 1L]] <- cbind(data.frame(country = cn, outcome = on, rung = rg, target = tg, stringsAsFactors = FALSE), r)
      cat(sprintf("  %-12s %-12s %-8s %-5s units %3d (single-cluster %3.0f%%) | ceiling honest %.2f | cluster-as-geography %.2f | cluster share %3.0f%%%s\n",
                  cn, on, rg, tg, r$n_units, 100 * r$share_single_cluster, r$ceiling_vc, r$ceiling_within_vc, 100 * r$cluster_share_of_ceiling, if (r$singular) " [singular]" else ""))
    }
  } }
R <- bind_rows(rows); if (!nrow(R)) stop("no cells fitted")

# ── joins: published split-half ceiling, CE-01 cluster split, achieved skill ──
HB <- tryCatch(read.csv(file.path(OUTDIR, "headroom_by_cell.csv")), error = function(e) NULL)
if (!is.null(HB)) R <- left_join(R, HB[, c("country", "outcome", "r_max_emp")], by = c("country", "outcome"))
CE <- tryCatch(read.csv(file.path(OUTDIR, "ceiling_cluster_split.csv")), error = function(e) NULL)
if (!is.null(CE)) R <- left_join(R, CE[, c("country", "outcome", "ceiling_within", "ceiling_cluster")] |> rename(ce01_within = ceiling_within, ce01_cluster = ceiling_cluster), by = c("country", "outcome"))
BM <- tryCatch(read.csv(file.path(OUTDIR, "benchmarks_v2_cells.csv")), error = function(e) NULL)
if (!is.null(BM)) { ach <- BM |> filter(arm == "domain_index", estimand == "infill") |> group_by(country, outcome, target) |> summarise(achieved_spearman = mean(spearman, na.rm = TRUE), .groups = "drop")
  R <- left_join(R, ach, by = c("country", "outcome", "target")) }
if (!"achieved_spearman" %in% names(R)) R$achieved_spearman <- NA_real_
if (!"r_max_emp" %in% names(R)) R$r_max_emp <- NA_real_
R$headroom_ratio <- ifelse(is.finite(R$achieved_spearman) & R$ceiling_vc > 0, R$achieved_spearman / R$ceiling_vc, NA_real_)
R <- R |> mutate(across(c(s2_region, s2_district, s2_cluster, s2_resid, icc_region, icc_district, icc_cluster, ceiling_vc, ceiling_within_vc, cluster_share_of_ceiling, headroom_ratio), ~ round(.x, 4)))
write.csv(R, file.path(OUTDIR, "variance_components_ceiling.csv"), row.names = FALSE)

cat("\n===== VC-01: variance-components ceiling (all units) =====\n")
A2 <- R[R$rung == "admin2", ]
cat("\n-- validation: cluster-as-geography VC ceiling vs published split-half r_max_emp (prevalence, Admin-2) --\n")
v <- A2[A2$target == "prev" & is.finite(A2$r_max_emp), ]
if (nrow(v) > 2) cat(sprintf("  cells %d | mean VC-within %.3f vs r_max_emp %.3f | Pearson across cells %.2f | mean abs diff %.3f\n",
            nrow(v), mean(v$ceiling_within_vc), mean(v$r_max_emp), suppressWarnings(stats::cor(v$ceiling_within_vc, v$r_max_emp)), mean(abs(v$ceiling_within_vc - v$r_max_emp))))
cat("\n-- honest ceiling by country and target (Admin-2 rung) --\n")
print(as.data.frame(A2 |> group_by(country, target) |> summarise(cells = dplyr::n(), single_cluster = round(mean(share_single_cluster), 2),
  ceiling_within = round(mean(ceiling_within_vc), 3), ceiling_honest = round(mean(ceiling_vc), 3), cluster_share = round(mean(cluster_share_of_ceiling, na.rm = TRUE), 3),
  icc_cluster = round(mean(icc_cluster), 3), icc_district = round(mean(icc_district), 3), icc_region = round(mean(icc_region), 3), .groups = "drop")), row.names = FALSE)
cat("\n-- per cell, prevalence, Admin-2 --\n")
print(as.data.frame(A2[A2$target == "prev", c("country", "outcome", "n_units", "share_single_cluster", "r_max_emp", "ceiling_within_vc", "ceiling_vc", "cluster_share_of_ceiling", "achieved_spearman", "headroom_ratio", "singular")]), row.names = FALSE)
cat("\n-- per cell, level, Admin-2 --\n")
print(as.data.frame(A2[A2$target == "level", c("country", "outcome", "n_units", "ceiling_within_vc", "ceiling_vc", "cluster_share_of_ceiling", "achieved_spearman", "headroom_ratio", "singular")]), row.names = FALSE)
MW <- R[R$rung == "district", ]
if (nrow(MW)) { cat("\n-- Malawi at its district rung (27 districts, 3 regions) --\n")
  print(as.data.frame(MW[, c("outcome", "target", "n_units", "mean_clusters_per_unit", "ceiling_within_vc", "ceiling_vc", "cluster_share_of_ceiling", "achieved_spearman", "singular")]), row.names = FALSE) }
pv <- A2[A2$target == "prev", ]
cat(sprintf("\noverall (Admin-2, prevalence): honest ceiling mean %.3f | cluster-as-geography %.3f | cluster share %.0f%% | cells with honest ceiling < 0.30: %d of %d\n",
            mean(pv$ceiling_vc), mean(pv$ceiling_within_vc), 100 * mean(pv$cluster_share_of_ceiling, na.rm = TRUE), sum(pv$ceiling_vc < 0.30), nrow(pv)))
cat("\nDONE\n")
