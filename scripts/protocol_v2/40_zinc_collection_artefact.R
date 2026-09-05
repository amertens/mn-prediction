# =============================================================================
# scripts/protocol_v2/40_zinc_collection_artefact.R   [ZN-02]
#
# IS MALAWI'S ZINC GEOGRAPHY A COLLECTION ARTEFACT?
#
# ZN-01 found zinc the one outcome with a high split-half ceiling and no
# covariate signal. VC-01 then showed that for child zinc the honest ceiling
# (cluster effects removed) is 0.27 on prevalence and 0.00 on the biomarker
# level, against 0.80 with cluster counted as geography: two thirds to all of
# what looked like zinc geography is something clusters share. Serum zinc is
# the biomarker most sensitive to collection conditions -- it falls after a
# meal and through the day, and drifts with time to centrifugation -- so the
# obvious candidate is when the team drew blood. The Malawi file carries
# time_blood_draw (1 = morning, 2 = afternoon) and an interview date for two
# thirds of respondents. This asks:
#   1  individual level: does serum zinc depend on time of draw and calendar
#      date, net of age, sex, district and cluster?
#   2  district level: does the district's share of afternoon draws track its
#      zinc-deficiency prevalence?
#   3  variance components: how much between-district and between-cluster
#      variance in zinc survives adjustment for time of draw and date?
#
#   Rscript scripts/protocol_v2/40_zinc_collection_artefact.R
# -> results/tables/protocol_v2/zinc_collection_artefact.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(lme4)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"
m <- readRDS("data/IPD/Malawi/clean_malawi_mn_data.RDS")
cat("rows", nrow(m), "| zn_gdl finite", sum(is.finite(m$zn_gdl)), "| zinc_def table:"); print(table(m$zinc_def, useNA = "ifany"))
cat("sex values:", paste(names(table(m$sex)), table(m$sex), collapse = " "), "| age_year range", range(m$age_year, na.rm = TRUE), "| age_month finite", sum(is.finite(m$age_month)), "\n")
cat("time_blood_draw x group:\n"); print(table(m$time_blood_draw, ifelse(is.finite(m$psc_agecat), "psc", ifelse(is.finite(m$women_agecat), "women", ifelse(is.finite(m$sac_agecat), "sac", "other"))), useNA = "ifany"))
d <- m[is.finite(m$zn_gdl) & m$zn_gdl > 0, ]
d$group <- ifelse(is.finite(d$psc_agecat) | (is.finite(d$age_month) & d$age_month >= 6 & d$age_month < 60), "child",
           ifelse(is.finite(d$women_agecat), "women", ifelse(is.finite(d$sac_agecat), "sac", "other")))
d$lzn <- log(d$zn_gdl); d$pm <- as.integer(d$time_blood_draw == 2); d$date <- as.Date(d$date_interview)
d$doy <- as.numeric(d$date - as.Date("2015-12-01")); d$cluster <- as.character(d$cluster); d$district <- as.character(d$Admin1); d$ta <- paste(d$Admin1, d$Admin2)
d$sexf <- as.factor(d$sex); d$age <- ifelse(d$group == "child", d$age_month / 12, d$age_year)
cat("analysis rows by group:", paste(names(table(d$group)), table(d$group), collapse = " "), "| dated", sum(!is.na(d$date)), "\n")
rows <- list()
for (g in c("child", "sac", "women")) { x <- d[d$group == g, ]; if (nrow(x) < 100) next
  cat(sprintf("\n===== %s (n = %d, clusters %d, districts %d, TAs %d) =====\n", g, nrow(x), dplyr::n_distinct(x$cluster), dplyr::n_distinct(x$district), dplyr::n_distinct(x$ta)))
  # covariates only where they vary and are mostly observed (women: one sex; age missing for some groups)
  x$age[!is.finite(x$age)] <- NA; adj <- c("pm", if (mean(is.finite(x$age)) > 0.5) "age", if (nlevels(droplevels(x$sexf)) > 1) "sexf")
  x <- x[is.finite(x$age) | !("age" %in% adj), ]
  fx <- function(extra, dat) lmer(stats::as.formula(paste("lzn ~", paste(c(adj, extra), collapse = " + "), "+ (1 | district) + (1 | ta) + (1 | cluster)")), data = dat, REML = TRUE, control = lmerControl(check.conv.singular = "ignore"))
  # 1 individual-level: time of draw and date
  f0 <- lmer(lzn ~ 1 + (1 | district) + (1 | ta) + (1 | cluster), data = x, REML = TRUE, control = lmerControl(check.conv.singular = "ignore"))
  f1 <- fx(character(0), x)
  xd <- x[!is.na(x$doy), ]
  has_date <- nrow(xd) >= 100 && dplyr::n_distinct(xd$cluster) < nrow(xd)
  f2 <- if (has_date) fx(c("doy", "I(doy^2)"), xd) else f1
  co <- summary(f1)$coefficients; cat("time of draw (afternoon vs morning) on log zinc:", sprintf("%+.3f (t = %.1f)", co["pm", 1], co["pm", 3]), "=> ", sprintf("%+.1f%%", 100 * (exp(co["pm", 1]) - 1)), "\n")
  co2 <- summary(f2)$coefficients
  if (has_date) cat("with date (n =", nrow(xd), "): afternoon", sprintf("%+.3f (t = %.1f)", co2["pm", 1], co2["pm", 3]), "| date linear", sprintf("%+.4f/day (t = %.1f)", co2["doy", 1], co2["doy", 3]), "| quadratic t", sprintf("%.1f", co2["I(doy^2)", 3]), "\n") else { cat("no usable interview dates for this group (", nrow(xd), "dated rows )\n"); co2 <- rbind(co2, doy = c(NA, NA, NA)); if (!"I(doy^2)" %in% rownames(co2)) co2 <- rbind(co2, `I(doy^2)` = c(NA, NA, NA)) }
  vc <- function(f) { v <- as.data.frame(VarCorr(f)); stats::setNames(v$vcov, v$grp) }
  v0 <- vc(f0); v1 <- vc(f1); v2 <- vc(f2)
  cat("variance components (district | TA | cluster | resid):\n")
  cat(sprintf("  unadjusted            %.4f | %.4f | %.4f | %.4f\n", v0["district"], v0["ta"], v0["cluster"], v0["Residual"]))
  cat(sprintf("  + time of draw, age   %.4f | %.4f | %.4f | %.4f\n", v1["district"], v1["ta"], v1["cluster"], v1["Residual"]))
  cat(sprintf("  + date (dated subset) %.4f | %.4f | %.4f | %.4f\n", v2["district"], v2["ta"], v2["cluster"], v2["Residual"]))
  # 2 district-level: share of afternoon draws vs prevalence
  agg <- x |> group_by(district) |> summarise(n = dplyr::n(), prev = mean(zinc_def == 1, na.rm = TRUE), share_pm = mean(pm, na.rm = TRUE), mean_lzn = mean(lzn), date_med = as.numeric(stats::median(doy, na.rm = TRUE)), .groups = "drop")
  aggt <- x |> group_by(ta) |> summarise(n = dplyr::n(), prev = mean(zinc_def == 1, na.rm = TRUE), share_pm = mean(pm, na.rm = TRUE), mean_lzn = mean(lzn), .groups = "drop") |> filter(n >= 5)
  r1 <- suppressWarnings(stats::cor(agg$prev, agg$share_pm, method = "spearman")); r2 <- suppressWarnings(stats::cor(agg$mean_lzn, agg$share_pm, method = "spearman"))
  r3 <- suppressWarnings(stats::cor(aggt$prev, aggt$share_pm, method = "spearman")); r4 <- suppressWarnings(stats::cor(agg$prev, agg$date_med, method = "spearman"))
  cat(sprintf("district (n = %d): Spearman(prevalence, share afternoon) %+.2f | (mean log zinc, share afternoon) %+.2f | (prevalence, median date) %+.2f\n", nrow(agg), r1, r2, r4))
  cat(sprintf("TA (n = %d with >= 5): Spearman(prevalence, share afternoon) %+.2f | share afternoon range across districts %.2f-%.2f\n", nrow(aggt), r3, min(agg$share_pm), max(agg$share_pm)))
  # 3 how much district variance survives, as a share of the unadjusted
  rows[[length(rows) + 1L]] <- data.frame(group = g, n = nrow(x), n_dated = nrow(xd), beta_pm = co["pm", 1], t_pm = co["pm", 3], pct_pm = 100 * (exp(co["pm", 1]) - 1), beta_date = co2["doy", 1], t_date = co2["doy", 3],
    v_district_raw = v0["district"], v_district_adj = v1["district"], v_district_adj_date = v2["district"], v_ta_raw = v0["ta"], v_ta_adj = v1["ta"], v_cluster_raw = v0["cluster"], v_cluster_adj = v1["cluster"], v_cluster_adj_date = v2["cluster"], v_resid_raw = v0["Residual"],
    rho_prev_sharepm_district = r1, rho_lzn_sharepm_district = r2, rho_prev_sharepm_ta = r3, rho_prev_date_district = r4, share_pm_min = min(agg$share_pm), share_pm_max = max(agg$share_pm), stringsAsFactors = FALSE)
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "zinc_collection_artefact.csv"), row.names = FALSE)
cat("\n===== ZN-02 summary =====\n")
print(as.data.frame(R |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |> select(group, n, pct_pm, t_pm, t_date, v_district_raw, v_district_adj, v_district_adj_date, v_cluster_raw, v_cluster_adj_date, rho_prev_sharepm_district, rho_prev_sharepm_ta, share_pm_min, share_pm_max)), row.names = FALSE)
cat("\nDONE\n")
