# =============================================================================
# explore/scripts/15_seasonal_offset.R   [probe SO-01]
#
# IS THE CROSS-COUNTRY LEVEL OFFSET PARTLY SEASONAL?
#
# RV-01 and LT-01 established that the between-country share of biomarker level
# variance is large for iron (0.71-0.80) and small for B12 (0.14), and AS-01
# ruled out the assay (the same VitMin ELISA in all four surveys). The offset is
# recorded as unexplained. This asks whether it tracks WHEN each survey ran.
#
# Mechanism: three of the four surveys sit in the dry / post-harvest season and
# Malawi sits in the lean season (FW-01). Post-harvest status should be better
# than lean-season status for nutrients that track recent intake. If the country
# offsets order by seasonal position, part of what looks like a survey artefact
# is real seasonal variation in the population - which is a different problem
# with a different fix (report the season, do not adjust it away).
#
# SEASONAL POSITION, from the Earth Engine monthly extraction (probe TM-01b):
# for each country, the NDVI anomaly at the fieldwork month, and the number of
# months from the fieldwork month to that country's NDVI peak (its growing
# season maximum). Both are computed from the cluster monthly stack, averaged
# over the country's clusters.
#
# POWER. Four countries, three for B12 and folate. This design cannot support a
# test; it can only show whether the ordering is consistent with the mechanism.
# Reported as descriptive, with the correlation given for completeness and the
# n printed beside it.
#
#   Rscript explore/scripts/15_seasonal_offset.R
# -> explore/out/15_seasonal_offset.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

E <- exp_load()
TG <- E$TG

G <- read.csv(file.path(EXP_ROOT, "explore/out/10_gee_cluster_monthly.csv"),
              stringsAsFactors = FALSE)
G$cal_month <- as.integer(substr(G$year_month, 6, 7))

# ── seasonal position of each country's fieldwork ───────────────────────────
# NDVI climatology by calendar month, per country, from the 36-month stack
nd <- G[G$layer == "ndvi", ]
clim <- aggregate(value ~ country + cal_month, data = nd, FUN = mean, na.rm = TRUE)
peak <- do.call(rbind, lapply(split(clim, clim$country), function(d)
  data.frame(country = d$country[1], peak_month = d$cal_month[which.max(d$value)],
             stringsAsFactors = FALSE)))

# fieldwork month = the month of lag 0 for that country's clusters
fw <- nd[nd$lag_months == 0, ]
fwm <- do.call(rbind, lapply(split(fw, fw$country), function(d) {
  mm <- as.integer(names(sort(table(d$cal_month), decreasing = TRUE))[1])
  # NDVI at fieldwork relative to that country's own annual range
  cc <- clim[clim$country == d$country[1], ]
  rng <- range(cc$value, na.rm = TRUE)
  at <- cc$value[cc$cal_month == mm]
  data.frame(country = d$country[1], fieldwork_month = mm,
             ndvi_at_fieldwork_rel = (at - rng[1]) / max(rng[2] - rng[1], 1e-9),
             stringsAsFactors = FALSE)
}))
S <- merge(fwm, peak, by = "country")
# circular months from the growing-season peak to the fieldwork month
S$months_since_peak <- (S$fieldwork_month - S$peak_month) %% 12
message("seasonal position per country:")
print(S, row.names = FALSE)

# ── country offsets per outcome ─────────────────────────────────────────────
rows <- list()
for (on in unique(TG$outcome)) {
  t <- TG[TG$outcome == on & is.finite(TG$y_level), ]
  if (dplyr::n_distinct(t$country) < 3) next
  mu <- tapply(t$y_level, t$country, mean)
  sdw <- tapply(t$y_level, t$country, stats::sd)
  d <- data.frame(outcome = on, country = names(mu),
                  country_mean_level = as.numeric(mu),
                  within_sd = as.numeric(sdw), stringsAsFactors = FALSE)
  d$offset <- d$country_mean_level - mean(d$country_mean_level)
  rows[[on]] <- merge(d, S, by = "country")
}
O <- dplyr::bind_rows(rows); exp_write(O, "15_seasonal_offset")

cat("\n== country offsets against seasonal position ==\n")
cat("   (y_level is NEGATED log biomarker: higher = more deficient)\n\n")
for (on in unique(O$outcome)) {
  d <- O[O$outcome == on, ]
  d <- d[order(d$months_since_peak), ]
  cat(sprintf("-- %s (%d countries) --\n", on, nrow(d)))
  print(d[, c("country", "fieldwork_month", "months_since_peak",
              "ndvi_at_fieldwork_rel", "offset")], row.names = FALSE, digits = 3)
  if (nrow(d) >= 3) {
    r1 <- suppressWarnings(stats::cor(d$months_since_peak, d$offset, method = "spearman"))
    r2 <- suppressWarnings(stats::cor(d$ndvi_at_fieldwork_rel, d$offset, method = "spearman"))
    cat(sprintf("   spearman(months since peak, offset) = %+.2f | spearman(NDVI at fieldwork, offset) = %+.2f   [n = %d, descriptive only]\n\n",
                r1, r2, nrow(d)))
  }
}

cat("== pooled across outcomes (still only 4 countries of information) ==\n")
pool <- aggregate(offset ~ country + months_since_peak + ndvi_at_fieldwork_rel,
                  data = O, FUN = mean)
print(pool[order(pool$months_since_peak), ], row.names = FALSE, digits = 3)
cat(sprintf("\nspearman(months since peak, mean offset) = %+.2f over %d countries\n",
            suppressWarnings(stats::cor(pool$months_since_peak, pool$offset,
                                        method = "spearman")), nrow(pool)))
cat("\nFOUR COUNTRIES. This cannot be a test. Read the ordering, not the number.\n")
