# =============================================================================
# explore/scripts/19_cv_bands.R   [probe CV-01]
#
# The Bayesian SAE preprint flags estimates by coefficient of variation using
# the standard survey-statistics bands: CV < 16.6% publish unrestricted,
# 16.6-33.3% publish with caution, > 33.3% unreliable. This project's dashboard
# reports interval width and coverage separately and has no such gate.
#
# Applied to the DIRECT survey estimates the project publishes per district,
# using the measured effective n already in targets_v2 (DE-01), so it says how
# many district figures would clear an externally recognised standard.
#
#   Rscript explore/scripts/19_cv_bands.R -> explore/out/19_cv_bands.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
E <- exp_load(); TG <- E$TG
t <- TG[is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
p <- pmin(pmax(t$y_prev, 1e-6), 1 - 1e-6)
se <- sqrt(p * (1 - p) / pmax(t$n_eff, 1))
t$cv <- 100 * se / p
t$band <- cut(t$cv, c(-Inf, 16.6, 33.3, Inf),
              labels = c("unrestricted", "caution", "unreliable"))
exp_write(t[, c("country","outcome","Admin1","Admin2","y_prev","n_eff","cv","band")],
          "19_cv_bands")

cat("\n== share of published district estimates in each CV band ==\n")
tb <- as.data.frame.matrix(table(paste(t$country, t$outcome), t$band))
tb$n <- rowSums(tb)
for (v in c("unrestricted","caution","unreliable")) tb[[paste0("pct_", v)]] <- round(100*tb[[v]]/tb$n)
print(tb[order(-tb$pct_unreliable), c("n","pct_unrestricted","pct_caution","pct_unreliable")])
cat("\n== overall, and by country ==\n")
cat(sprintf("ALL district-outcome estimates: %d | unrestricted %.0f%% | caution %.0f%% | UNRELIABLE %.0f%%\n",
    nrow(t), 100*mean(t$band=="unrestricted"), 100*mean(t$band=="caution"),
    100*mean(t$band=="unreliable")))
for (c0 in sort(unique(t$country))) {
  s <- t[t$country==c0,]
  cat(sprintf("  %-13s n=%4d  unreliable %.0f%%  median CV %.0f%%\n",
      c0, nrow(s), 100*mean(s$band=="unreliable"), median(s$cv)))
}
