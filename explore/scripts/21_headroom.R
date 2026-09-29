# =============================================================================
# explore/scripts/21_headroom.R   [probe HR-02]
# How much of the attainable district signal is the model already getting?
# r_max = sqrt(reliability); r_share = achieved / r_max.
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
R <- read.csv(file.path(EXP_OUT, "16_reliability_all.csv"), stringsAsFactors = FALSE)
R$country[R$country == "Sierra Leone"] <- "SierraLeone"
BM <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/benchmarks_v2_cells.csv"),
               stringsAsFactors = FALSE)
b <- BM[BM$estimand == "infill" & BM$target == "level" & BM$arm == "domain_index",
        c("country","outcome","spearman")]
names(b)[3] <- "achieved"
m <- merge(R[, c("country","outcome","r_full_respondent","districts")], b,
           by = c("country","outcome"))
m$r_max <- sqrt(pmax(m$r_full_respondent, 0))
m$r_share <- ifelse(m$r_max > 0.05, m$achieved / m$r_max, NA)
m$headroom <- m$r_max - m$achieved
m <- m[order(-m$headroom), ]
cat("r_max = sqrt(split-half reliability) is the most a predictor could reach.\n")
cat("r_share = achieved / r_max. headroom = what is left.\n\n")
print(m[, c("country","outcome","districts","r_full_respondent","r_max","achieved","r_share","headroom")],
      row.names = FALSE, digits = 3)
cat(sprintf("\nmedian r_share %.2f | median headroom %.3f | cells %d\n",
    median(m$r_share, na.rm=TRUE), median(m$headroom, na.rm=TRUE), nrow(m)))
cat(sprintf("cells already above 80%% of attainable: %d | below 50%%: %d\n",
    sum(m$r_share > 0.8, na.rm=TRUE), sum(m$r_share < 0.5, na.rm=TRUE)))
cat("\n== is the remaining headroom explained by the NUMBER of districts? ==\n")
ct <- suppressWarnings(cor.test(m$districts, m$r_share, method="spearman"))
cat(sprintf("spearman(n districts, r_share) = %+.3f  p = %.3f\n", ct$estimate, ct$p.value))
agg <- aggregate(cbind(r_share, headroom, r_max, achieved) ~ districts, data=m,
                 FUN=function(z) round(mean(z, na.rm=TRUE),3))
print(agg, row.names = FALSE)
exp_write(m, "21_headroom")
