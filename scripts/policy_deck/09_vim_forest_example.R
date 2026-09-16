# =============================================================================
# scripts/policy_deck/09_vim_forest_example.R   [VF-01, 2026-09-16]
#
# A VARIABLE-IMPORTANCE FOREST PLOT FOR ONE COUNTRY AND ONE OUTCOME
#
# The index is linear in the rank-normalised predictors (R/protocol_v2_importance.R),
# so every column has a standardised weight beta_std in a fit. This script fits
# the index on one country's surveyed districts, then refits it on B bootstrap
# resamples of those districts (domain PCs re-oriented on each resample, weights
# re-learned), and draws the twenty columns with the largest median |weight| as
# a forest plot: point = median over resamples, bar = 5th to 95th percentile,
# hollow marker = the weight in the full fit. The level target is the NEGATED log
# concentration (scripts/protocol_v2/01_build_targets_v2.R: higher = worse status),
# so a positive weight means districts with more of the predictor have worse status.
#
# Cell: by default the best in-fill index cell on the biomarker level in
# results/tables/protocol_v2/benchmarks_v2_cells.csv (RR-11: Malawi women's B12,
# Spearman 0.70; The Gambia child vitamin A is 0.69). Override with
#   VIM_COUNTRY=Gambia VIM_OUTCOME=child_vitA VIM_TARGET=level VIM_B=300
# Tiers follow V2_PREDICTOR_TIERS (default here: open,survey_public, the headline).
#
#   Rscript scripts/policy_deck/09_vim_forest_example.R
# -> results/figures/policy_deck/fig12_vim_forest_<country>_<outcome>.png
#    results/tables/policy_deck/vim_forest_<country>_<outcome>.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R"); source("R/predictor_plain_names.R")
P2 <- "results/tables/protocol_v2"; FDIR <- "results/figures/policy_deck"; TDIR <- "results/tables/policy_deck"
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE); dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)
B <- as.integer(Sys.getenv("VIM_B", "300")); TOP <- 20L; set.seed(20260916L)

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)

# ── which cell ────────────────────────────────────────────────────────────────
TARGET <- Sys.getenv("VIM_TARGET", "level")
if (nzchar(Sys.getenv("VIM_COUNTRY"))) { CN <- Sys.getenv("VIM_COUNTRY"); ON <- Sys.getenv("VIM_OUTCOME") } else {
  BC <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"), stringsAsFactors = FALSE)
  best <- BC |> filter(estimand == "infill", arm == "domain_index", target == TARGET) |> arrange(desc(spearman)) |> slice(1)
  CN <- best$country; ON <- best$outcome
  cat(sprintf("best in-fill cell on %s: %s %s, Spearman %.3f\n", TARGET, CN, ON, best$spearman))
}
ycol <- if (TARGET == "prev") "y_prev" else "y_level"; wcol <- if (TARGET == "prev") "n_eff" else "n_eff_cont"
t <- TG[TG$country == CN & TG$outcome == ON, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
m <- inner_join(t, S[S$country == CN, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
X <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); y <- if (TARGET == "prev") .v2_logit(m[[ycol]]) else m[[ycol]]
n <- nrow(X); cat(sprintf("%s %s: %d districts, %d predictors after per-country preparation\n", CN, ON, n, ncol(X)))

# ── full fit and bootstrap refits ─────────────────────────────────────────────
D0 <- domain_representation_v2(X, domain_of)
full <- index_importance_v2(seq_len(n), y, X, D0)$columns
boot <- sapply(seq_len(B), function(b) {
  tr <- sample.int(n, n, replace = TRUE)
  Db <- domain_representation_v2(X, domain_of, sign_rows = tr)
  im <- tryCatch(index_importance_v2(tr, y, X, Db)$columns, error = function(e) NULL)
  if (is.null(im)) return(rep(NA_real_, ncol(X)))
  im$beta_std[match(colnames(X), im$column)]
})
rownames(boot) <- colnames(X)
summ <- data.frame(column = colnames(X), full = full$beta_std[match(colnames(X), full$column)],
                   median = apply(boot, 1, median, na.rm = TRUE), lo = apply(boot, 1, quantile, 0.05, na.rm = TRUE), hi = apply(boot, 1, quantile, 0.95, na.rm = TRUE),
                   sign_agree = apply(boot, 1, function(v) mean(sign(v) == sign(median(v, na.rm = TRUE)), na.rm = TRUE)), stringsAsFactors = FALSE)
summ$domain <- unname(domain_of[summ$column]); summ$source <- MD$source[match(summ$column, MD$column)]
summ$plain <- unname(PLAIN[summ$column]); summ$label <- ifelse(is.na(summ$plain), clean_code(summ$column), summ$plain)   # codes stay in the CSV
summ <- summ |> arrange(desc(abs(median))) |> mutate(rank = row_number())
write.csv(summ, file.path(TDIR, sprintf("vim_forest_%s_%s.csv", CN, ON)), row.names = FALSE)

# ── figure ────────────────────────────────────────────────────────────────────
olab <- c(child_vitA = "child vitamin A", women_vitA = "women's vitamin A", child_iron = "child iron", women_iron = "women's iron", women_folate = "women's folate", women_b12 = "women's B12")
d <- summ |> slice_head(n = TOP) |> mutate(label = factor(label, levels = rev(label)))
pal <- c("#1f4e79", "#6baed6", "#8c510a", "#33a02c", "#b15928", "#6a3d9a", "#e08214", "#0F7B8A", "#cab2d6", "#fdbf6f", "#d7191c", "#9e9ac8", "#7fcdbb", "#fdae61")
doms <- sort(unique(d$domain)); names(pal) <- NULL; cols <- stats::setNames(rep(pal, length.out = length(doms)), doms)
g <- ggplot(d, aes(median, label, colour = domain)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, orientation = "y") +
  geom_point(size = 3) + geom_point(aes(x = full), shape = 1, size = 3.2, stroke = 0.9) +
  scale_colour_manual(values = cols) +
  labs(x = sprintf("Standardised weight in the index
(median and 5th to 95th percentile over %d bootstrap refits of %d districts; hollow marker = full fit)", B, n), y = NULL, colour = NULL,
       caption = sprintf("%s, %s, %s. Positive: districts with more of the predictor have worse status (%s).",
                         sub("SierraLeone", "Sierra Leone", CN), olab[[ON]], if (TARGET == "prev") "deficiency prevalence" else "biomarker level", if (TARGET == "prev") "higher prevalence" else "a lower concentration; the level target is the negated log concentration")) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", legend.text = element_text(size = 9), legend.key.height = grid::unit(10, "pt"), panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size = 10)) +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE))
f <- file.path(FDIR, sprintf("fig12_vim_forest_%s_%s.png", CN, ON))
ggsave(f, g, width = 12, height = 7.8, dpi = 200, bg = "white")
cat(sprintf("\n%s\n", f))
print(as.data.frame(d |> transmute(rank, column, domain, median = round(median, 3), lo = round(lo, 3), hi = round(hi, 3), sign_agree = round(sign_agree, 2))), row.names = FALSE)
cat("DONE\n")
