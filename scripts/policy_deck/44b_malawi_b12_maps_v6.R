# =============================================================================
# scripts/policy_deck/44b_malawi_b12_maps_v6.R
#
# Copy of scripts/policy_deck/22_malawi_b12_maps.R for the v6 MNF15 appendix
# slide "Causal driver or surrogate marker? Fish, anemia and B12 in Malawi",
# redrawn on the post-fix targets_v2.csv (27 September rebuild). The original
# is not edited. Differences from script 22:
#   * output: results/figures/mnf15_v6/a6_malawi_b12_surrogate.png (no rho file);
#   * panel titles and theme of the full-talk chunk that drew the old slide
#     (docs/slides/MN-proxy-full-talk-2026-09.qmd, chunk b12-surrogate);
#   * the deployment panel uses the headline predictor tiers (open +
#     public survey microdata, national constants dropped), as the dashboard
#     has since 27 September; script 22 used every column on disk, DHS
#     included. Both fits are computed and their agreement printed;
#   * the two correlations are also computed on the pre-fix targets and
#     printed old -> new. The pre-fix targets_v2.csv is in the snapshot taken
#     before the rebuild (C:/Users/andre/mn-prediction-snapshots/
#     pre_rebuild_2026-09-27); the copy in protocol_v2_pre_RR11_20260927 was
#     made after script 01 had rewritten the targets, so it is already post-fix.
#
# Four panels, each the within-country percentile of its variable: the
# survey's B12 deficiency, the deployment prediction for every district, the
# modelled anaemia surface (IHME, women), and the share of households eating
# fish in the last week (IHS4). rho = Spearman with B12 deficiency across the
# surveyed districts.
#
#   Rscript scripts/policy_deck/44b_malawi_b12_maps_v6.R
# -> results/figures/mnf15_v6/a6_malawi_b12_surrogate.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R"); if (!exists("is_water_admin2", mode = "function")) source("R/admin2_key_hygiene.R")
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)

TG  <- read.csv("results/tables/protocol_v2/targets_v2.csv")
PRE_TG <- "C:/Users/andre/mn-prediction-snapshots/pre_rebuild_2026-09-27/results/tables/protocol_v2/targets_v2.csv"
TGo <- if (file.exists(PRE_TG)) read.csv(PRE_TG) else TG[0, ]   # without the snapshot the old values print as NA
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS_ALL <- intersect(MD$column, names(S)); domain_of <- stats::setNames(MD$domain, MD$column)
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
PREDS <- drop_near_outcome_v2(PREDS_ALL, MD)
B <- readRDS("dashboard/data/admin2_boundaries.rds")[["malawi"]]

s <- S[S$country == "Malawi", c("Admin1", "Admin2", "ihme_allanemia", "hces_any_fish", PREDS_ALL)]
s <- s[!is_water_admin2(s$Admin2), ]
chk(!anyDuplicated(paste(s$Admin1, s$Admin2)) && !anyDuplicated(paste(B$Admin1, B$Admin2)), "one row per district")
cell <- function(tg) tg[tg$country == "Malawi" & tg$outcome == "women_b12" & is.finite(tg$y_prev), c("Admin1", "Admin2", "y_prev")]
m  <- dplyr::left_join(s, cell(TG), by = c("Admin1", "Admin2"))
mo <- dplyr::left_join(s[, c("Admin1", "Admin2", "ihme_allanemia", "hces_any_fish")], cell(TGo), by = c("Admin1", "Admin2"))

# the deployment fit: the index on all surveyed districts, applied to every district
tr <- which(is.finite(m$y_prev)); Y <- .v2_logit(m$y_prev)
deploy <- function(cols) {
  Xr <- prep_predictors_v2(as.matrix(m[, cols]))
  D  <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
  .v2_expit(ARMS_V2[["domain_index"]](tr, seq_len(nrow(m)), Y, NULL, D, list(Admin1 = m$Admin1, y_nat = Y)))
}
m$pred <- deploy(PREDS)
pred_all <- deploy(PREDS_ALL)
cat(sprintf("surveyed districts: %d of %d (pre-fix: %d)\n", length(tr), nrow(m), sum(is.finite(mo$y_prev))))
cat(sprintf("deployment fit, headline tiers vs every column (script 22): Spearman %.3f over all %d districts\n",
            cor(m$pred, pred_all, method = "spearman"), nrow(m)))

sp <- function(d, v) cor(d$y_prev, d[[v]], method = "spearman", use = "complete.obs")
r1 <- sp(m, "ihme_allanemia"); r2 <- sp(m, "hces_any_fish")
cat(sprintf("rho(B12 deficiency, modelled anaemia): old %.2f -> new %.2f (n = %d)\n", sp(mo, "ihme_allanemia"), r1, sum(is.finite(m$y_prev) & is.finite(m$ihme_allanemia))))
cat(sprintf("rho(B12 deficiency, fish eaten):       old %.2f -> new %.2f (n = %d)\n", sp(mo, "hces_any_fish"), r2, sum(is.finite(m$y_prev) & is.finite(m$hces_any_fish))))
cat(sprintf("rho(B12 deficiency, deployment fit, surveyed districts, in-sample): %.2f\n", sp(m, "pred")))

g <- dplyr::left_join(B, m[, c("Admin1", "Admin2", "y_prev", "pred", "ihme_allanemia", "hces_any_fish")], by = c("Admin1", "Admin2"))
chk(abs(cor(g$y_prev, g$ihme_allanemia, method = "spearman", use = "complete.obs") - r1) < 1e-12, "rho on the map layer equals rho on the district table")
cat(sprintf("polygons drawn: %d, with a prediction: %d\n", nrow(g), sum(is.finite(g$pred))))

pct <- function(x) 100 * (rank(x, na.last = "keep") - 1) / (sum(is.finite(x)) - 1)   # within-country percentile, one shared scale
labs4 <- c(y_prev = "Women's B12 deficiency\nsurvey 2015-16\n ", pred = "Model-predicted deficiency\nevery district, fitted on\nthe surveyed ones",
           ihme_allanemia = sprintf("Modelled anaemia\nin women (IHME)\nrho with B12 def. = %.2f", r1),
           hces_any_fish = sprintf("Households eating\nfish last week (IHS4)\nrho with B12 def. = %.2f", r2))
long <- do.call(rbind, lapply(names(labs4), function(v) data.frame(g[, c("Admin1", "Admin2")], value = pct(g[[v]]), panel = labs4[[v]])))
long$panel <- factor(long$panel, levels = labs4)
p <- ggplot(sf::st_as_sf(long)) + geom_sf(aes(fill = value), colour = "white", linewidth = 0.1) + facet_wrap(~ panel, ncol = 4) +
  scale_fill_viridis_c(na.value = "grey92", name = "Percentile among Malawi districts (100 = highest)", limits = c(0, 100)) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(1.6, "cm"), strip.text = element_text(face = "bold", size = 10, lineheight = 1.05),
        strip.clip = "off", panel.spacing.x = grid::unit(1.0, "cm"))
ggsave(file.path(OUT, "a6_malawi_b12_surrogate.png"), p, width = 12, height = 5.2, dpi = 220, bg = "white")   # 10.61 x 4.6 box
cat("wrote a6_malawi_b12_surrogate.png\nDONE\n")
