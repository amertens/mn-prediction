# =============================================================================
# scripts/policy_deck/22_malawi_b12_maps.R
#
# The Malawi B12 surrogate-marker maps as a standalone figure for the 3-minute
# short-oral deck. Same construction as the full talk's `b12-surrogate` chunk
# (docs/slides/MN-proxy-full-talk-2026-09.qmd), lifted out so the lightning deck
# does not have to re-run a Quarto chunk to get it.
#
# Four panels, each as a within-country percentile so they share one scale:
# the survey's B12 deficiency, the deployment prediction for every district,
# the modelled anaemia surface, and the share of households eating fish.
# The point: anaemia carries NEGATIVE weight for B12 because it marks the
# fish-eating lakeshore. A place-marker, not a mechanism.
#
#   Rscript scripts/policy_deck/22_malawi_b12_maps.R
# -> results/figures/mnf15/figM_malawi_b12_surrogate.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUT <- "results/figures/mnf15"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TG <- read.csv("results/tables/protocol_v2/targets_v2.csv")
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S)); domain_of <- stats::setNames(MD$domain, MD$column)
B <- readRDS("dashboard/data/admin2_boundaries.rds")[["malawi"]]

t <- TG[TG$country == "Malawi" & TG$outcome == "women_b12" & is.finite(TG$y_prev),
        c("Admin1", "Admin2", "y_prev")]
s <- S[S$country == "Malawi", c("Admin1", "Admin2", "ihme_allanemia", "hces_any_fish", PREDS)]
m <- dplyr::left_join(s, t, by = c("Admin1", "Admin2"))

tr <- which(is.finite(m$y_prev))
Xr <- prep_predictors_v2(as.matrix(m[, PREDS])); Y <- .v2_logit(m$y_prev)
D  <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
m$pred <- .v2_expit(ARMS_V2[["domain_index"]](tr, seq_len(nrow(m)), Y, NULL, D,
                                              list(Admin1 = m$Admin1, y_nat = Y)))
cat("surveyed units:", length(tr), "of", nrow(m), "\n")

g <- dplyr::left_join(B, m[, c("Admin1", "Admin2", "y_prev", "pred", "ihme_allanemia", "hces_any_fish")],
                      by = c("Admin1", "Admin2"))
r1 <- cor(g$y_prev, g$ihme_allanemia, method = "spearman", use = "complete.obs")
r2 <- cor(g$y_prev, g$hces_any_fish,  method = "spearman", use = "complete.obs")
cat(sprintf("rho(B12 deficiency, modelled anaemia) = %.2f\nrho(B12 deficiency, fish eaten) = %.2f\n", r1, r2))

pct <- function(x) 100 * (rank(x, na.last = "keep") - 1) / (sum(is.finite(x)) - 1)
labs4 <- c(
  y_prev         = "Survey:\nB12 deficiency",
  pred           = "Model:\nevery district",
  ihme_allanemia = sprintf("Modelled anaemia\nrho = %.2f", r1),
  hces_any_fish  = sprintf("Fish eaten\nrho = %.2f", r2))
long <- do.call(rbind, lapply(names(labs4), function(v)
  data.frame(g[, c("Admin1", "Admin2")], value = pct(g[[v]]), panel = labs4[[v]])))
long$panel <- factor(long$panel, levels = labs4)

p <- ggplot(sf::st_as_sf(long)) +
  geom_sf(aes(fill = value), colour = "white", linewidth = 0.08) +
  facet_wrap(~ panel, ncol = 4) +
  scale_fill_viridis_c(na.value = "grey92", limits = c(0, 100),
                       name = "Percentile among Malawi districts (100 = highest)") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(1.9, "cm"),
        strip.text = element_text(face = "bold", size = 13, lineheight = 1.15),
        strip.clip = "off", panel.spacing.x = grid::unit(0.2, "cm"),
        plot.margin = margin(6, 10, 4, 10))
ggsave(file.path(OUT, "figM_malawi_b12_surrogate.png"), p, width = 11.0, height = 5.4, dpi = 200, bg = "white")
cat("wrote figM_malawi_b12_surrogate.png\n")
writeLines(sprintf("%.2f\t%.2f", r1, r2), file.path(OUT, "figM_malawi_rho.txt"))
