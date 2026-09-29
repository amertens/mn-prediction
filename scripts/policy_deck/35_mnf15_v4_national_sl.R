# =============================================================================
# scripts/policy_deck/35_mnf15_v4_national_sl.R
#
# "But can proxy models recover national-level prevalence estimates?" redrawn
# with the SuperLearner (mean, ridge, lasso, elastic net, random forest) as the
# national model in BOTH panels, in place of ridge (left) and ridge + forest
# (right) of script 31's v4_national.png, and at the corrected survey years
# (Andrew, 28 September). Left: our four countries' national vitamin A
# prevalence, each predicted with the country held out of the WHO VMNIS panel.
# Right: every VMNIS panel, each country left out in turn, against the average
# of the other countries.
#
#   Rscript scripts/policy_deck/35_mnf15_v4_national_sl.R
# <- results/tables/national_levels_sl.csv, national_vmnis_loco_sl.csv  (scripts/covariates/19b)
# -> results/figures/mnf15_v4/v4_national_sl.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15_v4"
WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
th <- function(base = 15) theme_minimal(base_size = base) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

NL <- read.csv("results/tables/national_levels_sl.csv") |>
  mutate(ctry = recode(country, gambia = "The Gambia", ghana = "Ghana", sierraleone = "Sierra Leone", malawi = "Malawi"),
         pop = ifelse(grepl("^child", outcome), "children", "women"), lab = paste0(ctry, ", ", pop))
chk(nrow(NL) == 8 && all(is.finite(NL$sl_pp)), "eight vitamin A national cells with a SuperLearner level")
ordr <- NL |> arrange(pop, survey_pp) |> pull(lab)
nl <- NL |> select(lab, Survey = survey_pp, `Predicted from national indicators (SuperLearner)` = sl_pp) |>
  pivot_longer(-lab, names_to = "what", values_to = "v") |> mutate(lab = factor(lab, levels = ordr))
pn1 <- ggplot(nl, aes(v, lab)) +
  geom_line(aes(group = lab), colour = "grey70", linewidth = 1) +
  geom_point(aes(colour = what), size = 5) +
  geom_text(data = NL |> mutate(lab = factor(lab, levels = ordr)), aes(x = pmax(survey_pp, sl_pp), y = lab, label = sprintf("%+.0f", sl_pp - survey_pp)),
            hjust = -0.6, size = 4.6, colour = GREY) +
  scale_colour_manual(values = c(Survey = WARM, `Predicted from national indicators (SuperLearner)` = NAVY), name = NULL) +
  scale_x_continuous(limits = c(0, max(46, ceiling(max(nl$v) / 5) * 5 + 5))) +
  labs(title = "Vitamin A deficiency, each country predicted\nwithout its own survey", x = "National prevalence (%)", y = NULL) +
  th(15) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 15),
                 axis.text.y = element_text(size = 13.5, colour = INK))

VL <- read.csv("results/tables/national_vmnis_loco_sl.csv") |> mutate(model = ifelse(is.na(model) | model == "", "null", model)) |>
  filter(model %in% c("null", "sl")) |>
  mutate(panel = sub(" \\| Preschool-age children", ", children", sub(" \\| Non-pregnant women \\(NPW\\)", ", women", panel)),
         panel = sprintf("%s\n(%d countries)", panel, n_countries),
         who = factor(c(null = "Average of the other countries", sl = "SuperLearner on national indicators")[model],
                      levels = c("Average of the other countries", "SuperLearner on national indicators")))
chk(sum(VL$model == "sl") == 5, "five VMNIS panels with a SuperLearner row")
pn2 <- ggplot(VL, aes(mae_pp, panel, fill = who)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f", mae_pp)), position = position_dodge(width = 0.8), hjust = -0.3, size = 4.2, colour = GREY) +
  scale_fill_manual(values = c("grey70", NAVY), name = NULL) + scale_x_continuous(limits = c(0, 23)) +
  labs(title = "WHO VMNIS national surveys, each country\nleft out in turn", x = "Typical error in national prevalence (percentage points)", y = NULL) +
  th(15) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 15),
                 axis.text.y = element_text(size = 12.5, colour = INK), panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "v4_national_sl.png"),
       pn1 + pn2 + plot_layout(widths = c(1, 1.15)) +
         plot_annotation(caption = "Levels do not transport: that is why the model's ranking is anchored to a measured national figure (a small national sample, about 40 blood draws per nutrient).",
                         theme = theme(plot.caption = element_text(size = 13, colour = INK, hjust = 0))),
       width = 13, height = 5.8, dpi = 220, bg = "white")
cat("wrote v4_national_sl.png\n"); print(NL[, c("lab", "year", "survey_pp", "ridge_pp", "sl_pp", "null_pp")], row.names = FALSE)
print(as.data.frame(VL |> select(panel, model, mae_pp, spearman)), row.names = FALSE)
