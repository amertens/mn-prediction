# =============================================================================
# scripts/policy_deck/37_mnf15_v5_national_all.R
#
# "But can proxy models recover national-level prevalence estimates?" for the v5
# MNF15 talk: script 35's figure with EVERY outcome the national model covers in
# the left panel (Andrew, 28 Sep: "show all outcomes, not just vitamin A; drop the
# text at the bottom").
#
# Left: our four countries, each predicted from national indicators with the
# country held out of the WHO VMNIS panel (SuperLearner: mean, ridge, lasso,
# elastic net, random forest; folds grouped by country).
#   vitamin A (children, women): national_levels_sl.csv, at each survey's year,
#     against the survey's own national prevalence (targets_v2.csv)
#   folate, B12 (women) and zinc (children): the leave-one-country-out predictions
#     for our countries' own VMNIS rows (national_vmnis_loco_sl_pred.csv), against
#     the VMNIS value. Iron has no VMNIS national panel in the national track.
# Right: every VMNIS panel, each country left out in turn, SuperLearner against
# the average of the other countries.
#
#   Rscript scripts/policy_deck/37_mnf15_v5_national_all.R
# -> results/figures/mnf15_v5/v5_national_all.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15_v5"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
th <- function(base = 15) theme_minimal(base_size = base) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
CT <- c(gambia = "The Gambia", ghana = "Ghana", sierraleone = "Sierra Leone", malawi = "Malawi",
        GMB = "The Gambia", GHA = "Ghana", SLE = "Sierra Leone", MWI = "Malawi")

va <- read.csv("results/tables/national_levels_sl.csv") |>
  transmute(nut = ifelse(grepl("^child", outcome), "Vitamin A, children", "Vitamin A, women"), ctry = CT[country], year,
            survey = survey_pp, model = sl_pp)
chk(nrow(va) == 8 && all(is.finite(va$model)), "eight vitamin A cells with a SuperLearner level")
pr <- read.csv("results/tables/national_vmnis_loco_sl_pred.csv") |>
  filter(iso3c %in% c("GMB", "GHA", "SLE", "MWI"), !grepl("^Vitamin A", panel), is.finite(sl), is.finite(observed)) |>
  group_by(panel, iso3c) |> slice_max(year, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(nut = sub(" \\| Preschool-age children", ", children", sub(" \\| Non-pregnant women \\(NPW\\)", ", women", panel)),
            ctry = CT[iso3c], year, survey = 100 * observed, model = 100 * sl)
NUT <- c("Vitamin A, children" = "vitamin A, children", "Vitamin A, women" = "vitamin A, women", "Folate, women" = "folate, women",
         "Vitamin B12, women" = "B12, women", "Zinc, children" = "zinc, children")
nl <- bind_rows(va, pr) |> mutate(lab = sprintf("%s, %s", ctry, NUT[nut]))
chk(!anyNA(nl$lab), "every outcome has a label")
cat("left panel:\n"); print(as.data.frame(nl), row.names = FALSE)
nl <- nl |> arrange(nut, survey) |> mutate(lab = factor(lab, levels = rev(unique(lab))))
long <- nl |> select(lab, nut, Survey = survey, `Predicted from national indicators` = model) |> pivot_longer(c(Survey, `Predicted from national indicators`), names_to = "what", values_to = "v")
pn1 <- ggplot(long, aes(v, lab)) +
  geom_line(aes(group = lab), colour = "grey70", linewidth = 1) +
  geom_point(aes(colour = what), size = 4.2) +
  geom_text(data = nl, aes(x = pmax(survey, model), y = lab, label = sprintf("%+.0f", model - survey)), hjust = -0.5, size = 4, colour = GREY) +
  scale_colour_manual(values = c(Survey = WARM, `Predicted from national indicators` = NAVY), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  labs(title = "Each country predicted without its own survey", x = "National prevalence of deficiency (%)", y = NULL) +
  th(14) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 14.5),
                 axis.text.y = element_text(size = 11.5, colour = INK), panel.grid.major.y = element_blank())
VL <- read.csv("results/tables/national_vmnis_loco_sl.csv") |> mutate(model = ifelse(is.na(model) | model == "", "null", model)) |>
  filter(model %in% c("null", "sl")) |>
  mutate(panel = sub(" \\| Preschool-age children", ", children", sub(" \\| Non-pregnant women \\(NPW\\)", ", women", panel)),
         panel = sprintf("%s\n(%d countries)", panel, n_countries),
         who = factor(c(null = "Average of the other countries", sl = "Predicted from national indicators")[model],
                      levels = c("Average of the other countries", "Predicted from national indicators")))
chk(sum(VL$model == "sl") == 5, "five VMNIS panels with a SuperLearner row")
pn2 <- ggplot(VL, aes(mae_pp, panel, fill = who)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f", mae_pp)), position = position_dodge(width = 0.8), hjust = -0.3, size = 4.2, colour = GREY) +
  scale_fill_manual(values = c("grey70", NAVY), name = NULL) + scale_x_continuous(limits = c(0, 23)) +
  labs(title = "WHO national surveys, each country left out", x = "Typical error in national prevalence (points)", y = NULL) +
  th(14) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 14.5),
                 axis.text.y = element_text(size = 12, colour = INK), panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "v5_national_all.png"), pn1 + pn2 + plot_layout(widths = c(1.15, 1)), width = 13, height = 6.2, dpi = 220, bg = "white")
cat("wrote v5_national_all.png\n")
cat(sprintf("left panel: %d country-outcomes; median absolute miss %.1f points, largest %.1f\n", nrow(nl), median(abs(nl$model - nl$survey)), max(abs(nl$model - nl$survey))))
print(as.data.frame(VL |> select(panel, model, mae_pp, spearman)), row.names = FALSE)
