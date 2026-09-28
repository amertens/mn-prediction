# =============================================================================
# scripts/protocol_v2/72_plot_person_level_honest.R   [IL-02e]
# Draws the two person-level figures from il02_honest_person_level.csv (71):
#   Brier skill for the deficiency flags (the January figure, updated and scored out of fold)
#   MSE skill for the log concentrations
# Palette: the dataviz reference categorical slots 1-4 in fixed order (validated:
# adjacent CVD dE >= 9.1, normal-vision >= 22.9 on the light surface); slots 3-4
# sit below 3:1 contrast, so every row carries its set name as text. The
# ceiling is a bound, not a model, and is drawn in neutral ink.
#
#   Rscript scripts/protocol_v2/72_plot_person_level_honest.R
# -> results/figures/il02_person_level_brier_honest.png       all four model families
# -> results/figures/il02_person_level_mse_honest.png
# -> results/figures/il02_person_level_brier_proxy_only.png   index and proxies only (no survey answers)
# -> results/figures/il02_person_level_mse_proxy_only.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
M <- read.csv("results/tables/protocol_v2/il02_honest_person_level.csv", stringsAsFactors = FALSE)
SRC <- read.csv("results/tables/protocol_v2/il02_honest_survey_columns.csv", stringsAsFactors = FALSE) |> distinct(outcome, concentration_source)
FIGDIR <- "results/figures"; dir.create(FIGDIR, showWarnings = FALSE)
COUNTRY <- M$country[1]

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8983"; GRID <- "#e9e8e4"; SURF <- "#fcfcfb"
MODEL_SETS <- c(index = "Index", survey = "Survey only", proxies = "Proxies", both = "Survey + proxies")
COLS <- c("Index" = "#2a78d6", "Survey only" = "#eb6834", "Proxies" = "#1baf7a", "Survey + proxies" = "#eda100")

prep <- function(tp, nutrient_names, ceiling_label, keep = names(MODEL_SETS)) {
  d <- M[M$type == tp & M$set %in% c(keep, "ceiling"), ]
  sets <- c(MODEL_SETS[keep], ceiling = ceiling_label)
  d$population <- ifelse(grepl("^child", d$outcome), "Children", "Women")
  d$nutrient <- factor(nutrient_names[sub("^(child|women)_", "", d$outcome)], levels = unique(nutrient_names))
  d$set_lab <- factor(sets[d$set], levels = sets)
  lev <- as.vector(t(outer(c("Children", "Women"), sets, function(p, s) paste0(p, " \u2013 ", s))))
  d$row_lab <- factor(paste0(d$population, " \u2013 ", d$set_lab), levels = rev(lev))
  d |> mutate(across(c(skill, lo, hi), ~ 100 * .x))
}

draw <- function(d, ceiling_label, xlab, title, subtitle, caption, file, height = 7) {
  # colour follows the model, never its position: the proxy-only figure reuses the full figure's hues
  cols <- c(COLS, setNames(INK2, ceiling_label)); shapes <- c(setNames(rep(16, 4), names(COLS)), setNames(5, ceiling_label))
  xr <- range(c(d$lo, d$hi, 0), na.rm = TRUE); xr <- xr + c(-1, 1) * 0.04 * diff(xr)
  p <- ggplot(d, aes(x = skill, y = row_lab, colour = set_lab)) +
    geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.5) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.7, orientation = "y") +
    geom_point(aes(shape = set_lab), size = 2.8, stroke = 1.1, fill = SURF) +
    facet_wrap(~ nutrient, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    scale_shape_manual(values = shapes, name = NULL, drop = TRUE) +
    scale_x_continuous(limits = xr, breaks = scales::breaks_pretty(6), labels = function(x) paste0(x, "%")) +
    labs(x = xlab, y = NULL, title = title, subtitle = subtitle, caption = caption) +
    theme_minimal(base_size = 12) +
    theme(plot.background = element_rect(fill = SURF, colour = NA), panel.background = element_rect(fill = SURF, colour = NA),
          text = element_text(colour = INK), axis.text = element_text(colour = INK2), axis.text.y = element_text(size = 10.5),
          strip.text = element_text(face = "bold", size = 13, hjust = 0, colour = INK),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(colour = INK2, size = 11),
          plot.caption = element_text(colour = INK2, size = 9, hjust = 0), plot.title.position = "plot", plot.caption.position = "plot",
          legend.position = "bottom", legend.text = element_text(colour = INK2, size = 10.5), panel.spacing.x = unit(1.4, "lines"))
  ggsave(file, p, width = 11.5, height = height, dpi = 200, bg = SURF)
  invisible(p)
}

SUB <- paste0("Every model is scored on districts it never saw (5 folds blocked by district, averaged over 5 fold draws). ",
              "Bars: 95% cluster-bootstrap intervals.")
wv <- M[M$type == "bin" & M$outcome == "women_vitA" & M$set == "index", ]
cases_wva <- if (nrow(wv)) round(wv$n * wv$mean_outcome) else NA

NUT_B <- c(vitA = "Vitamin A deficiency", iron = "Iron deficiency", folate = "Folate deficiency", b12 = "Vitamin B12 deficiency")
NUT_C <- c(vitA = "Retinol-binding protein", iron = "Ferritin", folate = "Folate", b12 = "Vitamin B12")
ceil_b <- "Ceiling: true district rate"; ceil_c <- "Ceiling: true district mean"
CAP_CEIL_B <- "Ceiling: the between-district share of person-level variance. A model that knew every district's true prevalence would reach it; no model built on district data can pass it."
CAP_CEIL_C <- "Ceiling: the between-district share of person-level variance. A model that knew every district's true mean would reach it; no model built on district data can pass it."
XB <- "Improvement over predicting the survey prevalence for everyone (Brier skill score)"
XC <- "Reduction in squared error vs predicting the survey mean for everyone (log concentration)"
plain_src <- function(s) dplyr::case_when(grepl("RBP", s) ~ "retinol-binding protein BRINDA-adjusted",
                                          grepl("FerrAdjThurn", s) ~ "ferritin Thurnham-adjusted",
                                          grepl("Ferr", s) ~ "ferritin", grepl("Folate", s, ignore.case = TRUE) ~ "serum folate",
                                          grepl("B12", s) ~ "serum B12", TRUE ~ s)
CAP_SRC <- paste0("Concentrations on the log scale: ", paste(unique(plain_src(SRC$concentration_source)), collapse = "; "), ".")
CAP_JAN <- paste0("Women's vitamin A: ", cases_wva, " cases. The January 2026 version of this figure (20\u201367%) scored the models on their own training data.")

# all four model families
draw(prep("bin", NUT_B, ceil_b), ceil_b, XB, paste0("Predicting which individuals are deficient, ", COUNTRY, " 2017 (current data and models)"), SUB,
     paste0(CAP_CEIL_B, "\n", CAP_JAN), file.path(FIGDIR, "il02_person_level_brier_honest.png"))
draw(prep("cont", NUT_C, ceil_c), ceil_c, XC, paste0("Predicting individual biomarker concentrations, ", COUNTRY, " 2017 (current data and models)"), SUB,
     paste0(CAP_SRC, "\n", CAP_CEIL_C), file.path(FIGDIR, "il02_person_level_mse_honest.png"))
# proxy-only: the district-level models the project deploys (no respondent's own survey answers)
PX <- c("index", "proxies")
draw(prep("bin", NUT_B, ceil_b, PX), ceil_b, XB, paste0("Predicting which individuals are deficient from district proxies alone, ", COUNTRY, " 2017"), SUB,
     paste0(CAP_CEIL_B, "\nIndex: the PCA domain index, calibrated to respondents. Proxies: SuperLearner on the district's domain components. Women's vitamin A: ", cases_wva, " cases."),
     file.path(FIGDIR, "il02_person_level_brier_proxy_only.png"), height = 5.6)
draw(prep("cont", NUT_C, ceil_c, PX), ceil_c, XC, paste0("Predicting individual biomarker concentrations from district proxies alone, ", COUNTRY, " 2017"), SUB,
     paste0(CAP_SRC, "\n", CAP_CEIL_C), file.path(FIGDIR, "il02_person_level_mse_proxy_only.png"), height = 5.6)
cat("figures written\n")
