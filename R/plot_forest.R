# =============================================================================
# R/plot_forest.R
#
# Forest plots of OUT-OF-SAMPLE accuracy, for the funder brief, the manuscript
# and the slide deck. One generic panel builder plus four thin wrappers, each
# reading a results table the pipeline already writes -- no recomputation, so a
# figure can never disagree with the table it came from.
#
#   forest_national()   national prevalence error vs the survey's own 95% CI
#   forest_targeting()  programme reach at the top 20% of districts
#   forest_transport()  cross-country district rank r, bootstrap CI
#   forest_within()     within-country district rank r, bootstrap CI
#   table_anchor_arms() the admin-1 anchoring comparison (a table, not a plot)
#
# DESIGN NOTES
#   - Significance is encoded in COLOUR, not in a separate column, so the eye
#     picks out what survived a permutation null without reading numbers.
#   - Every panel draws its own reference line (0 for correlations, 1.0x for
#     lift, 0 pp for prevalence error) because "better" means something
#     different in each and a shared axis convention would mislead.
#   - Rows are ordered by effect size, not alphabetically: a forest plot's job
#     is to show the spread.
#   - Colours come from one palette so the deck, the report and the artifact
#     agree. ACCENT is the clinical teal used in the published brief.
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

# Local definition so this file can be sourced on its own (for a qmd or a
# slide deck) without pulling in the whole R/ directory.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

FOREST_PAL <- c(
  accent = "#0F6E66",  # significant / adds value
  null   = "#9AA6AE",  # not significant
  neg    = "#9E3A2C",  # outside the interval
  ref    = "#9C6206",  # reference line
  ink    = "#13171B",
  ink2   = "#48535E",
  rule   = "#D9E0E5"
)

OUTCOME_LABEL <- c(
  child_vitA = "Child vitamin A", women_vitA = "Women vitamin A",
  child_iron = "Child iron",      women_iron = "Women iron",
  child_zinc = "Child zinc",      women_zinc = "Women zinc",
  women_folate = "Women folate",  women_b12  = "Women B12")

.lab_outcome <- function(x) ifelse(x %in% names(OUTCOME_LABEL),
                                   OUTCOME_LABEL[x], x)
.cap_country <- function(x) {
  m <- c(gambia = "Gambia", ghana = "Ghana",
         sierraleone = "Sierra Leone", malawi = "Malawi")
  ifelse(tolower(x) %in% names(m), m[tolower(x)], x)
}

theme_forest <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = FOREST_PAL[["rule"]], linewidth = .3),
      axis.title.y       = element_blank(),
      axis.text.y        = element_text(colour = FOREST_PAL[["ink2"]], hjust = 1),
      axis.title.x       = element_text(colour = FOREST_PAL[["ink2"]],
                                        margin = margin(t = 8)),
      plot.title         = element_text(face = "bold", size = base_size + 3,
                                        colour = FOREST_PAL[["ink"]]),
      plot.subtitle      = element_text(colour = FOREST_PAL[["ink2"]],
                                        size = base_size - .5,
                                        margin = margin(b = 10)),
      plot.caption       = element_text(colour = FOREST_PAL[["ink2"]],
                                        size = base_size - 2, hjust = 0),
      legend.position    = "none",
      plot.margin        = margin(10, 14, 8, 8))
}

#' Generic forest panel.
#'
#' @param d data with columns: label, est, lo, hi (lo/hi may be NA), flag
#'   (TRUE = highlight), and optionally `bad` (TRUE = draw in the warning
#'   colour, used for "outside the interval")
#' @param ref reference line position
#' @param xlab,title,subtitle,caption text
forest_panel <- function(d, ref = 0, xlab = NULL, title = NULL,
                         subtitle = NULL, caption = NULL) {
  d <- d %>%
    mutate(label = factor(label, levels = rev(label)),
           col = dplyr::case_when(
             isTRUE(.data$bad) ~ FOREST_PAL[["neg"]],
             .data$flag        ~ FOREST_PAL[["accent"]],
             TRUE              ~ FOREST_PAL[["null"]]))
  p <- ggplot(d, aes(x = est, y = label, colour = I(col))) +
    geom_vline(xintercept = ref, colour = FOREST_PAL[["ref"]],
               linetype = "22", linewidth = .5)
  if (any(is.finite(d$lo) & is.finite(d$hi)))
    # geom_errorbar(orientation = "y"), not geom_errorbarh(): the latter is
    # deprecated in ggplot2 4.0 and warns on every panel.
    p <- p + geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                           width = .28, linewidth = .55, na.rm = TRUE,
                           alpha = .75)
  # Wrap: an unwrapped subtitle runs off the canvas at these widths.
  wrap <- function(x, n) if (is.null(x)) NULL else
    paste(strwrap(x, width = n), collapse = "
")
  p +
    geom_point(size = 2.1) +
    labs(x = xlab, title = wrap(title, 78), subtitle = wrap(subtitle, 96),
         caption = wrap(caption, 110)) +
    theme_forest()
}

# ── Panel 1: national prevalence recovery ───────────────────────────────────
forest_national <- function(path = here::here("results", "tables",
                                              "national_estimates_all.csv")) {
  n <- readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(err = (pred_prev - obs_prev) * 100,
           lo  = (obs_lo  - obs_prev) * 100,
           hi  = (obs_hi  - obs_prev) * 100,
           inside = err >= lo & err <= hi,
           label = paste(country, "·", .lab_outcome(outcome))) %>%
    arrange(desc(abs(err)))
  forest_panel(
    n %>% transmute(label, est = err, lo, hi, flag = inside, bad = !inside),
    ref = 0, xlab = "Predicted − surveyed national prevalence (pp)",
    title = "The model reproduces national prevalence within survey error",
    subtitle = sprintf(
      "Bar = the survey's own 95%% CI. %d of %d predictions fall inside it; mean absolute error %.2f pp.",
      sum(n$inside), nrow(n), mean(abs(n$err))),
    caption = "Source: results/tables/national_estimates_all.csv")
}

# ── Panel 2: targeting value ────────────────────────────────────────────────
forest_targeting <- function(path = here::here("results", "tables", "corrected",
                                               "decision_value.csv")) {
  d <- readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(label = paste(country, "·", .lab_outcome(outcome)),
           flag = lift_vs_no_targeting > 1) %>%
    arrange(lift_vs_no_targeting)
  forest_panel(
    d %>% transmute(label, est = lift_vs_no_targeting,
                    lo = NA_real_, hi = NA_real_, flag, bad = FALSE),
    ref = 1, xlab = "Cases reached vs no targeting (×)",
    title = "Ranking districts beats spending the budget evenly",
    subtitle = sprintf(
      "Reach at the top 20%% of districts. Above 1.0× in %d of %d combinations; best %.2f×.",
      sum(d$flag), nrow(d), max(d$lift_vs_no_targeting)),
    caption = paste("Point estimates — this metric carries no interval.",
                    "Source: results/tables/corrected/decision_value.csv"))
}

# ── Panel 3: cross-country transport ────────────────────────────────────────
forest_transport <- function(path = here::here("results", "tables",
                                    "transportability_area_loco_metrics.csv")) {
  d <- readr::read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(pearson_r)) %>%
    mutate(label = paste0(.cap_country(held_out), " held out · ",
                          .lab_outcome(outcome)),
           flag = pearson_perm_p < 0.05) %>%
    arrange(pearson_r)
  forest_panel(
    d %>% transmute(label, est = pearson_r, lo = pearson_ci_lo,
                    hi = pearson_ci_hi, flag, bad = FALSE),
    ref = 0, xlab = "Correlation with the held-out country's district estimates",
    title = "District rankings survive crossing a border",
    subtitle = sprintf(
      "Trained on three countries, tested on a fourth with no local biomarker data. %d of %d significant against a permutation null.",
      sum(d$flag), nrow(d)),
    caption = "1,000-replicate bootstrap intervals. Source: results/tables/transportability_area_loco_metrics.csv")
}

# ── Panel 4: within-country district ranking ────────────────────────────────
forest_within <- function(path = here::here("results", "tables", "corrected",
                                            "area_recipe_within.csv")) {
  d <- readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(label = paste(country, "·", .lab_outcome(outcome)),
           flag = pearson_perm_p < 0.05) %>%
    arrange(pearson_r)
  forest_panel(
    d %>% transmute(label, est = pearson_r, lo = pearson_ci_lo,
                    hi = pearson_ci_hi, flag, bad = FALSE),
    ref = 0, xlab = "Correlation with surveyed district prevalence",
    title = "Where district-level signal is strongest",
    subtitle = sprintf(
      "Block cross-validated within country. %d of %d significant; strongest in iron and in the largest district samples.",
      sum(d$flag), nrow(d)),
    caption = "1,000-replicate bootstrap intervals. Source: results/tables/corrected/area_recipe_within.csv")
}

# ── Panel 5: the anchoring comparison (a table) ─────────────────────────────
table_anchor_arms <- function(path = here::here("results", "tables",
                                                "admin1_arms.csv")) {
  nm <- c("admin-2 fit + ADMIN-1 benchmark (hard)"   = "District model anchored to regional totals",
          "admin-2 fit + ADMIN-1 benchmark (shrunk)" = "District model, shrunken regional anchor",
          "ADMIN-1 fit -> admin-2 (pooled)"          = "Fit at region, extrapolated to district",
          "admin-2 fit + national benchmark"         = "District model anchored to national total",
          "admin-2 fit (LORO), unbenchmarked"        = "District model alone, no anchor")
  readr::read_csv(path, show_col_types = FALSE) %>%
    group_by(arm) %>%
    summarise(Cells = n(), `Mean r` = round(mean(pearson_r, na.rm = TRUE), 3),
              `MAE (pp)` = round(mean(mae_pp, na.rm = TRUE), 2),
              `Mean |bias| (pp)` = round(mean(abs(bias_pp), na.rm = TRUE), 2),
              .groups = "drop") %>%
    mutate(Approach = ifelse(arm %in% names(nm), nm[arm], arm)) %>%
    arrange(desc(`Mean r`)) %>%
    select(Approach, Cells, `Mean r`, `MAE (pp)`, `Mean |bias| (pp)`)
}

#' Write all four panels to disk (used by the qmd and the deck).
save_forest_figures <- function(dir = here::here("results", "figures", "forest"),
                                width = 9, height = NULL, dpi = 200) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  panels <- list(national = forest_national(), targeting = forest_targeting(),
                 transport = forest_transport(), within = forest_within())
  for (nm in names(panels)) {
    p <- panels[[nm]]
    h <- height %||% max(3.2, 0.24 * nrow(p$data) + 2.1)
    ggsave(file.path(dir, paste0("forest_", nm, ".png")), p,
           width = width, height = h, dpi = dpi, bg = "white")
  }
  cat(sprintf("[forest] wrote %d panels to %s\n", length(panels), dir))
  invisible(panels)
}
