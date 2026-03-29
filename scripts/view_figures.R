# =============================================================================
# scripts/view_figures.R
#
# Load and display saved figure objects from the results report.
# Run interactively in RStudio to resize figures in the Plots pane.
#
# Usage:
#   source("scripts/view_figures.R")
#   # Then call any figure by name, e.g.:
#   show_fig("forest_auc")
#   show_fig("forest_prauc")
#   show_fig("auc_heatmap")
#   show_fig("national_ci_forest")
#
#   # Or browse all figures interactively:
#   browse_figures()
#
# Prerequisites:
#   Render docs/results_report.qmd first to generate the figures RDS file.
# =============================================================================

library(here)
library(ggplot2)

# Load saved figures
fig_file <- here("results", "figures", "report_figures.rds")
if (!file.exists(fig_file)) {
  stop("Figure file not found. Render docs/results_report.qmd first:\n  ",
       "quarto::quarto_render('docs/results_report.qmd')")
}

figs <- readRDS(fig_file)
cat(sprintf("Loaded %d figures:\n", length(figs)))
for (nm in names(figs)) cat(sprintf("  - %s\n", nm))


#' Show a single figure by name
#'
#' @param name Character name of the figure (e.g., "forest_auc")
#' @param width Optional width override (inches)
#' @param height Optional height override (inches)
show_fig <- function(name, width = NULL, height = NULL) {
  if (!name %in% names(figs)) {
    cat(sprintf("Figure '%s' not found. Available: %s\n",
                name, paste(names(figs), collapse = ", ")))
    return(invisible(NULL))
  }
  p <- figs[[name]]
  if (!is.null(width) || !is.null(height)) {
    # For RStudio, set device size
    if (!is.null(width) && !is.null(height)) {
      grDevices::dev.new(width = width, height = height)
    }
  }
  print(p)
  invisible(p)
}


#' Browse all figures interactively
#'
#' Prints each figure one by one. Press Enter to advance.
browse_figures <- function() {
  for (i in seq_along(figs)) {
    nm <- names(figs)[i]
    cat(sprintf("\n=== [%d/%d] %s ===\n", i, length(figs), nm))
    print(figs[[nm]])
    if (i < length(figs)) {
      readline(prompt = "Press Enter for next figure (or Esc to stop)... ")
    }
  }
  cat("\nDone.\n")
}


#' Modify a figure's theme for presentation
#'
#' @param name Figure name
#' @param base_size Base font size (default 16 for presentations)
#' @param ... Additional theme() arguments
presentation_fig <- function(name, base_size = 16, ...) {
  p <- figs[[name]]
  if (is.null(p)) {
    cat(sprintf("Figure '%s' not found.\n", name))
    return(invisible(NULL))
  }
  p + theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(size = base_size + 2, face = "bold"),
          strip.text = element_text(size = base_size - 1, face = "bold"),
          ...)
}


#' Save a figure to a file
#'
#' @param name Figure name
#' @param filename Output filename (default: results/figures/{name}.png)
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
save_fig <- function(name, filename = NULL, width = 10, height = 7, dpi = 300) {
  p <- figs[[name]]
  if (is.null(p)) {
    cat(sprintf("Figure '%s' not found.\n", name))
    return(invisible(NULL))
  }
  if (is.null(filename)) {
    filename <- here("results", "figures", paste0(name, ".png"))
  }
  ggsave(filename, p, width = width, height = height, dpi = dpi)
  cat(sprintf("Saved: %s\n", filename))
}


cat("\nReady. Try: show_fig('forest_auc') or browse_figures()\n")
