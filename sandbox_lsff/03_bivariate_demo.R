# Demo of the Tang et al. bivariate-map device, re-pointed at this project's
# targeting question. Two axes worth crossing:
#   (a) predicted prevalence x interval width  -> where to send the next survey
#   (b) predicted prevalence x child population -> where the burden actually is
# Layout follows their Fig 1: small multiples, ONE shared legend.
suppressPackageStartupMessages({library(dplyr); library(sf); library(ggplot2); library(cowplot)})
root <- ".."
source(file.path(root, "R/viz_bivariate.R"))

OUTCOME <- "child_iron"          # the slide-21 worked example
bnd  <- readRDS(file.path(root, "dashboard/data/admin2_boundaries.rds"))
fh   <- readRDS(file.path(root, "dashboard/data/admin2_fh_predictions.rds"))
pop  <- readRDS(file.path(root, "dashboard/data/admin2_population.rds"))
key  <- c(Gambia = "gambia", Ghana = "ghana", `Sierra Leone` = "sierraleone",
          Malawi = "malawi")

# bare map panel (no inset legend -- one shared legend for the whole figure)
panel <- function(cn, v, x_col, y_col) {
  d <- merge(bnd[[key[[cn]]]], v, by = "Admin2", all.x = TRUE)
  d$bi <- bivariate_class(d[[x_col]], d[[y_col]])
  # NB: missing must not be a light grey -- the palette's own low-low cell is
  # #e8e8e8. Draw no-estimate areas as white with a visible outline instead.
  d$fill <- unname(BIVARIATE_PAL_3[d$bi])
  miss <- is.na(d$fill)
  d$fill[miss] <- "white"
  d$edge <- ifelse(miss, "grey55", "white")
  n_miss <- sum(miss)
  ggplot(d) +
    geom_sf(aes(fill = fill, colour = edge), linewidth = 0.12) +
    scale_fill_identity() + scale_colour_identity() +
    labs(title = cn,
         subtitle = if (n_miss) paste0(n_miss, " areas without an estimate (white)")
                    else "all areas estimated") +
    theme_void(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, colour = "grey45", hjust = 0.5),
          plot.margin   = margin(4, 4, 4, 4))
}

build <- function(x_col, y_col, xlab, ylab, title, note, outfile) {
  ps <- lapply(names(key), function(cn) {
    v <- fh %>% filter(country == cn, outcome == OUTCOME) %>%
      left_join(pop %>% filter(country == cn) %>% select(Admin2, pop_child),
                by = "Admin2") %>%
      select(Admin2, pred_prev, ci_width, pop_child)
    panel(cn, v, x_col, y_col)
  })
  grid <- plot_grid(plotlist = ps, nrow = 2, ncol = 2)
  leg  <- bivariate_legend(xlab, ylab)
  head <- ggdraw() +
    draw_label(title, fontface = "bold", size = 13, x = 0.01, hjust = 0, y = 0.68) +
    draw_label(note, size = 9, colour = "grey35", x = 0.01, hjust = 0, y = 0.26)
  right <- plot_grid(NULL, leg, NULL, ncol = 1, rel_heights = c(1, 1.1, 1))
  body <- plot_grid(grid, right, nrow = 1, rel_widths = c(1, 0.24))
  ggsave(file.path(root, "results/figures", outfile),
         plot_grid(head, body, ncol = 1, rel_heights = c(0.10, 1)),
         width = 11, height = 10, dpi = 200, bg = "white")
}

dir.create(file.path(root, "results/figures"), showWarnings = FALSE, recursive = TRUE)
build("pred_prev", "ci_width", "predicted prevalence", "interval width",
      paste0("Where the next survey should go \u2014 ", OUTCOME, ", Fay-Herriot"),
      paste("Dark blue = high predicted prevalence AND a wide interval: bad, and we do not know how bad.",
            "Terciles within country."),
      "bivariate_prevalence_x_uncertainty.png")
build("pred_prev", "pop_child", "predicted prevalence", "child population",
      paste0("Where the burden is \u2014 ", OUTCOME),
      paste("Dark blue = high rate AND many children. A high rate in a small district is not the same target.",
            "Terciles within country."),
      "bivariate_prevalence_x_burden.png")
cat("wrote both figures\n")
