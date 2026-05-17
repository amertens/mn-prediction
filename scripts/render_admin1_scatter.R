# =============================================================================
# scripts/render_admin1_scatter.R
#
# Combine per-country plot_admin1_scatter_<country>_<outcome> targets into
# multi-country aggregated scatter PNGs (one per outcome) matching the
# legacy admin1_scatter_{pop}_{nut}.png naming.
#
# Inputs:  _targets_full/objects/plot_admin1_scatter_<country>_<outcome>
# Outputs: results/figures/admin1_scatter_{children,women}_{vitA,iron}.png
#
# Naming: filenames use "children" not "child" to match the pre-existing
# (orphaned) committed figures the funder deck references.
# =============================================================================

library(targets)
library(here)
library(ggplot2)
library(patchwork)

TARGETS_STORE <- here("_targets_full")
fig_dir <- here("results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

countries <- c("gambia", "ghana", "sierraleone", "malawi")
country_labels <- c(gambia = "Gambia", ghana = "Ghana",
                    sierraleone = "Sierra Leone", malawi = "Malawi")

# Map of outcome tag (target name) -> filename pop_nut tag
outcomes <- list(
  child_vitA  = "children_vitA",
  child_iron  = "children_iron",
  women_vitA  = "women_vitA",
  women_iron  = "women_iron"
)

for (oc in names(outcomes)) {
  fn_tag <- outcomes[[oc]]
  cat(sprintf("\n=== %s -> admin1_scatter_%s.png ===\n", oc, fn_tag))

  panels <- list()
  for (cn in countries) {
    tname <- paste0("plot_admin1_scatter_", cn, "_", oc)
    p <- tryCatch(
      tar_read_raw(tname, store = TARGETS_STORE),
      error = function(e) NULL
    )
    if (is.null(p) || !inherits(p, "ggplot")) {
      cat(sprintf("  [skip] %s: %s\n", cn,
                  if (is.null(p)) "no target" else "not a ggplot"))
      next
    }
    # Re-title for combined layout: small country label, no redundant outcome
    p <- p +
      labs(title = country_labels[[cn]], subtitle = NULL) +
      theme(plot.title = element_text(size = 11, face = "bold"))
    panels[[cn]] <- p
  }

  if (length(panels) == 0) {
    cat(sprintf("  [skip] no panels for %s\n", oc))
    next
  }

  combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("Admin-1 predicted vs. observed prevalence — %s",
                      switch(oc,
                             child_vitA = "Vitamin A (children)",
                             child_iron = "Iron deficiency (children)",
                             women_vitA = "Vitamin A (women)",
                             women_iron = "Iron deficiency (women)")),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  outfile <- file.path(fig_dir, sprintf("admin1_scatter_%s.png", fn_tag))
  ggsave(outfile, combined,
         width = 10, height = 8, dpi = 150, bg = "white")
  cat(sprintf("  Saved %s (%d countries)\n", basename(outfile), length(panels)))
}

cat("\nDone.\n")
