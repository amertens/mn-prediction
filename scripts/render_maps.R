# =============================================================================
# scripts/render_maps.R
#
# Pre-render all Admin-1 and Admin-2 maps as publication-quality PNG files.
# Maps use a SINGLE shared legend per outcome across all countries.
# Output: results/figures/admin1_maps_*.png, results/figures/admin2_maps_*.png
# =============================================================================

library(targets)
library(here)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(viridis)
library(patchwork)

TARGETS_STORE <- here("_targets_full")

countries <- c("gambia", "ghana", "sierraleone", "malawi")
outcomes  <- c("child_vitA", "women_vitA", "child_iron", "women_iron",
               "women_folate", "women_b12", "child_zinc", "women_zinc")

country_labels <- c(gambia = "Gambia", ghana = "Ghana",
                     sierraleone = "Sierra Leone", malawi = "Malawi")
outcome_labels <- c(
  child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
  child_iron = "Iron deficiency (children)", women_iron = "Iron deficiency (women)",
  women_folate = "Folate deficiency (women)", women_b12 = "B12 deficiency (women)",
  child_zinc = "Zinc deficiency (children)", women_zinc = "Zinc deficiency (women)")

gadm_codes <- c(gambia = "GMB", ghana = "GHA", sierraleone = "SLE", malawi = "MWI")

load_target <- function(name) {
  tryCatch(targets::tar_read_raw(name, store = TARGETS_STORE),
           error = function(e) NULL)
}

fig_dir <- here("results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# ── Load and simplify GADM polygons ──────────────────────────────────────
cat("Loading GADM polygons...\n")
gadm_cache <- here("results", "figures", "gadm_cache.rds")
if (file.exists(gadm_cache)) {
  cached <- readRDS(gadm_cache)
  poly_a1 <- cached$poly_a1
  poly_a2 <- cached$poly_a2
} else {
  poly_a1 <- poly_a2 <- list()
  for (ctry in countries) {
    code <- gadm_codes[ctry]
    poly_a1[[ctry]] <- tryCatch({
      g <- geodata::gadm(country = code, level = 1, path = here("data", "gadm"))
      st_as_sf(g) |> select(NAME_1) |> rename(Admin1 = NAME_1) |>
        st_transform(4326) |> st_simplify(dTolerance = 0.005, preserveTopology = TRUE)
    }, error = function(e) NULL)
    poly_a2[[ctry]] <- tryCatch({
      g <- geodata::gadm(country = code, level = 2, path = here("data", "gadm"))
      st_as_sf(g) |> select(NAME_1, NAME_2) |>
        rename(Admin1 = NAME_1, Admin2 = NAME_2) |>
        st_transform(4326) |> st_simplify(dTolerance = 0.005, preserveTopology = TRUE)
    }, error = function(e) NULL)
  }
  saveRDS(list(poly_a1 = poly_a1, poly_a2 = poly_a2), gadm_cache)
}


# ── Load all Admin-1 data ────────────────────────────────────────────────
cat("Loading Admin-1 data...\n")
a1_list <- list()
for (ctry in countries) {
  for (oc in outcomes) {
    obj <- load_target(paste0("admin1_prev_", ctry, "_", oc))
    if (!is.null(obj) && is.data.frame(obj) && nrow(obj) > 0) {
      if (!"country" %in% names(obj)) obj$country <- country_labels[ctry]
      if (!"country_code" %in% names(obj)) obj$country_code <- ctry
      # Standardize column names
      if ("sl_prev" %in% names(obj) && !"pred_prevalence" %in% names(obj))
        obj <- obj |> rename(pred_prevalence = sl_prev)
      if ("obs_prev" %in% names(obj) && !"obs_prevalence" %in% names(obj))
        obj <- obj |> rename(obs_prevalence = obs_prev)
      a1_list[[paste0(ctry, "_", oc)]] <- obj
    }
  }
}
a1_all <- if (length(a1_list) > 0) bind_rows(a1_list) else NULL
cat(sprintf("  Loaded %d Admin-1 data frames, %d total rows\n",
            length(a1_list), if (!is.null(a1_all)) nrow(a1_all) else 0))


# ── Load all Admin-2 data ────────────────────────────────────────────────
cat("Loading Admin-2 data...\n")
a2_sl <- a2_svy <- list()
for (ctry in countries) {
  for (oc in outcomes) {
    key <- paste0(ctry, "_", oc)
    sl_d <- load_target(paste0("admin2_sl_", ctry, "_", oc))
    if (!is.null(sl_d) && is.data.frame(sl_d) && nrow(sl_d) > 0) a2_sl[[key]] <- sl_d
    svy_d <- load_target(paste0("svy_admin2_", ctry, "_", oc))
    if (!is.null(svy_d) && is.data.frame(svy_d) && nrow(svy_d) > 0) a2_svy[[key]] <- svy_d
  }
}


# ── Compute global limits per outcome ────────────────────────────────────
cat("Computing global colour scale limits...\n")
a1_global_lims <- list()
if (!is.null(a1_all)) {
  for (oc in outcomes) {
    d <- a1_all |> filter(outcome == oc)
    vals <- c(d$obs_prevalence, d$pred_prevalence)
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if (length(vals) > 0) a1_global_lims[[oc]] <- range(vals)
  }
}

a2_global_lims <- list()
for (oc in outcomes) {
  vals <- c()
  for (ctry in countries) {
    key <- paste0(ctry, "_", oc)
    if (!is.null(a2_sl[[key]])) vals <- c(vals, a2_sl[[key]]$sl_prev)
    if (!is.null(a2_svy[[key]])) vals <- c(vals, a2_svy[[key]]$svy_prev)
  }
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (length(vals) > 0) a2_global_lims[[oc]] <- range(vals)
}


# ── Helper: single map panel (no individual legend) ──────────────────────
map_panel <- function(poly, data, fill_col, title, limits, join_col = "Admin1") {
  if (is.null(poly) || is.null(data)) return(NULL)
  map_data <- poly |> left_join(data, by = join_col)
  ggplot(map_data) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "grey70", linewidth = 0.15) +
    scale_fill_viridis_c(
      name = "Prevalence",
      labels = percent_format(accuracy = 1),
      limits = limits,
      na.value = "grey90",
      option = "plasma",
      guide = "none"  # suppress per-panel legend
    ) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
}


# ── Shared legend as a standalone grob ───────────────────────────────────
make_legend_bar <- function(limits, title = "Prevalence") {
  # Create a dummy plot just to extract the legend
  dummy_df <- data.frame(x = 1:2, y = 1:2, fill = limits)
  p <- ggplot(dummy_df, aes(x, y, fill = fill)) +
    geom_point() +
    scale_fill_viridis_c(
      name = title,
      labels = percent_format(accuracy = 1),
      limits = limits,
      option = "plasma"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 15, barheight = 0.8,
      title.position = "top", title.hjust = 0.5
    )) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  cowplot::get_legend(p)
}


# ── RENDER ADMIN-1 MAPS ─────────────────────────────────────────────────
cat("\n=== Rendering Admin-1 maps ===\n")

for (oc in outcomes) {
  lims <- a1_global_lims[[oc]]
  if (is.null(lims)) next

  panels <- list()
  for (ctry in countries) {
    ctry_label <- country_labels[ctry]
    ctry_poly <- poly_a1[[ctry]]
    if (is.null(ctry_poly)) next

    key <- paste0(ctry, "_", oc)
    if (is.null(a1_all)) next
    d <- a1_all |> filter(country_code == ctry | country == ctry_label, outcome == oc)
    if (nrow(d) == 0) next

    has_obs  <- "obs_prevalence" %in% names(d) && any(!is.na(d$obs_prevalence))
    has_pred <- "pred_prevalence" %in% names(d) && any(!is.na(d$pred_prevalence))

    if (has_obs) {
      panels[[paste0(ctry_label, "\nObserved")]] <-
        map_panel(ctry_poly, d, "obs_prevalence",
                  paste0(ctry_label, " — Observed"), lims)
    }
    if (has_pred) {
      panels[[paste0(ctry_label, "\nPredicted")]] <-
        map_panel(ctry_poly, d, "pred_prevalence",
                  paste0(ctry_label, " — Predicted"), lims)
    }
  }

  if (length(panels) == 0) next

  # Compose: panels in a grid + single shared legend at bottom
  n_panels <- length(panels)
  ncols <- min(n_panels, 8)
  combined <- wrap_plots(panels, ncol = ncols) +
    plot_annotation(
      title = paste0("Admin-1 Prevalence: ", outcome_labels[oc]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  # Add shared legend
  legend <- make_legend_bar(lims)
  final <- cowplot::plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.06))

  outfile <- file.path(fig_dir, paste0("admin1_maps_", oc, ".png"))
  ggsave(outfile, final, width = min(n_panels * 2.5, 20), height = 5,
         dpi = 200, bg = "white")
  cat(sprintf("  Saved %s (%d panels)\n", basename(outfile), n_panels))
}


# ── RENDER ADMIN-2 MAPS ─────────────────────────────────────────────────
cat("\n=== Rendering Admin-2 maps ===\n")

for (oc in outcomes) {
  lims <- a2_global_lims[[oc]]
  if (is.null(lims)) next

  panels <- list()
  for (ctry in countries) {
    ctry_label <- country_labels[ctry]
    ctry_poly <- poly_a2[[ctry]]
    if (is.null(ctry_poly)) next

    key <- paste0(ctry, "_", oc)
    sl_data <- a2_sl[[key]]
    svy_data <- a2_svy[[key]]

    if (!is.null(svy_data)) {
      svy_sub <- svy_data |> select(Admin2, svy_prev)
      panels[[paste0(ctry_label, "\nSurvey")]] <-
        map_panel(ctry_poly, svy_sub, "svy_prev",
                  paste0(ctry_label, " — Survey"), lims, "Admin2")
    }
    if (!is.null(sl_data)) {
      panels[[paste0(ctry_label, "\nPredicted")]] <-
        map_panel(ctry_poly, sl_data, "sl_prev",
                  paste0(ctry_label, " — Predicted"), lims, "Admin2")
    }
  }

  if (length(panels) == 0) next

  n_panels <- length(panels)
  ncols <- min(n_panels, 8)
  combined <- wrap_plots(panels, ncol = ncols) +
    plot_annotation(
      title = paste0("Admin-2 Prevalence: ", outcome_labels[oc]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  legend <- make_legend_bar(lims)
  final <- cowplot::plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.06))

  outfile <- file.path(fig_dir, paste0("admin2_maps_", oc, ".png"))
  ggsave(outfile, final, width = min(n_panels * 2.5, 20), height = 5,
         dpi = 200, bg = "white")
  cat(sprintf("  Saved %s (%d panels)\n", basename(outfile), n_panels))
}


# ── RENDER CONFORMAL CI WIDTH MAPS (Admin-1) ─────────────────────────────
cat("\n=== Rendering Admin-1 CI width maps ===\n")

for (oc in outcomes) {
  panels <- list()
  ci_vals_all <- c()

  # First pass: collect all CI widths for global limits
  for (ctry in countries) {
    key <- paste0(ctry, "_", oc)
    b <- load_target(paste0("conformal_ci_", ctry, "_", oc))
    if (is.null(b) || is.null(b$admin1_ci)) next
    ci <- b$admin1_ci
    if ("ci_lo" %in% names(ci) && "ci_hi" %in% names(ci)) {
      ci$ci_width <- ci$ci_hi - ci$ci_lo
      ci_vals_all <- c(ci_vals_all, ci$ci_width[is.finite(ci$ci_width)])
    }
  }

  if (length(ci_vals_all) == 0) next
  ci_lims <- range(ci_vals_all, na.rm = TRUE)

  # Second pass: build panels
  for (ctry in countries) {
    ctry_label <- country_labels[ctry]
    ctry_poly <- poly_a1[[ctry]]
    if (is.null(ctry_poly)) next

    b <- load_target(paste0("conformal_ci_", ctry, "_", oc))
    if (is.null(b) || is.null(b$admin1_ci)) next
    ci <- b$admin1_ci
    if (!"ci_width" %in% names(ci)) ci$ci_width <- ci$ci_hi - ci$ci_lo

    map_data <- ctry_poly |> left_join(ci, by = "Admin1")
    p <- ggplot(map_data) +
      geom_sf(aes(fill = ci_width), colour = "grey70", linewidth = 0.15) +
      scale_fill_viridis_c(
        name = "CI Width",
        labels = percent_format(accuracy = 1),
        limits = ci_lims,
        option = "inferno", direction = -1,
        na.value = "grey90",
        guide = "none"
      ) +
      labs(title = ctry_label) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    panels[[ctry_label]] <- p
  }

  if (length(panels) == 0) next

  combined <- wrap_plots(panels, ncol = length(panels)) +
    plot_annotation(
      title = paste0("Conformal PI Width: ", outcome_labels[oc]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  # Legend for CI width
  dummy_df <- data.frame(x = 1:2, y = 1:2, fill = ci_lims)
  p_leg <- ggplot(dummy_df, aes(x, y, fill = fill)) + geom_point() +
    scale_fill_viridis_c(name = "CI Width", labels = percent_format(accuracy = 1),
                          limits = ci_lims, option = "inferno", direction = -1) +
    guides(fill = guide_colorbar(barwidth = 15, barheight = 0.8,
                                  title.position = "top", title.hjust = 0.5)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 10))
  legend <- cowplot::get_legend(p_leg)

  final <- cowplot::plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.08))

  outfile <- file.path(fig_dir, paste0("admin1_ci_width_", oc, ".png"))
  ggsave(outfile, final, width = length(panels) * 3.5, height = 4.5,
         dpi = 200, bg = "white")
  cat(sprintf("  Saved %s (%d panels)\n", basename(outfile), length(panels)))
}

cat("\nDone. All maps saved to results/figures/\n")
