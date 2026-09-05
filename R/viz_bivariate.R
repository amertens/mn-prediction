#' Bivariate choropleth maps for Admin-2 estimates
#'
#' Adapted from the design used by Tang et al. (Nature Food 2026,
#' doi:10.1038/s43016-026-01412-2, Figs 2-3), who map fortification-vehicle
#' reach against consumption quantity on a two-axis colour grid. Here the
#' same device answers a targeting question this project already asks: which
#' Admin-2 areas are both badly off AND something-else (uncertain, populous).
#'
#' Implemented without the `biscale` package -- it is a tercile x tercile
#' cross with a fixed nine-colour palette, and adding a dependency for that
#' is not worth it.
#'
#' @name viz_bivariate
NULL

# Joshua Stevens' classic 3x3 bivariate palette. Names are "<x>-<y>".
BIVARIATE_PAL_3 <- c(
  "1-1" = "#e8e8e8", "2-1" = "#ace4e4", "3-1" = "#5ac8c8",
  "1-2" = "#dfb0d6", "2-2" = "#a5add3", "3-2" = "#5698b9",
  "1-3" = "#be64ac", "2-3" = "#8c62aa", "3-3" = "#3b4994"
)

#' Assign each row to a bivariate class
#'
#' @param x,y Numeric vectors of equal length.
#' @param dim Integer, number of bins per axis (default 3).
#' @param group Optional grouping vector; when supplied, bins are computed
#'   *within* each group. Default behaviour for this project is to bin within
#'   country, because prevalence levels do not transport across countries
#'   (see `docs/findings/` on the LOCO level offset) -- so a cross-country
#'   tercile would just re-draw the country boundaries.
#' @return Character vector of "<x>-<y>" class labels, NA where either input
#'   is missing.
bivariate_class <- function(x, y, dim = 3, group = NULL) {
  stopifnot(length(x) == length(y), dim >= 2)
  if (is.null(group)) group <- rep("all", length(x))

  bin <- function(v, g) {
    out <- rep(NA_integer_, length(v))
    for (lev in unique(g)) {
      i <- which(g == lev & is.finite(v))
      if (!length(i)) next
      qs <- stats::quantile(v[i], probs = seq(0, 1, length.out = dim + 1),
                            na.rm = TRUE, type = 7)
      qs <- unique(qs)
      if (length(qs) < 2) { out[i] <- 1L; next }   # constant within group
      out[i] <- as.integer(cut(v[i], breaks = qs, include.lowest = TRUE,
                               labels = FALSE))
    }
    out
  }

  xb <- bin(x, group); yb <- bin(y, group)
  ifelse(is.na(xb) | is.na(yb), NA_character_, paste0(xb, "-", yb))
}

#' Build the square legend that explains a bivariate map
#'
#' @param xlab,ylab Axis labels (short -- they sit under a 3x3 grid).
#' @param dim Bins per axis; only 3 has a bundled palette.
#' @param pal Named colour vector, defaults to `BIVARIATE_PAL_3`.
bivariate_legend <- function(xlab, ylab, dim = 3, pal = BIVARIATE_PAL_3) {
  grid <- expand.grid(x = seq_len(dim), y = seq_len(dim))
  grid$fill <- pal[paste0(grid$x, "-", grid$y)]

  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = paste0(xlab, " \u2192"), y = paste0(ylab, " \u2192")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 7, margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(size = 7, angle = 90,
                                           margin = ggplot2::margin(r = 2))
    )
}

#' Bivariate Admin-2 choropleth for one country
#'
#' @param sf_admin2 An `sf` object of Admin-2 polygons with an `Admin2` column.
#' @param values data.frame with `Admin2` plus the two numeric columns.
#' @param x_col,y_col Column names in `values`.
#' @param xlab,ylab Legend axis labels.
#' @param title,subtitle Passed to `ggplot2::labs()`.
#' @param dim Bins per axis.
#' @return A `cowplot` object (map with the legend inset bottom-left).
bivariate_admin2_map <- function(sf_admin2, values, x_col, y_col,
                                 xlab = x_col, ylab = y_col,
                                 title = NULL, subtitle = NULL, dim = 3) {
  if (!requireNamespace("cowplot", quietly = TRUE))
    stop("bivariate_admin2_map() needs the 'cowplot' package")

  d <- merge(sf_admin2, values, by = "Admin2", all.x = TRUE)
  d$bi_class <- bivariate_class(d[[x_col]], d[[y_col]], dim = dim)
  d$fill <- unname(BIVARIATE_PAL_3[d$bi_class])
  d$fill[is.na(d$fill)] <- "grey92"          # no estimate for this area

  map <- ggplot2::ggplot(d) +
    ggplot2::geom_sf(ggplot2::aes(fill = fill), colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                   plot.subtitle = ggplot2::element_text(size = 8.5,
                                                         colour = "grey30"))

  cowplot::ggdraw() +
    cowplot::draw_plot(map) +
    cowplot::draw_plot(bivariate_legend(xlab, ylab, dim),
                       x = 0.00, y = 0.02, width = 0.30, height = 0.30)
}
