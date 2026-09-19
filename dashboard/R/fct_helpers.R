# =============================================================================
# Helper functions used across modules
# =============================================================================

fmt_pct <- function(x, digits = 1) ifelse(is.finite(x), sprintf(paste0("%.", digits, "f%%"), x * 100), "—")
fmt_count <- function(x) ifelse(is.finite(x), formatC(round(x), format = "d", big.mark = ","), "—")
fmt_num <- function(x, d = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "—")

methods_note <- function(...) {
  htmltools::div(
    class = "methods-note",
    style = paste("font-size: 0.85em; color: #555; background: #f8f9fa; border-left: 3px solid #2c7bb6;",
                  "padding: 0.75em 1em; margin-top: 0.5em; margin-bottom: 1em; line-height: 1.5;"),
    htmltools::tags$strong("How to read this: "), ...)
}
empty_state <- function(msg) div(class = "alert alert-secondary", style = "font-size:0.9em;", msg)

# ── Level skill band (IS-01) ─────────────────────────────────────────────────
# The planning prevalence is the ranking mapped onto the outcome scale with a
# spread of rho x the survey's between-district spread, rho being the index's
# nested out-of-sample correlation in that country x outcome (idx_national$
# rho_train). Where rho is ~0 the map is flat at the national figure: the
# model has no level information there and only the ranking should be read.
# Bands are for reading, not thresholds of the analysis.
level_skill <- function(rho) {
  band <- ifelse(!is.finite(rho), "unknown",
          ifelse(rho < 0.10, "none", ifelse(rho < 0.30, "weak", ifelse(rho < 0.50, "moderate", "good"))))
  col <- c(none = "#b2182b", weak = "#e08214", moderate = "#4393c3", good = "#1a9850", unknown = "#999999")[band]
  list(band = band, colour = unname(col), rho = rho)
}
level_skill_badge <- function(rho) {
  ls <- level_skill(rho)
  htmltools::tags$span(style = sprintf("display:inline-block; padding:1px 7px; border-radius:9px; color:white; background:%s; font-size:0.85em;", ls$colour),
                       sprintf("level skill: %s (ρ = %s)", ls$band, fmt_num(rho, 2)))
}
level_skill_text <- function(rho) {
  switch(level_skill(rho)$band,
    none = "The model has no out-of-sample level skill for this outcome here: every district's planning prevalence sits at the national figure. Read the ranking, not the percentages.",
    weak = "Weak level skill: the planning prevalence varies little between districts and its spread is mostly the national figure. Read the ranking first.",
    moderate = "Moderate level skill: the planning prevalence carries about a third to a half of the true between-district spread.",
    good = "Good level skill: the planning prevalence carries half or more of the true between-district spread.",
    "Level skill not available for this cell.")
}

# GADM ships inland water as Admin-2 polygons (Lake Malawi, eight features).
# The prediction table drops them; the boundaries must too, or the map paints
# grey lakes and counts them as "no data" districts.
WATER_PATTERN <- "(^|\\b)(lake|lac|water ?body|waterbody|reservoir|lagoon)(\\b|$)"
is_water <- function(x) grepl(WATER_PATTERN, x, ignore.case = TRUE)

.key <- function(a1, a2) paste(trimws(as.character(a1)), trimws(as.character(a2)), sep = "|")

#' The joined district sf for one country and outcome: boundaries plus the
#' deployment ranking, the survey's estimate where one exists, and population.
get_country_admin2 <- function(ck, oc) {
  bnd <- admin2_bnds[[ck]]
  if (is.null(bnd)) return(NULL)
  bnd <- bnd[!is_water(bnd$Admin2), ]
  d <- idx_districts[idx_districts$country_key == ck & idx_districts$outcome == oc, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[!duplicated(.key(d$Admin1, d$Admin2)), ]
  i <- match(.key(bnd$Admin1, bnd$Admin2), .key(d$Admin1, d$Admin2))
  keep <- setdiff(names(d), c("Admin1", "Admin2", "country", "country_key", "outcome"))
  for (k in keep) bnd[[k]] <- d[[k]][i]
  bnd$country_key <- ck; bnd$outcome <- oc
  bnd$who_class[is.na(bnd$who_class)] <- "No data"
  bnd$worst_fifth <- is.finite(bnd$rank_worst) & bnd$rank_worst <= ceiling(bnd$n_districts / 5)
  bnd
}

#' Regions: population-weighted aggregates of the district table.
get_country_admin1 <- function(ck, oc) {
  a2 <- get_country_admin2(ck, oc)
  bnd1 <- admin1_bnds[[ck]]
  if (is.null(a2) || is.null(bnd1)) return(NULL)
  df <- sf::st_drop_geometry(a2)
  df <- df[!is.na(df$Admin1), ]
  pw <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) return(NA_real_); sum(x[ok] * w[ok]) / sum(w[ok]) }
  agg <- do.call(rbind, lapply(split(df, df$Admin1), function(g) {
    w <- g$population
    data.frame(Admin1 = g$Admin1[1],
               priority = pw(g$priority, w),
               prev_anchored = pw(g$prev_anchored, w),
               survey_prev = pw(g$survey_prev, w),
               p_worst_fifth = pw(g$p_worst_fifth, w),
               n_districts = nrow(g), n_surveyed = sum(g$surveyed, na.rm = TRUE),
               n_worst_fifth = sum(g$worst_fifth, na.rm = TRUE),
               population = sum(w, na.rm = TRUE),
               people_affected = sum(g$people_affected, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  bnd1$Admin1 <- trimws(bnd1$Admin1)
  out <- merge(bnd1, agg, by = "Admin1", all.x = TRUE, sort = FALSE)
  out$who_class <- "No data"
  th <- meta$who_thresholds[[oc]]
  if (!is.null(th)) {
    cl <- cut(out$prev_anchored, breaks = c(-Inf, sort(as.numeric(th)), Inf),
              labels = c("Low", "Mild", "Moderate", "Severe"), right = FALSE)
    out$who_class <- ifelse(is.na(cl), "No data", as.character(cl))
  }
  out$Admin2 <- NA_character_; out$surveyed <- out$n_surveyed > 0
  out$rank_worst <- rank(-out$priority, ties.method = "min", na.last = "keep")
  out
}

#' Predictor label for display: the plain-language name where one exists, else a cleaned code.
pred_label <- function(cols) {
  if (is.null(CAT) || !"label" %in% names(CAT$variables)) return(cols)
  v <- CAT$variables
  l <- v$label[match(cols, v$column)]
  ifelse(is.na(l), cols, l)
}
#' Labels made unique for a categorical axis: two columns can share a definition.
unique_labels <- function(labels, cols) {
  dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  ifelse(dup, paste0(labels, " (", cols, ")"), labels)
}
pred_source <- function(cols) {
  if (is.null(CAT)) return(rep(NA_character_, length(cols)))
  CAT$variables$source_label[match(cols, CAT$variables$column)]
}
pred_domain <- function(cols) {
  if (is.null(CAT)) return(rep(NA_character_, length(cols)))
  CAT$variables$domain[match(cols, CAT$variables$column)]
}

#' Exact decomposition of a district's score into predictor contributions.
#' The index is linear in the rank-normal predictors, so contribution_j =
#' beta_j (x_j - mean_j over training districts), on the logit scale of the
#' prediction. They sum to the district's deviation from the training mean.
decompose_district <- function(ck, oc, admin1, admin2) {
  fit <- idx_fits[[paste(ck, oc)]]; xr <- idx_xr[[ck]]
  if (is.null(fit) || is.null(xr)) return(NULL)
  row <- xr[match(.key(admin1, admin2), rownames(xr)), , drop = TRUE]
  if (all(is.na(row))) return(NULL)
  cols <- intersect(names(fit$beta), names(row))
  contrib <- fit$beta[cols] * (row[cols] - fit$mu[cols])
  out <- data.frame(column = cols, value = as.numeric(row[cols]), contribution = as.numeric(contrib),
                    label = pred_label(cols), source = pred_source(cols), domain = pred_domain(cols),
                    stringsAsFactors = FALSE)
  out <- out[order(-abs(out$contribution)), ]
  attr(out, "total") <- sum(contrib, na.rm = TRUE)
  out
}

#' Mean over cells with a 95% interval across cells.
cell_ci <- function(d, val, by) {
  d |> group_by(across(all_of(by))) |>
    summarise(n = sum(is.finite(.data[[val]])), est = mean(.data[[val]], na.rm = TRUE),
              se = sd(.data[[val]], na.rm = TRUE) / sqrt(pmax(n, 1)),
              lo = est - 1.96 * se, hi = est + 1.96 * se,
              pos = sum(.data[[val]] > 0, na.rm = TRUE), .groups = "drop")
}

#' Horizontal dot-and-interval plot (plotly).
forest_plotly <- function(d, y, null = NA, xlab = "Ranking accuracy (mean over cells, 95% interval)",
                          colour = PROXY_COL, height = NULL) {
  d <- as.data.frame(d)
  d <- d[is.finite(d$est), , drop = FALSE]
  d[[y]] <- factor(d[[y]], levels = rev(unique(d[[y]])))
  p <- plot_ly(d, height = height) |>
    add_segments(x = ~lo, xend = ~hi, y = ~.data[[y]], yend = ~.data[[y]],
                 line = list(color = "#9a9a9a", width = 2), showlegend = FALSE, hoverinfo = "none") |>
    add_markers(x = ~est, y = ~.data[[y]], marker = list(color = colour, size = 11),
                text = ~sprintf("%s<br>%.2f (%.2f to %.2f), %d cells", .data[[y]], est, lo, hi, n),
                hoverinfo = "text", showlegend = FALSE)
  shapes <- if (is.finite(null)) list(list(type = "line", x0 = null, x1 = null, y0 = 0, y1 = 1, yref = "paper",
                                           line = list(color = "#666", dash = "dash"))) else list()
  p |> layout(xaxis = list(title = xlab, zeroline = FALSE), yaxis = list(title = ""), shapes = shapes,
              margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
}
