# =============================================================================
# dashboard/data-raw/07_build_country_briefs.R
#
# One-page HTML brief per country, pre-rendered from the same bundles the app
# reads (a shinyapps dyno should serve files, not render reports). Each brief:
# the thesis line, the headline outcome's priority map, the worst-fifth list
# with how firmly each district is placed, every outcome's national figure and
# level skill, and the standing caveats. Written to dashboard/briefs/ and
# served by the Map explorer's download button.
#
#   Rscript dashboard/data-raw/07_build_country_briefs.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf); library(htmltools); library(here)})
setwd(here::here("dashboard"))
source("global.R")
dir.create("briefs", showWarnings = FALSE)

b64_png <- function(plot, w = 7.2, h = 5.2) {
  f <- tempfile(fileext = ".png")
  ggsave(f, plot, width = w, height = h, dpi = 130, bg = "white")
  on.exit(unlink(f))
  sprintf("data:image/png;base64,%s", base64enc::base64encode(f))
}

for (ck in names(meta$countries)) {
  cn <- meta$countries[[ck]]
  ocs <- outcomes_for(ck)
  headline_oc <- if ("child_vitA" %in% ocs) "child_vitA" else ocs[[1]]
  df <- get_country_admin2(ck, headline_oc)
  nat <- idx_national[idx_national$country_key == ck, ]
  k <- ceiling(g1(nat$n_districts[nat$outcome == headline_oc]) / 5)
  worst <- sf::st_drop_geometry(df) |> filter(is.finite(rank_worst)) |> arrange(rank_worst) |> head(k)

  map_plot <- ggplot(df) +
    geom_sf(aes(fill = priority), colour = "white", linewidth = 0.15) +
    geom_sf(data = df[isTRUE_v <- !is.na(df$surveyed) & df$surveyed, ], fill = NA, colour = "#1a1a1a", linewidth = 0.35) +
    scale_fill_gradientn(colours = c("#FFFFCC", "#FD8D3C", "#B10026"), limits = c(0, 100),
                         name = "Priority score\n(100 = ranked worst)", na.value = "#d9d9d9") +
    labs(caption = sprintf("%s, %s. Dark outline: districts the %s survey reached.", cn, outcome_short[[headline_oc]], meta$survey_years[[ck]])) +
    theme_void(base_size = 12) + theme(legend.position = "right", plot.caption = element_text(size = 8, colour = "grey35"))

  oc_rows <- lapply(ocs, function(oc) {
    nr <- nat[nat$outcome == oc, ]
    tags$tr(tags$td(outcome_short[[oc]]),
            tags$td(fmt_pct(g1(nr$national_prev))),
            tags$td(sprintf("%d of %d", g1(nr$n_surveyed), g1(nr$n_districts))),
            tags$td(skill_word[[level_skill(g1(nr$rho_train))$band]]))
  })
  worst_rows <- lapply(seq_len(nrow(worst)), function(i) {
    r <- worst[i, ]
    tags$tr(tags$td(sprintf("%d", r$rank_worst)), tags$td(r$Admin2), tags$td(r$Admin1),
            tags$td(if (is.finite(r$rank_lo)) sprintf("%d to %d", round(r$rank_lo), round(r$rank_hi)) else "—"),
            tags$td(if (is.finite(r$p_worst_fifth)) fmt_pct(r$p_worst_fifth, 0) else "—"))
  })

  page <- tags$html(tags$head(tags$meta(charset = "utf-8"),
    tags$title(sprintf("%s - micronutrient district brief", cn)),
    tags$style(HTML(paste("body{font-family:'Segoe UI',Arial,sans-serif;max-width:820px;margin:24px auto;color:#222;line-height:1.45;padding:0 16px;}",
      "h1{font-size:1.5em;color:#0F7B8A;margin-bottom:2px;} h2{font-size:1.05em;margin:18px 0 6px;}",
      ".thesis{color:#555;font-style:italic;margin-top:0;}",
      "table{border-collapse:collapse;width:100%;font-size:0.9em;} td,th{border-bottom:1px solid #e3e3e3;padding:4px 8px;text-align:left;}",
      ".note{font-size:0.8em;color:#666;background:#f6f8f9;border-left:3px solid #0F7B8A;padding:8px 12px;margin-top:14px;}",
      "img{max-width:100%;}")))),
    tags$body(
      h1(sprintf("%s: which districts to look at first", cn)),
      p(class = "thesis", "Modelling extends a survey's reach. It does not replace the survey."),
      p(sprintf(paste("Every district of %s is ranked by how likely it is to have high micronutrient deficiency, using",
                      "public data only (satellite imagery, climate, soil, crops, food prices and public household surveys),",
                      "and checked against the %s national biomarker survey. The ranking is the main result; any percentage",
                      "is a planning estimate based on the survey's national figure."), cn, meta$survey_years[[ck]])),
      h2(sprintf("The map: %s", outcome_short[[headline_oc]])),
      tags$img(src = b64_png(map_plot), alt = "district priority map"),
      h2(sprintf("The worst-ranked fifth (%d districts), and how firmly each is placed", nrow(worst))),
      tags$table(tags$thead(tags$tr(tags$th("Rank"), tags$th("District"), tags$th("Region"),
                                    tags$th("Rank range when re-estimated"), tags$th("Placed in worst fifth (share of test runs)"))),
                 tags$tbody(worst_rows)),
      h2("Every outcome the survey measured"),
      tags$table(tags$thead(tags$tr(tags$th("Outcome"), tags$th("National prevalence (survey)"),
                                    tags$th("Districts surveyed"), tags$th("Reliability of district percentages"))),
                 tags$tbody(oc_rows)),
      div(class = "note",
        p(strong("How much to trust this. "),
          sprintf(paste("Inside a surveyed country the model's ranking matches the survey's at %s, compared with %s for the",
                        "survey's own regional averages and at most %s for a random ranking. The rank range shows how far a",
                        "district's rank moves when the model is re-estimated; it is not a confidence interval. The ranking",
                        "shows which districts are likely worse off, not how severe the problem is, and a percentage needs",
                        "the national survey figure."),
                  f2(Q$infill), f2(Q$infill_jk), f2(Q$null_d))),
        p(strong("For the next survey. "), "The dashboard's Plan a survey page suggests which districts to visit for a",
          " chosen aim, and shows that a national sample 5% the size of a full survey, combined with the ranking, gives",
          " rough district estimates. The national prevalence itself still needs a probability sample.")),
      p(style = "font-size:0.75em;color:#999;margin-top:16px;",
        sprintf("Generated %s from the live dashboard's data build (%s). Interactive version: the Micronutrient Burden dashboard.",
                format(Sys.time(), "%Y-%m-%d"), PROTOCOL_LABEL))))

  out <- sprintf("briefs/brief_%s.html", ck)
  writeLines(as.character(page), out)
  cat(sprintf("  wrote %s (%s, %d outcomes)\n", out, cn, length(ocs)))
}
cat("done\n")
