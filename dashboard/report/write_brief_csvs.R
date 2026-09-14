# Companion CSVs for the country briefs, straight from the tables the brief
# prints, so a number can be checked against its source rather than read off a
# chart. Sourced by render_reports.R; also runnable on its own from the repo root:
#   Rscript dashboard/report/write_brief_csvs.R            # all countries
#   Rscript dashboard/report/write_brief_csvs.R ghana      # one

write_brief_csvs <- function(ck, out_dir = here::here("dashboard", "report", "out")) {
  owd <- setwd(here::here("dashboard")); on.exit(setwd(owd), add = TRUE)
  suppressPackageStartupMessages(source("global.R"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  clab <- unname(meta$countries[[ck]])
  ocs <- outcomes_for(ck)
  dis <- do.call(rbind, lapply(ocs, function(oc) {
    d <- sf::st_drop_geometry(get_country_admin2(ck, oc))
    data.frame(country = clab, outcome = oc, outcome_label = unname(meta$outcome_labels[oc]),
               district = d$Admin2, region = d$Admin1,
               rank_worst = d$rank_worst, n_districts = d$n_districts, priority_score = d$priority,
               chance_worst_fifth = d$p_worst_fifth,
               planning_prevalence = d$prev_anchored, who_class = d$who_class,
               surveyed = d$surveyed, survey_prevalence = d$survey_prev,
               survey_lo = d$survey_lo, survey_hi = d$survey_hi,
               survey_respondents = d$n_resp, survey_clusters = d$n_clusters,
               population = d$population, people_affected = d$people_affected,
               stringsAsFactors = FALSE)
  }))
  nat <- idx_national[idx_national$country_key == ck, ]
  nat <- data.frame(country = clab, outcome = nat$outcome, outcome_label = unname(meta$outcome_labels[nat$outcome]),
                    national_prevalence_survey = nat$national_prev, districts_surveyed = nat$n_surveyed,
                    districts = nat$n_districts, stringsAsFactors = FALSE)
  utils::write.csv(dis, file.path(out_dir, sprintf("%s_districts.csv", ck)), row.names = FALSE)
  utils::write.csv(nat, file.path(out_dir, sprintf("%s_national.csv", ck)), row.names = FALSE)
  cat(sprintf("    csv: %d district rows, %d national rows\n", nrow(dis), nrow(nat)))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  for (ck in if (length(args)) args else c("gambia", "ghana", "sierraleone", "malawi")) {
    cat(sprintf("== %s ==\n", ck)); write_brief_csvs(ck)
  }
}
