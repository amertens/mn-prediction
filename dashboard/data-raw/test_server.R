# Server-side checks, run from the repo root:
#   Rscript dashboard/data-raw/test_server.R
# smoke_test.R proves the UI constructs and the joins work; this exercises the
# reactives that decide what a user sees, module by module.

owd <- setwd(here::here("dashboard"))
on.exit(setwd(owd), add = TRUE)
source("global.R")
library(shiny)

fails <- character(0)
note <- function(ok, msg) { cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg)); if (!ok) fails <<- c(fails, msg) }
renders <- function(x) !is.null(x) && nzchar(as.character(x$html %||% x))

cat("Map explorer\n")
testServer(mod_map_explorer_server, args = list(id = "map"), {
  session$setInputs(country = "ghana", outcome = "child_vitA", admin_level = "admin2", layer = "priority")
  d <- map_data()
  note(!is.null(d) && nrow(d) == 260, sprintf("Ghana joins %d districts", if (is.null(d)) 0 else nrow(d)))
  note(sum(is.finite(d$priority)) == 260, "every Ghana district has a priority score")
  note(sum(d$surveyed, na.rm = TRUE) == 75, sprintf("%d surveyed districts carry a survey estimate", sum(d$surveyed, na.rm = TRUE)))
  note(renders(output$headline), "headline renders")
  h <- paste(as.character(output$headline$html), collapse = "")
  note(grepl("level skill: (none|weak|moderate|good)", h), "headline carries the level-skill badge (IS-01)")
  session$setInputs(layer = "prev_anchored")
  note(grepl("Level skill", output$caption), "caption names the level skill on the planning-prevalence layer")
  session$setInputs(layer = "priority", admin_level = "admin1")
  note(nrow(map_data()) == 16, sprintf("Ghana aggregates to %d regions", nrow(map_data())))
  session$setInputs(country = "malawi", outcome = "child_zinc", admin_level = "admin2")
  note(nrow(map_data()) > 200, "Malawi zinc draws")
  session$setInputs(country = "ghana", outcome = "child_vitA", layer = "rank_width", fade_unstable = TRUE)
  d <- map_data()
  note(sum(is.finite(d$rank_width)) > 200, sprintf("stability range joined for %d Ghana districts", sum(is.finite(d$rank_width))))
  session$setInputs(layer = "p_moderate_plus")
  note(sum(is.finite(d$p_moderate_plus)) > 200, "WHO exceedance joined")
  session$setInputs(country = "malawi", outcome = "child_selenium", layer = "priority")
  note(nrow(map_data()) > 200, "Malawi selenium draws")
})

cat("\nDistrict profiles\n")
testServer(mod_district_server, args = list(id = "district"), {
  session$setInputs(country = "ghana", district = "Northern|Tamale", outcome = "child_vitA")
  r <- rows()
  note(nrow(r) == 6, sprintf("Tamale has %d outcomes", nrow(r)))
  note(renders(output$summary), "summary renders")
  note(all(r$level_skill %in% c("none", "weak", "moderate", "good", "unknown")), "every outcome row carries a level-skill band")
})

cat("\nStart here\n")
testServer(mod_start_here_server, args = list(id = "start", go_to = function(...) invisible(NULL)), {
  note(renders(output$hero), "hero renders")
  note(renders(output$example), "worked example renders")
  note(renders(output$checklist), "checklist renders")
})

cat("\nCatalogue\n")
testServer(mod_catalogue_server, args = list(id = "catalogue"), {
  session$setInputs(by = "domain", pick = character(0), outcome = "child_iron", only_defined = FALSE, only_composite = FALSE, only_travel = FALSE)
  f <- filtered()
  note(nrow(f) == nrow(CAT$variables), sprintf("all %d predictors listed", nrow(f)))
  note(sum(is.finite(f$weight)) > 300, sprintf("%d carry a weight for child iron", sum(is.finite(f$weight))))
  session$setInputs(only_travel = TRUE)
  note(all(filtered()$climate_soil), "climate and soil filter works")
})

cat("\nTrust, targeting, survey design, Cote d'Ivoire\n")
testServer(mod_trust_server, args = list(id = "trust"), {
  session$setInputs(estimand = "country", target = "level", ceil_target = "prev")
  note(renders(output$test_text), "trust text renders")
})
testServer(mod_targeting_server, args = list(id = "targeting"), {
  session$setInputs(outcome = "child_vitA")
  note(renders(output$calibration_text), "calibration text renders")
})
testServer(mod_survey_design_server, args = list(id = "planning"), {
  session$setInputs(rank_from = "prev", set = "climate_soil", fraction = 0.05, size_country = "ghana",
                    country = "ghana", outcome = "child_vitA", k = 15, preset = "shrink")
  d <- dat()
  note(nrow(d) > 0, sprintf("%d design rows", nrow(d)))
  note(renders(output$at_fraction), "design table renders")
  pd <- plan_data()
  note(!is.null(pd) && sum(pd$selected) == 15, sprintf("planner proposes %d visits", if (is.null(pd)) 0 else sum(pd$selected)))
  note(all(is.finite(pd$plan_score)), "every district carries a planner score")
})
testServer(mod_civ_server, args = list(id = "civ"), {
  session$setInputs(outcome = "child_iron", layer = "priority")
  d <- dat()
  note(nrow(d) == 33, sprintf("Cote d'Ivoire joins %d districts", nrow(d)))
  note(sum(is.finite(d$p_worst3rd)) == 33, "rank uncertainty joined for children's iron")
  session$setInputs(outcome = "women_b12")
  note(sum(is.finite(dat()$p_worst3rd)) == 33, "rank uncertainty joined for women's B12 (CV-01 all outcomes)")
})

cat("\nImportance explorer, external check, roadmap\n")
testServer(mod_importance_server, args = list(id = "importance"), {
  session$setInputs(sc_scope = "pooled", sc_outcome = "child_iron", sc_target = "level", sc_domain = "", sc_country = "")
  d <- search_data()
  note(nrow(d) > 250, sprintf("importance explorer returns %d weights (pooled, child iron)", nrow(d)))
  note(sum(is.finite(d$lo)) > 200, "pooled weights carry resampling ranges")
  session$setInputs(sc_scope = "loco", sc_country = "Ghana")
  note(nrow(search_data()) > 250, "held-out-country fits searchable")
})
testServer(mod_external_server, args = list(id = "external"), {
  session$setInputs(target = "level")
  note(!is.null(output$cells), "external-validation chart renders")
})
testServer(mod_roadmap_server, args = list(id = "roadmap", go_to = function(...) invisible(NULL)), {
  note(!is.null(output$headroom), "headroom chart renders")
  note(!is.null(output$curve), "learning curve renders")
})

if (length(fails)) { cat("\nFAILURES:\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n"); quit(status = 1) }
cat("\nAll server checks passed.\n")
