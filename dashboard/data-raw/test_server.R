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
  session$setInputs(rank_from = "prev", set = "climate_soil", fraction = 0.05)
  d <- dat()
  note(nrow(d) > 0, sprintf("%d design rows", nrow(d)))
  note(renders(output$at_fraction), "design table renders")
})
testServer(mod_civ_server, args = list(id = "civ"), {
  session$setInputs(outcome = "child_iron", layer = "priority")
  d <- dat()
  note(nrow(d) == 33, sprintf("Cote d'Ivoire joins %d districts", nrow(d)))
  note(sum(is.finite(d$p_worst3rd)) == 33, "rank uncertainty joined for children's iron")
})

if (length(fails)) { cat("\nFAILURES:\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n"); quit(status = 1) }
cat("\nAll server checks passed.\n")
