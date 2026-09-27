# Smoke test, run from the repo root:
#   Rscript dashboard/data-raw/smoke_test.R
# Loads global.R the way the app does, walks every country x outcome through
# the map helpers at both levels, checks the district decomposition is exact,
# and constructs the app UI. The app has no other regression check before
# deployment.

owd <- setwd(here::here("dashboard"))
on.exit(setwd(owd), add = TRUE)
source("global.R")

fails <- character(0)
note <- function(ok, msg) { cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg)); if (!ok) fails <<- c(fails, msg) }

cat("Bundles\n")
note(nrow(idx_districts) > 0, sprintf("admin2_index: %d district rows, %d cells", nrow(idx_districts), length(idx_fits)))
note(all(is.finite(idx_national$national_prev)), "every cell has a national anchor")
note(all(is.finite(idx_national$anchor_shift)), "every cell's anchor converged")
note(!any(is_water(idx_districts$Admin2)), "no water polygons in the ranking table")
k <- paste(idx_districts$country, idx_districts$outcome, idx_districts$Admin1, idx_districts$Admin2)
note(!any(duplicated(k)), "district keys unique within each cell")
note(!is.null(CIV) && nrow(CIV$ranking) > 0, "Cote d'Ivoire bundle present")
note(length(EV) > 5, sprintf("protocol evidence: %d tables", length(EV) - 1))
note(!is.null(CAT) && nrow(CAT$variables) > 400, sprintf("catalogue: %d predictors", if (is.null(CAT)) 0 else nrow(CAT$variables)))
note(all(is.finite(c(Q$infill, Q$tr, Q$cs, Q$null_d, Q$cap_index, Q$ar_a1, Q$ceiling_prev))), "headline numbers computed")
note(length(UE) > 0 && !is.null(UE$cells) && nrow(UE$cells) > 3000,
     sprintf("stability ensembles: %d district rows", if (length(UE) && !is.null(UE$cells)) nrow(UE$cells) else 0))
note(!is.null(UE$pooled_weights) && nrow(UE$pooled_weights) > 1000, "pooled weight ranges present")
note(!is.null(EV$xv_pooled) && nrow(EV$xv_pooled) >= 5, "external validation tables present")
note(!is.null(EV$headroom) && !is.null(EV$rank_coverage), "headroom and rank-coverage tables present")
note(!is.null(EV$planner_summary) && nrow(EV$planner_summary) > 8, "survey-planner validation present")
note(!is.null(EV$importance_all) && nrow(EV$importance_all) > 20000, "searchable importance table present")
note("child_selenium" %in% idx_districts$outcome, "Malawi selenium ranked")
note(all(is.finite(c(Q$xv_level, Q$xv_off, Q$strong_mean, Q$stab_cov, Q$plan_spread))), "new headline numbers computed")
note("prev_cal_lo" %in% names(idx_districts) && sum(is.finite(idx_districts$prev_cal_lo)) > 3000,
     "calibrated prevalence bands on the districts (CP-01)")
note(is.finite(Q$cal_cov) && Q$cal_cov > 0.8 && Q$cal_cov < 1,
     sprintf("conformal LOO coverage computed (%.2f)", Q$cal_cov))
note(all(idx_districts$prev_cal_hi >= idx_districts$prev_cal_lo, na.rm = TRUE), "calibrated bands well-ordered")

cat("\nMap helpers\n")
n_ok <- 0L
for (ck in names(meta$countries)) {
  for (oc in outcomes_for(ck)) {
    res <- tryCatch({
      a2 <- get_country_admin2(ck, oc); stopifnot(!is.null(a2), nrow(a2) > 0, any(is.finite(a2$priority)))
      a1 <- get_country_admin1(ck, oc); stopifnot(!is.null(a1), nrow(a1) > 0)
      # the decomposition must reproduce the score exactly
      r <- sf::st_drop_geometry(a2); r <- r[is.finite(r$score_logit), ][1, ]
      dec <- decompose_district(ck, oc, r$Admin1, r$Admin2); stopifnot(!is.null(dec))
      fit <- idx_fits[[paste(ck, oc)]]
      mean_tr <- fit$intercept + sum(fit$beta * fit$mu)
      gap <- abs(attr(dec, "total") - (r$score_logit - mean_tr))
      stopifnot(gap < 1e-6)
      TRUE
    }, error = function(e) conditionMessage(e))
    if (isTRUE(res)) n_ok <- n_ok + 1L else fails <- c(fails, sprintf("%s / %s: %s", ck, oc, res))
  }
}
cat(sprintf("  %d country x outcome combinations through both levels and the decomposition\n", n_ok))

cat("\nApp UI\n")
res <- tryCatch({ eval(parse("app.R")[[2]]); TRUE },
                error = function(e) { fails <<- c(fails, sprintf("app.R: %s", conditionMessage(e))); FALSE })
note(isTRUE(res), "app.R UI constructs")

if (length(fails)) { cat("\nFAILURES:\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n"); quit(status = 1) }
cat("\nAll smoke checks passed.\n")
