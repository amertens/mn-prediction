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
