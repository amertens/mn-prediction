# Smoke test — run from the repo root:
#   Rscript dashboard/data-raw/smoke_test.R
# Exercises the data bundle the way the app does: loads global.R, walks every
# country x outcome through the map helpers for every prediction layer, and
# constructs both app UIs. Restored from archive/ after the 2026-08 data refresh
# (the app has no other regression check before deployment).

owd <- setwd(here::here("dashboard"))
on.exit(setwd(owd), add = TRUE)
source("global.R")

stopifnot(nrow(admin2_pred) > 0, nrow(admin2_pop) > 0, length(admin2_bnds) == 4)

layers <- list(
  sl     = admin2_pred,
  area   = admin2_area_pred,
  fh     = admin2_fh_pred,
  bym2   = admin2_bym2_pred,
  recipe = admin2_recipe_pred
)
layers <- layers[!vapply(layers, is.null, logical(1))]
cat("Prediction layers present:", paste(names(layers), collapse = ", "), "\n\n")

fails <- character(0)

# ── Admin-2 key hygiene ───────────────────────────────────────────────────
# GADM ships inland water as Admin-2 polygons and repeats Admin-2 names, and
# both reach the map: the area/recipe layers are built from the polygon-ordered
# covariate frame, which keeps them on purpose. A water row here means Lake
# Malawi is being painted with a deficiency prevalence; a duplicated key means
# get_country_admin2()'s !duplicated() drops one polygon's prediction and paints
# the survivor's value on both. Assert both are gone.
source(here::here("R", "admin2_key_hygiene.R"))
cat("Admin-2 key hygiene:\n")
for (lname in names(layers)) {
  d <- layers[[lname]]
  w <- unique(d$Admin2[is_water_admin2(as.character(d$Admin2))])
  if (length(w))
    fails <- c(fails, sprintf("%s layer contains water polygons: %s",
                              lname, paste(w, collapse = ", ")))
  k <- paste(d$country, d$outcome, d$Admin2)
  if (any(duplicated(k))) {
    dk <- unique(d$Admin2[duplicated(k)])
    fails <- c(fails, sprintf("%s layer has %d duplicated country/outcome/Admin2 key(s): %s",
                              lname, sum(duplicated(k)), paste(utils::head(dk, 6), collapse = ", ")))
  }
  cat(sprintf("  %-7s water=%d duplicate-keys=%d\n", lname, length(w), sum(duplicated(k))))
}
for (nm in c("boundaries", "population")) {
  obj <- if (nm == "boundaries") do.call(rbind, lapply(admin2_bnds, function(b)
             data.frame(Admin2 = b$Admin2))) else admin2_pop
  w <- unique(obj$Admin2[is_water_admin2(as.character(obj$Admin2))])
  if (length(w))
    fails <- c(fails, sprintf("%s contains water polygons: %s", nm, paste(w, collapse = ", ")))
  cat(sprintf("  %-7s water=%d\n", nm, length(w)))
}
cat("\n")

for (lname in names(layers)) {
  pd <- layers[[lname]]
  n_ok <- 0L
  for (ctry in names(meta$countries)) {
    clab <- meta$countries[[ctry]]
    for (oc in unique(pd$outcome[pd$country == clab])) {
      res <- tryCatch({
        df <- get_country_admin2(ctry, oc, admin2_bnds, pd, admin2_pop)
        stopifnot(!is.null(df), nrow(df) > 0)
        natl <- national_aggregate(df)
        df1 <- get_country_admin1(ctry, oc, admin1_bnds, admin2_bnds, pd, admin2_pop)
        stopifnot(nrow(df1) > 0)
        list(ok = TRUE, natl = natl$pred_prev_natl)
      }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
      if (isTRUE(res$ok)) {
        n_ok <- n_ok + 1L
      } else {
        fails <- c(fails, sprintf("%s / %s / %s: %s", lname, ctry, oc, res$msg))
      }
    }
  }
  cat(sprintf("  %-7s %3d country x outcome combos OK\n", lname, n_ok))
}

cat("\nConstructing app UIs...\n")
for (entry in c("app.R", "app_public.R")) {
  res <- tryCatch({ eval(parse(entry)[[2]]); TRUE },
                  error = function(e) { fails <<- c(fails, sprintf("%s: %s", entry, conditionMessage(e))); FALSE })
  cat(sprintf("  %-14s %s\n", entry, if (isTRUE(res)) "sourced OK" else "FAILED"))
}

if (length(fails)) {
  cat("\nFAILURES:\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("\nAll smoke checks passed.\n")
