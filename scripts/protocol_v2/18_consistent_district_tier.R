# =============================================================================
# scripts/protocol_v2/18_consistent_district_tier.R   [RUN 6]
#
# THE SAME ADMINISTRATIVE RUNG IN ALL FOUR COUNTRIES
#
# THE PROBLEM THIS FIXES
# ----------------------
# Every analysis so far has taken "Admin-1" or "Admin-2" uniformly across
# countries. Those labels do not denote the same thing:
#
#   country        Admin-1 in our data     Admin-2 in our data
#   Ghana          16 regions              260 districts
#   Gambia          6 divisions             37 districts
#   Sierra Leone    4 provinces             14 districts
#   Malawi         28 DISTRICTS           243 TRADITIONAL AUTHORITIES
#
# Malawi sits one rung finer at both tiers. GADM confirms it: gadm41_MWI_1 is
# the 28 districts, not the 3 official regions, which GADM does not carry. So
# the "Admin-2" set that carries most of our district-level results is a MIXTURE
# of districts (Ghana, Gambia, Sierra Leone) and sub-district Traditional
# Authorities (Malawi) -- and Malawi contributes 87 of its 206 units, the
# largest single block, at the finest and noisiest resolution in the set
# (median 9 biomarker measurements per unit against Gambia's 23).
#
# This run puts all four countries on the DISTRICT rung:
#
#   Ghana, Gambia, Sierra Leone   Admin-2 as-is
#   Malawi                        Admin-2 (TAs) aggregated up to Admin-1
#
# giving 146 genuinely comparable units. If the district-level numbers move,
# the Admin-1 vs Admin-2 gap reported earlier (0.50 vs 0.28) was partly a
# mislabelling artefact rather than pure aggregation gain, and the NCE's
# district-level claims need restating on this set.
#
# WHAT IS NOT RUN, AND WHY
# ------------------------
# No region-based arm (region_mean_jk) and no region estimand. At the district
# rung Malawi has no enclosing region in our data -- its 3 official regions are
# absent from GADM -- so a regional survey baseline cannot be computed for it
# without inventing one. Rather than hardcode a district-to-region lookup for
# one country and not the others, this run covers the two estimands that need
# no enclosing region and that carry the NCE's claims:
#   in-fill    5-fold over districts, replicated
#   transport  leave-one-country-out
#
#   Rscript scripts/protocol_v2/18_consistent_district_tier.R
# -> results/tables/protocol_v2/consistent_district_tier.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
REPS   <- as.integer(Sys.getenv("R6_REPS", "10"))
set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

# Malawi is the only country whose district rung is our "Admin1".
# Membership is tested with %in%, NOT [[ ]]: `c(Malawi = TRUE)[["Ghana"]]`
# throws "subscript out of bounds", and inside the tryCatch that wraps the cell
# builder it silently became "no cell for Ghana" -- the first run of this script
# reported Malawi alone and looked like a data problem rather than a typo.
DISTRICT_IS_ADMIN1 <- c("Malawi")

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

#' Build one country-outcome cell at the DISTRICT rung.
#' For Malawi that means collapsing Traditional Authorities into districts,
#' precision-weighting the outcome by effective n and averaging predictors.
build_district <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]
  if (!all(c(ycol, wcol) %in% names(t))) return(NULL)
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]
  if (!nrow(t)) return(NULL)
  collapse <- cn %in% DISTRICT_IS_ADMIN1

  if (collapse) {
    a <- t |> group_by(Admin1) |>
      summarise(y = stats::weighted.mean(.data[[ycol]], .data[[wcol]]),
                w = sum(.data[[wcol]]), .groups = "drop") |> filter(is.finite(y))
    x <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |>
      summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |>
      summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE),
                .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1") |> inner_join(c1, by = "Admin1")
    unit <- m$Admin1
  } else {
    a <- t[, c("Admin1", "Admin2", ycol, wcol)]
    names(a)[3:4] <- c("y", "w")
    m <- a |>
      inner_join(S[S$country == cn, c("Admin1", "Admin2", PREDS)],
                 by = c("Admin1", "Admin2")) |>
      inner_join(CENT[CENT$country == cn, c("Admin1","Admin2","lon","lat")],
                 by = c("Admin1", "Admin2"))
    unit <- m$Admin2
  }
  keep <- is.finite(m$lon) & is.finite(m$y)
  unit <- unit[keep]                 # subset BEFORE m, or the lengths diverge
  m <- m[keep, ]
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), unit = unit,
       y_nat = m$y, y_mod = if (target == "prev") .v2_logit(m$y) else m$y,
       X = Xr, w = m$w, lon = m$lon, lat = m$lat, collapsed = collapse)
}

rows <- list()
for (target in c("level", "prev")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      z <- tryCatch(build_district(cn, on, target), error = function(e) NULL)
      if (!is.null(z)) cl[[cn]] <- z
    }
    if (!length(cl)) next

    # ---- A. in-fill, replicated, within each country -----------------------
    for (cn in names(cl)) {
      z <- cl[[cn]]
      aux <- list(lon = z$lon, lat = z$lat, Admin1 = z$unit, y_nat = z$y_nat)
      D <- domain_representation_v2(z$X, domain_of)
      for (a in c("null_train_mean", "spatial", "domain_index",
                  "domain_enet", "spatial_plus_domain")) {
        for (r in seq_len(REPS)) {
          folds <- make_folds_v2("kfold_district", z$n, k = 5, rep_id = r)
          pred <- rep(NA_real_, z$n)
          for (f in unique(folds)) {
            te <- which(folds == f); tr <- which(folds != f)
            if (length(tr) < 12) next
            p <- tryCatch(ARMS_V2[[a]](tr, te, z$y_mod, z$X, D, aux),
                          error = function(e) rep(NA_real_, length(te)))
            if (length(p) == length(te)) pred[te] <- p
          }
          s <- score_v2(z$y_nat, pred, z$w,
                        scale = if (target == "prev") "prev" else "level")
          rows[[length(rows) + 1L]] <- data.frame(
            country = cn, outcome = on, target = target, estimand = "infill",
            arm = a, rep = r, n_units = z$n, collapsed = z$collapsed,
            spearman = s$spearman, wmae = s$wmae, stringsAsFactors = FALSE)
        }
      }
    }

    # ---- B. transport, leave-one-country-out --------------------------------
    if (length(cl) < 3) { cat("tier done", target, on, "(no LOCO)\n"); next }
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next
    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat))
    wv   <- unlist(lapply(cl, function(z) z$w))
    aux  <- list(lon = unlist(lapply(cl, function(z) z$lon)),
                 lat = unlist(lapply(cl, function(z) z$lat)),
                 Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$unit))),
                 y_nat = Y)
    folds <- as.integer(factor(ctry))
    for (a in c("domain_index", "domain_enet")) {
      pred <- rep(NA_real_, length(Y))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 20) next
        Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        p <- tryCatch(ARMS_V2[[a]](tr, te, Y, Xm, Dm, aux),
                      error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p
      }
      for (cn in unique(ctry)) {
        k <- which(ctry == cn)
        s <- score_v2(ynat[k], pred[k], wv[k],
                      scale = if (target == "prev") "prev" else "level")
        rows[[length(rows) + 1L]] <- data.frame(
          country = cn, outcome = on, target = target, estimand = "country",
          arm = a, rep = 1L, n_units = length(k),
          collapsed = cn %in% DISTRICT_IS_ADMIN1,
          spearman = s$spearman, wmae = NA_real_, stringsAsFactors = FALSE)
      }
    }
    cat("tier done", target, on, "\n")
  }
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "consistent_district_tier.csv"), row.names = FALSE)

CELL <- R |> group_by(country, outcome, target, estimand, arm, n_units) |>
  summarise(spearman = median(spearman, na.rm = TRUE),
            wmae = median(wmae, na.rm = TRUE), .groups = "drop")

cat("\n===== RUN 6: ALL FOUR COUNTRIES ON THE DISTRICT RUNG =====\n")
cat("units per country:\n")
print(as.data.frame(CELL |> distinct(country, n_units) |>
                      group_by(country) |>
                      summarise(units = max(n_units), .groups = "drop")),
      row.names = FALSE)

for (es in c("infill", "country")) {
  d <- CELL[CELL$estimand == es, ]
  if (!nrow(d)) next
  cat("\n--- ", es, " ---\n", sep = "")
  s <- d |> group_by(target, arm) |>
    summarise(cells = dplyr::n(),
              mean_rho = round(mean(spearman, na.rm = TRUE), 3),
              median_rho = round(median(spearman, na.rm = TRUE), 3),
              cells_positive = sum(spearman > 0, na.rm = TRUE),
              .groups = "drop") |> arrange(target, desc(mean_rho))
  print(as.data.frame(s), row.names = FALSE)
}

cat("\n--- Malawi only: does collapsing TAs to districts change it? ---\n")
mw <- CELL[CELL$country == "Malawi", ]
print(as.data.frame(mw |> group_by(estimand, target, arm) |>
        summarise(rho = round(median(spearman, na.rm = TRUE), 3),
                  .groups = "drop")), row.names = FALSE)
cat("\n(compare against benchmarks_v2_summary.csv, where Malawi's 87 units are\n")
cat(" Traditional Authorities rather than its 27 districts.)\n")
cat("\nDONE\n")
