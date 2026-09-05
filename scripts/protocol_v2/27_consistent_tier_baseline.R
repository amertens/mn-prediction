# =============================================================================
# scripts/protocol_v2/27_consistent_tier_baseline.R   [R6-02]
#
# THE SURVEY BASELINE ON THE CONSISTENT DISTRICT RUNG
#
# Run 6 (script 18) put all four countries on the same administrative rung --
# Malawi's 243 Traditional Authorities collapsed to its 27 districts -- and
# found the model numbers barely moved. It could not run the survey baseline
# (region_mean_jk) or burden capture, because at the district rung Malawi has
# no enclosing region in our data: GADM carries its 28 districts as level 1
# and omits the 3 official regions. This script supplies that lookup (verified
# against the district names in the data; it stops if any district is
# unmapped) and completes the comparison the NCE actually makes -- model vs
# the survey's own regional averages, on Spearman and burden captured -- on
# 146 genuinely comparable units.
#
# A second reason to care: with Malawi's real regions (9-13 districts each),
# its jackknifed regional mean is no longer a three-unit average, so the CF-01
# small-group artefact is smaller here than on the mixed rung.
#
#   Rscript scripts/protocol_v2/27_consistent_tier_baseline.R
# -> results/tables/protocol_v2/consistent_tier_baseline.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; REPS <- as.integer(Sys.getenv("R6_REPS", "10")); TOPFRAC <- 0.20
set.seed(20260903L)

MALAWI_REGION <- c(
  Chitipa = "Northern", Karonga = "Northern", Likoma = "Northern", Mzimba = "Northern",
  `Nkhata Bay` = "Northern", Rumphi = "Northern",
  Dedza = "Central", Dowa = "Central", Kasungu = "Central", Lilongwe = "Central", Mchinji = "Central",
  Nkhotakota = "Central", Ntcheu = "Central", Ntchisi = "Central", Salima = "Central",
  Balaka = "Southern", Blantyre = "Southern", Chikwawa = "Southern", Chiradzulu = "Southern",
  Machinga = "Southern", Mangochi = "Southern", Mulanje = "Southern", Mwanza = "Southern",
  Neno = "Southern", Nsanje = "Southern", Phalombe = "Southern", Thyolo = "Southern", Zomba = "Southern")

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds"); POP <- readRDS("dashboard/data/admin2_population.rds")
POP$country <- gsub(" ", "", POP$country)  # FIX 2026-09-04: the file spells "Sierra Leone" with a space; without this the join silently dropped the country
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2],
             stringsAsFactors = FALSE) }))

mw <- unique(TG$Admin1[TG$country == "Malawi"])
unm <- setdiff(mw, names(MALAWI_REGION))
if (length(unm)) stop("Malawi districts not in the region lookup: ", paste(unm, collapse = ", "))
cat("Malawi districts mapped:", length(mw), "of", length(mw), "| regions:", paste(names(table(MALAWI_REGION[mw])), table(MALAWI_REGION[mw]), collapse = ", "), "\n")

build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]; if (!nrow(t)) return(NULL)
  pc <- pop_for(on); pp <- POP[POP$country == cn, c("Admin2", pc)]; names(pp)[2] <- "pop"
  t <- left_join(t, pp, by = "Admin2"); t <- t[is.finite(t$pop) & t$pop > 0, ]; if (!nrow(t)) return(NULL)
  if (cn == "Malawi") {
    a <- t |> group_by(Admin1) |> summarise(y = stats::weighted.mean(y_prev, n_eff), w = sum(n_eff), pop = sum(pop), .groups = "drop")
    x <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1") |> inner_join(c1, by = "Admin1")
    unit <- m$Admin1; region <- unname(MALAWI_REGION[m$Admin1])
  } else {
    m <- t[, c("Admin1", "Admin2", "y_prev", "n_eff", "pop")]; names(m)[3:4] <- c("y", "w")
    m <- m |> inner_join(S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
      inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
    unit <- m$Admin2; region <- m$Admin1
  }
  keep <- is.finite(m$lon) & is.finite(m$y); m <- m[keep, ]; unit <- unit[keep]; region <- region[keep]
  if (nrow(m) < 12 || dplyr::n_distinct(region) < 2) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(n = nrow(m), y = m$y, w = m$w, pop = m$pop, region = region, X = Xr, D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = region, y_nat = m$y))
}
capture <- function(y, pop, s) { ok <- is.finite(y) & is.finite(pop) & is.finite(s); if (sum(ok) < 5) return(NA_real_)
  y <- y[ok]; pop <- pop[ok]; s <- s[ok]; b <- y * pop; k <- max(1L, round(TOPFRAC * length(y)))
  sel <- order(s, decreasing = TRUE)[seq_len(k)]; sum(b[sel]) / sum(b) }

rows <- list(); cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
ARMS <- c("null_train_mean", "region_mean_jk", "spatial", "domain_index", "spatial_plus_domain")
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) { cat("  build error", cn, on, conditionMessage(e), "\n"); NULL }); if (is.null(cl)) next
  ymod <- .v2_logit(cl$y)
  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    for (a in ARMS) {
      pred <- rep(NA_real_, cl$n)
      for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
        p <- tryCatch(ARMS_V2[[a]](tr, te, ymod, cl$X, cl$D, cl$aux), error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p }
      s <- score_v2(cl$y, pred, cl$w, scale = "prev")
      rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, rep = r, arm = a, n_units = cl$n,
        n_regions = dplyr::n_distinct(cl$region), spearman = s$spearman, capture = capture(cl$y, cl$pop, pred), stringsAsFactors = FALSE)
    }
  }
  cat("done", cn, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "consistent_tier_baseline.csv"), row.names = FALSE)
CELL <- R |> group_by(country, outcome, arm, n_units, n_regions) |> summarise(spearman = median(spearman, na.rm = TRUE), capture = mean(capture, na.rm = TRUE), .groups = "drop")
cat("\n===== R6-02: consistent district rung (146 units), in-fill, prevalence =====\n")
print(as.data.frame(CELL |> group_by(arm) |> summarise(cells = dplyr::n(), spearman = round(mean(spearman, na.rm = TRUE), 3),
  capture = round(mean(capture, na.rm = TRUE), 3), .groups = "drop") |> arrange(desc(capture))), row.names = FALSE)
for (m in c("spearman", "capture")) {
  W <- pivot_wider(CELL[, c("country", "outcome", "arm", m)], names_from = arm, values_from = all_of(m))
  for (a in c("domain_index", "spatial_plus_domain")) { d <- W[[a]] - W$region_mean_jk
    cat(sprintf("  %-8s %-20s vs region_mean_jk better in %2d of %2d | median %+.3f\n", m, a, sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE))) }
}
cat("\n--- Malawi on its real regions (Northern/Central/Southern) ---\n")
print(as.data.frame(CELL[CELL$country == "Malawi", c("outcome", "arm", "n_regions", "spearman", "capture")] |> mutate(across(c(spearman, capture), ~ round(.x, 3)))), row.names = FALSE)
cat("\nDONE\n")
