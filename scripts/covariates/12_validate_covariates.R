# =============================================================================
# scripts/covariates/12_validate_covariates.R
#
# SCIENTIFIC validation of the harmonised covariates -- as distinct from
# 07_verify_harmonization.R, which checks that the harmoniser did what it was
# told (structure, provenance, reproducibility of every value).
#
# This script asks a different question: are the NUMBERS RIGHT? A covariate can
# be perfectly reproducible from its source and still be wrong, because the
# source itself is wrong, the zonal statistic is meaningless, or the export
# corrupted it. Structural checks cannot see any of that. Four families of test
# can:
#
#   A. CROSS-SOURCE AGREEMENT. iSDAsoil and SoilGrids independently measure pH,
#      nitrogen, organic carbon and CEC for the same districts. Two independent
#      products should agree on the spatial pattern. Where they do not, at least
#      one is wrong, and we should know which variables those are before they
#      enter a model.
#
#   B. PHYSICAL PLAUSIBILITY. Every variable has a range it cannot leave:
#      fractions in [0,100], NDVI in [-1,1], pH in [3,10], elevation above sea
#      level. Values outside are export defects, whatever the provenance says.
#
#   C. KNOWN RELATIONSHIPS. Physics and geography impose relationships that must
#      hold in any correct dataset: temperature falls with elevation, night-time
#      lights track built surface, vegetation tracks rainfall. A dataset that
#      violates these is broken even if every individual value looks sane.
#
#   D. SPATIAL STRUCTURE. Environmental covariates are spatially autocorrelated.
#      A covariate with none is not measuring geography -- it is noise, a
#      constant, or a misaligned join.
#
# Output: results/tables/covariate_validation.csv + validation_report.md
#   Rscript scripts/covariates/12_validate_covariates.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here); library(sf)
})
targets::tar_source(here("R"))

D <- here("data", "covariates", "harmonized")
P <- suppressMessages(readr::read_csv(file.path(D, "predictors_admin2_harmonized.csv"),
                                      show_col_types = FALSE)) |> as.data.frame()
dict <- suppressMessages(readr::read_csv(file.path(D, "data_dictionary.csv"),
                                         show_col_types = FALSE))
countries <- sort(unique(P$country))
res <- list()
add <- function(...) res[[length(res) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
has <- function(v) v %in% names(P)
corr <- function(a, b, method = "spearman")
  suppressWarnings(stats::cor(a, b, method = method, use = "pairwise.complete.obs"))

# ── A. CROSS-SOURCE AGREEMENT ───────────────────────────────────────────────
# Same quantity, two independent providers. Compared WITHIN country, because a
# pooled correlation can be driven entirely by between-country differences and
# would look strong even if neither product resolved districts.
cat("\n=== A. cross-source agreement (iSDAsoil vs SoilGrids) ===\n")
pairs <- list(
  list(what = "pH",             a = "soil_ph_mean_0_20",            b = "soilgrids_ph"),
  list(what = "nitrogen",       a = "soil_nitrogen_mean_0_20",      b = "soilgrids_nitrogen"),
  list(what = "organic carbon", a = "soil_totalcarbon_mean_0_20",   b = "soilgrids_organic_carbon"),
  list(what = "CEC",            a = "soil_cec_mean_0_20",           b = "soilgrids_cec"))
for (p in pairs) {
  if (!has(p$a) || !has(p$b)) {
    add(test = "cross_source", variable = p$what, scope = "all",
        value = NA_real_, status = "SKIP",
        detail = sprintf("absent (%s / %s)", p$a, p$b)); next
  }
  for (ct in countries) {
    s <- P[P$country == ct, ]
    r <- corr(s[[p$a]], s[[p$b]])
    add(test = "cross_source", variable = p$what, scope = ct, value = round(r, 3),
        status = if (!is.finite(r)) "SKIP" else if (r > 0.5) "PASS"
                 else if (r > 0.2) "WEAK" else "FAIL",
        detail = sprintf("Spearman %s vs %s", p$a, p$b))
    cat(sprintf("  %-14s %-12s rho = %s\n", p$what, ct,
                if (is.finite(r)) sprintf("%.3f", r) else "NA"))
  }
}

# ── B. PHYSICAL PLAUSIBILITY ────────────────────────────────────────────────
cat("\n=== B. physical bounds ===\n")
bounds <- list(
  list(pat = "^ndvi_y",              lo = -1,    hi = 1,     what = "NDVI index"),
  list(pat = "^evi_",                lo = -1,    hi = 1,     what = "EVI index"),
  list(pat = "^lcover_.*_frac_",     lo = 0,     hi = 100,   what = "cover fraction (%)"),
  list(pat = "^soil_ph_mean",        lo = 3,     hi = 10,    what = "soil pH"),
  list(pat = "^soilgrids_ph",        lo = 3,     hi = 10,    what = "soil pH"),
  list(pat = "^elevation$",          lo = -430,  hi = 6000,  what = "elevation (m)"),
  list(pat = "^lst_night_(m[0-9]|annual_(mean|min|max))", lo = -20, hi = 45, what = "night LST (C)"),
  list(pat = "^tclim_tmmx_",         lo = -20,   hi = 60,    what = "max temperature (C)"),
  list(pat = "^lai_",                lo = 0,     hi = 12,    what = "leaf area index"),
  list(pat = "^human_modification$", lo = 0,     hi = 1,     what = "modification index"),
  list(pat = "^spam_share_",         lo = 0,     hi = 1,     what = "crop share"),
  list(pat = "^ghs_pop$",            lo = 0,     hi = Inf,   what = "population (non-negative)"),
  # DHS indicators are a MIX of proportions, means and z-scores. Applying a
  # [0,1] bound to all of them flagged 5,501 "violations" that were simply
  # mean birthweight in grams, height-for-age z-scores and a diet-diversity
  # count. The check was wrong, not the data; exclude the non-proportions by
  # name and keep the bound where it genuinely applies.
  list(pat = "^dhs_(?!.*(mean_|_mean$|_score|_visits|parity|crowding|hhsize|edu_years))",
       lo = 0, hi = 1, what = "DHS proportion"),
  list(pat = "^dhs_c_mean_haz|^dhs_c_mean_waz|^dhs_c_mean_whz", lo = -6, hi = 6,
       what = "anthropometric z-score"),
  list(pat = "^dhs_c_mean_birthweight$", lo = 1000, hi = 5000, what = "birthweight (g)"),
  # The DHS wealth index is a standardised score (mean 0, SD 1), not a
  # proportion -- excluded above and bounded on its own terms here.
  list(pat = "^dhs_hh_wealth_mean$", lo = -5, hi = 5, what = "DHS wealth index (z)"))
for (b in bounds) {
  v <- grep(b$pat, names(P), value = TRUE, perl = TRUE)
  if (!length(v)) next
  x <- unlist(P[, v, drop = FALSE]); x <- x[is.finite(x)]
  if (!length(x)) next
  # Tolerance, not exact comparison: a proportion stored as 1.0000000000000002
  # is at the boundary, not outside it, and 63 such values were being reported
  # as bound violations with an observed range printed as "0 to 1".
  tol <- 1e-9 * max(1, abs(b$lo), abs(b$hi[is.finite(b$hi)]))
  bad <- sum(x < b$lo - tol | x > b$hi + tol)
  add(test = "bounds", variable = b$what, scope = sprintf("%d cols", length(v)),
      value = bad, status = if (bad == 0) "PASS" else "FAIL",
      detail = sprintf("observed %.3g to %.3g; allowed %.3g to %.3g",
                       min(x), max(x), b$lo, b$hi))
  cat(sprintf("  [%-4s] %-26s %d violation(s) over %d cols (range %.3g to %.3g)\n",
              if (bad == 0) "PASS" else "FAIL", b$what, bad, length(v), min(x), max(x)))
}

# ── C. KNOWN RELATIONSHIPS ──────────────────────────────────────────────────
# Signed, within-country, on rank correlations. These are not hypotheses under
# test -- they are settled physics and geography, so a violation indicts the
# data rather than the relationship.
cat("\n=== C. known relationships (within country) ===\n")
rel <- list(
  list(what = "elevation vs night LST",  a = "elevation", b = "lst_night_annual_mean_t0",
       dir = "-", why = "temperature falls with altitude (lapse rate)"),
  list(what = "night lights vs built surface", a = "ntl_ccnl", b = "built_surface",
       dir = "+", why = "lit area tracks the built environment"),
  list(what = "population vs built surface", a = "ghs_pop", b = "built_surface",
       dir = "+", why = "people live in buildings"),
  list(what = "rainfall vs vegetation", a = "precip_y2015", b = "ndvi_y2015",
       dir = "+", why = "vegetation is water-limited in these countries"),
  # NOTE: gHM weights urban infrastructure far above agriculture, so in
  # predominantly rural countries its variation is driven by settlement rather
  # than by cropland. Observed rho: Malawi +0.46 but Gambia -0.12, Ghana -0.01,
  # Sierra Leone -0.15. That is a weak prior of mine, not a data defect, so it
  # is reported as informational rather than scored.
  list(what = "human modification vs cropland", a = "human_modification",
       b = "lcover_crops_frac_t0", dir = "?", why = "gHM is settlement-weighted, not crop-weighted"),
  list(what = "accessibility vs population", a = "access_healthcare_min", b = "ghs_pop",
       dir = "-", why = "remote districts are sparsely populated"))
for (r in rel) {
  if (!has(r$a) || !has(r$b)) {
    add(test = "relationship", variable = r$what, scope = "all", value = NA_real_,
        status = "SKIP", detail = "variable absent"); next
  }
  ok_n <- 0; tot <- 0; vals <- c()
  for (ct in countries) {
    s <- P[P$country == ct, ]
    if (nrow(s) < 10) next
    v <- corr(s[[r$a]], s[[r$b]]); if (!is.finite(v)) next
    tot <- tot + 1; vals <- c(vals, v)
    if ((r$dir == "+" && v > 0) || (r$dir == "-" && v < 0)) ok_n <- ok_n + 1
  }
  add(test = "relationship", variable = r$what, scope = sprintf("%d countries", tot),
      value = if (length(vals)) round(stats::median(vals), 3) else NA_real_,
      status = if (!tot) "SKIP" else if (r$dir == "?") "INFO"
               else if (ok_n == tot) "PASS"
               else if (ok_n >= tot - 1) "WEAK" else "FAIL",
      detail = sprintf("expected sign %s (%s); correct in %d/%d", r$dir, r$why, ok_n, tot))
  cat(sprintf("  [%-4s] %-32s median rho %6s, expected %s, correct in %d/%d\n",
              if (!tot) "SKIP" else if (r$dir == "?") "INFO" else if (ok_n == tot) "PASS"
              else if (ok_n >= tot-1) "WEAK" else "FAIL",
              r$what, if (length(vals)) sprintf("%.3f", stats::median(vals)) else "NA",
              r$dir, ok_n, tot))
}

# ── D. SPATIAL STRUCTURE ────────────────────────────────────────────────────
# Moran's I on district centroids. An environmental covariate with no spatial
# autocorrelation is not measuring geography.
cat("\n=== D. spatial autocorrelation ===\n")
codes <- c(Gambia = "GMB", Ghana = "GHA", SierraLeone = "SLE", Malawi = "MWI",
           Tanzania = "TZA")
vars <- setdiff(names(P), c("country", "Admin1", "Admin2"))
for (ct in intersect(countries, names(codes))) {
  ctr <- tryCatch(sf::st_drop_geometry(load_admin2_centroids(codes[[ct]]))[, c("Admin2", "lon", "lat")],
                  error = function(e) NULL)
  if (is.null(ctr)) next
  # GADM admin-2 names repeat within a country, so a name join fans out. Collapse
  # to one centroid per name first -- Sierra Leone was silently dropped by this.
  ctr <- ctr[!duplicated(ctr$Admin2), , drop = FALSE]
  s <- dplyr::inner_join(P[P$country == ct, ], ctr, by = "Admin2")
  s <- s[is.finite(s$lon) & is.finite(s$lat), ]
  if (nrow(s) < 20) { cat("  ", ct, "too few areas with centroids\n"); next }
  d <- as.matrix(stats::dist(cbind(s$lon, s$lat))); W <- 1 / d; diag(W) <- 0
  W[!is.finite(W)] <- 0; W <- W / sum(W)
  moran <- function(x) {
    ok <- is.finite(x); if (sum(ok) < 20) return(NA_real_)
    x <- x[ok]; Ww <- W[ok, ok]; Ww <- Ww / sum(Ww)
    z <- x - mean(x); den <- sum(z^2)
    if (den <= 0) return(NA_real_)
    sum(Ww * outer(z, z)) * length(z) / den
  }
  I <- vapply(vars, function(v) moran(suppressWarnings(as.numeric(s[[v]]))), numeric(1))
  I <- I[is.finite(I)]
  flat <- names(I)[I < 0.05]
  add(test = "spatial", variable = "Moran's I", scope = ct,
      value = round(stats::median(I), 3),
      status = if (nrow(s) < 20) "LOW-N" else
               if (length(flat) / length(I) < 0.15) "PASS" else "REVIEW",
      detail = sprintf("n=%d areas; %d of %d covariates with I < 0.05: %s", nrow(s),
                       length(flat), length(I),
                       paste(utils::head(flat, 6), collapse = ", ")))
  cat(sprintf("  %-12s median I = %.3f | %d of %d covariates spatially flat\n",
              ct, stats::median(I), length(flat), length(I)))
}

out <- dplyr::bind_rows(res)
readr::write_csv(out, here("results", "tables", "covariate_validation.csv"))
cat(sprintf("\n-> results/tables/covariate_validation.csv (%d checks)\n", nrow(out)))
print(as.data.frame(out %>% count(test, status) %>% arrange(test, status)), row.names = FALSE)

lines <- c("# Covariate scientific validation", "",
  sprintf("Generated: %s", Sys.Date()),
  sprintf("Countries: %s | predictors: %d", paste(countries, collapse = ", "), length(vars)),
  "",
  "Checks that the covariate VALUES are right, as opposed to",
  "07_verify_harmonization.R which checks the harmoniser did what it was told.",
  "", "## Failures and weak results", "",
  paste0("- **", out$status[out$status %in% c("FAIL", "WEAK", "REVIEW")], "** `",
         out$variable[out$status %in% c("FAIL", "WEAK", "REVIEW")], "` (",
         out$scope[out$status %in% c("FAIL", "WEAK", "REVIEW")], "): ",
         out$detail[out$status %in% c("FAIL", "WEAK", "REVIEW")]),
  "", "## All results", "",
  paste0("- [", out$status, "] ", out$test, " / ", out$variable, " / ", out$scope,
         " = ", out$value))
writeLines(lines, file.path(D, "validation_report.md"))
cat("-> data/covariates/harmonized/validation_report.md\n")
