# =============================================================================
# sandbox/20_upstream_data_fixes.R
#
# Applies the four upstream data fixes recommended after sandbox/19:
#
#   FIX 1 — B12 Ghana outlier cleaning
#     gw_B12_pmol_L has mean 12,139 with sd 50,012 — implausibly high
#     (normal range is ~150-700 pmol/L). Cap at biologically plausible
#     upper bound (1500 pmol/L) and report how many values were
#     winsorised.
#
#   FIX 2 — Folate biomarker harmonisation
#     Malawi's `rbcf` is red-blood-cell folate, biologically distinct
#     from serum folate (Ghana, SL). Drop Malawi from the folate LOCO
#     pool; report leave-one-out on the remaining 2 countries
#     (Ghana <-> Sierra Leone) as a sensitivity.
#
#   FIX 3 — Vitamin A Malawi shift correction
#     Malawi's `rbp` has mean 1.16 vs 0.83-0.88 for the other three
#     countries (30% systematic shift, likely an inflammation-
#     adjustment difference). Shift Malawi to match the pooled mean
#     of the other three.
#
#   FIX 4 — Sierra Leone polygons
#     Investigate admin-3 (chiefdom) availability. If chief polygons
#     exist in GADM, retest LOCO with finer resolution; otherwise
#     document why SL is structurally underpowered (n_test = 14).
#
# Each fix: rerun spatial+top-8-SoilGrids LOCO with the fix applied,
# compare to unfixed baseline.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "config.R"))

countries_lower <- c("gambia","ghana","sierraleone","malawi")
country_labels  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")
cc_list <- get_country_configs()

# Helper: load merged_<country> data
load_merged <- function() {
  setNames(lapply(countries_lower, function(c)
    tryCatch(targets::tar_read_raw(paste0("merged_", c), store = STORE),
              error = function(e) NULL)),
    names(country_labels))
}

# Helper: per-admin-2 mean of a biomarker column
admin2_mean <- function(merged, col) {
  if (is.null(merged) || !col %in% colnames(merged)) return(NULL)
  vals <- suppressWarnings(as.numeric(unclass(merged[[col]])))
  ok <- is.finite(vals)
  if (sum(ok) < 3) return(NULL)
  df <- data.frame(Admin2 = merged$Admin2[ok], value = vals[ok])
  df |> dplyr::group_by(Admin2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
}

# Helper: build pooled per outcome
build_pooled <- function(per_country_bio) {
  per_country <- list()
  for (n in names(per_country_bio)) {
    bio <- per_country_bio[[n]]
    if (is.null(bio) || nrow(bio) < 3) next
    gee <- tryCatch(targets::tar_read_raw(
      paste0("gee_admin2_", country_labels[[n]]), store = STORE),
      error = function(e) NULL)
    if (is.null(gee)) next
    df <- dplyr::inner_join(bio, gee, by = "Admin2")
    df$country <- n
    per_country[[n]] <- df
  }
  if (length(per_country) < 2) return(NULL)
  gee_vars <- Reduce(intersect, lapply(per_country, function(x)
    grep("^gee_", colnames(x), value = TRUE)))
  pooled <- dplyr::bind_rows(lapply(per_country, function(x)
    x[, c("country","Admin2","mean_bio","n_obs", gee_vars)]))
  svy_list <- setNames(lapply(country_labels[names(per_country)], function(c)
    tryCatch(targets::tar_read_raw(
      paste0("svy_admin2_", c, "_child_vitA"), store = STORE),
      error = function(e) NULL)),
    names(per_country))
  ok <- !vapply(svy_list, is.null, logical(1))
  pooled <- add_admin2_centroids(pooled, cc_list, svy_list[ok])
  list(pooled_data = pooled, gee_vars = gee_vars,
        country_names = names(per_country))
}

# Helper: LOCO with spatial + top-8 soil features predicting mean_bio
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS
spatial_soil_predict <- function(train, test) {
  if (!all(c("lon","lat") %in% colnames(train))) return(NULL)
  avail <- intersect(TOP8, colnames(train))
  med <- vapply(as.data.frame(train[, avail, drop = FALSE]),
                 median, numeric(1), na.rm = TRUE)
  imp_tr <- as.data.frame(train[, avail, drop = FALSE])
  imp_te <- as.data.frame(test [, avail, drop = FALSE])
  for (v in avail) {
    imp_tr[[v]][!is.finite(imp_tr[[v]])] <- med[v]
    imp_te[[v]][!is.finite(imp_te[[v]])] <- med[v]
  }
  fit_df <- data.frame(Y = train$mean_bio, lon = train$lon, lat = train$lat,
                        imp_tr, check.names = FALSE)
  test_df <- data.frame(   lon = test$lon,  lat = test$lat,
                            imp_te, check.names = FALSE)
  rhs <- paste(c("s(lon, lat, k = 20, bs = 'tp')", avail), collapse = " + ")
  form <- stats::as.formula(paste("Y ~", rhs))
  fit <- tryCatch(mgcv::gam(form, data = fit_df, method = "REML",
                              weights = pmax(train$n_obs, 1)),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  as.numeric(mgcv::predict.gam(fit, newdata = test_df))
}

run_loco <- function(pooled, label = "") {
  pd <- pooled$pooled_data
  countries <- pooled$country_names
  res <- list()
  for (held in countries) {
    train <- pd[pd$country != held, , drop = FALSE]
    test  <- pd[pd$country == held, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    p <- spatial_soil_predict(train, test)
    if (is.null(p) || length(p) != nrow(test)) next
    ok <- is.finite(p) & is.finite(test$mean_bio)
    if (sum(ok) < 3) next
    rp <- suppressWarnings(cor(test$mean_bio[ok], p[ok]))
    rs <- suppressWarnings(cor(test$mean_bio[ok], p[ok], method = "spearman"))
    res[[held]] <- data.frame(label = label, held_out = held,
                                n_test = sum(ok),
                                pearson_r = round(rp, 3),
                                spearman_r = round(rs, 3))
  }
  if (length(res) == 0) data.frame()
  else dplyr::bind_rows(res)
}

merged <- load_merged()

# ───────────────────────────────────────────────────────────────────────────
cat("\n##### FIX 1: B12 Ghana outlier cleaning #####\n")
# Inspect Ghana B12 distribution
m_gh <- merged$Ghana
b12_raw <- suppressWarnings(as.numeric(unclass(m_gh$gw_B12_pmol_L)))
b12_raw <- b12_raw[is.finite(b12_raw)]
cat(sprintf("Ghana gw_B12_pmol_L: n=%d, min=%.0f, p25=%.0f, median=%.0f, p75=%.0f, p95=%.0f, p99=%.0f, max=%.0f, mean=%.1f, sd=%.1f\n",
            length(b12_raw),
            min(b12_raw), quantile(b12_raw, 0.25),
            median(b12_raw), quantile(b12_raw, 0.75),
            quantile(b12_raw, 0.95), quantile(b12_raw, 0.99),
            max(b12_raw), mean(b12_raw), stats::sd(b12_raw)))
# Cap at 1500 pmol/L (generous upper bound; normal range up to ~700)
B12_CAP <- 1500
n_above <- sum(b12_raw > B12_CAP, na.rm = TRUE)
cat(sprintf("Capping at %d pmol/L would affect %d / %d values (%.1f%%)\n",
            B12_CAP, n_above, length(b12_raw), 100 * n_above / length(b12_raw)))

# Build per-country bio with FIXED Ghana
bio_b12_raw <- list(
  Ghana       = admin2_mean(merged$Ghana,        "gw_B12_pmol_L"),
  SierraLeone = admin2_mean(merged$SierraLeone,  "gw_B12"),
  Malawi      = admin2_mean(merged$Malawi,       "vitb12")
)
# Cap Ghana B12 BEFORE aggregation
m_gh_fixed <- m_gh
m_gh_fixed$gw_B12_pmol_L <- pmin(as.numeric(unclass(m_gh_fixed$gw_B12_pmol_L)),
                                   B12_CAP)
bio_b12_fixed <- list(
  Ghana       = admin2_mean(m_gh_fixed,          "gw_B12_pmol_L"),
  SierraLeone = admin2_mean(merged$SierraLeone,  "gw_B12"),
  Malawi      = admin2_mean(merged$Malawi,       "vitb12")
)
cat("\nGhana admin-2 mean B12 (raw vs capped):\n")
cat(sprintf("  raw:    mean %.1f, sd %.1f\n",
            mean(bio_b12_raw$Ghana$mean_bio), stats::sd(bio_b12_raw$Ghana$mean_bio)))
cat(sprintf("  capped: mean %.1f, sd %.1f\n",
            mean(bio_b12_fixed$Ghana$mean_bio), stats::sd(bio_b12_fixed$Ghana$mean_bio)))

pooled_b12_raw   <- build_pooled(bio_b12_raw)
pooled_b12_fixed <- build_pooled(bio_b12_fixed)
res_b12_raw   <- if (!is.null(pooled_b12_raw))   run_loco(pooled_b12_raw,   "b12_raw")   else data.frame()
res_b12_fixed <- if (!is.null(pooled_b12_fixed)) run_loco(pooled_b12_fixed, "b12_capped") else data.frame()
cat("\nB12 LOCO (continuous mean, Spearman r):\n")
print(rbind(res_b12_raw, res_b12_fixed), row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────
cat("\n##### FIX 2: Folate biomarker harmonisation #####\n")
# Drop Malawi (rbcf is RBC folate, not serum folate)
bio_fol_all <- list(
  Ghana       = admin2_mean(merged$Ghana,        "gw_Folate_nmol_L"),
  SierraLeone = admin2_mean(merged$SierraLeone,  "gw_wFolate"),
  Malawi      = admin2_mean(merged$Malawi,       "rbcf")
)
bio_fol_serum <- bio_fol_all[c("Ghana", "SierraLeone")]
cat("Folate columns + biomarker type:\n")
cat("  Ghana       gw_Folate_nmol_L  (serum folate, nmol/L)\n")
cat("  SierraLeone gw_wFolate         (serum folate, units unclear)\n")
cat("  Malawi      rbcf               (RED-BLOOD-CELL folate -- DROP)\n\n")
# Also check serum folate availability in Ghana — look for outliers
gh_fol <- suppressWarnings(as.numeric(unclass(merged$Ghana$gw_Folate_nmol_L)))
gh_fol <- gh_fol[is.finite(gh_fol)]
cat(sprintf("Ghana serum folate: median %.1f, p95 %.0f, p99 %.0f, max %.0f\n",
            median(gh_fol), quantile(gh_fol, 0.95),
            quantile(gh_fol, 0.99), max(gh_fol)))
# Cap Ghana folate outliers too (generous: 50 nmol/L is well above human range)
FOL_CAP <- 50
n_above_fol <- sum(gh_fol > FOL_CAP, na.rm = TRUE)
cat(sprintf("Capping at %g nmol/L affects %d / %d (%.1f%%)\n",
            FOL_CAP, n_above_fol, length(gh_fol),
            100 * n_above_fol / length(gh_fol)))

m_gh_fol_fixed <- merged$Ghana
m_gh_fol_fixed$gw_Folate_nmol_L <-
  pmin(as.numeric(unclass(m_gh_fol_fixed$gw_Folate_nmol_L)), FOL_CAP)
bio_fol_serum_fixed <- list(
  Ghana       = admin2_mean(m_gh_fol_fixed,      "gw_Folate_nmol_L"),
  SierraLeone = admin2_mean(merged$SierraLeone,  "gw_wFolate")
)
# With only 2 countries, "LOCO" = leave-one-out on the pair.
pooled_fol_serum <- build_pooled(bio_fol_serum_fixed)
res_fol <- if (!is.null(pooled_fol_serum) &&
                length(pooled_fol_serum$country_names) >= 2)
              run_loco(pooled_fol_serum, "folate_serum_2country") else data.frame()
cat("\nFolate LOCO (2-country, serum only, with Ghana cap):\n")
print(res_fol, row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────
cat("\n##### FIX 3: Vitamin A Malawi shift correction #####\n")
bio_vita_raw <- list(
  Gambia      = admin2_mean(merged$Gambia,      "gw_cRBP"),
  Ghana       = admin2_mean(merged$Ghana,       "gw_cRBP"),
  SierraLeone = admin2_mean(merged$SierraLeone, "gw_cRBP"),
  Malawi      = admin2_mean(merged$Malawi,      "rbp")
)
# Compute means
country_means_raw <- sapply(bio_vita_raw, function(x)
  if (!is.null(x)) mean(x$mean_bio) else NA_real_)
cat("Raw country mean RBP:\n")
print(round(country_means_raw, 3))

# Compute pooled mean of (Gambia, Ghana, SL) and shift Malawi to match
reference_mean <- mean(country_means_raw[c("Gambia","Ghana","SierraLeone")],
                        na.rm = TRUE)
malawi_mean <- country_means_raw[["Malawi"]]
shift <- reference_mean - malawi_mean
cat(sprintf("\nReference mean (G/Gh/SL): %.3f\nMalawi mean:              %.3f\nShift to apply to Malawi:%+.3f\n",
            reference_mean, malawi_mean, shift))

bio_vita_fixed <- bio_vita_raw
bio_vita_fixed$Malawi$mean_bio <- bio_vita_fixed$Malawi$mean_bio + shift
country_means_fixed <- sapply(bio_vita_fixed, function(x)
  if (!is.null(x)) mean(x$mean_bio) else NA_real_)
cat("\nShifted country mean RBP:\n")
print(round(country_means_fixed, 3))

pooled_vita_raw   <- build_pooled(bio_vita_raw)
pooled_vita_fixed <- build_pooled(bio_vita_fixed)
res_vita_raw   <- if (!is.null(pooled_vita_raw))   run_loco(pooled_vita_raw,   "vitA_raw")   else data.frame()
res_vita_fixed <- if (!is.null(pooled_vita_fixed)) run_loco(pooled_vita_fixed, "vitA_shifted") else data.frame()
cat("\nvitA LOCO (continuous mean RBP, Spearman r):\n")
print(rbind(res_vita_raw, res_vita_fixed), row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────
cat("\n##### FIX 4: Sierra Leone admin-3 polygons #####\n")
sl_gadm_a3 <- tryCatch(
  geodata::gadm("SLE", level = 3, path = here::here("data", "gadm")),
  error = function(e) NULL)
if (!is.null(sl_gadm_a3)) {
  cat(sprintf("Sierra Leone GADM level-3: %d polygons\n", length(sl_gadm_a3)))
  sf_a3 <- sf::st_as_sf(sl_gadm_a3)
  cat("Sample admin-3 names:\n")
  print(head(sf_a3$NAME_3, 12))
} else {
  cat("GADM level-3 unavailable for Sierra Leone.\n")
}
# Check whether merged_sierraleone has admin-3 column
sl_merged <- merged$SierraLeone
a3_cols <- grep("admin3|adm3|chiefdom|cluster", colnames(sl_merged),
                 ignore.case = TRUE, value = TRUE)
cat(sprintf("\nSierra Leone merged: %d possible admin-3 / cluster columns:\n",
            length(a3_cols)))
print(head(a3_cols, 10))

cat("\nVerdict for FIX 4: ")
if (!is.null(sl_gadm_a3) && length(sl_gadm_a3) > 14) {
  cat(sprintf("GADM admin-3 has %d polygons (vs %d at admin-2).\n",
              length(sl_gadm_a3), 14))
  cat("Building DHS admin-3 aggregates would require re-running cluster\n")
  cat("-> admin-3 spatial join. Outside the scope of this script; deferred\n")
  cat("as a recommendation: add 'sierraleone_admin3' as a separate pooled\n")
  cat("track if SL is the dominant uncertainty source for LOCO claims.\n")
} else {
  cat("No admin-3 expansion available -> Sierra Leone is structurally\n")
  cat("underpowered at n_test = 14 and should be reported as such, not\n")
  cat("used as the headline transportability test case.\n")
}

# ───────────────────────────────────────────────────────────────────────────
cat("\n##### CONSOLIDATED RESULTS #####\n")
all <- dplyr::bind_rows(
  if (nrow(res_b12_raw)   > 0) cbind(outcome="women_b12",  res_b12_raw)   else NULL,
  if (nrow(res_b12_fixed) > 0) cbind(outcome="women_b12",  res_b12_fixed) else NULL,
  if (nrow(res_fol)       > 0) cbind(outcome="women_folate", res_fol)     else NULL,
  if (nrow(res_vita_raw)  > 0) cbind(outcome="women_vitA", res_vita_raw)  else NULL,
  if (nrow(res_vita_fixed)> 0) cbind(outcome="women_vitA", res_vita_fixed)else NULL
)
readr::write_csv(all, file.path(SANDBOX_RESULTS, "20_upstream_fixes.csv"))
cat(sprintf("\nwritten: %s\n",
            file.path(SANDBOX_RESULTS, "20_upstream_fixes.csv")))

# Side-by-side per outcome
cat("\n===== before / after Spearman r per fix =====\n")
for (oc in unique(all$outcome)) {
  cat(sprintf("\n--- %s ---\n", oc))
  print(all[all$outcome == oc, c("label","held_out","n_test","spearman_r")],
        row.names = FALSE)
}
