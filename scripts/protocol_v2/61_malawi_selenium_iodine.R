# =============================================================================
# scripts/protocol_v2/61_malawi_selenium_iodine.R   [MW-SE / MW-IO]
#
# MALAWI-ONLY OUTCOMES, WITHIN COUNTRY: selenium and iodine
#
# The MNS 2015-16 measured plasma selenium in every group and urinary iodine
# in women and school-age children (the RA's literature set, 2026-09-15:
# Phiri et al. 2019 on selenium; Likoswe / Gebremedhin on zinc adjustment).
# Neither is measured in the other three surveys, so, like iodine in Gambia
# and Sierra Leone (IO-01, script 47), they have no leave-one-country-out test
# and are scored with the protocol's in-country arms only.
#
# Why selenium matters beyond one more cell: it is the one outcome in the
# panel with a known geochemical driver (soil pH / soil selenium; Phiri et al.
# showed plasma Se tracks soil type across Malawi), so it is the cleanest
# soil -> biomarker test case of the transport logic the project relies on.
#
# Targets (Malawi's Admin-2 = Traditional Authorities, >= 5 respondents):
#   child_selenium  share of preschool children with plasma Se < 84.6 ug/L;
#                   level = -mean(log Se)  (negated: higher = worse)
#   women_selenium  same for non-pregnant women
#   women_iodine    share of women with UIC < 100 ug/L; level = -mean(log UIC)
# Arms: null, jackknifed regional (district) mean, spatial smoother, domain
# index, spatial + domain; estimands: in-fill (5-fold x 10 draws over TAs) and
# leave-one-district-out. Plus the descriptive that motivates the outcome:
# the district-level correlation of Se status with soil pH.
#
#   Rscript -e "source('scripts/protocol_v2/61_malawi_selenium_iodine.R')"
# -> results/tables/protocol_v2/malawi_selenium_iodine_targets.csv
# -> results/tables/protocol_v2/malawi_selenium_iodine_in_country.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; REPS <- 10L; set.seed(20260915L)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
b <- BND[["malawi"]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
CENT <- data.frame(Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(unclass(x)))
wm <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
wmed <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (sum(ok) < 3) return(NA_real_); o <- order(x[ok]); cw <- cumsum(w[ok][o]) / sum(w[ok]); x[ok][o][which(cw >= 0.5)[1]] }

# ── individual-level data: the merged file, with the derived flags ───────────
cc <- get_country_configs()[["Malawi"]]
m <- readRDS(cc$data_path)
m <- derive_malawi_binary(m, "sel", "sel_def", cc$outcomes$child_selenium$cutoff)
m <- derive_malawi_binary(m, "iod", "iod_def", cc$outcomes$women_iodine$cutoff)
m$w <- num(m[[cc$weight_col]]); m$w[!is.finite(m$w) | m$w <= 0] <- NA
m$Admin1 <- as.character(m$Admin1); m$Admin2 <- as.character(m$Admin2)
CELLS <- list(
  child_selenium = list(pop = "preschool children", cont = "sel", bin = "sel_def"),
  women_selenium = list(pop = "women",              cont = "sel", bin = "sel_def"),
  women_iodine   = list(pop = "women",              cont = "iod", bin = "iod_def"))
for (nm in names(CELLS)) { ce <- CELLS[[nm]]; d <- m[m$population == ce$pop, ]
  cat(sprintf("%-15s n = %4d with %s | weighted prevalence %.1f%% | median %s %.1f\n", nm, sum(is.finite(d[[ce$bin]])), ce$cont,
              100 * wm(d[[ce$bin]], d$w), ce$cont, wmed(num(d[[ce$cont]]), d$w))) }

# ── district (TA) targets ────────────────────────────────────────────────────
tg_rows <- list()
for (nm in names(CELLS)) { ce <- CELLS[[nm]]; d <- m[m$population == ce$pop, ]
  x <- num(d[[ce$cont]]); x[is.finite(x) & x <= 0] <- NA
  t <- data.frame(Admin1 = d$Admin1, Admin2 = d$Admin2, y = num(d[[ce$bin]]), x = x, w = d$w) |>
    group_by(Admin1, Admin2) |>
    summarise(n = sum(is.finite(y)), n_eff = kish_n_v2(w[is.finite(y)]), p_low = wm(y, w),
              median_biomarker = wmed(x, w), mean_log = wm(log(x), w), .groups = "drop") |>
    filter(n >= 5) |> mutate(outcome = nm, .before = 1)
  tg_rows[[nm]] <- t
  cat(sprintf("%-15s %3d TAs with >= 5 respondents in %2d districts\n", nm, nrow(t), n_distinct(t$Admin1))) }
TGT <- bind_rows(tg_rows)
write.csv(TGT, file.path(OUTDIR, "malawi_selenium_iodine_targets.csv"), row.names = FALSE)

# ── descriptive: selenium status against soil chemistry across TAs ───────────
rows <- list()
soil_cols <- intersect(c("soilgrids_ph", "soilgrids_organic_carbon", "soilgrids_clay", "soilgrids_cec"), names(S))
for (nm in c("child_selenium", "women_selenium")) { t <- TGT[TGT$outcome == nm, ] |> inner_join(S[S$country == "Malawi", c("Admin1", "Admin2", soil_cols)], by = c("Admin1", "Admin2"))
  for (sc in soil_cols) { r <- suppressWarnings(stats::cor(t[[sc]], t$mean_log, method = "spearman", use = "complete.obs"))
    cat(sprintf("%-15s Spearman(%s, mean log Se) %+.2f over %d TAs\n", nm, sc, r, sum(is.finite(t[[sc]]) & is.finite(t$mean_log))))
    rows[[length(rows) + 1L]] <- data.frame(outcome = nm, analysis = "descriptive", target = "mean_log", arm = sc, estimand = "district_correlation", spearman = r, n_units = nrow(t), stringsAsFactors = FALSE) } }

# ── protocol arms (in-country) ───────────────────────────────────────────────
for (nm in names(CELLS)) { t <- TGT[TGT$outcome == nm, ]
  mm <- inner_join(t, S[S$country == "Malawi", c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |> inner_join(CENT, by = c("Admin1", "Admin2")); mm <- mm[is.finite(mm$lon), ]
  Xr <- prep_predictors_v2(as.matrix(mm[, PREDS, drop = FALSE])); n <- nrow(mm); if (n < 8) next
  targets <- list(p_low = list(y = mm$p_low, ymod = .v2_logit(pmin(pmax(mm$p_low, 0.005), 0.995)), scale = "prev"),
                  level = list(y = -mm$mean_log, ymod = -mm$mean_log, scale = "level"))   # negated: higher = worse
  for (tg in names(targets)) { Y <- targets[[tg]]$ymod; ynat <- targets[[tg]]$y; ok <- is.finite(Y); aux <- list(lon = mm$lon, lat = mm$lat, Admin1 = mm$Admin1, y_nat = Y)
    for (est in c("infill", "region")) { if (est == "infill" && n < 15) next; if (est == "region" && dplyr::n_distinct(mm$Admin1) < 3) next
      for (arm in c("null_train_mean", "region_mean_jk", "spatial", "domain_index", "spatial_plus_domain")) { if (arm == "region_mean_jk" && est == "region") next
        sp <- c()
        for (r in seq_len(if (est == "infill") REPS else 1L)) { folds <- if (est == "infill") make_folds_v2("kfold_district", n, k = 5, rep_id = r) else as.integer(factor(mm$Admin1)); pred <- rep(NA_real_, n)
          for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f & ok); if (length(tr) < 8) next
            D <- domain_representation_v2(Xr, domain_of, sign_rows = tr); p <- tryCatch(ARMS_V2[[arm]](tr, te, Y, Xr, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) pred[te] <- p }
          sp <- c(sp, score_v2(ynat[ok], pred[ok], mm$n_eff[ok], scale = targets[[tg]]$scale)$spearman) }
        rows[[length(rows) + 1L]] <- data.frame(outcome = nm, analysis = "protocol", target = tg, arm = arm, estimand = est, spearman = median(sp, na.rm = TRUE), n_units = n, stringsAsFactors = FALSE)
        cat(sprintf("%-15s %-6s %-7s %-20s Spearman %+.3f (n = %d TAs)\n", nm, tg, est, arm, median(sp, na.rm = TRUE), n)) } } } }
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "malawi_selenium_iodine_in_country.csv"), row.names = FALSE)
cat("\n===== MW-SE / MW-IO summary =====\n"); print(as.data.frame(R |> mutate(spearman = round(spearman, 3))), row.names = FALSE); cat("\nDONE\n")
