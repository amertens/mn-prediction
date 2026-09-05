# =============================================================================
# scripts/protocol_v2/47_iodine_in_country.R   [IO-01]
#
# IODINE, WITHIN COUNTRY: THE ONE OUTCOME THAT CANNOT BE VALIDATED ACROSS BORDERS
#
# Iodine status was in the January 2026 outcome table (urinary iodine, Sierra
# Leone) and left the modelling outcomes because a one-country outcome has
# no leave-one-country-out test. Two surveys carry it:
#   Gambia 2018   continuous urinary iodine concentration (gw_wUIC, ug/L,
#                 1,285 women) and household salt iodine (gw_salt_iod_conc_ppm)
#   Sierra Leone  UIC in six WHO epidemiological categories for non-lactating
#   2013          women (gw_UICatNonLact: < 20, 20-49, 50-99, 100-199,
#                 200-299, >= 300 ug/L; category labels are assumed from the
#                 WHO scheme and the category frequencies, the codebook label
#                 set could not be extracted), with its own UIC weight
#                 (gw_wStatWtUIC), and household salt iodine (gw_SaltIod ppm)
# District targets: prevalence of UIC < 100 ug/L (insufficient intake,
# non-pregnant women) and, for Gambia, the district median UIC. Then the
# protocol's in-fill arms where the district count allows it (Gambia, 30
# districts; 5-fold x 10 draws) and leave-one-region-out in both countries;
# and the descriptive that matters for a programme: does district salt
# iodisation track district iodine status?
#
#   Rscript scripts/protocol_v2/47_iodine_in_country.R
# -> results/tables/protocol_v2/iodine_in_country.csv, iodine_targets.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; STORE <- "_targets_full"; REPS <- 10L; set.seed(20260904L)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD); domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
cent <- function(lc) { b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b)))); data.frame(Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
wm <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
wmed <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (sum(ok) < 3) return(NA_real_); o <- order(x[ok]); cw <- cumsum(w[ok][o]) / sum(w[ok]); x[ok][o][which(cw >= 0.5)[1]] }

# ── individual-level extraction ──────────────────────────────────────────────
ind <- list()
d <- tar_read_raw("outcome_data_gambia_women_vitA", store = STORE)$data
uic <- num(d$gw_wUIC); w <- num(d$gw_svy_weight); preg <- if ("gw_wPregnant" %in% names(d)) num(d$gw_wPregnant) else rep(0, nrow(d))
ind$Gambia <- data.frame(country = "Gambia", Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2), uic = uic, low = as.numeric(uic < ifelse(is.finite(preg) & preg == 1, 150, 100)), w = w, salt_ppm = num(d$gw_salt_iod_conc_ppm), salt_iodised = num(d$gw_hSaltIodYN), stringsAsFactors = FALSE)
cat(sprintf("Gambia: %d women with UIC | median UIC %.0f ug/L | UIC < 100: %.1f%% | households with iodised salt %.0f%%\n", sum(is.finite(uic)), stats::median(uic, na.rm = TRUE), 100 * mean(uic < 100, na.rm = TRUE), 100 * mean(ind$Gambia$salt_iodised == 1, na.rm = TRUE)))
d <- tar_read_raw("outcome_data_sierraleone_women_vitA", store = STORE)$data
cat_nl <- num(d$gw_UICatNonLact); w <- num(d$gw_wStatWtUIC); if (all(!is.finite(w))) w <- num(d$gw_svy_weight)
cat("Sierra Leone UIC category frequencies (assumed WHO bands <20 / 20-49 / 50-99 / 100-199 / 200-299 / >=300):", paste(names(table(cat_nl)), table(cat_nl), collapse = " "), "\n")
ind$SierraLeone <- data.frame(country = "SierraLeone", Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2), uic = NA_real_, low = ifelse(is.finite(cat_nl), as.numeric(cat_nl <= 3), NA_real_), w = w, salt_ppm = num(d$gw_SaltIod), salt_iodised = ifelse(is.finite(num(d$gw_SaltIodAdeq)), as.numeric(num(d$gw_SaltIodAdeq) == 1), NA_real_), stringsAsFactors = FALSE)
cat(sprintf("Sierra Leone: %d non-lactating women with a UIC category | UIC < 100: %.1f%% | households with adequately iodised salt %.0f%%\n", sum(is.finite(cat_nl)), 100 * mean(ind$SierraLeone$low, na.rm = TRUE), 100 * mean(ind$SierraLeone$salt_iodised == 1, na.rm = TRUE)))
IND <- bind_rows(ind)
# ── district targets ─────────────────────────────────────────────────────────
TGT <- IND |> group_by(country, Admin1, Admin2) |> summarise(n = sum(is.finite(low)), n_eff = kish_n_v2(w[is.finite(low)]), p_low = wm(low, w), median_uic = wmed(uic, w), mean_log_uic = wm(log(uic), w),
  salt_ppm = wmed(salt_ppm, w), share_iodised = wm(salt_iodised, w), .groups = "drop") |> filter(n >= 5)
write.csv(TGT, file.path(OUTDIR, "iodine_targets.csv"), row.names = FALSE)
cat(sprintf("\ndistrict targets: Gambia %d districts, Sierra Leone %d districts (>= 5 women)\n", sum(TGT$country == "Gambia"), sum(TGT$country == "SierraLeone")))

# ── descriptive: salt iodisation vs iodine status across districts ───────────
rows <- list()
for (cn in unique(TGT$country)) { t <- TGT[TGT$country == cn, ]
  r1 <- suppressWarnings(stats::cor(t$share_iodised, t$p_low, method = "spearman", use = "complete.obs")); r2 <- suppressWarnings(stats::cor(t$salt_ppm, t$p_low, method = "spearman", use = "complete.obs"))
  r3 <- if (cn == "Gambia") suppressWarnings(stats::cor(t$salt_ppm, t$median_uic, method = "spearman", use = "complete.obs")) else NA_real_
  cat(sprintf("%-12s Spearman(share iodised salt, UIC<100) %+.2f | (salt ppm, UIC<100) %+.2f | (salt ppm, median UIC) %s | districts %d\n", cn, r1, r2, if (is.finite(r3)) sprintf("%+.2f", r3) else "n/a", nrow(t)))
  rows[[length(rows) + 1L]] <- data.frame(country = cn, analysis = "descriptive", target = "p_low", arm = "salt_share_iodised", estimand = "district_correlation", spearman = r1, n_units = nrow(t), stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <- data.frame(country = cn, analysis = "descriptive", target = "p_low", arm = "salt_ppm", estimand = "district_correlation", spearman = r2, n_units = nrow(t), stringsAsFactors = FALSE) }

# ── protocol arms ────────────────────────────────────────────────────────────
lc_of <- c(Gambia = "gambia", SierraLeone = "sierraleone")
for (cn in unique(TGT$country)) { t <- TGT[TGT$country == cn, ]
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |> inner_join(cent(lc_of[[cn]]), by = c("Admin1", "Admin2")); m <- m[is.finite(m$lon), ]
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); n <- nrow(m); if (n < 8) next
  targets <- list(p_low = list(y = m$p_low, ymod = .v2_logit(pmin(pmax(m$p_low, 0.005), 0.995)), scale = "prev"))
  if (cn == "Gambia") targets$log_uic <- list(y = -m$mean_log_uic, ymod = -m$mean_log_uic, scale = "level")   # negated so higher = worse, as the other level targets
  for (tg in names(targets)) { Y <- targets[[tg]]$ymod; ynat <- targets[[tg]]$y; ok <- is.finite(Y); aux <- list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = Y)
    for (est in c("infill", "region")) { if (est == "infill" && n < 15) next; if (est == "region" && dplyr::n_distinct(m$Admin1) < 3) next
      for (arm in c("null_train_mean", "region_mean_jk", "spatial", "domain_index", "spatial_plus_domain")) { if (arm == "region_mean_jk" && est == "region") next
        sp <- c()
        for (r in seq_len(if (est == "infill") REPS else 1L)) { folds <- if (est == "infill") make_folds_v2("kfold_district", n, k = 5, rep_id = r) else as.integer(factor(m$Admin1)); pred <- rep(NA_real_, n)
          for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f & ok); if (length(tr) < 8) next
            D <- domain_representation_v2(Xr, domain_of, sign_rows = tr); p <- tryCatch(ARMS_V2[[arm]](tr, te, Y, Xr, D, aux), error = function(e) rep(NA_real_, length(te))); if (length(p) == length(te)) pred[te] <- p }
          sp <- c(sp, score_v2(ynat[ok], pred[ok], m$n_eff[ok], scale = targets[[tg]]$scale)$spearman) }
        rows[[length(rows) + 1L]] <- data.frame(country = cn, analysis = "protocol", target = tg, arm = arm, estimand = est, spearman = median(sp, na.rm = TRUE), n_units = n, stringsAsFactors = FALSE)
        cat(sprintf("%-12s %-8s %-7s %-20s Spearman %+.3f (n = %d)\n", cn, tg, est, arm, median(sp, na.rm = TRUE), n)) } } } }
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "iodine_in_country.csv"), row.names = FALSE)
cat("\n===== IO-01 summary =====\n"); print(as.data.frame(R |> mutate(spearman = round(spearman, 3))), row.names = FALSE); cat("\nDONE\n")
