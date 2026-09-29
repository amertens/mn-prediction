# =============================================================================
# scripts/policy_deck/44_mnf15_v6_importance_appendix.R
#
# Three appendix slides of the v6 MNF15 talk redrawn on the post-fix tables
# (survey outcome fixes and full rebuild of 27-28 September 2026,
# logs/rr14_survey_fixes.log). The fourth, the Malawi B12 maps, is drawn by
# 44b_malawi_b12_maps_v6.R.
#
#   a6_top20_child_iron.png    "Top 20 individual predictors for child iron"
#       The twenty largest standardised weights in the pooled four-country
#       index, child iron, biomarker level. Label = share of the index's
#       variance; * = same sign in each country's own fit.
#       <- results/tables/protocol_v2/index_importance_columns.csv (scope pooled)
#   a6_malawi_b12_weights.png  "What are the weights behind the best ranking?"
#       Malawi, women's B12, biomarker level: the index refitted on B = 300
#       bootstrap resamples of the surveyed districts (domain components
#       re-oriented and weights re-learned each time); point = median weight,
#       bar = 5th to 95th percentile, hollow = full fit. This is the computation
#       of scripts/policy_deck/09_vim_forest_example.R, rerun here on the
#       post-fix targets_v2.csv and written only to the figure folder.
#   (text slide)               "What drives the transport models?"
#       Recurrence counts in the top twenty of the 22 leave-one-country-out
#       fits on the biomarker level, printed old -> new.
#       <- index_importance_columns.csv (scope loco)
#
# Every number is also computed on the pre-fix copy of the same tables
# (results/tables/protocol_v2_pre_RR11_20260927/) and printed old -> new. The
# importance and benchmark tables there predate the rebuild and reproduce the
# old slides' numbers exactly; its targets_v2.csv does not (it was copied after
# script 01 rewrote it), so nothing here reads it.
# Predictor tiers: open,survey_public, the headline set (script 57's post-fix
# run used the same; logs/rr11_57_index_importance.log).
#
#   Rscript scripts/policy_deck/44_mnf15_v6_importance_appendix.R
#   (optional A6_CACHE=<path.rds> caches the 300 refits between runs)
# -> results/figures/mnf15_v6/a6_top20_child_iron.png, a6_malawi_b12_weights.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R"); source("R/predictor_plain_names.R")
P2  <- "results/tables/protocol_v2"; PRE <- "results/tables/protocol_v2_pre_RR11_20260927"
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
rd  <- function(...) read.csv(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
hdr <- function(x) cat(sprintf("\n==== %s ====\n", x))

# plain names: R/predictor_plain_names.R first, then a few for columns that reach
# these two figures without one (kept here so no shared file is edited)
EXTRA <- c(mics_handwash_soap = "Handwashing place with soap (MICS)", hces_any_asf = "Households eating animal-source food (budget survey)",
           mics_c_mean_haz = "Mean child height-for-age (MICS)", wdist_any_km_mean = "Distance to any water")
lab_of <- function(cols) { p <- unname(PLAIN[cols]); e <- unname(EXTRA[cols]); out <- ifelse(!is.na(p), p, e)
  if (any(is.na(out))) { warning("no plain name for: ", paste(cols[is.na(out)], collapse = ", ")); out[is.na(out)] <- clean_code(cols[is.na(out)]) }
  out }

# one domain colour map for both figures; the first thirteen are the old child-iron
# figure's colours (its palette in sorted-domain order), so that slide is unchanged in look
DOMCOL <- c(
  "Agricultural production, land use" = "#1f4e79", "Built environment" = "#6baed6", "Child anthropometry" = "#8c510a",
  "Ecosystem productivity/greenness" = "#33a02c", "Education, employment, SES" = "#b15928", "Food prices and supply" = "#6a3d9a",
  "Household assets and characteristics" = "#e08214", "Infant and young child feeding" = "#0F7B8A",
  "Infection and inflammation burden" = "#cab2d6", "Livestock density" = "#fdbf6f", "Soil characteristics" = "#d7191c",
  "Water and coast proximity" = "#9e9ac8", "Water and sanitation" = "#7fcdbb",
  "Anaemia and haemoglobin" = "#e7298a", "Household diet and consumption (HCES)" = "#274C77",
  "Market prices (RTFP)" = "#8e6fbf", "Food fortification and supplementation" = "#a6d854")
SPARE <- c("#666666", "#bf812d", "#80cdc1", "#fb9a99", "#b2df8a")
dom_cols <- function(doms) { doms <- sort(unique(doms)); miss <- setdiff(doms, names(DOMCOL))
  if (length(miss)) { warning("no colour for domain(s): ", paste(miss, collapse = "; ")); DOMCOL <- c(DOMCOL, stats::setNames(SPARE[seq_along(miss)], miss)) }
  DOMCOL[doms] }

IC  <- rd(P2, "index_importance_columns.csv")
ICo <- rd(PRE, "index_importance_columns.csv")

# =============================================================================
# 1. Top 20 for child iron (pooled four-country fit, biomarker level)
# =============================================================================
hdr("Top 20 individual predictors for child iron")
top20 <- function(ic) ic |> filter(scope == "pooled", fit == "all", outcome == "child_iron", target == "level") |>
  mutate(n_cols = n()) |> arrange(rank) |> filter(rank <= 20) |>
  mutate(star = incountry_sign_agree == incountry_fits, loco_all = loco_sign_agree == loco_fits)
t_new <- top20(IC); t_old <- top20(ICo)
chk(nrow(t_new) == 20, "20 rows in the post-fix pooled child-iron fit")
IT <- rd(P2, "index_importance_top.csv") |> filter(outcome == "child_iron", target == "level", rank <= 20) |> arrange(rank)
chk(identical(IT$column, t_new$column) && max(abs(IT$beta_std - t_new$beta_std)) < 1e-9, "index_importance_top.csv agrees with the columns table")
cmp <- full_join(t_old |> transmute(column, rank_old = rank, w_old = round(beta_std, 2), share_old = round(100 * share, 1), star_old = star),
                 t_new |> transmute(column, rank_new = rank, w_new = round(beta_std, 2), share_new = round(100 * share, 1), star_new = star), by = "column") |>
  arrange(rank_new)
print(as.data.frame(cmp), row.names = FALSE)
cat(sprintf("columns in the pooled fit: old %d, new %d\n", t_old$n_cols[1], t_new$n_cols[1]))
cat(sprintf("share held by the top 20: old %.1f%%, new %.1f%%\n", 100 * sum(t_old$share), 100 * sum(t_new$share)))
cat(sprintf("starred (same sign in all four own-country fits): old %d, new %d of 20\n", sum(t_old$star), sum(t_new$star)))
cat(sprintf("same sign in every LOCO fit: old %d, new %d of 20\n", sum(t_old$loco_all), sum(t_new$loco_all)))
cat(sprintf("same twenty columns as before: %s\n", setequal(t_old$column, t_new$column)))

d1 <- t_new |> mutate(name = lab_of(column), lab = sprintf("%.1f%%%s", 100 * share, ifelse(star, "  *", "")),
                      name = factor(name, levels = rev(name)))
p1 <- ggplot(d1, aes(beta_std, name, fill = domain)) + geom_col(width = 0.72) + geom_vline(xintercept = 0, colour = "grey40") +
  geom_text(aes(label = lab, x = ifelse(beta_std > 0, beta_std + 0.1, beta_std - 0.1), hjust = ifelse(beta_std > 0, 0, 1)), size = 3.3, colour = "grey30") +
  scale_fill_manual(values = dom_cols(d1$domain)) + scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(x = "Standardised weight in the pooled index (left: district ranks better; right: worse)\nLabel: share of the index carried; * = same direction in each country's own fit", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 9.5)) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))
ggsave(file.path(OUT, "a6_top20_child_iron.png"), p1, width = 11, height = 6.4, dpi = 220, bg = "white")   # 7.91 x 4.6 box
cat("wrote a6_top20_child_iron.png\n")

# =============================================================================
# 2. What drives the transport models? (leave-one-country-out, level, top 20)
# =============================================================================
hdr("What drives the transport models? (LOCO top-20 recurrence)")
rec <- function(ic) ic |> filter(scope == "loco", target == "level", rank <= 20) |> group_by(column, domain) |>
  summarise(n_fits = n(), n_outcomes = n_distinct(outcome), n_pos = sum(beta_std > 0), n_neg = sum(beta_std < 0), .groups = "drop") |>
  mutate(unanimous = n_pos == 0 | n_neg == 0) |> arrange(desc(n_fits))
nfits <- function(ic) nrow(distinct(filter(ic, scope == "loco", target == "level"), fit, outcome))
R_new <- rec(IC); R_old <- rec(ICo)
cat(sprintf("held-out fits: old %d, new %d\n", nfits(ICo), nfits(IC)))
NAMED <- c(lcover_grass_frac_t0 = "grassland cover", lcover_crops_frac_t0 = "cropland cover", glw_ruminant_share = "ruminant share",
           glw_cattle_km2 = "cattle density", glw_sheep_km2 = "sheep density", spam_share_cereals = "cereal share of cropland",
           npp_npp_t0 = "net plant productivity", npp_gpp_t0 = "gross plant productivity",
           ihme_wastingprevalence = "modelled wasting", ihme_underweightprevalence = "modelled underweight", ihme_allanemia = "modelled anaemia (any)",
           ihme_mildanemia = "modelled mild anaemia", ihme_moderateanemia = "modelled moderate anaemia", ihme_severeanemia = "modelled severe anaemia",
           ihme_overweightprevalence = "modelled overweight", fprice_staple_rel = "relative staple price",
           wdist_coast_km_mean = "distance to coast (mean)", wdist_coast_km_min = "distance to coast (nearest)", wpop_dependency_ratio = "dependency ratio",
           map_sy_pf_mortality_rate = "malaria mortality", map_sy_pf_reproductive_number = "malaria reproductive number", map_sy_pf_parasite_rate = "malaria parasite rate")
nm <- data.frame(column = names(NAMED), what = unname(NAMED)) |>
  left_join(R_old |> transmute(column, fits_old = n_fits, outc_old = n_outcomes, sign_old = ifelse(n_pos > 0 & n_neg == 0, "+", ifelse(n_neg > 0 & n_pos == 0, "-", "mixed"))), by = "column") |>
  left_join(R_new |> transmute(column, fits_new = n_fits, outc_new = n_outcomes, sign_new = ifelse(n_pos > 0 & n_neg == 0, "+", ifelse(n_neg > 0 & n_pos == 0, "-", "mixed"))), by = "column")
print(nm, row.names = FALSE)
u5 <- function(R) sprintf("%d of the %d columns in five or more fits keep one sign", sum(R$unanimous[R$n_fits >= 5]), sum(R$n_fits >= 5))
cat("old:", u5(R_old), "\nnew:", u5(R_new), "\n")
cat("\npost-fix: columns in the top twenty of seven or more held-out fits\n")
print(as.data.frame(R_new |> filter(n_fits >= 7)), row.names = FALSE)
# household-survey items (MICS, HCES): how often they recur
sv <- R_new |> filter(grepl("^(mics|hces)_", column))
cat(sprintf("\nhousehold-survey items in any top twenty: %d columns; in five or more fits: %d; most: %s (%d)\n",
            nrow(sv), sum(sv$n_fits >= 5), sv$column[1], sv$n_fits[1]))
# the landscape gradient in every outcome?
LAND <- c("lcover_grass_frac_t0", "lcover_crops_frac_t0", "glw_ruminant_share", "glw_cattle_km2", "glw_sheep_km2", "spam_share_cereals", "npp_npp_t0", "npp_gpp_t0")
lo <- IC |> filter(scope == "loco", target == "level", rank <= 20) |> group_by(outcome) |>
  summarise(fits = n_distinct(fit), landscape_hits = sum(column %in% LAND), fits_with_landscape = n_distinct(fit[column %in% LAND]), .groups = "drop")
print(as.data.frame(lo), row.names = FALSE)

# =============================================================================
# 3. Malawi, women's B12: weights with 5th-95th percentile over 300 refits
# =============================================================================
hdr("Malawi women's B12: weights over bootstrap refits")
BC <- rd(P2, "benchmarks_v2_cells.csv"); BCo <- rd(PRE, "benchmarks_v2_cells.csv")
bestcells <- function(bc) bc |> filter(estimand == "infill", arm == "domain_index", target == "level") |> arrange(desc(spearman)) |>
  transmute(country, outcome, spearman = round(spearman, 3), rank = row_number())
cat("in-fill, index, level, old:\n"); print(head(bestcells(BCo), 4), row.names = FALSE)
cat("in-fill, index, level, new:\n"); print(head(bestcells(BC), 4), row.names = FALSE)

TG <- rd(P2, "targets_v2.csv")
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
CN <- "Malawi"; ON <- "women_b12"; B <- 300L
t <- TG[TG$country == CN & TG$outcome == ON, ]; t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
m <- inner_join(t, S[S$country == CN, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
X <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); y <- m$y_level; n <- nrow(X)
cat(sprintf("%s %s: %d surveyed districts, %d predictors after per-country preparation\n", CN, ON, n, ncol(X)))

D0 <- domain_representation_v2(X, domain_of)
full <- index_importance_v2(seq_len(n), y, X, D0)$columns
ref <- IC |> filter(scope == "country", fit == CN, outcome == ON, target == "level")
chk(nrow(ref) == ncol(X) && max(abs(ref$beta_std[match(full$column, ref$column)] - full$beta_std)) < 1e-8,
    "full fit reproduces the post-fix index_importance_columns.csv (scope country, Malawi women_b12, level)")
cat("full fit matches index_importance_columns.csv (scope country) exactly\n")

cache <- Sys.getenv("A6_CACHE", "")
if (nzchar(cache) && file.exists(cache)) { boot <- readRDS(cache); cat("refits read from", cache, "\n") } else {
  set.seed(20260916L)   # script 09's seed
  t0 <- Sys.time()
  boot <- sapply(seq_len(B), function(b) {
    tr <- sample.int(n, n, replace = TRUE)
    Db <- domain_representation_v2(X, domain_of, sign_rows = tr)
    im <- tryCatch(index_importance_v2(tr, y, X, Db)$columns, error = function(e) NULL)
    if (is.null(im)) return(rep(NA_real_, ncol(X)))
    im$beta_std[match(colnames(X), im$column)]
  })
  rownames(boot) <- colnames(X)
  cat(sprintf("%d refits in %.1f min\n", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (nzchar(cache)) saveRDS(boot, cache)
}
chk(ncol(boot) == B, "B refits")
cat(sprintf("failed refits: %d\n", sum(colSums(is.finite(boot)) == 0)))
summ <- data.frame(column = colnames(X), full = full$beta_std[match(colnames(X), full$column)],
                   median = apply(boot, 1, median, na.rm = TRUE), lo = apply(boot, 1, quantile, 0.05, na.rm = TRUE), hi = apply(boot, 1, quantile, 0.95, na.rm = TRUE),
                   sign_agree = apply(boot, 1, function(v) mean(sign(v) == sign(median(v, na.rm = TRUE)), na.rm = TRUE)), stringsAsFactors = FALSE)
summ$domain <- unname(domain_of[summ$column])
summ <- summ |> arrange(desc(abs(median))) |> mutate(rank = row_number())
d2 <- summ |> slice_head(n = 20) |> mutate(label = lab_of(column))

VO <- read.csv("results/tables/policy_deck/vim_forest_Malawi_women_b12.csv", stringsAsFactors = FALSE)   # the pre-fix run (16 Sep)
cmp2 <- full_join(VO |> filter(rank <= 20) |> transmute(column, rank_old = rank, med_old = round(median, 2), agree_old = round(sign_agree, 3)),
                  d2 |> transmute(column, rank_new = rank, med_new = round(median, 2), lo = round(lo, 2), hi = round(hi, 2), full_new = round(full, 2), agree_new = round(sign_agree, 3), domain), by = "column") |>
  arrange(rank_new)
print(as.data.frame(cmp2), row.names = FALSE)
cat(sprintf("top-20 overlap with the pre-fix figure: %d of 20\n", sum(d2$column %in% VO$column[VO$rank <= 20])))
cat(sprintf("lowest sign agreement in the top 20: old %.3f, new %.3f\n", min(VO$sign_agree[VO$rank <= 20]), min(d2$sign_agree)))
cat(sprintf("intervals crossing zero in the top 20: %d\n", sum(d2$lo < 0 & d2$hi > 0)))

d2 <- d2 |> mutate(label = factor(label, levels = rev(label)))
p2 <- ggplot(d2, aes(median, label)) + geom_vline(xintercept = 0, colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", colour = "grey60", linewidth = 0.8) +
  geom_point(aes(colour = domain), size = 3.2) + geom_point(aes(x = full), shape = 1, size = 3.2, colour = "grey20") +
  scale_colour_manual(values = dom_cols(d2$domain)) +
  labs(x = sprintf("Weight in the district index (median and 5th to 95th percentile over %d refits; hollow = full fit; right = worse B12 status)", B), y = NULL, colour = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.location = "plot", legend.text = element_text(size = 10), legend.key.height = grid::unit(12, "pt"),
        panel.grid.minor = element_blank(), axis.title.x = element_text(size = 10.5)) +
  guides(colour = guide_legend(ncol = 4, byrow = TRUE))
ggsave(file.path(OUT, "a6_malawi_b12_weights.png"), p2, width = 12, height = 5.6, dpi = 220, bg = "white")   # 9.86 x 4.6 box
cat("wrote a6_malawi_b12_weights.png\nDONE\n")
