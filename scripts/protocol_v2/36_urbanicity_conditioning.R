# =============================================================================
# scripts/protocol_v2/36_urbanicity_conditioning.R   [UR-01]
#
# IS THE TRANSPORTED SIGNAL MORE THAN AN URBAN-RURAL GRADIENT?
#
# The associations that replicate across the four countries -- legume
# cultivation, cattle, soil chemistry, night-time temperature, all tracking
# MORE deficiency -- describe where subsistence agriculture is. The obvious
# rival explanation for a climate + soil index that transports is that it is
# an urbanicity map with extra steps: cities sit on particular soils and
# climates, and deficiency is lower in cities. If so, one free urbanicity
# raster is the honest model and the agro-ecology reading is decoration.
#
# Urbanicity composite (all remotely sensed, no survey input): night lights
# (ntl_ccnl), population density, built surface, urban land-cover fraction,
# and travel time to healthcare (negated), each rank-normalised within
# country and averaged. Then, leave-one-country-out at both tiers:
#   rho_index         transported Spearman of the climate+soil (and full) index
#   rho_urb_only      an index built from urbanicity alone (sign learned on
#                     the training countries)
#   rho_index_resid   the index residualised on urbanicity using TRAINING-
#                     country coefficients -- usable at prediction time
#   partial_index     partial Spearman of outcome and index given urbanicity
#                     in the held-out country -- the diagnostic
#   rho_index_plus_urb  index and urbanicity together
# If partial_index ~ rho_index the signal is not urbanicity; if it collapses
# to zero while rho_urb_only ~ rho_index, the product is an urban-rural map.
#
#   Rscript scripts/protocol_v2/36_urbanicity_conditioning.R
# -> results/tables/protocol_v2/urbanicity_conditioning.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; MIN_TRAIN <- 20L; set.seed(20260904L)
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
CS <- c("Climate and weather", "Soil characteristics")
URB <- c(ntl_ccnl = 1, popdens_y2015 = 1, wpop_log_density_survey_year = 1, built_surface = 1, lcover_urban_frac_t0 = 1, access_healthcare_min = -1)
cat("urbanicity columns present:", paste(intersect(names(URB), PREDS), collapse = ", "), "\n")

build <- function(cn, on, tier) {
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont) & t$n_eff_cont > 0, ]
  if (nrow(t) < 12) return(NULL)
  if (tier == "admin1") {
    a <- t |> group_by(Admin1) |> summarise(y = stats::weighted.mean(y_level, n_eff_cont), yp = stats::weighted.mean(y_prev, n_eff_cont), w = sum(n_eff_cont), .groups = "drop")
    x <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1"); region <- m$Admin1
  } else {
    m <- t[, c("Admin1", "Admin2", "y_level", "y_prev", "n_eff_cont")]; names(m)[3:5] <- c("y", "yp", "w")
    m <- m |> inner_join(S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")); region <- m$Admin1
  }
  m <- m[is.finite(m$y), ]; if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  uc <- intersect(names(URB), colnames(Xr)); if (length(uc) < 2) return(NULL)
  urb <- rowMeans(sweep(Xr[, uc, drop = FALSE], 2, URB[uc], "*")); urb <- as.numeric(scale(urb))
  list(country = cn, n = nrow(m), y = m$y, yp = m$yp, X = Xr, w = m$w, urb = urb, region = region, n_urb_cols = length(uc))
}
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 4 || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) NA_real_ else suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman")) }
partial_sp <- function(y, p, u) { ok <- is.finite(y) & is.finite(p) & is.finite(u); if (sum(ok) < 5) return(NA_real_)
  ry <- stats::resid(stats::lm(rank(y[ok]) ~ rank(u[ok]))); rp <- stats::resid(stats::lm(rank(p[ok]) ~ rank(u[ok])))
  if (stats::sd(ry) == 0 || stats::sd(rp) == 0) NA_real_ else suppressWarnings(stats::cor(ry, rp)) }

rows <- list()
for (tier in c("admin1", "admin2")) for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on, tier), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  yfun <- function(z) if (target == "level") z$y else .v2_logit(z$yp)
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(yfun(z))))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) if (target == "level") z$y else z$yp))
  urb <- unlist(lapply(cl, function(z) z$urb)); all_rows <- seq_along(Y)
  aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$region))), y_nat = Y)
  for (h in names(cl)) { te <- which(ctry == h); tr <- which(ctry != h); if (length(tr) < MIN_TRAIN) next
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr); dom <- col_domain(colnames(Dm))
    D_cs <- Dm[, dom %in% CS, drop = FALSE]; U <- matrix(urb, ncol = 1, dimnames = list(NULL, "urb"))
    pr <- function(D) { p <- tryCatch(ARMS_V2[["domain_index"]](tr, all_rows, Y, NULL, D, aux), error = function(e) NULL); if (is.null(p) || length(p) != length(all_rows)) rep(NA_real_, length(all_rows)) else p }
    p_full <- pr(Dm); p_cs <- pr(D_cs); p_urb <- pr(U); p_cs_u <- pr(cbind(D_cs, U)); p_full_u <- pr(cbind(Dm, U))
    res_tr <- function(p) { if (!all(is.finite(p[tr]))) return(rep(NA_real_, length(p))); b <- stats::coef(stats::lm(p[tr] ~ urb[tr])); p - (b[1] + b[2] * urb) }
    y_te <- ynat[te]; u_te <- urb[te]
    rows[[length(rows) + 1L]] <- data.frame(tier = tier, target = target, outcome = on, heldout = h, n_units = length(te),
      rho_full = sp(y_te, p_full[te]), rho_cs = sp(y_te, p_cs[te]), rho_urb_only = sp(y_te, p_urb[te]),
      rho_full_resid = sp(y_te, res_tr(p_full)[te]), rho_cs_resid = sp(y_te, res_tr(p_cs)[te]),
      partial_full = partial_sp(y_te, p_full[te], u_te), partial_cs = partial_sp(y_te, p_cs[te], u_te),
      rho_full_plus_urb = sp(y_te, p_full_u[te]), rho_cs_plus_urb = sp(y_te, p_cs_u[te]),
      cor_y_urb = sp(y_te, u_te), cor_cs_urb = sp(p_cs[te], u_te), cor_full_urb = sp(p_full[te], u_te), stringsAsFactors = FALSE)
  }
  cat("done", tier, target, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "urbanicity_conditioning.csv"), row.names = FALSE)
cat("\n===== UR-01: transported signal with and without urbanicity =====\n")
ARMS <- c("rho_cs", "partial_cs", "rho_cs_resid", "rho_cs_plus_urb", "rho_full", "partial_full", "rho_full_resid", "rho_full_plus_urb", "rho_urb_only")
for (tier in unique(R$tier)) for (tg in unique(R$target)) { d <- R[R$tier == tier & R$target == tg, ]
  cat(sprintf("\n-- %s | %s target | %d cells --\n", tier, tg, nrow(d)))
  s <- data.frame(arm = ARMS, mean_rho = round(vapply(ARMS, function(a) mean(d[[a]], na.rm = TRUE), 0), 3),
                  median_rho = round(vapply(ARMS, function(a) median(d[[a]], na.rm = TRUE), 0), 3),
                  positive = vapply(ARMS, function(a) sum(d[[a]] > 0, na.rm = TRUE), 0L), n = vapply(ARMS, function(a) sum(is.finite(d[[a]])), 0L))
  print(s, row.names = FALSE)
  cat(sprintf("  climate+soil: partial - raw  median %+.3f (partial >= raw in %d of %d) | resid - raw median %+.3f | urb-only - cs median %+.3f\n",
              median(d$partial_cs - d$rho_cs, na.rm = TRUE), sum(d$partial_cs >= d$rho_cs, na.rm = TRUE), sum(is.finite(d$partial_cs - d$rho_cs)),
              median(d$rho_cs_resid - d$rho_cs, na.rm = TRUE), median(d$rho_urb_only - d$rho_cs, na.rm = TRUE)))
  cat(sprintf("  how urban is the index: median Spearman(cs index, urbanicity) %+.2f | (full index, urbanicity) %+.2f | (outcome, urbanicity) %+.2f\n",
              median(d$cor_cs_urb, na.rm = TRUE), median(d$cor_full_urb, na.rm = TRUE), median(d$cor_y_urb, na.rm = TRUE)))
}
cat("\n-- per cell, admin2 level: raw vs partial (climate+soil) --\n")
print(as.data.frame(R[R$tier == "admin2" & R$target == "level", c("outcome", "heldout", "n_units", "rho_cs", "partial_cs", "rho_cs_resid", "rho_urb_only", "cor_cs_urb", "cor_y_urb")] |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nDONE\n")
