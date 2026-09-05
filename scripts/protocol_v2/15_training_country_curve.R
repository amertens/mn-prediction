# =============================================================================
# scripts/protocol_v2/15_training_country_curve.R   [GAP 1]
#
# HOW MUCH DOES EACH ADDITIONAL TRAINING COUNTRY BUY?
#
# The NCE's first proposed activity is to add national survey datasets
# (Ethiopia, Pakistan). That ask is currently supported by an assertion. With
# four countries in hand we can MEASURE the slope: hold one country out, then
# train on every subset of the remaining three of size 1, 2 and 3, and see how
# transported ranking accuracy moves with the number of training countries.
#
#   held-out h  x  subsets of size 1 (3), 2 (3), 3 (1)  =  7 fits per h
#   4 held-out countries                                = 28 fits per outcome
#
# A positive slope converts "fund us to add countries" into "each added country
# buys this much". A flat slope is equally worth knowing before submission,
# because it would mean the ask needs a different justification.
#
# PROTOCOL. Identical to the country estimand in 02b_merge_and_loco.R:
#   - outcome z-scored WITHIN country, so this is a ranking claim only;
#   - predictors rank-normalised within country, pooled on common columns;
#   - domain PCs oriented from the TRAINING countries only (sign_rows = tr),
#     the per-country orientation hazard documented in PROTOCOL_V2.md.
#
#   Rscript scripts/protocol_v2/15_training_country_curve.R
# -> results/tables/protocol_v2/training_country_curve.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
MIN_TRAIN <- 20L
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

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  if (!all(c(ycol, wcol) %in% names(t))) return(NULL)
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1","Admin2","lon","lat")],
               by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  yn <- m[[ycol]]
  list(country = cn, n = nrow(m), y_nat = yn,
       y_mod = if (target == "prev") .v2_logit(yn) else yn,
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}

rows <- list()
for (target in c("level", "prev")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
      if (!is.null(z)) cl[[cn]] <- z
    }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next

    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat))
    wv   <- unlist(lapply(cl, function(z) z$w))
    aux  <- list(lon = unlist(lapply(cl, function(z) z$lon)),
                 lat = unlist(lapply(cl, function(z) z$lat)),
                 Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))),
                 y_nat = Y)

    for (h in names(cl)) {
      pool <- setdiff(names(cl), h)
      te <- which(ctry == h)
      for (k in seq_along(pool)) {
        subs <- utils::combn(pool, k, simplify = FALSE)
        for (sb in subs) {
          tr <- which(ctry %in% sb)
          if (length(tr) < MIN_TRAIN) next
          Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
          for (a in c("domain_index", "domain_enet")) {
            p <- tryCatch(ARMS_V2[[a]](tr, te, Y, Xm, Dm, aux),
                          error = function(e) rep(NA_real_, length(te)))
            if (length(p) != length(te)) next
            s <- score_v2(ynat[te], p, wv[te],
                          scale = if (target == "prev") "prev" else "level")
            rows[[length(rows) + 1L]] <- data.frame(
              target = target, outcome = on, heldout = h, arm = a,
              n_train_countries = k, train_set = paste(sb, collapse = "+"),
              n_train_areas = length(tr), n_areas = length(te),
              spearman = s$spearman, topk = s$topk, stringsAsFactors = FALSE)
          }
        }
      }
    }
    cat("curve done", target, on, "\n")
  }
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "training_country_curve.csv"), row.names = FALSE)

cat("\n===== TRANSPORT vs NUMBER OF TRAINING COUNTRIES =====\n")
for (tg in unique(R$target)) for (a in unique(R$arm)) {
  d <- R[R$target == tg & R$arm == a, ]
  if (!nrow(d)) next
  s <- d |> group_by(n_train_countries) |>
    summarise(fits = dplyr::n(),
              mean_rho = round(mean(spearman, na.rm = TRUE), 3),
              median_rho = round(median(spearman, na.rm = TRUE), 3),
              pct_positive = round(100 * mean(spearman > 0, na.rm = TRUE)),
              .groups = "drop")
  cat("\n--", tg, "/", a, "--\n"); print(as.data.frame(s), row.names = FALSE)
  if (nrow(s) > 1) {
    fit <- lm(spearman ~ n_train_countries, data = d)
    cat(sprintf("slope %+.4f per added country (p = %.3f)\n",
                coef(fit)[2], summary(fit)$coefficients[2, 4]))
  }
}
cat("\nDONE\n")
