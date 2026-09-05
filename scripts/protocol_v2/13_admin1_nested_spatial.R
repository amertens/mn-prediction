# =============================================================================
# scripts/protocol_v2/13_admin1_nested_spatial.R
#
# DOES ADDING COVARIATES TO A SPATIAL SMOOTHER HELP AT ADMIN-1?
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# Two gaps in what we have:
#
#   1. The protocol-v2 nested comparison (spatial vs spatial + covariates, where
#      the covariates are fitted to the SPATIAL MODEL'S RESIDUALS) has only ever
#      been run at ADMIN-2. It found no benefit. Whether that holds at ADMIN-1,
#      where a smoother has far fewer neighbours to interpolate between, is a
#      different question and is the one asked here.
#
#   2. scripts/covariates/13_resolution_comparison.R does compare admin-1 and
#      admin-2, but its covariate arm is the older NNLS ensemble on SCREENED RAW
#      COLUMNS. It has never been run with the domain-specific dimension
#      reduction (PCs per domain to 80% of variance) that protocol v2 made the
#      default. So the admin-1 result on record does not test the current model.
#
# WHY GAMBIA AND SIERRA LEONE MAY STILL NOT REPORT
# ------------------------------------------------
# It is NOT only the evaluation metric. A correlation on 4-6 held-out points is
# indeed uninterpretable, and MAE would be well defined -- but the protocol-v2
# learners themselves refuse to fit at small n:
#
#     .v2_spatial_fit : nrow(dtr) < 12          -> returns the training mean
#     .v2_enet        : nrow(Xtr) < 12          -> returns the training mean
#
# Under leave-one-out that means Gambia (6 regions -> 5 training) and Sierra
# Leone (4 -> 3) collapse BOTH arms to a constant, so spatial, spatial+domain
# and the null become the same number. The comparison is vacuous rather than
# noisy. This script reports that degeneracy explicitly instead of hiding it
# behind an "estimable" flag, and marks any cell where the arms are identical.
#
# SCORING
# Spearman is reported where there are at least 8 units, and MAE and a paired
# per-unit |error| sign test are reported ALWAYS, because those are defined at
# any n and are what a small-n country can contribute to.
#
#   Rscript scripts/protocol_v2/13_admin1_nested_spatial.R
# -> results/tables/protocol_v2/admin1_nested_spatial.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
MIN_R  <- 8L          # below this, report MAE only -- a correlation is noise

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

# Admin-1 centroids, from the Admin-2 polygon centroids
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
})) |>
  group_by(country, Admin1) |>
  summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE),
            .groups = "drop")

#' Aggregate an outcome and its predictors from Admin-2 up to Admin-1.
#' Outcome is precision-weighted by effective n; predictors are simple means,
#' matching how the signal probes aggregate.
build_a1 <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]
  if (!nrow(t)) return(NULL)
  a1 <- t |> group_by(Admin1) |>
    summarise(y = stats::weighted.mean(y_prev, n_eff),
              n_eff = sum(n_eff), .groups = "drop") |>
    filter(is.finite(y))
  sc <- S[S$country == cn, c("Admin1", PREDS)]
  x1 <- sc |> group_by(Admin1) |>
    summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  m <- a1 |> inner_join(x1, by = "Admin1") |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "lon", "lat")], by = "Admin1")
  m <- m[is.finite(m$lon), ]
  if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 2) return(NULL)
  list(n = nrow(m), y = m$y, Admin1 = m$Admin1,
       X = Xr, D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1))
}

rows <- list()
cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)

for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build_a1(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  ymod <- .v2_logit(cl$y)
  n <- cl$n

  arms <- c("null_train_mean", "spatial", "spatial_plus_domain", "domain_index")
  pred <- lapply(arms, function(a) rep(NA_real_, n))
  names(pred) <- arms

  # exhaustive leave-one-region-out
  for (k in seq_len(n)) {
    te <- k; tr <- setdiff(seq_len(n), k)
    for (a in arms) {
      p <- tryCatch(ARMS_V2[[a]](tr, te, ymod, cl$X, cl$D, cl$aux),
                    error = function(e) NA_real_)
      if (length(p) == 1L && is.finite(p)) pred[[a]][k] <- p
    }
  }

  obs <- cl$y
  err <- lapply(pred, function(p) abs(obs - .v2_expit(p)))

  # Did the learners actually fit, or did the n<12 guards collapse them?
  degenerate <- isTRUE(all.equal(pred$spatial, pred$spatial_plus_domain)) &&
                isTRUE(all.equal(pred$spatial, pred$null_train_mean))

  for (a in arms) {
    p <- pred[[a]]; ok <- is.finite(p) & is.finite(obs)
    rows[[length(rows) + 1L]] <- data.frame(
      country = cn, outcome = on, n_units = n, arm = a,
      spearman = if (sum(ok) >= MIN_R && stats::sd(p[ok]) > 0)
        round(suppressWarnings(stats::cor(obs[ok], p[ok], method = "spearman")), 3)
        else NA_real_,
      mae_pp = if (any(ok)) round(100 * mean(err[[a]][ok], na.rm = TRUE), 2) else NA_real_,
      degenerate = degenerate, stringsAsFactors = FALSE)
  }

  # the nested question: does adding covariates to the spatial fit help?
  d <- err$spatial - err$spatial_plus_domain      # >0 = covariates helped
  ok <- is.finite(d)
  rows[[length(rows) + 1L]] <- data.frame(
    country = cn, outcome = on, n_units = n, arm = "DIFF spatial+dom vs spatial",
    spearman = NA_real_,
    mae_pp = if (any(ok)) round(100 * mean(d[ok]), 2) else NA_real_,
    degenerate = degenerate, stringsAsFactors = FALSE)
  cat(sprintf("%-12s %-13s n=%2d %s\n", cn, on, n,
              if (degenerate) "DEGENERATE (n<12 guards)" else "fitted"))
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "admin1_nested_spatial.csv"), row.names = FALSE)

cat("\n===== ADMIN-1: does adding covariates to a spatial smoother help? =====\n")
cat("(mae_pp on the DIFF row is mean |err_spatial| - |err_spatial+domain|;\n")
cat(" POSITIVE = covariates helped. Degenerate cells are excluded.)\n\n")
D <- R[R$arm == "DIFF spatial+dom vs spatial" & !R$degenerate, ]
if (nrow(D)) {
  print(as.data.frame(D[order(-D$mae_pp), c("country","outcome","n_units","mae_pp")]),
        row.names = FALSE)
  cat(sprintf("\ncells %d | covariates helped in %d | mean %+.3f pp | median %+.3f pp\n",
              nrow(D), sum(D$mae_pp > 0, na.rm = TRUE),
              mean(D$mae_pp, na.rm = TRUE), median(D$mae_pp, na.rm = TRUE)))
  nz <- D$mae_pp[is.finite(D$mae_pp) & D$mae_pp != 0]
  if (length(nz))
    cat(sprintf("sign test p = %.3f\n", binom.test(sum(nz > 0), length(nz))$p.value))
}
cat("\n--- degenerate cells (learners returned the training mean) ---\n")
dg <- unique(R[R$degenerate, c("country","outcome","n_units")])
if (nrow(dg)) print(as.data.frame(dg), row.names = FALSE) else cat("  none\n")
cat("\nDONE\n")
