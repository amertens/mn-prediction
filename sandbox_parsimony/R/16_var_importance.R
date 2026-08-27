# =============================================================================
# sandbox_parsimony/R/16_var_importance.R
#
# Which of the 16 a-priori constructs is actually carrying the signal?
#
# HOW THE 16 WERE CHOSEN (so the importance table is read correctly):
# they were written down from the causal story before any performance number
# was looked at -- one variable per distinct mechanism, no two measuring the
# same thing:
#
#   diet quality / agro-ecology  NDVI, precipitation, soil organic carbon,
#                                soil zinc, soil iron, soil nitrogen
#   market access / poverty      nightlights, population density, built
#                                surface, settlement degree, global human
#                                modification, elevation
#   inflammation & blood loss    Pf parasite rate, ITN use
#   inherited anaemia            HbS allele frequency, G6PD allele frequency
#
# The constraint was availability in all five countries. That is NOT a neutral
# constraint -- the five-country intersection has already lost travel time to
# cities, land surface temperature, TerraClimate, cropland fraction and
# vegetation dynamics (see FINDINGS.md section 4), all of which have a better
# claim on the diet-quality pathway than soil chemistry does. So this list is
# the best available set, not the best conceivable one.
#
# Importance is permutation importance, averaged over CV folds and repeats, and
# reported RELATIVE to the best variable in each cell so cells with different
# signal strength can be averaged.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr); library(ranger)})

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
audit <- read.csv("sandbox_parsimony/out/noise_audit.csv", stringsAsFactors = FALSE)

# short, readable labels for the constructs
LAB <- c(
  gee_ndvi = "NDVI (vegetation)", gee_elevation = "elevation",
  gee_popdensity = "population density", gee_ccnl = "nightlights",
  gee_globalhumanmodification = "human modification",
  gee_ghsl_smod = "settlement degree", gee_wsf = "built surface",
  gee_trmm = "precipitation",
  gee_soiltotalcarbon_mean_0_20 = "soil organic carbon",
  gee_soilzinc_mean_0_20 = "soil zinc", gee_soiliron_mean_0_20 = "soil iron",
  gee_soilnitrogen_mean_0_20 = "soil nitrogen",
  lon = "longitude", lat = "latitude")
lab_of <- function(v) {
  out <- LAB[v]
  out[is.na(out)] <- sub("^MAP_[0-9]+_", "", v[is.na(out)])
  unname(out)
}

rows <- list()
for (oc in names(pooled_all)) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  rmax <- audit$r_max_d15[audit$outcome == oc]
  cur <- curated_vars(P$predictors)
  for (ctry in P$countries) {
    rm_c <- audit$r_max_d15[audit$outcome == oc & audit$country == ctry]
    # only cells where there is signal to attribute
    if (!length(rm_c) || !is.finite(rm_c) || rm_c < 0.35) next
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < 25) next

    pp <- prep_X(dat, dat, P$predictors)
    sel <- intersect(cur, pp$vars)
    X <- cbind(pp$Xtr[, sel, drop = FALSE], lon = dat$lon, lat = dat$lat)
    df <- data.frame(y = .logit(dat$svy_prev), X, check.names = TRUE)
    imp <- matrix(NA_real_, nrow = 10, ncol = ncol(X),
                  dimnames = list(NULL, colnames(X)))
    for (i in 1:10) {
      rf <- tryCatch(ranger::ranger(y ~ ., data = df, num.trees = 1000,
                                    min.node.size = 5, importance = "permutation",
                                    case.weights = pmax(dat$n_svy, 1), seed = i),
                     error = function(e) NULL)
      if (!is.null(rf)) imp[i, names(rf$variable.importance)] <- rf$variable.importance
    }
    mi <- colMeans(imp, na.rm = TRUE)
    mi[!is.finite(mi)] <- 0
    rows[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, variable = names(mi),
      imp = as.numeric(mi),
      rel = as.numeric(mi) / max(mi[mi > 0], na.rm = TRUE),
      stringsAsFactors = FALSE)
    message(sprintf("%-14s %-12s done", oc, ctry))
  }
}

im <- bind_rows(rows)
write.csv(im, "sandbox_parsimony/out/var_importance.csv", row.names = FALSE)

cat("\n=== Permutation importance, averaged over signal-bearing cells ===\n")
cat("rel = importance relative to the strongest variable in that cell.\n")
cat("share_top5 = share of cells where the variable lands in that cell's top 5.\n\n")
s <- im |> group_by(variable) |>
  summarise(cells = n(),
            mean_rel = round(mean(pmax(rel, 0)), 3),
            share_top5 = round(mean(rel >= sort(rel, decreasing = TRUE)[1] * 0), 3),
            .groups = "drop")
# recompute share_top5 properly, per cell
top5 <- im |> group_by(outcome, country) |>
  mutate(rank = rank(-rel, ties.method = "min")) |> ungroup() |>
  group_by(variable) |> summarise(share_top5 = round(mean(rank <= 5), 2), .groups = "drop")
s <- s |> select(-share_top5) |> left_join(top5, by = "variable") |>
  mutate(construct = lab_of(variable)) |>
  select(construct, variable, cells, mean_rel, share_top5) |>
  arrange(desc(mean_rel))
print(as.data.frame(s), row.names = FALSE)

cat("\n=== Top 3 variables per signal-bearing cell ===\n")
t3 <- im |> group_by(outcome, country) |> slice_max(rel, n = 3) |>
  summarise(top3 = paste(lab_of(variable), collapse = ", "), .groups = "drop")
print(as.data.frame(t3), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/var_importance.csv")
