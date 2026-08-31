# =============================================================================
# scripts/covariates/08_admin1_arms.R
#
# Evaluate the two admin-1 additions before the full pipeline rebuild:
#
#   ARM A  admin-1 BENCHMARKING  -- keep fitting at admin-2, but anchor each
#          region's predictions to its (shrunk) design-based regional total
#          instead of anchoring the whole country to one national number
#   ARM B  admin-1 DOWNSCALING   -- fit at admin-1 pooled across countries and
#          predict admin-2 (the Planetary Prediction Engine design)
#
# Both are scored against the admin-2 direct estimates with the noise ceiling
# r_max shown alongside, because at admin-2 that ceiling is the binding
# constraint on what any number here can mean.
#
# Output: results/tables/admin1_arms.csv
#   Rscript scripts/covariates/08_admin1_arms.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
targets::tar_source("R/")

STORE <- here("_targets_full")
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()
rows <- list()

score <- function(svy, pred_df, col, country, outcome, arm, rmax) {
  # Pair key where both sides carry Admin1. Joining on the district NAME alone
  # fans Malawi from 87 rows to 90, because six of its Admin-2 names occur in
  # more than one region; every other consumer of a GADM-derived covariate
  # table was migrated to admin2_join_by() and this one was missed.
  jb <- admin2_join_by(svy, pred_df)
  keep_s <- intersect(c("Admin1", "Admin2", "svy_prev", "n_svy"), names(svy))
  keep_p <- intersect(c("Admin1", "Admin2", col), names(pred_df))
  m <- dplyr::inner_join(svy[, keep_s, drop = FALSE],
                         pred_df[, keep_p, drop = FALSE], by = jb)
  names(m)[names(m) == col] <- "p"
  m <- m[is.finite(m$svy_prev) & is.finite(m$p), , drop = FALSE]
  if (nrow(m) < 5) return(NULL)
  r <- suppressWarnings(stats::cor(m$svy_prev, m$p))
  data.frame(country = country, outcome = outcome, arm = arm, n_admin2 = nrow(m),
             pearson_r = round(r, 3), r_max = round(rmax, 3),
             r_share = if (is.finite(rmax) && rmax > 0) round(r / rmax, 2) else NA_real_,
             mae_pp = round(100 * mean(abs(m$svy_prev - m$p)), 2),
             bias_pp = round(100 * mean(m$p - m$svy_prev), 2),
             stringsAsFactors = FALSE)
}

for (ctry in names(cfgs)) {
  cc <- cfgs[[ctry]]
  # Deduplicate on the PAIR, not the name. `!duplicated(key$Admin2)` discarded
  # the second district of every duplicate-named pair -- the same defect that
  # was found in the dashboard's get_country_admin2().
  key <- H[H$country == ctry, c("Admin1", "Admin2")]
  key <- key[!duplicated(key[, c("Admin1", "Admin2")]), ]
  if (!nrow(key)) next
  hc  <- H[H$country == ctry, ]
  # Carry Admin1 so the benchmark's population weights key on the pair too;
  # without it a duplicate-named district takes the other one's population.
  pop <- data.frame(Admin1 = hc$Admin1, Admin2 = hc$Admin2,
                    pop = if ("ghs_pop" %in% names(hc)) hc$ghs_pop else NA_real_)

  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    suffix <- paste0(tolower(ctry), "_", ocn)
    svy <- rd(paste0("svy_admin2_", suffix)); od <- rd(paste0("outcome_data_", suffix))
    if (is.null(svy) || is.null(od) || nrow(svy) < 5) next
    rel <- admin2_reliability(svy, deff = 1.5, boot = 0)
    rmax <- rel$r_max %||% NA_real_
    message("\n=== ", ctry, " / ", ocn, " (r_max ", round(rmax, 3), ") ===")

    # Baseline: the same estimator fitted at admin-2, but scored OUT OF FOLD by
    # leave-one-region-out. This matters -- an in-sample admin-2 fit would be
    # compared against a genuinely out-of-sample downscaling arm, and would win
    # for that reason alone rather than on merit. Region-block folds also match
    # the spatial-block structure the downscaling arm is validated under.
    tr <- dplyr::inner_join(svy[, intersect(c("Admin1", "Admin2", "svy_prev", "n_svy"),
                                            names(svy)), drop = FALSE],
                            hc, by = admin2_join_by(svy, hc))
    tr <- tr[is.finite(tr$svy_prev), , drop = FALSE]
    if (nrow(tr) < 10 || dplyr::n_distinct(tr$Admin1) < 3) next
    Xtr <- as.matrix(tr[, COVS, drop = FALSE])
    oof <- rep(NA_real_, nrow(tr))
    for (r in unique(tr$Admin1)) {
      i <- which(tr$Admin1 == r)
      m <- .ds_fit(Xtr[-i, , drop = FALSE], tr$svy_prev[-i], k_screen = 20L)
      p <- .ds_predict(m, Xtr[i, , drop = FALSE])
      if (!is.null(p)) oof[i] <- p
    }
    p2 <- data.frame(Admin2 = tr$Admin2, Admin1 = tr$Admin1, pred = oof)
    p2 <- p2[is.finite(p2$pred), , drop = FALSE]
    if (nrow(p2) < 5) next
    rows[[length(rows) + 1L]] <- score(svy, p2, "pred", ctry, ocn,
                                       "admin-2 fit (LORO), unbenchmarked", rmax)

    nat <- national_design_based(od, cc, oc)
    if (!is.null(nat)) {
      pn <- benchmark_admin2_table(p2, "pred", nat$prev, pop)
      rows[[length(rows) + 1L]] <- score(svy, pn, "pred", ctry, ocn,
                                         "admin-2 fit + national benchmark", rmax)
    }

    a1t <- admin1_design_based(od, cc, oc)
    if (!is.null(a1t) && !is.null(nat)) {
      # Shrunk targets: protects against a thin region being treated as an exact
      # anchor, at the cost of compressing genuine between-region variation.
      a1s <- shrink_admin1_targets(a1t, nat$prev)
      pa <- benchmark_admin2_to_admin1(p2, "pred", a1s, national = nat$prev,
                                       admin1_map = key, pop = pop)
      rows[[length(rows) + 1L]] <- score(svy, pa, "pred", ctry, ocn,
                                         "admin-2 fit + ADMIN-1 benchmark (shrunk)", rmax)
      # Hard targets: which side of that trade-off wins is an empirical
      # question, and the admin-1 ceiling (r_max ~ 0.66) is high enough that it
      # is not obvious shrinkage is needed at all.
      a1h <- a1t; a1h$target <- a1h$prev
      ph <- benchmark_admin2_to_admin1(p2, "pred", a1h, national = nat$prev,
                                       admin1_map = key, pop = pop)
      rows[[length(rows) + 1L]] <- score(svy, ph, "pred", ctry, ocn,
                                         "admin-2 fit + ADMIN-1 benchmark (hard)", rmax)
    }
  }
}

# ── ARM B: pooled admin-1 downscaling, one model per outcome across countries ─
message("\n\n########## ARM B: pooled admin-1 -> admin-2 downscaling ##########")
outcomes <- unique(unlist(lapply(cfgs, function(c) names(c$outcomes))))
for (ocn in outcomes) {
  a1_list <- list(); a2_list <- list(); svy_list <- list()
  for (ctry in names(cfgs)) {
    cc <- cfgs[[ctry]]; if (!ocn %in% names(cc$outcomes)) next
    oc <- cc$outcomes[[ocn]]
    svy <- rd(paste0("svy_admin2_", tolower(ctry), "_", ocn))
    od  <- rd(paste0("outcome_data_", tolower(ctry), "_", ocn))
    if (is.null(svy) || is.null(od)) next
    hc <- H[H$country == ctry, ]
    a1t <- admin1_design_based(od, cc, oc); if (is.null(a1t)) next
    wts <- data.frame(Admin1 = hc$Admin1, Admin2 = hc$Admin2,
                      w = if ("ghs_pop" %in% names(hc)) hc$ghs_pop else 1)
    agg <- aggregate_covariates_to_admin1(hc[, c("Admin1", "Admin2", COVS)], wts)
    a1 <- dplyr::inner_join(a1t, agg, by = "Admin1")
    if (!nrow(a1)) next
    a1$country <- ctry
    a1_list[[ctry]] <- a1
    a2_list[[ctry]] <- hc
    svy_list[[ctry]] <- svy
  }
  if (length(a1_list) < 2) next
  A1 <- dplyr::bind_rows(a1_list); A2 <- dplyr::bind_rows(a2_list)
  message("\n=== ", ocn, ": ", nrow(A1), " pooled admin-1 units from ",
          length(a1_list), " countries ===")
  ds <- fit_downscale_admin1(A1, A2, covs = COVS, k_screen = 20L, min_units = 30L)
  if (is.null(ds)) next
  for (ctry in names(svy_list)) {
    pr <- ds$pred_admin2[ds$pred_admin2$country == ctry, ]
    rel <- admin2_reliability(svy_list[[ctry]], deff = 1.5, boot = 0)
    rows[[length(rows) + 1L]] <- score(svy_list[[ctry]], pr, "downscale_prev",
                                       ctry, ocn, "ADMIN-1 fit -> admin-2 (pooled)",
                                       rel$r_max %||% NA_real_)
  }
}

out <- dplyr::bind_rows(rows)
readr::write_csv(out, here("results", "tables", "admin1_arms.csv"))
message("\n-> results/tables/admin1_arms.csv (", nrow(out), " rows)")

print(as.data.frame(out %>% group_by(arm) %>%
  summarise(cells = dplyr::n(),
            median_r = round(stats::median(pearson_r, na.rm = TRUE), 3),
            median_mae_pp = round(stats::median(mae_pp, na.rm = TRUE), 2),
            median_abs_bias_pp = round(stats::median(abs(bias_pp), na.rm = TRUE), 2),
            .groups = "drop") %>% arrange(median_mae_pp)), row.names = FALSE)
