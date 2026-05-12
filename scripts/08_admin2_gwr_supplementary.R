# =============================================================================
# scripts/08_admin2_gwr_supplementary.R
#
# Supplementary spatial analysis: Moran's I on Admin2 SL residuals, plus
# Geographically Weighted Regression (GWR) of survey-weighted Admin2
# prevalence on aggregated proxy predictors. Loops over countries x outcomes.
#
# Policy (set after smoke-testing on Gambia/child_iron, n=17, R^2 +0.007):
#   * Moran's I residual diagnostic runs ALWAYS (cheap, defensible).
#   * GWR runs only when:
#       - country is large enough (n_admin2 >= MIN_N_GWR), AND
#       - residual Moran's I p_value <= MORANS_GATE_P (i.e. there is
#         spatial structure for GWR to capture).
#     Otherwise GWR is skipped and only Moran's I is written.
#   * Gambia (~17 effective Admin2s) is excluded from GWR by the size gate.
#
# Inputs
#   - results/shared/admin2_predictions_<country>_<outcome>.csv
#     (produced by scripts/06_admin2_predictions_map.R / export script)
#   - data/IPD/<country>/<country>_merged_dataset.rds (individual-level,
#     used to aggregate proxy predictors to Admin2)
#   - GADM level-2 polygons (cached in data/gadm/)
#
# Outputs
#   results/admin2/gwr/
#     <country>_<outcome>_morans_i.csv          (residual Moran's I)
#     <country>_<outcome>_gwr_summary.csv       (bandwidth, AIC, F-test)
#     <country>_<outcome>_gwr_local_coefs.csv   (per-Admin2 betas + t)
#     figures/<country>_<outcome>_local_coef_<var>.png
#     figures/<country>_<outcome>_local_r2.png
#
# Usage
#   Rscript scripts/08_admin2_gwr_supplementary.R [country] [outcome]
#   (no args -> all countries x outcomes)
#
# Proxy-only constraint: gw_-prefixed survey variables are excluded from
# candidate predictors (see R/gwr_analysis.R).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(geodata)
  library(viridis)
  library(patchwork)
  library(spdep)
  library(GWmodel)
  library(car)
})

source(here::here("R", "gwr_analysis.R"))


# ── Gating thresholds ─────────────────────────────────────────────────────────
MIN_N_GWR      <- 50    # require >=50 Admin2 with both survey + SL predictions
MORANS_GATE_P  <- 0.10  # require residual spatial autocorrelation p <= 0.10
FORCE_GWR      <- FALSE # set TRUE to bypass gates (e.g. for sensitivity runs)


# ── Parse args ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
arg_country <- if (length(args) >= 1) args[[1]] else NA_character_
arg_outcome <- if (length(args) >= 2) args[[2]] else NA_character_
if (length(args) >= 3 && tolower(args[[3]]) == "force") FORCE_GWR <- TRUE


# ── Country / outcome registry (mirrors R/config.R, kept local for portability)
country_reg <- list(
  Gambia       = list(slug = "gambia",       gadm = "GMB",
                      data = here::here("data", "IPD", "Gambia",
                                        "Gambia_merged_dataset.rds")),
  Ghana        = list(slug = "ghana",        gadm = "GHA",
                      data = here::here("data", "IPD", "Ghana",
                                        "Ghana_merged_dataset.rds")),
  `Sierra Leone` = list(slug = "sierraleone", gadm = "SLE",
                        data = here::here("data", "IPD", "Sierra Leone",
                                          "Sierra_Leone_merged_dataset.rds")),
  Malawi       = list(slug = "malawi",       gadm = "MWI",
                      data = here::here("data", "IPD", "Malawi",
                                        "Malawi_merged_dataset.rds"))
)

outcomes_by_country <- list(
  Gambia       = c("child_vitA", "women_vitA", "child_iron", "women_iron"),
  Ghana        = c("child_vitA", "women_vitA", "child_iron", "women_iron",
                   "women_b12", "women_folate"),
  `Sierra Leone` = c("child_vitA", "women_vitA", "child_iron", "women_iron",
                   "women_b12", "women_folate"),
  Malawi       = c("child_vitA", "women_vitA", "child_iron", "women_iron",
                   "women_b12", "women_folate", "child_zinc", "women_zinc")
)


# ── Output paths ──────────────────────────────────────────────────────────────
out_root <- here::here("results", "admin2", "gwr")
out_fig  <- file.path(out_root, "figures")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)


# ── Helper: enumerate proxy predictor candidates from individual data ─────────
enumerate_proxy_candidates <- function(ind_df) {
  cn <- colnames(ind_df)
  # Drop survey-derived (gw_) per project memory; admin labels; weights/IDs.
  drop_pat <- "^(gw_|hh_|cluster|psu|strata|weight|svy_)"
  drop_meta <- c("Admin1", "Admin2", "country", "outcome")
  cand <- cn[!grepl(drop_pat, cn, ignore.case = TRUE)]
  cand <- setdiff(cand, drop_meta)
  cand <- cand[vapply(cand, function(v) is.numeric(ind_df[[v]]), logical(1))]
  cand
}


# ── Main loop ────────────────────────────────────────────────────────────────
process_one <- function(country, outcome) {
  reg <- country_reg[[country]]
  pred_csv <- here::here("results", "shared",
                         sprintf("admin2_predictions_%s_%s.csv",
                                 reg$slug, outcome))
  if (!file.exists(pred_csv)) {
    cat(sprintf("[skip] missing %s\n", pred_csv))
    return(invisible(NULL))
  }

  cat(sprintf("\n[gwr] %s | %s\n", country, outcome))

  # ── Stage 1: Moran's I (always) ──────────────────────────────────────────
  morans_only <- tryCatch(
    run_admin2_gwr_pipeline(
      admin2_pred_csv = pred_csv,
      ind_data        = NULL,
      gadm_code       = reg$gadm,
      skip_gwr        = TRUE
    ),
    error = function(e) {
      cat(sprintf("  ERROR (Moran's): %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(morans_only)) return(invisible(NULL))

  m <- morans_only$morans
  if (!is.na(m$morans_i)) {
    cat(sprintf("  Moran's I: %.3f (p = %.3f, n = %d)\n",
                m$morans_i, m$p_value, m$n))
  } else {
    cat(sprintf("  Moran's I: NA (%s)\n", m$note))
  }

  # ── Stage 2: gate GWR ────────────────────────────────────────────────────
  gate_n     <- !is.na(m$n)        && m$n        >= MIN_N_GWR
  gate_pval  <- !is.na(m$p_value)  && m$p_value  <= MORANS_GATE_P
  run_gwr    <- FORCE_GWR || (gate_n && gate_pval)

  if (!run_gwr) {
    reasons <- character(0)
    if (!gate_n)    reasons <- c(reasons, sprintf("n=%s < %d", m$n, MIN_N_GWR))
    if (!gate_pval) reasons <- c(reasons,
      sprintf("Moran's p=%s > %.2f (no spatial structure to capture)",
              format(round(m$p_value, 3)), MORANS_GATE_P))
    cat(sprintf("  GWR skipped: %s\n", paste(reasons, collapse = "; ")))
    res <- morans_only
  } else {
    ind <- if (file.exists(reg$data)) {
      tryCatch(readRDS(reg$data), error = function(e) NULL)
    } else NULL
    candidates <- if (!is.null(ind)) enumerate_proxy_candidates(ind) else character(0)
    if (length(candidates) > 0) {
      cat(sprintf("  %d candidate proxy predictors\n", length(candidates)))
    }
    res <- tryCatch(
      run_admin2_gwr_pipeline(
        admin2_pred_csv      = pred_csv,
        ind_data             = ind,
        gadm_code            = reg$gadm,
        candidate_predictors = candidates,
        k_neighbors          = 5,
        top_k                = 6
      ),
      error = function(e) {
        cat(sprintf("  ERROR (GWR): %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) res <- morans_only
  }

  # Write Moran's I
  morans <- cbind(country = country, outcome = outcome, res$morans)
  morans_path <- file.path(out_root,
    sprintf("%s_%s_morans_i.csv", reg$slug, outcome))
  utils::write.csv(morans, morans_path, row.names = FALSE)
  cat(sprintf("  wrote %s\n", morans_path))

  if (is.null(res$gwr_summary)) {
    cat("  GWR not run (insufficient predictors or individual data).\n")
    return(invisible(NULL))
  }

  # Write GWR tables
  gwr_sum <- cbind(country = country, outcome = outcome, res$gwr_summary)
  utils::write.csv(gwr_sum,
    file.path(out_root, sprintf("%s_%s_gwr_summary.csv", reg$slug, outcome)),
    row.names = FALSE)

  utils::write.csv(res$gwr_local_coefs,
    file.path(out_root, sprintf("%s_%s_gwr_local_coefs.csv", reg$slug, outcome)),
    row.names = FALSE)

  # Choropleths: local R^2 + per-predictor local beta
  polys <- res$polys %>%
    dplyr::left_join(res$gwr_local_coefs, by = c("Admin1", "Admin2"))

  p_r2 <- ggplot(polys) +
    geom_sf(aes(fill = local_R2), colour = "white", linewidth = 0.1) +
    scale_fill_viridis_c(na.value = "grey90", name = expression(R^2)) +
    theme_void() +
    ggtitle(sprintf("%s — %s — Local R² (GWR)", country, outcome))
  ggsave(file.path(out_fig,
                   sprintf("%s_%s_local_r2.png", reg$slug, outcome)),
         p_r2, width = 7, height = 5, dpi = 150)

  for (v in res$selected_predictors) {
    bcol <- paste0("beta_", v)
    pcol <- paste0("padj_", v)
    if (!bcol %in% colnames(polys)) next
    polys$.signif <- if (pcol %in% colnames(polys)) {
      ifelse(!is.na(polys[[pcol]]) & polys[[pcol]] < 0.05, "p<0.05", "ns")
    } else "n/a"

    p_b <- ggplot(polys) +
      geom_sf(aes(fill = .data[[bcol]]), colour = "white", linewidth = 0.1) +
      scale_fill_gradient2(name = expression(beta), na.value = "grey90") +
      geom_sf(data = polys[polys$.signif == "p<0.05", ],
              fill = NA, colour = "black", linewidth = 0.4) +
      theme_void() +
      ggtitle(sprintf("%s — %s — Local β for %s", country, outcome, v),
              subtitle = "outlined: BH-adjusted p < 0.05")
    ggsave(file.path(out_fig,
                     sprintf("%s_%s_local_coef_%s.png", reg$slug, outcome, v)),
           p_b, width = 7, height = 5, dpi = 150)
  }

  cat(sprintf("  GWR done: bw=%s, AIC delta=%.1f, R2 global=%.3f -> GWR=%.3f\n",
              format(res$gwr_summary$bandwidth),
              res$gwr_summary$aic_delta,
              res$gwr_summary$r2_global,
              res$gwr_summary$r2_gwr))
  invisible(res)
}


# ── Build job list ────────────────────────────────────────────────────────────
jobs <- list()
countries <- if (!is.na(arg_country)) arg_country else names(country_reg)
for (cn in countries) {
  outs <- if (!is.na(arg_outcome)) arg_outcome else outcomes_by_country[[cn]]
  for (oc in outs) {
    jobs[[length(jobs) + 1]] <- list(country = cn, outcome = oc)
  }
}

cat(sprintf("[gwr] %d job(s) queued.\n", length(jobs)))
for (j in jobs) process_one(j$country, j$outcome)

cat("\n[gwr] done.\n")
